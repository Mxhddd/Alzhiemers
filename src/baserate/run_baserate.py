"""Module 2 entry point: fit variants A, B, C to autopsy stage counts.

    python -m src.baserate.run_baserate                       (real data; gated)
    python -m src.baserate.run_baserate --dry-run-simulated   (simulated; no gates)

Real-data mode refuses to start unless Module 2's validation passed on the
current code and prereg/baserate.md is frozen. The dry run never reads
data/: it simulates counts with the real structure (stage scheme, decade
bands) from a known parameter set, writes to results/dryrun/ with every
output stamped SIMULATED, and exists to test the pipeline end to end.

What Module 2 tests (prereg/baserate.md, revised 2026-10-13):
  * criterion 2 — which variant (A: LC first, B: independent cortical
    entry, C: multi-hit) best explains the stage distribution, by k-fold
    cross-validation over whole age bands;
  * falsification — each variant's predicted incidence of the dementia-
    threshold stage against published incidence it never saw.
It does NOT test amyloid gating (criterion 1 was dropped): the stage counts
carry no per-brain amyloid status, so the population amyloid curve A(t) is
just a function of age. The age-adjusted amyloid HR from variant A is
reported as "descriptive, not a test of gating"; the brief's literal
variant A is a sensitivity analysis labelled "confounded with age".
``require_gating_disabled`` refuses to run if the config ever re-enables
gating as a test, and the validation suite checks that it does.

Every stage (primary fit per variant, CV per variant, each sensitivity
refit) is checkpointed under results/dryrun/checkpoints/ (or
results/checkpoints/baserate/) and the tables are rewritten as soon as it
finishes. A run stopped part-way (e.g. by a job time limit) loses at most
the stage in progress; rerunning resumes. Checkpoints are keyed on the
config, the module's code fingerprint and the run settings, so stale ones
are never reused. ``--max-minutes`` stops cleanly (exit code 3) before
starting a stage once the budget is spent.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import sys
import time
from dataclasses import dataclass, field, replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.baserate.data import (StageCounts, jansen_bands, load_amyloid_phase_counts, load_braak_counts,
                               load_incidence)
from src.baserate.fit import fit_bayesian, fit_mle, group_probabilities
from src.baserate.jaxlik import make_incidence_batch
from src.baserate.model import AGated, AGatedAge, BIndependent, ChainSpec, CMultiHit, StepFunction
from src.baserate.predict import compare_models, predict_incidence, stage_prevalence_table
from src.baserate.solve import error_bound, stage_probabilities_ode
from src.common.config import load_config
from src.common.paths import CONFIG_DIR, FIGURES_DIR, REPO_ROOT, RESULTS_DIR, TABLES_DIR
from src.common.plotting import INK, apply_theme, save_figure, series_colors
from src.common.prereg import require_preregistration
from src.common.provenance import sha256_file
from src.common.runlog import start_run
from src.common.validation import code_fingerprint, require_validated

CONFIG = CONFIG_DIR / "baserate.yaml"
DESCRIPTIVE = "descriptive, not a test of gating"
CONFOUNDED = "confounded with age"
INCOMPLETE_EXIT = 3


class GatingReenabledError(RuntimeError):
    pass


def require_gating_disabled(cfg) -> None:
    """Criterion 1 (amyloid gating) was dropped on 2026-10-13. Refuse to run
    if the configuration ever turns it back on."""
    flag = cfg.get("amyloid.gating_is_a_test")
    if flag is not False:
        raise GatingReenabledError(
            f"amyloid.gating_is_a_test is {flag!r}: Module 2 cannot test amyloid gating (stage-by-age counts "
            "carry no per-brain amyloid status; prereg/baserate.md change log 2026-10-13). Re-enabling it "
            "needs per-person amyloid x tau-stage data, a new model and a new pre-registration.")


@dataclass(frozen=True)
class Scenario:
    name: str
    threshold: str
    amyloid: str          # "primary" | "sensitivity"
    eval_mode: str
    c_hits: int
    open_age: float
    refit: tuple[str, ...]      # variants refitted in this scenario (others reuse the primary fit)
    variants: tuple[str, ...] = field(default=())   # variants reported in this scenario (default: primary ones)


@dataclass
class Inputs:
    data: list[StageCounts]
    amyloid: dict[str, StepFunction]
    incidence: pd.DataFrame          # band, eval_age, observed_min, observed_max (open band at primary age)
    incidence_open_lower: float | None
    files: list[Path] = field(default_factory=list)   # real input files (checkpoint key)


# ---------------------------------------------------------------- inputs ----
BRAAK_TOTALS = {  # quantity in braak2011_abstract_totals.csv -> stages it sums over
    "total brains examined (ages 1-100)": None,
    "AT8-immunonegative after restaging (no abnormal tau)": ["0"],
    "subcortical stages a-c (LC-predominant tau; no cortical tau)": ["a", "b", "c"],
    "stage 1a (transentorhinal neurites)": ["1a"],
    "stage 1b (transentorhinal pyramidal-cell pretangles)": ["1b"],
}


def check_braak_totals(data: list[StageCounts], stages: list[str], run, amyloid=None) -> None:
    """The extracted stage table must reproduce the abstract's totals exactly,
    the amyloid table must give the abstract's plaque total, and every band
    total must agree between the two tables (pre-registration §3); any
    mismatch stops the analysis. `amyloid`: load_amyloid_phase_counts rows
    with min_phase = 1."""
    from src.common.extracted import load_extracted

    tot = load_extracted(REPO_ROOT / "data/extracted/braak2011_abstract_totals.csv",
                         provided_by="extracted from the Braak 2011 abstract", run=run)
    counts = np.sum([g.counts for g in data], axis=0)
    problems = []
    for q, st in BRAAK_TOTALS.items():
        want = int(float(tot.loc[tot["quantity"] == q, "value"].iloc[0]))
        got = int(counts.sum()) if st is None else int(sum(counts[stages.index(x)] for x in st))
        if got != want:
            problems.append(f"{q}: table {got} vs abstract {want}")
    if amyloid is not None:
        plaques = int(float(tot.loc[(tot["quantity"] == "cases with beta-amyloid plaques")
                                    & (tot["unit"] == "cases"), "value"].iloc[0]))
        if sum(r[3] for r in amyloid) != plaques:
            problems.append(f"plaque cases: amyloid table {sum(r[3] for r in amyloid)} vs abstract {plaques}")
        band_n = {(g.age_lower, g.age_upper): g.n for g in data}
        for lo, hi, _, _, n in amyloid:
            if band_n.get((lo, hi)) != n:
                problems.append(f"band {lo}-{hi}: stage table {band_n.get((lo, hi))} vs amyloid table {n} brains")
    if problems:
        raise ValueError("Braak 2011 stage table does not reproduce the abstract: " + "; ".join(problems))
    run.info("Braak 2011 stage table reproduces all abstract totals")


def real_inputs(cfg, run) -> Inputs:
    stages = [str(s) for s in cfg.get("data.stage_order")]
    inc = cfg.get("data.upper_inclusive")
    files = [REPO_ROOT / cfg.get(k) for k in ("data.braak_counts", "data.amyloid_prevalence",
                                              "data.amyloid_prevalence_sensitivity", "data.dementia_incidence")]
    data = load_braak_counts(files[0], stages, upper_inclusive=inc, run=run)
    min_phase = int(cfg.get("data.amyloid_positive_min_phase"))
    amy_rows = load_amyloid_phase_counts(files[1], upper_inclusive=inc, min_phase=min_phase, run=run)
    check_braak_totals(data, stages, run, amyloid=load_amyloid_phase_counts(files[1], upper_inclusive=inc,
                                                                            min_phase=1, run=None))
    primary = StepFunction([(lo, hi, p) for lo, hi, p, _, _ in amy_rows])
    sens = splice_amyloid(primary, jansen_bands(files[2], run=run))
    incidence = load_incidence(files[3], open_band_age=cfg.get("data.open_band_evaluation_age"), run=run)
    return Inputs(data, {"primary": primary, "sensitivity": sens}, incidence, _open_lower(files[3]),
                  files + [REPO_ROOT / "data/extracted/braak2011_abstract_totals.csv"])


def splice_amyloid(primary: StepFunction, inner: list[tuple[float, float, float]]) -> StepFunction:
    """Use `inner` bands where they exist (Jansen 2015 covers 47.5-92.5) and the
    primary curve outside that range, cutting primary bands that straddle it."""
    lo, hi = min(b[0] for b in inner), max(b[1] for b in inner)
    out = list(inner)
    for a, b, v in primary.bands:
        if b <= lo or a >= hi:
            out.append((a, b, v))
            continue
        if a < lo:
            out.append((a, lo, v))
        if b > hi:
            out.append((hi, b, v))
    return StepFunction(sorted(out))


def _open_lower(path) -> float | None:
    df = pd.read_csv(path, dtype=str)
    o = df[df["age_upper"] == "open"]
    return float(o["age_lower"].iloc[0]) if len(o) else None


SIM_TRUTH = {  # SIMULATION SETTINGS (not data): log-rate at age 50 and log-slope per decade per edge
    "alpha[0]": -2.5, "beta[0]": 0.0, "alpha[1]": -2.5, "beta[1]": 0.0, "alpha[2]": -3.0, "beta[2]": 0.1,
    # LC -> cortex: rises with age (0.1 per decade, like its neighbours) and has NO amyloid effect,
    # so the dry run shows the literal variant A manufacturing a hazard ratio from age alone.
    "log_mu0": -4.0, "beta_mu": 0.1, "log_hr_amyloid": 0.0,
    "alpha[4]": -2.5, "beta[4]": 0.1, "alpha[5]": -3.0, "beta[5]": 0.2, "alpha[6]": -3.0, "beta[6]": 0.3,
    "alpha[7]": -3.5, "beta[7]": 0.5, "alpha[8]": -4.0, "beta[8]": 0.6, "alpha[9]": -4.5, "beta[9]": 0.8,
    "alpha[10]": -4.0, "beta[10]": 0.5}
SIM_AMYLOID = [(1, 41, 0.0), (41, 51, 0.05), (51, 61, 0.15), (61, 71, 0.3), (71, 81, 0.45), (81, 91, 0.6),
               (91, 101, 0.7)]
SIM_AMYLOID_SENS = [(1, 47.5, 0.0), (47.5, 67.5, 0.1), (67.5, 87.5, 0.3), (87.5, 101, 0.45)]


def simulated_inputs(cfg, spec: ChainSpec, step: float, start: float, centre: float, scale: float,
                     truth_settings: dict | None = None, seed: int | None = None) -> tuple[Inputs, np.ndarray]:
    """Synthetic world with the real structure. All numbers here are
    SIMULATION SETTINGS, not data, and never leave the dry run."""
    from sim.baserate_sim import simulate_counts

    amy, amy_sens = StepFunction(SIM_AMYLOID), StepFunction(SIM_AMYLOID_SENS)
    truth_model = AGatedAge(spec, amy, centre, scale)
    settings = dict(SIM_TRUTH, **(truth_settings or {}))
    theta = np.array([settings[n] for n in truth_model.param_names])
    bands = [(lo, lo + 10) for lo in range(1, 101, 10)]
    data = simulate_counts(truth_model, theta, bands, 230, step=step, start_age=start,
                           seed=cfg.get("seed") if seed is None else seed)
    ages = np.array([67.5, 72.5, 77.5, 82.5, float(cfg.get("data.open_band_evaluation_age"))])
    thr = list(spec.stages).index(str(cfg.get("data.dementia_threshold_stage")))
    true_inc = 1000 * predict_incidence(truth_model, theta, ages, threshold=thr, step=step, start_age=start)
    incidence = pd.DataFrame({"band": ["65-69", "70-74", "75-79", "80-84", "85-open"], "eval_age": ages,
                              "n_estimates": 4, "observed_min": 0.8 * true_inc, "observed_max": 1.2 * true_inc})
    return Inputs(data, {"primary": amy, "sensitivity": amy_sens}, incidence, 85.0), theta


# ------------------------------------------------------------- analysis ----
def build_model(spec, name: str, sc: Scenario, inputs: Inputs, centre, scale):
    if name == "A_gated_age":
        return AGatedAge(spec, inputs.amyloid[sc.amyloid], centre, scale)
    if name == "A_gated_literal":
        return AGated(spec, inputs.amyloid[sc.amyloid], centre, scale)
    if name == "B_independent":
        return BIndependent(spec, centre, scale)
    if name == "C_multihit":
        return CMultiHit(spec, sc.c_hits, centre, scale)
    raise ValueError(f"unknown variant {name}")


def build_models(spec, sc: Scenario, inputs: Inputs, centre, scale, names):
    return {v: build_model(spec, v, sc, inputs, centre, scale) for v in names}


def audit(model, theta, data, *, step, start, eval_mode, n_points) -> tuple[float, float]:
    """max |p_h - p_ODE| over evaluation ages/stages, and |LL_h - LL_ODE|."""
    P_h = group_probabilities(model, theta, data, step=step, start_age=start, eval_mode=eval_mode, n_points=n_points)
    pts = [g.eval_points(eval_mode, n_points) for g in data]
    ages = np.concatenate([a for a, _ in pts])
    P_ode_all = stage_probabilities_ode(model, theta, ages, start_age=start)
    P_ode, i = [], 0
    for a, w in pts:
        P_ode.append(w @ P_ode_all[i:i + len(a)])
        i += len(a)
    P_ode = np.array(P_ode)

    def ll(P):
        return sum(float(np.sum(g.counts[g.counts > 0] * np.log(np.maximum(p[g.counts > 0], 1e-300))))
                   for g, p in zip(data, P))

    return float(np.abs(P_h - P_ode).max()), abs(ll(P_h) - ll(P_ode))


def fit_audited(model, data, *, cfg, step, start, eval_mode, n_points, seed, run, warm=None, cache=None, tag=None):
    """MLE, then the pre-registered step audit; halve h and refit on failure.
    `warm` (e.g. the primary fit) is added as one extra start. `cache`
    (a SubCheckpoints) saves each attempt under `tag`, so a long run resumes
    per attempt; every attempt is an independent fit, so a cached one is the
    same computation as a fresh one."""
    for _attempt in range(4):
        label = f"{tag}/step{step:g}"
        rec = cache.get(label) if cache is not None and tag else None
        if rec is None:
            fit = fit_mle(model, data, step=step, start_age=start, cfg_mle=cfg.get("mle"), seed=seed,
                          eval_mode=eval_mode, n_points=n_points, extra_starts=None if warm is None else [warm])
            dp, dll = audit(model, fit.theta, data, step=step, start=start, eval_mode=eval_mode, n_points=n_points)
            end = max(a.max() for a, _ in (g.eval_points(eval_mode, n_points) for g in data))
            bound = step**2 * error_bound(model, fit.theta, start_age=start, end_age=end)
            if cache is not None and tag:
                cache.put(label, {"fit": fit, "dp": dp, "dll": dll, "bound": bound})
        else:
            fit, dp, dll, bound = rec["fit"], rec["dp"], rec["dll"], rec["bound"]
        ok = (dp < cfg.get("solver.max_prob_error") and dll < cfg.get("solver.max_loglik_error")
              and bound <= cfg.get("solver.max_prob_error"))
        run.info(f"{model.name}: step {step:g} y: bound h^2*Phi = {bound:.2e}, Radau audit max|dp| = {dp:.2e}, "
                 f"|dLL| = {dll:.2e} -> {'ok' if ok else 'FAIL: halving the step'}")
        if ok:
            return fit, step, dp, dll, bound
        step /= 2
    raise RuntimeError(f"{model.name}: step audit still failing after 3 halvings")


def falsification(model, theta, inputs: Inputs, *, threshold: int, open_age: float, step, start, tol,
                  draws: np.ndarray | None = None, level: float = 0.95):
    """Pre-registered incidence test (prereg §6, revised 2026-10-13).

    For each CFAS band the published envelope is [min, max] of the four
    rates; the acceptance region is that envelope widened by the factor
    e^tol (= 2) on both sides. Predicted incidence is computed at the fitted
    parameters (point) and for each parameter draw (posterior draws for the
    primary fit, Laplace draws N(theta_hat, cov) for MLE-only fits).
      * A band FAILS if the whole `level` interval of the predictions lies
        outside the acceptance region (distance of its nearest end > tol);
      * a band is UNINFORMATIVE if its interval is wider than the whole
        acceptance region (it could not have failed);
      * verdict: FAIL if any band fails; else UNINFORMATIVE if every band is
        uninformative; else PASS.
    The point-estimate rule of the first draft (fail if the point prediction
    is > tol outside) is reported alongside as `point_rule_verdict`. Without
    draws (e.g. a non-positive-definite Hessian) the interval collapses to
    the point and this is recorded."""
    inc = inputs.incidence.copy()
    if inputs.incidence_open_lower is not None:
        inc.loc[inc["band"].str.endswith("open"), "eval_age"] = open_age
    ages = inc["eval_age"].to_numpy()

    def pred(th):
        return 1000 * predict_incidence(model, th, ages, threshold=threshold, step=step, start_age=start)

    point = pred(np.asarray(theta, float))
    sims, failed_draws = [], 0
    if draws is not None and len(draws):
        # compiled + vectorised over draws (equal to `pred` to 1e-10; unit-tested)
        batch = make_incidence_batch(model, ages, threshold=threshold, step=step, start_age=start)
        V = 1000 * np.asarray(batch(np.asarray(draws, float)))
        ok = np.all(np.isfinite(V) & (V >= 0), axis=1)
        sims, failed_draws = list(V[ok]), int((~ok).sum())
    if sims:
        q_lo, q_hi = np.quantile(np.array(sims), [(1 - level) / 2, (1 + level) / 2], axis=0)
    else:
        q_lo, q_hi = point.copy(), point.copy()
    lo, hi = inc["observed_min"].to_numpy(), inc["observed_max"].to_numpy()

    def outside(a_lo, a_hi):
        return np.where(a_lo > hi, np.log(a_lo / hi), np.where(a_hi < lo, np.log(lo / a_hi), 0.0))

    d_point, d_int = outside(point, point), outside(q_lo, q_hi)
    region = np.log(hi * np.exp(tol)) - np.log(lo / np.exp(tol))
    uninformative = (np.log(q_hi) - np.log(q_lo)) > region
    inc["predicted_per_1000py"] = point
    inc["pred_lower"], inc["pred_upper"] = q_lo, q_hi
    inc["interval_from"] = f"{len(sims)} parameter draws" if sims else "point estimate only (no draws)"
    inc["failed_draws"] = failed_draws
    inc["log_distance_outside_envelope"] = d_point
    inc["interval_log_distance_outside_envelope"] = d_int
    inc["band_uninformative"] = uninformative
    inc["point_rule_verdict"] = "PASS" if np.all(d_point <= tol) else "FAIL"
    verdict = ("FAIL" if np.any(d_int > tol) else "UNINFORMATIVE" if np.all(uninformative) else "PASS")
    return inc, verdict


def laplace_draws(fit_theta, cov, n: int, seed: int) -> np.ndarray | None:
    if cov is None:
        return None
    return np.random.default_rng(seed).multivariate_normal(np.asarray(fit_theta, float), np.asarray(cov, float), size=n)


def cv_verdict_3way(cv_df: pd.DataFrame, margin: float, a_name: str) -> str:
    """Pre-registered three-way reading of the CV comparison (prereg §2,
    criterion 2): A = the LC-first variant `a_name`."""
    ll = dict(zip(cv_df["model"], cv_df["cv_heldout_loglik"]))
    a = ll[a_name]
    others = {k: v for k, v in ll.items() if k != a_name}
    if all(a - v > margin for v in others.values()):
        return "supports (A beats B and C by more than the margin)"
    if any(v - a > margin for v in others.values()):
        return "AGAINST (" + ", ".join(k for k, v in others.items() if v - a > margin) + " beats A by more than the margin)"
    return "inconclusive (A not separated from the alternatives by the margin)"


def expectation_comparison(interval: tuple[float, float], expected: float, exp_lo: float, exp_hi: float, *,
                           on_bound: bool = False, max_width: float | None = None) -> str:
    """Descriptive comparison of a fitted log HR interval with the published
    expectation (manuscript/derivations/module2_gamma.md). Not a test.
    An estimate on the parameter bound, or an interval wider than `max_width`
    on the log scale (config gamma_check.uninformative_log_width), is labelled
    as such instead of being compared: "published value inside the interval"
    would otherwise read as agreement when the data say nothing."""
    lo, hi = interval
    if on_bound:
        return "not identified: estimate on the parameter bound"
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return "not identified: no finite standard error"
    if max_width is not None and hi - lo > max_width:
        return f"uninformative: interval spans more than a {np.exp(max_width):.0f}-fold range of HR"
    if lo <= expected <= hi:
        return "published value inside the interval"
    if hi < exp_lo or lo > exp_hi:
        return "intervals do not overlap"
    return "published value outside the interval; intervals overlap"


# --------------------------------------------------------- checkpoints ----
def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return [_jsonable(v) for v in x.tolist()]
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


class Checkpoints:
    """One JSON file per finished stage, valid only for the key it was written
    under (config + code fingerprint + run settings [+ input file hashes])."""

    def __init__(self, directory: Path, key: str, run):
        self.dir = Path(directory)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.key = key
        self.run = run

    def _path(self, stage: str) -> Path:
        return self.dir / (stage.replace("/", "__") + ".json")

    def get(self, stage: str):
        p = self._path(stage)
        if not p.exists():
            return None
        rec = json.loads(p.read_text(encoding="utf-8"))
        if rec.get("key") != self.key:
            self.run.warning(f"checkpoint {p.name} was written under different code/config/settings: ignored")
            return None
        self.run.info(f"stage {stage}: reusing checkpoint {p.name} (finished {rec.get('finished_utc')})")
        return rec["data"]

    def put(self, stage: str, data: dict) -> dict:
        """Write atomically (a killed run never leaves a half-written file) and
        return the data exactly as a later run will read it back."""
        rec = {"key": self.key, "stage": stage, "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "data": _jsonable(data)}
        text = json.dumps(rec, indent=1)
        tmp = self._path(stage).with_suffix(".tmp")
        tmp.write_text(text, encoding="utf-8")
        tmp.replace(self._path(stage))
        return json.loads(text)["data"]


class SubCheckpoints:
    """Finer-grained saves inside a stage (a step-size attempt, a NUTS chain,
    a cross-validation fold), so a stage longer than the machine's lifetime
    still makes progress. Same key rule as Checkpoints; pickled because the
    items are fit objects and arrays."""

    def __init__(self, directory: Path, key: str, run):
        self.dir = Path(directory)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.key = key
        self.run = run

    def _path(self, label: str) -> Path:
        return self.dir / (label.replace("/", "__") + ".pkl")

    def get(self, label: str):
        p = self._path(label)
        if not p.exists():
            return None
        with p.open("rb") as fh:
            rec = pickle.load(fh)
        if rec.get("key") != self.key:
            self.run.warning(f"sub-checkpoint {p.name} was written under different code/config/settings: ignored")
            return None
        self.run.info(f"reusing sub-checkpoint {label}")
        return rec["data"]

    def put(self, label: str, data) -> None:
        tmp = self._path(label).with_suffix(".tmp")
        with tmp.open("wb") as fh:
            pickle.dump({"key": self.key, "label": label, "data": data,
                         "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, fh)
        tmp.replace(self._path(label))
        self.run.info(f"sub-checkpoint written: {label}")

    def pair(self, prefix: str):
        """(get, put) functions for code that takes a plain callable pair."""
        return (lambda i: self.get(f"{prefix}/{i}"), lambda i, d: self.put(f"{prefix}/{i}", d))


def checkpoint_key(cfg_path: Path, settings: dict, input_files: list[Path]) -> str:
    h = hashlib.sha256()
    h.update(sha256_file(cfg_path).encode())
    h.update(code_fingerprint("baserate").encode())
    h.update(sha256_file(REPO_ROOT / "sim" / "baserate_sim.py").encode())
    h.update(json.dumps(settings, sort_keys=True).encode())
    for f in sorted(input_files):
        h.update(f"{f.name}:{sha256_file(f)}".encode())
    return h.hexdigest()


class OutOfTime(Exception):
    pass


# ----------------------------------------------------------------- main ----
def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run-simulated", action="store_true")
    ap.add_argument("--dry-run-nuts", type=int, nargs=3, metavar=("TUNE", "DRAWS", "CHAINS"),
                    help="reduced NUTS settings for timing the dry run")
    ap.add_argument("--dry-run-mle-starts", type=int, help="fewer random MLE starts in the dry run")
    ap.add_argument("--max-minutes", type=float, help="do not START a new stage after this many minutes")
    args = ap.parse_args(argv)
    dry = args.dry_run_simulated
    cfg = load_config(CONFIG)
    require_gating_disabled(cfg)
    prereg = None
    if (args.dry_run_nuts or args.dry_run_mle_starts) and not dry:
        raise SystemExit("reduced settings are only allowed with --dry-run-simulated")
    if not dry:
        require_validated("baserate")
        prereg = require_preregistration("baserate")
    base = RESULTS_DIR / "dryrun" if dry else RESULTS_DIR
    tables = base / "tables" if dry else TABLES_DIR
    figures = base / "figures" if dry else FIGURES_DIR
    ck_dir = RESULTS_DIR / "dryrun" / "checkpoints" if dry else RESULTS_DIR / "checkpoints" / "baserate"
    t_start = time.time()

    with start_run(module="baserate_dryrun" if dry else "baserate", config_path=CONFIG, seed=cfg.get("seed"),
                   packages=["numpy", "scipy", "pandas", "numpyro", "arviz", "jax"], preregistration=prereg) as run:
        stages = [str(s) for s in cfg.get("data.stage_order")]
        spec = ChainSpec.from_labels(stages, cfg.get("data.first_cortical_stage"))
        step0, start = cfg.get("solver.step_years"), cfg.get("solver.start_age")
        centre, scale = cfg.get("solver.age_centre"), cfg.get("solver.age_scale")
        npts = cfg.get("solver.band_average_points")
        seed = cfg.get("seed")
        if dry:
            run.mark_simulated()
            inputs, truth = simulated_inputs(cfg, spec, step0, start, centre, scale)
        else:
            inputs, truth = real_inputs(cfg, run), None
        run.info(f"{len(inputs.data)} age-band groups, {sum(g.n for g in inputs.data)} brains, {len(stages)} stages")
        bayes = dict(cfg.get("bayes"))
        if dry and args.dry_run_nuts:
            bayes.update(tune=args.dry_run_nuts[0], draws=args.dry_run_nuts[1], chains=args.dry_run_nuts[2])
            run.warning(f"dry run: NUTS reduced to tune={bayes['tune']}, draws={bayes['draws']}, chains={bayes['chains']}")
        if dry and args.dry_run_mle_starts:
            cfg._data["mle"]["n_starts"] = args.dry_run_mle_starts
            cv_r = min(args.dry_run_mle_starts, cfg.get("cross_validation.n_random_starts"))
            cfg._data["cross_validation"]["n_random_starts"] = cv_r
            run.warning(f"dry run: MLE random starts reduced to {args.dry_run_mle_starts} ({cv_r} per CV fold)")
        ck = Checkpoints(ck_dir, checkpoint_key(CONFIG, {"dry": dry, "nuts": args.dry_run_nuts,
                                                         "mle_starts": args.dry_run_mle_starts},
                                                inputs.files), run)
        sub = SubCheckpoints(ck_dir / "sub", ck.key, run)

        variants = list(cfg.get("variants.fit"))
        a_name = cfg.get("variants.lc_first")
        literal = cfg.get("variants.literal_a_sensitivity")
        primary = Scenario("primary", cfg.get("data.dementia_threshold_stage"), "primary", cfg.get("solver.evaluate_at"),
                           cfg.get("variants.c_hits"), cfg.get("data.open_band_evaluation_age"), tuple(variants),
                           tuple(variants))
        scenarios = [
            replace(primary, name="threshold_IV", threshold=cfg.get("data.dementia_threshold_stage_sensitivity"), refit=()),
            replace(primary, name="open_band_87.5", open_age=cfg.get("data.open_band_evaluation_age_sensitivity"), refit=()),
            replace(primary, name="band_average", eval_mode=cfg.get("solver.evaluate_at_sensitivity"), refit=tuple(variants)),
            replace(primary, name="c_hits_3", c_hits=cfg.get("variants.c_hits_sensitivity"), refit=("C_multihit",),
                    variants=("C_multihit",)),
            replace(primary, name="amyloid_jansen2015", amyloid="sensitivity", refit=(a_name,), variants=(a_name,)),
            replace(primary, name="literal_A_confounded_with_age", refit=(literal,), variants=(literal,)),
        ]
        tol = cfg.get("falsification.max_abs_log_ratio")
        n_draws = cfg.get("falsification.n_draws")
        g_exp = (cfg.get("gamma_check.expected_log_hr"), cfg.get("gamma_check.expected_log_hr_ci_lower"),
                 cfg.get("gamma_check.expected_log_hr_ci_upper"))
        g_width = cfg.get("gamma_check.uninformative_log_width")
        run.info(f"falsification tolerance |log ratio| <= {tol}; amyloid HR is reported as '{DESCRIPTIVE}' "
                 f"(published reference log HR {g_exp[0]} [{g_exp[1]}, {g_exp[2]}], {cfg.get('gamma_check.source')})")

        models = build_models(spec, primary, inputs, centre, scale, variants)
        prim: dict[str, dict] = {}
        cvrec: dict[str, dict] = {}
        cvsens: dict[str, dict] = {}
        sens: dict[str, dict] = {}
        flags = {"cv_done": False}

        def stage(name, fn):
            rec = ck.get(name)
            if rec is not None:
                return rec
            if args.max_minutes is not None and (time.time() - t_start) / 60 > args.max_minutes:
                raise OutOfTime(name)
            t0 = time.time()
            run.info(f"stage {name}: started")
            rec = fn()
            rec["seconds"] = time.time() - t0
            rec = ck.put(name, rec)
            run.info(f"stage {name}: finished in {rec['seconds']:.0f} s; checkpoint written")
            return rec

        # ---------------------------------------------- stage functions -----
        def primary_fit(name):
            m = models[name]
            fit, step, dp, dll, bound = fit_audited(m, inputs.data, cfg=cfg, step=step0, start=start,
                                                    eval_mode=primary.eval_mode, n_points=npts, seed=seed, run=run,
                                                    cache=sub, tag=f"primary/{name}/mle")
            for msg in fit.messages:
                if msg.startswith(("WARNING", "NOT")):
                    run.warning(f"{name}: {msg}")
            run.info(f"{name}: MLE logLik {fit.loglik:.2f}, AIC {fit.aic:.2f}, {fit.n_at_best}/{fit.converged_starts} "
                     f"starts at best, Hessian {'ok' if fit.hessian_ok else 'NOT positive definite'}")
            _, summ, diag = fit_bayesian(m, inputs.data, step=step, start_age=start, cfg_bayes=bayes, seed=seed,
                                         init_theta=fit.theta, eval_mode=primary.eval_mode, n_points=npts,
                                         chain_cache=sub.pair(f"primary/{name}/nuts/chain"))
            clean = diag["max_rhat"] <= bayes["max_rhat"] and diag["divergences"] <= bayes["max_divergences"]
            (run.info if clean else run.warning)(
                f"{name} NUTS: max R-hat {diag['max_rhat']:.3f}, min bulk ESS {diag['min_ess_bulk']:.0f}, "
                f"divergences {diag['divergences']}, mean leapfrog {diag['mean_leapfrog_steps']:.1f}, "
                f"{diag['chain_method']} chains" + ("" if clean else " — POSTERIOR NOT CLEAN: reported with this caveat"))
            lo_col, hi_col = [c for c in summ.columns if c.startswith("hdi_")]
            post = {p: {"mean": summ.loc[p.replace("[", "_").replace("]", ""), "mean"],
                        "ci_lower": summ.loc[p.replace("[", "_").replace("]", ""), lo_col],
                        "ci_upper": summ.loc[p.replace("[", "_").replace("]", ""), hi_col],
                        "r_hat": summ.loc[p.replace("[", "_").replace("]", ""), "r_hat"],
                        "ess_bulk": summ.loc[p.replace("[", "_").replace("]", ""), "ess_bulk"]} for p in m.param_names}
            all_draws = diag.pop("draws")
            sel = all_draws[np.linspace(0, len(all_draws) - 1, min(n_draws, len(all_draws))).astype(int)]
            tab, verdict = falsification(m, fit.theta, inputs, threshold=stages.index(primary.threshold),
                                         open_age=primary.open_age, step=step, start=start, tol=tol, draws=sel)
            (run.info if verdict == "PASS" else run.warning)(
                f"falsification {name} (primary): {verdict} (point-estimate rule: {tab['point_rule_verdict'].iloc[0]})")
            return {"fit": _fit_dict(fit), "step": step, "audit_max_dp": dp, "audit_dll": dll, "bound": bound,
                    "posterior": post, "nuts": {**diag, "posterior_clean": clean},
                    "falsification": tab.to_dict("records"), "verdict": verdict}

        def cv_stage(name, scheme):
            # the step the primary fit's audit settled on (prereg §4), not the default
            r = compare_models([models[name]], inputs.data, k=cfg.get("cross_validation.k"), step=prim[name]["step"],
                               start_age=start, cfg_mle=cfg.get("mle"), seed=seed, eval_mode=primary.eval_mode,
                               n_points=npts, warm_starts={name: np.array(prim[name]["fit"]["theta"])},
                               n_random=cfg.get("cross_validation.n_random_starts"), scheme=scheme,
                               split_seed=cfg.get("cross_validation.split_seed"),
                               fold_cache=sub.pair("cv"))[0]
            (run.info if r["unverified_folds"] == 0 else run.warning)(
                f"cross-validation ({scheme}) {name}: held-out logLik {r['cv_heldout_loglik']:.2f}; "
                f"folds without a verified fit: {r['unverified_folds']} of {r['k']}")
            return {"model": name, "cv_heldout_loglik": r["cv_heldout_loglik"], "k": r["k"], "scheme": scheme,
                    "unverified_folds": r["unverified_folds"], "per_fold": r["per_fold"]}

        def sens_stage(sc, name):
            m = build_model(spec, name, sc, inputs, centre, scale)
            if name in sc.refit:
                warm = prim[name]["fit"]["theta"] if name in prim and len(prim[name]["fit"]["theta"]) == len(
                    m.param_names) else None
                fit, step, _, _, _ = fit_audited(m, inputs.data, cfg=cfg, step=step0, start=start,
                                                 eval_mode=sc.eval_mode, n_points=npts, seed=seed, run=run,
                                                 warm=None if warm is None else np.array(warm),
                                                 cache=sub, tag=f"sensitivity/{sc.name}/{name}/mle")
                fd = _fit_dict(fit)
            else:
                fd, step = prim[name]["fit"], prim[name]["step"]
            if name in sc.refit:
                draws = laplace_draws(fd["theta"], fd["cov"], n_draws, seed)
            else:
                draws = laplace_draws(prim[name]["fit"]["theta"], prim[name]["fit"]["cov"], n_draws, seed)
            tab, verdict = falsification(m, np.array(fd["theta"]), inputs, threshold=stages.index(sc.threshold),
                                         open_age=sc.open_age, step=step, start=start, tol=tol, draws=draws)
            run.info(f"sensitivity {sc.name} / {name}: falsification {verdict} "
                     f"(point-estimate rule: {tab['point_rule_verdict'].iloc[0]})")
            return {"scenario": sc.name, "model": name, "refitted": name in sc.refit, "fit": fd, "step": step,
                    "falsification": tab.to_dict("records"), "verdict": verdict}

        # ------------------------------------------------------ tables ------
        def flush_tables():
            params, diags, fals, amy_rows, sens_rows = [], [], [], [], []
            for name, rec in prim.items():
                f = rec["fit"]
                for i, p in enumerate(f["param_names"]):
                    note = (DESCRIPTIVE if p == "log_hr_amyloid" and name == a_name else
                            CONFOUNDED if p == "log_hr_amyloid" else "")
                    params.append({"model": name, "parameter": p, "note": note, "mle": f["theta"][i],
                                   "mle_se": np.nan if f["se"][i] is None else f["se"][i], "posterior_mean": rec["posterior"][p]["mean"],
                                   "ci_lower": rec["posterior"][p]["ci_lower"], "ci_upper": rec["posterior"][p]["ci_upper"],
                                   "r_hat": rec["posterior"][p]["r_hat"], "ess_bulk": rec["posterior"][p]["ess_bulk"],
                                   **({"truth_simulated": truth[i]} if dry and name == "A_gated_age" else {})})
                    if p == "log_hr_amyloid":
                        ci = (rec["posterior"][p]["ci_lower"], rec["posterior"][p]["ci_upper"])
                        amy_rows.append({"label": DESCRIPTIVE if name == a_name else CONFOUNDED, "scenario": "primary",
                                         "model": name, "interval_from": "NUTS 95% HDI",
                                         "log_hr_estimate": rec["posterior"][p]["mean"], "ci_lower": ci[0],
                                         "ci_upper": ci[1], "hr_estimate": np.exp(rec["posterior"][p]["mean"]),
                                         "published_log_hr": g_exp[0],
                                         "descriptive_comparison": expectation_comparison(
                                             ci, *g_exp, on_bound=p in f["on_bound"], max_width=g_width)})
                diags.append({"model": name, "step_years": rec["step"], "error_bound_h2phi": rec["bound"],
                              "audit_max_dp": rec["audit_max_dp"], "audit_dll": rec["audit_dll"],
                              "mle_loglik": f["loglik"], "aic": f["aic"], "hessian_ok": f["hessian_ok"],
                              "on_bound": ";".join(f["on_bound"]), "starts_at_best": f["n_at_best"],
                              "converged_starts": f["converged_starts"], "stage_seconds": rec.get("seconds"),
                              **rec["nuts"]})
                fals += [dict(r, scenario="primary", model=name, verdict=rec["verdict"]) for r in rec["falsification"]]
            for rec in sens.values():
                f = rec["fit"]
                row = {"scenario": rec["scenario"], "model": rec["model"], "refitted": rec["refitted"],
                       "mle_loglik": f["loglik"], "hessian_ok": f["hessian_ok"], "falsification": rec["verdict"]}
                if "log_hr_amyloid" in f["param_names"]:
                    i = f["param_names"].index("log_hr_amyloid")
                    se = f["se"][i] if f["se"][i] is not None else np.nan
                    ci = (f["theta"][i] - 1.96 * se, f["theta"][i] + 1.96 * se)
                    lab = DESCRIPTIVE if rec["model"] == a_name else CONFOUNDED
                    row.update(amyloid_label=lab, log_hr_mle=f["theta"][i], log_hr_ci_lower=ci[0], log_hr_ci_upper=ci[1])
                    amy_rows.append({"label": lab, "scenario": rec["scenario"], "model": rec["model"],
                                     "interval_from": "MLE +/- 1.96 SE", "log_hr_estimate": f["theta"][i],
                                     "ci_lower": ci[0], "ci_upper": ci[1], "hr_estimate": np.exp(f["theta"][i]),
                                     "published_log_hr": g_exp[0],
                                     "descriptive_comparison": expectation_comparison(
                                         ci, *g_exp, on_bound="log_hr_amyloid" in f["on_bound"], max_width=g_width)})
                sens_rows.append(row)
                fals += [dict(r, scenario=rec["scenario"], model=rec["model"], verdict=rec["verdict"])
                         for r in rec["falsification"]]
            out = [(params, "baserate_parameters.csv"), (diags, "baserate_diagnostics.csv"),
                   (amy_rows, "baserate_amyloid_hr_descriptive.csv"), (sens_rows, "baserate_sensitivity.csv"),
                   (fals, "baserate_incidence_falsification.csv")]
            crit = []
            if flags["cv_done"]:
                cv_df = pd.DataFrame([{k: v for k, v in r.items() if k != "per_fold"} for r in cvrec.values()]
                                     ).sort_values("cv_heldout_loglik", ascending=False)
                cv_rows = cv_df.to_dict("records")
                crit.append({"criterion": "2: variant comparison (5-fold CV over brains within age bands)",
                             "verdict": cv_verdict_3way(cv_df, cfg.get("cross_validation.margin"), a_name)})
                if len(cvsens) == len(cvrec):
                    sdf = pd.DataFrame([{k: v for k, v in r.items() if k != "per_fold"} for r in cvsens.values()])
                    cv_rows += sdf.to_dict("records")
                    crit.append({"criterion": "2, sensitivity: whole-age-band CV",
                                 "verdict": cv_verdict_3way(sdf, cfg.get("cross_validation.margin"), a_name)})
                out.append((cv_rows, "baserate_cv.csv"))
            for name, rec in prim.items():
                crit.append({"criterion": f"falsification: {name} (primary)", "verdict": rec["verdict"]})
            crit.append({"criterion": "1: amyloid gating", "verdict": "NOT TESTED (dropped 2026-10-13: no per-brain "
                                                                     "amyloid status in stage-by-age counts)"})
            out.append((crit, "baserate_criteria.csv"))
            for rows, fname in out:
                if rows:
                    run.write_table(pd.DataFrame(rows), tables / fname)

        # ------------------------------------------------------- run --------
        try:
            for name in variants:
                prim[name] = stage(f"primary/{name}", lambda name=name: primary_fit(name))
                flush_tables()
            scheme = cfg.get("cross_validation.scheme")
            for name in variants:
                cvrec[name] = stage(f"cv/{name}", lambda name=name: cv_stage(name, scheme))
            flags["cv_done"] = True
            flush_tables()
            cv_df = pd.DataFrame([{k: v for k, v in r.items() if k != "per_fold"} for r in cvrec.values()])
            run.info(f"PRE-REGISTERED CRITERION 2 (variant comparison): "
                     f"{cv_verdict_3way(cv_df, cfg.get('cross_validation.margin'), a_name)}")
            for sc in scenarios:
                for name in sc.variants:
                    sens[f"{sc.name}/{name}"] = stage(f"sensitivity/{sc.name}/{name}",
                                                      lambda sc=sc, name=name: sens_stage(sc, name))
                    flush_tables()
            s_scheme = cfg.get("cross_validation.scheme_sensitivity")
            for name in variants:
                cvsens[name] = stage(f"sensitivity/cv_{s_scheme}/{name}", lambda name=name: cv_stage(name, s_scheme))
                flush_tables()
        except OutOfTime as e:
            flush_tables()
            run.warning(f"time budget of {args.max_minutes} min spent: stopping before stage {e}; "
                        "rerun the same command to continue from the checkpoints")
            print(f"INCOMPLETE: stopped before stage {e}; rerun to continue", file=sys.stderr)
            return INCOMPLETE_EXIT

        # ------------------------------------------------ final outputs -----
        prev = pd.DataFrame([dict(r, model=n) for n, m in models.items()
                             for r in stage_prevalence_table(m, np.array(prim[n]["fit"]["theta"]), inputs.data,
                                                             step=prim[n]["step"], start_age=start,
                                                             eval_mode=primary.eval_mode, n_points=npts)])
        prev["stage"] = [stages[k] for k in prev["stage_index"]]
        run.write_table(prev, tables / "baserate_prevalence_fit.csv")
        _figures(prim, sens, prev, models, stages, primary, tol, figures, run, a_name, literal, g_exp, dry)
        run.info("ALL STAGES COMPLETE")
    return 0


def _fit_dict(fit) -> dict:
    return {"param_names": fit.param_names, "theta": fit.theta, "se": fit.se, "cov": fit.cov,
            "loglik": fit.loglik, "aic": fit.aic,
            "hessian_ok": fit.hessian_ok, "on_bound": fit.on_bound, "n_at_best": fit.n_at_best,
            "converged_starts": fit.converged_starts, "n_starts": fit.n_starts, "messages": fit.messages}


def _figures(prim, sens, prev, models, stages, primary, tol, figures, run, a_name, literal, g_exp, dry=False):
    apply_theme()
    cols = series_colors(len(models))
    fal = pd.DataFrame([dict(r, model=n, verdict=rec["verdict"]) for n, rec in prim.items() for r in rec["falsification"]])
    fig, ax = plt.subplots(figsize=(6.5, 4))
    env = fal.drop_duplicates("band")
    ax.vlines(env["eval_age"], env["observed_min"], env["observed_max"], color=INK["primary"], lw=3,
              label="SIMULATED target range (truth x0.8 to x1.2)" if dry else "published range (sex x CFAS wave)")
    ax.vlines(env["eval_age"], env["observed_min"] / np.exp(tol), env["observed_max"] * np.exp(tol),
              color=INK["muted"], lw=1, label="pre-set tolerance (x2)")
    for j, ((name, g), c) in enumerate(zip(fal.groupby("model", sort=False), cols)):
        x = g["eval_age"] + 0.6 * (j - 1)                     # small offsets so the intervals do not overlap
        ax.errorbar(x, g["predicted_per_1000py"], yerr=[g["predicted_per_1000py"] - g["pred_lower"],
                                                        g["pred_upper"] - g["predicted_per_1000py"]],
                    fmt="-o", ms=3, lw=1, capsize=2, color=c, label=f"{name}: {g['verdict'].iloc[0]} (95% interval)")
    ax.set_yscale("log")
    ax.set_xlabel("Age (years)")
    ax.set_ylabel(f"Incidence of Braak {primary.threshold}+ / dementia per 1,000 PY")
    ax.legend(fontsize=7)
    save_figure(fig, "baserate_incidence_falsification", run, figures_dir=figures)

    fig, axes = plt.subplots(1, len(models), figsize=(4 * len(models), 3.6), sharey=True)
    for ax, name in zip(np.atleast_1d(axes), models):
        d = prev[prev["model"] == name]
        for k, c in zip(range(len(stages)), plt.cm.viridis(np.linspace(0, 0.95, len(stages)))):
            s = d[d["stage_index"] == k]
            mid = (s["age_lower"] + s["age_upper"]) / 2
            ax.plot(mid, s["observed_prop"], "o", ms=2.5, color=c)
            ax.plot(mid, s["predicted_prop"], "-", lw=1, color=c, label=stages[k])
        ax.set_title(name, fontsize=8)
        ax.set_xlabel("Age (years)")
    np.atleast_1d(axes)[0].set_ylabel("Proportion at stage (dots observed, lines fitted)")
    np.atleast_1d(axes)[-1].legend(fontsize=6, ncol=2)
    save_figure(fig, "baserate_prevalence_fit", run, figures_dir=figures)

    # Amyloid HR: every estimate labelled; the reference band is the published value (not a criterion).
    # The x axis is fixed to log HR in [-XLIM, XLIM] (HR 0.02 to 55) so one unidentified estimate cannot
    # squash the rest; anything beyond is drawn to the edge with an arrow and its row says why.
    XLIM = 4.0
    rows = []
    p = prim[a_name]["posterior"]["log_hr_amyloid"]
    rows.append((f"{a_name} (primary, NUTS 95% HDI)\n{DESCRIPTIVE}", p["mean"], p["ci_lower"], p["ci_upper"],
                 "log_hr_amyloid" in prim[a_name]["fit"]["on_bound"]))
    for rec in sens.values():
        f = rec["fit"]
        if "log_hr_amyloid" in f["param_names"]:
            i = f["param_names"].index("log_hr_amyloid")
            se = f["se"][i] if f["se"][i] is not None else np.nan
            lab = DESCRIPTIVE if rec["model"] == a_name else CONFOUNDED
            rows.append((f"{rec['model']} [{rec['scenario']}] (MLE ± 1.96 SE)\n{lab}", f["theta"][i],
                         f["theta"][i] - 1.96 * se, f["theta"][i] + 1.96 * se, "log_hr_amyloid" in f["on_bound"]))
    fig, ax = plt.subplots(figsize=(6.5, 0.9 + 0.75 * len(rows)))
    ax.axvspan(g_exp[1], g_exp[2], color=INK["muted"], alpha=0.2, lw=0, label="published 95% CI (Ossenkoppele 2022)")
    ax.axvline(0, color=INK["primary"], lw=0.8)
    for y, (lab, est, lo, hi, on_bound) in enumerate(rows):
        c = INK["primary"] if CONFOUNDED not in lab else INK["muted"]
        if on_bound or not (np.isfinite(lo) and np.isfinite(hi)):
            ax.annotate("not identified (estimate on parameter bound)", (0, y), xytext=(0, -9),
                        textcoords="offset points", ha="center", fontsize=6, color=INK["muted"])
            continue
        lo_c, hi_c = max(lo, -XLIM), min(hi, XLIM)
        ax.plot([lo_c, hi_c], [y, y], color=c, lw=2)
        for edge, clipped in ((-XLIM, lo < -XLIM), (XLIM, hi > XLIM)):
            if clipped:
                ax.annotate("", (edge, y), xytext=(edge * 0.9, y), arrowprops=dict(arrowstyle="->", color=c, lw=1.5))
        if -XLIM <= est <= XLIM:
            ax.plot([est], [y], "o", color=c)
    ax.set_xlim(-XLIM * 1.05, XLIM * 1.05)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("log amyloid hazard ratio on the LC → cortex step (arrows: interval continues off scale)")
    ax.set_title(f"Amyloid HR — {DESCRIPTIVE}", fontsize=9)
    ax.legend(fontsize=7, loc="lower right")
    save_figure(fig, "baserate_amyloid_hr_descriptive", run, figures_dir=figures)


if __name__ == "__main__":
    sys.exit(main())
