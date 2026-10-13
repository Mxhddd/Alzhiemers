"""Module 2 entry point: fit variants A, B, C to autopsy stage counts.

    python -m src.baserate.run_baserate                       (real data; gated)
    python -m src.baserate.run_baserate --dry-run-simulated   (simulated; no gates)

Real-data mode refuses to start unless Module 2's validation passed on the
current code and prereg/baserate.md is frozen. The dry run never reads
data/: it simulates counts with the real structure (stage scheme, decade
bands) from a known parameter set, writes to results/dryrun/ with every
output stamped SIMULATED, and exists to test the pipeline end to end and
measure run time before the real run.

Primary analysis (config/baserate.yaml, all fixed before any fit):
  per variant: multi-start MLE with exact gradients; step-error audit against
  a Radau solution (step halved and refitted if it fails); NUTS posterior
  with R-hat / ESS / divergences; k-fold CV by age band; falsification
  against published incidence the model never saw; gamma vs the derivation.
Sensitivity analyses (MLE only, pre-specified): dementia threshold IV;
open incidence band at 87.5; band-average evaluation; 3 hits for C;
Jansen 2015 amyloid curve for A.
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass, replace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.baserate.data import StageCounts, jansen_bands, load_amyloid_bands, load_braak_counts, load_incidence
from src.baserate.fit import fit_bayesian, fit_mle, group_probabilities
from src.baserate.model import AGated, BIndependent, ChainSpec, CMultiHit, StepFunction
from src.baserate.predict import compare_models, predict_incidence, stage_prevalence_table
from src.baserate.solve import error_bound, stage_probabilities_ode
from src.common.config import load_config
from src.common.paths import CONFIG_DIR, FIGURES_DIR, REPO_ROOT, RESULTS_DIR, TABLES_DIR
from src.common.plotting import INK, apply_theme, save_figure, series_colors
from src.common.prereg import require_preregistration
from src.common.runlog import start_run
from src.common.validation import require_validated

CONFIG = CONFIG_DIR / "baserate.yaml"


@dataclass(frozen=True)
class Scenario:
    name: str
    threshold: str
    amyloid: str          # "primary" | "sensitivity"
    eval_mode: str
    c_hits: int
    open_age: float
    refit: tuple[str, ...]  # variants refitted in this scenario (others reuse the primary fit)


@dataclass
class Inputs:
    data: list[StageCounts]
    amyloid: dict[str, StepFunction]
    incidence: pd.DataFrame          # band, eval_age, observed_min, observed_max (open band at primary age)
    incidence_open_lower: float | None


# ---------------------------------------------------------------- inputs ----
BRAAK_TOTALS = {  # quantity in braak2011_abstract_totals.csv -> stages it sums over
    "total brains examined (ages 1-100)": None,
    "AT8-immunonegative after restaging (no abnormal tau)": ["0"],
    "subcortical stages a-c (LC-predominant tau; no cortical tau)": ["a", "b", "c"],
    "stage 1a (transentorhinal neurites)": ["1a"],
    "stage 1b (transentorhinal pyramidal-cell pretangles)": ["1b"],
}


def check_braak_totals(data: list[StageCounts], stages: list[str], run) -> None:
    """The extracted stage table must reproduce the abstract's totals exactly
    (pre-registration §3); any mismatch stops the analysis."""
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
    if problems:
        raise ValueError("Braak 2011 stage table does not reproduce the abstract: " + "; ".join(problems))
    run.info("Braak 2011 stage table reproduces all abstract totals")


def real_inputs(cfg, run) -> Inputs:
    stages = [str(s) for s in cfg.get("data.stage_order")]
    inc = cfg.get("data.upper_inclusive")
    data = load_braak_counts(REPO_ROOT / cfg.get("data.braak_counts"), stages, upper_inclusive=inc, run=run)
    check_braak_totals(data, stages, run)
    primary = StepFunction(load_amyloid_bands(REPO_ROOT / cfg.get("data.amyloid_prevalence"),
                                              upper_inclusive=inc, run=run))
    jb = jansen_bands(REPO_ROOT / cfg.get("data.amyloid_prevalence_sensitivity"), run=run)
    sens = splice_amyloid(primary, jb)
    incidence = load_incidence(REPO_ROOT / cfg.get("data.dementia_incidence"),
                               open_band_age=cfg.get("data.open_band_evaluation_age"), run=run)
    return Inputs(data, {"primary": primary, "sensitivity": sens}, incidence,
                  _open_lower(REPO_ROOT / cfg.get("data.dementia_incidence")))


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


def simulated_inputs(cfg, spec: ChainSpec, step: float, start: float, centre: float, scale: float) -> tuple[Inputs, np.ndarray]:
    """Synthetic world with the real structure. All numbers here are
    SIMULATION SETTINGS, not data, and never leave the dry run."""
    from sim.baserate_sim import simulate_counts

    amy = StepFunction([(1, 41, 0.0), (41, 51, 0.05), (51, 61, 0.15), (61, 71, 0.3), (71, 81, 0.45),
                        (81, 91, 0.6), (91, 101, 0.7)])
    amy_sens = StepFunction([(1, 47.5, 0.0), (47.5, 67.5, 0.1), (67.5, 87.5, 0.3), (87.5, 101, 0.45)])
    truth_model = AGated(spec, amy, centre, scale)
    # log-rate at age 50 and log-slope per decade for each edge, in edge order
    settings = {"alpha[0]": -2.5, "beta[0]": 0.0, "alpha[1]": -2.5, "beta[1]": 0.0, "alpha[2]": -3.0, "beta[2]": 0.1,
                "log_mu0": -4.0, "log_hr_amyloid": np.log(2.0), "alpha[4]": -2.5, "beta[4]": 0.1,
                "alpha[5]": -3.0, "beta[5]": 0.2, "alpha[6]": -3.0, "beta[6]": 0.3, "alpha[7]": -3.5, "beta[7]": 0.5,
                "alpha[8]": -4.0, "beta[8]": 0.6, "alpha[9]": -4.5, "beta[9]": 0.8, "alpha[10]": -4.0, "beta[10]": 0.5}
    theta = np.array([settings[n] for n in truth_model.param_names])
    bands = [(lo, lo + 10) for lo in range(1, 101, 10)]
    data = simulate_counts(truth_model, theta, bands, 230, step=step, start_age=start, seed=cfg.get("seed"))
    ages = np.array([67.5, 72.5, 77.5, 82.5, float(cfg.get("data.open_band_evaluation_age"))])
    thr = list(spec.stages).index(str(cfg.get("data.dementia_threshold_stage")))
    true_inc = 1000 * predict_incidence(truth_model, theta, ages, threshold=thr, step=step, start_age=start)
    incidence = pd.DataFrame({"band": ["65-69", "70-74", "75-79", "80-84", "85-open"], "eval_age": ages,
                              "n_estimates": 4, "observed_min": 0.8 * true_inc, "observed_max": 1.2 * true_inc})
    return Inputs(data, {"primary": amy, "sensitivity": amy_sens}, incidence, 85.0), theta


# ------------------------------------------------------------- analysis ----
def build_models(spec, sc: Scenario, inputs: Inputs, centre, scale, names):
    out = {}
    for v in names:
        if v == "A_gated":
            out[v] = AGated(spec, inputs.amyloid[sc.amyloid], centre, scale)
        elif v == "B_independent":
            out[v] = BIndependent(spec, centre, scale)
        elif v == "C_multihit":
            out[v] = CMultiHit(spec, sc.c_hits, centre, scale)
        else:
            raise ValueError(f"unknown variant {v}")
    return out


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
    ll = lambda P: sum(float(np.sum(g.counts[g.counts > 0] * np.log(np.maximum(p[g.counts > 0], 1e-300))))  # noqa: E731
                       for g, p in zip(data, P))
    return float(np.abs(P_h - P_ode).max()), abs(ll(P_h) - ll(P_ode))


def fit_audited(model, data, *, cfg, step, start, eval_mode, n_points, seed, run, warm=None):
    """MLE, then the pre-registered step audit; halve h and refit on failure.
    `warm` (e.g. the primary fit) is added as one extra start."""
    for attempt in range(4):
        fit = fit_mle(model, data, step=step, start_age=start, cfg_mle=cfg.get("mle"), seed=seed,
                      eval_mode=eval_mode, n_points=n_points, extra_starts=None if warm is None else [warm])
        dp, dll = audit(model, fit.theta, data, step=step, start=start, eval_mode=eval_mode, n_points=n_points)
        end = max(a.max() for a, _ in (g.eval_points(eval_mode, n_points) for g in data))
        bound = step**2 * error_bound(model, fit.theta, start_age=start, end_age=end)
        ok = (dp < cfg.get("solver.max_prob_error") and dll < cfg.get("solver.max_loglik_error")
              and bound <= cfg.get("solver.max_prob_error"))
        run.info(f"{model.name}: step {step:g} y: bound h^2*Phi = {bound:.2e}, Radau audit max|dp| = {dp:.2e}, "
                 f"|dLL| = {dll:.2e} -> {'ok' if ok else 'FAIL: halving the step'}")
        if ok:
            return fit, step, dp, dll
        step /= 2
    raise RuntimeError(f"{model.name}: step audit still failing after 3 halvings")


def falsification(model, theta, inputs: Inputs, *, threshold: int, open_age: float, step, start, tol):
    inc = inputs.incidence.copy()
    if inputs.incidence_open_lower is not None:
        inc.loc[inc["band"].str.endswith("open"), "eval_age"] = open_age
    pred = 1000 * predict_incidence(model, theta, inc["eval_age"].to_numpy(), threshold=threshold, step=step,
                                    start_age=start)
    lo, hi = inc["observed_min"].to_numpy(), inc["observed_max"].to_numpy()
    dist = np.where(pred > hi, np.log(pred / hi), np.where(pred < lo, np.log(lo / pred), 0.0))
    inc["predicted_per_1000py"] = pred
    inc["log_distance_outside_envelope"] = dist
    verdict = "PASS" if np.all(dist <= tol) else "FAIL"
    return inc, verdict


def gating_verdict(interval: tuple[float, float], exp_lo: float) -> str:
    """Pre-registered three-way reading of the amyloid log HR interval (prereg §2.1)."""
    lo, hi = interval
    if lo > 0:
        return "supports gating (interval above 0)"
    if hi < exp_lo:
        return "AGAINST gating (rules out the smallest published-compatible effect)"
    return "inconclusive (interval includes values on both sides; underpowered)"


def cv_verdict_3way(cv_df: pd.DataFrame, margin: float) -> str:
    """Pre-registered three-way reading of the CV comparison (prereg §2.2)."""
    ll = dict(zip(cv_df["model"], cv_df["cv_heldout_loglik"]))
    a = ll["A_gated"]
    others = {k: v for k, v in ll.items() if k != "A_gated"}
    if all(a - v > margin for v in others.values()):
        return "supports (A beats B and C by more than the margin)"
    if any(v - a > margin for v in others.values()):
        return "AGAINST (" + ", ".join(k for k, v in others.items() if v - a > margin) + " beats A by more than the margin)"
    return "inconclusive (A not separated from the alternatives by the margin)"


def gamma_verdict(interval: tuple[float, float], expected: float, exp_lo: float, exp_hi: float) -> str:
    lo, hi = interval
    if lo <= expected <= hi:
        return "consistent"
    if hi < exp_lo or lo > exp_hi:
        return "DISAGREES (intervals do not overlap)"
    return "tension (expected value outside the interval; intervals overlap)"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run-simulated", action="store_true")
    ap.add_argument("--dry-run-nuts", type=int, nargs=3, metavar=("TUNE", "DRAWS", "CHAINS"),
                    help="reduced NUTS settings for timing the dry run")
    ap.add_argument("--dry-run-mle-starts", type=int, help="fewer random MLE starts in the dry run")
    args = ap.parse_args(argv)
    dry = args.dry_run_simulated
    cfg = load_config(CONFIG)
    prereg = None
    if (args.dry_run_nuts or args.dry_run_mle_starts) and not dry:
        raise SystemExit("reduced settings are only allowed with --dry-run-simulated")
    if not dry:
        require_validated("baserate")
        prereg = require_preregistration("baserate")
    tables = RESULTS_DIR / "dryrun" / "tables" if dry else TABLES_DIR
    figures = RESULTS_DIR / "dryrun" / "figures" if dry else FIGURES_DIR

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
            inputs = real_inputs(cfg, run)
        run.info(f"{len(inputs.data)} age-band groups, {sum(g.n for g in inputs.data)} brains, {len(stages)} stages")
        bayes = dict(cfg.get("bayes"))
        if dry and args.dry_run_nuts:
            bayes.update(tune=args.dry_run_nuts[0], draws=args.dry_run_nuts[1], chains=args.dry_run_nuts[2])
            run.warning(f"dry run: NUTS reduced to tune={bayes['tune']}, draws={bayes['draws']}, chains={bayes['chains']}")

        if dry and args.dry_run_mle_starts:
            cfg._data["mle"]["n_starts"] = args.dry_run_mle_starts
            run.warning(f"dry run: MLE random starts reduced to {args.dry_run_mle_starts}")
        variants = list(cfg.get("variants.fit"))
        primary = Scenario("primary", cfg.get("data.dementia_threshold_stage"), "primary", cfg.get("solver.evaluate_at"),
                           cfg.get("variants.c_hits"), cfg.get("data.open_band_evaluation_age"), tuple(variants))
        scenarios = [
            replace(primary, name="threshold_IV", threshold=cfg.get("data.dementia_threshold_stage_sensitivity"), refit=()),
            replace(primary, name="open_band_87.5", open_age=cfg.get("data.open_band_evaluation_age_sensitivity"), refit=()),
            replace(primary, name="band_average", eval_mode=cfg.get("solver.evaluate_at_sensitivity"), refit=tuple(variants)),
            replace(primary, name="c_hits_3", c_hits=cfg.get("variants.c_hits_sensitivity"), refit=("C_multihit",)),
            replace(primary, name="amyloid_jansen2015", amyloid="sensitivity", refit=("A_gated",)),
        ]
        tol = cfg.get("falsification.max_abs_log_ratio")
        g_exp = (cfg.get("gamma_check.expected_log_hr"), cfg.get("gamma_check.expected_log_hr_ci_lower"),
                 cfg.get("gamma_check.expected_log_hr_ci_upper"))
        run.info(f"pre-stated expectation: log HR_A = {g_exp[0]} [{g_exp[1]}, {g_exp[2]}] "
                 f"({cfg.get('gamma_check.source')}); falsification tolerance |log ratio| <= {tol}")

        # ------------------------------------------------ primary ----------
        models = build_models(spec, primary, inputs, centre, scale, variants)
        fits, steps, params, diags, fals, gammas = {}, {}, [], [], [], []
        sens_rows, cv_df = [], None

        def flush_tables():
            """Write every table that has content so far (overwriting earlier
            versions), so a run stopped part-way still leaves its completed
            stages on disk; the run log's footer keeps the final hashes."""
            out = [(pd.DataFrame(params), "baserate_parameters.csv"), (pd.DataFrame(diags), "baserate_diagnostics.csv"),
                   (pd.DataFrame(gammas), "baserate_gamma_check.csv"), (pd.DataFrame(sens_rows), "baserate_sensitivity.csv")]
            if fals:
                out.append((pd.concat(fals, ignore_index=True), "baserate_incidence_falsification.csv"))
            if cv_df is not None:
                out.append((cv_df, "baserate_cv.csv"))
            for df, fname in out:
                if len(df):
                    run.write_table(df, tables / fname)
        for name, m in models.items():
            t0 = time.time()
            fit, step, dp, dll = fit_audited(m, inputs.data, cfg=cfg, step=step0, start=start,
                                             eval_mode=primary.eval_mode, n_points=npts, seed=seed, run=run)
            fits[name], steps[name] = fit, step
            for msg in fit.messages:
                if msg.startswith(("WARNING", "NOT")):
                    run.warning(f"{name}: {msg}")
            run.info(f"{name}: MLE logLik {fit.loglik:.2f}, AIC {fit.aic:.2f}, {fit.n_at_best}/{fit.converged_starts} "
                     f"starts at best, Hessian {'ok' if fit.hessian_ok else 'NOT positive definite'}, "
                     f"{time.time() - t0:.0f} s")
            t1 = time.time()
            idata, summ, diag = fit_bayesian(m, inputs.data, step=step, start_age=start, cfg_bayes=bayes, seed=seed,
                                             init_theta=fit.theta, eval_mode=primary.eval_mode, n_points=npts)
            clean = diag["max_rhat"] <= bayes["max_rhat"] and diag["divergences"] <= bayes["max_divergences"]
            (run.info if clean else run.warning)(
                f"{name} NUTS: max R-hat {diag['max_rhat']:.3f}, min bulk ESS {diag['min_ess_bulk']:.0f}, "
                f"divergences {diag['divergences']}, mean leapfrog {diag['mean_leapfrog_steps']:.1f}, "
                f"{diag['chain_method']} chains, {time.time() - t1:.0f} s"
                + ("" if clean else " — POSTERIOR NOT CLEAN: reported with this caveat"))
            diags.append({"model": name, "step_years": step, "audit_max_dp": dp, "audit_dll": dll,
                          "mle_loglik": fit.loglik, "aic": fit.aic, "hessian_ok": fit.hessian_ok,
                          "on_bound": ";".join(fit.on_bound), "starts_at_best": fit.n_at_best,
                          **diag, "posterior_clean": clean})
            lo_col, hi_col = [c for c in summ.columns if c.startswith("hdi_")]
            for i, pname in enumerate(m.param_names):
                r = summ.loc[pname.replace("[", "_").replace("]", "")]
                params.append({"model": name, "parameter": pname, "mle": fit.theta[i], "mle_se": fit.se[i],
                               "posterior_mean": r["mean"], "ci_lower": r[lo_col], "ci_upper": r[hi_col],
                               "r_hat": r["r_hat"], "ess_bulk": r["ess_bulk"],
                               **({"truth_simulated": truth[i]} if dry and name == "A_gated" else {})})
                if pname == "log_hr_amyloid":
                    v = gamma_verdict((r[lo_col], r[hi_col]), *g_exp)
                    mu0 = np.exp(fit.theta[m.param_names.index("log_mu0")])
                    gv = gating_verdict((r[lo_col], r[hi_col]), g_exp[1])
                    gammas.append({"scenario": "primary", "interval_from": "NUTS 95% HDI", "log_hr_estimate": r["mean"],
                                   "ci_lower": r[lo_col], "ci_upper": r[hi_col], "hr_estimate": np.exp(r["mean"]),
                                   "implied_gamma_per_year_at_mle_mu0": mu0 * (np.exp(r["mean"]) - 1),
                                   "expected_log_hr": g_exp[0], "verdict": v, "gating_verdict": gv})
                    (run.warning if v.startswith("DIS") else run.info)(f"gamma check (primary): {v}")
                    run.info(f"PRE-REGISTERED CRITERION 1 (amyloid gating): {gv}")
            thr = stages.index(primary.threshold)
            tab, verdict = falsification(m, fit.theta, inputs, threshold=thr, open_age=primary.open_age, step=step,
                                         start=start, tol=tol)
            tab.insert(0, "model", name); tab.insert(0, "scenario", "primary"); tab["verdict"] = verdict
            fals.append(tab)
            (run.info if verdict == "PASS" else run.warning)(f"falsification {name} (primary): {verdict}")
            flush_tables()

        cv = compare_models(list(models.values()), inputs.data, k=cfg.get("cross_validation.k"), step=step0,
                            start_age=start, cfg_mle=cfg.get("mle"), seed=seed, eval_mode=primary.eval_mode,
                            n_points=npts, warm_starts={n: f.theta for n, f in fits.items()},
                            n_random=cfg.get("cross_validation.n_random_starts"))
        cv_df = pd.DataFrame([{k: v for k, v in r.items() if k != "per_fold"} for r in cv]).sort_values(
            "cv_heldout_loglik", ascending=False)
        best, second = cv_df.iloc[0], cv_df.iloc[1]
        margin = cfg.get("cross_validation.margin")
        cv_verdict = (f"{best['model']} best by {best['cv_heldout_loglik'] - second['cv_heldout_loglik']:.2f} "
                      f"(margin {margin}): " + ("clear" if best["cv_heldout_loglik"] - second["cv_heldout_loglik"] > margin
                                                else "NOT distinguishable from the runner-up"))
        run.info(f"cross-validation: {cv_verdict}")
        run.info(f"PRE-REGISTERED CRITERION 2 (variant comparison): {cv_verdict_3way(cv_df, margin)}")
        flush_tables()

        # --------------------------------------------- sensitivity ---------
        for sc in scenarios:
            sm = build_models(spec, sc, inputs, centre, scale, variants)
            for name, m in sm.items():
                if name in sc.refit:
                    fit, step, _, _ = fit_audited(m, inputs.data, cfg=cfg, step=step0, start=start,
                                                  eval_mode=sc.eval_mode, n_points=npts, seed=seed, run=run,
                                                  warm=fits[name].theta if len(m.param_names) == len(fits[name].theta) else None)
                else:
                    fit, step = fits[name], steps[name]
                thr = stages.index(sc.threshold)
                tab, verdict = falsification(m, fit.theta, inputs, threshold=thr, open_age=sc.open_age, step=step,
                                             start=start, tol=tol)
                tab.insert(0, "model", name); tab.insert(0, "scenario", sc.name); tab["verdict"] = verdict
                fals.append(tab)
                row = {"scenario": sc.name, "model": name, "refitted": name in sc.refit, "mle_loglik": fit.loglik,
                       "falsification": verdict}
                if "log_hr_amyloid" in m.param_names:
                    i = m.param_names.index("log_hr_amyloid")
                    ci = (fit.theta[i] - 1.96 * fit.se[i], fit.theta[i] + 1.96 * fit.se[i])
                    v = gamma_verdict(ci, *g_exp)
                    row.update(log_hr_mle=fit.theta[i], log_hr_ci_lower=ci[0], log_hr_ci_upper=ci[1], gamma_check=v,
                               gating_verdict=gating_verdict(ci, g_exp[1]))
                    gammas.append({"scenario": sc.name, "interval_from": "MLE +/- 1.96 SE", "log_hr_estimate": fit.theta[i],
                                   "ci_lower": ci[0], "ci_upper": ci[1], "hr_estimate": np.exp(fit.theta[i]),
                                   "implied_gamma_per_year_at_mle_mu0":
                                       np.exp(fit.theta[m.param_names.index("log_mu0")]) * (np.exp(fit.theta[i]) - 1),
                                   "expected_log_hr": g_exp[0], "verdict": v})
                sens_rows.append(row)
                run.info(f"sensitivity {sc.name} / {name}: falsification {verdict}"
                         + (f", gamma {row['gamma_check']}" if "gamma_check" in row else ""))
                flush_tables()

        prev = pd.DataFrame([dict(r, model=n) for n, m in models.items()
                             for r in stage_prevalence_table(m, fits[n].theta, inputs.data, step=steps[n], start_age=start,
                                                             eval_mode=primary.eval_mode, n_points=npts)])
        prev["stage"] = [stages[k] for k in prev["stage_index"]]
        fal_df = pd.concat(fals, ignore_index=True)
        flush_tables()
        run.write_table(prev, tables / "baserate_prevalence_fit.csv")

        apply_theme()
        cols = series_colors(len(models))
        fig, ax = plt.subplots(figsize=(6.5, 4))
        p0 = fal_df[fal_df["scenario"] == "primary"]
        env = p0.drop_duplicates("band")
        ax.vlines(env["eval_age"], env["observed_min"], env["observed_max"], color=INK["primary"], lw=3,
                  label="published range (sex x CFAS wave)")
        ax.vlines(env["eval_age"], env["observed_min"] / np.exp(tol), env["observed_max"] * np.exp(tol),
                  color=INK["muted"], lw=1, label="pre-set tolerance (x2)")
        for (name, g), c in zip(p0.groupby("model", sort=False), cols):
            ax.plot(g["eval_age"], g["predicted_per_1000py"], "-o", ms=3, color=c,
                    label=f"{name}: {g['verdict'].iloc[0]}")
        ax.set_yscale("log")
        ax.set_xlabel("Age (years)")
        ax.set_ylabel(f"Incidence of Braak {primary.threshold}+ / dementia per 1,000 PY")
        ax.legend(fontsize=7)
        save_figure(fig, "baserate_incidence_falsification", run, figures_dir=figures)

        fig, axes = plt.subplots(1, len(models), figsize=(4 * len(models), 3.6), sharey=True)
        for ax, (name, m) in zip(np.atleast_1d(axes), models.items()):
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
