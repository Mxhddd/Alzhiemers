"""SIMULATED: refit chosen power-study replicates with the analysis's own multi-start settings.

Diagnostic for the FIRST cross-validation design (whole age bands held out),
whose results are archived in results/dryrun/power/whole_band_design/. It found
the stage-only-in-the-held-out-band failure that led to the 2026-10-14 revision.

    python -m sim.baserate_power_refit "LC-first, no amyloid, c->1a slope 0" 19

The power study (sim/baserate_power.py) fits each replicate with one warm-started
L-BFGS-B run per fit, for speed. The pre-registered analysis instead uses
mle.n_starts random starts for the full data and the full-data MLE plus
cross_validation.n_random_starts random starts per fold, each verified by the
Newton decrement (src/baserate/fit.py::fit_mle). This script rebuilds one
replicate's data exactly (same truth, same seed) and refits it both ways, so a
surprising power-study verdict can be traced to the data or to the shortcut.
Writes results/dryrun/power/refit_<truth index>_<rep>.json. Reports both
verdicts; nothing in the power-study table is changed.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from sim.baserate_power import BANDS, N_PER_BAND, OUT, POWER_STEP
from sim.baserate_sim import simulate_counts
from src.baserate.data import StageCounts
from src.baserate.fit import fit_mle
from src.baserate.jaxlik import make_loglik_family
from src.baserate.model import AGatedAge, BIndependent, ChainSpec, CMultiHit
from src.baserate.predict import age_band_folds
from src.baserate.run_baserate import SIM_TRUTH, cv_verdict_3way, simulated_inputs
from src.common.config import load_config
from src.common.paths import CONFIG_DIR


def main(truth_name: str, rep: int) -> int:
    cfg = load_config(CONFIG_DIR / "baserate.yaml")
    spec = ChainSpec.from_labels([str(s) for s in cfg.get("data.stage_order")], cfg.get("data.first_cortical_stage"))
    step, start = POWER_STEP, cfg.get("solver.start_age")
    centre, scale = cfg.get("solver.age_centre"), cfg.get("solver.age_scale")
    k, margin = cfg.get("cross_validation.k"), cfg.get("cross_validation.margin")
    mle = dict(cfg.get("mle"))
    inp, _ = simulated_inputs(cfg, spec, step, start, centre, scale)
    amy = inp.amyloid["primary"]
    A, B, C = AGatedAge(spec, amy, centre, scale), BIndependent(spec, centre, scale), CMultiHit(spec, 2, centre, scale)
    models = {"A_gated_age": A, "B_independent": B, "C_multihit": C}

    def theta_for(m, extra):
        d = dict(SIM_TRUTH, **extra)
        return np.array([d[n] for n in m.param_names])

    # identical to sim/baserate_power.py (order matters: it sets the seed)
    truths = {
        "LC-first, no amyloid, c->1a slope 0": (A, theta_for(A, {"beta_mu": 0.0, "log_hr_amyloid": 0.0})),
        "LC-first, no amyloid, c->1a slope 0.1": (A, theta_for(A, {"beta_mu": 0.1, "log_hr_amyloid": 0.0})),
        "LC-first, no amyloid, c->1a slope 0.3": (A, theta_for(A, {"beta_mu": 0.3, "log_hr_amyloid": 0.0})),
        "B independent": (B, theta_for(B, {"a_nu": SIM_TRUTH["log_mu0"] + np.log(1.05), "b_nu": 0.1})),
        "C multi-hit (2 hits)": (C, theta_for(C, {"log_lambda": -3.5})),
    }
    ti = list(truths).index(truth_name)
    tm, th = truths[truth_name]
    data = simulate_counts(tm, th, BANDS, N_PER_BAND, step=step, start_age=start, seed=100_000 + 1_000 * ti + rep)
    counts = np.array([g.counts for g in data], dtype=float)
    template = [StageCounts(float(lo), float(hi), "SIMULATED", np.zeros(len(spec.stages))) for lo, hi in BANDS]
    folds = np.array(age_band_folds(template, k))
    fam = {n: make_loglik_family(m, template, step=step, start_age=start) for n, m in models.items()}
    ones = np.ones(len(BANDS))

    def bound(name, mask):
        ll, vg, h = fam[name]
        return (lambda t: ll(t, counts, mask), lambda t: vg(t, counts, mask), lambda t: h(t, counts, mask))

    out = {"truth": truth_name, "rep": rep, "step": step, "models": {}}
    for name, m in models.items():
        full = fit_mle(m, [], step=step, start_age=start, cfg_mle=mle, seed=cfg.get("seed"),
                       compiled=bound(name, ones))                      # mle.n_starts random starts
        held = 0.0
        for f in range(k):
            train = (folds != f).astype(float)
            ff = fit_mle(m, [], step=step, start_age=start, cfg_mle=mle, seed=cfg.get("seed") + f,
                         extra_starts=[full.theta], n_random=cfg.get("cross_validation.n_random_starts"),
                         compiled=bound(name, train))
            held += float(fam[name][0](ff.theta, counts, 1.0 - train))
        out["models"][name] = {"full_loglik": full.loglik, "cv_heldout_loglik": held,
                               "converged_starts": full.converged_starts, "n_at_best": full.n_at_best}
        print(f"{name}: full {full.loglik:.2f} (converged starts {full.converged_starts}), held-out {held:.2f}", flush=True)
    cv_df = pd.DataFrame({"model": list(models), "cv_heldout_loglik": [out["models"][n]["cv_heldout_loglik"] for n in models]})
    out["verdict_multistart"] = cv_verdict_3way(cv_df, margin, "A_gated_age")
    reps = pd.concat([pd.read_csv(p) for p in (OUT / "whole_band_design").glob("baserate_power_reps_shard*.csv")])
    old = reps[(reps["truth"] == truth_name) & (reps["rep"] == rep)].iloc[0]
    out["verdict_power_study"] = old["verdict"]
    out["power_study_cv"] = {n: float(old[f"cv_{n}"]) for n in models}
    (OUT / "whole_band_design" / f"refit_{ti}_{rep}.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"power study: {old['verdict']} | multi-start refit: {out['verdict_multistart']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], int(sys.argv[2])))
