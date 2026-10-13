"""SIMULATED power study for Module 2's pre-registered amyloid-gating criterion.

In the dry run's synthetic world (12 Braak-2011 stages, 10 decade bands, 230
brains per band; settings in src/baserate/run_baserate.py::simulated_inputs),
how often does the three-way reading of log HR_A come out supports / against /
inconclusive when the true HR is (a) 1 (no amyloid effect) and (b) 2.43 (the
published value used for the expectation)? Uses MLE Wald intervals (the
sensitivity-analysis estimator), warm-started at the truth plus one random
start. Indicative only: power depends on the true stage distribution, which
this synthetic world only imitates.

    python -m sim.baserate_power --reps 30
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

from sim.baserate_sim import simulate_counts
from src.baserate.fit import fit_mle
from src.baserate.model import ChainSpec
from src.baserate.run_baserate import Scenario, build_models, gating_verdict, simulated_inputs
from src.common.config import load_config
from src.common.paths import CONFIG_DIR, RESULTS_DIR
from src.common.runlog import start_run


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, required=True)
    args = ap.parse_args(argv)
    cfg = load_config(CONFIG_DIR / "baserate.yaml")
    spec = ChainSpec.from_labels([str(s) for s in cfg.get("data.stage_order")], cfg.get("data.first_cortical_stage"))
    step, start = cfg.get("solver.step_years"), cfg.get("solver.start_age")
    inp, truth = simulated_inputs(cfg, spec, step, start, cfg.get("solver.age_centre"), cfg.get("solver.age_scale"))
    model = build_models(spec, Scenario("p", "V", "primary", "midpoint", 2, 90, ()), inp,
                         cfg.get("solver.age_centre"), cfg.get("solver.age_scale"), ["A_gated"])["A_gated"]
    i = model.param_names.index("log_hr_amyloid")
    exp_lo = cfg.get("gamma_check.expected_log_hr_ci_lower")
    bands = [(g.age_lower, g.age_upper) for g in inp.data]
    mle_cfg = dict(cfg.get("mle"))
    rows = []
    with start_run(module="baserate_power_sim", config_path=CONFIG_DIR / "baserate.yaml", seed=cfg.get("seed"),
                   packages=["numpy", "scipy", "jax"]) as run:
        run.mark_simulated()
        for true_log_hr in (0.0, float(np.log(2.43))):
            th = truth.copy(); th[i] = true_log_hr
            for r in range(args.reps):
                data = simulate_counts(model, th, bands, 230, step=step, start_age=start, seed=10_000 + r)
                fit = fit_mle(model, data, step=step, start_age=start, cfg_mle=mle_cfg, seed=r,
                              extra_starts=[th], n_random=1)
                ci = (fit.theta[i] - 1.96 * fit.se[i], fit.theta[i] + 1.96 * fit.se[i])
                rows.append({"true_hr": float(np.exp(true_log_hr)), "rep": r, "log_hr_mle": fit.theta[i],
                             "ci_lower": ci[0], "ci_upper": ci[1], "hessian_ok": fit.hessian_ok,
                             "verdict": gating_verdict(ci, exp_lo).split(" (")[0]})
            sub = pd.DataFrame([x for x in rows if np.isclose(x["true_hr"], np.exp(true_log_hr))])
            run.info(f"true HR {np.exp(true_log_hr):.2f}: " + ", ".join(
                f"{k} {v / len(sub):.0%}" for k, v in sub["verdict"].value_counts().items()))
        run.write_table(pd.DataFrame(rows), RESULTS_DIR / "dryrun" / "tables" / "baserate_power_sim.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
