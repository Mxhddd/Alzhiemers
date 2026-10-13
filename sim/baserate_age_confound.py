"""SIMULATED: is variant A's amyloid HR confounded with age?
Truth: sequential chain, NO amyloid effect, c->1a rate exp(alpha3 + beta3*u).
Fit A_gated to the EXPECTED counts (230 x true probabilities per decade band),
which gives the large-sample (pseudo-true) log HR_A without Monte Carlo noise.

    python -m sim.baserate_age_confound 0.0 0.1 0.3

Results and interpretation: manuscript/notes/module2_identifiability.md.
"""
import sys
import time

import numpy as np

from src.baserate.data import StageCounts
from src.baserate.fit import fit_mle
from src.baserate.model import ChainSpec, Sequential
from src.baserate.run_baserate import Scenario, build_models, simulated_inputs
from src.baserate.solve import stage_probabilities
from src.common.config import load_config
from src.common.paths import CONFIG_DIR

cfg = load_config(CONFIG_DIR / "baserate.yaml")
spec = ChainSpec.from_labels([str(s) for s in cfg.get("data.stage_order")], cfg.get("data.first_cortical_stage"))
step, start, C, S = 0.25, 1, 50, 10
inp, truthA = simulated_inputs(cfg, spec, step, start, C, S)
A = build_models(spec, Scenario("p", "V", "primary", "midpoint", 2, 90, ()), inp, C, S, ["A_gated"])["A_gated"]
seq = Sequential(spec, C, S)
dA = A.unpack(truthA)
bands = [(g.age_lower, g.age_upper) for g in inp.data]
mle = dict(cfg.get("mle"))
for beta3 in [float(x) for x in sys.argv[1:]]:
    th = []
    for nm in seq.param_names:
        if nm == "alpha[3]":
            th.append(dA["log_mu0"] + np.log(1.05))   # same c->1a rate at age 50 as the A truth
        elif nm == "beta[3]":
            th.append(beta3)
        else:
            th.append(dA[nm])
    th = np.array(th)
    P = stage_probabilities(seq, th, [(lo + hi) / 2 for lo, hi in bands], step=step, start_age=start)
    data = [StageCounts(float(lo), float(hi), "SIMULATED-EXPECTED", 230 * p / p.sum()) for (lo, hi), p in zip(bands, P)]
    t0 = time.time()
    fit = fit_mle(A, data, step=step, start_age=start, cfg_mle=mle, seed=1, n_random=2,
                  extra_starts=[np.where(np.array(A.param_names) == "log_hr_amyloid", 0.0, truthA)])
    i = A.param_names.index("log_hr_amyloid")
    print(f"truth: no amyloid effect, c->1a log-slope {beta3}/decade -> A_gated pseudo-true log HR_A = "
          f"{fit.theta[i]:.3f} (HR {np.exp(fit.theta[i]):.2f}), SE at n=2300 {fit.se[i]:.3f}, "
          f"95% CI [{fit.theta[i]-1.96*fit.se[i]:.2f}, {fit.theta[i]+1.96*fit.se[i]:.2f}], {time.time()-t0:.0f} s", flush=True)
