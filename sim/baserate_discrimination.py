"""SIMULATED: can ~2,300 brains in 10 decade bands discriminate A/B/C at all?
For each truth (A with HR 2 = the dry run's world; B; C), generate EXPECTED
counts (230 x true probabilities) and fit all three variants. With expected
counts, the full-data log-likelihood gap between the truth and a wrong variant
equals n x KL divergence: the AVERAGE evidence a study of this size would get
(no Monte Carlo noise). Gaps < 2 mean the pre-registered margin cannot be met
on average even in-sample.

    python -m sim.baserate_discrimination

Results and interpretation: manuscript/notes/module2_identifiability.md.
"""
import time

import numpy as np

from src.baserate.data import StageCounts
from src.baserate.fit import fit_mle
from src.baserate.model import AGated, BIndependent, ChainSpec, CMultiHit
from src.baserate.run_baserate import simulated_inputs
from src.baserate.solve import stage_probabilities
from src.common.config import load_config
from src.common.paths import CONFIG_DIR

cfg = load_config(CONFIG_DIR / "baserate.yaml")
spec = ChainSpec.from_labels([str(s) for s in cfg.get("data.stage_order")], cfg.get("data.first_cortical_stage"))
step, start, Cn, S = 0.25, 1, 50, 10
inp, thA = simulated_inputs(cfg, spec, step, start, Cn, S)
A = AGated(spec, inp.amyloid["primary"], Cn, S); B = BIndependent(spec, Cn, S); C = CMultiHit(spec, 2, Cn, S)
dA = A.unpack(thA)
def theta_for(m, extra):
    return np.array([extra[n] if n in extra else dA[n] for n in m.param_names])
truths = {
    "A_gated(HR 2)": (A, thA),
    "B_independent": (B, theta_for(B, {"a_nu": dA["log_mu0"] + np.log(1.05), "b_nu": 0.1})),
    "C_multihit(2 hits)": (C, theta_for(C, {"log_lambda": -3.5})),
}
bands = [(g.age_lower, g.age_upper) for g in inp.data]
mle = dict(cfg.get("mle"))
for tname, (tm, th) in truths.items():
    P = stage_probabilities(tm, th, [(lo + hi) / 2 for lo, hi in bands], step=step, start_age=start)
    data = [StageCounts(float(lo), float(hi), "SIMULATED-EXPECTED", 230 * p / p.sum()) for (lo, hi), p in zip(bands, P)]
    res = {}
    for m in (A, B, C):
        t0 = time.time()
        ws = [th] if m is tm else None
        f = fit_mle(m, data, step=step, start_age=start, cfg_mle=mle, seed=3, n_random=3, extra_starts=ws)
        res[m.name] = f.loglik
        print(f"  truth {tname}: fit {m.name} logLik {f.loglik:.3f} ({time.time()-t0:.0f} s)", flush=True)
    best = res[tm.name]
    print(f"truth {tname}: logLik gap (truth minus wrong variant) = "
          + ", ".join(f"{k}: {best - v:.2f}" for k, v in res.items() if k != tm.name), flush=True)
