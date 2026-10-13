"""Independent machine check of manuscript/derivations/module2_gamma.md.

  1. symbolic: the two-type mixture gives mu0 + gamma A with
     gamma = mu0 (HR - 1), and the age-adjusted form m(t)(1 + (HR - 1)A(t));
  2. the numbers: ln 2.43, ln 1.38, ln 4.26, the implied SE, gamma/mu0;
  3. claim A3 (attenuation) numerically: an INDIVIDUAL-level model in which
     each person's own amyloid status multiplies their LC -> cortex rate by
     HR, fitted with the population-level form (which uses only the
     population prevalence A(t)) on the expected counts. A3 predicts the
     fitted HR is pulled towards 1.

    uv run python manuscript/derivations/checks/module2_gamma.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
results = []


def check(name, ok, detail=""):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}{(': ' + detail) if detail else ''}")


# 1. symbolic -------------------------------------------------------------------
mu0, g, A, HR, m = sp.symbols("mu0 gamma A HR m", positive=True)
mix = (1 - A) * mu0 + A * (mu0 + g)
check("mixture (1-A) mu0 + A (mu0 + gamma) = mu0 + gamma A", sp.simplify(mix - (mu0 + g * A)) == 0)
sol = sp.solve(sp.Eq((mu0 + g) / mu0, HR), g)
check("HR = (mu0 + gamma)/mu0  <=>  gamma = mu0 (HR - 1)", len(sol) == 1 and sp.simplify(sol[0] - mu0 * (HR - 1)) == 0,
      str(sol))
check("age-adjusted form: (1-A) m + A HR m = m (1 + (HR - 1) A)",
      sp.simplify((1 - A) * m + A * HR * m - m * (1 + (HR - 1) * A)) == 0)

# 2. numbers ----------------------------------------------------------------------
lhr, llo, lhi = math.log(2.43), math.log(1.38), math.log(4.26)
se = (lhi - llo) / (2 * 1.959963984540054)
check("log HR = 0.888 (0.322-1.449), SE 0.288", (round(lhr, 3), round(llo, 3), round(lhi, 3), round(se, 3)) ==
      (0.888, 0.322, 1.449, 0.288), f"{lhr:.4f} ({llo:.4f}-{lhi:.4f}), SE {se:.4f}")
check("gamma/mu0 = HR - 1 = 1.43 (0.38-3.26)", (round(2.43 - 1, 2), round(1.38 - 1, 2), round(4.26 - 1, 2)) ==
      (1.43, 0.38, 3.26))
# published CI symmetric on the log scale? (a check of the extracted values, not of the derivation)
mid = (llo + lhi) / 2
check("published CI is symmetric on the log scale around ln 2.43 (to rounding)", abs(mid - lhr) < 0.01,
      f"midpoint {mid:.4f} vs {lhr:.4f}")

# 3. A3: attenuation ----------------------------------------------------------------
# Individual model on a 3-stage chain (none -> LC -> cortex) x amyloid (neg -> pos):
# hidden state (stage, amyloid). Amyloid onset hazard a(t) chosen so that population
# prevalence rises with age; LC -> cortex rate mu0 * HR^(amyloid). Expected stage
# counts by decade band are then fitted by the population form
# mu0 (1 + (HR_pop - 1) A(t)), A(t) = the SAME model's population amyloid prevalence.
from src.baserate.data import StageCounts  # noqa: E402
from src.baserate.fit import fit_mle  # noqa: E402
from src.baserate.model import AGated, ChainSpec, StepFunction  # noqa: E402


def individual_probs(ages, hr, mu0=0.02, r_lc=0.08, onset=lambda t: 0.002 * math.exp(0.07 * (t - 40)), step=0.05):
    # states: 0 (none,-) 1 (none,+) 2 (LC,-) 3 (LC,+) 4 (cortex,-) 5 (cortex,+)
    p = np.zeros(6); p[0] = 1.0
    out, t, k = {}, 0.0, 0
    targets = sorted(ages)
    amy = {}
    while k < len(targets):
        nxt = min(t + step, targets[k])
        tm = 0.5 * (t + nxt)
        a = onset(tm)
        Q = np.zeros((6, 6))
        for s in (0, 2, 4):
            Q[s, s + 1] += a                      # amyloid onset (independent of stage)
        Q[0, 2] += r_lc; Q[1, 3] += r_lc          # none -> LC
        Q[2, 4] += mu0; Q[3, 5] += mu0 * hr       # LC -> cortex, x HR if amyloid-positive
        Q[np.diag_indices(6)] = -Q.sum(1)
        p = p @ expm(Q * (nxt - t))
        t = nxt
        if abs(t - targets[k]) < 1e-12:
            out[targets[k]] = np.array([p[0] + p[1], p[2] + p[3], p[4] + p[5]])
            amy[targets[k]] = p[1] + p[3] + p[5]
            k += 1
    return out, amy


bands = [(lo, lo + 10) for lo in range(0, 100, 10)]
mids = [lo + 5.0 for lo, _ in bands]
cfg_mle = {"n_starts": 4, "method": "l-bfgs-b", "maxiter": 3000, "newton_steps": 20, "max_newton_decrement": 1e-3,
           "init": {"log_rate": [-6.0, -2.0], "slope": [-0.5, 1.0], "log_hr": [-0.5, 0.5]},
           "bounds": {"log_rate": [-20.0, 2.0], "slope": [-2.0, 2.0], "log_hr": [-5.0, 5.0]}}
for hr_true in (1.0, 2.43, 5.0):
    probs, amy = individual_probs(mids, hr_true)
    # population prevalence A(t) at a fine grid, as the step function the population model uses
    fine = list(np.arange(0.5, 100, 1.0))
    _, amy_f = individual_probs(fine, hr_true)
    A_step = StepFunction([(t - 0.5, t + 0.5, float(amy_f[t])) for t in fine])
    data = [StageCounts(float(lo), float(hi), "EXPECTED", 100000 * probs[lo + 5.0]) for lo, hi in bands]
    model = AGated(ChainSpec(("none", "LC", "cortex"), 2), A_step, 50, 10)
    start = np.array([math.log(0.08), 0.0, math.log(0.02), math.log(hr_true)])
    fit = fit_mle(model, data, step=0.5, start_age=0, cfg_mle=cfg_mle, seed=1, extra_starts=[start], n_random=2)
    hr_fit = math.exp(fit.theta[model.param_names.index("log_hr_amyloid")])
    if hr_true == 1.0:
        check("A3 control: with HR = 1 the population form recovers HR = 1", abs(hr_fit - 1) < 0.02, f"fitted {hr_fit:.3f}")
    else:
        check(f"A3: individual HR {hr_true} is attenuated towards 1 by the population form", 1 < hr_fit < hr_true,
              f"fitted population-form HR {hr_fit:.3f}")

print(f"\n{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)
