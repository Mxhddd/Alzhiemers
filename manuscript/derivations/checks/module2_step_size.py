"""Independent machine check of manuscript/derivations/module2_step_size.md.

Re-derives each step by a different route from the one in the document:
  1. the two Magnus integrals, symbolically (sympy);
  2. the local error (★) from the EXACT Dyson (time-ordered) series of the
     row-form propagator, in exact rational arithmetic with random rational
     generator matrices: U - exp(hQ0) - (★) must have no terms below h^4;
  3. the commutator closed form for a sequential chain, symbolically;
  4. the a-priori envelope arithmetic;
  5. numerically, that the actual global error of the stepped solution is
     below h^2 Phi(theta) at the dry-run parameter set, for two step sizes.

    uv run python manuscript/derivations/checks/module2_step_size.py
Prints PASS/FAIL per check and exits non-zero on any failure.
"""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
results = []


def check(name, ok, detail=""):
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}{(': ' + detail) if detail else ''}")


# 1. Magnus integrals ---------------------------------------------------------
h, t1, t2, tau = sp.symbols("h t1 t2 tau", positive=True)
inner = sp.integrate(t2 - t1, (t2, -h / 2, t1))
outer = sp.simplify(sp.integrate(inner, (t1, -h / 2, h / 2)))
check("double integral of (tau2 - tau1) over tau2 < tau1 equals -h^3/6", sp.simplify(outer + h**3 / 6) == 0, str(outer))
om1 = sp.integrate(tau**2 / 2, (tau, -h / 2, h / 2))
check("integral of tau^2/2 over [-h/2, h/2] equals h^3/24", sp.simplify(om1 - h**3 / 24) == 0, str(om1))
check("integral of tau over [-h/2, h/2] vanishes", sp.integrate(tau, (tau, -h / 2, h / 2)) == 0)


# 2. Dyson series vs (★) --------------------------------------------------------
def rand_generator(seed, n=3):
    import random
    rnd = random.Random(seed)
    M = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                M[i, j] = sp.Rational(rnd.randint(-5, 9), rnd.randint(2, 7))
        M[i, i] = -sum(M[i, j] for j in range(n) if j != i)
    return M


s = sp.symbols("s")
ok_all = True
for seed in (1, 2, 3):
    Q0, Q1, Q2 = rand_generator(seed), rand_generator(seed + 10), rand_generator(seed + 20)
    Q = lambda x: Q0 + x * Q1 + x**2 / 2 * Q2  # noqa: E731  Q(m + x), x in [-h/2, h/2]
    a = -h / 2
    # Row form U' = U Q, U(a) = I: U(t) = I + int_a^t U(s) Q(s) ds. Iterate to third order.
    I3 = sp.eye(3)
    U1 = (Q(s)).applyfunc(lambda e: sp.integrate(e, (s, a, t1)))                     # int_a^t1 Q
    U2 = (U1.subs(t1, s) * Q(s)).applyfunc(lambda e: sp.integrate(e, (s, a, t1)))    # int_a^t1 U1(s) Q(s)
    U3 = (U2.subs(t1, s) * Q(s)).applyfunc(lambda e: sp.integrate(e, (s, a, t1)))
    U = (I3 + U1 + U2 + U3).subs(t1, h / 2)
    expm_series = I3 + h * Q0 + (h * Q0) ** 2 / 2 + (h * Q0) ** 3 / 6
    star = h**3 / 24 * Q2 + h**3 / 12 * (Q0 * Q1 - Q1 * Q0)
    R = (U - expm_series - star).applyfunc(lambda e: sp.expand(e))
    low = [sp.Poly(e, h).coeff_monomial(h**p) for e in R for p in range(0, 4)]
    ok = all(c == 0 for c in low)
    # negative control: a wrong commutator coefficient (1/6 instead of 1/12) must leave an h^3 term
    wrong = (U - expm_series - h**3 / 24 * Q2 - h**3 / 6 * (Q0 * Q1 - Q1 * Q0)).applyfunc(sp.expand)
    neg = any(sp.Poly(e, h).coeff_monomial(h**3) != 0 for e in wrong)
    ok_all &= ok and neg
check("Dyson series: U - exp(hQ0) - [(h^3/24)Q'' + (h^3/12)[Q,Q']] has no h^0..h^3 terms (3 random rational "
      "generators), and a wrong coefficient 1/6 leaves an h^3 term (negative control)", ok_all)

# 3. Commutator closed form ---------------------------------------------------
K = 4
tt = sp.symbols("t")
a_s = sp.symbols("a0:3"); b_s = sp.symbols("b0:3")
r = [sp.exp(a_s[k] + b_s[k] * tt) for k in range(K - 1)]
Qs = sp.zeros(K, K)
for k in range(K - 1):
    Qs[k, k + 1] = r[k]
    Qs[k, k] = -r[k]
Qd = Qs.diff(tt)
C = (Qs * Qd - Qd * Qs).applyfunc(sp.simplify)
want = sp.zeros(K, K)
for k in range(K - 2):
    want[k, k + 1] = r[k] * r[k + 1] * (b_s[k] - b_s[k + 1])
    want[k, k + 2] = r[k] * r[k + 1] * (b_s[k + 1] - b_s[k])
check("commutator [Q, Q'] of a sequential chain has exactly the entries in §3", (C - want).applyfunc(sp.simplify) == sp.zeros(K, K))

# 4. Envelope arithmetic --------------------------------------------------------
eps = Fraction(2, 10) * Fraction(15, 100) ** 2 / 12 + Fraction(2, 10) ** 2 * Fraction(30, 100) / 6
phi = 100 * eps
h_max = (Fraction(1, 1000) / phi) ** 0.5
check("envelope: eps = 2.375e-3 per year, Phi = 0.2375, h <= 0.065 y", eps == Fraction(2375, 10**6) and
      phi == Fraction(2375, 10**4) and abs(float(h_max) - 0.0649) < 5e-4, f"eps={float(eps)}, Phi={float(phi)}, h={float(h_max):.4f}")
check("realistic rule: Phi = 0.0099 gives h <= 0.32 y", abs((1e-3 / 0.0099) ** 0.5 - 0.318) < 1e-3,
      f"h = {(1e-3 / 0.0099) ** 0.5:.3f}")

# 5. Global error vs h^2 Phi at the dry-run parameters ---------------------------
import numpy as np  # noqa: E402

from src.baserate.model import AGatedAge, ChainSpec, StepFunction  # noqa: E402
from src.baserate.run_baserate import SIM_AMYLOID, SIM_TRUTH  # noqa: E402
from src.baserate.solve import error_bound, stage_probabilities, stage_probabilities_ode  # noqa: E402

spec = ChainSpec.from_labels(["0", "a", "b", "c", "1a", "1b", "I", "II", "III", "IV", "V", "VI"], "1a")
m = AGatedAge(spec, StepFunction(SIM_AMYLOID), 50, 10)
th = np.array([SIM_TRUTH[n] for n in m.param_names])
ages = np.arange(6, 100, 5.0)
ref = stage_probabilities_ode(m, th, ages, start_age=1)
phi_th = error_bound(m, th, start_age=1, end_age=ages.max())
for hh in (0.5, 0.25):
    err = float(np.abs(stage_probabilities(m, th, ages, step=hh, start_age=1) - ref).max())
    check(f"global error at h={hh} is below h^2 Phi(theta)", err <= hh**2 * phi_th,
          f"max|p_h - p| = {err:.2e} vs h^2 Phi = {hh**2 * phi_th:.2e} (Phi = {phi_th:.4f})")
row_sums = []
for t in (30.0, 60.0, 85.0):
    from src.baserate.model import build_generator  # noqa: E402
    d = 1e-3
    Qp, Q0n, Qm = (build_generator(m, t + d, th), build_generator(m, t, th), build_generator(m, t - d, th))
    Q2n = (Qp - 2 * Q0n + Qm) / d**2
    rates = np.array([Q0n[k, k + 1] for k in range(len(spec.stages) - 1)])
    bpy = np.array([0.0] * len(rates))
    for k, n in enumerate([e for e in m.edges]):
        bpy[k] = th[n.p1] / 10 if n.kind != "const" else 0.0
    row_sums.append(np.allclose(np.abs(Q2n).sum(axis=1)[:-1], 2 * bpy**2 * rates, rtol=1e-4, atol=1e-12))
check("numerically, row sums of |Q''| equal 2 b_k^2 r_k (dry-run parameters, 3 ages, inside amyloid bands)",
      all(row_sums))

print(f"\n{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)
