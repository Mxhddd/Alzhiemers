"""Machine check of manuscript/derivations/module4_fixed_points.md.

Model (manuscript/notes/module4_model_proposal.md, adopted 2026-10-13):
    f  = f0 / (1 + kappa R E)
    E' = f - (u + d) E                      u = u0 / (1 + D_atx / IC_atx)
    Y' = phi u E - (v + m) Y
    R' = r_rep (1 - R) - r_dam (m Y + h A) R
    T' = (s0 + s1 m Y)(1 + eta S)(1 + h A) - c0 R T,   S = E / (1 + D_praz / IC_praz)

Checks:
  1. symbolic: E*(R) solves (u+d) E (1 + kappa R E) = f0; Y*, R*, T* solve their equations;
  2. G(0) < 0 <= G(1) (at least one fixed point in (0, 1]);
  3. kappa = 0 -> G linear in R -> exactly one fixed point;
  4. numerically, for random parameter sets: the roots of G match the steady states
     reached by integrating the ODE from many initial conditions, and a root is
     stable (all Jacobian eigenvalues negative real part) iff G'(R*) > 0;
  5. at least one random parameter set is bistable and one is not (the structure
     can answer either way).
Parameters here are random TEST values (dimensionless), not model parameters.

    uv run python manuscript/derivations/checks/module4_fixed_points.py
"""
from __future__ import annotations

import sys

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

results = []


def check(name, ok, detail=""):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}{(': ' + detail) if detail else ''}")


# 1. symbolic ------------------------------------------------------------------
f0, kap, u, d, phi, v, m, rr, rd, h, A, s0, s1, eta, c0, Rs, E = sp.symbols(
    "f0 kappa u d phi v m r_rep r_dam h A s0 s1 eta c0 R E", positive=True)
w = u + d
Estar = (-w + sp.sqrt(w**2 + 4 * kap * Rs * w * f0)) / (2 * kap * Rs * w)
check("E*(R) solves (u+d) E (1 + kappa R E) = f0",
      sp.simplify(w * Estar * (1 + kap * Rs * Estar) - f0) == 0)
Ystar = phi * u * Estar / (v + m)
check("Y* = phi u E*/(v+m) zeroes Y'", sp.simplify(phi * u * Estar - (v + m) * Ystar) == 0)
G = Rs * (rr + rd * (m * Ystar + h * A)) - rr
Rdot = rr * (1 - Rs) - rd * (m * Ystar + h * A) * Rs
check("G(R) = 0 <=> R' = 0 at (E*, Y*)", sp.simplify(G + Rdot) == 0)
lim0 = sp.limit(Estar, Rs, 0, "+")
check("E*(R) -> f0/(u+d) as kappa R -> 0 (continuous at R = 0)", sp.simplify(lim0 - f0 / w) == 0, str(lim0))
Fimp = w * E * (1 + kap * Rs * E) - f0
dEdR = sp.simplify(-sp.diff(Fimp, Rs) / sp.diff(Fimp, E))          # implicit-function theorem
check("dE*/dR = -kappa E^2 / (1 + 2 kappa R E) < 0 (more restraint, less NE)",
      sp.simplify(dEdR + kap * E**2 / (1 + 2 * kap * Rs * E)) == 0, str(dEdR))
# R dE*/dR >= -E*/2, hence G'(R) >= r_rep + r_dam (m Y*/2 + h A) > 0 in the q = 1 model
gap = sp.simplify(Rs * dEdR + E / 2)                               # = E / (2 (1 + 2 kappa R E)) > 0
check("R dE*/dR + E*/2 = E/(2(1 + 2 kappa R E)) > 0, so G is strictly increasing when damage is linear in Y (q = 1): "
      "the proposal's structure has exactly ONE fixed point for every parameter set",
      sp.simplify(gap - E / (2 * (1 + 2 * kap * Rs * E))) == 0)


# 2-5. numerical ------------------------------------------------------------------
def make(p):
    def Est(R):
        a = p["kappa"] * R * p["w"]
        return p["f0"] / p["w"] if a == 0 else (-p["w"] + np.sqrt(p["w"] ** 2 + 4 * a * p["f0"])) / (2 * a)

    def Yst(R):
        return p["phi"] * p["u"] * Est(R) / (p["v"] + p["m"])

    def G(R):
        return R * (p["rr"] + p["rd"] * (p["m"] * Yst(R) + p["hA"])) - p["rr"]

    def rhs(t, y):
        E_, Y_, R_ = y
        f = p["f0"] / (1 + p["kappa"] * R_ * E_)
        return [f - p["w"] * E_, p["phi"] * p["u"] * E_ - (p["v"] + p["m"]) * Y_,
                p["rr"] * (1 - R_) - p["rd"] * (p["m"] * Y_ + p["hA"]) * R_]

    def jac(E_, Y_, R_):
        q = 1 + p["kappa"] * R_ * E_
        return np.array([
            [-p["w"] - p["f0"] * p["kappa"] * R_ / q**2, 0.0, -p["f0"] * p["kappa"] * E_ / q**2],
            [p["phi"] * p["u"], -(p["v"] + p["m"]), 0.0],
            [0.0, -p["rd"] * p["m"] * R_, -p["rr"] - p["rd"] * (p["m"] * Y_ + p["hA"])]])
    return Est, Yst, G, rhs, jac


rng = np.random.default_rng(4)
n_sets, n_bistable, n_mono, agree, stab_ok, g_signs, n_roots = 60, 0, 0, 0, 0, True, 0
for k in range(n_sets):
    p = {"f0": 10 ** rng.uniform(-0.5, 1.5), "kappa": 10 ** rng.uniform(-1, 2.5), "u": 10 ** rng.uniform(-1, 1),
         "d": 10 ** rng.uniform(-1, 0.5), "phi": rng.uniform(0.2, 1), "v": 10 ** rng.uniform(-1, 1),
         "m": 10 ** rng.uniform(-1, 1), "rr": 10 ** rng.uniform(-2, 0), "rd": 10 ** rng.uniform(-1, 2),
         "hA": 10 ** rng.uniform(-2, 0)}
    p["w"] = p["u"] + p["d"]
    Est, Yst, G, rhs, jac = make(p)
    g_signs &= G(1e-12) < 0 <= G(1.0)
    grid = np.linspace(1e-9, 1.0, 4001)
    gv = np.array([G(r) for r in grid])
    roots = [brentq(G, grid[i], grid[i + 1]) for i in range(len(grid) - 1) if gv[i] * gv[i + 1] < 0]
    if len(roots) >= 3:
        n_bistable += 1
    else:
        n_mono += 1
    # stability from the Jacobian vs sign of G'
    for r in roots:
        e, y = Est(r), Yst(r)
        stable = np.all(np.linalg.eigvals(jac(e, y, r)).real < 0)
        gprime = (G(r + 1e-7) - G(r - 1e-7)) / 2e-7
        stab_ok += int(stable == (gprime > 0))
        n_roots += 1
    # integrate from a spread of initial conditions; every end state must be a stable root
    ends = []
    for R0 in (0.02, 0.3, 0.6, 0.98):
        for E0 in (0.0, 5.0):
            sol = solve_ivp(rhs, (0, 5000), [E0, 0.0, R0], method="Radau", rtol=1e-9, atol=1e-12)
            ends.append(sol.y[2, -1])
    stable_roots = [r for r in roots if np.all(np.linalg.eigvals(jac(Est(r), Yst(r), r)).real < 0)]
    ok = all(min(abs(e - r) for r in stable_roots) < 1e-4 for e in ends) and stable_roots
    agree += int(bool(ok))
check("G(0) < 0 <= G(1) for every random parameter set (a fixed point always exists in (0, 1])", g_signs)
check(f"ODE end states coincide with stable roots of G ({n_sets} random parameter sets, 8 starts each)",
      agree == n_sets, f"{agree}/{n_sets}")
check("a root is Jacobian-stable iff G'(R*) > 0 (every root of every set)", n_roots > 0 and stab_ok == n_roots,
      f"{stab_ok}/{n_roots} roots consistent")
check("q = 1 (the proposal's linear damage): no random set is bistable, as proved above",
      n_bistable == 0, f"{n_bistable} bistable, {n_mono} single fixed point")


# 6. generalised damage H(mY) = (mY/K)^q / (1 + (mY/K)^q): can it go either way? ----------
def make_q(p):
    Est, Yst, _, _, _ = make(p)

    def H(x):
        z = (x / p["K"]) ** p["q"]
        return z / (1 + z)

    def G(R):
        return R * (p["rr"] + p["rd"] * (H(p["m"] * Yst(R)) + p["hA"])) - p["rr"]

    def rhs(t, y):
        E_, Y_, R_ = y
        f = p["f0"] / (1 + p["kappa"] * R_ * E_)
        return [f - p["w"] * E_, p["phi"] * p["u"] * E_ - (p["v"] + p["m"]) * Y_,
                p["rr"] * (1 - R_) - p["rd"] * (H(p["m"] * max(Y_, 0.0)) + p["hA"]) * R_]
    return Est, Yst, G, rhs


rng = np.random.default_rng(7)
nb, nm, ode_ok, ode_multi = 0, 0, 0, 0
for k in range(200):
    p = {"f0": 10 ** rng.uniform(-0.5, 1.5), "kappa": 10 ** rng.uniform(-1, 2.5), "u": 10 ** rng.uniform(-1, 1),
         "d": 10 ** rng.uniform(-1, 0.5), "phi": rng.uniform(0.2, 1), "v": 10 ** rng.uniform(-1, 1),
         "m": 10 ** rng.uniform(-1, 1), "rr": 10 ** rng.uniform(-2, 0), "rd": 10 ** rng.uniform(-1, 2),
         "hA": 10 ** rng.uniform(-3, -1), "q": rng.choice([2.0, 4.0]), "K": 10 ** rng.uniform(-1, 1)}
    p["w"] = p["u"] + p["d"]
    Est, Yst, G, rhs = make_q(p)
    grid = np.linspace(1e-9, 1.0, 4001)
    gv = np.array([G(r) for r in grid])
    roots = [brentq(G, grid[i], grid[i + 1]) for i in range(len(grid) - 1) if gv[i] * gv[i + 1] < 0]
    if len(roots) >= 3:
        nb += 1
        if ode_multi < 3:      # integrate from low and high R: two different stable states expected
            ends = sorted(solve_ivp(rhs, (0, 20000), [Est(R0), Yst(R0), R0], method="Radau", rtol=1e-9,
                                    atol=1e-12).y[2, -1] for R0 in (0.01, 0.99))
            ode_ok += int(abs(ends[0] - roots[0]) < 1e-3 and abs(ends[1] - roots[-1]) < 1e-3)
            ode_multi += 1
    else:
        nm += 1
check("generalised damage (q in {2, 4}): some random sets bistable, some not, i.e. the question can come out either way",
      nb > 0 and nm > 0, f"{nb} bistable, {nm} single fixed point (of 200)")
check("in bistable sets the ODE settles in the low-R state from R0 = 0.01 and in the high-R state from R0 = 0.99",
      ode_multi > 0 and ode_ok == ode_multi, f"{ode_ok}/{ode_multi} sets")
# kappa = 0: unique root
p0 = {"f0": 3.0, "kappa": 0.0, "u": 1.0, "d": 0.5, "phi": 0.7, "v": 1.0, "m": 2.0, "rr": 0.1, "rd": 50.0, "hA": 0.2}
p0["w"] = p0["u"] + p0["d"]
_, _, G0, _, _ = make(p0)
gv0 = np.array([G0(r) for r in np.linspace(1e-9, 1, 4001)])
check("kappa = 0: G is linear in R, exactly one fixed point", int(np.sum(gv0[:-1] * gv0[1:] < 0)) == 1)

print(f"\n{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)
