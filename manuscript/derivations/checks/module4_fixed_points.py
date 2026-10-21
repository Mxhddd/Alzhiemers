"""Machine check of manuscript/derivations/module4_fixed_points.md (and the
model of manuscript/derivations/module4_lc_model.md).

MACHINE-CHECKED IS NOT HUMAN-CHECKED. This script re-derives the algebra with
sympy and tests every claim numerically. It cannot tell whether the model is a
good description of the biology.

Model (dimensionless; doses D in units of IC50):
    f  = 1 / (1 + kappa R E)
    E' = f - w E,           w = u + (1 - u0),  u = u0 / (1 + D_atx)
    Y' = phi u E - (v + m) Y
    R' = r_rep (1 - R) - r_dam (H(m Y) + hA) R,   H(x) = (x/K)^q / (1 + (x/K)^q)
    T' = c0 [ (1 + sigma m Y)(1 + eta s E)(1 + hA) - c R T ],
         s = 1/(1 + D_praz),  c = 1 - eps + eps s

Part A is symbolic. Part B is numerical and self-contained (its own right-hand
side; random TEST parameter sets drawn from a deliberately wide distribution,
not the registered parameter box, with early-exit existence searches so that no
registered fraction is estimated here). Part C checks that the implementation
in src/dynamics/lc_model.py agrees with Part B's independent code.

    uv run python manuscript/derivations/checks/module4_fixed_points.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

results = []


def check(name, ok, detail=""):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}{(': ' + detail) if detail else ''}")


# =========================================================================== A. symbolic
kap, w, R, E, phi, u, v, m, rr, rd, hA, K, q, H, x = sp.symbols(
    "kappa w R E phi u v m r_rep r_dam hA K q H x", positive=True)
sig, eta, eps, s, c, c0, d, u0, D = sp.symbols("sigma eta epsilon s c c0 d u0 D", positive=True)

# A1. E*(R)
Estar = 2 / (w + sp.sqrt(w**2 + 4 * kap * R * w))
check("A1  E*(R) = 2/(w + sqrt(w^2 + 4 kappa R w)) solves w E (1 + kappa R E) = 1",
      sp.simplify(sp.radsimp(w * Estar * (1 + kap * R * Estar) - 1)) == 0)

# A2. dE*/dR and its elasticity (implicit-function theorem)
F = 1 / (1 + kap * R * E) - w * E
dEdR = -sp.diff(F, R) / sp.diff(F, E)
dEdR_fp = sp.simplify(dEdR.subs(w, 1 / (E * (1 + kap * R * E))))   # use w = f/E at the fixed point
check("A2  at a fixed point dE*/dR = -kappa E^2 / (1 + 2 kappa R E)",
      sp.simplify(dEdR_fp + kap * E**2 / (1 + 2 * kap * R * E)) == 0, str(dEdR_fp))
el = sp.simplify(-R * dEdR_fp / E)
check("A2b elasticity eps_E = -R E*'/E* = kappa R E/(1 + 2 kappa R E), and 1/2 - eps_E = 1/(2(1 + 2 kappa R E)) > 0",
      sp.simplify(el - kap * R * E / (1 + 2 * kap * R * E)) == 0
      and sp.simplify(sp.Rational(1, 2) - el - 1 / (2 * (1 + 2 * kap * R * E))) == 0)

# A3. G(R) and its derivative: G'(R) = r_rep + r_dam [H + hA - eps_E q H (1 - H)]
Rf = sp.Function("Rf")
Hill = lambda z: (z / K)**q / (1 + (z / K)**q)  # noqa: E731
xR = m * phi * u * Estar / (v + m)
Gs = R * (rr + rd * (Hill(xR) + hA)) - rr
Rdot = rr * (1 - R) - rd * (Hill(xR) + hA) * R
check("A3  G(R) = 0 <=> R' = 0 at (E*(R), Y*(R))", sp.simplify(Gs + Rdot) == 0)
Gp_sym = sp.diff(Gs, R)
Hx = Hill(xR)
Ex = Estar
Gp_formula = rr + rd * (Hx + hA - (kap * R * Ex / (1 + 2 * kap * R * Ex)) * q * Hx * (1 - Hx))
rng = np.random.default_rng(11)


def rat(z):
    """Exact rational test value (so only the algebra, not float rounding, is compared)."""
    return sp.Rational(str(round(float(z), 12)))


worst = 0.0
for _ in range(40):
    vals = {kap: 10**rng.uniform(-1, 2), w: rng.uniform(0.2, 2), R: rng.uniform(0.01, 1), phi: rng.uniform(0.1, 1),
            u: rng.uniform(0.05, 0.9), v: 10**rng.uniform(-1, 1), m: 10**rng.uniform(-1, 1), rr: 10**rng.uniform(-3, -1),
            rd: 10**rng.uniform(-3, 0), hA: rng.uniform(0, 1), K: 10**rng.uniform(-2, 0), q: rng.uniform(1, 4)}
    vals = {k_: rat(v_) for k_, v_ in vals.items()}
    a, b = sp.N(Gp_sym.subs(vals), 40), sp.N(Gp_formula.subs(vals), 40)
    worst = max(worst, float(abs(a - b) / max(1e-30, abs(a))))
check("A3b G'(R) = r_rep + r_dam [H + hA - eps_E q H(1 - H)] (sympy derivative vs formula, 40 random rational points, 40 digits)",
      worst < 1e-20, f"max relative difference {worst:.1e}")

# A4. characteristic polynomial of the (E, Y, R) block and Routh-Hurwitz
a11, a13, a21, a22, a32, a33, lam = sp.symbols("a11 a13 a21 a22 a32 a33 lambda", positive=True)
J = sp.Matrix([[-a11, 0, -a13], [a21, -a22, 0], [0, -a32, -a33]])
p3 = sp.expand((lam * sp.eye(3) - J).det())
target = sp.expand((lam + a11) * (lam + a22) * (lam + a33) - a13 * a21 * a32)
check("A4  det(lambda I - J) = (lambda + a11)(lambda + a22)(lambda + a33) - a13 a21 a32 (one positive feedback loop E->Y->R->E)",
      sp.simplify(p3 - target) == 0)
c1 = a11 + a22 + a33
c2 = a11 * a22 + a11 * a33 + a22 * a33
c3 = a11 * a22 * a33 - a13 * a21 * a32
hur = sp.expand(c1 * c2 - c3)
pos = sp.expand(a11**2 * a22 + a11**2 * a33 + a11 * a22**2 + a22**2 * a33 + a11 * a33**2 + a22 * a33**2
                + 2 * a11 * a22 * a33 + a13 * a21 * a32)
check("A4b c1 c2 - c3 = (a sum of positive terms) > 0 always: the Hurwitz determinant never vanishes, so NO Hopf bifurcation; "
      "stability <=> c3 > 0", sp.simplify(hur - pos) == 0)
# c3 = a11 a22 G'(R): with E*' = -a13/a11, Y*' = (a21/a22) E*', G' = a33 + a32 Y*'
Gp_chain = a33 + a32 * (a21 / a22) * (-a13 / a11)
check("A4c c3 = a11 a22 G'(R*): a fixed point is stable iff G'(R*) > 0 (and a saddle with one positive eigenvalue iff G' < 0)",
      sp.simplify(c3 - a11 * a22 * Gp_chain) == 0)

# A5. q <= 2 => G' > 0 (monostable). For eps_E in [0, 1/2], q in [1, 2], H in [0, 1]:
#     H - eps_E q H (1 - H) >= H - H (1 - H) = H^2 >= 0.
e_ = sp.symbols("e_", nonnegative=True)
lhs = H - e_ * q * H * (1 - H)
worst_case = lhs.subs({e_: sp.Rational(1, 2), q: 2})
check("A5  at the worst case eps_E = 1/2, q = 2: H - eps_E q H (1 - H) = H^2 >= 0, and the expression decreases in eps_E and q, "
      "so q <= 2 => G' > r_rep > 0: exactly ONE fixed point for every parameter set with q <= 2",
      sp.simplify(worst_case - H**2) == 0 and sp.simplify(sp.diff(lhs, e_) + q * H * (1 - H)) == 0
      and sp.simplify(sp.diff(lhs, q) + e_ * H * (1 - H)) == 0)

# A6. necessary condition for bistability: G' < 0 somewhere needs
#     r_rep/r_dam + hA < max_H H (q (1 - H)/2 - 1) = (q - 2)^2 / (8 q)   (q > 2)
phiH = H * (q * (1 - H) / 2 - 1)
Hopt = sp.solve(sp.diff(phiH, H), H)[0]
check("A6  max over H of H(q(1 - H)/2 - 1) is (q - 2)^2/(8q) at H = 1/2 - 1/q: bistability at hA needs q > 2 and "
      "r_rep/r_dam + hA < (q - 2)^2/(8q) (<= 1/8 for q <= 4)",
      sp.simplify(Hopt - (sp.Rational(1, 2) - 1 / q)) == 0
      and sp.simplify(phiH.subs(H, Hopt) - (q - 2)**2 / (8 * q)) == 0
      and sp.simplify(((q - 2)**2 / (8 * q)).subs(q, 4) - sp.Rational(1, 8)) == 0)

# A7. Phi(R): G = r_dam R (hA - Phi(R)); Phi'(R) = [-rho/R + eps_E q H (1 - H)]/R
rho = rr / rd
PhiR = rho * (1 - R) / R - Hill(xR)
check("A7  G(R) = r_dam R (hA - Phi(R)) with Phi(R) = (r_rep/r_dam)(1 - R)/R - H(x*(R))",
      sp.simplify(Gs - rd * R * (hA - PhiR)) == 0)
Phip_formula = (-rho / R + (kap * R * Ex / (1 + 2 * kap * R * Ex)) * q * Hx * (1 - Hx)) / R
worst = 0.0
for _ in range(30):
    vals = {kap: 10**rng.uniform(-1, 2), w: rng.uniform(0.2, 2), R: rng.uniform(0.01, 1), phi: rng.uniform(0.1, 1),
            u: rng.uniform(0.05, 0.9), v: 10**rng.uniform(-1, 1), m: 10**rng.uniform(-1, 1), rr: 10**rng.uniform(-3, -1),
            rd: 10**rng.uniform(-3, 0), K: 10**rng.uniform(-2, 0), q: rng.uniform(1, 4)}
    vals = {k_: rat(v_) for k_, v_ in vals.items()}
    a, b = sp.N(sp.diff(PhiR, R).subs(vals), 40), sp.N(Phip_formula.subs(vals), 40)
    worst = max(worst, float(abs(a - b) / max(1e-30, abs(a))))
check("A7b Phi'(R) = [-r_rep/(r_dam R) + eps_E q H(1 - H)]/R (30 random rational points, 40 digits)", worst < 1e-20, f"{worst:.1e}")

# A8. readouts do not depend on c0, and depend on r_rep, r_dam only through their ratio
Tstar = (1 + sig * x) * (1 + eta * s * E) * (1 + hA) / (c * R)
k = sp.symbols("k", positive=True)
check("A8  T* = (1 + sigma x)(1 + eta s E)(1 + hA)/(c R) contains no c0; G(R; k r_rep, k r_dam) = k G(R): "
      "the fixed points depend on r_rep and r_dam only through r_rep/r_dam",
      c0 not in Tstar.free_symbols and sp.simplify(Gs.subs({rr: k * rr, rd: k * rd}, simultaneous=True) - k * Gs) == 0)

# A9. prazosin: E, Y, R do not depend on D_praz, so T*(D)/T*(0) = (1 + eta s E)/(1 + eta E) * 1/(1 - eps + eps s)
ratio = (1 + eta * s * E) / ((1 + eta * E) * (1 - eps + eps * s))
dds = sp.simplify(sp.diff(ratio, s))
num = sp.factor(sp.numer(sp.together(dds)))
check("A9  prazosin: d ratio/ds has the sign of eta E (1 - eps) - eps (s = 1/(1 + D) falls with dose): the dose-response is "
      "MONOTONE (never U-shaped); protective iff eta E (1 - eps) > eps, harmful iff eta E (1 - eps) < eps",
      sp.simplify(num / (eta * E * (1 - eps) - eps)).is_positive or
      sp.simplify(num / (eta * E * (1 - eps) - eps)).free_symbols <= {eta, E} and
      sp.simplify(dds * (1 + eta * E) * (1 - eps + eps * s)**2 - (eta * E * (1 - eps) - eps)) == 0,
      f"d ratio/ds = {dds}")
check("A9b prazosin at D = 1 (s = 1/2): ratio = (1 + eta E/2) / ((1 + eta E)(1 - eps/2))",
      sp.simplify(ratio.subs(s, sp.Rational(1, 2)) - (1 + eta * E / 2) / ((1 + eta * E) * (1 - eps / 2))) == 0)

# A10. atomoxetine, open loop (kappa = 0 and r_dam = 0, so R = 1): with h = 1/w (rises with dose)
#      T*(h) = (1 + sigma mu phi (1 - d h))(1 + eta h)(1 + hA)/c is a concave quadratic in h:
#      the minimum over any dose range is at an end, so NO U-shape is possible.
hh, mu = sp.symbols("h mu", positive=True)
T_open = (1 + sig * mu * phi * (1 - d * hh)) * (1 + eta * hh) * (1 + hA) / c
check("A10 open loop (kappa = r_dam = 0): uE = u/w = 1 - d h, so T*(h) has d^2T/dh^2 = -2 sigma mu phi d eta (1 + hA)/c < 0: "
      "concave, hence atomoxetine can be harmful at an intermediate dose (inverted U) but NEVER U-shaped without the feedbacks",
      sp.simplify(sp.diff(T_open, hh, 2) + 2 * sig * mu * phi * d * eta * (1 + hA) / c) == 0)

# A11. the validation fixed point: kappa = 4, u0 = 1/2 (w = 1 without drug, u = 1/2), phi = 1, v = m = 1,
#      K = 1/8, q = 3, r_rep = r_dam = 1/20, hA = 1/2, sigma = 8, eta = 2, no drug (s = c = 1)
V = {kap: 4, w: 1, u: sp.Rational(1, 2), phi: 1, v: 1, m: 1, K: sp.Rational(1, 8), q: 3, rr: sp.Rational(1, 20),
     rd: sp.Rational(1, 20), hA: sp.Rational(1, 2)}
E_v = sp.nsimplify(Estar.subs(V).subs(R, sp.Rational(1, 2)))
x_v = sp.nsimplify(xR.subs(V).subs(R, sp.Rational(1, 2)))
G_v = sp.nsimplify(Gs.subs(V).subs(R, sp.Rational(1, 2)))
T_v = sp.nsimplify(Tstar.subs({sig: 8, eta: 2, s: 1, c: 1, hA: sp.Rational(1, 2), x: x_v, E: E_v, R: sp.Rational(1, 2)}))
Y_v = x_v / 1
Hv = Hill(x_v).subs({K: sp.Rational(1, 8), q: 3})
A = {a11: 1 + 4 * sp.Rational(1, 2) / (1 + 4 * sp.Rational(1, 2) * E_v)**2, a13: 4 * E_v / (1 + 4 * sp.Rational(1, 2) * E_v)**2,
     a21: sp.Rational(1, 2), a22: 2, a32: sp.Rational(1, 20) * sp.Rational(1, 2) * 1 * 3 * Hv * (1 - Hv) / x_v,
     a33: sp.Rational(1, 20) + sp.Rational(1, 20) * (Hv + sp.Rational(1, 2))}
c3v = sp.nsimplify(c3.subs(A))
bound = sp.Rational(1, 1) + sp.Rational(1, 2)       # r_rep/r_dam + hA
check("A11 validation point: (E, Y, R, T) = (1/2, 1/8, 1/2, 12) exactly; c3 = 21/80 > 0 (stable); "
      "r_rep/r_dam + hA = 3/2 > (q - 2)^2/(8q) = 1/24, so it is the ONLY fixed point",
      E_v == sp.Rational(1, 2) and Y_v == sp.Rational(1, 8) and G_v == 0 and T_v == 12 and c3v == sp.Rational(21, 80)
      and bound > sp.Rational(1, 24), f"E={E_v}, Y={Y_v}, G(1/2)={G_v}, T={T_v}, c3={c3v}")


# =========================================================================== B. numerical (self-contained)
def P_derived(p):
    uu = p["u0"] / (1 + p.get("D_atx", 0.0))
    ss = 1 / (1 + p.get("D_praz", 0.0))
    return uu, uu + 1 - p["u0"], ss, 1 - p["eps"] + p["eps"] * ss


def H_(z, p):
    t = (max(z, 0.0) / p["K"]) ** p["q"]
    return t / (1 + t)


def f_rhs(t, y, p):
    E_, Y_, R_, T_ = y
    uu, ww, ss, cc = P_derived(p)
    fr = 1 / (1 + p["kappa"] * max(R_, 0) * max(E_, 0))
    xx = p["m"] * Y_
    return [fr - ww * E_, p["phi"] * uu * E_ - (p["v"] + p["m"]) * Y_,
            p["r_rep"] * (1 - R_) - p["r_dam"] * (H_(xx, p) + p["hA"]) * R_,
            p["c0"] * ((1 + p["sigma"] * xx) * (1 + p["eta"] * ss * E_) * (1 + p["hA"]) - cc * R_ * T_)]


def fd_jac(y, p, h=1e-7):
    Jm = np.empty((4, 4))
    for j in range(4):
        e = np.zeros(4); e[j] = h * max(1, abs(y[j]))
        Jm[:, j] = (np.array(f_rhs(0, y + e, p)) - np.array(f_rhs(0, y - e, p))) / (2 * e[j])
    return Jm


def G_num(Rv, p):
    uu, ww, _, _ = P_derived(p)
    Ev = 2 / (ww + np.sqrt(ww * ww + 4 * p["kappa"] * Rv * ww))
    xv = p["m"] * p["phi"] * uu * Ev / (p["v"] + p["m"])
    return Rv * (p["r_rep"] + p["r_dam"] * (H_(xv, p) + p["hA"])) - p["r_rep"]


def fp_at(Rv, p):
    uu, ww, ss, cc = P_derived(p)
    Ev = 2 / (ww + np.sqrt(ww * ww + 4 * p["kappa"] * Rv * ww))
    Yv = p["phi"] * uu * Ev / (p["v"] + p["m"])
    xv = p["m"] * Yv
    return np.array([Ev, Yv, Rv, (1 + p["sigma"] * xv) * (1 + p["eta"] * ss * Ev) * (1 + p["hA"]) / (cc * Rv)])


def roots_num(p, n=20001):
    """Independent root finder: a dense uniform-in-log + linear grid, no derivative refinement."""
    grid = np.unique(np.concatenate([np.logspace(-9, -2, 3000), np.linspace(1e-2, 1, n)]))
    gv = np.array([G_num(r, p) for r in grid])
    return sorted([float(r) for r, g in zip(grid, gv) if g == 0.0] +
                  [brentq(G_num, grid[i], grid[i + 1], args=(p,), xtol=1e-15) for i in range(len(grid) - 1)
                   if gv[i] * gv[i + 1] < 0])


def wide_sample(rng, **over):
    p = dict(kappa=10**rng.uniform(-2, 3), u0=rng.uniform(0.05, 0.95), phi=rng.uniform(0.05, 1), v=10**rng.uniform(-2, 2),
             m=10**rng.uniform(-2, 2), r_rep=10**rng.uniform(-4, 0), r_dam=10**rng.uniform(-4, 1), K=10**rng.uniform(-3, 1),
             q=rng.uniform(1, 6), sigma=10**rng.uniform(-2, 3), eta=10**rng.uniform(-3, 2), eps=rng.uniform(0, 0.8),
             hA=rng.uniform(0, 1.5), c0=10**rng.uniform(-3, -1))
    p.update(over)
    return p


def find_bistable(rng, tries=20000):
    """Early-exit search for a parameter set with three fixed points (biased towards the region allowed by A6)."""
    for _ in range(tries):
        p = wide_sample(rng, q=rng.uniform(3, 6), kappa=10**rng.uniform(1, 3), r_rep=10**rng.uniform(-4, -2),
                        r_dam=10**rng.uniform(-1, 1), hA=rng.uniform(0, 0.1))
        if len(roots_num(p, 4001)) >= 3:
            return p
    return None


rng = np.random.default_rng(2026)
# B1. analytic fixed points = 4-D Newton fixed points; stability G' > 0 <=> Jacobian stable
n_sets, n_fp, ok_match, ok_stab, ok_exist = 40, 0, 0, 0, True
sets = [wide_sample(rng) for _ in range(n_sets - 4)]
for _ in range(4):
    pb = find_bistable(rng)
    if pb is not None:
        sets.append(pb)
for p in sets:
    rts = roots_num(p)
    ok_exist &= len(rts) >= 1 and G_num(1e-12, p) < 0 <= G_num(1.0, p)
    for r in rts:
        y = fp_at(r, p)
        sol = root(lambda z: f_rhs(0, z, p), y * (1 + 1e-3), method="hybr", tol=1e-13)
        ok_match += int(np.max(np.abs(sol.x - y) / np.maximum(1, np.abs(y))) < 1e-7)
        gp = (G_num(r + 1e-8, p) - G_num(r - 1e-8, p)) / 2e-8
        eig = np.linalg.eigvals(fd_jac(y, p))
        ok_stab += int((np.max(eig.real) < 0) == (gp > 0))
        n_fp += 1
check("B1  G(0+) < 0 <= G(1) in every set: at least one fixed point always exists in (0, 1]", ok_exist)
check(f"B1b each root of G is a fixed point of the 4-D system (Newton from a perturbed start), {len(sets)} sets "
      f"({len(sets) - (n_sets - 4)} of them bistable)", ok_match == n_fp, f"{ok_match}/{n_fp}")
check("B1c Jacobian-stable <=> G'(R*) > 0 at every fixed point", ok_stab == n_fp, f"{ok_stab}/{n_fp}")

# B2. from the healthy start the ODE goes MONOTONICALLY (E, Y, T up; R down) to the fixed point with the largest R*
ok_lim, ok_mono, n_b = 0, 0, 0
for p in sets[:12] + sets[n_sets - 4:]:
    rts = roots_num(p)
    top = fp_at(max(rts), p)
    uu, ww, ss, cc = P_derived(p)
    gp = (G_num(max(rts) + 1e-8, p) - G_num(max(rts) - 1e-8, p)) / 2e-8
    slow = min(p["c0"] * cc * max(rts), abs(np.max(np.linalg.eigvals(fd_jac(top, p)).real)))
    sol = solve_ivp(f_rhs, (0, 60 / slow), [0, 0, 1, 0], args=(p,), method="Radau", rtol=1e-10, atol=1e-13,
                    dense_output=True)
    ok_lim += int(np.max(np.abs(sol.y[:, -1] - top) / np.maximum(1e-9, np.abs(top))) < 1e-5)
    dy = np.diff(sol.y, axis=1)
    tol = 1e-9 * np.maximum(1, np.abs(sol.y[:, 1:]))
    ok_mono += int(np.all(dy[[0, 1, 3]] >= -tol[[0, 1, 3]]) and np.all(dy[2] <= tol[2]))
    n_b += 1
check("B2  the trajectory from (0, 0, 1, 0) converges to the fixed point with the LARGEST R* (incl. bistable sets)",
      ok_lim == n_b, f"{ok_lim}/{n_b}")
check("B2b and it is monotone: E, Y, T never decrease and R never increases (cooperative system, derivation §4)",
      ok_mono == n_b, f"{ok_mono}/{n_b}")

# B3. in bistable sets the low-R fixed point is reached from a damaged start
pb = sets[-1]
rts = roots_num(pb)
low = fp_at(min(rts), pb)
sol = solve_ivp(f_rhs, (0, 1e6), [low[0], low[1], min(rts) * 0.9, 0], args=(pb,), method="Radau", rtol=1e-10, atol=1e-13)
check("B3  bistable set: three roots, the middle one a saddle (G' < 0), and a damaged start goes to the low-R state",
      len(rts) == 3 and (G_num(rts[1] + 1e-8, pb) - G_num(rts[1] - 1e-8, pb)) < 0
      and abs(sol.y[2, -1] - min(rts)) < 1e-5, f"roots {np.round(rts, 6)}")

# B4. the necessary condition of A6 holds in every bistable set found
viol = 0
for p in sets[n_sets - 4:]:
    viol += int(not (p["q"] > 2 and p["r_rep"] / p["r_dam"] + p["hA"] < (p["q"] - 2)**2 / (8 * p["q"])))
check("B4  every bistable set found satisfies q > 2 and r_rep/r_dam + hA < (q - 2)^2/(8q)", viol == 0 and len(sets) > n_sets - 4)
# and q <= 2 never gives three roots (random search)
n_q2 = 0
for _ in range(300):
    p = wide_sample(rng, q=rng.uniform(1, 2), kappa=10**rng.uniform(1, 3), r_rep=10**rng.uniform(-4, -2),
                    r_dam=10**rng.uniform(-1, 1), hA=rng.uniform(0, 0.1))
    n_q2 += int(len(roots_num(p, 4001)) >= 3)
check("B4b q <= 2: no set with more than one fixed point (300 sets drawn where bistability is otherwise easiest)", n_q2 == 0)

# B5. the bistable window in hA from Phi(R) agrees with a fine scan of hA
p = dict(sets[-1])
grid = np.unique(np.concatenate([np.logspace(-9, -2, 3000), np.linspace(1e-2, 1, 20001)]))
uu, ww, _, _ = P_derived(p)
Eg = 2 / (ww + np.sqrt(ww * ww + 4 * p["kappa"] * grid * ww))
xg = p["m"] * p["phi"] * uu * Eg / (p["v"] + p["m"])
Phi_g = (p["r_rep"] / p["r_dam"]) * (1 - grid) / grid - np.array([H_(z, p) for z in xg])
inc = np.where(np.diff(Phi_g) > 0)[0]
w_lo, w_hi = Phi_g[inc[0]], Phi_g[inc[-1] + 1]
scan = np.linspace(w_lo - 0.2 * (w_hi - w_lo), w_hi + 0.2 * (w_hi - w_lo), 201)   # hA < 0 allowed: algebra only
three = [h for h in scan if len(roots_num(dict(p, hA=h), 4001)) >= 3]
check("B5  the hA window with three fixed points equals (local min of Phi, local max of Phi) to scan resolution",
      len(three) > 0 and abs(min(three) - w_lo) <= scan[1] - scan[0] and abs(max(three) - w_hi) <= scan[1] - scan[0],
      f"Phi window ({w_lo:.5f}, {w_hi:.5f}); scan ({min(three):.5f}, {max(three):.5f})")


# B6. each question can come out either way (early-exit existence searches, wide TEST distribution)
def ratio_num(p, key, dose):
    t0 = fp_at(max(roots_num(p, 4001)), p)[3]
    q_ = dict(p, **{key: dose})
    return fp_at(max(roots_num(q_, 4001)), q_)[3] / t0


def exists(pred, tries=3000):
    for _ in range(tries):
        p = wide_sample(rng)
        if pred(p):
            return True
    return False


check("B6a atomoxetine can be protective (ratio < 0.95 at D = 1)", exists(lambda p: ratio_num(p, "D_atx", 1.0) < 0.95))
check("B6b atomoxetine can be harmful (ratio > 1.05 at D = 1)", exists(lambda p: ratio_num(p, "D_atx", 1.0) > 1.05))
check("B6c prazosin can be protective", exists(lambda p: ratio_num(p, "D_praz", 1.0) < 0.95))
check("B6d prazosin can be harmful (needs eps > 0, A9)", exists(lambda p: ratio_num(p, "D_praz", 1.0) > 1.05))
dgrid = np.r_[0, np.logspace(-2, 1, 41)]


def ushape(p):
    y = np.array([fp_at(max(roots_num(dict(p, D_atx=dd), 2001)), dict(p, D_atx=dd))[3] for dd in dgrid])
    i = int(np.argmin(y))
    return 0 < i < len(y) - 1 and y[0] - y[i] > 0.05 * y[i] and y[-1] - y[i] > 0.05 * y[i]


# The U-shape is rare: an early-exit search over the wide test distribution (seed 99, outside this script)
# first met one at draw 1,681 (an inverted U at draw 2). That set is checked here, so existence is shown
# deterministically. It lies outside the registered box (eta = 16 > 10).
P_U = {"kappa": 11.376078531577555, "u0": 0.7237681857061803, "phi": 0.7564661794999585, "v": 7.133389055662706,
       "m": 0.7973803965681399, "r_rep": 0.006536740773393628, "r_dam": 0.020060467193409752, "K": 0.0204807465537441,
       "q": 3.6869063593778226, "sigma": 1.2835260133821838, "eta": 16.034179557899343, "eps": 0.2494542512028148,
       "hA": 0.8707037023800883, "c0": 0.01}
check("B6e atomoxetine CAN be U-shaped (a found example; needs the R feedback, A10)", ushape(P_U))


def inv_u(p):
    y = np.array([fp_at(max(roots_num(dict(p, D_atx=dd), 2001)), dict(p, D_atx=dd))[3] for dd in dgrid])
    i = int(np.argmax(y))
    return 0 < i < len(y) - 1 and y[i] - y[0] > 0.05 * y[i] and y[i] - y[-1] > 0.05 * y[i]


check("B6f atomoxetine can be harmful at an intermediate dose only (inverted U)", exists(inv_u, tries=200))
p_open = dict(P_U, kappa=0.0, r_dam=0.0)
check("B6g the same set with the feedbacks removed (kappa = r_dam = 0) is NOT U-shaped (A10)", not ushape(p_open))

# =========================================================================== C. implementation agrees
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from src.dynamics import lc_model as L  # noqa: E402

ok_rhs = ok_jac = ok_roots = ok_T = 0
for p in sets:
    y = np.array([rng.uniform(0, 1.5), rng.uniform(0, 0.5), rng.uniform(0.05, 1), rng.uniform(0, 20)])
    ok_rhs += int(np.allclose(L.rhs(0, y, p), f_rhs(0, y, p), rtol=1e-12, atol=1e-14))
    ok_jac += int(np.allclose(L.jac(0, y, p), fd_jac(y, p), rtol=1e-5, atol=1e-7))
    a, b = L.roots_of_G(p), np.array(roots_num(p))
    ok_roots += int(len(a) == len(b) and np.allclose(a, b, rtol=1e-9, atol=1e-12))
    ok_T += int(np.isclose(L.T_star(p), fp_at(max(b), p)[3], rtol=1e-9))
n = len(sets)
check("C1  src/dynamics/lc_model.rhs equals this script's right-hand side", ok_rhs == n, f"{ok_rhs}/{n}")
check("C2  src/dynamics/lc_model.jac equals the finite-difference Jacobian", ok_jac == n, f"{ok_jac}/{n}")
check("C3  src/dynamics/lc_model.roots_of_G finds the same roots as this script's independent dense-grid search",
      ok_roots == n, f"{ok_roots}/{n}")
check("C4  src/dynamics/lc_model.T_star (the registered readout) equals T* at the largest root", ok_T == n, f"{ok_T}/{n}")
pv = dict(kappa=4.0, u0=0.5, phi=1.0, v=1.0, m=1.0, K=0.125, q=3.0, r_rep=0.05, r_dam=0.05, hA=0.5, sigma=8.0,
          eta=2.0, eps=0.25, c0=0.01)
fps = L.fixed_points(pv)
check("C5  implementation at the validation point: one fixed point (1/2, 1/8, 1/2, 12), stable",
      len(fps) == 1 and np.allclose(fps[0]["y"], [0.5, 0.125, 0.5, 12.0], atol=1e-12) and fps[0]["stable"])

print(f"\n{sum(results)}/{len(results)} checks passed (machine-checked, NOT human-checked)")
sys.exit(0 if all(results) else 1)
