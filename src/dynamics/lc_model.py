"""The LC rate-controller ODE model (Module 4).

Equations, units and the meaning of every term: manuscript/derivations/module4_lc_model.md.
Fixed points and stability (derived by Claude, machine-checked, NOT human-checked):
manuscript/derivations/module4_fixed_points.md, checked by
manuscript/derivations/checks/module4_fixed_points.py.

Dimensionless model. Time is in units of the baseline extracellular
noradrenaline (NE) turnover time 1/w0; E and Y in units of the unrestrained,
drug-free extracellular NE level; R is a fraction in [0, 1]; T is in units of
the tau burden of a fully regulated (R = 1) system with no second hit and no
NE-dependent seeding or harm.

    f  = 1 / (1 + kappa R E)                             firing, alpha2A restraint
    E' = f - w E,  w = u + (1 - u0),  u = u0 / (1 + D_atx)
    Y' = phi u E - (v + m) Y
    R' = r_rep (1 - R) - r_dam (H(m Y) + hA) R,   H(x) = (x/K)^q / (1 + (x/K)^q)
    T' = c0 [ (1 + sigma m Y)(1 + eta S)(1 + hA) - c R T ]
         S = E / (1 + D_praz),  c = 1 - eps + eps / (1 + D_praz)

Doses D_atx and D_praz are in units of each drug's IC50 at its target (NET,
alpha1). A parameter dict without a dose key means that drug is absent
(dose 0); drugs are applied only through ``drugs()``.

The pre-registered readout is the steady-state tau burden T* reached from the
healthy start (E, Y, R, T) = (0, 0, 1, 0). The system is cooperative in the
order (E, Y, -R, T), so that trajectory converges monotonically to the fixed
point with the LARGEST R* (derivation §4); ``healthy_steady_state`` computes it
from the roots of the scalar function G(R) and is cross-checked against ODE
integration (unit tests, derivation check, and on a random subset in the run).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from src.dynamics.model import Drug, ODEModel

STATES = ("E", "Y", "R", "T")
PARAMETERS = ("kappa", "u0", "phi", "v", "m", "r_rep", "r_dam", "K", "q", "sigma", "eta", "eps", "hA", "c0")
HEALTHY_START = (0.0, 0.0, 1.0, 0.0)


# ---------------------------------------------------------------------------------------------- pieces
def _dose(p: dict, key: str) -> float:
    d = float(p.get(key, 0.0))      # absent key = drug not given
    if d < 0:
        raise ValueError(f"{key} = {d} < 0")
    return d


def derived(p: dict) -> dict:
    """Quantities that depend on the drugs: NET reuptake rate u, total
    clearance w, alpha1 transmission share s = S/E, clearance factor c."""
    missing = [k for k in PARAMETERS if k not in p]
    if missing:
        raise KeyError(f"missing parameters: {missing}")
    if not 0 < p["u0"] < 1:
        raise ValueError("u0 must lie in (0, 1): it is a share of baseline NE clearance")
    datx, dpraz = _dose(p, "D_atx"), _dose(p, "D_praz")
    u = p["u0"] / (1.0 + datx)
    w = u + (1.0 - p["u0"])
    s = 1.0 / (1.0 + dpraz)
    c = 1.0 - p["eps"] + p["eps"] * s
    return {"u": u, "w": w, "s": s, "c": c}


def hill(x, K: float, q: float):
    z = (np.maximum(x, 0.0) / K) ** q
    return z / (1.0 + z)


# ---------------------------------------------------------------------------------------------- ODE
def rhs(t: float, y: np.ndarray, p: dict) -> np.ndarray:
    E, Y, R, T = y
    d = derived(p)
    kre = p["kappa"] * max(R, 0.0) * max(E, 0.0)     # identical to kappa R E on the physical domain
    f = 1.0 / (1.0 + kre)
    x = p["m"] * Y
    return np.array([
        f - d["w"] * E,
        p["phi"] * d["u"] * E - (p["v"] + p["m"]) * Y,
        p["r_rep"] * (1.0 - R) - p["r_dam"] * (hill(x, p["K"], p["q"]) + p["hA"]) * R,
        p["c0"] * ((1.0 + p["sigma"] * x) * (1.0 + p["eta"] * d["s"] * E) * (1.0 + p["hA"]) - d["c"] * R * T),
    ])


def hill_prime(x, K: float, q: float):
    """dH/dx; for x > 0 equal to q H (1 - H) / x."""
    x = np.maximum(x, 0.0)
    z = (x / K) ** q
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(x > 0, q * z / (x * (1.0 + z) ** 2), 0.0 if q > 1 else (1.0 / K if q == 1 else np.inf))
    return out


def jac(t: float, y: np.ndarray, p: dict) -> np.ndarray:
    E, Y, R, T = y
    d = derived(p)
    kre = p["kappa"] * R * E
    g = 1.0 / (1.0 + kre) ** 2
    x = p["m"] * Y
    H, Hp = hill(x, p["K"], p["q"]), hill_prime(x, p["K"], p["q"])
    prod_seed = 1.0 + p["sigma"] * x
    prod_harm = 1.0 + p["eta"] * d["s"] * E
    return np.array([
        [-d["w"] - p["kappa"] * R * g, 0.0, -p["kappa"] * E * g, 0.0],
        [p["phi"] * d["u"], -(p["v"] + p["m"]), 0.0, 0.0],
        [0.0, -p["r_dam"] * R * p["m"] * Hp, -p["r_rep"] - p["r_dam"] * (H + p["hA"]), 0.0],
        [p["c0"] * prod_seed * p["eta"] * d["s"] * (1 + p["hA"]),
         p["c0"] * p["sigma"] * p["m"] * prod_harm * (1 + p["hA"]),
         -p["c0"] * d["c"] * T, -p["c0"] * d["c"] * R],
    ])


def build_model() -> ODEModel:
    return ODEModel("lc_rate_controller", STATES, rhs, nonnegative=STATES, conserved={}, jac=jac)


def initial_state(p: dict) -> np.ndarray:
    """The healthy start: no NE yet, full regulation, no cortical tau."""
    return np.array(HEALTHY_START, float)


# ---------------------------------------------------------------------------------------------- drugs
def _with(p: dict, key: str, dose: float) -> dict:
    if dose < 0:
        raise ValueError("dose must be >= 0")
    return dict(p, **{key: float(dose)})


def drugs() -> dict[str, Drug]:
    """Atomoxetine blocks NET (u -> u0/(1 + D)); prazosin blocks alpha1
    (S -> E/(1 + D), and the alpha1 share eps of clearance is reduced)."""
    return {"atomoxetine": Drug("atomoxetine", lambda p, d: _with(p, "D_atx", d)),
            "prazosin": Drug("prazosin", lambda p, d: _with(p, "D_praz", d))}


def outcome(traj) -> float:
    """Tau burden at the end of a trajectory (lower = better). The registered
    readout is its steady-state limit, ``healthy_steady_state(p)['T']``."""
    return float(traj.y[STATES.index("T"), -1])


# ---------------------------------------------------------------------------------------------- degenerate case
def coupling_params() -> list[str]:
    """Setting these to zero decouples the four equations (derivation §6)."""
    return ["kappa", "phi", "r_dam", "sigma", "eta"]


def uncoupled_solution(t, p: dict, y0) -> np.ndarray:
    """Closed form with kappa = phi = r_dam = sigma = eta = 0:
    E = 1/w + (E0 - 1/w) e^{-w t}; Y = Y0 e^{-(v+m) t}; R = 1 - (1 - R0) e^{-r_rep t};
    T solves T' = c0 [(1 + hA) - c R(t) T]: closed form when R0 = 1, otherwise
    the exact solution written as one quadrature."""
    for k in coupling_params():
        if p[k] != 0:
            raise ValueError(f"uncoupled_solution needs {k} = 0")
    t = np.asarray(t, float)
    d = derived(p)
    E0, Y0, R0, T0 = map(float, y0)
    w, vm, rr, c0, c, P = d["w"], p["v"] + p["m"], p["r_rep"], p["c0"], d["c"], 1.0 + p["hA"]
    E = 1.0 / w + (E0 - 1.0 / w) * np.exp(-w * t)
    Y = Y0 * np.exp(-vm * t)
    R = 1.0 - (1.0 - R0) * np.exp(-rr * t)
    k = c0 * c
    if R0 == 1.0:
        T = P / c + (T0 - P / c) * np.exp(-k * t)
    else:
        def intR(s):   # int_0^s R
            return s - (1.0 - R0) * (1.0 - np.exp(-rr * s)) / rr
        # T(t) = T0 e^{-k I(t)} + c0 P int_0^t e^{k (I(s) - I(t))} ds,  I(t) = int_0^t R
        T = np.array([T0 * np.exp(-k * intR(ti)) +
                      c0 * P * quad(lambda s, ti=ti: np.exp(k * (intR(s) - intR(ti))), 0.0, ti,
                                    epsabs=1e-13, epsrel=1e-12, limit=200)[0] for ti in t])
    return np.vstack([E, Y, R, T])


# ---------------------------------------------------------------------------------------------- fixed points
def E_star(R, p: dict):
    """Positive root of w E (1 + kappa R E) = 1 (rationalised: no cancellation)."""
    w = derived(p)["w"]
    R = np.asarray(R, float)
    return 2.0 / (w + np.sqrt(w * w + 4.0 * p["kappa"] * R * w))


def x_star(R, p: dict):
    """DOPEGAL load x = m Y* at the fixed point, as a function of R."""
    u = derived(p)["u"]
    return p["m"] * p["phi"] * u * E_star(R, p) / (p["v"] + p["m"])


def G(R, p: dict):
    """Fixed points <-> roots of G in (0, 1] (derivation §2)."""
    R = np.asarray(R, float)
    return R * (p["r_rep"] + p["r_dam"] * (hill(x_star(R, p), p["K"], p["q"]) + p["hA"])) - p["r_rep"]


def elasticity(R, p: dict):
    """eps_E = -R E*'(R)/E*(R) = kappa R E / (1 + 2 kappa R E), in [0, 1/2)."""
    kre = p["kappa"] * np.asarray(R, float) * E_star(R, p)
    return kre / (1.0 + 2.0 * kre)


def G_prime(R, p: dict):
    """G'(R) = r_rep + r_dam [H + hA - eps_E q H (1 - H)] (derivation §3)."""
    H = hill(x_star(R, p), p["K"], p["q"])
    return p["r_rep"] + p["r_dam"] * (H + p["hA"] - elasticity(R, p) * p["q"] * H * (1.0 - H))


def phi_curve(R, p: dict):
    """Phi(R) = (r_rep/r_dam)(1 - R)/R - H(x*(R)); fixed points at hA are the
    solutions of Phi(R) = hA (derivation §5)."""
    R = np.asarray(R, float)
    return (p["r_rep"] / p["r_dam"]) * (1.0 - R) / R - hill(x_star(R, p), p["K"], p["q"])


def phi_prime(R, p: dict):
    """Phi'(R) = [ -rho/R + eps_E q H (1 - H) ] / R, rho = r_rep / r_dam."""
    R = np.asarray(R, float)
    H = hill(x_star(R, p), p["K"], p["q"])
    return (-(p["r_rep"] / p["r_dam"]) / R + elasticity(R, p) * p["q"] * H * (1.0 - H)) / R


def r_grid(log_min: float = -9.0, log_switch: float = -2.0, n_log: int = 400, n_lin: int = 3600) -> np.ndarray:
    return np.unique(np.concatenate([np.logspace(log_min, log_switch, n_log),
                                     np.linspace(10.0 ** log_switch, 1.0, n_lin)]))


def _refine_extrema(fprime, grid: np.ndarray) -> np.ndarray:
    """Insert the zeros of fprime (located by sign changes on the grid) so that
    the function is monotone between consecutive points of the returned grid."""
    d = fprime(grid)
    pts = [brentq(fprime, grid[i], grid[i + 1], xtol=1e-15, rtol=1e-14)
           for i in range(len(grid) - 1) if d[i] * d[i + 1] < 0]
    return np.unique(np.concatenate([grid, pts]))


def roots_of_G(p: dict, grid: np.ndarray | None = None) -> np.ndarray:
    """All roots of G in (0, 1], ascending. G(0+) = -r_rep < 0 <= G(1)."""
    g = _refine_extrema(lambda r: G_prime(r, p), r_grid() if grid is None else grid)
    v = G(g, p)
    roots = [float(r) for r, val in zip(g, v) if val == 0.0]
    roots += [brentq(lambda r: float(G(r, p)), g[i], g[i + 1], xtol=1e-15, rtol=1e-14)
              for i in range(len(g) - 1) if v[i] * v[i + 1] < 0]
    roots = np.unique(np.round(roots, 14))
    if len(roots) == 0:
        raise RuntimeError(f"no root of G found in (0, 1] (impossible by derivation §2): {p}")
    return roots


def fixed_point_at(R: float, p: dict) -> np.ndarray:
    d = derived(p)
    E = float(E_star(R, p))
    Y = p["phi"] * d["u"] * E / (p["v"] + p["m"])
    x = p["m"] * Y
    T = (1.0 + p["sigma"] * x) * (1.0 + p["eta"] * d["s"] * E) * (1.0 + p["hA"]) / (d["c"] * R)
    return np.array([E, Y, R, T])


def fixed_points(p: dict) -> list[dict]:
    """Every fixed point with its analytic stability: stable iff G'(R*) > 0
    (derivation §3; the T direction always has eigenvalue -c0 c R* < 0)."""
    out = []
    for r in roots_of_G(p):
        gp = float(G_prime(r, p))
        out.append({"y": fixed_point_at(r, p), "G_prime": gp, "stable": gp > 0})
    return out


def healthy_steady_state(p: dict) -> dict:
    """The fixed point reached from the healthy start: the largest root of G
    (derivation §4)."""
    fps = fixed_points(p)
    top = fps[-1]
    return {**dict(zip(STATES, top["y"])), "stable": top["stable"], "n_fixed": len(fps),
            "n_stable": sum(f["stable"] for f in fps)}


def T_star(p: dict) -> float:
    return float(healthy_steady_state(p)["T"])


# ---------------------------------------------------------------------------------------------- readouts
def dose_response(p: dict, drug: str, doses) -> np.ndarray:
    dr = drugs()[drug]
    return np.array([T_star(dr.apply(p, float(d))) for d in doses])


def bistable_window(p: dict, lo: float = 0.0, hi: float = 1.0, grid: np.ndarray | None = None) -> dict:
    """The set of second-hit values hA in [lo, hi] with >= 2 stable fixed
    points (derivation §5). Fixed points at hA are the crossings of Phi(R) = hA
    on (0, 1]. Phi is monotone between consecutive points of the refined grid,
    so the number of crossings is constant between the breakpoints {Phi at its
    extrema, Phi(1)}; it is counted exactly at the midpoint of each interval
    between breakpoints. Three or more crossings = two or more stable states
    (crossings alternate stable / saddle, starting stable at the largest R).
    Returns the windows (intersected with [lo, hi]) and their total length."""
    g = _refine_extrema(lambda r: phi_prime(r, p), r_grid() if grid is None else grid)
    ph = phi_curve(g, p)
    interior = np.where(np.sign(np.diff(ph[:-1])) != np.sign(np.diff(ph[1:])))[0] + 1
    brk = np.unique(np.concatenate([[lo, hi], ph[interior], [ph[-1]]]))
    brk = brk[(brk >= lo) & (brk <= hi)]
    windows = []
    for a, b in zip(brk[:-1], brk[1:]):
        h = 0.5 * (a + b)
        n_cross = int(np.sum((ph[:-1] - h) * (ph[1:] - h) < 0))
        if n_cross >= 3:
            if windows and windows[-1][1] == a:
                windows[-1] = (windows[-1][0], float(b))
            else:
                windows.append((float(a), float(b)))
    return {"windows": windows, "length": float(sum(b - a for a, b in windows)), "bistable": bool(windows)}


def is_inverted_u(outcome, *, min_relative_rise: float) -> bool:
    """Mirror image of the U-shape test: an intermediate dose is WORST. The
    maximum is at an interior dose and exceeds the outcome at both the lowest
    and the highest dose by more than `min_relative_rise` of its value."""
    y = np.asarray(outcome, float)
    if len(y) < 3 or not np.all(np.isfinite(y)):
        raise ValueError("dose-response needs >= 3 finite points")
    i = int(np.argmax(y))
    if i in (0, len(y) - 1):
        return False
    top = abs(y[i]) if y[i] != 0 else 1.0
    return bool((y[i] - y[0]) > min_relative_rise * top and (y[i] - y[-1]) > min_relative_rise * top)


def readouts(p: dict, *, reference_dose: float, protective_ratio: float, harmful_ratio: float,
             dose_grid=None, min_relative_rise: float | None = None, bifurcation_range=(0.0, 1.0)) -> dict:
    """All pre-registered readouts for one parameter set (prereg §4)."""
    from src.dynamics.sensitivity import is_u_shaped

    base = healthy_steady_state(p)
    T0 = float(base["T"])
    dr = drugs()
    r_atx = T_star(dr["atomoxetine"].apply(p, reference_dose)) / T0
    r_praz = T_star(dr["prazosin"].apply(p, reference_dose)) / T0
    out = {"T0": T0, "R0": float(base["R"]), "ratio_atomoxetine": r_atx, "ratio_prazosin": r_praz,
           "protective_atomoxetine": r_atx < protective_ratio, "protective_prazosin": r_praz < protective_ratio,
           "harmful_atomoxetine": r_atx > harmful_ratio, "harmful_prazosin": r_praz > harmful_ratio}
    out["both_protective"] = out["protective_atomoxetine"] and out["protective_prazosin"]
    if dose_grid is not None:
        for name in ("atomoxetine", "prazosin"):
            y = dose_response(p, name, dose_grid)
            out[f"u_shape_{name}"] = is_u_shaped(np.asarray(dose_grid), y, min_relative_rise=min_relative_rise)
            out[f"inverted_u_{name}"] = is_inverted_u(y, min_relative_rise=min_relative_rise)
    bw = bistable_window(p, *bifurcation_range)
    out["bistable"] = bw["bistable"]
    out["bistable_window_length"] = bw["length"]
    return out
