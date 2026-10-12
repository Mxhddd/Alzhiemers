"""Integration with built-in checks: non-negativity, conservation, and
solver independence (two stiff solvers, tolerances tightened 10x)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from src.dynamics.model import ODEModel


@dataclass
class Trajectory:
    t: np.ndarray
    y: np.ndarray  # shape (n_states, n_times)
    method: str
    rtol: float
    atol: float


def simulate(model: ODEModel, p: dict, y0, t_span, t_eval, *, method: str, rtol: float, atol: float,
             neg_tol: float = 1e-9, conservation_tol: float = 1e-6) -> Trajectory:
    if method not in ("Radau", "BDF", "LSODA"):
        raise ValueError("use a stiff solver: Radau, BDF or LSODA")
    sol = solve_ivp(lambda t, y: model.f(t, y, p), t_span, np.asarray(y0, float), method=method,
                    t_eval=t_eval, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"{model.name} / {method}: {sol.message}")
    for name in model.nonnegative:
        i = model.state_names.index(name)
        if sol.y[i].min() < -neg_tol:
            raise AssertionError(f"{model.name}: state {name} went negative ({sol.y[i].min():.3g})")
    for name, g in model.conserved.items():
        vals = np.array([g(sol.y[:, k]) for k in range(sol.y.shape[1])])
        drift = np.abs(vals - vals[0]).max()
        if drift > conservation_tol * max(1.0, abs(vals[0])):
            raise AssertionError(f"{model.name}: conserved quantity {name} drifted by {drift:.3g}")
    return Trajectory(sol.t, sol.y, method, rtol, atol)


def solver_independence(model: ODEModel, p: dict, y0, t_span, t_eval, *, rtol: float, atol: float,
                        methods=("Radau", "BDF")) -> dict:
    """Two stiff solvers at the working tolerance and with both tolerances
    tightened 10x (brief §5). Reports the cross-solver max relative
    difference at each level. Results are numerics-independent when the
    tight-level difference is small AND it shrank on tightening (the
    solvers are converging to the same answer)."""
    m1, m2 = methods
    runs = {(m, k): simulate(model, p, y0, t_span, t_eval, method=m, rtol=rtol / k, atol=atol / k)
            for m in methods for k in (1, 10)}

    def rel(a, b):
        scale = np.maximum(np.abs(b), atol)
        return float(np.max(np.abs(a - b) / scale))

    working = rel(runs[(m2, 1)].y, runs[(m1, 1)].y)
    tight = rel(runs[(m2, 10)].y, runs[(m1, 10)].y)
    return {"cross_solver_working": working, "cross_solver_tight": tight,
            "converging": tight <= working, "methods": methods}
