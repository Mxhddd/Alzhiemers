"""Stage probabilities p(t) for a StageModel.

``stage_probabilities`` — the analysis method (brief §5): piecewise matrix
exponential, p(t + h) = p(t) expm(Q(t + h/2) h) — the exponential midpoint
rule. Its error is derived in manuscript/derivations/module2_step_size.md.

The integration grid (``integration_grid``) is shared with the JAX
likelihood: between consecutive knots — the start age, every evaluation age
and every breakpoint of A(t) — the interval is cut into equal steps no
longer than ``step``. No step straddles a discontinuity of A(t) (a jump
inside a step would make the local error O(h) instead of O(h^3)).

``stage_probabilities_ode`` — a high-accuracy ODE reference (Radau,
rtol 1e-10), used only to audit the step error at fitted parameters.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from src.baserate.model import StageModel, build_generator

SUM_TOL = 1e-9
NEG_TOL = 1e-12


@dataclass(frozen=True)
class Grid:
    t_mid: np.ndarray       # midpoint of each step
    h: np.ndarray           # length of each step
    record_after: np.ndarray  # for each requested age (input order): number of steps completed


def integration_grid(ages, *, step: float, start_age: float, breakpoints=()) -> Grid:
    if step is None or not step > 0:
        raise ValueError("step must be a positive number of years")
    ages = np.asarray(ages, dtype=float)
    if (ages < start_age).any():
        raise ValueError("ages before start_age")
    end = float(ages.max()) if len(ages) else start_age
    knots = {float(start_age), *map(float, ages)}
    knots |= {b for b in breakpoints if start_age < b < end}
    knots = sorted(knots)
    t_mid, h, ends = [], [], []
    for a, b in zip(knots, knots[1:]):
        n = int(np.ceil((b - a) / step - 1e-12))
        edges = np.linspace(a, b, n + 1)
        t_mid += list((edges[:-1] + edges[1:]) / 2)
        h += list(np.diff(edges))
        ends.append((b, len(h)))
    after = {start_age: 0, **{b: k for b, k in ends}}
    return Grid(np.array(t_mid), np.array(h), np.array([after[float(a)] for a in ages], dtype=int))


def aggregation_matrix(model: StageModel) -> np.ndarray:
    """Hidden-state probabilities -> observed-stage probabilities."""
    M = np.zeros((model.n_hidden, model.spec.K))
    M[np.arange(model.n_hidden), model.observed_of_hidden] = 1.0
    return M


def _check(p: np.ndarray, ages: np.ndarray) -> np.ndarray:
    if (p < -NEG_TOL).any():
        raise AssertionError(f"negative stage probability {p.min():.3g}")
    s = p.sum(axis=1)
    bad = np.abs(s - 1) > SUM_TOL
    if bad.any():
        raise AssertionError(f"stage probabilities do not sum to 1 at ages {ages[bad][:5]}")
    return np.clip(p, 0.0, None)  # only values within NEG_TOL of zero are touched


def hidden_probabilities(model: StageModel, theta: np.ndarray, grid: Grid) -> np.ndarray:
    """Hidden-state probabilities after each step (row k = after k+1 steps)."""
    p = model.initial()
    out = np.empty((len(grid.h), model.n_hidden))
    for k, (tm, h) in enumerate(zip(grid.t_mid, grid.h)):
        p = p @ expm(build_generator(model, tm, theta) * h)
        out[k] = p
    return out


def stage_probabilities(model: StageModel, theta: np.ndarray, ages, *, step: float,
                        start_age: float) -> np.ndarray:
    """Observed-stage probabilities at each age (rows follow the input order)."""
    ages = np.asarray(ages, dtype=float)
    grid = integration_grid(ages, step=step, start_age=start_age, breakpoints=model.breakpoints)
    traj = hidden_probabilities(model, theta, grid)
    p0 = model.initial()
    hidden = np.array([traj[k - 1] if k > 0 else p0 for k in grid.record_after])
    return _check(hidden @ aggregation_matrix(model), ages)


def stage_probabilities_ode(model: StageModel, theta: np.ndarray, ages, *, start_age: float) -> np.ndarray:
    """Radau reference, integrated piecewise between A(t) breakpoints so the
    solver never steps across a discontinuity."""
    ages = np.asarray(ages, dtype=float)
    targets = np.sort(np.unique(ages))
    knots = sorted({start_age, *targets, *[b for b in model.breakpoints if start_age < b < targets.max()]})
    p = model.initial()
    lookup = {start_age: p}
    for a, b in zip(knots, knots[1:]):
        sol = solve_ivp(lambda t, y: y @ build_generator(model, t, theta), (a, b), p, method="Radau",
                        rtol=1e-10, atol=1e-13)
        if not sol.success:
            raise RuntimeError(f"ODE reference failed: {sol.message}")
        p = sol.y[:, -1]
        lookup[b] = p
    hidden = np.array([lookup[a] for a in ages])
    return _check(hidden @ aggregation_matrix(model), ages)


def max_step_error(model: StageModel, theta: np.ndarray, ages, *, step: float, start_age: float) -> float:
    """Largest absolute difference between the stepped and reference solutions —
    logged with every fit so the step choice is audited at the fitted parameters."""
    a = stage_probabilities(model, theta, ages, step=step, start_age=start_age)
    b = stage_probabilities_ode(model, theta, ages, start_age=start_age)
    return float(np.abs(a - b).max())


def error_bound(model: StageModel, theta: np.ndarray, *, start_age: float, end_age: float,
                dt: float = 0.25, d: float = 1e-3) -> float:
    """Phi(theta) = integral of ||Q''/24 + [Q, Q']/12||_inf over [start, end]
    (manuscript/derivations/module2_step_size.md §3): the global error of the
    exponential-midpoint rule is <= h^2 Phi. Q' and Q'' by central differences
    within amyloid bands (points within 2d of a breakpoint are skipped: the
    grid puts a knot there, so no step crosses it)."""
    bps = np.array(model.breakpoints, dtype=float)
    total = 0.0
    for t in np.arange(start_age + dt / 2, end_age, dt):
        if len(bps) and np.min(np.abs(bps - t)) < 2 * d:
            continue
        Q0 = build_generator(model, t, theta)
        Qp, Qm = build_generator(model, t + d, theta), build_generator(model, t - d, theta)
        Q1, Q2 = (Qp - Qm) / (2 * d), (Qp - 2 * Q0 + Qm) / d**2
        E = Q2 / 24 + (Q0 @ Q1 - Q1 @ Q0) / 12
        total += float(np.abs(E).sum(axis=1).max()) * dt
    return total
