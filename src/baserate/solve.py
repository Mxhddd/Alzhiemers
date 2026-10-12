"""Stage probabilities p(t) for a StageModel.

``stage_probabilities`` — the analysis method from the brief: piecewise
matrix exponential, p(t + h) = p(t) expm(Q(t + h/2) h), with step h supplied
by the researcher (Q(t1) and Q(t2) do not commute, so this is an
approximation whose error the researcher bounds by hand).

``stage_probabilities_ode`` — a high-accuracy ODE reference (Radau,
rtol 1e-10) used only to measure the step error of the method above.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from src.baserate.model import StageModel, build_generator

SUM_TOL = 1e-9
NEG_TOL = 1e-12


def _observed(model: StageModel, p_hidden: np.ndarray) -> np.ndarray:
    out = np.zeros((p_hidden.shape[0], model.spec.K))
    for h, k in enumerate(model.observed_of_hidden):
        out[:, k] += p_hidden[:, h]
    return out


def _check(p: np.ndarray, ages: np.ndarray) -> np.ndarray:
    if (p < -NEG_TOL).any():
        raise AssertionError(f"negative stage probability {p.min():.3g}")
    s = p.sum(axis=1)
    bad = np.abs(s - 1) > SUM_TOL
    if bad.any():
        raise AssertionError(f"stage probabilities do not sum to 1 at ages {ages[bad][:5]}")
    return np.clip(p, 0.0, None)  # only values within NEG_TOL of zero are touched


def stage_probabilities(model: StageModel, theta: np.ndarray, ages, *, step: float,
                        start_age: float) -> np.ndarray:
    """Observed-stage probabilities at each age (rows follow the input order)."""
    if step is None or not step > 0:
        raise ValueError("step must be a positive number of years (researcher-supplied)")
    ages = np.asarray(ages, dtype=float)
    if (ages < start_age).any():
        raise ValueError("ages before start_age")
    order = np.argsort(ages)
    p = model.initial()
    t = float(start_age)
    out = np.empty((len(ages), model.n_hidden))
    for idx in order:
        target = ages[idx]
        while t < target - 1e-12:
            h = min(step, target - t)
            Q = build_generator(model, t + h / 2, theta)
            p = p @ expm(Q * h)
            t += h
        out[idx] = p
    return _check(_observed(model, out), ages)


def stage_probabilities_ode(model: StageModel, theta: np.ndarray, ages, *, start_age: float) -> np.ndarray:
    ages = np.asarray(ages, dtype=float)
    sol = solve_ivp(lambda t, p: p @ build_generator(model, t, theta), (start_age, float(ages.max())),
                    model.initial(), method="Radau", t_eval=np.sort(np.unique(ages)),
                    rtol=1e-10, atol=1e-13)
    if not sol.success:
        raise RuntimeError(f"ODE reference failed: {sol.message}")
    lookup = {t: sol.y[:, i] for i, t in enumerate(sol.t)}
    hidden = np.array([lookup[a] for a in ages])
    return _check(_observed(model, hidden), ages)


def max_step_error(model: StageModel, theta: np.ndarray, ages, *, step: float, start_age: float) -> float:
    """Largest absolute difference between the stepped and reference solutions —
    logged with every fit so the researcher's step choice can be audited."""
    a = stage_probabilities(model, theta, ages, step=step, start_age=start_age)
    b = stage_probabilities_ode(model, theta, ages, start_age=start_age)
    return float(np.abs(a - b).max())
