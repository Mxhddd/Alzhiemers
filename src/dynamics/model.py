"""ODE model interface for Module 4.

A model is a pure right-hand side ``rhs(t, y, p) -> dy/dt`` plus metadata.
The machinery (simulate, fixed points, bifurcation, sensitivity) works on
any ``ODEModel``; it is validated on reference systems with analytic answers
(tests/validation/test_dynamics.py). The LC model itself is the researcher's
— see src/dynamics/lc_model.py and manuscript/notes/module4_model_proposal.md.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

Rhs = Callable[[float, np.ndarray, dict], np.ndarray]


@dataclass(frozen=True)
class Drug:
    """A drug as a dose-dependent change to named parameters:
    ``apply(p, dose) -> new parameter dict`` (pure)."""

    name: str
    apply: Callable[[dict, float], dict]


@dataclass(frozen=True)
class ODEModel:
    name: str
    state_names: tuple[str, ...]
    rhs: Rhs
    # States that must stay >= 0 (concentrations, densities).
    nonnegative: tuple[str, ...] = ()
    # Conserved quantities: name -> function(y) that must stay constant.
    conserved: dict[str, Callable[[np.ndarray], float]] = field(default_factory=dict)
    # Optional analytic Jacobian jac(t, y, p); finite differences otherwise.
    jac: Callable[[float, np.ndarray, dict], np.ndarray] | None = None

    def f(self, t: float, y: np.ndarray, p: dict) -> np.ndarray:
        out = np.asarray(self.rhs(t, np.asarray(y, float), p), dtype=float)
        if out.shape != (len(self.state_names),):
            raise ValueError(f"{self.name}: rhs returned shape {out.shape}, expected ({len(self.state_names)},)")
        return out

    def jacobian(self, y: np.ndarray, p: dict, t: float = 0.0, h: float = 1e-7) -> np.ndarray:
        if self.jac is not None:
            return np.asarray(self.jac(t, y, p), dtype=float)
        n = len(y)
        J = np.empty((n, n))
        for j in range(n):
            step = h * max(1.0, abs(y[j]))
            e = np.zeros(n); e[j] = step
            J[:, j] = (self.f(t, y + e, p) - self.f(t, y - e, p)) / (2 * step)
        return J
