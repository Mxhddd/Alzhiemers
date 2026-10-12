"""Fixed points (roots of rhs) and their stability from Jacobian eigenvalues."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import root

from src.dynamics.model import ODEModel


@dataclass
class FixedPoint:
    y: np.ndarray
    eigenvalues: np.ndarray
    stable: bool
    kind: str


def classify(eig: np.ndarray, tol: float = 1e-9) -> tuple[bool, str]:
    re = eig.real
    if np.any(np.abs(re) <= tol):
        return False, "non-hyperbolic (eigenvalue with ~zero real part: linearisation inconclusive)"
    if np.all(re < 0):
        return True, "stable focus" if np.any(np.abs(eig.imag) > tol) else "stable node"
    if np.all(re > 0):
        return False, "unstable focus" if np.any(np.abs(eig.imag) > tol) else "unstable node"
    return False, "saddle"


def find_fixed_points(model: ODEModel, p: dict, guesses: np.ndarray, *, tol: float = 1e-10,
                      dedupe_tol: float = 1e-6, require_nonnegative: bool = True) -> list[FixedPoint]:
    """Newton-type root finding from every guess; duplicates merged. Roots
    with negative entries in non-negative states are discarded only when
    `require_nonnegative` (they are unphysical), and the count is returned in
    the log by the caller via len(guesses) vs results."""
    found: list[np.ndarray] = []
    for g in np.atleast_2d(guesses):
        sol = root(lambda y: model.f(0.0, y, p), g, method="hybr", tol=tol)
        if not sol.success or np.max(np.abs(model.f(0.0, sol.x, p))) > 1e-8:
            continue
        if require_nonnegative:
            idx = [model.state_names.index(s) for s in model.nonnegative]
            if np.any(sol.x[idx] < -1e-9):
                continue
        if not any(np.max(np.abs(sol.x - f)) < dedupe_tol * max(1.0, np.max(np.abs(f))) for f in found):
            found.append(sol.x)
    out = []
    for y in sorted(found, key=lambda v: tuple(v)):
        eig = np.linalg.eigvals(model.jacobian(y, p))
        stable, kind = classify(eig)
        out.append(FixedPoint(y=y, eigenvalues=eig, stable=stable, kind=kind))
    return out


def grid_guesses(bounds: dict[str, tuple[float, float]], state_names, n_per_dim: int) -> np.ndarray:
    axes = [np.linspace(*bounds[s], n_per_dim) for s in state_names]
    return np.array(np.meshgrid(*axes, indexing="ij")).reshape(len(state_names), -1).T
