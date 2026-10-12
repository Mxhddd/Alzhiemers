"""Network diffusion of tau on a structural connectome (Raj et al. 2012).

    dx/dt = -beta L x + alpha x
solved through the eigendecomposition L = V diag(lambda) V^T:
    x(t) = e^{alpha t} V diag(e^{-beta lambda t}) V^T x0.
beta and t enter only as tau = beta t, which is what is fitted. Uniform
growth alpha rescales every region equally, so it cannot change a spatial
correlation; it is carried for completeness and has no effect on seed
comparisons (stated in the output).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def laplacian(W: np.ndarray, kind: str) -> np.ndarray:
    W = np.asarray(W, float)
    if W.ndim != 2 or W.shape[0] != W.shape[1]:
        raise ValueError("connectome must be a square matrix")
    if not np.allclose(W, W.T):
        raise ValueError("connectome must be symmetric")
    if np.any(W < 0) or np.any(np.diag(W) != 0):
        raise ValueError("connectome must be non-negative with a zero diagonal")
    d = W.sum(axis=1)
    if np.any(d == 0):
        raise ValueError(f"{int(np.sum(d == 0))} isolated region(s): diffusion is undefined there")
    if kind == "unnormalized":
        return np.diag(d) - W
    if kind == "normalized":
        s = 1 / np.sqrt(d)
        return np.eye(len(W)) - s[:, None] * W * s[None, :]
    raise ValueError("laplacian kind must be 'unnormalized' or 'normalized'")


@dataclass
class Diffuser:
    eigvals: np.ndarray
    eigvecs: np.ndarray

    @classmethod
    def from_laplacian(cls, L: np.ndarray) -> "Diffuser":
        lam, V = np.linalg.eigh(L)
        if lam.min() < -1e-9:
            raise ValueError("Laplacian has a negative eigenvalue: check the connectome")
        return cls(np.clip(lam, 0, None), V)

    def evolve(self, x0: np.ndarray, tau: float, alpha_t: float = 0.0) -> np.ndarray:
        V = self.eigvecs
        return np.exp(alpha_t) * (V @ (np.exp(-self.eigvals * tau) * (V.T @ x0)))


def seed_vector(n: int, seed_idx: list[int]) -> np.ndarray:
    x = np.zeros(n)
    x[seed_idx] = 1.0 / len(seed_idx)
    return x


def best_fit(diff: Diffuser, x0: np.ndarray, observed: np.ndarray, eval_idx: np.ndarray,
             tau_grid: np.ndarray) -> tuple[float, float, np.ndarray]:
    """Pearson r between prediction and observed over eval_idx, maximised over tau."""
    rs = np.empty(len(tau_grid))
    for k, tau in enumerate(tau_grid):
        x = diff.evolve(x0, tau)[eval_idx]
        rs[k] = np.corrcoef(x, observed)[0, 1] if np.std(x) > 0 else np.nan
    if np.all(np.isnan(rs)):
        raise ValueError("prediction is constant over the evaluated regions at every tau")
    k = int(np.nanargmax(rs))
    return float(rs[k]), float(tau_grid[k]), rs


def crossover_tau(diff: Diffuser, x0_a: np.ndarray, x0_b: np.ndarray, eval_idx: np.ndarray,
                  tau_grid: np.ndarray, threshold: float) -> float | None:
    """First tau at which the two seeds' predicted patterns correlate above
    `threshold` over eval_idx — after it, the data cannot tell the seeds apart
    (high-eigenvalue modes have decayed). None if never reached in the grid."""
    for tau in tau_grid:
        a, b = diff.evolve(x0_a, tau)[eval_idx], diff.evolve(x0_b, tau)[eval_idx]
        if np.std(a) > 0 and np.std(b) > 0 and np.corrcoef(a, b)[0, 1] >= threshold:
            return float(tau)
    return None
