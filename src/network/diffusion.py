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


# --- Vectorised fitting (used by the analysis; equal to best_fit/crossover_tau,
# which the validation exercises, and checked against them in the unit tests).

def prediction_curves(diff: Diffuser, x0: np.ndarray, eval_idx: np.ndarray, tau_grid: np.ndarray) -> np.ndarray:
    """(n_tau, n_eval): the predicted pattern over eval_idx at every tau."""
    coef = diff.eigvecs.T @ x0
    decay = np.exp(-np.outer(tau_grid, diff.eigvals))
    return (decay * coef) @ diff.eigvecs[eval_idx].T


def _ranks(a: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata
    return rankdata(a, axis=-1)


def standardise_rows(P: np.ndarray, method: str) -> np.ndarray:
    """Rows centred and scaled so that a row-by-row dot product divided by the
    row length is Pearson r. method='spearman' ranks first. A constant row
    becomes NaN (its correlation is undefined)."""
    if method not in ("pearson", "spearman"):
        raise ValueError("correlation method must be 'pearson' or 'spearman'")
    P = np.atleast_2d(np.asarray(P, float))
    if method == "spearman":
        P = _ranks(P)
    C = P - P.mean(axis=1, keepdims=True)
    sd = np.sqrt((C ** 2).mean(axis=1, keepdims=True))
    with np.errstate(invalid="ignore", divide="ignore"):
        Z = np.where(sd > 1e-12 * (np.abs(P).max(axis=1, keepdims=True) + 1e-300), C / sd, np.nan)
    return Z


def corr_curves(Zpred: np.ndarray, Zmaps: np.ndarray) -> np.ndarray:
    """(n_maps, n_tau) correlations between standardised maps and standardised
    prediction curves (both from standardise_rows with the same method)."""
    return (np.atleast_2d(Zmaps) @ Zpred.T) / Zpred.shape[1]


def max_over_tau(R: np.ndarray, tau_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Best r and the tau where it occurs, per row of R (NaN columns ignored)."""
    if np.all(np.isnan(R)):
        raise ValueError("prediction is constant over the evaluated regions at every tau")
    filled = np.where(np.isnan(R), -np.inf, R)
    k = np.argmax(filled, axis=1)
    return filled[np.arange(len(R)), k], np.asarray(tau_grid)[k]


def crossover_from_curves(Pa: np.ndarray, Pb: np.ndarray, tau_grid: np.ndarray, threshold: float) -> float | None:
    """First tau at which two seeds' predicted patterns (rows of Pa, Pb)
    correlate (Pearson) at or above threshold; None if never in the grid."""
    Za, Zb = standardise_rows(Pa, "pearson"), standardise_rows(Pb, "pearson")
    r = (Za * Zb).mean(axis=1)
    hit = np.where(np.nan_to_num(r, nan=-np.inf) >= threshold)[0]
    return float(tau_grid[hit[0]]) if len(hit) else None
