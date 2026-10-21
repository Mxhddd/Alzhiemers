"""Group connectome assembly for Module 6 (prereg/network.md, section 4).

The per-subject matrices come from tractography run outside this repository
on HCP data (the steps and settings are fixed in config/network.yaml,
`connectome`). This file does only what follows: average across subjects,
optionally keep edges present in enough subjects, and build the
'uniform LC projection' variant used in a sensitivity analysis.
"""

from __future__ import annotations

import numpy as np


def _check_square_symmetric(W: np.ndarray, name: str) -> np.ndarray:
    W = np.asarray(W, float)
    if W.ndim != 2 or W.shape[0] != W.shape[1]:
        raise ValueError(f"{name}: must be a square matrix")
    if not np.allclose(W, W.T, rtol=1e-6, atol=1e-9):
        raise ValueError(f"{name}: must be symmetric")
    if np.any(W < 0):
        raise ValueError(f"{name}: negative weights")
    return W


def group_average(mats: list[np.ndarray], min_fraction: float | None) -> np.ndarray:
    """Mean edge weight across subjects (absent edges count as zero). With
    min_fraction, an edge present (non-zero) in fewer than that fraction of
    subjects is set to zero (consistency thresholding)."""
    if not mats:
        raise ValueError("no subject matrices")
    stack = np.stack([_check_square_symmetric(m, f"subject {i}") for i, m in enumerate(mats)])
    W = stack.mean(axis=0)
    if min_fraction is not None:
        if not 0 < min_fraction <= 1:
            raise ValueError("min_fraction must be in (0, 1]")
        present = (stack > 0).mean(axis=0)
        W = np.where(present >= min_fraction, W, 0.0)
    W = (W + W.T) / 2
    np.fill_diagonal(W, 0.0)
    return W


def uniform_projection(W: np.ndarray, node: int, targets: list[int]) -> np.ndarray:
    """Replace one node's connections by equal weights to `targets`, keeping
    its total connection weight. Every other edge is unchanged."""
    W = _check_square_symmetric(W, "connectome").copy()
    if node in targets:
        raise ValueError("a node cannot project to itself")
    if not targets:
        raise ValueError("no targets")
    total = W[node].sum()
    if total <= 0:
        raise ValueError("node has no connections to redistribute")
    W[node, :] = 0.0
    W[:, node] = 0.0
    W[node, targets] = total / len(targets)
    W[targets, node] = total / len(targets)
    return W
