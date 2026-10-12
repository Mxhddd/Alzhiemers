"""Spatial-permutation ("spin") nulls for brain maps (Alexander-Bloch et al.
2018, NeuroImage 178:540; parcel reassignment without duplicates as in Váša
et al. 2018, Cereb Cortex 28:281).

Brain maps are spatially autocorrelated, so two unrelated maps correlate far
above chance; shuffling regions ignores this and inflates false positives.
A spin rotates each hemisphere's spherical parcel centroids by a uniformly
random rotation (the right hemisphere gets the mirror-image rotation) and
reassigns parcels to the rotated positions, preserving spatial structure.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

MIRROR = np.diag([-1.0, 1.0, 1.0])


def random_rotation(rng: np.random.Generator) -> np.ndarray:
    """Uniform (Haar) random rotation in SO(3) via QR with sign correction."""
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def _check_sphere(c: np.ndarray, name: str) -> np.ndarray:
    c = np.asarray(c, float)
    if c.ndim != 2 or c.shape[1] != 3:
        raise ValueError(f"{name}: coordinates must be (n, 3)")
    norms = np.linalg.norm(c, axis=1)
    if np.any(norms == 0):
        raise ValueError(f"{name}: zero-length coordinate")
    return c / norms[:, None]


def spin_permutations(coords_lh: np.ndarray, coords_rh: np.ndarray, n_perm: int, *, seed: int,
                      method: str) -> np.ndarray:
    """(n_perm, n_lh + n_rh) index arrays: null map = original[perm]. Left
    parcels come first, right second, matching the concatenated map order."""
    lh, rh = _check_sphere(coords_lh, "left"), _check_sphere(coords_rh, "right")
    if method not in ("hungarian", "nearest"):
        raise ValueError("method must be 'hungarian' (no duplicates) or 'nearest' (Alexander-Bloch original)")
    rng = np.random.default_rng(seed)
    n_l = len(lh)
    out = np.empty((n_perm, n_l + len(rh)), dtype=int)
    for k in range(n_perm):
        R = random_rotation(rng)
        for coords, Rh, offset, sl in ((lh, R, 0, slice(0, n_l)), (rh, MIRROR @ R @ MIRROR, n_l, slice(n_l, None))):
            D = cdist(coords @ Rh.T, coords)  # rows: rotated parcels, cols: original positions
            if method == "hungarian":
                rows, cols = linear_sum_assignment(D)
                perm = np.empty(len(coords), dtype=int)
                perm[cols] = rows
            else:
                perm = np.argmin(D, axis=0)
            out[k, sl] = perm + offset
    return out


def naive_permutations(n: int, n_perm: int, *, seed: int) -> np.ndarray:
    """Plain shuffles — provided ONLY to demonstrate the inflated false-positive
    rate in validation. Analyses must not use them (brief §5)."""
    rng = np.random.default_rng(seed)
    return np.array([rng.permutation(n) for _ in range(n_perm)])


def perm_pvalue(observed: float, null: np.ndarray) -> float:
    """Two-sided, with the +1 correction so p is never 0."""
    null = np.asarray(null, float)
    return float((1 + np.sum(np.abs(null) >= abs(observed))) / (1 + len(null)))


def corr(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(x, y)[0, 1])
