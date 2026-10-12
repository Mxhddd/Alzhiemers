"""SIMULATED geometry and maps for Module 6 validation (visible generator;
never read by an analysis pipeline; nothing written to disk).

* ``fibonacci_sphere`` — near-uniform parcel centroids on a unit sphere
  (one sphere per hemisphere, as in FreeSurfer's sphere.reg).
* ``smooth_field`` — a spatially autocorrelated random map: a sum of random
  von Mises–Fisher bumps. Two calls with different seeds give INDEPENDENT
  maps that nonetheless correlate strongly by chance — the known null.
* ``random_connectome`` — a symmetric, connected, distance-dependent graph.
"""

from __future__ import annotations

import numpy as np


def fibonacci_sphere(n: int) -> np.ndarray:
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])


def smooth_field(coords: np.ndarray, *, n_bumps: int, kappa: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    centres = rng.normal(size=(n_bumps, 3))
    centres /= np.linalg.norm(centres, axis=1, keepdims=True)
    w = rng.normal(size=n_bumps)
    return np.exp(kappa * (coords @ centres.T - 1)) @ w


def random_connectome(coords: np.ndarray, *, scale: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    W = np.exp(-d / scale) * rng.uniform(0.5, 1.5, size=d.shape)
    W = (W + W.T) / 2
    np.fill_diagonal(W, 0)
    return W
