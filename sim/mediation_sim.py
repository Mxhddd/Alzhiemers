"""SIMULATED cohort for Module 3 validation. Visible generator (ground rule 2);
never read by an analysis pipeline; nothing written to disk.

DAG used: X -> T, X -> M, X -> Y, T -> M -> Y, T -> Y; optional correlated
errors between M and Y (corr = rho_true) to emulate unmeasured M-Y
confounding, so the sensitivity analysis can be checked against the truth.
"""

from __future__ import annotations

import numpy as np


def simulate_cohort(n: int, *, beta2: float, gamma: float, beta3: float, rho_true: float = 0.0,
                    binary_exposure: bool = True, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, 2))
    p = 1 / (1 + np.exp(-(0.3 * X[:, 0] - 0.2 * X[:, 1])))
    T = rng.binomial(1, p).astype(float) if binary_exposure else 0.5 * X[:, 0] + rng.normal(size=n)
    cov = np.array([[1.0, rho_true], [rho_true, 1.0]])
    e = rng.multivariate_normal([0, 0], cov, size=n)
    M = 0.2 + beta2 * T + X @ np.array([0.3, -0.1]) + e[:, 0]
    Y = 1.0 + beta3 * T + gamma * M + X @ np.array([-0.2, 0.4]) + e[:, 1]
    return {"T": T, "M": M, "Y": Y, "X": X, "true_acme": beta2 * gamma, "true_ade": beta3}
