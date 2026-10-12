"""Causal mediation for Module 3 (Killer 2): does LC integrity (M) mediate
the effect of a risk factor (T) on tau (Y)?

Linear-linear models without exposure-mediator interaction, covariates X:
    M = a2 + beta2 T + X xi2 + e2                    (mediator model)
    Y = a3 + beta3 T + gamma M + X xi3 + e3          (outcome model)
    Y = a1 + beta1 T + X xi1 + e1                    (total-effect model)
ACME = beta2 * gamma, ADE = beta3, total = ACME + ADE (= beta1 by OLS algebra),
proportion mediated = ACME / total.

Sensitivity to unmeasured mediator-outcome confounding (Imai, Keele &
Yamamoto 2010, Stat Sci 25:51), rho = corr(e2, e3):
    ACME(rho) = (beta2 sigma1 / sigma2) * (rho~ - rho sqrt((1 - rho~^2) / (1 - rho^2)))
with sigma1, sigma2 the residual SDs of the total-effect and mediator models
and rho~ = corr of their residuals. ACME(0) is the OLS product of
coefficients exactly, and ACME(rho~) = 0.

These estimates are causal only under sequential ignorability (no unmeasured
T-M, T-Y or M-Y confounding given X, and no M-Y confounder affected by T).
Observational data cannot establish that; the rho analysis shows how much
violation the result survives.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _ols(y: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return coef, y - X @ coef


@dataclass
class PointEstimates:
    acme: float
    ade: float
    total: float
    prop_mediated: float
    beta2: float
    gamma: float
    sigma1: float
    sigma2: float
    rho_tilde: float

    def acme_at(self, rho: np.ndarray) -> np.ndarray:
        rho = np.asarray(rho, float)
        if np.any(np.abs(rho) >= 1):
            raise ValueError("rho must lie strictly inside (-1, 1)")
        rt = self.rho_tilde
        return self.beta2 * self.sigma1 / self.sigma2 * (rt - rho * np.sqrt((1 - rt**2) / (1 - rho**2)))


def estimate(T: np.ndarray, M: np.ndarray, Y: np.ndarray, X: np.ndarray | None) -> PointEstimates:
    n = len(T)
    for name, v in (("T", T), ("M", M), ("Y", Y)):
        if len(v) != n or not np.all(np.isfinite(v)):
            raise ValueError(f"{name} has the wrong length or non-finite values")
    Z = np.ones((n, 1)) if X is None else np.column_stack([np.ones(n), X])
    if X is not None and not np.all(np.isfinite(X)):
        raise ValueError("covariates contain non-finite values")
    if np.std(T) == 0 or np.std(M) == 0:
        raise ValueError("exposure or mediator has no variation")
    c2, e2 = _ols(M, np.column_stack([Z, T]))
    c3, _ = _ols(Y, np.column_stack([Z, T, M]))
    c1, e1 = _ols(Y, np.column_stack([Z, T]))
    beta2, beta3, gamma = c2[-1], c3[-2], c3[-1]
    acme = beta2 * gamma
    total = acme + beta3
    s1, s2 = np.sqrt(np.mean(e1**2)), np.sqrt(np.mean(e2**2))
    rt = float(np.mean(e1 * e2) / (s1 * s2))
    return PointEstimates(acme=float(acme), ade=float(beta3), total=float(total),
                          prop_mediated=float(acme / total) if total != 0 else np.nan,
                          beta2=float(beta2), gamma=float(gamma), sigma1=float(s1), sigma2=float(s2),
                          rho_tilde=rt)


def bootstrap(T, M, Y, X, *, n_boot: int, seed: int) -> list[PointEstimates]:
    rng = np.random.default_rng(seed)
    n = len(T)
    out = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        out.append(estimate(T[i], M[i], Y[i], None if X is None else X[i]))
    return out


def percentile_ci(values: np.ndarray, level: float) -> tuple[float, float]:
    a = (1 - level) / 2
    lo, hi = np.quantile(np.asarray(values, float), [a, 1 - a])
    return float(lo), float(hi)


def summarise(point: PointEstimates, boots: list[PointEstimates], *, level: float, rho_grid: np.ndarray) -> dict:
    out = {"estimates": {}, "sensitivity": {}}
    for k in ("acme", "ade", "total", "prop_mediated"):
        vals = np.array([getattr(b, k) for b in boots])
        lo, hi = percentile_ci(vals[np.isfinite(vals)], level)
        out["estimates"][k] = {"estimate": getattr(point, k), "ci_lower": lo, "ci_upper": hi}
    curves = np.array([b.acme_at(rho_grid) for b in boots])
    a = (1 - level) / 2
    lo, hi = np.quantile(curves, [a, 1 - a], axis=0)
    est = point.acme_at(rho_grid)
    includes_zero = (lo <= 0) & (hi >= 0)
    out["sensitivity"] = {
        "rho": rho_grid, "acme": est, "ci_lower": lo, "ci_upper": hi,
        "rho_at_zero_acme": point.rho_tilde,
        # smallest |rho| at which the ACME interval includes 0 (None if never in the grid)
        "rho_ci_includes_zero": (float(rho_grid[includes_zero][np.argmin(np.abs(rho_grid[includes_zero]))])
                                 if includes_zero.any() else None),
    }
    return out


def killer2_verdict(summary: dict, predicted_sign: int) -> str:
    e = summary["estimates"]["acme"]
    excludes_zero = e["ci_lower"] > 0 or e["ci_upper"] < 0
    if excludes_zero and np.sign(e["estimate"]) == predicted_sign:
        return "ACME in the predicted direction, CI excludes 0 — Killer 2 NOT met for this exposure"
    if excludes_zero:
        return "ACME OPPOSITE to prediction, CI excludes 0 — Killer 2 MET (hypothesis fails)"
    return "ACME CI includes 0 — Killer 2 MET (no evidence of mediation by LC integrity)"
