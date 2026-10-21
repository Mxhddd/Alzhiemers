"""Cumulative-link (ordinal) models and ordinal-outcome mediation for Module 3.

Braak NFT stage (RADC `braaksc`, 0-VI) is an ordered outcome, so the
association models are cumulative-link regressions:

    P(Y <= k | x) = F(theta_k - x beta),   k = 0 .. J-2

with F the logistic CDF (proportional odds; exp(beta) is an odds ratio for a
higher stage) or the normal CDF (ordered probit; beta on a latent N(0, 1)
scale). Observations may be weighted (inverse-probability weights).

Mediation with an ordinal outcome (the conditional LC branch) uses an ordered
probit outcome model and a linear mediator model,

    M  = a2 + beta2 T + X xi2 + e2,          e2 ~ N(0, sigma2^2)
    Y* = beta3 T + gamma M + X xi3 + e3,     e3 ~ N(0, 1),  Y = k if theta_{k-1} < Y* <= theta_k

so that every potential outcome Y(t, M(t')) has a closed-form distribution
(Y* given X is normal) and the natural indirect / direct effects on the
expected-stage scale are computed exactly, without Monte Carlo error.

Sensitivity to unmeasured mediator-outcome confounding: rho = corr(e2, e3).
Writing e3 = (rho / sigma2) e2 + sqrt(1 - rho^2) u, the ordered probit fitted
to the observed data estimates (with s = sqrt(1 - rho^2))

    gamma_obs = (gamma + rho / sigma2) / s,     beta3_obs = (beta3 - rho beta2 / sigma2) / s,
    xi3_obs   = (xi3 - rho xi2 / sigma2) / s,   theta_obs = (theta + rho a2 / sigma2) / s,

so for any assumed rho the true outcome-model parameters are recovered by
inverting these relations; the indirect effect is zero exactly at
rho* = gamma_obs sigma2 / sqrt(1 + gamma_obs^2 sigma2^2). This is the
ordered-probit analogue of Imai, Keele & Yamamoto (2010, Stat Sci 25:51); it
is checked against simulated truth in tests/validation/test_mediation.py.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import optimize, special

LINKS = ("logit", "probit")


def _cdf_pdf(z: np.ndarray, link: str) -> tuple[np.ndarray, np.ndarray]:
    if link == "logit":
        F = special.expit(z)
        return F, F * (1 - F)
    if link == "probit":
        return special.ndtr(z), np.exp(-0.5 * z * z) / np.sqrt(2 * np.pi)
    raise ValueError(f"link must be one of {LINKS}")


def _interval_prob(upper: np.ndarray, lower: np.ndarray, link: str) -> np.ndarray:
    """F(upper) - F(lower), computed in the tail where it is accurate."""
    flip = lower > 0
    up = np.where(flip, -lower, upper)
    lo = np.where(flip, -upper, lower)
    Fu, _ = _cdf_pdf(up, link)
    Fl, _ = _cdf_pdf(lo, link)
    return np.maximum(Fu - Fl, 1e-300)


@dataclass
class OrdinalFit:
    link: str
    levels: np.ndarray            # observed outcome levels, ascending
    beta: np.ndarray              # slopes (no intercept)
    theta: np.ndarray             # J-1 increasing thresholds
    cov_beta: np.ndarray          # inverse observed information, slope block
    loglik: float
    n: int
    converged: bool
    names: list[str] = field(default_factory=list)

    def se(self) -> np.ndarray:
        return np.sqrt(np.diag(self.cov_beta))

    def coef(self, name: str) -> tuple[float, float]:
        i = self.names.index(name)
        return float(self.beta[i]), float(self.se()[i])

    def prob_ge(self, eta: np.ndarray, scale: float | np.ndarray = 1.0) -> np.ndarray:
        """P(Y >= level_j), j = 1..J-1, for linear predictors eta (n x J-1)."""
        F, _ = _cdf_pdf((np.asarray(eta)[:, None] - self.theta[None, :]) / np.asarray(scale).reshape(-1, 1), self.link)
        return F


def _unpack(par: np.ndarray, p: int) -> tuple[np.ndarray, np.ndarray]:
    beta = par[:p]
    raw = par[p:]
    theta = np.cumsum(np.concatenate([[raw[0]], np.exp(raw[1:])]))
    return beta, theta


def _negloglik_grad(par, X, yi, w, link, p):
    beta, theta = _unpack(par, p)
    eta = X @ beta
    J = len(theta) + 1
    th = np.concatenate([[-np.inf], theta, [np.inf]])
    upper = th[yi + 1] - eta
    lower = th[yi] - eta
    prob = _interval_prob(upper, lower, link)
    _, fu = _cdf_pdf(np.where(np.isfinite(upper), upper, 0.0), link)
    _, fl = _cdf_pdf(np.where(np.isfinite(lower), lower, 0.0), link)
    fu = np.where(np.isfinite(upper), fu, 0.0)
    fl = np.where(np.isfinite(lower), fl, 0.0)
    nll = -np.sum(w * np.log(prob))
    g_eta = w * (fu - fl) / prob                  # d nll / d eta
    g_beta = X.T @ g_eta
    g_theta = np.zeros(J - 1)
    up_idx = yi < J - 1                           # theta_{y} appears as upper bound
    np.add.at(g_theta, yi[up_idx], -(w * fu / prob)[up_idx])
    lo_idx = yi > 0                               # theta_{y-1} appears as lower bound
    np.add.at(g_theta, yi[lo_idx] - 1, (w * fl / prob)[lo_idx])
    raw = par[p:]
    g_raw = np.empty_like(raw)
    g_raw[0] = g_theta.sum()
    tail = np.cumsum(g_theta[::-1])[::-1]         # sum_{k>=m} g_theta_k
    g_raw[1:] = np.exp(raw[1:]) * tail[1:]
    return nll, np.concatenate([g_beta, g_raw])


def fit_ordinal(y, X, *, link: str = "logit", weights=None, names: list[str] | None = None,
                hessian: bool = True) -> OrdinalFit:
    """Maximum-likelihood cumulative-link model. y: ordered codes (any numeric
    levels); X: n x p covariates WITHOUT an intercept column."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    n, p = X.shape
    if len(y) != n or not (np.all(np.isfinite(y)) and np.all(np.isfinite(X))):
        raise ValueError("outcome or covariates have the wrong length or non-finite values")
    w = np.ones(n) if weights is None else np.asarray(weights, float)
    if np.any(w < 0) or not np.all(np.isfinite(w)):
        raise ValueError("weights must be finite and non-negative")
    levels = np.unique(y)
    if len(levels) < 2:
        raise ValueError("outcome has fewer than two observed levels")
    yi = np.searchsorted(levels, y)
    cum = np.cumsum(np.bincount(yi, weights=w, minlength=len(levels)))[:-1] / w.sum()
    cum = np.clip(cum, 1e-4, 1 - 1e-4)
    th0 = (special.logit(cum) if link == "logit" else special.ndtri(cum))
    th0 = np.maximum.accumulate(th0 + np.arange(len(th0)) * 1e-6)
    raw0 = np.concatenate([[th0[0]], np.log(np.maximum(np.diff(th0), 1e-4))])
    # centre columns for conditioning; slopes are unchanged, thresholds shift
    mu = np.average(X, axis=0, weights=w)
    Xc = X - mu
    par0 = np.concatenate([np.zeros(p), raw0])
    res = optimize.minimize(_negloglik_grad, par0, args=(Xc, yi, w, link, p), jac=True, method="BFGS",
                            options={"gtol": 1e-8, "maxiter": 2000})
    if not res.success:
        res2 = optimize.minimize(_negloglik_grad, res.x, args=(Xc, yi, w, link, p), jac=True,
                                 method="L-BFGS-B", options={"maxiter": 5000, "gtol": 1e-9})
        if res2.fun <= res.fun:
            res = res2
    beta, theta_c = _unpack(res.x, p)
    theta = theta_c + mu @ beta                   # back to uncentred X
    cov = np.full((p, p), np.nan)
    if hessian:
        H = _num_hessian(lambda q: _negloglik_grad(q, Xc, yi, w, link, p)[1], res.x)
        try:
            cov = np.linalg.inv(H)[:p, :p]
        except np.linalg.LinAlgError:
            pass
    gnorm = float(np.max(np.abs(res.jac))) if res.jac is not None else np.inf
    return OrdinalFit(link=link, levels=levels, beta=beta, theta=theta, cov_beta=cov, loglik=-float(res.fun),
                      n=n, converged=bool(res.success or gnorm < 1e-4), names=list(names or [f"x{i}" for i in range(p)]))


def _num_hessian(grad, x, eps=1e-5):
    k = len(x)
    H = np.empty((k, k))
    for i in range(k):
        h = eps * max(1.0, abs(x[i]))
        e = np.zeros(k)
        e[i] = h
        H[:, i] = (grad(x + e) - grad(x - e)) / (2 * h)
    return (H + H.T) / 2


# --------------------------------------------------------------------------
# Ordinal-outcome mediation (ordered probit outcome, linear mediator)
# --------------------------------------------------------------------------

@dataclass
class OrdinalMediation:
    a2: float
    beta2: float
    xi2: np.ndarray
    sigma2: float
    beta3: float
    gamma: float
    xi3: np.ndarray
    theta: np.ndarray
    levels: np.ndarray
    X: np.ndarray                 # covariates used for averaging (n x q, may be q = 0)

    def _true_params(self, rho: float):
        """Outcome-model parameters implied by an assumed error correlation rho."""
        if not -1 < rho < 1:
            raise ValueError("rho must lie strictly inside (-1, 1)")
        s = np.sqrt(1 - rho**2)
        k = rho / self.sigma2
        return (s * self.beta3 + k * self.beta2, s * self.gamma - k, s * self.xi3 + k * self.xi2,
                s * self.theta - k * self.a2)

    def expected_stage(self, t: float, t_m: float, rho: float = 0.0) -> float:
        """Sample-average E[Y(t, M(t_m))] on the stage scale (levels as coded)."""
        beta3, gamma, xi3, theta = self._true_params(rho)
        n = self.X.shape[0]
        xm = self.X @ self.xi2 if self.X.shape[1] else np.zeros(n)
        xy = self.X @ xi3 if self.X.shape[1] else np.zeros(n)
        mu = beta3 * t + gamma * (self.a2 + self.beta2 * t_m + xm) + xy
        sd = np.sqrt(1 + gamma**2 * self.sigma2**2 + 2 * gamma * rho * self.sigma2)
        pge = special.ndtr((mu[:, None] - theta[None, :]) / sd)
        steps = np.diff(self.levels)
        return float(np.mean(self.levels[0] + pge @ steps))

    def effects(self, t1: float, t0: float, rho: float = 0.0) -> dict:
        e11 = self.expected_stage(t1, t1, rho)
        e10 = self.expected_stage(t1, t0, rho)
        e00 = self.expected_stage(t0, t0, rho)
        nie, nde = e11 - e10, e10 - e00
        total = e11 - e00
        return {"nie": nie, "nde": nde, "total": total,
                "prop_mediated": nie / total if total != 0 else np.nan,
                "nie_latent_per_unit": self.beta2 * self._true_params(rho)[1]}

    @property
    def rho_star(self) -> float:
        """rho at which the indirect effect is exactly 0."""
        g = self.gamma * self.sigma2
        return float(g / np.sqrt(1 + g * g))


def fit_ordinal_mediation(T, M, Y, X=None) -> OrdinalMediation:
    T = np.asarray(T, float)
    M = np.asarray(M, float)
    n = len(T)
    Xm = np.zeros((n, 0)) if X is None else np.asarray(X, float).reshape(n, -1)
    for name, v in (("T", T), ("M", M), ("Y", np.asarray(Y, float))):
        if len(v) != n or not np.all(np.isfinite(v)):
            raise ValueError(f"{name} has the wrong length or non-finite values")
    if np.std(T) == 0 or np.std(M) == 0:
        raise ValueError("exposure or mediator has no variation")
    Z = np.column_stack([np.ones(n), T, Xm])
    c2, *_ = np.linalg.lstsq(Z, M, rcond=None)
    e2 = M - Z @ c2
    sigma2 = float(np.sqrt(np.mean(e2**2)))
    out = fit_ordinal(Y, np.column_stack([T, M, Xm]), link="probit", hessian=False)
    return OrdinalMediation(a2=float(c2[0]), beta2=float(c2[1]), xi2=c2[2:], sigma2=sigma2,
                            beta3=float(out.beta[0]), gamma=float(out.beta[1]), xi3=out.beta[2:],
                            theta=out.theta, levels=out.levels, X=Xm)


def bootstrap_ordinal_mediation(T, M, Y, X, *, t1: float, t0: float, n_boot: int, seed: int,
                                rho_grid: np.ndarray | None = None) -> dict:
    """Case-resampling bootstrap of the effects (and, optionally, of the NIE
    curve over rho). Returns arrays of replicate values."""
    rng = np.random.default_rng(seed)
    T, M, Y = (np.asarray(v, float) for v in (T, M, Y))
    n = len(T)
    X = None if X is None else np.asarray(X, float).reshape(n, -1)
    keys = ("nie", "nde", "total", "prop_mediated")
    out = {k: [] for k in keys}
    curves = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        if len(np.unique(Y[i])) < 2 or np.std(T[i]) == 0:
            continue
        fit = fit_ordinal_mediation(T[i], M[i], Y[i], None if X is None else X[i])
        e = fit.effects(t1, t0)
        for k in keys:
            out[k].append(e[k])
        if rho_grid is not None:
            curves.append([fit.effects(t1, t0, r)["nie"] for r in rho_grid])
    res = {k: np.array(v) for k, v in out.items()}
    if rho_grid is not None:
        res["nie_curve"] = np.array(curves)
    return res


def e_value(odds_ratio: float, common_outcome: bool = True) -> float:
    """E-value (VanderWeele & Ding 2017, Ann Intern Med 167:268) for an odds
    ratio; for a common outcome the OR is first converted with RR ~ sqrt(OR).
    Returns 1 when the ratio is at (or crosses) the null."""
    if not np.isfinite(odds_ratio) or odds_ratio <= 0:
        return np.nan
    rr = np.sqrt(odds_ratio) if common_outcome else odds_ratio
    rr = rr if rr >= 1 else 1 / rr
    return float(rr + np.sqrt(rr * (rr - 1)))
