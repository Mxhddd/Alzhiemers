"""Likelihood, maximum-likelihood and Bayesian fitting for Module 2.

The data are cross-sectional autopsy counts: in each age band, how many
brains were at each stage. The log-likelihood is multinomial:
    sum over groups g, stages s of  n_gs * log p_s(age_g)   (+ multinomial constant)
with p evaluated at the band midpoint (config: solver.evaluate_at).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln

from src.baserate.data import StageCounts
from src.baserate.model import InvalidParameters, StageModel
from src.baserate.solve import stage_probabilities


def log_likelihood(model: StageModel, theta: np.ndarray, data: list[StageCounts], *,
                   step: float, start_age: float) -> float:
    try:
        P = stage_probabilities(model, theta, [g.midpoint for g in data], step=step, start_age=start_age)
    except InvalidParameters:
        return -np.inf
    ll = 0.0
    for g, p in zip(data, P):
        pos = g.counts > 0
        if (p[pos] <= 0).any():
            return -np.inf
        ll += float(np.sum(g.counts[pos] * np.log(p[pos])))
        ll += float(gammaln(g.n + 1) - np.sum(gammaln(g.counts + 1)))
    return ll


def numerical_hessian(f, x: np.ndarray, h: float = 1e-4) -> np.ndarray:
    """Central-difference Hessian of a scalar function."""
    n = len(x)
    H = np.empty((n, n))
    for i in range(n):
        for j in range(i, n):
            e_i = np.eye(n)[i] * h
            e_j = np.eye(n)[j] * h
            H[i, j] = H[j, i] = (f(x + e_i + e_j) - f(x + e_i - e_j) - f(x - e_i + e_j) + f(x - e_i - e_j)) / (4 * h * h)
    return H


@dataclass
class MLEResult:
    model: str
    param_names: list[str]
    theta: np.ndarray
    se: np.ndarray            # NaN where the Hessian is not positive definite
    loglik: float
    n_params: int
    converged_starts: int
    n_starts: int
    optima_spread: float      # range of log-lik among converged starts within 2 units of the best
    n_at_best: int            # converged starts within 2 log-lik units of the best
    hessian_ok: bool
    messages: list[str]

    @property
    def aic(self) -> float:
        return 2 * self.n_params - 2 * self.loglik


def initial_values(model: StageModel, rng: np.random.Generator, cfg_mle: dict) -> np.ndarray:
    lo_a, hi_a = cfg_mle["init_alpha_range"]
    lo_b, hi_b = cfg_mle["init_beta_range"]
    x = {}
    for name in model.param_names:
        if name.startswith("beta") or name.startswith("b_"):
            x[name] = rng.uniform(lo_b, hi_b)
        elif name == "gamma":
            # Symmetric around zero (no favoured sign) and inside the valid
            # region mu0 + gamma * A >= 0 for every A in [0, 1].
            x[name] = rng.uniform(-0.9, 0.9) * np.exp(x["log_mu0"])
        else:
            x[name] = rng.uniform(lo_a, hi_a)
    return np.array([x[n] for n in model.param_names])


def fit_mle(model: StageModel, data: list[StageCounts], *, step: float, start_age: float,
            cfg_mle: dict, seed: int) -> MLEResult:
    rng = np.random.default_rng(seed)
    nll = lambda th: -log_likelihood(model, th, data, step=step, start_age=start_age)  # noqa: E731
    results = []
    messages = []
    for _ in range(cfg_mle["n_starts"]):
        x0 = initial_values(model, rng, cfg_mle)
        if not np.isfinite(nll(x0)):
            messages.append("start skipped: infinite likelihood at initial values")
            continue
        if cfg_mle["method"] != "nelder-mead+l-bfgs-b":
            raise ValueError(f"unknown mle.method {cfg_mle['method']!r}")
        # Stage 1: Nelder-Mead tolerates the -inf region (negative rates).
        r1 = minimize(nll, x0, method="Nelder-Mead",
                      options={"maxiter": cfg_mle["nm_maxiter"], "xatol": 1e-6, "fatol": 1e-8})
        # Stage 2: gradient polish from the Nelder-Mead point.
        r = minimize(nll, r1.x, method="L-BFGS-B")
        if not (np.isfinite(r.fun) and r.fun <= r1.fun + 1e-9):
            r = r1
        ok = bool(np.isfinite(r.fun) and (r.success or r1.success))
        messages.append(f"{'ok' if ok else 'NOT CONVERGED'}: NM {r1.message}; polish {r.message}")
        if ok:
            results.append(r)
    if not results:
        raise RuntimeError(f"{model.name}: no start converged ({messages[:3]})")
    best = min(results, key=lambda r: r.fun)
    near = [r.fun for r in results if r.fun - best.fun < 2]
    if len(near) < 2:
        messages.append(f"WARNING: only {len(near)} of {len(results)} converged starts reached the best optimum "
                        "(within 2 log-lik units): the optimum may not be global")
    H = numerical_hessian(nll, best.x)
    try:
        cov = np.linalg.inv(H)
        se = np.sqrt(np.diag(cov))
        ok = bool(np.all(np.linalg.eigvalsh(H) > 0) and np.all(np.isfinite(se)))
    except np.linalg.LinAlgError:
        se, ok = np.full(len(best.x), np.nan), False
    if not ok:
        se = np.full(len(best.x), np.nan)
    return MLEResult(model=model.name, param_names=model.param_names, theta=best.x, se=se,
                     loglik=-best.fun, n_params=len(best.x), converged_starts=len(results),
                     n_starts=cfg_mle["n_starts"], optima_spread=float(max(near) - min(near)), n_at_best=len(near),
                     hessian_ok=ok, messages=messages)


def fit_bayesian(model: StageModel, data: list[StageCounts], *, step: float, start_age: float,
                 cfg_bayes: dict, seed: int):
    """NUTS via PyMC with a black-box likelihood (finite-difference gradient).

    Returns (InferenceData, diagnostics) where diagnostics lists max R-hat,
    min bulk ESS and the number of divergent transitions — always reported."""
    import arviz as az
    import pymc as pm
    import pytensor.tensor as pt
    from pytensor.graph import Apply, Op

    def ll(th):
        return log_likelihood(model, np.asarray(th, dtype=float), data, step=step, start_age=start_age)

    def grad(th, h=1e-5):
        g = np.empty(len(th))
        for i in range(len(th)):
            e = np.zeros(len(th)); e[i] = h
            g[i] = (ll(th + e) - ll(th - e)) / (2 * h)
        return g

    class LogLikeGrad(Op):
        def make_node(self, theta):
            theta = pt.as_tensor_variable(theta)
            return Apply(self, [theta], [theta.type()])

        def perform(self, node, inputs, outputs):
            outputs[0][0] = grad(inputs[0])

    class LogLike(Op):
        def make_node(self, theta):
            theta = pt.as_tensor_variable(theta)
            return Apply(self, [theta], [pt.dscalar()])

        def perform(self, node, inputs, outputs):
            outputs[0][0] = np.asarray(ll(inputs[0]))

        def grad(self, inputs, g):
            return [g[0] * LogLikeGrad()(inputs[0])]

    names = model.param_names
    with pm.Model():
        comps = []
        for name in names:
            if name.startswith("beta") or name.startswith("b_"):
                pr = cfg_bayes["beta_prior"]
            elif name == "gamma":
                pr = cfg_bayes["gamma_prior"]
            else:
                pr = cfg_bayes["alpha_prior"]
            comps.append(pm.Normal(name.replace("[", "_").replace("]", ""), mu=pr["mu"], sigma=pr["sigma"]))
        theta = pt.stack(comps)
        pm.Potential("loglik", LogLike()(theta))
        idata = pm.sample(draws=cfg_bayes["draws"], tune=cfg_bayes["tune"], chains=cfg_bayes["chains"],
                          target_accept=cfg_bayes["target_accept"], random_seed=seed, progressbar=False,
                          cores=1)
    summ = az.summary(idata, hdi_prob=cfg_bayes["ci_level"])
    diagnostics = {
        "max_rhat": float(summ["r_hat"].max()),
        "min_ess_bulk": float(summ["ess_bulk"].min()),
        "divergences": int(idata.sample_stats["diverging"].sum()),
    }
    return idata, summ, diagnostics
