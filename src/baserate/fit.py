"""Likelihood, maximum-likelihood and Bayesian fitting for Module 2.

The data are cross-sectional autopsy counts: in each age band, how many
brains were at each stage. The log-likelihood is multinomial:
    sum over groups g, stages s of  n_gs * log p_s(age_g)   (+ multinomial constant)
with p evaluated at the band midpoint (config: solver.evaluate_at).

``log_likelihood`` is the numpy reference. Fitting uses the JAX version in
jaxlik.py (same edges, same grid, exact gradients); tests require agreement.

Parameter roles (decide initial values, bounds and priors):
  log-rate intercepts:  alpha[k], a_nu, log_mu0, log_lambda
  age slopes (per age_scale years):  beta[k], b_nu, beta_mu
  amyloid log hazard ratio:  log_hr_amyloid
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln

from src.baserate.data import StageCounts
from src.baserate.jaxlik import LOG_FLOOR, make_loglik
from src.baserate.model import InvalidParameters, StageModel
from src.baserate.solve import stage_probabilities


def role(name: str) -> str:
    if name.startswith("beta[") or name in ("b_nu", "beta_mu"):
        return "slope"
    if name == "log_hr_amyloid":
        return "log_hr"
    if name.startswith("alpha[") or name in ("a_nu", "log_mu0", "log_lambda"):
        return "log_rate"
    raise ValueError(f"no role defined for parameter {name}")


def group_probabilities(model: StageModel, theta: np.ndarray, data: list[StageCounts], *, step: float,
                        start_age: float, eval_mode: str = "midpoint", n_points: int = 10) -> np.ndarray:
    """Stage probabilities per observation group (G, K), at the band midpoint
    or averaged over the band, per eval_mode."""
    pts = [g.eval_points(eval_mode, n_points) for g in data]
    ages = np.concatenate([a for a, _ in pts])
    P = stage_probabilities(model, theta, ages, step=step, start_age=start_age)
    out, i = [], 0
    for a, w in pts:
        out.append(w @ P[i:i + len(a)])
        i += len(a)
    return np.array(out)


def log_likelihood(model: StageModel, theta: np.ndarray, data: list[StageCounts], *,
                   step: float, start_age: float, eval_mode: str = "midpoint", n_points: int = 10) -> float:
    """Numpy reference implementation (see jaxlik.py for the floor on p)."""
    try:
        P = group_probabilities(model, theta, data, step=step, start_age=start_age, eval_mode=eval_mode,
                                n_points=n_points)
    except InvalidParameters:
        return -np.inf
    ll = 0.0
    for g, p in zip(data, P):
        pos = g.counts > 0
        ll += float(np.sum(g.counts[pos] * np.log(np.maximum(p[pos], LOG_FLOOR))))
        ll += float(gammaln(g.n + 1) - np.sum(gammaln(g.counts + 1)))
    return ll


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
    on_bound: list[str]       # parameters whose estimate sits on the search box
    messages: list[str]
    cov: np.ndarray | None = None   # inverse observed information (None if the Hessian is not positive definite)
    verified: bool = True           # False: no start met the convergence test (allow_unverified fits only)

    @property
    def aic(self) -> float:
        return 2 * self.n_params - 2 * self.loglik


def _free(g: np.ndarray, x: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Parameters not held at a bound by the gradient (minimisation)."""
    at_lo = (x <= lo + 1e-10) & (g > 0)
    at_hi = (x >= hi - 1e-10) & (g < 0)
    return ~(at_lo | at_hi)


def newton_polish(f, hess, x, lo, hi, *, steps: int):
    """Damped Newton steps on the free parameters with the exact Hessian of
    the negative log-likelihood. Returns (x, f(x), Newton decrement, steps).
    The decrement 1/2 g' H^-1 g is the predicted remaining improvement in
    log-likelihood; inf if the Hessian of the free block is not positive
    definite (then the point is not a verified local optimum)."""
    fx, g = f(x)
    n = 0
    for n in range(steps + 1):
        free = _free(g, x, lo, hi)
        H = -np.asarray(hess(x), dtype=float)[np.ix_(free, free)]    # Hessian of the NEGATIVE log-lik
        try:
            np.linalg.cholesky(H)
        except np.linalg.LinAlgError:
            return x, fx, np.inf, n
        d = np.linalg.solve(H, g[free])
        dec = 0.5 * float(g[free] @ d)
        if dec <= 1e-12 or n == steps:
            return x, fx, dec, n
        t = 1.0
        while t > 1e-6:
            xn = x.copy()
            xn[free] = np.clip(x[free] - t * d, lo[free], hi[free])
            fn, gn = f(xn)
            if fn <= fx - 1e-4 * t * float(g[free] @ d):
                x, fx, g = xn, fn, gn
                break
            t /= 2
        else:
            return x, fx, dec, n                                    # no further progress possible
    return x, fx, dec, n


def bounds_for(model: StageModel, cfg_mle: dict) -> list[tuple[float, float]]:
    return [tuple(cfg_mle["bounds"][role(n)]) for n in model.param_names]


def initial_values(model: StageModel, rng: np.random.Generator, cfg_mle: dict) -> np.ndarray:
    return np.array([rng.uniform(*cfg_mle["init"][role(n)]) for n in model.param_names])


def fit_mle(model: StageModel, data: list[StageCounts], *, step: float, start_age: float,
            cfg_mle: dict, seed: int, eval_mode: str = "midpoint", n_points: int = 10,
            extra_starts: list[np.ndarray] | None = None, n_random: int | None = None,
            compiled: tuple | None = None, allow_unverified: bool = False) -> MLEResult:
    """Multi-start L-BFGS-B on exact (JAX) gradients inside the configured
    box; SEs from the exact Hessian. An optimum on the box edge is reported.
    allow_unverified: if no start passes the convergence test, return the
    best finite start (by log-likelihood) with verified=False instead of
    raising. Used only for cross-validation folds, where a training subset
    can leave a direction flat (prereg §5); the flag is reported.
    `compiled` = (loglik, value_and_grad, hessian) functions of theta only,
    already bound to `data` (e.g. from jaxlik.make_loglik_family); default:
    compile them here."""
    import jax

    if cfg_mle["method"] != "l-bfgs-b":
        raise ValueError(f"unknown mle.method {cfg_mle['method']!r}")
    if compiled is None:
        ll, vg = make_loglik(model, data, step=step, start_age=start_age, eval_mode=eval_mode, n_points=n_points)
        hess = jax.jit(jax.hessian(ll))
    else:
        ll, vg, hess = compiled

    def f(th):
        v, g = vg(th)
        v, g = float(v), np.asarray(g, dtype=float)
        if not np.isfinite(v) or not np.all(np.isfinite(g)):
            return 1e300, np.zeros_like(th)  # outside the finite region: reject the step
        return -v, -g

    rng = np.random.default_rng(seed)
    box = bounds_for(model, cfg_mle)
    lo_b = np.array([b[0] for b in box]); hi_b = np.array([b[1] for b in box])
    results, messages, unverified = [], [], []
    starts = list(extra_starts or []) + [initial_values(model, rng, cfg_mle)
                                         for _ in range(cfg_mle["n_starts"] if n_random is None else n_random)]
    for x0 in starts:
        x0 = np.clip(np.asarray(x0, float), lo_b, hi_b)
        r = minimize(f, x0, jac=True, method="L-BFGS-B", bounds=box,
                     options={"maxiter": cfg_mle["maxiter"], "ftol": 1e-15, "gtol": 1e-10})
        if not r.fun < 1e299:
            messages.append(f"NOT CONVERGED: {r.message} (non-finite objective)")
            continue
        x, fx, dec, n_newton = newton_polish(f, hess, r.x, lo_b, hi_b, steps=cfg_mle["newton_steps"])
        # Converged = the predicted remaining gain to the optimum (Newton
        # decrement) is negligible. The optimiser's own message is not trusted:
        # L-BFGS-B can report convergence far from an optimum.
        ok = bool(np.isfinite(dec) and dec <= cfg_mle["max_newton_decrement"])
        messages.append(f"{'ok' if ok else 'NOT CONVERGED'}: L-BFGS-B {r.nit} it ({r.message}); "
                        f"{n_newton} Newton steps; remaining gain {dec:.2e} log-lik")
        if ok:
            results.append(type(r)(x=x, fun=fx))
        elif np.isfinite(fx):
            unverified.append(type(r)(x=x, fun=fx))
    verified = bool(results)
    if not results:
        if not (allow_unverified and unverified):
            raise RuntimeError(f"{model.name}: no start converged ({messages[:3]})")
        results = unverified
        messages.append("WARNING: no start passed the convergence test; the best finite start is used (verified = False)")
    best = min(results, key=lambda r: r.fun)
    near = [r.fun for r in results if r.fun - best.fun < 2]
    if len(near) < 2:
        messages.append(f"WARNING: only {len(near)} of {len(results)} converged starts reached the best optimum "
                        "(within 2 log-lik units): the optimum may not be global")
    on_bound = [n for n, x, (lo, hi) in zip(model.param_names, best.x, box)
                if min(x - lo, hi - x) < 1e-6 * max(1.0, hi - lo)]
    if on_bound:
        messages.append(f"WARNING: estimate on the search-box edge for {on_bound}: not an interior optimum")
    H = -np.asarray(hess(best.x), dtype=float)
    cov = None
    try:
        eig = np.linalg.eigvalsh(H)
        ok = bool(np.all(eig > 0))
        if ok:
            cov = np.linalg.inv(H)
            se = np.sqrt(np.diag(cov))
        else:
            se = np.full(len(best.x), np.nan)
    except np.linalg.LinAlgError:
        ok, se = False, np.full(len(best.x), np.nan)
    return MLEResult(model=model.name, param_names=model.param_names, theta=best.x, se=se,
                     loglik=-best.fun, n_params=len(best.x), converged_starts=len(results),
                     n_starts=len(starts), optima_spread=float(max(near) - min(near)), n_at_best=len(near),
                     hessian_ok=ok, on_bound=on_bound, messages=messages, cov=cov, verified=verified)


def fit_bayesian(model: StageModel, data: list[StageCounts], *, step: float, start_age: float,
                 cfg_bayes: dict, seed: int, init_theta: np.ndarray, eval_mode: str = "midpoint",
                 n_points: int = 10, chain_cache=None):
    """NUTS posterior (NumPyro: the same NUTS algorithm as PyMC, compiled end
    to end in JAX together with the likelihood).

    Geometry: parameters differ in posterior scale by ~100x and slopes are
    correlated with intercepts, which drives a diagonally-adapted NUTS to its
    maximum tree depth. So the initial DENSE inverse mass matrix is the
    Laplace covariance (inverse Hessian of the log-posterior) at `init_theta`
    (the MLE), adaptation continues during warm-up, and chains start at
    over-dispersed points drawn from N(MLE, (k^2) * Laplace covariance) with
    k = cfg_bayes["init_overdispersion"] — dispersed starts keep R-hat able to
    detect multimodality. Chains run in parallel when enough devices exist.

    Returns (InferenceData, summary, diagnostics); diagnostics give max
    R-hat, min bulk ESS, divergences and mean leapfrog steps — always reported."""
    import arviz as az
    import jax
    import jax.numpy as jnp
    import numpyro
    import numpyro.distributions as dist
    from numpyro.infer import MCMC, NUTS

    ll, _ = make_loglik(model, data, step=step, start_age=start_age, eval_mode=eval_mode, n_points=n_points)
    prior_of = {"log_rate": cfg_bayes["log_rate_prior"], "slope": cfg_bayes["slope_prior"],
                "log_hr": cfg_bayes["log_hr_prior"]}
    names = [n.replace("[", "_").replace("]", "") for n in model.param_names]
    mu = jnp.array([prior_of[role(n)]["mu"] for n in model.param_names])
    sd = jnp.array([prior_of[role(n)]["sigma"] for n in model.param_names])

    def logpost(th):
        return ll(th) + jnp.sum(dist.Normal(mu, sd).log_prob(th))

    H = -np.asarray(jax.hessian(logpost)(jnp.asarray(init_theta)), dtype=float)
    try:
        cov = np.linalg.inv(H)
        np.linalg.cholesky(cov)
    except np.linalg.LinAlgError as e:
        raise RuntimeError(f"{model.name}: log-posterior Hessian at the MLE is not positive definite; "
                           "cannot build the NUTS mass matrix") from e

    def npmodel():
        # One vector site: its dense mass matrix is unambiguous (a dict keyed by
        # many separate sites must match NumPyro's internal site order).
        th = numpyro.sample("theta", dist.Normal(mu, sd))
        numpyro.factor("loglik", ll(th))

    chains = cfg_bayes["chains"]
    rng = np.random.default_rng(seed)
    k = cfg_bayes["init_overdispersion"]
    starts = rng.multivariate_normal(init_theta, (k ** 2) * cov, size=chains)
    kernel = NUTS(npmodel, target_accept_prob=cfg_bayes["target_accept"], dense_mass=True,
                  inverse_mass_matrix=jnp.asarray(cov))
    if chain_cache is None:
        method = "parallel" if jax.local_device_count() >= chains else "sequential"
        mc = MCMC(kernel, num_warmup=cfg_bayes["tune"], num_samples=cfg_bayes["draws"], num_chains=chains,
                  chain_method=method, progress_bar=False)
        mc.run(jax.random.PRNGKey(seed), init_params={"theta": jnp.asarray(starts)},
               extra_fields=("diverging", "num_steps"))
        draws = np.asarray(mc.get_samples(group_by_chain=True)["theta"])          # (chains, draws, P)
        extra = mc.get_extra_fields(group_by_chain=True)
    else:
        # One chain at a time, each saved when it finishes, so a long run
        # resumes per chain. Chain i gets the key NumPyro itself derives for
        # it (random.split of the run key) and the same start as above.
        method = "per-chain (checkpointed)"
        keys = jax.random.split(jax.random.PRNGKey(seed), chains)
        per = []
        for i in range(chains):
            rec = chain_cache[0](i)
            if rec is None:
                mc = MCMC(kernel, num_warmup=cfg_bayes["tune"], num_samples=cfg_bayes["draws"], num_chains=1,
                          progress_bar=False)
                mc.run(keys[i], init_params={"theta": jnp.asarray(starts[i])},
                       extra_fields=("diverging", "num_steps"))
                ex = mc.get_extra_fields()
                rec = {"draws": np.asarray(mc.get_samples()["theta"]),
                       "diverging": np.asarray(ex["diverging"]), "num_steps": np.asarray(ex["num_steps"])}
                chain_cache[1](i, rec)
            per.append(rec)
        draws = np.stack([np.asarray(r["draws"]) for r in per])
        extra = {"diverging": np.stack([np.asarray(r["diverging"]) for r in per]),
                 "num_steps": np.stack([np.asarray(r["num_steps"]) for r in per])}
    idata = az.from_dict(posterior={nm: draws[:, :, i] for i, nm in enumerate(names)},
                         sample_stats={"diverging": np.asarray(extra["diverging"]),
                                       "n_steps": np.asarray(extra["num_steps"])})
    summ = az.summary(idata, var_names=names, hdi_prob=cfg_bayes["ci_level"])
    diagnostics = {
        "max_rhat": float(summ["r_hat"].max()),
        "min_ess_bulk": float(summ["ess_bulk"].min()),
        "divergences": int(np.asarray(extra["diverging"]).sum()),
        "mean_leapfrog_steps": float(np.asarray(extra["num_steps"]).mean()),
        "chain_method": method,
    }
    diagnostics["draws"] = draws.reshape(-1, draws.shape[-1])        # (chains*draws, P), for predictions
    return idata, summ, diagnostics
