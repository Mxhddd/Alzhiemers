"""JAX implementation of the Module 2 log-likelihood, with exact gradients.

Evaluates the SAME declarative edge list and the SAME integration grid as the
numpy reference (src/baserate/model.py, solve.py); tests require the two to
agree to 1e-9. Used for fitting only because it is ~100x faster and gives
exact gradients (L-BFGS-B and NUTS); the numpy version stays the reference.

Log-likelihood (multinomial, with constant):
    sum_g [ log n_g! - sum_s log n_gs! + sum_s n_gs log p_s(age_g) ]
p is floored at 1e-300 inside the log in BOTH implementations, so an
impossible observation costs ~690 log-units per brain instead of -inf (which
breaks gradient optimisers); this changes nothing for any parameter set with
non-negligible likelihood.
"""

from __future__ import annotations

import numpy as np
from scipy.special import gammaln

from src.baserate.data import StageCounts
from src.baserate.model import StageModel
from src.baserate.solve import aggregation_matrix, integration_grid

LOG_FLOOR = 1e-300
# Squarings allowed in JAX's expm (default 16). Inside the configured search
# box (log-rates <= 2, slopes <= 2 per decade) the largest ||Q h|| needs ~15;
# 24 leaves margin (a test checks the box corner is finite). Higher caps cost
# run time (64 doubled it), lower ones return NaN at the box edge.
MAX_SQUARINGS = 24


def _jax():
    import jax

    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    from jax.scipy.linalg import expm

    return jax, jnp, expm


def make_loglik(model: StageModel, data: list[StageCounts], *, step: float, start_age: float,
                eval_mode: str = "midpoint", n_points: int = 10):
    """Return (loglik, value_and_grad) — jitted functions of theta."""
    jax, jnp, expm = _jax()
    pts = [g.eval_points(eval_mode, n_points) for g in data]
    ages = np.concatenate([a for a, _ in pts])
    W = np.zeros((len(data), len(ages)))                               # group <- evaluation-age weights
    i = 0
    for gidx, (a, w) in enumerate(pts):
        W[gidx, i:i + len(a)] = w
        i += len(a)
    grid = integration_grid(ages, step=step, start_age=start_age, breakpoints=model.breakpoints)
    n = model.n_hidden
    edges = model.edges
    ei = np.array([e.i for e in edges]); ej = np.array([e.j for e in edges])
    kind = np.array([{"loglinear": 0, "gated": 1, "const": 2}[e.kind] for e in edges])
    p0 = np.array([e.p0 for e in edges]); p1 = np.array([max(e.p1, 0) for e in edges])
    u = np.array([model.u(t) for t in grid.t_mid])                     # (S,)
    A = (np.array([model.amyloid(t) for t in grid.t_mid]) if model.amyloid is not None
         else np.zeros(len(grid.t_mid)))
    M = aggregation_matrix(model)
    counts = np.array([g.counts for g in data], dtype=float)           # (G, K)
    const = float(sum(gammaln(g.n + 1) - np.sum(gammaln(g.counts + 1)) for g in data))
    rec = grid.record_after                                            # (E,) per evaluation age
    init = model.initial()

    u_j, A_j, h_j = jnp.asarray(u), jnp.asarray(A), jnp.asarray(grid.h)
    kind_j, p0_j, p1_j = jnp.asarray(kind), jnp.asarray(p0), jnp.asarray(p1)

    def rates(theta, uk, ak):
        a, b = theta[p0_j], theta[p1_j]
        loglin = jnp.exp(a + b * uk)
        gated = jnp.exp(a) * (1.0 + (jnp.exp(b) - 1.0) * ak)
        const_r = jnp.exp(a)
        return jnp.where(kind_j == 0, loglin, jnp.where(kind_j == 1, gated, const_r))

    def generator(theta, uk, ak):
        Q = jnp.zeros((n, n)).at[ei, ej].add(rates(theta, uk, ak))
        return Q - jnp.diag(Q.sum(axis=1))

    def loglik(theta):
        def body(p, xs):
            uk, ak, hk = xs
            p = p @ expm(generator(theta, uk, ak) * hk, max_squarings=MAX_SQUARINGS)
            return p, p
        _, traj = jax.lax.scan(body, jnp.asarray(init), (u_j, A_j, h_j))
        traj = jnp.vstack([jnp.asarray(init)[None, :], traj])          # row k = after k steps
        P = jnp.asarray(W) @ (traj[jnp.asarray(rec)] @ jnp.asarray(M))  # (G, K)
        ll = jnp.sum(jnp.where(counts > 0, counts * jnp.log(jnp.maximum(P, LOG_FLOOR)), 0.0)) + const
        return jnp.where(jnp.isfinite(ll), ll, -jnp.inf)

    return jax.jit(loglik), jax.jit(jax.value_and_grad(loglik))
