"""Predictions and model comparison for Module 2."""

from __future__ import annotations

import numpy as np

from src.baserate.data import StageCounts
from src.baserate.fit import fit_mle, log_likelihood
from src.baserate.model import StageModel, build_generator
from src.baserate.solve import (
    hidden_probabilities,
    integration_grid,
)


def predict_incidence(model: StageModel, theta: np.ndarray, ages, *, threshold: int,
                      step: float, start_age: float) -> np.ndarray:
    """Hazard of first entering any stage >= `threshold` at each age:
    flux into the set / probability of being outside it.

    Stages never regress in these models, so the flux out of stages
    < threshold into stages >= threshold is the incidence numerator exactly.
    Hidden states (variant C) are handled through observed_of_hidden. Uses
    the same integration grid as the likelihood (breakpoints respected).
    """
    ages = np.asarray(ages, dtype=float)
    grid = integration_grid(ages, step=step, start_age=start_age, breakpoints=model.breakpoints)
    traj = hidden_probabilities(model, theta, grid)
    below = model.observed_of_hidden < threshold
    out = np.empty(len(ages))
    for idx, (age, k) in enumerate(zip(ages, grid.record_after)):
        p = traj[k - 1] if k > 0 else model.initial()
        Q = build_generator(model, age, theta)
        flux = float(p[below] @ Q[np.ix_(below, ~below)].sum(axis=1))
        at_risk = float(p[below].sum())
        if at_risk <= 0:
            raise ValueError(f"nobody at risk at age {age}")
        out[idx] = flux / at_risk
    return out


def stage_prevalence_table(model, theta, data: list[StageCounts], *, step, start_age, eval_mode="midpoint",
                           n_points=10):
    """Predicted vs observed stage proportions per age band."""
    from src.baserate.fit import group_probabilities
    P = group_probabilities(model, theta, data, step=step, start_age=start_age, eval_mode=eval_mode,
                            n_points=n_points)
    rows = []
    for g, p in zip(data, P):
        for k, (n, pk) in enumerate(zip(g.counts, p)):
            rows.append({"age_lower": g.age_lower, "age_upper": g.age_upper, "source": g.source,
                         "stage_index": k, "observed_n": int(n), "observed_prop": n / g.n,
                         "predicted_prop": float(pk), "band_n": g.n})
    return rows


def age_band_folds(data: list[StageCounts], k: int) -> list[int]:
    """Fold index per group. Distinct age bands are ranked by midpoint and dealt
    round-robin (rank % k), so every fold spans the age range — the test is of
    interpolation between bands. Groups sharing a band share a fold."""
    bands = sorted({(g.age_lower, g.age_upper) for g in data}, key=lambda b: (b[0] + b[1]) / 2)
    if len(bands) < k:
        raise ValueError(f"{len(bands)} distinct age bands but k = {k} folds")
    fold_of = {b: i % k for i, b in enumerate(bands)}
    return [fold_of[(g.age_lower, g.age_upper)] for g in data]


def compare_models(models: list[StageModel], data: list[StageCounts], *, k: int, step: float,
                   start_age: float, cfg_mle: dict, seed: int, eval_mode: str = "midpoint",
                   n_points: int = 10, warm_starts: dict | None = None, n_random: int | None = None) -> list[dict]:
    """k-fold cross-validated held-out log-likelihood per model (higher is better).
    warm_starts: {model name: theta} added as a start in every fold (the
    full-data MLE); n_random: number of additional random starts per fold."""
    folds = age_band_folds(data, k)
    out = []
    for m in models:
        held = 0.0
        per_fold = []
        for f in range(k):
            train = [g for g, fo in zip(data, folds) if fo != f]
            test = [g for g, fo in zip(data, folds) if fo == f]
            ws = [warm_starts[m.name]] if warm_starts and m.name in warm_starts else None
            fit = fit_mle(m, train, step=step, start_age=start_age, cfg_mle=cfg_mle, seed=seed + f,
                          eval_mode=eval_mode, n_points=n_points, extra_starts=ws, n_random=n_random)
            ll = log_likelihood(m, fit.theta, test, step=step, start_age=start_age, eval_mode=eval_mode,
                                n_points=n_points)
            per_fold.append(ll)
            held += ll
        out.append({"model": m.name, "cv_heldout_loglik": held, "per_fold": per_fold, "k": k})
    return out
