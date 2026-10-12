"""Sensitivity and robustness (brief §5, Module 4).

* ``sobol_indices`` — Saltelli sampling + Sobol analysis (SALib).
* ``is_u_shaped`` — the U-shape test on a dose–response curve.
* ``lhs_robustness`` — Latin hypercube over plausible parameter ranges; the
  fraction of samples in which the U-shape appears, with a Wilson 95% CI.
  If it is 3%, the caption says 3%.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy.stats import norm, qmc


def sobol_indices(f: Callable[[dict], float], ranges: dict[str, tuple[float, float]], *, n_base: int,
                  seed: int) -> dict:
    from SALib.analyze import sobol
    from SALib.sample import sobol as sobol_sample

    names = list(ranges)
    problem = {"num_vars": len(names), "names": names, "bounds": [list(ranges[n]) for n in names]}
    X = sobol_sample.sample(problem, n_base, calc_second_order=False, seed=seed)
    Y = np.array([f(dict(zip(names, x))) for x in X])
    if not np.all(np.isfinite(Y)):
        raise ValueError(f"{np.sum(~np.isfinite(Y))} non-finite model outputs in the Sobol sample")
    Si = sobol.analyze(problem, Y, calc_second_order=False, seed=seed)
    return {"names": names, "S1": Si["S1"], "S1_conf": Si["S1_conf"], "ST": Si["ST"], "ST_conf": Si["ST_conf"],
            "n_evaluations": len(Y)}


def is_u_shaped(doses: np.ndarray, outcome: np.ndarray, *, min_relative_rise: float) -> bool:
    """U-shaped (i.e. an intermediate dose is best — lower outcome = better):
    the minimum is at an interior dose, and the outcome at BOTH the lowest and
    the highest dose exceeds that minimum by more than `min_relative_rise`
    of its absolute value."""
    y = np.asarray(outcome, float)
    if len(y) < 3 or not np.all(np.isfinite(y)):
        raise ValueError("dose-response needs >= 3 finite points")
    i = int(np.argmin(y))
    if i in (0, len(y) - 1):
        return False
    floor = abs(y[i]) if y[i] != 0 else 1.0
    return bool((y[0] - y[i]) > min_relative_rise * floor and (y[-1] - y[i]) > min_relative_rise * floor)


def wilson_interval(k: int, n: int, level: float = 0.95) -> tuple[float, float]:
    z = norm.ppf(1 - (1 - level) / 2)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return float(c - h), float(c + h)


def lhs_robustness(curve: Callable[[dict], tuple[np.ndarray, np.ndarray]], ranges: dict[str, tuple[float, float]],
                   *, n_samples: int, seed: int, min_relative_rise: float) -> dict:
    names = list(ranges)
    sampler = qmc.LatinHypercube(d=len(names), seed=seed)
    U = qmc.scale(sampler.random(n_samples), [ranges[n][0] for n in names], [ranges[n][1] for n in names])
    flags = []
    for row in U:
        d, y = curve(dict(zip(names, row)))
        flags.append(is_u_shaped(d, y, min_relative_rise=min_relative_rise))
    k = int(np.sum(flags))
    lo, hi = wilson_interval(k, n_samples)
    return {"n_samples": n_samples, "n_u_shaped": k, "fraction": k / n_samples, "ci95": (lo, hi),
            "samples": U, "u_shaped": np.array(flags), "names": names}
