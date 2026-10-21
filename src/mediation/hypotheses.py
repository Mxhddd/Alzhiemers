"""Module 3 hypothesis tests on the RADC Standard dataset (see prereg/mediation.md).

H1  amyloid load -> higher Braak stage (proportional odds), given age at death,
    sex, education, APOE e4 and cohort. Family A.
H2  amyloid load -> cortical entry (Braak >= I versus 0), same adjustment:
    the split closest to the LC -> cortex step that silver-stain Braak
    staging can see. Family A.
H3  each social / sensory exposure, oriented to its predicted risk side ->
    higher Braak stage. Family B (four tests).
H4  (only if LC measures are delivered) natural indirect effect of each
    exposure on expected Braak stage through the LC measure. Family C.

Every association is a log odds ratio per SD of the exposure, oriented so that
a positive value is the predicted direction. Verdicts use Bonferroni-adjusted
intervals within each family.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from src.mediation.ordinal import e_value, fit_ordinal


def adjusted_z(level: float, m: int) -> float:
    return float(stats.norm.ppf(1 - (1 - level) / (2 * m)))


def design(df: pd.DataFrame, covariates: list[str]) -> tuple[np.ndarray, list[str], list[str]]:
    """Covariate matrix; constant columns (e.g. the cohort indicator in a
    MAP-only analysis) are dropped and reported, never silently kept."""
    keep, dropped = [], []
    for c in covariates:
        (keep if df[c].nunique() > 1 else dropped).append(c)
    return df[keep].to_numpy(float), keep, dropped


def association(df: pd.DataFrame, *, outcome: str, exposure: str, direction: int, covariates: list[str],
                link: str = "logit", split: float | None = None, weights: np.ndarray | None = None) -> dict:
    x = df[exposure].to_numpy(float)
    sd = float(np.std(x, ddof=1))
    if not sd > 0:
        raise ValueError(f"{exposure} has no variation in the analysis sample")
    z = direction * (x - x.mean()) / sd
    y = df[outcome].to_numpy(float)
    if split is not None:
        y = (y >= split).astype(float)
    Xc, kept, dropped = design(df, covariates)
    fit = fit_ordinal(y, np.column_stack([z, Xc]), link=link, weights=weights,
                      names=[exposure] + kept, hessian=weights is None)
    b, se = fit.coef(exposure) if weights is None else (float(fit.beta[0]), np.nan)
    levels, counts = np.unique(y, return_counts=True)
    return {"exposure": exposure, "log_or": b, "se": se, "n": int(len(y)), "sd": sd,
            "levels": {float(k): int(c) for k, c in zip(levels, counts)}, "covariates": kept,
            "dropped_constant": dropped, "converged": fit.converged, "fit": fit}


def verdict(log_or: float, se: float, *, m: int, level: float, equivalence: float, max_width: float) -> dict:
    """Pre-registered three-way reading (plus 'uninformative'), on the
    Bonferroni-adjusted interval for a family of m tests."""
    z95, za = adjusted_z(level, 1), adjusted_z(level, m)
    lo, hi = log_or - za * se, log_or + za * se
    if lo > 0:
        v = "supports"
    elif hi < 0:
        v = "AGAINST"
    elif -equivalence < lo and hi < equivalence:
        v = "AGAINST (practically null)"
    elif hi - lo > max_width:
        v = "uninformative"
    else:
        v = "inconclusive"
    return {"or": float(np.exp(log_or)), "ci95": (float(np.exp(log_or - z95 * se)), float(np.exp(log_or + z95 * se))),
            "ci_adj": (float(np.exp(lo)), float(np.exp(hi))), "adj_level": 1 - (1 - level) / m,
            "p": float(2 * stats.norm.sf(abs(log_or / se))), "verdict": v,
            "e_value": e_value(float(np.exp(log_or))),
            "e_value_ci": 1.0 if lo <= 0 else e_value(float(np.exp(lo)))}


def mediation_verdict(nie: float, boot: np.ndarray, *, m: int, level: float, predicted_sign: int = 1) -> dict:
    a = (1 - level) / m
    boot = boot[np.isfinite(boot)]
    lo, hi = np.quantile(boot, [a / 2, 1 - a / 2])
    lo95, hi95 = np.quantile(boot, [(1 - level) / 2, 1 - (1 - level) / 2])
    if (lo > 0 and predicted_sign > 0) or (hi < 0 and predicted_sign < 0):
        v = "NIE in the predicted direction: Killer 2 NOT met for this exposure"
    elif lo > 0 or hi < 0:
        v = "NIE OPPOSITE to prediction: Killer 2 MET (hypothesis fails) for this exposure"
    else:
        v = "NIE interval includes 0: Killer 2 MET (no evidence of mediation by the LC) for this exposure"
    return {"nie": float(nie), "ci95": (float(lo95), float(hi95)), "ci_adj": (float(lo), float(hi)),
            "n_boot_used": int(len(boot)), "verdict": v}


def rubin(estimates: list[float], ses: list[float]) -> tuple[float, float]:
    q = np.asarray(estimates)
    u = np.asarray(ses) ** 2
    m = len(q)
    t = u.mean() + (1 + 1 / m) * q.var(ddof=1)
    return float(q.mean()), float(np.sqrt(t))


def multiply_imputed_association(df: pd.DataFrame, *, columns: list[str], n_imputations: int, cycles: int,
                                 seed: int, **kw) -> dict:
    """Chained-equations imputation (statsmodels MICEData, predictive mean
    matching) over `columns` (outcome, exposure, covariates); the association
    model is fitted in each completed data set and pooled by Rubin's rules."""
    from statsmodels.imputation.mice import MICEData

    np.random.seed(seed)   # MICEData draws from numpy's global generator
    data = df[columns].astype(float).reset_index(drop=True)
    imp = MICEData(data)
    est, ses = [], []
    for _ in range(n_imputations):
        imp.update_all(cycles)
        done = imp.data.copy()
        done[kw["outcome"]] = np.round(done[kw["outcome"]])   # outcome stays on its ordinal codes
        r = association(done, **kw)
        est.append(r["log_or"])
        ses.append(r["se"])
    b, se = rubin(est, ses)
    return {"log_or": b, "se": se, "n": int(len(data)), "m": n_imputations}


def autopsy_weights(people: pd.DataFrame, *, predictors: list[str], observed: str, truncate_q: float) -> pd.Series:
    """Stabilised inverse-probability-of-autopsy weights among deceased
    participants: P(Braak observed | baseline predictors), logistic model,
    truncated at the `truncate_q` quantile."""
    import statsmodels.api as sm

    d = people.dropna(subset=predictors)
    X = sm.add_constant(d[predictors].astype(float), has_constant="add")
    r = d[observed].astype(float)
    p = sm.Logit(r, X).fit(disp=0).predict(X)
    w = r.mean() / p
    w = np.minimum(w, np.quantile(w[r == 1], truncate_q))
    return w[r == 1]
