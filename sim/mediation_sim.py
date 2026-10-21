"""SIMULATED cohorts for Module 3 validation and the power / identifiability
study. Visible generators (ground rule 2); never read by an analysis
pipeline; nothing written to disk by this module.

simulate_cohort       -- linear-linear mediation (X -> T, X -> M, X -> Y,
                         T -> M -> Y, T -> Y), optional corr(e_M, e_Y) = rho_true.
simulate_ordinal_mediation -- linear mediator, ordered-probit outcome, optional
                         corr(e_M, e_Y*) = rho_true; returns the TRUE natural
                         indirect / direct effects on the expected-stage scale,
                         computed by Monte Carlo over potential outcomes (an
                         independent check of the closed form in src/mediation).
simulate_radc_like    -- a synthetic cohort shaped like the RADC Standard
                         dataset (autopsied ROS + MAP participants): age at
                         death, sex, education, APOE e4, amyloid, Braak stage
                         0-6, the four social / sensory exposures (two of them
                         MAP-only) and an optional LC measure. EVERY number in
                         it is an ASSUMPTION for simulation, not a RADC value.
"""

from __future__ import annotations

import numpy as np
from scipy import special

BRAAK_LEVELS = np.arange(7)          # braaksc 0..6 (0, I..VI)


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


def simulate_ordinal_mediation(n: int, *, beta2: float, gamma: float, beta3: float, rho_true: float = 0.0,
                               sigma2: float = 1.0, thresholds=(-1.5, -0.8, -0.2, 0.4, 1.0, 1.7),
                               t1: float = 1.0, t0: float = 0.0, n_truth: int = 400_000, seed: int) -> dict:
    """Continuous exposure T ~ N(0.5 X1, 1); M linear; Y* = beta3 T + gamma M + X xi3 + e3,
    e3 ~ N(0, 1), corr(e2 / sigma2, e3) = rho_true; Y = number of thresholds below Y*."""
    rng = np.random.default_rng(seed)
    th = np.asarray(thresholds, float)
    xi2, xi3 = np.array([0.3, -0.1]), np.array([-0.2, 0.4])

    def draw(m):
        X = rng.normal(size=(m, 2))
        T = 0.5 * X[:, 0] + rng.normal(size=m)
        z = rng.multivariate_normal([0, 0], [[1, rho_true], [rho_true, 1]], size=m)
        return X, T, sigma2 * z[:, 0], z[:, 1]

    X, T, e2, e3 = draw(n)
    M = 0.2 + beta2 * T + X @ xi2 + e2
    Y = np.searchsorted(th, beta3 * T + gamma * M + X @ xi3 + e3).astype(float)
    # truth by Monte Carlo over potential outcomes (unit-level errors shared across worlds)
    Xt, _, e2t, e3t = draw(n_truth)

    def ey(t, tm):
        Mt = 0.2 + beta2 * tm + Xt @ xi2 + e2t
        return float(np.mean(np.searchsorted(th, beta3 * t + gamma * Mt + Xt @ xi3 + e3t)))

    e11, e10, e00 = ey(t1, t1), ey(t1, t0), ey(t0, t0)
    return {"T": T, "M": M, "Y": Y, "X": X, "true_nie": e11 - e10, "true_nde": e10 - e00,
            "true_total": e11 - e00, "true_nie_latent_per_unit": beta2 * gamma}


# ---------------------------------------------------------------------------
# RADC-shaped synthetic cohort (assumptions only)
# ---------------------------------------------------------------------------

# Assumed marginal Braak distribution among brains that have entered the
# cortex (stages I..VI). Shape only; swept in the power study through p0.
ASSUMED_CORTICAL_SHARES = np.array([0.08, 0.12, 0.25, 0.25, 0.18, 0.12])

# Exposures: name -> (MAP only?, risk direction: +1 if a HIGHER value is the
# predicted risk side, -1 if a LOWER value is). Directions follow the
# hypothesis, not data.
EXPOSURES = {"soc_net": (False, -1), "social_isolation": (True, +1),
             "late_life_soc_act": (True, -1), "vision": (False, +1)}


def _quantile_thresholds(latent: np.ndarray, shares: np.ndarray) -> np.ndarray:
    cum = np.cumsum(shares)[:-1]
    return np.quantile(latent, cum)


def simulate_radc_like(n: int, *, seed: int, p0: float = 0.05, map_share: float = 0.5,
                       log_or_amyloid_entry: float = 0.0, log_or_amyloid_spread: float = 0.0,
                       log_or_age_spread: float = 0.3, log_or_expo: float = 0.0,
                       lc_beta2: float = 0.0, lc_gamma: float = 0.0, with_lc: bool = False) -> dict:
    """Synthetic autopsy cohort. Braak stage is generated in two steps:
    cortical ENTRY (stage >= I) by a logistic model (amyloid log-OR per SD
    `log_or_amyloid_entry`, the LC -> cortex step's closest analogue in these
    data), then SPREAD among entered brains (stages I..VI) by a proportional-
    odds logistic model (amyloid log-OR per SD `log_or_amyloid_spread`).
    Each exposure, oriented to its risk direction and standardised, has
    log-OR `log_or_expo` per SD on both steps. With `with_lc`, an LC-damage
    measure M (z-score, higher = more damage) lies on the exposure path:
    M = lc_beta2 * exposure_risk_z + e, and M enters both steps with log-OR
    `lc_gamma` per SD. All values are assumptions."""
    rng = np.random.default_rng(seed)
    study = np.where(rng.random(n) < map_share, "MAP", "ROS")
    age = rng.normal(89.0, 6.0, n)
    male = rng.binomial(1, 0.32, n).astype(float)
    educ = np.clip(rng.normal(16.0, 3.5, n), 0, 30)
    apoe4 = rng.binomial(1, 0.25, n).astype(float)
    agez = (age - 89.0) / 6.0
    amyloid = 0.35 * agez + 0.6 * apoe4 + rng.normal(0, 1, n)
    amyloid = (amyloid - amyloid.mean()) / amyloid.std()

    risk = {}
    raw = {}
    for name, (map_only, direction) in EXPOSURES.items():
        z = rng.normal(0, 1, n) + 0.1 * agez
        risk[name] = (z - z.mean()) / z.std()                # standardised, oriented to risk
        val = direction * risk[name]
        if name == "soc_net":
            raw[name] = rng.poisson(np.exp(1.8 + 0.5 * val)).astype(float)
        elif name == "vision":
            raw[name] = np.clip(np.round(2.2 + 1.1 * val), 1, 7)
        else:                                                  # 1-5 scale scores
            raw[name] = np.clip(3.0 + 0.6 * val, 1, 5)
        if map_only:
            raw[name] = np.where(study == "MAP", raw[name], np.nan)

    expo_sum = sum(risk.values()) / np.sqrt(len(risk))         # each exposure contributes equally
    M = None
    lin = log_or_age_spread * agez + 0.4 * apoe4
    if with_lc:
        M = lc_beta2 * expo_sum + rng.normal(0, 1, n)
        lin = lin + lc_gamma * M
    lin_expo = log_or_expo * sum(risk.values())

    entry_lin = log_or_amyloid_entry * amyloid + lin + lin_expo
    alpha = special.logit(1 - p0) - np.mean(entry_lin)
    # calibrate the entry intercept so the overall share at stage 0 is close to p0
    for _ in range(30):
        cur = np.mean(special.expit(alpha + entry_lin))
        alpha += (special.logit(1 - p0) - special.logit(np.clip(cur, 1e-6, 1 - 1e-6)))
    entered = rng.random(n) < special.expit(alpha + entry_lin)
    spread_latent = log_or_amyloid_spread * amyloid + lin + lin_expo + rng.logistic(0, 1, n)
    th = _quantile_thresholds(spread_latent, ASSUMED_CORTICAL_SHARES)
    stage = np.where(entered, 1 + np.searchsorted(th, spread_latent), 0).astype(float)

    out = {"study": study, "age_death": age, "msex": male, "educ": educ, "apoe4": apoe4,
           "amyloid": amyloid, "braaksc": stage, **raw}
    if with_lc:
        out["lc"] = M
    return out


def radc_like_files(n: int, *, seed: int, autopsy_share: float = 0.8, alive_share: float = 0.3,
                    lc_name: str | None = None, **kw):
    """Two DataFrames shaped like the (assumed) delivered RADC files: a
    cross-sectional table and a longitudinal one (projid, fu_year). Built from
    simulate_radc_like; participants alive or not autopsied have no pathology.
    Column names are the public RADC names the config uses. Assumptions only."""
    import pandas as pd

    rng = np.random.default_rng(seed + 1)
    s = simulate_radc_like(n, seed=seed, with_lc=lc_name is not None, **kw)
    alive = rng.random(n) < alive_share
    autopsied = (~alive) & (rng.random(n) < autopsy_share)
    n_vis = rng.integers(1, 15, n)
    age_bl = s["age_death"] - n_vis - rng.uniform(0, 2, n)
    path = lambda x: np.where(autopsied, x, np.nan)          # noqa: E731
    genos = np.where(s["apoe4"] == 1, rng.choice([24, 34, 44], n, p=[0.05, 0.85, 0.10]),
                     rng.choice([22, 23, 33], n, p=[0.02, 0.18, 0.80]))
    amy_pos = s["amyloid"] > -0.6
    basic = pd.DataFrame({
        "projid": np.arange(1, n + 1), "study": s["study"], "braaksc": path(s["braaksc"]),
        "b_score": path(np.select([s["braaksc"] == 0, s["braaksc"] <= 2, s["braaksc"] <= 4], [0, 1, 2], 3)),
        "amylsqrt_est_8reg": path(np.maximum(0, 1.5 + s["amyloid"])),
        "a_score_st4": path(np.where(amy_pos, np.clip(np.round(1.5 + s["amyloid"]), 1, 3), 0)),
        "age_death": np.where(alive, np.nan, s["age_death"]), "age_bl": age_bl, "msex": s["msex"],
        "educ": np.round(s["educ"]), "apoe_genotype": genos.astype(float), "died": (~alive).astype(float)})
    if lc_name is not None:
        basic[lc_name] = path(s["lc"])
    rows = []
    for i in range(n):
        for v in range(int(n_vis[i])):
            r = {"projid": i + 1, "fu_year": v, "dcfdx": 1.0 if v < 0.7 * n_vis[i] else 2.0,
                 "cogn_global": rng.normal(0, 0.5)}
            for name in EXPOSURES:
                x = s[name][i]
                r[name] = np.nan if np.isnan(x) else (x if v == 0 else x + rng.normal(0, 0.3))
            rows.append(r)
    return basic, pd.DataFrame(rows)
