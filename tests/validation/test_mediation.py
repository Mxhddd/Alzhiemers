"""Module 3 validation: exact identities, recovery and coverage on simulated
cohorts (sim/mediation_sim.py), the rho sensitivity against a known
confounding correlation, the ordinal (Braak-stage) estimators against
simulated truth, and reproduction of two published mediation analyses from
the authors' public replication data (tests/validation/targets/mediation_published.yaml;
data: python -m src.mediation.fetch_validation_data)."""

import numpy as np
import pandas as pd
import pytest

from sim.mediation_sim import simulate_cohort, simulate_ordinal_mediation, simulate_radc_like
from src.common.config import load_config
from src.common.io import analysis_input
from src.common.paths import REPO_ROOT
from src.mediation.hypotheses import association
from src.mediation.mediation import bootstrap, estimate, killer2_verdict, percentile_ci, summarise
from src.mediation.ordinal import bootstrap_ordinal_mediation, fit_ordinal, fit_ordinal_mediation

pytestmark = pytest.mark.validation

TARGETS = REPO_ROOT / "tests/validation/targets/mediation_published.yaml"
FETCH = "public CC0 replication data: run `python -m src.mediation.fetch_validation_data`"


# ---- linear-linear estimator (used for the published checks and sensitivity S3) ----

def test_acme_at_rho_zero_is_the_product_of_ols_coefficients_and_vanishes_at_rho_tilde():
    d = simulate_cohort(2000, beta2=0.5, gamma=0.4, beta3=0.3, seed=1)
    pe = estimate(d["T"], d["M"], d["Y"], d["X"])
    assert pe.acme_at(np.array([0.0]))[0] == pytest.approx(pe.beta2 * pe.gamma, rel=1e-10)
    assert pe.acme_at(np.array([pe.rho_tilde]))[0] == pytest.approx(0.0, abs=1e-12)
    assert pe.total == pytest.approx(pe.acme + pe.ade, rel=1e-12)


def test_recovery_and_ci_coverage_when_sequential_ignorability_holds():
    covered, errs = 0, []
    reps = 150
    for r in range(reps):
        d = simulate_cohort(400, beta2=0.5, gamma=0.4, beta3=0.3, seed=100 + r)
        pe = estimate(d["T"], d["M"], d["Y"], d["X"])
        lo, hi = percentile_ci([b.acme for b in bootstrap(d["T"], d["M"], d["Y"], d["X"], n_boot=300, seed=r)], 0.95)
        covered += lo <= d["true_acme"] <= hi
        errs.append(pe.acme - d["true_acme"])
    assert abs(np.mean(errs)) < 0.01                   # unbiased
    assert 0.89 <= covered / reps <= 0.99, covered / reps  # nominal 95%


def test_sensitivity_recovers_truth_at_the_true_rho_while_naive_estimate_is_biased():
    d = simulate_cohort(40000, beta2=0.5, gamma=0.4, beta3=0.3, rho_true=0.4, seed=7)
    pe = estimate(d["T"], d["M"], d["Y"], d["X"])
    assert abs(pe.acme - d["true_acme"]) > 0.1           # confounding biases the naive ACME
    assert pe.acme_at(np.array([0.4]))[0] == pytest.approx(d["true_acme"], abs=0.02)


def test_a_null_mediator_triggers_killer2():
    d = simulate_cohort(3000, beta2=0.5, gamma=0.0, beta3=0.3, seed=9)   # LC integrity does not affect tau
    pe = estimate(d["T"], d["M"], d["Y"], d["X"])
    s = summarise(pe, bootstrap(d["T"], d["M"], d["Y"], d["X"], n_boot=500, seed=9), level=0.95,
                  rho_grid=np.arange(-0.9, 0.91, 0.05))
    assert "Killer 2 MET" in killer2_verdict(s, predicted_sign=1)


# ---- published reproductions (authors' replication data, CC0) ----

def _load(rel: str) -> pd.DataFrame:
    df = pd.read_csv(analysis_input(REPO_ROOT / rel, provided_by=FETCH))
    if df.isna().any().any():
        raise AssertionError(f"{rel} has missing values; the published analyses used complete data")
    return df


def test_reproduces_imai_keele_tingley_2010_table1_jobs_ii():
    t = load_config(TARGETS)
    p, inp = t.get("published"), t.get("inputs")
    df = _load(inp["data"])
    X = pd.get_dummies(df[inp["covariates"]], columns=inp["factors"], drop_first=True).astype(float).to_numpy()
    pe = estimate(df[inp["exposure"]].to_numpy(float), df[inp["mediator"]].to_numpy(float),
                  df[inp["outcome"]].to_numpy(float), X)
    assert len(df) == 1285
    assert abs(pe.acme - p["acme"]) <= p["acme_tolerance"], (pe.acme, p["acme"])
    assert p["acme_ci_lower"] <= pe.acme <= p["acme_ci_upper"]
    assert abs(pe.ade - p["ade"]) <= p["ade_tolerance"], (pe.ade, p["ade"])
    assert abs(pe.beta2 - p["effect_on_mediator"]) <= p["effect_on_mediator_tolerance"], pe.beta2


def test_reproduces_imai_keele_yamamoto_2010_sensitivity_rho_nelson_framing():
    s = load_config(TARGETS).get("sensitivity")
    df = _load(s["data"])
    pe = estimate(df[s["exposure"]].to_numpy(float), df[s["mediator"]].to_numpy(float),
                  df[s["outcome"]].to_numpy(float), None)
    assert abs(pe.acme - s["acme"]) <= s["acme_tolerance"], pe.acme
    assert abs(pe.total - s["total"]) <= s["acme_tolerance"], pe.total
    # published: the ACME keeps its sign "unless rho is less than -0.68" (rho grid step 0.01)
    assert abs(pe.rho_tilde - s["rho_at_zero_acme"]) <= s["rho_tolerance"], pe.rho_tilde
    grid = np.round(np.arange(-0.95, 0.951, 0.01), 2)
    curve = pe.acme_at(grid)
    assert np.all(curve[grid >= -0.67] < 0) and np.all(curve[grid <= -0.68] > 0)


# ---- ordinal (Braak stage) estimators against simulated truth ----

def test_proportional_odds_recovers_known_log_odds_ratio_with_nominal_coverage():
    covered, errs = 0, []
    reps = 120
    for r in range(reps):
        d = pd.DataFrame(simulate_radc_like(1000, seed=500 + r, p0=0.05, log_or_amyloid_entry=0.7,
                                            log_or_amyloid_spread=0.7))
        # entry and spread are separate logistic steps in the generator; the spread step
        # (stages I..VI among brains past entry) is exactly a proportional-odds model.
        dd = d[d.braaksc >= 1]
        X = np.column_stack([dd.amyloid, (dd.age_death - 89) / 6, dd.apoe4, dd.msex, dd.educ])
        f = fit_ordinal(dd.braaksc, X, link="logit")
        b, se = f.beta[0], f.se()[0]
        covered += abs(b - 0.7) <= 1.96 * se
        errs.append(b - 0.7)
    assert abs(np.mean(errs)) < 0.03, np.mean(errs)
    assert 0.89 <= covered / reps <= 0.99, covered / reps


def test_entry_split_test_sees_amyloid_only_when_it_acts_on_cortical_entry():
    def h2(entry, seed):
        d = pd.DataFrame(simulate_radc_like(3000, seed=seed, p0=0.15, log_or_amyloid_entry=entry,
                                            log_or_amyloid_spread=0.7))
        d["map"] = 0.0
        return association(d, outcome="braaksc", exposure="amyloid", direction=1,
                           covariates=["age_death", "msex", "educ", "apoe4"], split=1)
    gated, part = h2(0.7, 11), h2(0.0, 12)
    assert abs(gated["log_or"] - 0.7) < 3 * gated["se"]
    assert abs(part["log_or"]) < 3 * part["se"]


def test_ordinal_mediation_recovers_true_natural_effects_and_rho_correction():
    for rho in (0.0, 0.4):
        d = simulate_ordinal_mediation(20000, beta2=0.5, gamma=0.4, beta3=0.3, rho_true=rho, sigma2=1.3, seed=21)
        f = fit_ordinal_mediation(d["T"], d["M"], d["Y"], d["X"])
        e = f.effects(1.0, 0.0, rho)
        assert e["nie"] == pytest.approx(d["true_nie"], abs=0.02), (rho, e, d["true_nie"])
        assert e["nde"] == pytest.approx(d["true_nde"], abs=0.03), (rho, e, d["true_nde"])
        assert f.effects(1.0, 0.0, f.rho_star)["nie"] == pytest.approx(0.0, abs=1e-12)
    naive = f.effects(1.0, 0.0, 0.0)["nie"]
    assert abs(naive - d["true_nie"]) > 0.1          # rho = 0.4 biases the naive NIE


def test_ordinal_mediation_bootstrap_coverage():
    covered = 0
    reps = 60
    for r in range(reps):
        d = simulate_ordinal_mediation(500, beta2=0.5, gamma=0.4, beta3=0.3, seed=900 + r, n_truth=100_000)
        b = bootstrap_ordinal_mediation(d["T"], d["M"], d["Y"], d["X"], t1=1.0, t0=0.0, n_boot=200, seed=r)
        lo, hi = np.quantile(b["nie"], [0.025, 0.975])
        covered += lo <= d["true_nie"] <= hi
    assert 0.85 <= covered / reps <= 1.0, covered / reps
