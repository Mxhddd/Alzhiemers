"""Module 3 validation: exact identities, recovery and coverage on simulated
cohorts (sim/mediation_sim.py), the rho sensitivity against a known
confounding correlation, and reproduction of a published analysis
(researcher-supplied)."""

import numpy as np
import pandas as pd
import pytest

from sim.mediation_sim import simulate_cohort
from src.common.config import load_config
from src.common.io import analysis_input
from src.common.paths import REPO_ROOT
from src.mediation.mediation import bootstrap, estimate, killer2_verdict, percentile_ci, summarise

pytestmark = pytest.mark.validation


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


def test_reproduces_published_mediation_analysis():
    t = load_config(REPO_ROOT / "tests/validation/targets/mediation_published.yaml")
    p, inp = t.get("published"), t.get("inputs")
    for k in ("citation", "acme", "acme_ci_lower", "acme_ci_upper", "ade"):
        if p.get(k) is None:
            raise AssertionError(f"mediation_published.yaml: researcher must transcribe published.{k}")
    for k in ("data", "exposure", "mediator", "outcome", "covariates"):
        if inp.get(k) is None:
            raise AssertionError(f"mediation_published.yaml: researcher must supply inputs.{k}")
    df = pd.read_csv(analysis_input(REPO_ROOT / inp["data"], provided_by="researcher (validation dataset)"))
    if df.isna().any().any():
        raise AssertionError("validation dataset has missing values; the published analysis's handling must be replicated explicitly")
    X = df[inp["covariates"]].to_numpy(float) if inp["covariates"] else None
    pe = estimate(df[inp["exposure"]].to_numpy(float), df[inp["mediator"]].to_numpy(float),
                  df[inp["outcome"]].to_numpy(float), X)
    assert p["acme_ci_lower"] <= pe.acme <= p["acme_ci_upper"], (pe.acme, p)
    assert pe.ade == pytest.approx(p["ade"], abs=0.01 * max(1.0, abs(p["ade"])))
