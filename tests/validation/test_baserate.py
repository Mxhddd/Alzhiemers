"""Module 2 validation (PROJECT_BRIEF.md §5): (a) simulation recovery at three
sample sizes, (b) analytic two-state check to 8 decimal places, (c) stage
probabilities sum to 1 at every age. (d), the falsification test against
published incidence, needs researcher data and is part of the analysis run.

The step sizes here are TEST settings for the machinery, not the analysis
step, which the researcher supplies in config/baserate.yaml.
"""

import numpy as np
import pytest

from sim.baserate_sim import simulate_counts
from src.baserate.fit import fit_mle
from src.baserate.model import AGated, BIndependent, ChainSpec, CMultiHit, Sequential, piecewise_constant
from src.baserate.predict import predict_incidence
from src.baserate.solve import max_step_error, stage_probabilities

pytestmark = pytest.mark.validation

BANDS = [(a, a + 10) for a in range(0, 100, 10)]
CFG_MLE = {"n_starts": 6, "method": "l-bfgs-b", "maxiter": 2000, "newton_steps": 20, "max_newton_decrement": 1e-3,
           "init": {"log_rate": [-8.0, -2.0], "slope": [-0.5, 1.0], "log_hr": [-0.5, 0.5]},
           "bounds": {"log_rate": [-20.0, 3.0], "slope": [-3.0, 3.0], "log_hr": [-5.0, 5.0]}}
C, S = 50.0, 10.0  # age centre and scale used throughout these tests
SPEC3 = ChainSpec(("none", "LC", "cortex"), first_cortical=2)


def test_b_two_state_constant_rate_matches_analytic_to_8dp():
    lam = 0.037
    m = Sequential(ChainSpec(("0", "1"), first_cortical=1))
    ages = np.linspace(0, 120, 97)
    for step in (0.5, 3.7, 25.0):   # Q is constant, so every step is exact
        p = stage_probabilities(m, np.array([np.log(lam), 0.0]), ages, step=step, start_age=0)
        np.testing.assert_allclose(p[:, 1], 1 - np.exp(-lam * ages), atol=5e-9, rtol=0)
    inc = predict_incidence(m, np.array([np.log(lam), 0.0]), [20, 60], threshold=1, step=1.0, start_age=0)
    np.testing.assert_allclose(inc, lam, rtol=1e-12)


@pytest.mark.parametrize("model,theta", [
    (Sequential(SPEC3, C, S), [-4, 0.3, -7, 0.5]),
    (AGated(SPEC3, piecewise_constant([(0, 50, 0.0), (50, 70, 0.3), (70, 200, 0.6)]), C, S), [-4, 0.3, -6, 1.2]),
    (BIndependent(SPEC3, C, S), [-4, 0.3, -7, 0.5]),
    (CMultiHit(SPEC3, hits=3, age_centre=C, age_scale=S), [-4, 0.3, -3.5]),
])
def test_c_probabilities_sum_to_one_at_every_age(model, theta):
    ages = np.arange(0, 110.01, 0.5)
    p = stage_probabilities(model, np.array(theta, float), ages, step=1.0, start_age=0)
    np.testing.assert_allclose(p.sum(axis=1), 1.0, atol=1e-12)
    assert (p >= 0).all()


def _recover(model, theta, n, seed):
    data = simulate_counts(model, theta, BANDS, n, step=1.0, start_age=0, seed=seed)
    return fit_mle(model, data, step=1.0, start_age=0, cfg_mle=CFG_MLE, seed=seed)


def test_a_simulation_recovery_at_three_sample_sizes():
    model = Sequential(SPEC3, C, S)
    theta = np.array([-4.0, 0.3, -7.0, 0.5])
    errors = []
    for n in (100, 400, 1600):
        fit = _recover(model, theta, n, seed=1000 + n)
        assert fit.hessian_ok, f"n={n}: Hessian not positive definite"
        z = (fit.theta - theta) / fit.se
        assert np.all(np.abs(z) < 3.5), f"n={n}: estimates {fit.theta} vs truth {theta}, z = {z}"
        errors.append(np.sqrt(np.mean(((fit.theta - theta) / np.abs(theta)) ** 2)))
    assert errors[2] < errors[0], f"relative error did not shrink with n: {errors}"


def test_a_amyloid_hazard_ratio_is_recovered_in_the_gated_variant():
    amyloid = piecewise_constant([(0, 50, 0.0), (50, 70, 0.3), (70, 200, 0.6)])
    model = AGated(SPEC3, amyloid, C, S)
    theta = np.array([-4.0, 0.3, np.log(0.002), np.log(4.0)])
    fit = _recover(model, theta, 3000, seed=7)
    i = model.param_names.index("log_hr_amyloid")
    assert fit.hessian_ok and not fit.on_bound
    assert abs(fit.theta[i] - theta[i]) < 3 * fit.se[i], (fit.theta[i], fit.se[i])


def test_null_amyloid_effect_is_not_manufactured():
    """With amyloid having NO effect in the truth (HR = 1), the fit must not
    report one: the analysis can come out against the hypothesis."""
    amyloid = piecewise_constant([(0, 50, 0.0), (50, 70, 0.3), (70, 200, 0.6)])
    model = AGated(SPEC3, amyloid, C, S)
    theta = np.array([-4.0, 0.3, np.log(0.01), 0.0])
    fit = _recover(model, theta, 3000, seed=8)
    i = model.param_names.index("log_hr_amyloid")
    assert abs(fit.theta[i]) < 3 * fit.se[i]


def test_age_acceleration_is_not_read_as_amyloid_gating():
    """No amyloid effect in the truth, but the LC -> cortex rate rises with age
    (0.3 per decade). The fit must not report an amyloid effect. Added
    2026-10-13: the brief's variant A has no age term on the gated edge, so
    it attributes any age acceleration to the age-rising amyloid curve and
    this test is EXPECTED TO FAIL until the open issue in
    manuscript/notes/module2_identifiability.md is resolved. It keeps the
    validation gate shut on that issue."""
    amyloid = piecewise_constant([(0, 50, 0.0), (50, 70, 0.3), (70, 200, 0.6)])
    truth = Sequential(SPEC3, C, S)
    theta = np.array([-4.0, 0.3, np.log(0.01), 0.3])
    data = simulate_counts(truth, theta, BANDS, 3000, step=1.0, start_age=0, seed=8)
    model = AGated(SPEC3, amyloid, C, S)
    fit = fit_mle(model, data, step=1.0, start_age=0, cfg_mle=CFG_MLE, seed=8)
    i = model.param_names.index("log_hr_amyloid")
    assert abs(fit.theta[i]) < 3 * fit.se[i], (
        f"amyloid log HR {fit.theta[i]:.2f} (SE {fit.se[i]:.2f}) manufactured from age acceleration alone")


def test_multihit_incidence_scales_as_t_to_the_m_minus_1():
    m = CMultiHit(ChainSpec(("0", "1"), first_cortical=1), hits=3)
    ages = np.array([0.5, 1.0])
    inc = predict_incidence(m, np.array([np.log(0.01)]), ages, threshold=1, step=0.01, start_age=0)
    slope = np.log(inc[1] / inc[0]) / np.log(ages[1] / ages[0])
    assert slope == pytest.approx(2.0, abs=0.01)


def test_step_error_shrinks_with_step_like_a_second_order_method():
    model = Sequential(SPEC3, C, S)
    theta = np.array([-4.0, 0.6, -7.0, 0.8])
    ages = np.arange(5, 100, 10.0)
    e = [max_step_error(model, theta, ages, step=h, start_age=0) for h in (4.0, 2.0, 1.0)]
    assert e[0] > e[1] > e[2]
    assert 3.0 < e[1] / e[2] < 5.0, e   # midpoint rule: error ~ h^2
