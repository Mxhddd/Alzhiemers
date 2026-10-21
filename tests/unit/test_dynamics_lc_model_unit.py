"""Unit tests of the LC rate-controller model (src/dynamics/lc_model.py) and
the Module 4 analysis helpers. The validation gate is in
tests/validation/test_dynamics*.py; these tests are stricter companions."""

import numpy as np
import pytest

from src.common.config import load_config
from src.common.paths import CONFIG_DIR
from src.dynamics import lc_model as L
from src.dynamics.analysis import Space, dose_grid, evaluate
from src.dynamics.model import ODEModel
from src.dynamics.params import ParameterSet
from src.dynamics.sensitivity import label_fraction
from src.dynamics.simulate import simulate

PS = ParameterSet.from_yaml(CONFIG_DIR / "dynamics_params.yaml")
CFG = load_config(CONFIG_DIR / "dynamics.yaml")
P0 = PS.values()
MODEL = L.build_model()


def test_analytic_jacobian_matches_finite_differences():
    fd = ODEModel("fd", MODEL.state_names, MODEL.rhs)
    for y in ([0.3, 0.1, 0.7, 2.0], [1.2, 0.02, 0.05, 30.0]):
        for p in (P0, L.drugs()["atomoxetine"].apply(P0, 2.0), L.drugs()["prazosin"].apply(P0, 0.5)):
            np.testing.assert_allclose(MODEL.jacobian(np.array(y), p), fd.jacobian(np.array(y), p), rtol=1e-5, atol=1e-8)


def test_steady_state_readout_equals_the_end_of_the_ode_from_the_healthy_start():
    for p in (P0, L.drugs()["atomoxetine"].apply(P0, 1.0), dict(P0, q=4.0, kappa=50.0, hA=0.0)):
        ss = L.healthy_steady_state(p)
        tr = simulate(MODEL, p, L.initial_state(p), (0, 2e5), np.array([0, 2e5]), method="Radau", rtol=1e-10, atol=1e-13)
        np.testing.assert_allclose(tr.y[:, -1], [ss[s] for s in L.STATES], rtol=1e-6)


def test_degenerate_case_from_a_nontrivial_start_matches_the_closed_form():
    p = dict(P0, **{k: 0.0 for k in L.coupling_params()})
    y0 = [0.3, 0.4, 0.6, 5.0]                     # R0 != 1: T uses the quadrature form
    t = np.linspace(0, 400, 41)
    tr = simulate(MODEL, p, y0, (0, 400), t, method="Radau", rtol=1e-11, atol=1e-13)
    np.testing.assert_allclose(tr.y, L.uncoupled_solution(t, p, y0), rtol=1e-7, atol=1e-10)


def test_fixed_points_and_stability_agree_with_the_4d_root_finder():
    from src.dynamics.fixedpoints import find_fixed_points, grid_guesses
    guesses = grid_guesses(CFG.get("fixed_points.bounds"), L.STATES, 4)
    for p in (P0, dict(P0, q=4.0, kappa=50.0, hA=0.0, r_rep=0.001, r_dam=0.5, K=0.03)):
        ana = L.fixed_points(p)
        num = find_fixed_points(MODEL, p, np.vstack([guesses] + [f["y"] for f in ana]))
        assert len(num) == len(ana)
        for a, n in zip(sorted(ana, key=lambda f: f["y"][0]), sorted(num, key=lambda f: f.y[0])):
            np.testing.assert_allclose(n.y, a["y"], rtol=1e-7)
            assert n.stable == a["stable"]


def test_q_at_most_two_is_monostable_and_prazosin_is_monotone():
    rng = np.random.default_rng(0)
    sp = Space.from_params(PS)
    for _ in range(50):
        p = sp.to_params({k: rng.uniform(*r) for k, r in sp.sampling_ranges().items()})
        p["q"] = rng.uniform(1, 2)
        assert len(L.roots_of_G(p)) == 1
        y = L.dose_response(p, "prazosin", dose_grid(CFG.get("dose_response.dose_grid")))
        dy = np.diff(y)
        assert np.all(dy <= 1e-12 * y[0]) or np.all(dy >= -1e-12 * y[0])


def test_readouts_do_not_depend_on_c0():
    a = evaluate(P0, CFG, with_dose_curves=False)
    b = evaluate(dict(P0, c0=0.1), CFG, with_dose_curves=False)
    for k in ("T0", "ratio_atomoxetine", "ratio_prazosin", "bistable_window_length"):
        assert a[k] == pytest.approx(b[k], rel=1e-12)


def test_space_scales_and_sensitivity_variants():
    sp = Space.from_params(PS)
    assert sp.sampling_ranges()["kappa"] == pytest.approx((-1.0, 2.0))
    assert sp.to_params({"kappa": 1.0})["kappa"] == pytest.approx(10.0)
    half = sp.variant(CFG.get("sensitivity_analyses.SA1_central_half"))
    assert half.ranges["kappa"] == pytest.approx((10 ** -0.25, 10 ** 1.25))
    assert half.ranges["u0"] == pytest.approx((0.3, 0.7))
    no_a1 = sp.variant(CFG.get("sensitivity_analyses.SA2_no_alpha1_clearance"))
    assert "eps" not in no_a1.ranges and no_a1.to_params({})["eps"] == 0.0
    lin = sp.variant(CFG.get("sensitivity_analyses.SA3_eta_linear"))
    assert lin.scales["eta"] == "linear" and lin.ranges["eta"] == (0.0, 10.0)


def test_labels_follow_the_registered_rule():
    assert label_fraction(2, 1000, fragile_upper=0.1, robust_lower=0.5)["label"] == "fragile"
    assert label_fraction(600, 1000, fragile_upper=0.1, robust_lower=0.5)["label"] == "robust"
    assert label_fraction(300, 1000, fragile_upper=0.1, robust_lower=0.5)["label"] == "intermediate"
    assert label_fraction(0, 1000, fragile_upper=0.1, robust_lower=0.5)["fraction"] == 0.0


def test_r_grid_in_code_is_the_registered_grid():
    g = CFG.get("fixed_points.r_grid")
    np.testing.assert_array_equal(L.r_grid(), L.r_grid(g["log10_min"], g["log10_switch"], g["n_log"], g["n_lin"]))


def test_dose_grid_is_zero_plus_41_log_spaced_doses():
    d = dose_grid(CFG.get("dose_response.dose_grid"))
    assert d[0] == 0 and len(d) == 42 and d[1] == pytest.approx(0.01) and d[-1] == pytest.approx(10.0)
