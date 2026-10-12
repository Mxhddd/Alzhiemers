"""Module 4 machinery validation on reference systems with analytic answers.
(The LC-model-specific checks the brief lists are in test_dynamics_lc_model.py
and need the researcher's model.)"""

import numpy as np
import pytest

from src.dynamics.bifurcation import bistable_range, folds, scan
from src.dynamics.fixedpoints import classify, find_fixed_points
from src.dynamics.model import ODEModel
from src.dynamics.sensitivity import is_u_shaped, lhs_robustness, sobol_indices
from src.dynamics.simulate import simulate, solver_independence

pytestmark = pytest.mark.validation


def test_logistic_fixed_points_and_eigenvalues():
    m = ODEModel("logistic", ("x",), lambda t, y, p: np.array([p["r"] * y[0] * (1 - y[0] / p["K"])]),
                 nonnegative=("x",))
    fps = find_fixed_points(m, {"r": 0.7, "K": 3.0}, np.linspace(-0.5, 5, 12)[:, None])
    assert [round(f.y[0], 8) for f in fps] == [0.0, 3.0]
    assert fps[0].eigenvalues[0].real == pytest.approx(0.7, rel=1e-5) and not fps[0].stable
    assert fps[1].eigenvalues[0].real == pytest.approx(-0.7, rel=1e-5) and fps[1].stable


def test_cusp_folds_match_the_analytic_value():
    # x' = h + r x - x^3: saddle-node folds at h = +/- 2 (r/3)^(3/2).
    m = ODEModel("cusp", ("x",), lambda t, y, p: np.array([p["h"] + p["r"] * y[0] - y[0] ** 3]))
    rows = scan(m, {"r": 1.0, "h": 0.0}, "h", np.linspace(-1, 1, 401), np.linspace(-2, 2, 9)[:, None])
    exact = 2 * (1 / 3) ** 1.5
    locs = sorted(f["location"] for f in folds(rows))
    res = folds(rows)[0]["resolution"]
    assert len(locs) == 2
    assert locs[0] == pytest.approx(-exact, abs=res) and locs[1] == pytest.approx(exact, abs=res)
    lo, hi = bistable_range(rows)
    assert -exact - res <= lo and hi <= exact + res


def test_classification_of_linear_systems():
    assert classify(np.array([-1.0, 2.0]))[1] == "saddle"
    assert classify(np.array([-1 + 2j, -1 - 2j])) == (True, "stable focus")
    assert classify(np.array([0.0, -1.0]))[0] is False


ROBERTSON = ODEModel(
    "robertson", ("y1", "y2", "y3"),
    lambda t, y, p: np.array([-0.04 * y[0] + 1e4 * y[1] * y[2],
                              0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] ** 2,
                              3e7 * y[1] ** 2]),
    nonnegative=("y1", "y2", "y3"), conserved={"total": lambda y: y.sum()})


def test_stiff_benchmark_conservation_nonnegativity_and_solver_independence():
    t = np.logspace(-5, 4, 40)
    # atol must resolve y2 (~1e-5 down to ~1e-10): the standard benchmark uses ~1e-14.
    simulate(ROBERTSON, {}, [1, 0, 0], (0, 1e4), t, method="Radau", rtol=1e-8, atol=1e-14)
    rep = solver_independence(ROBERTSON, {}, [1, 0, 0], (0, 1e4), t, rtol=1e-8, atol=1e-14)
    assert rep["converging"] and rep["cross_solver_tight"] < 1e-3, rep


def test_conservation_violation_is_caught():
    leaky = ODEModel("leak", ("a", "b"), lambda t, y, p: np.array([-y[0], 0.5 * y[0]]),
                     conserved={"total": lambda y: y.sum()})
    with pytest.raises(AssertionError, match="drifted"):
        simulate(leaky, {}, [1, 0], (0, 5), np.linspace(0, 5, 11), method="Radau", rtol=1e-8, atol=1e-10)


def test_zero_coupling_reduces_to_exponentials():
    # x' = -a x ; y' = -b y + c x. With c = 0: x = x0 e^{-at}, y = y0 e^{-bt}.
    m = ODEModel("lin", ("x", "y"), lambda t, y, p: np.array([-p["a"] * y[0], -p["b"] * y[1] + p["c"] * y[0]]))
    t = np.linspace(0, 10, 50)
    tr = simulate(m, {"a": 0.3, "b": 1.1, "c": 0.0}, [2.0, 5.0], (0, 10), t, method="BDF", rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(tr.y[0], 2.0 * np.exp(-0.3 * t), rtol=1e-7, atol=1e-10)
    np.testing.assert_allclose(tr.y[1], 5.0 * np.exp(-1.1 * t), rtol=1e-7, atol=1e-10)


def test_sobol_matches_ishigami_analytic_indices():
    a, b = 7.0, 0.1
    f = lambda q: np.sin(q["x1"]) + a * np.sin(q["x2"]) ** 2 + b * q["x3"] ** 4 * np.sin(q["x1"])  # noqa: E731
    r = {k: (-np.pi, np.pi) for k in ("x1", "x2", "x3")}
    out = sobol_indices(f, r, n_base=4096, seed=1)
    V = a**2 / 8 + b * np.pi**4 / 5 + b**2 * np.pi**8 / 18 + 0.5
    S1 = np.array([0.5 * (1 + b * np.pi**4 / 5) ** 2, a**2 / 8, 0.0]) / V
    ST = np.array([0.5 * (1 + b * np.pi**4 / 5) ** 2 + 8 * b**2 * np.pi**8 / 225, a**2 / 8,
                   8 * b**2 * np.pi**8 / 225]) / V
    np.testing.assert_allclose(out["S1"], S1, atol=0.03)
    np.testing.assert_allclose(out["ST"], ST, atol=0.03)


def test_u_shape_detector():
    d = np.linspace(0, 1, 11)
    assert is_u_shaped(d, (d - 0.5) ** 2 + 1, min_relative_rise=0.05)
    assert not is_u_shaped(d, 1 - d, min_relative_rise=0.05)                # monotone
    assert not is_u_shaped(d, 1 + 0.001 * (d - 0.5) ** 2, min_relative_rise=0.05)  # ripple, not U
    assert not is_u_shaped(d, (d - 0.9) ** 2 + 1, min_relative_rise=0.05)   # high end barely rises


def test_lhs_fraction_recovers_a_known_region():
    d = np.linspace(0, 1, 11)
    curve = lambda q: (d, (d - 0.5) ** 2 + 1 if q["theta"] > 0.3 else 1 - d)  # noqa: E731  U iff theta > 0.3
    out = lhs_robustness(curve, {"theta": (0.0, 1.0)}, n_samples=500, seed=2, min_relative_rise=0.05)
    assert out["ci95"][0] <= 0.7 <= out["ci95"][1]
    assert out["fraction"] == pytest.approx(0.7, abs=0.01)   # LHS stratifies: near-exact
