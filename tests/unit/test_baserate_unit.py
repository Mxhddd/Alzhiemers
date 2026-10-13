import numpy as np
import pytest

from src.baserate.data import load_braak_counts
from src.baserate.fit import log_likelihood
from src.baserate.data import StageCounts
from src.baserate.model import (AGated, BIndependent, ChainSpec, CMultiHit, InvalidParameters, Sequential,
                                build_generator, piecewise_constant)
from src.baserate.predict import age_band_folds

SPEC = ChainSpec(("0", "LC", "I", "II"), first_cortical=2)
HEADER = "age_lower,age_upper,stage,n,source_id,locator,verified\n"


def _csv(tmp_path, body, monkeypatch):
    import src.common.extracted as ex
    monkeypatch.setattr(ex, "analysis_input", lambda p, provided_by: p)
    p = tmp_path / "braak.csv"
    p.write_text(HEADER + body)
    return p


def _rows(*rows):
    return "".join(f"{r},S1,Table 1,false\n" for r in rows)


def test_loader_reads_complete_groups_and_converts_inclusive_bands(tmp_path, monkeypatch):
    p = _csv(tmp_path, _rows("41,50,0,5", "41,50,LC,3", "41,50,I,0", "41,50,II,1"), monkeypatch)
    g = load_braak_counts(p, list(SPEC.stages), upper_inclusive=True)
    assert len(g) == 1 and list(g[0].counts) == [5, 3, 0, 1]
    assert (g[0].age_lower, g[0].age_upper, g[0].midpoint) == (41, 51, 46)


@pytest.mark.parametrize("rows,match", [
    (("40,50,0,5", "40,50,LC,3", "40,50,I,0"), "does not list stages"),     # II absent: not imputed as 0
    (("40,50,0,5", "40,50,LC,3", "40,50,I,0", "40,50,III,1"), "not in stage_order"),
    (("40,50,0,5", "40,50,0,5"), "duplicated"),
    (("40,50,0,5.5",), "non-negative integers"),
    (("50,40,0,5",), "age_upper <= age_lower"),
])
def test_loader_refuses_bad_transcriptions(tmp_path, monkeypatch, rows, match):
    with pytest.raises(ValueError, match=match):
        load_braak_counts(_csv(tmp_path, _rows(*rows), monkeypatch), list(SPEC.stages), upper_inclusive=False)


def test_loader_refuses_blank_cells_and_missing_provenance(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="empty cells"):
        load_braak_counts(_csv(tmp_path, "40,50,0,,S1,T1,false\n", monkeypatch), list(SPEC.stages),
                          upper_inclusive=False)
    p = tmp_path / "b.csv"
    p.write_text("age_lower,age_upper,stage,n,verified\n40,50,0,1,false\n")
    with pytest.raises(ValueError, match="source"):
        load_braak_counts(p, list(SPEC.stages), upper_inclusive=False)


def test_band_average_points_are_uniform_within_the_band():
    g = StageCounts(40, 50, "x", np.array([1, 1]))
    a, w = g.eval_points("band_average", 5)
    np.testing.assert_allclose(a, [41, 43, 45, 47, 49]) and np.testing.assert_allclose(w.sum(), 1)


def test_generator_rows_sum_to_zero_for_every_variant():
    am = piecewise_constant([(0, 200, 0.4)])
    rng = np.random.default_rng(0)
    for m in (Sequential(SPEC, 50), AGated(SPEC, am, 50), BIndependent(SPEC, 50), CMultiHit(SPEC, 3, 50)):
        th = rng.uniform(-6, -2, len(m.param_names))

        Q = build_generator(m, 63.0, th)
        assert Q.shape == (m.n_hidden, m.n_hidden)
        np.testing.assert_allclose(Q.sum(axis=1), 0, atol=1e-15)


def test_variant_structures():
    m = BIndependent(SPEC)
    Q = build_generator(m, 50.0, np.zeros(len(m.param_names)))
    assert Q[0, 2] > 0 and Q[1, 2] > 0          # cortex reachable from every pre-cortical stage
    a = AGated(SPEC, piecewise_constant([(0, 200, 0.5)]))
    Qa = build_generator(a, 50.0, np.zeros(len(a.param_names)))
    assert Qa[0, 2] == 0                         # gated: cortex only via the LC stage
    c = CMultiHit(SPEC, hits=2)
    assert c.n_hidden == 2 * 2 + 2 and list(c.observed_of_hidden) == [0, 0, 1, 1, 2, 3]


def test_gated_rate_is_positive_for_any_log_hr_and_equals_mu0_times_hr_when_all_amyloid_positive():
    a = AGated(SPEC, piecewise_constant([(0, 200, 1.0)]))
    th = np.zeros(len(a.param_names))
    th[a.param_names.index("log_mu0")] = np.log(0.01)
    for log_hr in (-4.0, 0.0, 1.5):
        th[a.param_names.index("log_hr_amyloid")] = log_hr
        Q = build_generator(a, 60.0, th)
        assert Q[1, 2] == pytest.approx(0.01 * np.exp(log_hr), rel=1e-12)   # A = 1: rate = mu0 * HR
    half = AGated(SPEC, piecewise_constant([(0, 200, 0.5)]))
    th[a.param_names.index("log_hr_amyloid")] = np.log(3.0)
    assert build_generator(half, 60.0, th)[1, 2] == pytest.approx(0.01 * (1 + 2.0 * 0.5))   # mu0 + gamma A


def test_overflowing_rate_is_invalid_and_likelihood_is_minus_infinity():
    m = Sequential(SPEC)
    th = np.zeros(len(m.param_names))
    th[0], th[1] = 5.0, 50.0
    with pytest.raises(InvalidParameters):
        build_generator(m, 100.0, th)
    data = [StageCounts(90, 100, "x", np.array([1, 1, 1, 1]))]
    assert log_likelihood(m, th, data, step=1.0, start_age=0) == -np.inf


def test_amyloid_step_function_never_extrapolates():
    f = piecewise_constant([(40, 60, 0.1), (60, 80, 0.3)])
    assert f(59.9) == 0.1 and f(60) == 0.3
    with pytest.raises(ValueError, match="refusing to extrapolate"):
        f(85)
    with pytest.raises(ValueError, match="overlapping"):
        piecewise_constant([(40, 60, 0.1), (55, 80, 0.3)])


def test_folds_are_by_whole_age_band_and_span_the_range():
    data = [StageCounts(lo, lo + 10, s, np.array([1, 1])) for lo in range(40, 100, 10) for s in ("a", "b")]
    folds = age_band_folds(data, 3)
    assert folds[0] == folds[1]                  # two sources, same band, same fold
    assert sorted(set(folds)) == [0, 1, 2]
    with pytest.raises(ValueError, match="folds"):
        age_band_folds(data, 7)


@pytest.mark.slow
def test_bayesian_fit_runs_and_reports_diagnostics():
    """Small NUTS run on simulated two-state data: diagnostics are reported and
    the posterior mean is near the truth (lambda = 0.02 per year)."""
    from sim.baserate_sim import simulate_counts
    from src.baserate.fit import fit_bayesian

    m = Sequential(ChainSpec(("0", "1"), first_cortical=1), 50, 10)
    truth = np.array([np.log(0.02), 0.0])
    data = simulate_counts(m, truth, [(a, a + 10) for a in range(0, 100, 10)], 500, step=1.0, start_age=0, seed=3)
    cfg = {"log_rate_prior": {"mu": -5.0, "sigma": 3.0}, "slope_prior": {"mu": 0.0, "sigma": 1.0},
           "log_hr_prior": {"mu": 0.0, "sigma": 1.0}, "draws": 200, "tune": 200, "chains": 2,
           "init_overdispersion": 2.0,
           "target_accept": 0.9, "ci_level": 0.95}
    idata, summ, diag = fit_bayesian(m, data, step=1.0, start_age=0, cfg_bayes=cfg, seed=1,
                                     init_theta=np.array([np.log(0.02), 0.0]))
    assert {"max_rhat", "min_ess_bulk", "divergences", "mean_leapfrog_steps"} <= set(diag)
    assert diag["max_rhat"] < 1.05 and diag["divergences"] == 0
    assert abs(summ.loc["alpha_0", "mean"] - np.log(0.02)) < 0.1
    assert abs(summ.loc["beta_0", "mean"]) < 0.1   # per decade


@pytest.mark.parametrize("mode", ["midpoint", "band_average"])
@pytest.mark.parametrize("which", ["sequential", "A", "B", "C"])
def test_jax_likelihood_and_gradient_agree_with_the_numpy_reference(which, mode):
    from src.baserate.jaxlik import make_loglik

    am = piecewise_constant([(0, 40, 0.0), (40, 60, 0.2), (60, 200, 0.6)])
    m = {"sequential": Sequential(SPEC, 50, 10), "A": AGated(SPEC, am, 50, 10),
         "B": BIndependent(SPEC, 50, 10), "C": CMultiHit(SPEC, 2, 50, 10)}[which]
    data = [StageCounts(lo, lo + 10, "x", np.array([30, 10, 5, 2])) for lo in range(30, 100, 10)]
    rng = np.random.default_rng(1)
    ll_j, vg = make_loglik(m, data, step=0.5, start_age=0, eval_mode=mode, n_points=4)
    for _ in range(3):
        th = rng.uniform(-4, -2, len(m.param_names))
        for i, n in enumerate(m.param_names):
            if n.startswith("beta") or n == "b_nu":
                th[i] = rng.uniform(-0.3, 0.6)
            if n == "log_hr_amyloid":
                th[i] = rng.uniform(-1, 1)
        ref = log_likelihood(m, th, data, step=0.5, start_age=0, eval_mode=mode, n_points=4)
        val, grad = vg(th)
        assert float(val) == pytest.approx(ref, rel=1e-9, abs=1e-9)
        h = 1e-6
        fd = np.array([(log_likelihood(m, th + h * e, data, step=0.5, start_age=0, eval_mode=mode, n_points=4)
                        - log_likelihood(m, th - h * e, data, step=0.5, start_age=0, eval_mode=mode, n_points=4))
                       / (2 * h)
                       for e in np.eye(len(th))])
        np.testing.assert_allclose(np.asarray(grad), fd, rtol=1e-4, atol=1e-4)


def test_grid_puts_breakpoints_and_ages_on_knots_and_never_exceeds_the_step():
    from src.baserate.solve import integration_grid

    g = integration_grid([45.0, 55.0], step=0.3, start_age=0, breakpoints=[40.0, 50.0, 120.0])
    edges = np.concatenate([[0], np.cumsum(g.h)])
    for knot in (40.0, 45.0, 50.0, 55.0):
        assert np.min(np.abs(edges - knot)) < 1e-9
    assert g.h.max() <= 0.3 + 1e-12 and edges[-1] == pytest.approx(55.0)
    assert list(g.record_after) == [int(np.argmin(np.abs(edges - 45.0))), len(g.h)]


def test_error_bound_matches_the_closed_form_for_a_sequential_chain():
    from src.baserate.solve import error_bound

    m = Sequential(SPEC, 50, 10)
    th = np.array([-2.0, 0.5, -2.5, -0.3, -3.0, 1.0])
    got = error_bound(m, th, start_age=0, end_age=100, dt=0.05)
    t = np.arange(0.025, 100, 0.05)
    a, b = th[0::2], th[1::2] / 10.0                     # per-year slopes
    r = np.exp(a[None, :] + (th[1::2][None, :]) * ((t[:, None] - 50) / 10))
    rows = np.column_stack([r[:, :-1] * b[:-1] ** 2 / 12 + r[:, :-1] * r[:, 1:] * np.abs(b[:-1] - b[1:]) / 6,
                            r[:, -1] * b[-1] ** 2 / 12])
    np.testing.assert_allclose(got, rows.max(axis=1).sum() * 0.05, rtol=1e-3)


def test_jax_likelihood_is_finite_at_the_corner_of_the_search_box():
    from src.baserate.jaxlik import make_loglik

    m = BIndependent(SPEC, 50, 10)
    data = [StageCounts(91, 101, "x", np.array([1, 1, 1, 1]))]
    ll, _ = make_loglik(m, data, step=0.25, start_age=1)
    x = np.full(len(m.param_names), 2.0)                # log-rates and slopes at the upper box edge
    assert np.isfinite(float(ll(x)))


def test_amyloid_splice_uses_the_inner_curve_inside_and_the_primary_outside():
    from src.baserate.run_baserate import splice_amyloid

    primary = piecewise_constant([(1, 41, 0.0), (41, 51, 0.1), (51, 91, 0.4), (91, 101, 0.7)])
    s = splice_amyloid(primary, [(47.5, 52.5, 0.2), (52.5, 92.5, 0.5)])
    assert [s(t) for t in (30, 45, 50, 60, 92, 95)] == [0.0, 0.1, 0.2, 0.5, 0.5, 0.7]
    assert s.breakpoints[0] == 1 and s.breakpoints[-1] == 101


def test_preregistered_three_way_readings():
    import pandas as pd
    from src.baserate.run_baserate import cv_verdict_3way, gating_verdict

    assert gating_verdict((0.1, 1.2), 0.322).startswith("supports")
    assert gating_verdict((-0.9, 0.2), 0.322).startswith("AGAINST")
    assert gating_verdict((-0.07, 1.06), 0.322).startswith("inconclusive")   # the dry-run case: not "against"
    cv = lambda a, b, c: pd.DataFrame({"model": ["A_gated", "B_independent", "C_multihit"],  # noqa: E731
                                       "cv_heldout_loglik": [a, b, c]})
    assert cv_verdict_3way(cv(-100, -105, -110), 2).startswith("supports")
    assert cv_verdict_3way(cv(-100, -95, -110), 2).startswith("AGAINST")
    assert cv_verdict_3way(cv(-100, -101, -110), 2).startswith("inconclusive")
