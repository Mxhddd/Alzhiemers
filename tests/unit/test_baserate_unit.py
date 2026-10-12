import numpy as np
import pytest

from src.baserate.data import load_braak_counts
from src.baserate.fit import log_likelihood
from src.baserate.data import StageCounts
from src.baserate.model import (AGated, BIndependent, ChainSpec, CMultiHit, InvalidParameters, Sequential,
                                build_generator, piecewise_constant)
from src.baserate.predict import age_band_folds

SPEC = ChainSpec(("0", "LC", "I", "II"), first_cortical=2)
HEADER = "age_lower,age_upper,stage,n,source\n"


def _csv(tmp_path, body):
    p = tmp_path / "braak.csv"
    p.write_text(HEADER + body)
    return p


def test_loader_reads_complete_groups(tmp_path):
    p = _csv(tmp_path, "40,50,0,5,S1\n40,50,LC,3,S1\n40,50,I,0,S1\n40,50,II,1,S1\n")
    g = load_braak_counts(p, list(SPEC.stages))
    assert len(g) == 1 and list(g[0].counts) == [5, 3, 0, 1] and g[0].midpoint == 45


@pytest.mark.parametrize("body,match", [
    ("40,50,0,5,S1\n40,50,LC,3,S1\n40,50,I,0,S1\n", "does not list stages"),          # II absent: not imputed as 0
    ("40,50,0,5,S1\n40,50,LC,3,S1\n40,50,I,0,S1\n40,50,III,1,S1\n", "not in stage_order"),
    ("40,50,0,5,S1\n40,50,0,5,S1\n", "duplicated"),
    ("40,50,0,5.5,S1\n", "non-negative integers"),
    ("50,40,0,5,S1\n", "age_upper <= age_lower"),
    ("40,50,0,,S1\n", "empty cells"),
])
def test_loader_refuses_bad_transcriptions(tmp_path, body, match):
    with pytest.raises(ValueError, match=match):
        load_braak_counts(_csv(tmp_path, body), list(SPEC.stages))


def test_loader_requires_exact_columns(tmp_path):
    p = tmp_path / "b.csv"
    p.write_text("age_lo,age_upper,stage,n,source\n")
    with pytest.raises(ValueError, match="columns must be exactly"):
        load_braak_counts(p, list(SPEC.stages))


def test_generator_rows_sum_to_zero_for_every_variant():
    am = piecewise_constant([(0, 200, 0.4)])
    rng = np.random.default_rng(0)
    for m in (Sequential(SPEC, 50), AGated(SPEC, am, 50), BIndependent(SPEC, 50), CMultiHit(SPEC, 3, 50)):
        th = rng.uniform(-6, -2, len(m.param_names))
        if "gamma" in m.param_names:
            th[m.param_names.index("gamma")] = 0.001
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


def test_negative_gated_rate_is_invalid_not_clipped():
    a = AGated(SPEC, piecewise_constant([(0, 200, 1.0)]))
    th = np.zeros(len(a.param_names))
    th[a.param_names.index("log_mu0")] = np.log(0.01)
    th[a.param_names.index("gamma")] = -0.02     # mu0 + gamma < 0
    with pytest.raises(InvalidParameters):
        build_generator(a, 60.0, th)
    data = [StageCounts(40, 50, "x", np.array([1, 1, 1, 1]))]
    assert log_likelihood(a, th, data, step=1.0, start_age=0) == -np.inf


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
    cfg = {"alpha_prior": {"mu": -5.0, "sigma": 3.0}, "beta_prior": {"mu": 0.0, "sigma": 1.0},
           "gamma_prior": {"mu": 0.0, "sigma": 1.0}, "draws": 100, "tune": 100, "chains": 2,
           "target_accept": 0.9, "ci_level": 0.95}
    idata, summ, diag = fit_bayesian(m, data, step=1.0, start_age=0, cfg_bayes=cfg, seed=1)
    assert set(diag) == {"max_rhat", "min_ess_bulk", "divergences"}
    assert abs(summ.loc["alpha_0", "mean"] - np.log(0.02)) < 0.1
    assert abs(summ.loc["beta_0", "mean"]) < 0.1   # per decade
