from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from src.mediation.mediation import estimate, killer2_verdict


def test_bad_inputs_raise():
    T = np.array([0.0, 1, 0, 1]); M = np.array([1.0, 2, 1, 2]); Y = np.array([1.0, 2, 3, 4])
    with pytest.raises(ValueError, match="non-finite"):
        estimate(T, M, np.array([1.0, np.nan, 3, 4]), None)
    with pytest.raises(ValueError, match="no variation"):
        estimate(np.zeros(4), M, Y, None)
    with pytest.raises(ValueError, match="rho"):
        estimate(T, M + np.array([0.1, -0.1, 0.2, 0.0]), Y, None).acme_at(np.array([1.0]))


@pytest.mark.parametrize("acme,lo,hi,expect", [
    (0.2, 0.1, 0.3, "NOT met"),
    (-0.2, -0.3, -0.1, "OPPOSITE"),
    (0.05, -0.02, 0.12, "includes 0"),
])
def test_killer2_verdicts(acme, lo, hi, expect):
    s = {"estimates": {"acme": {"estimate": acme, "ci_lower": lo, "ci_upper": hi}}}
    assert expect in killer2_verdict(s, predicted_sign=1)


# ---- ordinal models ----

def test_ordinal_fit_matches_statsmodels_for_both_links():
    from statsmodels.miscmodels.ordinal_model import OrderedModel

    from src.mediation.ordinal import fit_ordinal
    rng = np.random.default_rng(3)
    X = rng.normal(size=(800, 3))
    y = np.digitize(X @ np.array([0.5, -0.3, 0.0]) + 3 + rng.logistic(size=800), [0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
    for link in ("logit", "probit"):
        f = fit_ordinal(y, X, link=link)
        m = OrderedModel(y, X, distr=link).fit(method="bfgs", disp=False, maxiter=2000)
        assert np.allclose(f.beta, m.params[:3], atol=2e-4)
        assert np.allclose(f.se(), m.bse[:3], rtol=2e-3)
        assert f.loglik == pytest.approx(m.llf, abs=1e-4)


def test_weighted_fit_equals_duplicated_rows():
    from src.mediation.ordinal import fit_ordinal
    rng = np.random.default_rng(1)
    X = rng.normal(size=(300, 2))
    y = np.digitize(X[:, 0] + rng.logistic(size=300), [-1, 0, 1])
    w = rng.integers(1, 4, 300).astype(float)
    a = fit_ordinal(y, X, weights=w, hessian=False)
    b = fit_ordinal(np.repeat(y, w.astype(int)), np.repeat(X, w.astype(int), axis=0), hessian=False)
    assert np.allclose(a.beta, b.beta, atol=1e-5) and np.allclose(a.theta, b.theta, atol=1e-5)


def test_e_value_known_values():
    from src.mediation.ordinal import e_value
    assert e_value(1.0) == pytest.approx(1.0)
    assert e_value(4.0, common_outcome=False) == pytest.approx(4 + np.sqrt(12))
    assert e_value(4.0) == pytest.approx(2 + np.sqrt(2))       # sqrt(OR) for a common outcome
    assert e_value(0.25) == pytest.approx(e_value(4.0))


@pytest.mark.parametrize("b,se,expect", [
    (0.5, 0.1, "supports"),
    (-0.5, 0.1, "AGAINST"),
    (0.0, 0.02, "AGAINST (practically null)"),
    (0.05, 0.05, "inconclusive"),
    (0.0, 0.5, "uninformative"),
])
def test_preregistered_verdict_rule(b, se, expect):
    from src.mediation.hypotheses import verdict
    v = verdict(b, se, m=4, level=0.95, equivalence=0.10536, max_width=1.386)
    assert v["verdict"] == expect
    assert v["adj_level"] == pytest.approx(1 - 0.05 / 4)


# ---- RADC loader ----

def test_apoe4_carrier_coding():
    from src.mediation.radc import apoe4_carrier
    s = apoe4_carrier(pd.Series([22, 23, 24, 33, 34, 44, np.nan]), [24, 34, 44])
    assert s.tolist()[:6] == [0, 0, 1, 0, 1, 1] and np.isnan(s.iloc[6])


def test_exposure_summaries_earliest_mean_last_and_visit_age():
    from src.mediation.radc import exposure_summaries
    long = pd.DataFrame({"projid": [1, 1, 1, 2], "fu_year": [2, 0, 1, 0], "x": [3.0, np.nan, 1.0, 5.0],
                         "dcfdx": [2, 1, 1, 1]})
    s = exposure_summaries(long, id_col="projid", visit_col="fu_year", var="x", age_col=None,
                           age_bl=pd.Series({1: 80.0, 2: 70.0}), dx_col="dcfdx")
    assert s.loc[1, "x__earliest"] == 1.0 and s.loc[1, "x__visit"] == 1 and s.loc[1, "x__age"] == 81.0
    assert s.loc[1, "x__mean"] == 2.0 and s.loc[1, "x__last"] == 3.0 and s.loc[1, "x__dx"] == 1
    assert s.loc[2, "x__earliest"] == 5.0


def test_real_path_fails_loudly_and_names_the_provider():
    from src.common.config import load_config
    from src.common.io import MissingInputError
    from src.common.paths import CONFIG_DIR
    from src.mediation.radc import load_radc
    with pytest.raises(MissingInputError, match="requested, not received"):
        load_radc(load_config(CONFIG_DIR / "mediation.yaml"))


def _synthetic_cfg(monkeypatch, *, lc: bool, drop_col: str | None = None, missing_share: float = 0.0):
    import src.mediation.radc as radc
    from sim.mediation_sim import radc_like_files
    from src.common.config import Config
    from src.common.paths import CONFIG_DIR
    basic, long = radc_like_files(500, seed=2, log_or_amyloid_spread=0.7, log_or_expo=0.1,
                                  lc_name="lc_neuronal_density" if lc else None, lc_beta2=-0.5, lc_gamma=-0.4)
    if drop_col:
        basic = basic.drop(columns=[drop_col])
    if missing_share:
        rng = np.random.default_rng(0)
        basic.loc[rng.random(len(basic)) < missing_share, "educ"] = np.nan
    monkeypatch.setattr(radc, "_read", lambda p, provided_by: basic if "basic" in p else long)
    d = yaml.safe_load((CONFIG_DIR / "mediation.yaml").read_text())
    d["autopsy_weights"]["n_boot"] = 5
    d["mediation"]["n_boot"] = 20
    d["linear"]["n_boot"] = 20
    d["missing"]["mi_imputations"] = 2
    d["missing"]["mi_cycles"] = 2
    if lc:
        d["lc"]["variables"] = [{"name": "lc_neuronal_density", "damage_direction": -1}]
    return Config(d, Path("synthetic.yaml"))


class _Run:
    def __init__(self):
        self.drops, self.msgs = [], []

    def drop(self, **kw):
        self.drops.append(kw)

    def info(self, m):
        self.msgs.append(m)

    warning = info


def test_absent_codebook_column_stops_the_run(monkeypatch):
    from src.mediation.radc import CodebookError, load_radc
    cfg = _synthetic_cfg(monkeypatch, lc=False, drop_col="amylsqrt_est_8reg")
    with pytest.raises(CodebookError, match="amylsqrt_est_8reg"):
        load_radc(cfg)


@pytest.mark.filterwarnings("ignore")
def test_pipeline_on_synthetic_data_without_lc_measures(monkeypatch):
    from src.mediation.radc import load_radc
    from src.mediation.run_mediation import Analysis
    cfg = _synthetic_cfg(monkeypatch, lc=False, missing_share=0.1)
    people, info = load_radc(cfg)
    run = _Run()
    a = Analysis(cfg, run, people)
    a.amyloid_tests()
    a.social_tests()
    rows = pd.DataFrame(a.rows)
    prim = rows[rows.confirmatory]
    assert prim.test.tolist() == ["H1", "H2", "H3 soc_net", "H3 social_isolation", "H3 late_life_soc_act",
                                  "H3 vision"]
    assert set(rows.analysis) >= {"primary", "S1", "S2", "S3", "S4-mean", "S4-last", "S5", "S6", "S7", "S8",
                                  "S9", "S10-ROS", "S10-MAP", "descriptive"}
    assert run.drops, "incomplete participants must be counted and logged"
    # the MAP-only exposures drop the (constant) cohort indicator, and say so
    assert (rows[rows.test == "H3 social_isolation"].dropped_constant.str.contains("map")).all()
    assert a.lc_branch() is None
    assert any("cannot test mediation through the LC" in m for m in run.msgs)
    assert set(a.robustness().columns) == {"test", "verdict", "robust", "analyses_checked"}


@pytest.mark.filterwarnings("ignore")
def test_lc_branch_runs_only_when_configured_and_gives_a_killer2_reading(monkeypatch):
    from src.mediation.radc import load_radc
    from src.mediation.run_mediation import Analysis
    cfg = _synthetic_cfg(monkeypatch, lc=True)
    people, _ = load_radc(cfg)
    a = Analysis(cfg, _Run(), people)
    res = a.lc_branch()
    assert len(res) == 4 and "Killer 2" in res.attrs["overall"]
    assert {"nie", "nde", "rho_star", "nie_ciadj_lower", "S3_linear_acme"} <= set(res.columns)
