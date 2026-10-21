import numpy as np
import pytest

from src.network.diffusion import Diffuser, laplacian, seed_vector
from src.network.spin import perm_pvalue, spin_permutations


def test_connectome_validation():
    with pytest.raises(ValueError, match="symmetric"):
        laplacian(np.array([[0, 1], [0, 0.0]]), "normalized")
    with pytest.raises(ValueError, match="isolated"):
        laplacian(np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0.0]]), "unnormalized")
    with pytest.raises(ValueError, match="non-negative"):
        laplacian(np.array([[0, -1], [-1, 0.0]]), "normalized")
    with pytest.raises(ValueError, match="kind"):
        laplacian(np.array([[0, 1], [1, 0.0]]), "random_walk")


def test_normalized_laplacian_spectrum_lies_in_0_2():
    rng = np.random.default_rng(0)
    W = rng.uniform(size=(8, 8)); W = (W + W.T) / 2; np.fill_diagonal(W, 0)
    lam = Diffuser.from_laplacian(laplacian(W, "normalized")).eigvals
    assert lam.min() == pytest.approx(0, abs=1e-10) and lam.max() <= 2 + 1e-10


def test_seed_vector_and_pvalue():
    assert seed_vector(4, [1, 2]).tolist() == [0, 0.5, 0.5, 0]
    assert perm_pvalue(0.9, np.zeros(99)) == pytest.approx(0.01)   # never exactly 0


def test_spin_requires_valid_method_and_coordinates():
    c = np.eye(3)
    with pytest.raises(ValueError, match="method"):
        spin_permutations(c, c, 2, seed=0, method="shuffle")
    with pytest.raises(ValueError, match="zero-length"):
        spin_permutations(np.zeros((3, 3)), c, 2, seed=0, method="hungarian")


# --- Vectorised fitting, verdict, connectome assembly, analysis, ADNI step ---

import pandas as pd

from sim.network_sim import fibonacci_sphere, random_connectome, smooth_field
from src.network.adni_tau import attach_nearest, build_group_maps, cortical_regions, suvr_columns, w_scores
from src.network.analysis import analyse, seed_ranking
from src.network.connectome import group_average, uniform_projection
from src.network.diffusion import (best_fit, corr_curves, crossover_from_curves, crossover_tau, max_over_tau,
                                   prediction_curves, standardise_rows)
from src.network import verdict as V


def _world(n_hemi=20, seed=0):
    """Two hemispheres on spheres, a distance-dependent cortex, 2 'LC' nodes
    wired to every cortical node of their side with random weights."""
    c = fibonacci_sphere(n_hemi)
    xyz = np.vstack([c + [-1.5, 0, 0], c + [1.5, 0, 0]])
    Wc = random_connectome(xyz, scale=0.6, seed=seed)
    n = 2 * n_hemi + 2
    W = np.zeros((n, n))
    W[:2 * n_hemi, :2 * n_hemi] = Wc
    rng = np.random.default_rng(seed + 1)
    for k, side in enumerate((range(0, n_hemi), range(n_hemi, 2 * n_hemi))):
        w = rng.lognormal(0, 1, n_hemi) * Wc.sum() / Wc.shape[0] / n_hemi
        W[2 * n_hemi + k, list(side)] = w
        W[list(side), 2 * n_hemi + k] = w
    names = [f"ctx-lh-r{i}" for i in range(n_hemi)] + [f"ctx-rh-r{i}" for i in range(n_hemi)] + ["LC-lh", "LC-rh"]
    return W, names, c


def test_vectorised_curves_equal_best_fit_and_crossover():
    W, names, _ = _world()
    d = Diffuser.from_laplacian(laplacian(W, "normalized"))
    ev = np.arange(40)
    taus = np.geomspace(0.01, 50, 120)
    obs = np.random.default_rng(1).normal(size=40)
    x_lc, x_ec = seed_vector(42, [40, 41]), seed_vector(42, [3, 23])
    r, t, curve = best_fit(d, x_lc, obs, ev, taus)
    P = prediction_curves(d, x_lc, ev, taus)
    R = corr_curves(standardise_rows(P, "pearson"), standardise_rows(obs, "pearson"))
    np.testing.assert_allclose(R[0], curve, atol=1e-10)
    rb, tb = max_over_tau(R, taus)
    assert rb[0] == pytest.approx(r, abs=1e-10) and tb[0] == t
    Pe = prediction_curves(d, x_ec, ev, taus)
    assert crossover_from_curves(P, Pe, taus, 0.95) == crossover_tau(d, x_lc, x_ec, ev, taus, 0.95)


def test_spearman_is_invariant_to_monotone_transforms_and_constant_rows_are_nan():
    P = np.array([[1.0, 2.0, 3.0, 5.0], [2.0, 2.0, 2.0, 2.0]])
    Z = standardise_rows(P, "spearman")
    assert np.isnan(Z[1]).all()
    a = corr_curves(Z, standardise_rows(np.array([0.1, 0.5, 0.2, 9.0]), "spearman"))
    b = corr_curves(Z, standardise_rows(np.exp(np.array([0.1, 0.5, 0.2, 9.0])), "spearman"))
    np.testing.assert_allclose(a[:, 0], b[:, 0])
    with pytest.raises(ValueError, match="correlation method"):
        standardise_rows(P, "kendall")


def _cmp(**kw):
    base = dict(r_lc=0.6, r_ec=0.4, tau_lc=0.5, tau_ec=0.5, p_lc=0.01, p_lc_better=0.01, p_ec_better=0.99,
                crossover=5.0)
    base.update(kw)
    return V.SeedComparison(**base)


def test_verdict_rules_in_order():
    a = 0.05
    assert V.verdict(_cmp(), a) == V.SUPPORTS
    assert V.verdict(_cmp(tau_lc=6, tau_ec=7), a) == V.UNINFORMATIVE
    assert V.verdict(_cmp(r_lc=0.4, r_ec=0.4, p_ec_better=0.5), a) == V.AGAINST       # a tie is against
    assert V.verdict(_cmp(r_lc=0.2, p_ec_better=0.01), a) == V.AGAINST_BEYOND_NULL
    assert V.verdict(_cmp(p_lc_better=0.2), a) == V.INCONCLUSIVE
    assert V.verdict(_cmp(p_lc=0.2), a) == V.INCONCLUSIVE
    assert V.verdict(_cmp(tau_lc=6), a) == V.INCONCLUSIVE                            # LC fits only past crossover
    assert V.verdict(_cmp(crossover=None, tau_lc=40), a) == V.SUPPORTS


def test_spin_p_values_and_robustness_rule():
    null = np.arange(99.0)
    assert V.spin_p_upper(98.5, null) == pytest.approx(1 / 100)
    assert V.spin_p_lower(-1, null) == pytest.approx(1 / 100)
    sens = {"S1": V.SUPPORTS, "S2": V.SUPPORTS, "S4": V.AGAINST}
    assert V.robust(V.SUPPORTS, sens, ["S1", "S2"])
    assert not V.robust(V.SUPPORTS, sens, ["S1", "S4"])
    assert not V.robust(V.SUPPORTS, sens, ["S1", "S9"])                              # missing -> not robust
    assert V.robust(V.AGAINST, {"S1": V.AGAINST_BEYOND_NULL}, ["S1"])


def test_group_average_threshold_and_uniform_projection():
    A = np.array([[0, 2, 0], [2, 0, 1], [0, 1, 0.0]])
    B = np.array([[0, 4, 3], [4, 0, 0], [3, 0, 0.0]])
    np.testing.assert_allclose(group_average([A, B], None), [[0, 3, 1.5], [3, 0, 0.5], [1.5, 0.5, 0]])
    np.testing.assert_allclose(group_average([A, B], 1.0), [[0, 3, 0], [3, 0, 0], [0, 0, 0]])
    W = np.array([[0, 1, 2, 3], [1, 0, 0, 0], [2, 0, 0, 1], [3, 0, 1, 0.0]])
    U = uniform_projection(W, 0, [1, 2])
    assert U[0].sum() == pytest.approx(W[0].sum()) and U[0, 3] == 0 and U[0, 1] == U[0, 2] == 3
    np.testing.assert_allclose(U, U.T)
    assert U[2, 3] == 1                                                           # other edges unchanged


def _analysis_inputs(true_seed, noise, rep):
    W, names, c = _world(seed=rep)
    perms = spin_permutations(c, c, 200, seed=7, method="hungarian")
    d = Diffuser.from_laplacian(laplacian(W, "normalized"))
    x0 = seed_vector(42, [40, 41] if true_seed == "LC" else [3, 23])
    sig = d.evolve(x0, 1.0)[:40]
    sm = np.concatenate([smooth_field(c, n_bumps=4, kappa=4.0, seed=100 + rep),
                         smooth_field(c, n_bumps=4, kappa=4.0, seed=200 + rep)])
    obs = sig / sig.std() + noise * sm / sm.std()
    seeds = {"LC": ["LC-lh", "LC-rh"], "entorhinal": ["ctx-lh-r3", "ctx-rh-r3"]}
    return W, names, obs, perms, seeds


@pytest.mark.parametrize("true_seed, expected", [("LC", "supports"), ("EC", "against")])
def test_analysis_recovers_the_true_seed(true_seed, expected):
    W, names, obs, perms, seeds = _analysis_inputs(true_seed, 0.3, 0)
    res = analyse(W, names, np.arange(40), obs, perms, seeds, ["ctx-lh-r3", "ctx-rh-r3"],
                  taus=np.geomspace(0.01, 50, 100), laplacian_kind="normalized", method="pearson",
                  crossover_threshold=0.95, alpha=0.05)
    assert V.direction(res.verdict) == expected
    assert res.fits["LC"].n_eval == res.fits["entorhinal"].n_eval == 38


def test_analysis_rejects_mismatched_inputs_and_ranking_is_sorted():
    W, names, obs, perms, seeds = _analysis_inputs("LC", 0.3, 0)
    with pytest.raises(ValueError, match="same length"):
        analyse(W, names, np.arange(39), obs, perms, seeds, [], taus=np.geomspace(0.01, 50, 20),
                laplacian_kind="normalized", method="pearson", crossover_threshold=0.95, alpha=0.05)
    rk = seed_ranking(W, names, np.arange(40), obs, {"LC": seeds["LC"], "r5": ["ctx-lh-r5", "ctx-rh-r5"]}, [],
                      taus=np.geomspace(0.01, 50, 20), laplacian_kind="normalized", method="pearson")
    assert [r for _, r in rk] == sorted([r for _, r in rk], reverse=True)


def test_adni_column_names_and_nearest_match():
    cols = suvr_columns("CTX_{hemi}_{region}_SUVR")
    assert cols["ctx-lh-entorhinal"] == "CTX_LH_ENTORHINAL_SUVR" and len(cortical_regions()) == 68
    base = pd.DataFrame({"RID": [1, 2], "SCANDATE": pd.to_datetime(["2020-01-01", "2020-01-01"])})
    other = pd.DataFrame({"RID": [1, 1, 2], "D": ["2019-06-01", "2020-03-01", "2018-01-01"], "v": [0, 1, 1]})
    out = attach_nearest(base, other, id_col="RID", base_date="SCANDATE", other_date="D", value_col="v",
                         out_col="amy", max_days=365)
    assert out["amy"].tolist()[0] == 1 and np.isnan(out["amy"].tolist()[1])     # nearest; outside window


def test_w_scores_are_standardised_residuals():
    rng = np.random.default_rng(0)
    ref = pd.DataFrame({"age": rng.uniform(60, 80, 50), "sex": rng.integers(0, 2, 50).astype(float)})
    ref["x"] = 1 + 0.01 * ref["age"] + 0.1 * ref["sex"] + rng.normal(0, 0.05, 50)
    W = w_scores(ref, ref, ["x"], ["age", "sex"])
    assert W["x"].mean() == pytest.approx(0, abs=1e-10)
    assert (W["x"] ** 2).sum() / (50 - 3) == pytest.approx(1.0)


def _adni_tables(n=120, seed=0):
    rng = np.random.default_rng(seed)
    cols = suvr_columns("CTX_{hemi}_{region}_SUVR")
    tau = pd.DataFrame({"RID": np.arange(n), "SCANDATE": "2020-06-01", "TRACER": "FTP", "qc_flag": 2})
    for c in cols.values():
        tau[c] = 1 + rng.normal(0, 0.05, n)
    tau.loc[0, "TRACER"] = "MK6240"
    tau.loc[1, "qc_flag"] = 0
    amy = pd.DataFrame({"RID": np.arange(n), "SCANDATE": "2020-03-01", "AMYLOID_STATUS": np.arange(n) % 2})
    amy.loc[2, "SCANDATE"] = "2015-01-01"                                         # outside the window
    dx = pd.DataFrame({"RID": np.arange(n), "EXAMDATE": "2020-06-15", "DIAGNOSIS": np.where(np.arange(n) % 4 < 3, 1, 2)})
    dem = pd.DataFrame({"RID": np.arange(n), "PTDOBYY": rng.integers(1940, 1960, n), "PTGENDER": rng.integers(1, 3, n)})
    pos_cu = (np.arange(n) % 2 == 1) & (np.arange(n) % 4 < 3)
    tau.loc[pos_cu, cols["ctx-lh-entorhinal"]] += 0.3                             # a planted effect
    return tau, amy, dx, dem


def test_adni_group_maps_count_exclusions_and_find_a_planted_effect():
    import yaml
    a = yaml.safe_load(open("config/network.yaml"))["adni"]
    a = {**a, "min_group_n": 10}
    maps, counts = build_group_maps(*_adni_tables(), a)
    assert counts["tau_rows_other_tracer"] == 1 and counts["tau_rows_failed_qc"] == 1
    assert counts["excluded_no_amyloid_scan_in_window"] == 1
    m = maps["apos_cu"]
    assert m.idxmax() == "ctx-lh-entorhinal" and m["ctx-lh-entorhinal"] > 3
    assert "apos_mci" in maps
    with pytest.raises(ValueError, match="reference group"):
        build_group_maps(*_adni_tables(), {**a, "min_group_n": 1000})
