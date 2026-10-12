"""Module 6 validation. The gate (brief §6, Session 7): the spin test is
validated against a known null — independent but spatially smooth maps on a
synthetic sphere must give ~5% false positives under the spin test, while
naive region shuffling is shown to be inflated."""

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.stats import binom

from sim.network_sim import fibonacci_sphere, random_connectome, smooth_field
from src.network.diffusion import Diffuser, best_fit, crossover_tau, laplacian, seed_vector
from src.network.spin import corr, naive_permutations, perm_pvalue, random_rotation, spin_permutations

pytestmark = pytest.mark.validation

N_HEMI = 100
COORDS = fibonacci_sphere(N_HEMI)


def _maps(seed):
    return np.concatenate([smooth_field(COORDS, n_bumps=6, kappa=4.0, seed=seed),
                           smooth_field(COORDS, n_bumps=6, kappa=4.0, seed=seed + 50_000)])


@pytest.fixture(scope="module")
def spins():
    return spin_permutations(COORDS, COORDS, 500, seed=1, method="hungarian")


def test_spin_test_holds_its_false_positive_rate_on_a_known_null_and_naive_shuffling_does_not(spins):
    pairs, alpha = 200, 0.05
    naive = naive_permutations(2 * N_HEMI, 500, seed=2)
    fp_spin = fp_naive = 0
    for k in range(pairs):
        x, y = _maps(10 * k + 1), _maps(10 * k + 2)        # independent by construction
        r = corr(x, y)
        fp_spin += perm_pvalue(r, np.array([corr(x[p], y) for p in spins])) < alpha
        fp_naive += perm_pvalue(r, np.array([corr(x[p], y) for p in naive])) < alpha
    upper = binom.ppf(0.99, pairs, alpha) / pairs
    print(f"false-positive rate: spin {fp_spin / pairs:.3f}, naive {fp_naive / pairs:.3f} (limit {upper:.3f})")
    assert fp_spin / pairs <= upper
    assert fp_naive / pairs > 2 * alpha                 # the trap the brief warns about is real here


def test_spin_test_still_detects_a_real_association(spins):
    rng = np.random.default_rng(3)
    hits = 0
    for k in range(20):
        x = _maps(1000 + k)
        y = x + 0.5 * np.std(x) * rng.normal(size=x.size)
        hits += perm_pvalue(corr(x, y), np.array([corr(x[p], y) for p in spins])) < 0.05
    assert hits >= 18


def test_rotations_are_proper_and_spins_are_permutations(spins):
    R = random_rotation(np.random.default_rng(0))
    np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-12)
    assert np.linalg.det(R) == pytest.approx(1.0)
    for p in spins[:20]:
        assert sorted(p[:N_HEMI]) == list(range(N_HEMI))           # left stays left, no duplicates
        assert sorted(p[N_HEMI:]) == list(range(N_HEMI, 2 * N_HEMI))


def test_ring_laplacian_has_the_analytic_spectrum():
    n = 12
    W = np.zeros((n, n))
    for i in range(n):
        W[i, (i + 1) % n] = W[(i + 1) % n, i] = 1
    lam = np.sort(Diffuser.from_laplacian(laplacian(W, "unnormalized")).eigvals)
    np.testing.assert_allclose(lam, np.sort(2 - 2 * np.cos(2 * np.pi * np.arange(n) / n)), atol=1e-12)


def test_eigen_solution_matches_matrix_exponential_conserves_mass_and_equilibrates():
    W = random_connectome(fibonacci_sphere(30), scale=0.5, seed=4)
    L = laplacian(W, "unnormalized")
    d = Diffuser.from_laplacian(L)
    x0 = seed_vector(30, [3])
    np.testing.assert_allclose(d.evolve(x0, 0.7), expm(-0.7 * L) @ x0, atol=1e-12)
    assert d.evolve(x0, 5.0).sum() == pytest.approx(1.0, rel=1e-12)
    np.testing.assert_allclose(d.evolve(x0, 1e4), np.full(30, 1 / 30), atol=1e-9)


def test_best_fit_recovers_the_true_seed_and_time_and_seeds_merge_late():
    coords = fibonacci_sphere(60)
    d = Diffuser.from_laplacian(laplacian(random_connectome(coords, scale=0.4, seed=5), "unnormalized"))
    taus = np.linspace(0.01, 3, 300)
    ev = np.arange(60)
    observed = d.evolve(seed_vector(60, [7]), 0.8)
    r_true, tau_true, _ = best_fit(d, seed_vector(60, [7]), observed, ev, taus)
    r_other, _, _ = best_fit(d, seed_vector(60, [40]), observed, ev, taus)
    assert r_true == pytest.approx(1.0, abs=1e-9) and tau_true == pytest.approx(0.8, abs=0.011)
    assert r_other < r_true
    cross = crossover_tau(d, seed_vector(60, [7]), seed_vector(60, [40]), ev, np.linspace(0.01, 50, 2000), 0.95)
    assert cross is not None and cross > 0.8
