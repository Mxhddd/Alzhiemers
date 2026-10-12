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
