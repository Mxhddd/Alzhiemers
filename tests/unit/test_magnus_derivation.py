"""Numerical check of manuscript/derivations/module2_step_size.md, §2-3:
the exponential-midpoint local error equals (h^3/24) Q'' + (h^3/12)[Q, Q']
to leading order, and the sequential-chain commutator has the closed form
2 r_k r_{k+1} |b_k - b_{k+1}| per row. (Checks algebra; not a substitute for
the researcher's hand check.)"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

rng = np.random.default_rng(0)
K = 5
A_ = rng.uniform(-3, -1.5, K - 1)
B_ = rng.uniform(-0.15, 0.15, K - 1)


def _Q(t, n=0):
    r = np.exp(A_ + B_ * t) * B_**n
    M = np.zeros((K, K))
    M[np.arange(K - 1), np.arange(1, K)] = r
    M[np.diag_indices(K)] = -M.sum(1)
    return M


def _one_step(t0, h):
    m = t0 + h / 2
    exact = solve_ivp(lambda t, y: (y.reshape(K, K) @ _Q(t)).ravel(), (t0, t0 + h), np.eye(K).ravel(),
                      method="DOP853", rtol=1e-13, atol=1e-15).y[:, -1].reshape(K, K)
    err = exact - expm(h * _Q(m))
    pred = h**3 / 24 * _Q(m, 2) + h**3 / 12 * (_Q(m) @ _Q(m, 1) - _Q(m, 1) @ _Q(m))
    return np.abs(err).max(), np.abs(err - pred).max()


def test_leading_local_error_term():
    e = [_one_step(10.0, h) for h in (1.0, 0.5, 0.25)]
    for (err1, res1), (err2, res2) in zip(e, e[1:]):
        assert 6.5 < err1 / err2 < 9.5      # local error O(h^3)
        assert 13 < res1 / res2 < 19        # residual after (*) is O(h^4): (*) is the leading term


def test_commutator_closed_form():
    m = 3.0
    r = np.exp(A_ + B_ * m)
    C = _Q(m) @ _Q(m, 1) - _Q(m, 1) @ _Q(m)
    closed = [2 * r[k] * r[k + 1] * abs(B_[k] - B_[k + 1]) for k in range(K - 2)] + [0.0, 0.0]
    np.testing.assert_allclose(np.abs(C).sum(axis=1), closed, rtol=1e-12, atol=1e-15)
