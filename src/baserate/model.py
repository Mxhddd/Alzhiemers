"""Continuous-time Markov chain over tau stages (Module 2).

dp/dt = p(t) Q(t), with progression rates q_{k,k+1}(t) = exp(alpha_k + beta_k (t - t_c) / s).
Age is centred at t_c and scaled by s (config: solver.age_centre, age_scale):
the same model as the brief's exp(alpha + beta t), with alpha re-expressed as
the log-rate at age t_c and beta per s years. Uncentred, unscaled alpha and
beta are almost collinear and on very different scales, which stalls both the
optimiser and NUTS.
Three variants share one interface (``StageModel``) so none is favoured by
construction:

* **A_gated** — the LC -> cortex step (stage c-1 -> c) has rate
  mu(t) = mu0 + gamma * A(t), A(t) = amyloid-positivity probability at age t.
  gamma is unconstrained in sign: the fit can say amyloid does not help.
* **B_independent** — cortical tau arises unrelated to LC stage: every
  pre-cortical stage k < c enters stage c at the same rate nu(t) =
  exp(a_nu + b_nu t), so entry does not require passing through the LC stages.
* **C_multihit** — cortical entry needs m independent hits accumulating from
  birth at constant rate lambda each (hidden counter), independent of LC
  stage; early incidence scales as t^(m-1) (Armitage–Doll).

Stage c is ``first_cortical``; stages are indexed 0..K-1 in the
researcher-supplied order. Hidden states (variant C) map onto observed stages
through ``observed_of_hidden``.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np


class InvalidParameters(ValueError):
    """Parameters give a negative rate; the likelihood is -inf there."""


@dataclass(frozen=True)
class ChainSpec:
    stages: tuple[str, ...]
    first_cortical: int  # index c of the first cortical stage

    def __post_init__(self):
        if len(self.stages) < 2:
            raise ValueError("need at least two stages")
        if len(set(self.stages)) != len(self.stages):
            raise ValueError("stage labels must be unique")
        if not 1 <= self.first_cortical <= len(self.stages) - 1:
            raise ValueError("first_cortical must index a stage after the first")

    @classmethod
    def from_labels(cls, stages: list[str], first_cortical: str) -> "ChainSpec":
        stages = [str(s) for s in stages]
        if str(first_cortical) not in stages:
            raise ValueError(f"first cortical stage {first_cortical!r} not in stage_order")
        return cls(tuple(stages), stages.index(str(first_cortical)))

    @property
    def K(self) -> int:
        return len(self.stages)


def _rate(a: float, b: float, t: float) -> float:
    with np.errstate(over="raise"):
        try:
            return float(np.exp(a + b * t))
        except FloatingPointError as e:
            raise InvalidParameters(f"rate overflow: exp({a} + {b} * {t})") from e


class StageModel(ABC):
    name: str

    def __init__(self, spec: ChainSpec, age_centre: float = 0.0, age_scale: float = 1.0):
        self.spec = spec
        self.age_centre = float(age_centre)
        if not age_scale > 0:
            raise ValueError("age_scale must be positive")
        self.age_scale = float(age_scale)

    def _u(self, t: float) -> float:
        """Rescaled age entering every rate: (t - centre) / scale."""
        return (t - self.age_centre) / self.age_scale

    @property
    @abstractmethod
    def param_names(self) -> list[str]: ...

    @property
    def n_hidden(self) -> int:
        return self.spec.K

    @property
    def observed_of_hidden(self) -> np.ndarray:
        return np.arange(self.spec.K)

    def initial(self) -> np.ndarray:
        p = np.zeros(self.n_hidden)
        p[0] = 1.0
        return p

    @abstractmethod
    def transitions(self, t: float, theta: np.ndarray) -> list[tuple[int, int, float]]:
        """(from, to, rate) triples among hidden states at age t."""

    def unpack(self, theta: np.ndarray) -> dict[str, float]:
        theta = np.asarray(theta, dtype=float)
        if theta.shape != (len(self.param_names),):
            raise ValueError(f"{self.name} expects {len(self.param_names)} parameters, got {theta.shape}")
        return dict(zip(self.param_names, theta))

    def _seq(self, k: int) -> tuple[str, str]:
        return f"alpha[{k}]", f"beta[{k}]"


def build_generator(model: StageModel, t: float, theta: np.ndarray) -> np.ndarray:
    """Q(t). Off-diagonals are rates; rows sum to zero (asserted)."""
    Q = np.zeros((model.n_hidden, model.n_hidden))
    for i, j, r in model.transitions(t, theta):
        if not np.isfinite(r) or r < 0:
            raise InvalidParameters(f"{model.name}: rate {i}->{j} at age {t} is {r}")
        Q[i, j] += r
    Q[np.diag_indices_from(Q)] = -Q.sum(axis=1)
    if not np.allclose(Q.sum(axis=1), 0.0, atol=1e-12 * max(1.0, np.abs(Q).max())):
        raise AssertionError("generator rows do not sum to zero")
    return Q


class Sequential(StageModel):
    """Plain sequential chain; the shared backbone and the analytic test case."""

    name = "sequential"

    @property
    def param_names(self) -> list[str]:
        return [n for k in range(self.spec.K - 1) for n in self._seq(k)]

    def transitions(self, t, theta):
        p = self.unpack(theta)
        return [(k, k + 1, _rate(p[f"alpha[{k}]"], p[f"beta[{k}]"], self._u(t))) for k in range(self.spec.K - 1)]


class AGated(StageModel):
    name = "A_gated"

    def __init__(self, spec: ChainSpec, amyloid: Callable[[float], float], age_centre: float = 0.0,
                 age_scale: float = 1.0):
        super().__init__(spec, age_centre, age_scale)
        self.amyloid = amyloid

    @property
    def param_names(self) -> list[str]:
        c = self.spec.first_cortical
        names = []
        for k in range(self.spec.K - 1):
            names += ["log_mu0", "gamma"] if k == c - 1 else list(self._seq(k))
        return names

    def transitions(self, t, theta):
        p = self.unpack(theta)
        c = self.spec.first_cortical
        out = []
        for k in range(self.spec.K - 1):
            if k == c - 1:
                a = self.amyloid(t)
                if not 0.0 <= a <= 1.0:
                    raise ValueError(f"amyloid probability {a} at age {t} outside [0, 1]")
                out.append((k, k + 1, np.exp(p["log_mu0"]) + p["gamma"] * a))
            else:
                out.append((k, k + 1, _rate(p[f"alpha[{k}]"], p[f"beta[{k}]"], self._u(t))))
        return out


class BIndependent(StageModel):
    name = "B_independent"

    @property
    def param_names(self) -> list[str]:
        c = self.spec.first_cortical
        names = [n for k in range(c - 1) for n in self._seq(k)]
        names += ["a_nu", "b_nu"]
        names += [n for k in range(c, self.spec.K - 1) for n in self._seq(k)]
        return names

    def transitions(self, t, theta):
        p = self.unpack(theta)
        c = self.spec.first_cortical
        out = [(k, k + 1, _rate(p[f"alpha[{k}]"], p[f"beta[{k}]"], self._u(t))) for k in range(c - 1)]
        nu = _rate(p["a_nu"], p["b_nu"], self._u(t))
        out += [(k, c, nu) for k in range(c)]
        out += [(k, k + 1, _rate(p[f"alpha[{k}]"], p[f"beta[{k}]"], self._u(t))) for k in range(c, self.spec.K - 1)]
        return out


class CMultiHit(StageModel):
    """Hidden states: (LC stage k < c, hits h < m) laid out k*m + h, then the
    cortical stages c..K-1."""

    name = "C_multihit"

    def __init__(self, spec: ChainSpec, hits: int, age_centre: float = 0.0, age_scale: float = 1.0):
        super().__init__(spec, age_centre, age_scale)
        if not isinstance(hits, int) or hits < 1:
            raise ValueError("hits must be a positive integer")
        self.m = hits

    @property
    def param_names(self) -> list[str]:
        c = self.spec.first_cortical
        names = [n for k in range(c - 1) for n in self._seq(k)]
        names += ["log_lambda"]
        names += [n for k in range(c, self.spec.K - 1) for n in self._seq(k)]
        return names

    @property
    def n_hidden(self) -> int:
        c = self.spec.first_cortical
        return c * self.m + (self.spec.K - c)

    def _h(self, k: int, h: int) -> int:
        return k * self.m + h

    def _cort(self, k: int) -> int:
        return self.spec.first_cortical * self.m + (k - self.spec.first_cortical)

    @property
    def observed_of_hidden(self) -> np.ndarray:
        c = self.spec.first_cortical
        obs = [k for k in range(c) for _ in range(self.m)]
        return np.array(obs + list(range(c, self.spec.K)))

    def transitions(self, t, theta):
        p = self.unpack(theta)
        c = self.spec.first_cortical
        lam = float(np.exp(p["log_lambda"]))
        out = []
        for k in range(c):
            for h in range(self.m):
                if k < c - 1:
                    out.append((self._h(k, h), self._h(k + 1, h), _rate(p[f"alpha[{k}]"], p[f"beta[{k}]"], self._u(t))))
                target = self._h(k, h + 1) if h < self.m - 1 else self._cort(c)
                out.append((self._h(k, h), target, lam))
        for k in range(c, self.spec.K - 1):
            out.append((self._cort(k), self._cort(k + 1), _rate(p[f"alpha[{k}]"], p[f"beta[{k}]"], self._u(t))))
        return out


def piecewise_constant(bands: list[tuple[float, float, float]]) -> Callable[[float], float]:
    """Step function from (lower, upper, value) bands — no interpolation.
    Ages outside every band raise rather than being extrapolated."""
    bands = sorted(bands)
    for (l1, u1, _), (l2, _, _) in zip(bands, bands[1:]):
        if l2 < u1:
            raise ValueError(f"overlapping bands [{l1},{u1}) and [{l2},...)")

    def f(t: float) -> float:
        for lo, hi, v in bands:
            if lo <= t < hi:
                return float(v)
        raise ValueError(f"age {t} is outside every supplied band; refusing to extrapolate")

    return f
