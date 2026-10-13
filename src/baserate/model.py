"""Continuous-time Markov chain over tau stages (Module 2).

dp/dt = p(t) Q(t), with progression rates q_{k,k+1}(t) = exp(alpha_k + beta_k u),
u = (t - t_c) / s. Age is centred at t_c and scaled by s (config:
solver.age_centre, age_scale): the same model as the brief's
exp(alpha + beta t), with alpha re-expressed as the log-rate at age t_c and
beta per s years (uncentred, unscaled alpha and beta are almost collinear and
on very different scales, which stalls both the optimiser and NUTS).

Three variants share one interface (``StageModel``) so none is favoured by
construction:

* **A_gated_age** (primary A) — LC first: cortical tau (stage c) is entered
  only from the last LC stage (c-1), at rate
      mu(t) = exp(log_mu0 + beta_mu u) * (1 + (HR_A - 1) A(t)),
  A(t) = population amyloid-positivity probability at age t. The age term
  is what makes this the age-adjusted form: without it every age-related
  acceleration of the step is attributed to amyloid, because A(t) rises
  with age (manuscript/notes/module2_identifiability.md). HR_A from this
  model is DESCRIPTIVE ONLY, not a test of gating: stage-by-age counts carry
  no per-brain amyloid status (prereg/baserate.md change log).
* **A_gated_literal** (sensitivity only, "confounded with age") — the
  brief's literal form mu(t) = mu0 (1 + (HR_A - 1) A(t)) with no age term,
  i.e. mu0 + gamma A(t) with gamma = mu0 (HR_A - 1)
  (manuscript/derivations/module2_gamma.md). Any real log_hr_amyloid gives
  a positive rate.
* **B_independent** — cortical tau arises unrelated to LC stage: every
  pre-cortical stage k < c enters stage c at the same rate
  nu(t) = exp(a_nu + b_nu u), so entry does not require the LC stages.
* **C_multihit** — cortical entry needs m independent hits accumulating from
  birth at constant rate lambda each (hidden counter), independent of LC
  stage; early incidence scales as t^(m-1) (Armitage–Doll).

Each model is described ONCE, declaratively, as a list of ``Edge``s (from,
to, rate formula, parameter indices). ``build_generator`` (numpy, the
reference) and ``src/baserate/jaxlik.py`` (JAX, used for fitting) both
evaluate that same list, and a test requires them to agree to 1e-9.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np


class InvalidParameters(ValueError):
    """Parameters give a non-finite rate; the likelihood is -inf there."""


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


@dataclass(frozen=True)
class Edge:
    """One transition i -> j between hidden states.

    kind "loglinear":       rate = exp(theta[p0] + theta[p1] * u(t))
    kind "gated":           rate = exp(theta[p0]) * (1 + (exp(theta[p1]) - 1) * A(t))
    kind "gated_loglinear": rate = exp(theta[p0] + theta[p1] * u(t)) * (1 + (exp(theta[p2]) - 1) * A(t))
    kind "const":           rate = exp(theta[p0])            (p1 unused, -1)
    """

    i: int
    j: int
    kind: str
    p0: int
    p1: int
    p2: int = -1


class StepFunction:
    """Piecewise-constant function of age from (lower, upper, value) bands.
    No interpolation; ages outside every band raise. ``breakpoints`` are
    forced onto the integration grid so no step straddles a jump."""

    def __init__(self, bands: list[tuple[float, float, float]]):
        bands = sorted((float(a), float(b), float(v)) for a, b, v in bands)
        for (l1, u1, _), (l2, _, _) in zip(bands, bands[1:]):
            if l2 < u1:
                raise ValueError(f"overlapping bands [{l1},{u1}) and [{l2},...)")
        for lo, hi, v in bands:
            if not lo < hi:
                raise ValueError(f"empty band [{lo},{hi})")
        self.bands = bands

    def __call__(self, t: float) -> float:
        for lo, hi, v in self.bands:
            if lo <= t < hi:
                return v
        raise ValueError(f"age {t} is outside every supplied band; refusing to extrapolate")

    @property
    def breakpoints(self) -> list[float]:
        return sorted({x for lo, hi, _ in self.bands for x in (lo, hi)})


def piecewise_constant(bands) -> StepFunction:
    return StepFunction(bands)


class StageModel(ABC):
    name: str

    def __init__(self, spec: ChainSpec, age_centre: float = 0.0, age_scale: float = 1.0):
        self.spec = spec
        self.age_centre = float(age_centre)
        if not age_scale > 0:
            raise ValueError("age_scale must be positive")
        self.age_scale = float(age_scale)
        self._names: list[str] = []
        self._edges: list[Edge] = []
        self._build()

    # -- declarative description -------------------------------------------
    def _p(self, name: str) -> int:
        self._names.append(name)
        return len(self._names) - 1

    def _seq_edge(self, i: int, j: int, k: int) -> None:
        a, b = self._p(f"alpha[{k}]"), self._p(f"beta[{k}]")
        self._edges.append(Edge(i, j, "loglinear", a, b))

    @abstractmethod
    def _build(self) -> None: ...

    @property
    def param_names(self) -> list[str]:
        return list(self._names)

    @property
    def edges(self) -> list[Edge]:
        return list(self._edges)

    @property
    def n_hidden(self) -> int:
        return self.spec.K

    @property
    def observed_of_hidden(self) -> np.ndarray:
        return np.arange(self.spec.K)

    @property
    def amyloid(self) -> StepFunction | None:
        return None

    @property
    def breakpoints(self) -> list[float]:
        return self.amyloid.breakpoints if self.amyloid is not None else []

    def u(self, t: float) -> float:
        """Rescaled age entering every log-linear rate: (t - centre) / scale."""
        return (t - self.age_centre) / self.age_scale

    def initial(self) -> np.ndarray:
        p = np.zeros(self.n_hidden)
        p[0] = 1.0
        return p

    def unpack(self, theta: np.ndarray) -> dict[str, float]:
        theta = np.asarray(theta, dtype=float)
        if theta.shape != (len(self._names),):
            raise ValueError(f"{self.name} expects {len(self._names)} parameters, got {theta.shape}")
        return dict(zip(self._names, theta))

    # -- numpy evaluation (reference implementation) -------------------------
    def transitions(self, t: float, theta: np.ndarray) -> list[tuple[int, int, float]]:
        theta = np.asarray(theta, dtype=float)
        if theta.shape != (len(self._names),):
            raise ValueError(f"{self.name} expects {len(self._names)} parameters, got {theta.shape}")
        u = self.u(t)
        a = self.amyloid(t) if self.amyloid is not None else None
        out = []
        with np.errstate(over="ignore", invalid="ignore"):
            for e in self._edges:
                if e.kind == "loglinear":
                    r = np.exp(theta[e.p0] + theta[e.p1] * u)
                elif e.kind in ("gated", "gated_loglinear"):
                    if not 0.0 <= a <= 1.0:
                        raise ValueError(f"amyloid probability {a} at age {t} outside [0, 1]")
                    if e.kind == "gated":
                        r = np.exp(theta[e.p0]) * (1.0 + (np.exp(theta[e.p1]) - 1.0) * a)
                    else:
                        r = np.exp(theta[e.p0] + theta[e.p1] * u) * (1.0 + (np.exp(theta[e.p2]) - 1.0) * a)
                elif e.kind == "const":
                    r = np.exp(theta[e.p0])
                else:
                    raise ValueError(f"unknown edge kind {e.kind}")
                if not np.isfinite(r):
                    raise InvalidParameters(f"{self.name}: rate {e.i}->{e.j} at age {t} is not finite")
                out.append((e.i, e.j, float(r)))
        return out


def build_generator(model: StageModel, t: float, theta: np.ndarray) -> np.ndarray:
    """Q(t). Off-diagonals are rates; rows sum to zero (asserted)."""
    Q = np.zeros((model.n_hidden, model.n_hidden))
    for i, j, r in model.transitions(t, theta):
        if r < 0:
            raise InvalidParameters(f"{model.name}: negative rate {i}->{j} at age {t}")
        Q[i, j] += r
    Q[np.diag_indices_from(Q)] = -Q.sum(axis=1)
    if not np.allclose(Q.sum(axis=1), 0.0, atol=1e-12 * max(1.0, np.abs(Q).max())):
        raise AssertionError("generator rows do not sum to zero")
    return Q


class Sequential(StageModel):
    """Plain sequential chain; the shared backbone and the analytic test case."""

    name = "sequential"

    def _build(self):
        for k in range(self.spec.K - 1):
            self._seq_edge(k, k + 1, k)


class AGated(StageModel):
    """The brief's literal variant A (no age term on the gated step).
    Sensitivity analysis only: confounded with age."""

    name = "A_gated_literal"

    def __init__(self, spec: ChainSpec, amyloid: StepFunction, age_centre: float = 0.0, age_scale: float = 1.0):
        if not isinstance(amyloid, StepFunction):
            raise TypeError("amyloid must be a StepFunction (piecewise constant, with breakpoints)")
        self._amyloid = amyloid
        super().__init__(spec, age_centre, age_scale)

    @property
    def amyloid(self) -> StepFunction:
        return self._amyloid

    def _build(self):
        c = self.spec.first_cortical
        for k in range(self.spec.K - 1):
            if k == c - 1:
                m, h = self._p("log_mu0"), self._p("log_hr_amyloid")
                self._edges.append(Edge(k, k + 1, "gated", m, h))
            else:
                self._seq_edge(k, k + 1, k)


class AGatedAge(AGated):
    """Primary variant A: the gated step also has its own age slope."""

    name = "A_gated_age"

    def _build(self):
        c = self.spec.first_cortical
        for k in range(self.spec.K - 1):
            if k == c - 1:
                m, b, h = self._p("log_mu0"), self._p("beta_mu"), self._p("log_hr_amyloid")
                self._edges.append(Edge(k, k + 1, "gated_loglinear", m, b, h))
            else:
                self._seq_edge(k, k + 1, k)


class BIndependent(StageModel):
    name = "B_independent"

    def _build(self):
        c = self.spec.first_cortical
        for k in range(c - 1):
            self._seq_edge(k, k + 1, k)
        a, b = self._p("a_nu"), self._p("b_nu")
        for k in range(c):
            self._edges.append(Edge(k, c, "loglinear", a, b))
        for k in range(c, self.spec.K - 1):
            self._seq_edge(k, k + 1, k)


class CMultiHit(StageModel):
    """Hidden states: (LC stage k < c, hits h < m) laid out k*m + h, then the
    cortical stages c..K-1."""

    name = "C_multihit"

    def __init__(self, spec: ChainSpec, hits: int, age_centre: float = 0.0, age_scale: float = 1.0):
        if not isinstance(hits, int) or hits < 1:
            raise ValueError("hits must be a positive integer")
        self.m = hits
        super().__init__(spec, age_centre, age_scale)

    def _h(self, k: int, h: int) -> int:
        return k * self.m + h

    def _cort(self, k: int) -> int:
        return self.spec.first_cortical * self.m + (k - self.spec.first_cortical)

    def _build(self):
        c = self.spec.first_cortical
        seq = {}
        for k in range(c - 1):
            seq[k] = (self._p(f"alpha[{k}]"), self._p(f"beta[{k}]"))
        lam = self._p("log_lambda")
        for k in range(c):
            for h in range(self.m):
                if k < c - 1:
                    self._edges.append(Edge(self._h(k, h), self._h(k + 1, h), "loglinear", *seq[k]))
                target = self._h(k, h + 1) if h < self.m - 1 else self._cort(c)
                self._edges.append(Edge(self._h(k, h), target, "const", lam, -1))
        for k in range(c, self.spec.K - 1):
            self._seq_edge(self._cort(k), self._cort(k + 1), k)

    @property
    def n_hidden(self) -> int:
        c = self.spec.first_cortical
        return c * self.m + (self.spec.K - c)

    @property
    def observed_of_hidden(self) -> np.ndarray:
        c = self.spec.first_cortical
        obs = [k for k in range(c) for _ in range(self.m)]
        return np.array(obs + list(range(c, self.spec.K)))
