# Module 2 — step size for the piecewise matrix exponential

> **Status: derived by Claude. NOT independently checked.** The researcher
> has not yet verified this derivation. The algebra in §2–3 was checked
> numerically (script reproduced in §6), which catches algebra slips but is
> not a substitute for a hand check. Results that depend on it: the choice
> `solver.step_years` in `config/baserate.yaml`. Check 2–3 of the
> a-posteriori audit in §5 do **not** depend on this derivation, so a mistake here can cost
> run time but cannot silently corrupt a reported fit.

## 1. Setting

Row-vector Kolmogorov forward equation for stage probabilities

    p'(t) = p(t) Q(t),    p(t0) given.

The exact one-step propagator is U(t0, t0+h) with p(t0+h) = p(t0) U,
U' = U Q(t). The method used is the **exponential midpoint rule**

    p(t0 + h) ≈ p(t0) exp(h Q(m)),    m = t0 + h/2.

It is exact when Q is constant, and more generally when Q(t1) and Q(t2)
commute for all t1, t2. The question is how large h may be when they do not.

## 2. Local error (Magnus expansion)

Write Q(m + τ) = Q0 + τ Q1 + τ²/2 Q2 + O(τ³) with Q0 = Q(m), Q1 = Q'(m),
Q2 = Q''(m), τ ∈ [−h/2, h/2].

For the column form y' = A y (A = Qᵀ) the Magnus expansion gives
U_col = exp(Ω), Ω = Ω1 + Ω2 + …, with

    Ω1 = ∫ A(τ) dτ                         = h A0 + (h³/24) A2 + O(h⁵)
    Ω2 = ½ ∫_{τ1>τ2} [A(τ1), A(τ2)] dτ2 dτ1.

For Ω2, [A0 + τ1A1, A0 + τ2A1] = (τ2 − τ1)[A0, A1] + O(τ²), and

    ∫_{−h/2}^{h/2} dτ1 ∫_{−h/2}^{τ1} (τ2 − τ1) dτ2 = −∫ (τ1 + h/2)²/2 dτ1 = −h³/6,

so Ω2 = −(h³/12)[A0, A1] + O(h⁵). Transposing back to row form
(A = Qᵀ, [A0, A1]ᵀ = [Q1, Q0]):

    U = exp( h Q0 + (h³/24) Q2 + (h³/12) [Q0, Q1] + O(h⁵) ).

The midpoint rule drops the two h³ terms, so its **local error** is

    U − exp(h Q0) = (h³/24) Q2 + (h³/12) [Q0, Q1] + O(h⁴).          (★)

(Numerical check, §6: removing (★) from the true one-step error leaves a
residual that falls 16-fold per halving of h, i.e. (★) is the exact leading
term.)

## 3. Global error

Every propagator exp(h Q) is a stochastic matrix, so in the row-vector
1-norm ‖p M‖₁ ≤ ‖p‖₁ ‖M‖∞ (‖·‖∞ = largest absolute row sum) and errors made
at one step are not amplified by later steps. Summing (★) over the N = T/h
steps on [0, T]:

    max |p_h(t) − p(t)| ≤ ‖p_h(t) − p(t)‖₁ ≤ h² ∫_0^T ε(t) dt + O(h³),
    ε(t) = ‖ Q''(t)/24 + [Q(t), Q'(t)]/12 ‖∞.

So the method is second order in h; its constant is Φ = ∫ ε dt.

### Closed form for the sequential chain

With rates r_k(t) = exp(α_k + β_k (t − t_c)/s) on the edges k → k+1, write
b_k = β_k / s (per-year log-slope). Then Q' has rates b_k r_k and Q'' has
rates b_k² r_k. Direct multiplication of the bidiagonal matrices gives

* ‖Q''‖∞ = max_k 2 b_k² r_k,
* [Q, Q'] has entries r_k r_{k+1}(b_k − b_{k+1}) at (k, k+1) and
  r_k r_{k+1}(b_{k+1} − b_k) at (k, k+2), so ‖[Q, Q']‖∞ = max_k 2 r_k r_{k+1}|b_k − b_{k+1}|.

(If all b_k are equal the commutator vanishes: Q(t) is then a scalar
function times a fixed matrix, and all Q(t) commute.) Hence

    ε(t) ≤ max_k [ r_k b_k² / 12 + r_k r_{k+1} |b_k − b_{k+1}| / 6 ].

Variants A–C add edges (B: k → c for every k < c; C: constant-rate hit
edges; A: a gated edge whose rate is constant within each amyloid band, so
its Q' is zero inside a band). Their ε has the same form plus terms of the
same order; the sequential bound is used for the a-priori choice and §5
covers every variant exactly.

## 4. A-priori choice of h

Before any fit (the parameters are unknown), bound ε over a plausible
envelope:

* rates r ≤ r_max = 0.2 per year (no stage with a mean sojourn shorter than
  5 years at any age);
* per-year log-slopes |b| ≤ b_max = 0.15 per year. For scale: the extracted
  CFAS I incidence for men rises from 10.3 (65–69) to 71.5 (85+) per 1,000
  person-years (data/extracted/dementia_incidence_cfas.csv, unverified),
  i.e. a log-slope of ln(71.5/10.3)/20 ≈ 0.097 per year; 0.15 leaves
  margin;
* adjacent slopes may differ by up to |b_k − b_{k+1}| ≤ 2 b_max = 0.30.

Then

    ε ≤ 0.2 × 0.15² / 12 + 0.2² × 0.30 / 6
      = 3.75 × 10⁻⁴ + 2.00 × 10⁻³ = 2.375 × 10⁻³ per year,
    Φ ≤ 100 years × 2.375 × 10⁻³ = 0.2375.

Target: the numerical error in any stage probability at most
δ_max = 10⁻³ (rationale: with ~2,300 brains in ~10 age bands, ~230 per
band, the binomial standard error of a stage proportion near 5% is about
0.014, so 10⁻³ is an order of magnitude below the sampling noise; §5
additionally bounds the effect on the log-likelihood, which is what matters
for small proportions). Solving h² Φ ≤ δ_max:

    h ≤ √(10⁻³ / 0.2375) = 0.065 years.

This worst-case envelope guarantees δ ≤ 10⁻³ for *any* parameters in it with
h = 0.05 years — but it is very pessimistic: it assumes every rate sits at
its maximum with maximally different slopes for a full century. Measured on
a realistic simulated 12-stage problem (the Module 2 dry run's parameter
set), the log-likelihood changes by only 0.0003 units between h = 0.05 and
h = 0.25 and by 0.0013 between 0.05 and 0.5, while h = 0.05 costs ~5–10×
the run time. So the bound is applied at realistic and then at the actual
fitted parameters instead of at the envelope's corner:

* Φ(θ) evaluated numerically at the simulated parameter set
  (`src/baserate/solve.py::error_bound`, which a unit test checks against the
  closed form above) is **Φ = 0.0099**, giving h ≤ √(10⁻³/0.0099) = 0.32 years.
* **Choice: h = 0.25 years.**
* **At every fitted estimate θ̂ the bound is recomputed**: h² Φ(θ̂) ≤ 10⁻³ is
  required (together with the Radau audit, §5); if it fails, h is halved and
  the fit repeated. The guarantee therefore holds at the parameters actually
  reported, which is where it matters.

Two further grid rules (implemented in `src/baserate/solve.py`):

* every evaluation age and every breakpoint of A(t) is a grid knot, so no
  step straddles a jump in A(t) — a jump inside a step would make that
  step's error O(h), not O(h³);
* between knots the interval is cut into equal steps no longer than h.

## 5. A-posteriori audit (pre-registered; does not rely on §2–4)

At every reported fit, the stepped stage probabilities are compared with a
Radau solution (rtol 10⁻¹⁰, integrated piecewise between breakpoints):

1. the bound h² Φ(θ̂) ≤ 10⁻³ (§4) at the fitted parameters,
2. max over ages and stages of |p_h − p_ODE| < 10⁻³, and
3. |LL(p_h) − LL(p_ODE)| < 0.1 log-likelihood units (negligible against the
   2-unit margins used for model comparison).

If any fails, h is halved and the fit repeated; the run log records both
numbers for every model.

## 6. Numerical check of §2–3 (script)

`tests/unit/test_magnus_derivation.py` runs this check on every `make test` (script
below, as first run): random 5-stage chain, rates
e^{U(−3,−1.5)}, slopes U(−0.15, 0.15).

    h=  2.0: |err|=1.732e-03  |err - predicted leading term|=1.139e-03
    h=  1.0: |err|=2.582e-04  |err - predicted leading term|=7.203e-05
    h=  0.5: |err|=3.535e-05  |err - predicted leading term|=4.519e-06
    h= 0.25: |err|=4.615e-06  |err - predicted leading term|=2.828e-07
    commutator row sums: [0.00110819 0.00115806 0.00026955 0. 0.]
    closed form         : [0.00110819 0.00115806 0.00026955 0. 0.]

The error falls ~8× per halving (local O(h³)); the residual after removing
(★) falls 16× (O(h⁴)), confirming (★) is the leading term; the commutator
closed form matches the matrix computation.

```python
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
rng = np.random.default_rng(0)
K = 5
a = rng.uniform(-3, -1.5, K - 1); b = rng.uniform(-0.15, 0.15, K - 1)
def Q(t):
    r = np.exp(a + b * t); M = np.zeros((K, K))
    M[np.arange(K-1), np.arange(1, K)] = r; M[np.diag_indices(K)] = -M.sum(1); return M
def dQ(t, n):
    r = np.exp(a + b * t) * b**n; M = np.zeros((K, K))
    M[np.arange(K-1), np.arange(1, K)] = r; M[np.diag_indices(K)] = -M.sum(1); return M
t0 = 10.0
for h in (2.0, 1.0, 0.5, 0.25):
    m = t0 + h / 2
    exact = solve_ivp(lambda t, y: (y.reshape(K, K) @ Q(t)).ravel(), (t0, t0 + h), np.eye(K).ravel(),
                      method="DOP853", rtol=1e-13, atol=1e-15).y[:, -1].reshape(K, K)
    approx = expm(h * Q(m))
    pred = h**3 / 24 * dQ(m, 2) + h**3 / 12 * (Q(m) @ dQ(m, 1) - dQ(m, 1) @ Q(m))
    err = exact - approx
    print(f"h={h:5}: |err|={np.abs(err).max():.3e}  |err - predicted leading term|={np.abs(err - pred).max():.3e}")
```
