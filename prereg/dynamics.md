# Pre-registration — Module 4: ODE dynamical model of LC regulation

**Status: COMPLETE DRAFT, ready for the researcher to review, upload to OSF
and freeze — not yet frozen.** Written by Claude under the standing policy.
Claude must not run `make freeze-prereg`. This module simulates; it reads
no patient data. Its conclusions are statements about a model, and every
parameter is an assumption (§3).

Implementation status: the equations below are not yet coded in
`src/dynamics/lc_model.py`, and `config/dynamics_params.yaml` is not yet
filled. That is the next Module 4 step, and it must match this document.

## Change log

| Date | Change |
|---|---|
| 2026-10-13 | Complete draft. Replaces `dynamics.DRAFT.md`. Adopts the model of `manuscript/notes/module4_model_proposal.md` with one change. The damage term on LC regulation becomes cooperative (Hill coefficient q, swept over 1–4). The machine-checked derivation `manuscript/derivations/module4_fixed_points.md` §3 proves that the proposal's linear damage term (q = 1) has exactly one fixed point for every parameter set. With it, the bistability question would be answered "no" by construction. |

## 1. Questions

1. Can one parameter set make **both** atomoxetine and prazosin protective?
2. Is the system **bistable** in the second-hit parameter h·A?
3. Is the atomoxetine dose–response **U-shaped**?

Each must be able to come out "no". §5 states what counts against.

## 2. Model (dimensionless)

    f  = 1 / (1 + κ R E)                               (firing, α2A restraint)
    E' = f − E,     with reuptake share u = u0 / (1 + D_atx)      (time unit: NE turnover)
    Y' = φ u E − (v + m) Y                             (cytosolic NE → MAO-A → DOPEGAL)
    R' = r_rep (1 − R) − r_dam (H(mY) + hA) R,   H(x) = (x/K)^q / (1 + (x/K)^q)
    T' = (s0 + s1 m Y)(1 + η S)(1 + hA) − c0 R T,  S = E / (1 + D_praz)

Doses D are in units of the drugs' IC values. The fixed points reduce to
one equation in R (derivation §2):

- a root of G(R) is stable iff G′(R) > 0;
- q = 1 is provably monostable;
- q > 1 allows bistability but does not force it (1 of 200 random
  parameter sets in the check).

## 3. Parameters

**All parameters are weak** (no published value for this abstract model).
Each is swept over a log-uniform range, linear for fractions. Point values
are never reported as estimates.

| Parameter | Meaning | Range |
|---|---|---|
| κ | α2A restraint strength | 0.1–100 |
| u0 | reuptake share of NE clearance | 0.1–0.9 (linear) |
| φ | share of reuptaken NE reaching the cytosol | 0.1–1 (linear) |
| v | VMAT2 repackaging rate | 0.1–10 |
| m | MAO-A rate | 0.1–10 |
| r_rep | repair rate of regulation | 0.001–0.1 |
| r_dam | damage rate | 0.001–1 |
| K | damage half-saturation (units of mY) | 0.1–10 |
| q | damage cooperativity | 1–4 (linear) |
| s0 | baseline seeding | 0.01–1 |
| s1 | DOPEGAL-dependent seeding | 0–10 (linear) |
| η | α1-dependent harm | 0–10 (linear) |
| c0 | regulation-dependent clearance | 0.001–0.1 |
| hA | second hit, the bifurcation parameter | 0–1 |

The ranges span two orders of magnitude around the scale set by NE
turnover. They are assumptions, chosen before any simulation.

## 4. Pre-specified analyses and readouts

- **Readout:** steady-state cortical tau burden T\* at the stable fixed
  point. Where there are two stable states, T\* is taken in the state
  reached from the healthy start (R0 = 1, E0 = Y0 = T0 = 0).
- **Protective:** T\*(drug at D = 1) / T\*(no drug) < 0.95.
- **U-shape:** over D_atx ∈ [0, 10] (41 log-spaced doses plus 0), both
  ends exceed the interior minimum of T\* by more than 5%.
- **Bistability:** at least two stable fixed points for some hA in the
  201-point scan over [0, 1].
- **Robustness:** a Latin-hypercube sample of 1,000 parameter sets over
  §3 gives the fraction with:
  - (i) both drugs protective;
  - (ii) an atomoxetine U-shape;
  - (iii) bistability.

  Each fraction is reported with a Wilson 95% CI **whatever it is**.
  Labels:
  - **fragile:** upper bound < 10%;
  - **robust:** lower bound ≥ 50%;
  - **intermediate:** otherwise.
- **Sensitivity:** Sobol indices (Saltelli, N = 1,024) for T\* and for the
  three outcomes.
- **Numerical validity:**
  - Results must agree between Radau and BDF with tolerances tightened
    10× (relative disagreement ≤ 10⁻³).
  - The analytic fixed point (derivation §2) must match the numerical one.
  - The degenerate case (κ = φ = r_dam = s1 = η = 0, A = 0) must reduce
    to its closed forms.

  A result failing these is not interpreted.

## 5. What counts against the LC rate-controller account (stated first)

- **(i) fragile:** the claim that one regulatory mechanism makes both
  atomoxetine and prazosin protective is not supported by this model
  class.
- **(iii) zero:** in no sampled parameter set does LC regulation show a
  tipping point. The "runaway" framing is not supported by this model
  class.
- **(ii) fragile:** a U-shaped dose–response is not a robust prediction.

These are statements about the model class, not about patients. Module 4
generates predictions to test elsewhere; it cannot confirm the hypothesis.

## 6. Limitations, stated in advance

- An abstract, dimensionless model.
- All parameters are assumptions.
- R is a single "regulation" variable.
- T is a scalar tau burden, with no spatial spread.
- The q > 1 damage form is a modelling choice made to keep bistability
  possible, not forced. Its biological basis (cooperative DOPEGAL / AEP
  damage) is a hypothesis.
