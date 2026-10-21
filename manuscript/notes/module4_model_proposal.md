# Module 4 — candidate model structure (PROPOSAL for the researcher)

> **Superseded 2026-10-21.** The implemented model is in
> `manuscript/derivations/module4_lc_model.md` (changes and reasons: change log of
> `prereg/dynamics.md`). Kept for the record.

> Written by Claude as a starting point. **No parameter values are given and
> none should be inferred from this note.** The researcher decides the
> equations, derives the analytic fixed point by hand (Session 5 gate) and
> writes `config/dynamics_params.yaml` with a source and confidence for every
> parameter. Nothing here is implemented in `src/dynamics/lc_model.py`.

## What the model has to be able to answer

1. Can one parameter set make **both** atomoxetine (NET inhibitor: raises
   extracellular NE) and prazosin (α1 antagonist) protective?
2. Is the system **bistable** in a "second-hit" parameter (e.g. amyloid load),
   i.e. a preserved-regulation state and a runaway state?
3. Is there a **U-shaped** dose–response (intermediate dose best)?

It must also be able to answer **no** to each. A structure that forces a
U-shape or bistability by construction would violate the brief.

## Candidate states (all non-negative, dimensionless after scaling)

| State | Meaning |
|---|---|
| E | extracellular noradrenaline in LC projection fields |
| Y | cytosolic noradrenaline in LC neurons (substrate for MAO-A → DOPEGAL) |
| R | LC regulatory integrity (α2A autoreceptor function), in [0, 1] |
| T | cortical tau burden |

Input: A(t), amyloid (a forcing or a constant second-hit parameter h).

## Candidate equations

Firing with α2A restraint (restraint weakens as R falls):

    f = f0 / (1 + kappa * R * E)

Noradrenaline release and handling (u = NET reuptake rate; atomoxetine at
dose D multiplies u by 1/(1 + D/IC_atx); v = VMAT2 packaging rate, reduced in
APOE4 carriers; m = MAO-A rate):

    dE/dt = f - u * E - dE * E
    dY/dt = phi * u * E - (v + m) * Y

DOPEGAL-driven seeding (AEP cleavage of tau) is proportional to m * Y.
α1 signalling S = E / (1 + D_praz/IC_praz) (prazosin blocks it), entering as a
harmful drive at high levels (neuroinflammation / excitotoxicity — to be
justified by the researcher or dropped).

Regulation repairs towards 1 and is damaged by local pathology (the feedback
that could create bistability):

    dR/dt = r_rep * (1 - R) - r_dam * (m * Y + h * A) * R

Cortical tau: seeding and α1-dependent spread versus clearance that depends
on intact regulation (NE oscillation amplitude in NREM sleep → glymphatic
clearance, proxied by R):

    dT/dt = (sigma0 + sigma1 * m * Y) * (1 + eta * S) * (1 + h * A) - c0 * R * T

## Why each question can come out either way

- **Atomoxetine** lowers u: less NE re-enters LC neurons (Y falls, so less
  DOPEGAL seeding — protective) but extracellular E rises (S rises — harmful
  through eta). Whether the net effect is protective, harmful or U-shaped
  depends on phi, m, eta and the IC — not fixed by structure.
- **Prazosin** lowers S only; protective iff eta > 0. If eta = 0 it does
  nothing — a possible finding.
- **Bistability** requires the R feedback (damage rising with Y, Y rising as
  R falls through firing) to be strong enough; for weak r_dam there is one
  fixed point. The bifurcation scan reports which.

## Degenerate case for validation (brief §5)

With coupling parameters kappa = phi = r_dam = sigma1 = eta = 0 and A = 0 the
system decouples: E relaxes exponentially to f0/(u + dE), Y decays as
exp(-(v + m) t), R relaxes as 1 - (1 - R0) exp(-r_rep t), and T is linear in
the remaining seeding with exponential clearance — closed forms the
researcher writes in `uncoupled_solution`.

## Parameters that would need a source

f0, kappa, u, dE, phi, v, m, IC_atx, IC_praz, eta, r_rep, r_dam, sigma0,
sigma1, c0, h. Most will be **weak** (swept over plausible ranges, never
reported as fitted values).
