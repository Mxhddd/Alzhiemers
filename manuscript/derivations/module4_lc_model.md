# Module 4 — the LC rate-controller model: equations, terms and assumptions

> **Status: written by Claude, 2026-10-21. The equations are implemented in
> `src/dynamics/lc_model.py`, and the implementation is MACHINE-CHECKED
> against an independent re-implementation (`checks/module4_fixed_points.py`,
> part C). NOT human-checked.** Every term below is a modelling assumption. No
> published number sets any parameter (`config/dynamics_params.yaml`).
> Fixed points and stability are derived in `module4_fixed_points.md`.

## 1. What the model is for

The model asks whether three claims of the rate-controller account can hold
in one mechanism: (i) atomoxetine and prazosin both protective, (ii) a
U-shaped atomoxetine dose–response, (iii) a tipping point (bistability) in LC
regulation as the second hit grows. It does **not** test whether LC regulation
controls tau. That is built in: tau clearance is proportional to R, and α2A
restraint is proportional to R. So the model takes the core claim as given and
can only test whether these three consequences follow from it.

## 2. State variables (all dimensionless, all ≥ 0)

| State | Meaning | Unit |
|---|---|---|
| E | extracellular noradrenaline (NE) in LC projection fields | unrestrained, drug-free steady level f_max/w0 |
| Y | cytosolic NE in LC terminals (substrate of MAO-A → DOPEGAL) | same as E |
| R | LC regulatory integrity (α2A autoreceptor function, and the NREM NE-oscillation amplitude that drives clearance), 0 ≤ R ≤ 1 | fraction of intact |
| T | cortical tau burden | burden of a fully regulated system with no second hit and no NE-dependent seeding or harm |

Time is measured in baseline extracellular-NE turnover times 1/w0. The real
gap between that time scale (minutes) and tau (years) is about 10⁶ or more.
The model's slow rates (r_rep, r_dam, c0 ≥ 10⁻³) are far less separated.
This is harmless for the registered readouts: they are steady states, and
fixed points and their stability do not depend on the time scales
(`module4_fixed_points.md` §2–3). Transients are not interpreted.

## 3. Equations

    f  = 1 / (1 + κ R E)
    E' = f − w E,                 w = u + (1 − u0),   u = u0 / (1 + D_atx)
    Y' = φ u E − (v + m) Y
    R' = r_rep (1 − R) − r_dam (H(mY) + hA) R,     H(x) = (x/K)^q / (1 + (x/K)^q)
    T' = c0 [ (1 + σ m Y)(1 + η S)(1 + hA) − c R T ]
         S = E / (1 + D_praz),    c = 1 − ε + ε / (1 + D_praz)

Healthy start: (E, Y, R, T) = (0, 0, 1, 0).

## 4. Every term

| Term | Biological meaning | Status |
|---|---|---|
| f = 1/(1 + κRE) | LC firing (release), restrained by NE acting on α2A autoreceptors. The restraint is proportional to both the agonist (E) and receptor integrity (R) | Assumption: a single hyperbolic restraint; no phasic/tonic split |
| −wE, w = u + (1 − u0) | Clearance of extracellular NE: NET reuptake u plus everything else (diffusion, uptake-2, metabolism), 1 − u0. Scaled so that w = 1 without drug | Assumption: linear clearance |
| u = u0/(1 + D_atx) | Atomoxetine blocks NET. D_atx is the dose in units of its IC50, so occupancy is D/(1 + D) | Assumption: simple competitive occupancy; no off-target effects |
| φ u E | Reuptaken NE that enters the cytosolic pool of LC terminals | Assumption |
| −(v + m) Y | Cytosolic NE is repackaged by VMAT2 (v; reduced by APOE4 in the hypothesis) or deaminated by MAO-A (m) | Assumption: first-order; no cytosolic leak from vesicles, no de novo synthesis into Y |
| x = mY | DOPEGAL production rate (MAO-A flux), the "DOPEGAL load". x = μφuE with μ = m/(v + m), so x < 1 always | Assumption: DOPEGAL is not tracked as a separate state |
| r_rep (1 − R) | Repair of LC regulation towards full integrity | Assumption |
| r_dam H(x) R | DOPEGAL-driven damage to regulation (via AEP activation, in the hypothesis), saturating, with cooperativity q | Assumption. q > 1 is needed for a tipping point to be possible at all (q ≤ 2 provably excludes it; `module4_fixed_points.md` §5). Its biological basis is a hypothesis |
| r_dam hA R | Second hit (amyloid) damaging regulation | Assumption: hA is a constant, not a time course |
| (1 + σx) | Tau seeding: baseline 1, plus DOPEGAL/AEP-cleaved tau (σ) | Assumption |
| (1 + ηS) | Harm from α1 transmission (e.g. neuroinflammation), S = α1-transmitted NE | Assumption: monotone harm, no α1 benefit except through c |
| (1 + hA) | Amyloid speeds tau production | Assumption |
| c R T | Clearance of tau. Glymphatic clearance in NREM sleep is driven by NE oscillations, proxied by regulation R. A share ε of it needs α1 transmission (NE-driven vasomotion), which prazosin blunts: c = 1 − ε + ε/(1 + D_praz) | Assumption. ε lets prazosin be harmful as well as helpful; with ε = 0 prazosin could never be harmful (sensitivity analysis SA2) |
| S = E/(1 + D_praz) | Prazosin blocks α1 receptors | Assumption: competitive occupancy |
| c0 | Tau turnover rate: only the time scale of T | Provably no effect on any steady-state readout |

The mechanisms motivating the structure come from the project brief (§2:
APOE4–VMAT2, MAO-A → DOPEGAL → AEP → tau; α2A restraint; NE oscillations
and glymphatic clearance). They set the *form* of the equations, never a
number. Published sources for these mechanisms are to be added by the
researcher. Claude did not extract any of them in this session, so none is
cited here.

## 5. Parameters

There are 14, all dimensionless, all assumptions (confidence `weak`). Each
is swept over a stated range, never reported as an estimate. The values,
ranges, scales and reasons are in `config/dynamics_params.yaml`. The
reference values are the centres of the ranges on their sampling scales. They
are used only for validation runs and illustrative figures.

**Non-identifiable combinations** (from the steady-state equations,
`module4_fixed_points.md` §2):
- v and m enter the steady states only through μ = m/(v + m);
- φ, μ and u0 enter the damage only through φμ/K;
- r_rep and r_dam enter only through ρ = r_rep/r_dam;
- c0 does not enter at all.

The sweeps vary them separately, because each is a distinct biological
lever (VMAT2, MAO-A, repair, damage).

## 6. What is not in the model (limitations)

- A single regulation variable R. α2A integrity and oscillation amplitude
  are not separated.
- No phasic/tonic firing, and no sleep–wake cycle. The three time scales of
  the brief (dimension 5) are not represented.
- A scalar tau burden, with no spatial spread (Module 6 handles spread).
- A constant second hit, not an amyloid time course.
- Atomoxetine's rise in tonic NE does not flatten oscillations (a plausible
  extra harm route that is left out).
- Drug doses are relative to IC50, with no pharmacokinetics or brain
  penetrance.
