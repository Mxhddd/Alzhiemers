# Module 4 — fixed points and stability of the LC model

> **Status: derived by Claude; MACHINE-CHECKED 2026-10-13
> (`checks/module4_fixed_points.py`, 13/13), NOT human-checked.** It was
> written for the model adopted in `prereg/dynamics.md`. Results that depend
> on it:
> - the analytic fixed point used by the Module 4 validation gate;
> - the finding (§3) that the proposal's original structure cannot be
>   bistable, which is why the pre-registered model generalises the
>   damage term.

## 1. Model (the T equation does not feed back)

    f  = f0 / (1 + κ R E)
    E' = f − w E,                     w = u + d,   u = u0 / (1 + D_atx/IC_atx)
    Y' = φ u E − (v + m) Y
    R' = r_rep (1 − R) − r_dam (H(mY) + hA) R,   H(x) = (x/K)^q / (1 + (x/K)^q)
    T' = (s0 + s1 m Y)(1 + η S)(1 + hA) − c0 R T,    S = E / (1 + D_praz/IC_praz)

q = 1 with large K is the proposal's linear damage term, rescaled
(r_dam · mY / K).

## 2. Fixed points reduce to one equation in R

- E′ = 0: w E (1 + κRE) = f0, a quadratic in E with one positive root:

      E*(R) = [−w + √(w² + 4κR w f0)] / (2κR w),    E*(R) → f0/w as κR → 0.

- Y′ = 0: Y*(R) = φ u E*(R) / (v + m).
- R′ = 0 with E*, Y* substituted:

      G(R) ≡ R [r_rep + r_dam (H(m Y*(R)) + hA)] − r_rep = 0.

- T′ = 0: T* = (s0 + s1 m Y*)(1 + η S*)(1 + hA) / (c0 R*).

Every fixed point corresponds to one root of G in (0, 1]. G(0) = −r_rep < 0
and G(1) = r_dam (H + hA) ≥ 0, so at least one fixed point always exists.

**Stability.** T′ has eigenvalue −c0 R* < 0. For the (E, Y, R) block, a root
is stable iff G′(R*) > 0. This was checked numerically against the Jacobian
eigenvalues at every root of 60 random parameter sets: all consistent. The
ODE, integrated from 8 initial conditions per set, always settles on a
stable root of G.

## 3. The proposal's structure (q = 1) cannot be bistable

By implicit differentiation, dE*/dR = −κE²/(1 + 2κRE) < 0: more restraint
means less extracellular NE. Hence

    R dE*/dR + E*/2 = E* / (2(1 + 2κRE*)) > 0,   i.e.  R dE*/dR > −E*/2,

and with damage linear in Y:

    G′(R) = r_rep + r_dam (m Y* + hA) + r_dam R m Y*′(R)
          > r_rep + r_dam (m Y*/2 + hA) > 0.

G is strictly increasing, so there is **exactly one fixed point for every
parameter set**. The firing restraint makes NE (and hence the DOPEGAL
load) rise only like R^(−1/2) as regulation falls. That is too weak to
make the damage feedback self-reinforcing. Checked symbolically, and
numerically on 60 random sets (0 bistable).

A model in which bistability is impossible by construction would answer
the brief's bistability question "no" before any analysis. The proposal
itself rules that out.

## 4. Generalised damage (adopted)

A cooperative damage term H with Hill coefficient q > 1 lets the damage
rate rise steeply once the DOPEGAL load passes K. G can then be
non-monotone and have three roots (two stable, one unstable), but it does
not have to.

Numerically, with q ∈ {2, 4} and K log-uniform in [0.1, 10], 1 of 200
random parameter sets was bistable and 199 had a single fixed point. In
the bistable set, the ODE started from R0 = 0.01 settles in the low-R
state, and from R0 = 0.99 in the high-R state. So the question can come
out either way, and bistability is rare in this random space rather than
built in.

## 5. Degenerate case (validation)

κ = φ = r_dam = s1 = η = 0 and A = 0 gives:

- E(t) → f0/w exponentially at rate w;
- Y decays at rate v + m;
- R = 1 − (1 − R0) e^(−r_rep t);
- T is linear seeding s0 with clearance c0 R(t).

With κ = 0 alone, G is linear in R and has one root (checked).
