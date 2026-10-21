# Module 4 — fixed points, stability and bifurcations of the LC model

> **Status: derived by Claude; MACHINE-CHECKED 2026-10-21
> (`checks/module4_fixed_points.py`, 39/39), NOT human-checked.**
> Machine-checked means that sympy re-derives the algebra and a numerical
> re-implementation tests every claim. It does not mean that a person has
> checked the derivation or judged the model. Rewritten for the model of
> `module4_lc_model.md`. The 2026-10-13 version covered an earlier scaling;
> its result "q = 1 is monostable" is kept and sharpened to q ≤ 2 (§5).
>
> Results that depend on this note:
> - the Module 4 validation gate (`tests/validation/targets/dynamics_analytic_fixed_point.yaml`, §7);
> - the registered readout ("largest root", §4);
> - the bistability window (§5).

Model and notation: `module4_lc_model.md` §3. Throughout, the doses are
held fixed, with u = u0/(1 + D_atx), w = u + 1 − u0, s = 1/(1 + D_praz) and
c = 1 − ε + εs.

## 1. T does not feed back

E, Y and R do not depend on T. The Jacobian is block lower-triangular, and T
contributes the eigenvalue −c0 c R* < 0. Prazosin enters only the T equation,
so it never changes E*, Y* or R*. [checks A8, A9]

## 2. Fixed points reduce to one equation in R

- E′ = 0: w E (1 + κRE) = 1, whose positive root is
  **E\*(R) = 2 / (w + √(w² + 4κRw))** (A1). E\* ≤ 1/w.
- Y′ = 0: Y\* = φuE\*/(v + m). The DOPEGAL load is
  x\* = mY\* = μφuE\*, with μ = m/(v + m). Since uE\* ≤ u/w < 1, **x\* < 1**.
- R′ = 0 ⇔ **G(R) ≡ R [r_rep + r_dam (H(x\*(R)) + hA)] − r_rep = 0** (A3).
- T′ = 0: **T\* = (1 + σx\*)(1 + ηsE\*)(1 + hA) / (c R\*)**. It contains no
  c0. G scales with (r_rep, r_dam) jointly, so the fixed points depend on them
  only through ρ = r_rep/r_dam (A8).

G(0+) = −r_rep < 0 and G(1) = r_dam(H + hA) ≥ 0, so **at least one fixed
point always exists in (0, 1]** (B1).

## 3. Stability: Jacobian, Routh–Hurwitz, no Hopf

The (E, Y, R) block is

    J = [ −a11    0    −a13 ]     a11 = w + κR/(1+κRE)²,  a13 = κE/(1+κRE)²
        [  a21  −a22    0   ]     a21 = φu,               a22 = v + m
        [   0   −a32  −a33  ]     a32 = r_dam R m H′(x),  a33 = r_rep + r_dam (H + hA)

with all a_ij > 0 (a32 ≥ 0). There is one feedback loop, E → Y → R → E. Its
sign is (+)(−)(−) = +, so it is a positive feedback. The characteristic
polynomial is (A4)

    (λ + a11)(λ + a22)(λ + a33) − a13 a21 a32 = λ³ + c1 λ² + c2 λ + c3.

Routh–Hurwitz for a cubic requires c1 > 0, c3 > 0 and c1c2 > c3. Here
c1 > 0 trivially, and

    c1 c2 − c3 = a11²a22 + a11²a33 + a11a22² + a22²a33 + a11a33² + a22a33² + 2a11a22a33 + a13a21a32 > 0

always (A4b). So **no Hopf bifurcation can occur**, and stability is decided
by c3 alone. Implicit differentiation gives E\*′ = −a13/a11 and
Y\*′ = (a21/a22)E\*′, hence (A4c)

    c3 = a11 a22 G′(R*).

**A fixed point is stable iff G′(R\*) > 0. If G′(R\*) < 0 it is a saddle
with exactly one positive eigenvalue.** This was checked against numerical
eigenvalues at every fixed point of 40 test sets, 4 of them bistable (B1c).

**Explicit G′.** At a fixed point, dE\*/dR = −κE²/(1 + 2κRE) (A2). So the
elasticity ε_E ≡ −R E\*′/E\* = κRE/(1 + 2κRE) lies in [0, ½) (A2b). Since
xH′(x) = qH(1 − H),

    G′(R) = r_rep + r_dam [ H + hA − ε_E q H (1 − H) ]     (A3b).

## 4. The readout: from the healthy start the system goes to the largest R\*

With R replaced by −R, every off-diagonal Jacobian entry of the full
(E, Y, −R, T) system is ≥ 0 (∂E′/∂R < 0, ∂Y′/∂E > 0, ∂R′/∂Y ≤ 0, and
∂T′/∂E, ∂T′/∂Y ≥ 0, ∂T′/∂R ≤ 0). The system is therefore **cooperative**
on the convex box E, Y, T ≥ 0, R ∈ [0, 1], and its flow preserves the order
(Kamke's theorem; Smith HL, *Monotone Dynamical Systems*, AMS 1995, ch. 3).

The healthy start x₀ = (0, 0, 1, 0) is the minimum of the box in this order,
and its velocity is ≥ 0 in the order. So the trajectory is monotone: E, Y and
T never fall, and R never rises. It is bounded, so it converges to an
equilibrium. Every equilibrium e satisfies x₀ ≤ e, so x(t) ≤ e for all t, and
the limit is the **least** equilibrium in the order. The equilibria are
totally ordered: E\* and Y\* fall as R rises, and T\* falls too. The least
equilibrium is therefore the one with the **largest R\***, and it is stable,
since G changes sign from − to + there.

The registered readout T\* is computed from that root (`lc_model.healthy_steady_state`).
B2 and B2b check this numerically (16/16 sets, including all bistable ones).
The registered run re-checks it on 50 Latin-hypercube samples.

## 5. Bistability: when it is possible, and the window in hA

**q ≤ 2 ⇒ exactly one fixed point.** With ε_E < ½ and q ≤ 2,
H − ε_E qH(1 − H) ≥ H − H(1 − H) = H² ≥ 0. So G′ > r_rep > 0 and G is
strictly increasing (A5; B4b: 0 of 300 sets). This extends the earlier
result for q = 1.

**Necessary condition for q > 2.** G′ < 0 somewhere needs
ρ + hA < H[ε_E q(1 − H) − 1] < H[q(1 − H)/2 − 1] ≤ (q − 2)²/(8q), the
maximum being at H = ½ − 1/q (A6). So

    bistability at a given hA requires  q > 2  and  r_rep/r_dam + hA < (q − 2)²/(8q)  (≤ 1/8 for q ≤ 4).

The condition also needs ε_E close to ½, which means strong α2A restraint
(κRE ≫ 1). Every bistable test set satisfied it (B4).

**The window.** Write G(R) = r_dam R (hA − Φ(R)), with
Φ(R) = ρ(1 − R)/R − H(x\*(R)) (A7) and
Φ′(R) = [−ρ/R + ε_E qH(1 − H)]/R (A7b). The fixed points at hA are the
solutions of Φ(R) = hA in (0, 1]. A root is stable iff Φ′ < 0 there.
Φ(0+) = +∞ and Φ(1) = −H(x\*(1)) ≤ 0. With one increasing stretch
[R_a, R_b], hA has three fixed points, two of them stable, iff

    max(Φ(R_a), Φ(1)) < hA < Φ(R_b).

The third root must lie at R ≤ 1, which is why Φ(1) appears. Φ(1) ≤ 0, so
on the registered range hA ∈ [0, 1] the window is (Φ(R_a), Φ(R_b)) ∩ [0, 1].
At its two ends the fixed points meet in **saddle-node (fold)
bifurcations**: G = G′ = 0, with transversality ∂G/∂hA = r_dam R > 0. Checked
against a 201-point scan of hA (B5). The implementation counts the crossings
exactly between the breakpoints of Φ (C4b).

## 6. Dose–responses

**Prazosin.** E\*, Y\* and R\* do not depend on D_praz, so
T\*(D)/T\*(0) = (1 + ηsE)/((1 + ηE)(1 − ε + εs)), with s = 1/(1 + D). Its
derivative in s has the sign of ηE(1 − ε) − ε (A9). So the prazosin
dose–response is **monotone, never U-shaped**. It is protective iff
ηE\*(1 − ε) > ε and harmful iff ηE\*(1 − ε) < ε. At D = 1 the ratio is
(1 + ηE/2)/((1 + ηE)(1 − ε/2)) (A9b).

**Atomoxetine, without feedback** (κ = r_dam = 0, so R ≡ 1): uE\* = u/w = 1 − (1 − u0)h,
with h = 1/w rising with dose. T\*(h) is then a quadratic in h with
d²T/dh² = −2σμφ(1 − u0)η(1 + hA)/c < 0 (A10). It is **concave, so the
minimum over any dose range lies at an end: no U-shape**. An inverted U
(worst at an intermediate dose) is possible.

**Atomoxetine, with feedback:** a U-shape is possible but needs the R
feedback. A found example (outside the registered box, η = 16) is U-shaped,
and the same set with κ = r_dam = 0 is not (B6e, B6g). In an early-exit
search over a wide test distribution, the first U-shape came at draw 1,681
and the first inverted U at draw 2. **A U-shape is therefore uncommon in
this model class.** This is foreknowledge, disclosed in the plan.

**Every question can come out either way** (B6a–f): atomoxetine can be
protective or harmful; prazosin can be protective or harmful (harm needs
ε > 0); bistability exists and is absent in other sets; and the atomoxetine
U-shape exists, alongside the inverted U.

## 7. The analytic fixed point used by the validation gate

Parameters: κ = 4, u0 = ½, φ = 1, v = m = 1, K = ⅛, q = 3,
r_rep = r_dam = 1/20, hA = ½, σ = 8, η = 2, no drug (w = 1, u = ½, s = c = 1).

1. Try R = ½. E solves E(1 + 2E) = 1, i.e. 2E² + E − 1 = 0, so **E = ½**.
2. **Y** = φuE/(v + m) = (1 · ½ · ½)/2 = **⅛**; x = mY = ⅛ = K, so H = ½ for any q.
3. G(½) = ½ [1/20 + (1/20)(½ + ½)] − 1/20 = 0, so **R = ½** is a root.
4. **T** = (1 + 8 · ⅛)(1 + 2 · ½)(1 + ½) / (1 · ½) = 2 · 2 · (3/2) · 2 = **12**.
5. Stability. a11 = 1 + 2/4 = 3/2; a13 = 2/4 = ½; a21 = ½; a22 = 2;
   H′(x) = qH(1 − H)/x = 6, so a32 = (1/20)(½)(1)(6) = 3/20;
   a33 = 1/20 + (1/20)(1) = 1/10.
   c3 = (3/2)(2)(1/10) − (½)(½)(3/20) = 3/10 − 3/80 = **21/80 > 0**: stable.
   The T eigenvalue is −c0/2.
6. Uniqueness: ρ + hA = 1 + ½ = 3/2 > (q − 2)²/(8q) = 1/24, so by §5 no
   other fixed point exists.

Checked exactly in sympy (A11) and against the implementation (C5). **This is
Claude's derivation, machine-checked; the brief's gate asks for the
researcher's own hand derivation, which these six steps are meant to make
quick.**
