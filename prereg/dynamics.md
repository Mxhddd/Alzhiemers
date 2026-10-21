# Can one regulatory mechanism make opposite noradrenergic drugs protective? A pre-registered simulation study of the locus coeruleus as rate controller

Mohammad Hassan, independent researcher · October 2026

This is my plan for a small mathematical model of the locus coeruleus (LC) as the rate controller of tau progression (Module 4 of my project). It is a simulation study: no patient data, measured outcomes or published effect sizes are used anywhere. I'm writing the model, the parameter ranges and the reading rules down before the registered simulations are run, so that what counts as support and what counts as failure can't move once I've seen the results. The equations are in `manuscript/derivations/module4_lc_model.md`, and the code and settings are in the project repository (`src/dynamics/`, `config/dynamics.yaml`, `config/dynamics_params.yaml`), at the commit where I register this file.

## Tools

I used Claude, Anthropic's AI model, through Claude Code to write the code, derive and machine-check the fixed points, and draft this text. I told it to decide methodological questions for me and explain each decision, which it did in the config files (every setting carries its reason) and in the change log below. By registering this plan I adopt those decisions as mine.

The fixed points and their stability were derived by Claude and checked by machine: sympy re-derives the algebra, and an independent numerical implementation tests every claim (39 of 39 checks pass, `manuscript/derivations/checks/module4_fixed_points.py`). Machine-checked is not the same as human-checked. I have not yet re-derived them by hand, and the validation gate's analytic fixed point is Claude's, not mine (section 5).

**What I already know about the model's behaviour (foreknowledge).** No registered fraction and no Sobol index has been computed. But model outputs have been looked at while the model was being built, and I list all of it here:

* The derivation shows, before any sampling, that the model needs a cooperative damage term (q > 2) to have a tipping point at all. It also needs r_rep/r_dam + hA < (q − 2)²/(8q), which is at most 1/8 in my ranges, and strong α2A restraint. So bistability can only occur in a small corner of the parameter space.
* The prazosin dose–response is monotone by construction, and the atomoxetine dose–response cannot be U-shaped without the regulation feedback.
* In an existence search over a wide test distribution (not my registered ranges), the first U-shaped atomoxetine curve appeared at draw 1,681, and the first inverted U at draw 2.
* A pilot screening (`src/dynamics/pilot.py`, its own seed 20261021, generated note `manuscript/notes/module4_pilot.md`) used Morris elementary effects over the registered ranges (750 model evaluations) and one-at-a-time sweeps around the reference point. In it, the bistable window had length zero at every one of the 750 points. The atomoxetine effect was most sensitive to the repair rate, the reuptake share and the seeding gain. The prazosin effect was most sensitive to η and ε. The dummy parameter c0 had no effect, as derived.

From this I expect bistability to be rare and the U-shape to be uncommon. I did not change the model, the ranges or the reading rules after the pilot. All of them were fixed before it, and section 2 says plainly that these outcomes count against the hypothesis.

## Change log

None of the changes below was made after a registered simulation. Earlier versions were working drafts written by Claude.

| Date | Change | Reason |
|---|---|---|
| 2026-10-13 | First complete version. It adopted the structure of `manuscript/notes/module4_model_proposal.md` with a cooperative damage term (Hill coefficient q). | The proposal's linear damage term provably has exactly one fixed point, so the bistability question would have been answered "no" by construction. |
| 2026-10-21 | **Non-reuptake clearance restored: E′ = f − wE with w = u + (1 − u0).** | In the earlier dimensionless form (E′ = f − E), atomoxetine could not raise extracellular NE, so its only effects were protective. It would have been protective and never U-shaped by construction, an analysis that could only come out one way. |
| 2026-10-21 | **New α1-dependent share ε of tau clearance (NE-driven vasomotion), range 0–0.5.** Sensitivity analysis SA2 sets ε = 0. | Without it, prazosin could be neutral but never harmful, so "prazosin protective" could not be contradicted except by η ≈ 0. |
| 2026-10-21 | **Tau is scaled by its baseline seeding (s0 removed; σ = s1/s0). c0 is kept only as a time scale**, and serves as a dummy in the sweeps. | s0 and c0 only rescale T and cancel from every ratio. Fewer redundant parameters make the sensitivity indices easier to read. |
| 2026-10-21 | **Ranges changed:** K 0.01–1 (was 0.1–10), log-uniform; η log-uniform 0.01–10 (was uniform 0–10); σ log-uniform 0.1–100. SA3 and SA4 rerun the earlier η and K ranges. | The DOPEGAL load is always below 1 in model units, so K above 1 leaves the damage switch almost off: 2 of the 3 decades of the earlier range would have made bistability nearly impossible. A uniform η on 0–10 puts 90% of samples at η > 1, which would make prazosin protective almost everywhere. |
| 2026-10-21 | **Bistability is read from the exact window of the scalar function Φ(R)** (derivation §5), not only from a 201-point scan. | A grid scan misses windows narrower than its spacing. The 201-point scan is kept as a cross-check and for the figure. |
| 2026-10-21 | **Readout computed from the largest root of G(R).** | The system is cooperative, so the healthy start provably converges to that fixed point (derivation §4). It is re-checked by ODE integration on 50 samples. |
| 2026-10-21 | **Harmful fractions, the inverted U, and verdict rules for harm added (section 2).** | They let the model contradict the Module 1 predictions, not only fail to support them. |
| 2026-10-21 | Pilot screening run, and foreknowledge disclosed (Tools). | — |

## 1. Questions

1. Can one parameter set make **both** atomoxetine (raises noradrenaline) and prazosin (blocks α1 signalling) protective?
2. Is the atomoxetine dose–response **U-shaped**?
3. Is LC regulation **bistable** in the second-hit parameter hA, that is, is there a tipping point?

The model cannot answer whether LC regulation controls tau. That is assumed: tau clearance and α2A restraint are both proportional to the regulation variable R. The model asks only whether these three consequences follow from it, and how much of the plausible parameter space they need.

## 2. What would count against the LC rate-controller account

I am stating these first, on purpose. Each fraction below is the share of 1,000 Latin-hypercube parameter sets (section 4) with the property, reported with its Wilson 95% CI, whatever it is. A fraction is labelled:

* **fragile** if the upper bound is below 10%;
* **robust** if the lower bound is at least 50%;
* **intermediate** otherwise.

| Quantity | Counts **against** | Counts **in favour** |
|---|---|---|
| (i) both drugs protective | fragile: one regulatory mechanism making both drugs protective is not a robust property of this model class | robust |
| (ii) atomoxetine U-shape | fragile: a U-shaped dose–response is not a robust prediction | robust |
| (iii) bistability for some hA in [0, 1] | fragile: a tipping point in LC regulation is not a robust property of this model class (and if no sampled set is bistable, I say so) | robust |
| atomoxetine harmful | robust: the model class predicts NET inhibition to be harmful, against my Module 1 prediction | — |
| prazosin harmful | robust: the model class predicts α1 blockade to be harmful, against my Module 1 prediction | — |

"Intermediate" means that the answer depends on parameters this model cannot pin down. I report it as that, not as support. A verdict is called robust to my choices only if its label is the same in sensitivity analyses SA1–SA4 (section 6).

Given what I already know (Tools), I expect (iii) to come out fragile, which counts against the tipping-point part of the account. These are statements about a model class, not about patients. Module 4 generates predictions to test elsewhere; it cannot confirm the hypothesis.

## 3. Data

None. No empirical data, published effect sizes or measured parameter values are used. **No published value sets any parameter.** The model is dimensionless: every parameter is a ratio of model rates or a model gain, and no published measurement maps onto one without further assumptions. All 14 parameters are therefore assumptions, labelled confidence "weak", each swept over a stated range and never reported as an estimate (`config/dynamics_params.yaml`). If a published value were ever used, it would first be extracted to `data/extracted/` with its source and verbatim text, and a test enforces this.

## 4. Model and methods

**Model** (dimensionless; time in units of the baseline noradrenaline turnover time; derivation in `manuscript/derivations/module4_lc_model.md`):

    f  = 1 / (1 + κ R E)                               firing, with α2A restraint
    E′ = f − w E,      w = u + (1 − u0),  u = u0 / (1 + D_atx)
    Y′ = φ u E − (v + m) Y                             cytosolic NE → MAO-A → DOPEGAL
    R′ = r_rep (1 − R) − r_dam (H(mY) + hA) R,   H(x) = (x/K)^q / (1 + (x/K)^q)
    T′ = c0 [ (1 + σ m Y)(1 + η S)(1 + hA) − c R T ],  S = E/(1 + D_praz),  c = 1 − ε + ε/(1 + D_praz)

Doses D are in units of each drug's IC50. The healthy start is (E, Y, R, T) = (0, 0, 1, 0).

**Parameter ranges** (log-uniform unless marked linear):

| Parameter | Meaning | Range |
|---|---|---|
| κ | α2A restraint gain | 0.1–100 |
| u0 | NET share of NE clearance | 0.1–0.9, linear |
| φ | share of reuptaken NE entering the cytosolic pool | 0.1–1, linear |
| v | VMAT2 re-packaging rate | 0.1–10 |
| m | MAO-A rate | 0.1–10 |
| r_rep | repair rate of regulation | 0.001–0.1 |
| r_dam | damage rate of regulation | 0.001–1 |
| K | half-saturation of DOPEGAL damage | 0.01–1 |
| q | cooperativity of DOPEGAL damage | 1–4, linear |
| σ | DOPEGAL/AEP-dependent seeding | 0.1–100 |
| η | α1-dependent harm | 0.01–10 |
| ε | α1-dependent share of tau clearance | 0–0.5, linear |
| hA | second hit (amyloid) | 0–1, linear |
| c0 | tau turnover rate (time scale only; a dummy) | 0.001–0.1 |

**Readout.** The steady-state tau burden T\* reached from the healthy start. The system is cooperative, so this is the fixed point with the largest R\* (derivation §4). It is computed from the roots of G(R) on a grid of 400 log-spaced points on [10⁻⁹, 10⁻²] plus 3,600 linear points on [10⁻², 1], with the extrema of G located exactly.

**Per parameter set:**

* **Drug effect:** the ratio T\*(drug at 1 IC50) / T\*(no drug). It is *protective* below 0.95 and *harmful* above 1.05.
* **U-shape:** over 0 plus 41 log-spaced doses from 0.01 to 10 IC50, the minimum of T\* lies at an interior dose and both ends exceed it by more than 5%. The mirror image (the maximum at an interior dose) is reported as an *inverted U*.
* **Bistability:** the set of hA in [0, 1] with two stable fixed points, computed exactly from Φ(R) (derivation §5). A set is bistable if that set is non-empty.

**Analyses:**

1. **Robustness:** a Latin hypercube of 1,000 parameter sets over the ranges, on each parameter's sampling scale, seed 20261019. Fractions with Wilson 95% CIs: both drugs protective, atomoxetine U-shape, bistable, each drug protective, each drug harmful, and atomoxetine inverted U.
2. **Sensitivity:** Sobol indices (Saltelli sampling, N = 1,024; first and total order) for log T\*, the log drug ratios, and the three yes/no outcomes as 0/1 indicators.
3. **Illustration only:** at the reference point (the centre of every range), the trajectory, the fixed points, the bifurcation diagram over hA, and the dose–responses. These values are not results.

## 5. Numerical validity (a result failing these is not interpreted)

* The analysis runs only if Module 4's validation passed on the current code: `make validate-dynamics`, currently 12 of 12. That includes the analytic fixed point (E, Y, R, T) = (1/2, 1/8, 1/2, 12) and its stability, derived by Claude, machine-checked, and not yet re-derived by me.
* At the reference point, Radau and BDF must agree within a relative 10⁻³ with tolerances tightened 10×. The analytic and numerical fixed points must also agree.
* For the first 50 Latin-hypercube sets, ODE integration from the healthy start (to 50 time constants of the slowest eigenvalue) must reproduce the analytic T\* within a relative 10⁻³.

## 6. Sensitivity analyses (all pre-specified)

The robustness fractions are recomputed with the same seed and size under:

* **SA1:** every range shrunk to its central half on its sampling scale.
* **SA2:** ε = 0 (no α1-dependent clearance; prazosin can then never be harmful).
* **SA3:** η uniform on [0, 10] (the earlier range).
* **SA4:** K log-uniform on [0.1, 10] (the earlier range).
* **SA5:** other effect thresholds, 0.90/1.10 and 0.99/1.01 (reported, not used to judge robustness).

## 7. Limitations, stated in advance

* An abstract, dimensionless model. Every parameter is an assumption, and every fraction depends on the ranges I chose. A fraction is a property of my plausible space, not a probability about the brain.
* The core claim (regulation controls tau) is assumed, not tested.
* R is a single regulation variable, and T is a scalar tau burden with no spatial spread. There is no sleep–wake cycle, and the second hit is constant.
* Time scales are compressed. Only steady states are interpreted, and they do not depend on the time scales.
* The cooperative damage term (q > 2 needed for a tipping point) is a modelling choice made so that bistability is possible, not forced. Its biological basis is a hypothesis.
* The derivations are machine-checked, not human-checked.
