# Pre-registration — Module 2: Markov base-rate model of tau stage progression

**Status: COMPLETE DRAFT for researcher review — not yet frozen.** Written by
Claude under the standing policy (Claude decides methodological questions and
states them; the researcher reviews, uploads to OSF and runs
`make freeze-prereg NAME=baserate OSF=<url>`). Nothing in Module 2 has been
fitted to real data. Configuration: `config/baserate.yaml` at the commit
that freezes this file.

> **DO NOT FREEZE — open issue (raised 2026-10-13, awaiting the researcher's
> decision).** Criterion 1 (§2) cannot reliably come out against the
> hypothesis. Variant A's LC → cortex rate μ0(1 + (HR_A − 1)A(t)) has no age
> term of its own, and A(t) rises with age, so ANY age-related acceleration
> of that transition is read as amyloid gating. Simulated check: with no
> amyloid effect at all and an age slope of 0.1 per decade on that edge,
> the fit gives HR_A = 2.22 (95% CI on the log scale 0.28–1.31), i.e.
> "supports gating" and "consistent" with the published 2.43. Details,
> further checks and the proposed fixes:
> `manuscript/notes/module2_identifiability.md`. `make freeze-prereg`
> refuses this file while this marker is present.

## 1. Question

Nearly everyone has locus coeruleus (LC) tau by mid-life, yet only a
minority develop Alzheimer's dementia. Which model of the transition from
LC-only tau to cortical tau best explains the age distribution of tau stages
at autopsy, and does it predict the published age-specific dementia
incidence it never sees?

## 2. What would count against the LC rate-controller hypothesis (stated first)

Each criterion has three pre-defined outcomes, because the simulated dry run
showed that with ~2,300 brains a true amyloid HR of 2 can still give a 95%
interval that includes 1. An interval that includes 1 is therefore NOT by
itself evidence against the hypothesis.

1. **Amyloid gating of the LC → cortex transition (variant A, log HR_A,
   95% credible interval):**
   * *supports gating* — the interval lies entirely above 0 (HR_A > 1);
   * *against gating* — the interval's upper end is below 0.322 = ln 1.38,
     i.e. the data rule out even the smallest effect compatible with the
     published hazard ratio (§7);
   * *inconclusive* — anything else; reported as uninformative, with the
     interval, never as support or refutation.
2. **Which variant explains the stage distribution (5-fold CV by age band):**
   * *supports* — A beats both B and C by > 2 held-out log-likelihood units;
   * *against* — B or C beats A by > 2 units;
   * *inconclusive* — otherwise.
3. **Falsification** (§6): every variant fails the incidence test → the
   framework itself is inadequate, whatever 1–2 say.

The hypothesis as stated is supported only if 1 and 2 both support and A
passes §6; it counts against it if 1 or 2 is "against".

## 3. Data (all extracted by Claude, `verified: false` until checked)

| Input | Source | File |
|---|---|---|
| Tau stage counts by age band (stages 0, a, b, c, 1a, 1b, I–VI) | Braak H, Thal DR, Ghebremedhin E, Del Tredici K. *J Neuropathol Exp Neurol* 2011;70:960–969 | `data/extracted/braak2011_stage_counts.csv` (**pending the full text**) |
| Amyloid-plaque prevalence by age band, same cohort (primary A(t)) | Braak et al. 2011 | `data/extracted/braak2011_amyloid_by_age.csv` (**pending**) |
| Amyloid-positivity prevalence, cognitively normal (sensitivity A(t)) | Jansen WJ et al. *JAMA* 2015;313:1924–1938, Table 2 | `data/extracted/amyloid_prevalence_jansen2015.csv` |
| Dementia incidence by age (falsification only) | Matthews FE et al. *Nat Commun* 2016;7:11398, Table 1 | `data/extracted/dementia_incidence_cfas.csv` |
| Amyloid hazard ratio (expected γ) | Ossenkoppele R et al. *Nat Med* 2022;28:2381–2387, Table 2 | `data/extracted/amyloid_hr_ossenkoppele2022.csv` |

Consistency check before fitting: the extracted Braak table must reproduce
the abstract totals (2,332 brains; 10 without AT8 tau; 58 at stages a–c;
38 at 1a; 236 at 1b; 1,031 with plaques). A mismatch stops the analysis.

## 4. Model

Continuous-time Markov chain, dp/dt = p Q(t), stages in Braak 2011 order
0 → a → b → c → 1a → 1b → I → … → VI; everyone at stage 0 at age 1 (youngest
age in the series). Progression rates exp(α_k + β_k (t − 50)/10).

* **A, gated:** c → 1a at rate μ0 + γA(t) = μ0(1 + (HR_A − 1)A(t)).
* **B, independent:** every pre-cortical stage enters 1a at a common rate
  exp(a + b(t − 50)/10).
* **C, multi-hit:** 1a is entered after m = 2 independent hits at constant
  rate (sensitivity m = 3).

Likelihood: multinomial counts per age band, stage probabilities at the band
midpoint. Numerical solution: exponential midpoint rule with step 0.25 years
(manuscript/derivations/module2_step_size.md). At every fitted estimate the
derived error bound h²Φ(θ̂) must be ≤ 10⁻³ and a Radau audit must give
max |Δp| < 10⁻³ and |ΔLL| < 0.1; otherwise the step is halved and the fit
repeated.

## 5. Estimation

* MLE: 20 random starts, L-BFGS-B with exact gradients inside a stated box
  (per-year rates e⁻²⁰–e², at most e²-fold change per decade, amyloid HR
  e⁻⁵–e⁵). Each start is finished with exact-Hessian Newton steps and counts
  as converged only if the predicted remaining log-likelihood gain (Newton
  decrement) is < 10⁻³ with a positive-definite Hessian (the optimiser's own
  message is not trusted). Estimates on the box
  edge, optima reached by a single start, or a non-positive-definite Hessian
  are reported as such.
* Bayesian (primary analysis): NUTS (NumPyro), 4 chains × 1,000 draws after
  1,000 warm-up, dense mass matrix initialised from the Laplace covariance at
  the MLE, chains started at over-dispersed points N(MLE, 4 × Laplace
  covariance); priors log-rate N(−5, 3), log-slope per decade N(0, 1),
  log HR_A N(0, 1) (centred on *no* amyloid effect). A posterior is reported
  without caveat only if max R-hat ≤ 1.01 and there are zero divergences.
* Model comparison: 5-fold cross-validation by whole age band (bands dealt
  round-robin by age; each fold started from the full-data MLE plus 4 random
  starts); margin 2 log-likelihood units.

## 6. Falsification test (fixed before any fit)

For each variant, the predicted incidence of entering Braak stage **V**
(per 1,000 person-years) is computed at the midpoint of each CFAS age band
(65–69, …, 80–84; 85+ evaluated at age 90). The reference for each band is
the interval between the lowest and highest of the four published rates
(men/women × CFAS I/II). **A variant fails if, in any band, its prediction
lies more than a factor of 2 outside that interval** (|log ratio to the
nearest end| > ln 2).

## 7. Expected amyloid effect (stated before any fit)

From Ossenkoppele 2022 (A+T− vs A−T−, progression to MCI, HR 2.43, 95% CI
1.38–4.26), under the proxy assumption that this HR applies to the LC →
cortex transition (manuscript/derivations/module2_gamma.md):
**log HR_A = 0.888 (0.322–1.449)**, i.e. γ = 1.43 μ0. Comparison:
consistent (0.888 inside the fitted 95% interval) / tension (outside, intervals
overlap) / disagrees (no overlap, flagged). The derivation also predicts
attenuation of the fitted HR_A towards 1 (A3), so a smaller fitted HR_A is
weak evidence against this expectation, but HR_A ≤ 1 is not explained by it.

## 8. Sensitivity analyses (MLE only, all pre-specified)

1. Dementia threshold stage IV instead of V.
2. Open incidence band evaluated at 87.5 instead of 90.
3. Band-average instead of midpoint evaluation (all variants refitted).
4. Variant C with 3 hits.
5. Variant A with the Jansen 2015 amyloid curve (Braak-cohort values outside
   ages 47.5–92.5).

A conclusion in §2 is called robust only if it holds in all five.

## 9. Known limitations, stated in advance

* Cross-sectional autopsy data: the model assumes stage does not change the
  chance of coming to autopsy at a given age (non-informative death).
* Stage V is a neuropathological proxy for clinical dementia; the incidence
  data come from a different population (England) and partly a different era.
* The expected γ rests on a proxy (clinical progression, not tau spread).
* All extracted values are unverified until the researcher checks them.
