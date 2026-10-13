# Pre-registration — Module 2: Markov base-rate model of tau stage progression

**Status: COMPLETE DRAFT, ready for the researcher to review, upload to OSF
and freeze — not yet frozen.** Written by Claude under the standing policy:
Claude decides methodological questions and states them; the researcher
reviews, uploads to OSF and runs `make freeze-prereg NAME=baserate OSF=<url>`.
Claude must not run the freeze: the agent that writes the criteria should not
be the one that locks them. Nothing in Module 2 has been fitted to real data.
Configuration: `config/baserate.yaml` at the commit that freezes this file.

## Change log

This document has never been frozen. Every change below is a revision of a
draft made before any real data were fitted, not a post-hoc change.

| Date | Change | Reason |
|---|---|---|
| 2026-10-13 | First complete draft: criteria 1 (amyloid gating), 2 (variant comparison) and the incidence falsification test. | — |
| 2026-10-13 | **Falsification test made uncertainty-aware (§6).** A variant now fails only if the whole 95% interval of its predicted incidence is outside the acceptance region, not just its point estimate. Bands too uncertain to fail are reported as uninformative, and the point-estimate verdict is still reported. The factor-of-2 tolerance and the envelope are unchanged. | In the simulated dry run, the point-estimate rule rejected the *true* model. Stage V is almost never observed before age 70, so the IV → V rate at young ages is estimated with very large error (log-rate SE 1.7), and a noisy point prediction in the youngest band (0.03 vs 0.14 per 1,000 person-years) failed the test. A test that rejects the truth through estimation noise is miscalibrated. Decided on simulated data only, before any real fit. |
| 2026-10-13 | **Criterion 1 (amyloid gating) dropped.** Variant A is now the age-adjusted LC-first model. Its amyloid hazard ratio is reported as *descriptive, not a test of gating*. The brief's literal variant A is kept only as a sensitivity analysis labelled *confounded with age*. Decided by the researcher (option (c) of `manuscript/notes/module2_identifiability.md`). | The autopsy stage counts carry no per-brain amyloid status, so the population amyloid curve A(t) enters only as a function of age. In simulation, with no amyloid effect at all, the literal variant A reported HR_A = 2.22 (log-scale 95% CI 0.28–1.31) when the LC → cortex step speeds up by 0.1 per decade, and 7.58 at 0.3 per decade. The criterion would have read ageing as gating and could essentially only come out in favour of the hypothesis. |

**Module 2 does not test amyloid gating.** The stage-by-age counts carry no
per-brain amyloid status, so no model fitted to them can separate "amyloid
speeds the LC → cortex transition" from "the transition speeds up with age".
The amyloid-gating question is deferred to data with per-person amyloid
status and tau stage, e.g. autopsy series cross-tabulating amyloid by Braak
stage by age, or the per-person cohort analyses in Modules 3 and 5.

## 1. Question

Nearly everyone has locus coeruleus (LC) tau by mid-life, yet only a
minority develop Alzheimer's dementia. Which model of the transition from
LC-only tau to cortical tau best explains the age distribution of tau stages
at autopsy? Does it predict the published age-specific dementia incidence,
which it never sees?

* **A, LC first:** cortical tau arises only from the last LC stage.
* **B, independent:** cortical tau arises from any pre-cortical stage,
  unrelated to LC stage.
* **C, multi-hit:** cortical entry needs m independent hits.

## 2. What would count against the LC-first model (stated first)

1. *(Dropped 2026-10-13 — see the change log. Amyloid gating is not tested
   here.)*
2. **Which variant explains the stage distribution (5-fold cross-validation
   by whole age band, held-out multinomial log-likelihood):**
   * *supports LC first* — A beats both B and C by more than 2 held-out
     log-likelihood units;
   * *against LC first* — B or C beats A by more than 2 units;
   * *inconclusive* — anything else; reported as uninformative.

   Known asymmetry, measured on simulated data before any real fit
   (`sim/baserate_power.py`; results in
   `manuscript/notes/module2_identifiability.md`): with ~2,300 brains, B or
   C worlds are much easier to tell from A than A worlds are from B. The
   power of each outcome is reported with the result.
3. **Falsification (§6):** a variant whose predicted incidence fails the test
   is rejected whatever criterion 2 says. If every variant fails, the
   framework itself is inadequate.

A result counts **against the LC-first model** if criterion 2 is "against"
or if A fails §6 while B or C passes. It counts **in favour** only if
criterion 2 "supports" and A passes §6.

## 3. Data (all extracted by Claude, `verified: false` until a human checks them)

| Input | Source | File |
|---|---|---|
| Tau stage counts by age band (stages 0, a, b, c, 1a, 1b, I–VI) | Braak H, Thal DR, Ghebremedhin E, Del Tredici K. *J Neuropathol Exp Neurol* 2011;70:960–969 | `data/extracted/braak2011_stage_counts.csv` (**pending: full text not openly available**) |
| Amyloid-plaque prevalence by age band, same cohort (A(t) for the descriptive HR) | Braak et al. 2011 | `data/extracted/braak2011_amyloid_by_age.csv` (**pending**) |
| Amyloid-positivity prevalence, cognitively normal (sensitivity A(t)) | Jansen WJ et al. *JAMA* 2015;313:1924–1938, Table 2 | `data/extracted/amyloid_prevalence_jansen2015.csv` |
| Dementia incidence by age (falsification only) | Matthews FE et al. *Nat Commun* 2016;7:11398, Table 1 | `data/extracted/dementia_incidence_cfas.csv` |
| Published amyloid hazard ratio (descriptive reference only) | Ossenkoppele R et al. *Nat Med* 2022;28:2381–2387, Table 2 | `data/extracted/amyloid_hr_ossenkoppele2022.csv` |

Consistency check before fitting: the extracted Braak table must reproduce
the abstract totals exactly. Those totals are 2,332 brains; 10 without AT8
tau; 58 at stages a–c; 38 at 1a; 236 at 1b; 1,031 with plaques. A mismatch
stops the analysis.

## 4. Model

Continuous-time Markov chain, dp/dt = p Q(t). Stages follow the Braak 2011
order 0 → a → b → c → 1a → 1b → I → … → VI. Everyone is at stage 0 at age 1,
the youngest age in the series. Progression rates are exp(α_k + β_k (t − 50)/10).

* **A, LC first (age-adjusted):** c → 1a at rate
  exp(α_μ + β_μ (t − 50)/10) · (1 + (HR_A − 1) A(t)). The age slope β_μ is
  what keeps ageing from being read as amyloid. HR_A is *descriptive, not a
  test of gating*. A nests the plain LC-first chain (HR_A = 1).
* **B, independent:** every pre-cortical stage enters 1a at a common rate
  exp(a + b(t − 50)/10).
* **C, multi-hit:** 1a is entered after m = 2 independent hits at constant
  rate (sensitivity: m = 3).

Likelihood: multinomial counts per age band, with stage probabilities at the
band midpoint. Numerical solution: exponential midpoint rule with step 0.25
years (manuscript/derivations/module2_step_size.md). At every fitted
estimate, the derived error bound h²Φ(θ̂) must be ≤ 10⁻³. A Radau audit must
also give max |Δp| < 10⁻³ and |ΔLL| < 0.1. Otherwise the step is halved and
the fit repeated.

## 5. Estimation

* **MLE.**
  * 20 random starts, L-BFGS-B with exact gradients, inside a stated box:
    per-year rates e⁻²⁰ to e², at most an e²-fold change per decade, and
    amyloid HR e⁻⁵ to e⁵.
  * Each start is finished with exact-Hessian Newton steps. It counts as
    converged only if the Newton decrement is < 10⁻³ and the Hessian is
    positive definite.
  * Estimates on the box edge, optima reached by a single start, or a
    non-positive-definite Hessian are reported as such.
* **Bayesian.**
  * NUTS (NumPyro): 4 chains × 1,000 draws after 1,000 warm-up.
  * Dense mass matrix initialised from the Laplace covariance at the MLE;
    chains start at over-dispersed points N(MLE, 4 × Laplace covariance).
  * Priors: log-rate N(−5, 3), log-slope per decade N(0, 1), log HR_A N(0, 1).
  * A posterior is reported without caveat only if max R-hat ≤ 1.01 and
    there are zero divergences.
* **Model comparison.**
  * 5-fold cross-validation by whole age band; bands are dealt round-robin
    by age.
  * Each fold starts from the full-data MLE plus 4 random starts.
  * Margin: 2 log-likelihood units.

## 6. Falsification test (fixed before any real fit)

For each variant, compute the predicted incidence of entering Braak stage
**V**, per 1,000 person-years. It is evaluated at the midpoint of each CFAS
age band (65–69, …, 80–84; 85+ is evaluated at age 90). The reference for
each band is the interval between the lowest and highest of the four
published rates (men/women × CFAS I/II). The *acceptance region* is that
interval widened by a factor of 2 on both sides.

Predictions are computed for 200 parameter draws: posterior draws for the
primary fits, and Laplace draws N(MLE, inverse observed information) for the
MLE-only sensitivity fits.

* **A variant fails if, in any band, the whole 95% interval of its predicted
  incidence lies outside the acceptance region**: the nearer end of the
  interval is more than a factor of 2 outside the published interval.
* A band whose 95% interval is wider than the whole acceptance region is
  reported as *uninformative* (it could not have failed). If every band is
  uninformative, the verdict is UNINFORMATIVE, not PASS.
* The point-estimate rule of the first draft is reported alongside: fail if
  the point prediction is more than a factor of 2 outside.

## 7. Amyloid hazard ratio (descriptive only)

Reported as *descriptive, not a test of gating*:

* the posterior for HR_A from the age-adjusted variant A;
* the same under the Jansen 2015 amyloid curve;
* next to the published reference from Ossenkoppele 2022 (A+T− vs A−T−,
  progression to MCI): HR 2.43 (95% CI 1.38–4.26), i.e. log HR 0.888
  (0.322–1.449) (manuscript/derivations/module2_gamma.md).

The literal variant A's HR is reported in the sensitivity table labelled
*confounded with age*. None of these numbers enters any verdict.

## 8. Sensitivity analyses (MLE only, all pre-specified)

1. Dementia threshold stage IV instead of V.
2. Open incidence band evaluated at 87.5 instead of 90.
3. Band-average instead of midpoint evaluation (all variants refitted).
4. Variant C with 3 hits.
5. Variant A with the Jansen 2015 amyloid curve (Braak-cohort values outside
   ages 47.5–92.5).
6. The brief's literal variant A, labelled *confounded with age*.

A conclusion in §2 is called robust only if it holds in 1–5 (6 is shown for
transparency, not used for robustness).

## 9. Known limitations, stated in advance

* Cross-sectional autopsy data: the model assumes stage does not change the
  chance of coming to autopsy at a given age (non-informative death).
* Stage V is a neuropathological proxy for clinical dementia. The incidence
  data come from a different population (England) and partly a different era.
* With ~2,300 brains, the variant comparison has limited power to support A
  even when A is true (§2). An inconclusive result is the expected outcome
  in that case and is not evidence against A.
* Amyloid gating cannot be tested with these data (change log).
* All extracted values are unverified until the researcher checks them.
