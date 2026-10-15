# How tau spreads from the locus coeruleus to the cortex: a pre-registered base-rate model

Mohammad Hassan, independent researcher · October 2026

This is my plan for fitting a simple model of tau progression to the 2,332 autopsy brains in Braak et al. (2011) (Module 2 of my project). I'm writing it down before anything is fitted to the real counts, so that what counts as support and what counts as failure can't move once I've seen the results. Nothing in this module has been fitted to real data yet. The code and settings are in the project repository (`config/baserate.yaml`), at the commit where I register this file.

## Tools

I used Claude, Anthropic's AI model, through Claude Code to write the code, extract the published numbers, run the simulations and draft this text. I told it to decide methodological questions for me and explain each decision, which it did in the config file and in the change log below. By registering this plan I adopt those decisions as mine, including the one I made myself: dropping criterion 1. The analysis code refuses to run on real data until this plan is frozen.

## Change log

This document has never been frozen. Every change below is a revision of a draft, made before any real data were fitted. None of them is a change made after seeing a result.

| Date | Change | Reason |
|---|---|---|
| 2026-10-13 | First complete draft: criterion 1 (amyloid gating), criterion 2 (variant comparison) and the incidence falsification test. | — |
| 2026-10-13 | **The falsification test now takes uncertainty into account (section 6).** A variant fails only if the whole 95% interval of its predicted incidence lies outside the acceptance region, not just its point estimate. Bands that are too uncertain to fail are reported as uninformative. The point-estimate verdict is still reported. The factor-of-2 tolerance and the envelope are unchanged. | In the simulated dry run, the point-estimate rule rejected the *true* model. Stage V is almost never seen before age 70, so the IV → V rate at young ages is estimated with very large error (log-rate SE 1.7). A noisy point prediction in the youngest band (0.03 against 0.14 per 1,000 person-years) failed the test. A test that rejects the truth through estimation noise is miscalibrated. I decided this on simulated data only. |
| 2026-10-13 | **Criterion 1 (amyloid gating) dropped.** Variant A is now the age-adjusted LC-first model. Its amyloid hazard ratio is reported as *descriptive, not a test of gating*. The brief's literal variant A is kept only as a sensitivity analysis, labelled *confounded with age*. I made this decision (option (c) of `manuscript/notes/module2_identifiability.md`). | The autopsy stage counts carry no per-brain amyloid status, so the population amyloid curve A(t) enters only as a function of age. In simulation with no amyloid effect at all, the literal variant A reported HR_A = 2.22 (log-scale 95% CI 0.28–1.31) when the LC → cortex step speeds up by 0.1 per decade, and 7.58 at 0.3 per decade. The criterion would have read ageing as amyloid gating, and could in practice only come out in favour of the hypothesis. |
| 2026-10-13 | **Labels for the descriptive amyloid hazard ratio (section 7).** An estimate sitting on its parameter bound is labelled *not identified*. An interval spanning more than a 100-fold range of hazard ratio (log width above 4.605, `gamma_check.uninformative_log_width`) is labelled *uninformative*, and is not compared with the published value. | In the dry run, the Jansen sensitivity estimate sat on its bound (log HR −5, interval ±47), yet the table said "published value inside the interval". That reads as agreement when the data said nothing. This changes labelling only. No verdict uses these numbers. |
| 2026-10-14 | **Criterion 2: cross-validation now runs over brains within age bands, not over whole age bands. Fold fits that cannot be verified no longer stop the run.** Whole-band cross-validation is kept as sensitivity analysis 7. | Found by the power study, on simulated data only. With whole bands held out, a stage seen in only one band (in the simulated world, stage VI only at ages 91–100) is missing from that fold's training data. Every model then gives it a probability near 10⁻¹² (about −730 log-likelihood units in one band), and the verdict is decided by numerical detail in that band. Refitting such a replicate with the full multi-start procedure flipped it from AGAINST to supports. In 2 of 4 refitted replicates, a fold fit had no start that passed the convergence test (a flat direction in the training subset), and the code stopped instead of "reporting as such" (section 5). With brains dealt to the folds within each band, every training set covers every band and, in practice, every stage. I read "cross-validation by age band" in the brief as *stratified by* age band. |
| 2026-10-14 | **Data section (section 3) filled in from the paper's own supplements.** The stage-by-age counts are in Supplemental Digital Content 1 (Supplementary Table 1). Amyloid phase by age is in SDC 2. Both are openly downloadable from the publisher. A(t) is defined as the share of brains in amyloid phases 1–4 (any plaque) in each age band. The consistency check before fitting now also requires the band totals to match Table 1, and the plaque total to match the abstract (1,031). | The full text and supplements became available (I supplied the PDF, and the supplement permalinks are open). "Any plaque" matches the abstract's own definition ("1,031 (44.2%) had β-amyloid plaques"). The values were extracted with `verified: false`, and 43/43 count-consistency checks pass (`data/extracted/braak2011_sum_checks.md`). No model has been fitted to these counts. |

**Module 2 does not test amyloid gating.** The stage-by-age counts carry no per-brain amyloid status, so no model fitted to them can tell apart "amyloid speeds up the LC → cortex transition" from "the transition speeds up with age". That question is left for data with per-person amyloid status and tau stage, for example an autopsy series that cross-tabulates amyloid by Braak stage by age, or the per-person cohort analyses in Modules 3 and 5.

## 1. Question

Nearly everyone has locus coeruleus (LC) tau by mid-life, yet only a minority develop Alzheimer's dementia. Which model of the step from LC-only tau to cortical tau best explains the age distribution of tau stages at autopsy? And does that model predict the published age-specific dementia incidence, which it never sees?

* **A, LC first:** cortical tau arises only from the last LC stage.
* **B, independent:** cortical tau arises from any pre-cortical stage, unrelated to LC stage.
* **C, multi-hit:** cortical entry needs m independent hits.

## 2. What would count against the LC-first model

I am stating these first, on purpose.

1. *Dropped on 2026-10-13; see the change log. Amyloid gating is not tested here.*
2. **Which variant best explains the stage distribution.** I use 5-fold cross-validation over brains within age bands, scored by the held-out multinomial log-likelihood.
   * *Supports LC first:* A beats both B and C by more than 2 held-out log-likelihood units.
   * *Against LC first:* B or C beats A by more than 2 units.
   * *Inconclusive:* anything else. It is reported as uninformative.

   I measured the operating characteristics on simulated data before any real fit (`sim/baserate_power.py`; results in `manuscript/notes/module2_identifiability.md`), and I will report them alongside the result. In 100 simulated data sets of the real size (20 per world; light single-start fits, step 0.5 years):
   * When LC-first is true, with an age slope on the LC → cortex step of 0, 0.1 or 0.3 per decade and no amyloid effect, the verdict was "supports" in 57 of 60 (95%; Wilson 95% CI 86–98%) and "against" in 1 of 60 (1.7%; 0.3–8.9%).
   * When B is true, the verdict was "against" in 20 of 20.
   * When C is true, it was "against" in 16 of 20 and inconclusive in 4. "Supports" never occurred when B or C was true (0 of 40; upper 95% bound 8.8%).

   These numbers describe simulated worlds, built around the dry-run parameter values. They are not guarantees for the real data.
3. **Falsification (section 6).** A variant whose predicted incidence fails that test is rejected whatever criterion 2 says. If every variant fails, the framework itself is inadequate.

A result counts **against the LC-first model** if criterion 2 is "against", or if A fails section 6 while B or C passes. It counts **in favour** only if criterion 2 is "supports" and A passes section 6.

## 3. Data

Everything below was extracted with the help of Claude and is marked `verified: false` until I have checked it against the source.

| Input | Source | File |
|---|---|---|
| Tau stage counts by age band (12 stages: 0, a, b, c, 1a, 1b, I–VI; 10 bands, ages 1–10 to 91–100; 2,332 brains) | Braak H, Thal DR, Ghebremedhin E, Del Tredici K. *J Neuropathol Exp Neurol* 2011;70:960–969, **Supplemental Digital Content 1** (Supplementary Table 1), http://links.lww.com/NEN/A273 | `data/extracted/braak2011_stage_counts.csv` |
| Amyloid phase 0–4 by age band, same cohort. A(t) = brains in phases 1–4 / all brains in the band (any plaque). Used only for the descriptive hazard ratio | Braak et al. 2011, **Supplemental Digital Content 2** (Supplementary Table 2), http://links.lww.com/NEN/A274 | `data/extracted/braak2011_sdc2_amyloid_by_age.csv` |
| Amyloid-positivity prevalence in cognitively normal people (sensitivity A(t)) | Jansen WJ et al. *JAMA* 2015;313:1924–1938, Table 2 | `data/extracted/amyloid_prevalence_jansen2015.csv` |
| Dementia incidence by age (falsification test only) | Matthews FE et al. *Nat Commun* 2016;7:11398, Table 1 | `data/extracted/dementia_incidence_cfas.csv` |
| Published amyloid hazard ratio (descriptive reference only) | Ossenkoppele R et al. *Nat Med* 2022;28:2381–2387, Table 2 | `data/extracted/amyloid_hr_ossenkoppele2022.csv` |

**Consistency check before fitting.** Any mismatch stops the analysis.

* The stage table must reproduce the abstract totals exactly: 2,332 brains; 10 without AT8 tau; 58 at stages a–c; 38 at 1a; 236 at 1b.
* The amyloid table must give 1,031 brains with plaques (2,332 minus phase 0).
* Every age band's total must agree between the stage table and the amyloid table.

These checks also passed at extraction, together with the row and column sums of the paper's Table 1, Table 2 and both supplements (`data/extracted/braak2011_sum_checks.md`).

**What these data cannot do.** Neither the paper nor its supplements cross amyloid with tau stage *and* age for each brain.

* Table 2 crosses amyloid with tau stage, but has no age.
* SDC 1 crosses tau stage with age.
* SDC 2 crosses amyloid with age.

That is why amyloid gating is not tested here (see the change log).

## 4. Model

A continuous-time Markov chain, dp/dt = p Q(t). The stages follow the Braak 2011 order 0 → a → b → c → 1a → 1b → I → … → VI. Everyone is at stage 0 at age 1, the youngest age in the series. Progression rates are exp(α_k + β_k (t − 50)/10).

* **A, LC first (age-adjusted):** the c → 1a step runs at rate exp(α_μ + β_μ (t − 50)/10) · (1 + (HR_A − 1) A(t)). The age slope β_μ is what stops ageing being read as amyloid. HR_A is *descriptive, not a test of gating*. A contains the plain LC-first chain as the special case HR_A = 1.
* **B, independent:** every pre-cortical stage enters 1a at a common rate exp(a + b(t − 50)/10).
* **C, multi-hit:** 1a is entered after m = 2 independent hits at a constant rate (sensitivity: m = 3).

The likelihood is multinomial, using counts per age band, with stage probabilities taken at the band midpoint. I solve the model numerically with an exponential midpoint rule and a step of 0.25 years (`manuscript/derivations/module2_step_size.md`). At every fitted estimate, the derived error bound h²Φ(θ̂) must be at most 10⁻³. A Radau audit must also give max |Δp| < 10⁻³ and |ΔLL| < 0.1. If either fails, the step is halved and the fit repeated.

## 5. Estimation

**Maximum likelihood**

* 20 random starts, L-BFGS-B with exact gradients, inside a stated box: per-year rates e⁻²⁰ to e², at most an e²-fold change per decade, and amyloid HR e⁻⁵ to e⁵.
* Each start is finished with exact-Hessian Newton steps. It counts as converged only if the Newton decrement is below 10⁻³ and the Hessian is positive definite.
* Estimates on the edge of the box, optima reached by a single start, and a Hessian that is not positive definite are all reported as such.

**Bayesian**

* NUTS (NumPyro): 4 chains of 1,000 draws after 1,000 warm-up.
* The dense mass matrix starts from the Laplace covariance at the maximum-likelihood estimate. Chains start at over-dispersed points N(MLE, 4 × Laplace covariance).
* Priors: log-rate N(−5, 3), log-slope per decade N(0, 1), log HR_A N(0, 1).
* Target acceptance probability 0.9.
* A posterior is reported without a caveat only if the largest R-hat is at most 1.01 and there are zero divergences.

**Model comparison**

* 5-fold cross-validation over brains, stratified by age band. Within every band, the brains are shuffled with a fixed seed (20261014) and dealt to the 5 folds in turn. Each fold is fitted to 4/5 of every band and scored on the held-out 1/5.
* Each fold starts from the full-data maximum-likelihood estimate plus 4 random starts.
* If no start in a fold passes the convergence test, the best start by training log-likelihood is used. That fit is counted and reported as unverified, and the run does not stop.
* Margin: 2 log-likelihood units.

## 6. Falsification test (fixed before any real fit)

For each variant I compute the predicted incidence of entering Braak stage **V**, per 1,000 person-years. It is evaluated at the midpoint of each CFAS age band (65–69, …, 80–84), and the 85+ band is evaluated at age 90. The reference for each band is the interval between the lowest and highest of the four published rates (men and women, CFAS I and II). The *acceptance region* is that interval widened by a factor of 2 on both sides.

Predictions are computed for 200 parameter draws. These are posterior draws for the primary fits, and Laplace draws N(MLE, inverse observed information) for the sensitivity fits that use the maximum-likelihood estimate only.

* **A variant fails if, in any band, the whole 95% interval of its predicted incidence lies outside the acceptance region.** In other words, the nearer end of the interval is more than a factor of 2 outside the published interval.
* A band whose 95% interval is wider than the whole acceptance region is reported as *uninformative*, because it could not have failed. If every band is uninformative, the verdict is UNINFORMATIVE, not PASS.
* The point-estimate rule from the first draft is reported alongside: a variant fails if its point prediction is more than a factor of 2 outside.

## 7. Amyloid hazard ratio (descriptive only)

This is reported as *descriptive, not a test of gating*:

* the posterior for HR_A from the age-adjusted variant A,
* the same result under the Jansen 2015 amyloid curve,
* both next to the published reference from Ossenkoppele 2022 (A+T− versus A−T−, progression to MCI): HR 2.43 (95% CI 1.38–4.26), that is, log HR 0.888 (0.322–1.449) (`manuscript/derivations/module2_gamma.md`).

The literal variant A's hazard ratio goes in the sensitivity table, labelled *confounded with age*. An estimate on its parameter bound is labelled *not identified*. An interval spanning more than a 100-fold range of hazard ratio is labelled *uninformative* (see the change log). None of these numbers enters any verdict.

## 8. Sensitivity analyses (maximum likelihood only; all pre-specified)

1. Dementia threshold at stage IV instead of V.
2. Open incidence band evaluated at 87.5 instead of 90.
3. Band-average instead of midpoint evaluation (all variants refitted).
4. Variant C with 3 hits.
5. Variant A with the Jansen 2015 amyloid curve (Braak-cohort values outside ages 47.5–92.5).
6. The brief's literal variant A, labelled *confounded with age*.
7. Criterion 2 with whole age bands held out in turn (the first design). It is reported with its known weakness (see the change log): a stage seen only in a held-out band gets near-zero probability from every model.

I will call a conclusion in section 2 robust only if it holds in analyses 1–5. Analyses 6 and 7 are shown for transparency and are not used to judge robustness.

## 9. Known limitations, stated in advance

* The data are cross-sectional autopsy data. The model assumes that stage does not change the chance of coming to autopsy at a given age (non-informative death).
* Stage V is a neuropathological stand-in for clinical dementia. The incidence data come from a different population (England) and partly from a different era.
* With about 2,300 brains, the variant comparison has limited power to support A even when A is true (section 2). An inconclusive result is the expected outcome in that case, and it is not evidence against A.
* Amyloid gating cannot be tested with these data (see the change log).
* Every extracted value is unverified until I have checked it against its source.
