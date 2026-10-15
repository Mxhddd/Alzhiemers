# Noradrenergic drug classes and dementia: a pre-registered meta-analysis

Mohammad Hassan, independent researcher · October 2026

This is my plan for pooling published studies of noradrenergic drugs and dementia (Module 5 of my project). I'm writing the rules down before any of the real pooling is done, so that what counts as support and what counts as failure can't shift once I've seen the results. The analysis code and settings are in the project repository (`config/metaanalysis.yaml`), at the commit where I register this file.

## Tools

I used Claude, Anthropic's AI model, through Claude Code to write the analysis code, to find and transcribe the published estimates, and to draft this text. Where I left a methodological setting to it, it chose one and gave its reason in the config file. By registering this plan I adopt those choices as mine.

Claude read the study estimates while transcribing them, so the extraction and the criteria were not done independently. The predictions and criteria below come from the target table in my project brief, which was written before the extraction, not from the extracted numbers. No pooled analysis of any noradrenergic class has been run. The one exception is anticholinergic drugs and dementia, a known association, which I pooled in section 5 to check that the machinery works.

## Change log

| Date | Change |
|---|---|
| 2026-10-13 | First complete draft. |
| 2026-10-14 | **A property of the verdict rule, measured before any real pooling. The method did not change.** While dry-running the pipeline on simulated null data (every true ratio set to 1), one contrast, the dihydropyridine (DHP) calcium-channel-blocker comparison, came out AGAINST at the 95% level. It was only inconclusive at the Bonferroni level. I wanted to know how often that happens by chance, so I ran a calibration simulation (`sim/meta_verdict_calibration.R`; 4,000 data sets per cell; 2–10 studies; τ = 0 or 0.1). With no true effect, REML with Knapp–Hartung gives a false AGAINST in **2.2–3.5%** of contrasts (nominal 2.5%), and an interval that excludes 1 on either side in 4.7–6.9%. The plain z interval gives 4.4% and 8.8% when heterogeneity is present. The ad hoc Knapp–Hartung variant gives 0–2.0%, which is too conservative when there are six studies or fewer. **I am keeping Knapp–Hartung as originally planned.** The consequence, stated now: with 4 predicted contrasts and no true effect, "any contrast AGAINST at 95%" will fire in roughly 10–13% of worlds (1 − (1 − 0.025…0.035)⁴). That is why the Bonferroni-level verdict is always reported next to it. This was measured on simulated data only. The real noradrenergic contrasts have not been pooled. |

## 1. Question

Across published cohorts, do noradrenergic drug classes show a consistent association with dementia incidence that scales with how well the drug reaches the brain, in the directions the locus coeruleus rate-controller hypothesis predicts?

## 2. Predictions, and what would count against them

I am stating these first, before the data section, on purpose.

The confirmatory evidence is the set of **active-comparator contrasts**. When two drugs treat the same condition, comparing them largely cancels out confounding by indication, which is the main weakness of comparing users with non-users. Users-versus-non-users comparisons are pooled and reported, but they do not enter any verdict.

Each pooled ratio is oriented as "drug class versus its comparator":

| Contrast | Target | Prediction (pooled ratio) |
|---|---|---|
| lipophilic vs hydrophilic β-blockers | ADRB1 | < 1 (lower ADRB1 activity is protective) |
| terazosin/doxazosin/alfuzosin vs tamsulosin | ADRA1A | < 1 |
| atomoxetine vs stimulants | SLC6A2 | < 1 |
| brain-penetrant DHP calcium-channel blockers vs amlodipine | CACNA1D | < 1 |
| dual orexin antagonists vs Z-drugs | HCRTR2 | none; exploratory in the brief, reported without a verdict |

For each contrast that has data:

* **Supports:** the 95% confidence interval excludes 1 on the predicted side.
* **AGAINST:** the 95% confidence interval excludes 1 on the other side.
* **Inconclusive:** anything else.

I also report the same three verdicts at a Bonferroni-adjusted level, across the contrasts that have both a prediction and data.

**The result counts against the hypothesis** if any predicted contrast is AGAINST at the 95% level. I will report alongside it whether that verdict survives the Bonferroni level.

**The result counts in its favour** only if every predicted contrast with data is "supports" and none is AGAINST.

If the contrasts are all inconclusive, I will say that there is no evidence either way.

**The known false-alarm rate of this rule** (see the change log, 2026-10-14): with no true effect, a single contrast comes out AGAINST in about 2.2–3.5% of data sets, and "any of the 4 predicted contrasts AGAINST" in about 10–13%. An AGAINST verdict that does not survive the Bonferroni level will be reported as such.

**A limitation I expect, and am stating now:** several contrasts have only one or two studies, and atomoxetine versus stimulants has none in open access. Many verdicts will be inconclusive or will rest on a single study. I will say so when it happens.

## 3. Data

The file is `data/extracted/meta_extraction.csv`.

* It has 37 rows, of which 32 are included.
* Every row records its source (a PMCID), where in the paper the number sits, the verbatim snippet, `verified: false`, and `machine_rechecked`. All 37 rows agree with a second, independent blind re-read.
* The citations, licences, and the searches that found nothing are in `data/extracted/meta_extraction.sources.md`.

**Search.** I searched Europe PMC and NCBI PMC, open-access full text only. No paywalled paper could be read. The searches were run on 2026-10-13 and are logged in `data/extracted/candidates/`. Restricting to open access is a selection rule, and studies behind paywalls are missing from the analysis.

**Selection rules.** Each study contributes one estimate per drug class and comparator type. The rules are coded in `src/metaanalysis/build_extraction.py`.

| Rule | When it applies | Estimate used |
|---|---|---|
| S1 | Always, unless S2–S5 apply | The primary, fully adjusted estimate |
| S2 | Several analyses with equal stated priority | The first one listed |
| S3 | Only per-drug estimates against a shared reference | The one with the narrowest log CI |
| S4 | Only duration or dose categories | The highest category versus none |
| S5 | Several follow-up windows | The longest |
| S6 | Choice of outcome | All-cause dementia, otherwise Alzheimer's disease |

**Exclusion rules.** An excluded row stays in the file with its reason and is counted as a drop.

* **E1:** dementia within Parkinson's disease. Populations restricted by the drug's indication (hypertension, atrial fibrillation) are kept.
* **E2:** no usable standard error. That means a zero-width confidence interval, or one that is asymmetric on the log scale by more than 15% of the half-width plus an allowance for rounding of the printed decimals.
* **E3:** the exposure contrast is confounded by indication in a way unlike the other studies (for example, antidepressant users against people who do not use antidepressants).
* **E4:** anticholinergic rows with an active comparator are kept out of the use-versus-non-use positive control.

**Risk of bias** has not been assessed. All seven ROBINS-I domains are marked `not_assessed`. Once ratings exist, I will add by amendment a sensitivity analysis that excludes studies rated serious or critical overall. The primary analysis does not depend on those ratings.

## 4. Methods

**Scale and pooling**

* Estimates are put on the log scale, and the standard error is worked out from the transcribed confidence interval.
* Rows printed with exposure and comparator the wrong way round (`invert_ratio`) are negated. The transcribed numbers themselves are never edited.
* Random effects: REML, with Hartung–Knapp–Sidik–Jonkman confidence intervals (called Knapp–Hartung elsewhere in this document) and 95% prediction intervals. Fisher-scoring step halving is used to help convergence.
* Pooling needs at least 2 studies. A single study is shown with its own transcribed interval.

**Grouping**

* Hazard ratios, risk ratios and incidence rate ratios are pooled together. Odds ratios are pooled separately. Pooling them together is a sensitivity analysis.
* Active-comparator rows are pooled only within their own contrast.

**Other analyses**

* Subgroups by contrast, comparator type and brain penetrance, each with a meta-regression test for differences between subgroups.
* Publication bias (Egger test and trim-and-fill) only when there are at least 10 studies.
* Leave-one-out only when there are at least 3 studies.

## 5. Checks that must pass first (`make validate-metaanalysis`)

1. **Reproduce a published pooled estimate.** Moraros 2017 (Depress Anxiety), five studies from its Table 3. With DerSimonian–Laird random effects I get 2.169 (1.411–3.334), and the paper reports 2.17 (1.41–3.33). The paper says only "random effects", and DerSimonian–Laird is the default in Stata's `metan`, so I used it. **Passed.**
2. **A three-study example I can compute by hand,** in exact rational arithmetic (`manuscript/derivations/module5_hand_three_study.md`). The machine has checked it. I have not yet checked it myself. **Passed** to 1e-9.
3. **Simulated homogeneous data** should give a median I² near 0. **Passed** (I² ≈ 0).
4. **A positive control:** anticholinergic use versus non-use, hazard ratio, risk ratio and incidence rate ratio group. The pooled ratio is 1.62 (1.26–2.08) from k = 7, harmful, with an interval that excludes 1. **Passed.** The odds-ratio group, reported alongside, gives 1.45 (1.19–1.76) from k = 5. The analysis run repeats this check first and stops if it fails.

## 6. Triangulation (descriptive)

For each Mendelian randomisation target, one figure shows three things on a single log scale, as the effect of lower target activity:

* the pooled active-comparator estimate,
* the Module 1 Mendelian randomisation estimate,
* the mechanistic prediction.

Concordance flags are drawn with them. Missing pieces appear as labelled gaps. Until Module 1 has been run, the Mendelian randomisation column is empty.

## 7. Multiplicity

Each prediction is a separate claim that can fail, read at the 95% level. I report Bonferroni-adjusted verdicts for every contrast. A conclusion that changes between the two levels will be called fragile.

## 8. Sensitivity analyses (pre-specified; `config/metaanalysis.yaml`, key `sensitivity`)

1. DerSimonian–Laird instead of REML.
2. Plain z-based confidence intervals instead of Knapp–Hartung.
3. Odds ratios pooled with the hazard, risk and rate ratios (the rare-outcome approximation).
4. Rows excluded under E1 (Parkinson's populations) put back in.

I will call a verdict robust only if it is the same in all four.

## 9. Limitations, stated in advance

* The studies are observational. Residual confounding remains, especially confounding by indication in the non-user comparisons, and so does protopathic bias close to diagnosis.
* Restricting to open access selects which studies are seen.
* Few studies per contrast means wide intervals or single-study estimates.
* Outcomes mix all-cause dementia and Alzheimer's disease (rule S6).
* Risk of bias has not been assessed.
* Every extracted value is unverified by a human until I have checked it against its source.
