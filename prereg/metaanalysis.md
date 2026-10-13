# Pre-registration — Module 5: meta-analysis and triangulation of noradrenergic drug classes

**Status: COMPLETE DRAFT, ready for the researcher to review, upload to OSF
and freeze — not yet frozen.** Written by Claude under the standing policy.
Claude must not run `make freeze-prereg`: the agent that writes the
criteria should not be the one that locks them. Configuration:
`config/metaanalysis.yaml` at the commit that freezes this file.

**What has and has not been seen.** The study-level estimates were
extracted by Claude, and extracting them means reading them. No pooled
analysis of the noradrenergic classes has been run: the run refuses
without this file frozen. The predictions and criteria below come from the
brief's target table, not from the extracted numbers. The anticholinergic
positive control has been pooled as part of validation (§5). It is a known
association, not a question this module asks.

## Change log

| Date | Change |
|---|---|
| 2026-10-13 | First complete draft. |

## 1. Question

Across published cohorts, do noradrenergic drug classes show a consistent,
brain-penetrance-graded association with dementia incidence, in the
directions the LC rate-controller hypothesis predicts?

## 2. Predictions and what counts against them (stated first)

The confirmatory evidence is the **active-comparator contrasts**. Comparing
two treatments for the same indication largely cancels confounding by
indication. Users-vs-non-users comparisons are pooled and reported but do
not enter the verdicts.

For each contrast, the pooled ratio is oriented as "drug class vs its
comparator":

| Contrast | Target | Prediction (pooled ratio) |
|---|---|---|
| lipophilic vs hydrophilic β-blockers | ADRB1 | < 1 (lower ADRB1 activity protective) |
| terazosin/doxazosin/alfuzosin vs tamsulosin | ADRA1A | < 1 |
| atomoxetine vs stimulants | SLC6A2 | < 1 |
| brain-penetrant DHP calcium-channel blockers vs amlodipine | CACNA1D | < 1 |
| dual orexin antagonists vs Z-drugs | HCRTR2 | none: exploratory in the brief, reported without a verdict |

For each contrast with data:

- **supports**: the 95% CI excludes 1 on the predicted side;
- **AGAINST**: the 95% CI excludes 1 on the other side;
- **inconclusive**: otherwise.

The same verdicts are also reported at a Bonferroni-adjusted level over
the contrasts that have both a prediction and data.

**The result counts against the hypothesis** if any predicted contrast is
AGAINST at the 95% level. Whether it survives the Bonferroni level is
reported alongside. **It counts in its favour** only if every predicted
contrast with data is "supports" and none is AGAINST. A set of
inconclusive contrasts is reported as no evidence either way.

Expected limitation, stated now: several contrasts have one or two studies
(atomoxetine vs stimulants has none in open access). Many verdicts will be
inconclusive or rest on a single study, and this will be said.

## 3. Data

`data/extracted/meta_extraction.csv`:

- 37 rows, 32 included;
- every row has its source (PMCID), locator, verbatim snippet,
  `verified: false` and `machine_rechecked` (all 37 agree with an
  independent blind re-read);
- citations, licences and the searches that found nothing are in
  `data/extracted/meta_extraction.sources.md`.

**Search.** Europe PMC and NCBI PMC, open-access full text only (no
paywalled papers could be read). Searches were run on 2026-10-13 and are
logged in `data/extracted/candidates/`. *Limitation:* restricting to open
access is a selection rule. Studies behind paywalls are missing.

**Selection rules** (one estimate per study per class and comparator type;
coded in `src/metaanalysis/build_extraction.py`):

| Rule | When it applies | Estimate used |
|---|---|---|
| S1 | Always, unless S2–S5 apply | The primary fully adjusted estimate |
| S2 | Several analyses with equal stated priority | The first listed |
| S3 | Only per-drug estimates against a shared reference | The narrowest log CI |
| S4 | Only duration or dose categories | The highest category vs none |
| S5 | Several follow-up windows | The longest |
| S6 | Outcome choice | All-cause dementia, else Alzheimer's disease |

**Exclusion rules** (each row is kept in the file with its reason and
logged as a drop):

- **E1:** dementia within Parkinson's disease. Populations restricted by
  the drug's indication (hypertension, atrial fibrillation) are kept.
- **E2:** no usable standard error. That means a zero-width CI, or a CI
  asymmetric on the log scale beyond 15% of the half-width plus the
  rounding allowance for the printed decimals.
- **E3:** the exposure contrast is confounded by indication unlike the
  other studies (antidepressant users vs non-users of antidepressants).
- **E4:** anticholinergic rows with an active comparator are kept out of
  the use-vs-non-use positive control.

**Risk of bias:** not assessed. All seven ROBINS-I domains are
`not_assessed`. A sensitivity analysis excluding studies rated
serious/critical overall will be added by amendment once ratings exist.
The primary analysis does not depend on them.

## 4. Methods

**Scale and pooling.**
- Estimates are put on the log scale. The SE is backed out of the
  transcribed CI.
- Rows printed with exposure and comparator reversed (`invert_ratio`) are
  negated. The transcribed numbers are never edited.
- Random effects: REML with Hartung–Knapp–Sidik–Jonkman CIs and 95%
  prediction intervals. Fisher-scoring step halving is used for convergence.
- At least 2 studies are needed to pool. A single study is shown with its
  own transcribed CI.

**Grouping.**
- HR, RR and IRR are pooled together. ORs are pooled separately (pooling
  them together is a sensitivity analysis).
- Active-comparator rows are pooled only within their contrast.

**Other analyses.**
- Subgroups by contrast, comparator type and brain penetrance, each with
  a meta-regression test of subgroup differences.
- Publication bias (Egger test and trim-and-fill) only with ≥ 10 studies.
- Leave-one-out only with ≥ 3 studies.

## 5. Gates that must pass first (validation, `make validate-metaanalysis`)

1. **Published pooled estimate reproduced:** Moraros 2017 (Depress
   Anxiety), five studies from Table 3.
   - DerSimonian–Laird random effects: ours 2.169 (1.411–3.334) vs
     published 2.17 (1.41–3.33).
   - The estimator is metan's default; the paper says only "random
     effects". PASS.
2. **Hand-computable three-study example** in exact rational arithmetic
   (`manuscript/derivations/module5_hand_three_study.md`; machine-checked,
   not human-checked). PASS to 1e-9.
3. **Simulated homogeneous data:** median I² ≈ 0. PASS.
4. **Positive control:** anticholinergic use vs non-use, HR/RR/IRR group.
   - Pooled ratio 1.62 (1.26–2.08), k = 7, harmful with a CI excluding 1. PASS.
   - The OR group, reported alongside: 1.45 (1.19–1.76), k = 5.
   - The analysis run repeats this check first and stops if it fails.

## 6. Triangulation (descriptive)

For each MR target, the figure shows on one log scale, as the effect of
*lower* target activity:

- the pooled active-comparator estimate;
- the Module 1 MR estimate;
- the mechanistic prediction.

Concordance flags are shown with them. Missing pieces are drawn as
labelled gaps. Until Module 1 runs, the MR column is empty.

## 7. Multiplicity

Each prediction is a separate falsifiable claim, read at 95%. The
Bonferroni-adjusted verdicts are reported for every contrast. A conclusion
that changes between the two levels is reported as fragile.

## 8. Sensitivity analyses (pre-specified; `config/metaanalysis.yaml` → `sensitivity`)

1. DerSimonian–Laird instead of REML.
2. Plain z-based CIs instead of Knapp–Hartung.
3. ORs pooled with HR/RR/IRR (rare-outcome approximation).
4. Rows excluded under E1 (Parkinson's populations) put back.

A verdict is called robust only if it is the same in all four.

## 9. Limitations, stated in advance

- Observational studies: residual confounding, especially by indication in
  non-user comparisons, and protopathic bias near diagnosis.
- Open-access restriction: selection of studies.
- Few studies per contrast: wide or single-study intervals.
- Outcomes mix all-cause dementia and Alzheimer's disease (rule S6).
- Risk of bias not assessed.
- All extracted values unverified by a human until checked.
