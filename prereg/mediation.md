# Amyloid, social and sensory risk factors and cortical tau stage at autopsy: a pre-registered analysis of the Rush ROS and MAP cohorts

Mohammad Hassan, independent researcher · October 2026

This is my plan for Module 3 of my project, the test of whether the locus coeruleus (LC) carries the effect of modifiable risk factors on tau pathology ("Killer 2"). I have requested the Standard dataset of the Religious Orders Study (ROS) and the Rush Memory and Aging Project (MAP) from the Rush Alzheimer's Disease Center (RADC). It has not arrived, and nobody working on this project has seen any of its values. I am writing the plan down now, so that what counts as support and what counts as failure is fixed before I can look. The code and settings are in the project repository (`config/mediation.yaml`, `src/mediation/`), at the commit where I register this file.

**What this dataset cannot do, stated first.** The Standard dataset's tau stage, `braaksc`, is the classic Braak neurofibrillary tangle stage, 0 to VI, from a silver stain, starting in the entorhinal cortex. It has no LC or subcortical pre-tangle stages (a–c, 1a, 1b). A brain at stage 0 may or may not have LC tau. So **these data cannot test the LC → cortex step itself.** They can test whether amyloid and age track cortical tau stage, and whether the social and sensory variables I chose track it. Mediation through the LC is tested only if RADC confirms that LC measures are in the data I receive (section 5.5); otherwise Killer 2 stays untested and I will say so.

## Tools

I used Claude, Anthropic's AI model, through Claude Code to write the code, read RADC's public documentation pages, extract the published values used to validate the code, run the simulations and draft this text. I told it to make the methodological choices for me and give the reason for each, which it did in the config file. By registering this plan I adopt those choices as mine. Claude did not request, download, open or analyse any RADC data, and used no RADC account. The analysis code refuses to run on real data until this plan is frozen, until the code's validation has passed, and until I have confirmed every variable name against the codebook that comes with the data.

## Change log

This document has never been frozen. Every change below was made before any RADC data existed on my side.

| Date | Change | Reason |
|---|---|---|
| 2026-10-12 | First version: a mediation module with the ρ sensitivity analysis, built on simulated data. | Project brief, session 6. |
| 2026-10-13 | First complete draft: hearing loss as the primary exposure, neuromelanin-sensitive LC MRI as the mediator, tau-PET as the outcome; the cohort left for me to name. | The brief asks for a cohort with LC measures, risk factors and tau; I did not have one. |
| 2026-10-21 | **Rewritten for the RADC ROS/MAP Standard dataset.** Outcome: Braak stage at autopsy (ordinal models). Exposures: amyloid load, and the four social and sensory variables I named (`soc_net`, `social_isolation`, `late_life_soc_act`, `vision`). Mediation through the LC becomes a conditional branch, decided by what RADC delivers, not by results. Hearing loss and LC MRI are dropped. | This is the dataset I requested. Its public documentation lists no LC measure and no hearing measure. |
| 2026-10-21 | **Validation now reproduces two published analyses** from the authors' public replication data (section 7), and checks the new ordinal estimators against simulated truth. | The published-reproduction check had been waiting for values; it now passes. |
| 2026-10-21 | **Operating characteristics measured on simulated data** (section 6). One property is stated here because it shapes how results can be read: in every simulated world, the "practically null" verdict was never reached. | The sample sizes I expect are too small for it. A null H3 will read "inconclusive". |

## 1. Questions

1. Does amyloid load track cortical tau stage in the same person, once age is accounted for? In my base-rate model (Module 2) amyloid and age could not be told apart, because the autopsy counts carried no per-person amyloid status. These data have it.
2. Does amyloid act at **cortical entry** (Braak stage I or above versus 0)? This is the split closest to the LC → cortex step that a silver-stain Braak stage can see. If amyloid gates the passage of tau from the LC into the cortex, it should show here. If cortical entry is driven by age alone, as the "primary age-related tauopathy" view holds, it should not.
3. Do smaller social networks, loneliness, less late-life social activity and poorer vision go with a higher Braak stage? The LC-first hypothesis says these act on tau through the LC, so they should at least go with more tau.
4. Only if LC measures are delivered: is part of that association carried by the LC (Killer 2)?

## 2. Hypotheses, and what would count against them

I am stating these before the data section, on purpose. Every association is a proportional-odds odds ratio per standard deviation (SD), oriented so that a value above 1 is the predicted direction. Verdicts use Bonferroni-adjusted intervals within each family: Family A is H1 and H2 (2 tests, 97.5% intervals); Family B is H3 for the four exposures (4 tests, 98.75% intervals). Unadjusted 95% intervals are reported alongside.

For each test:

* **Supports:** the adjusted interval lies above 1.
* **AGAINST:** the adjusted interval lies below 1.
* **AGAINST (practically null):** the adjusted interval lies entirely inside 0.90 to 1.11 (|log OR| < 0.10536, an odds change under 10% per SD, too small to carry a rate-controller claim).
* **Uninformative:** none of the above, and the adjusted interval spans more than a 4-fold range of odds ratio (log width above 1.386).
* **Inconclusive:** anything else.

**H1 (amyloid → Braak stage).** Higher amyloid load goes with a higher Braak stage, given age at death, sex, education, APOE ε4 carrier status and cohort. *Against:* AGAINST or practically null. That would mean amyloid does not track cortical tau in the same person, and it would count against the amyloid-gated model (Module 2, variant A) and for the independent model (variant B).
*This test can hardly fail, and it does not discriminate.* The amyloid cascade predicts the same association, and in simulation H1 came out "supports" in every world where amyloid speeds cortical spread, including a world where amyloid has no effect on cortical entry (section 6). Support for H1 is not evidence for LC-first.

**H2 (amyloid → cortical entry).** Higher amyloid load goes with being at Braak stage I or above rather than 0, with the same adjustment. *Against:* AGAINST or practically null: amyloid would not act at the step nearest the LC → cortex transition. If H1 supports and H2 does not, amyloid acts on spread within the cortex, not on entry; I will report that as evidence against amyloid gating of the LC → cortex step. This is the one confirmatory test here that discriminates.

**H3 (social and sensory exposures → Braak stage),** one test per exposure: smaller social network (`soc_net`), more loneliness (`social_isolation`), less late-life social activity (`late_life_soc_act`), and worse near vision (`vision`), each going with a higher Braak stage, with the same adjustment. *Against:* AGAINST for that exposure (the association is on the protective side). That would count against the premise that this factor raises tau pathology, which is the path Killer 2 needs. *Known property:* with the sample sizes I expect, the practically-null verdict is out of reach (section 6), so a factor with no effect will read "inconclusive", not "against". An inconclusive H3 is no evidence either way. A "supports" H3 does not show that the LC is involved: in simulation, a world where the whole effect runs through the LC and a world where none of it does gave the same H3 result (section 6).

**H4 (only if LC measures are delivered; Killer 2).** For each of the four exposures, the natural indirect effect (NIE) through the LC measure on expected Braak stage is positive: risk factor → more LC damage → higher Braak stage. Family C: 4 tests, Bonferroni-adjusted bootstrap intervals.

* *Killer 2 is met for an exposure* if the adjusted NIE interval includes 0 or lies below 0.
* *Killer 2 is met overall* if it is met for every exposure that can be estimated. It is not met if at least one exposure shows a predicted-direction NIE.
* A predicted-direction NIE whose ρ\* (the error correlation at which it vanishes, section 5.5) is below 0.1 in absolute value is reported as **fragile**. Fragile is not the same as support.
* If Killer 2 is met, the hypothesis is wrong in its current form, and the report says so in its first paragraph.

**Overall reading.** The result counts against LC-first if H2 is AGAINST or practically null, if H1 is, or (with LC measures) if Killer 2 is met. Without LC measures nothing in this module can count in favour of mediation through the LC.

## 3. Data

**Source.** RADC Standard dataset for ROS and MAP, requested by me from the RADC Research Resource Sharing Hub; not yet received. Both studies enrol older people without known dementia, follow them annually and include brain donation (https://www.radc.rush.edu/docs/parentStudyDesigns.htm). MARS, Clinical Core and Latino Core participants are excluded if delivered.

**Analysis sample.** All ROS and MAP participants with a Braak stage (`braaksc`), that is, autopsied participants. Each test uses the complete cases for its own variables (section 5.4). For `social_isolation` and `late_life_soc_act`, which ROS does not collect, the sample is MAP only.

**Sample size.** I do not know it. It is fixed by the cohort, not chosen by me. I did not consult RADC's Data Availability Report. Section 6 gives the operating characteristics at 500, 1000 and 1500 autopsied participants.

**Variables.** The Standard dataset page (https://www.radc.rush.edu/docs/var/standardDatasetVariables.htm, read 2026-10-21) says it offers 84 variables and lists 83. Names and codings below come from that page and from each variable's public page at https://www.radc.rush.edu/docs/var/detail.htm?variable=NAME (read 2026-10-21). Every name will be checked against the delivered codebook, names only, before any analysis (`data.codebook_confirmed`).

| Role | Variable | Definition on the public page | Status |
|---|---|---|---|
| Outcome | `braaksc` | Braak NFT stage 0–VI, Bielschowsky silver stain | confirmed (public page) |
| Outcome, sensitivity | `b_score` | Braak stage in 4 levels: 0, I–II, III–IV, V–VI | confirmed (public page) |
| Exposure H1/H2 | `amylsqrt_est_8reg` | amyloid density, mean over 8 regions of square-root % area | confirmed (public page) |
| Exposure, sensitivity | `a_score_st4` | Thal-based A score, 4 levels | confirmed (public page) |
| Exposure H3 | `soc_net` | social network size (number seen at least monthly); ROS and MAP, different questions in ROS; longitudinal | confirmed (public page) |
| Exposure H3 | `social_isolation` | perceived social isolation (loneliness), 5 items, 1–5; MAP only; longitudinal | confirmed (public page) |
| Exposure H3 | `late_life_soc_act` | late-life social activity, 6 items, 1–5; MAP only; longitudinal | confirmed (public page) |
| Exposure H3 | `vision` | near-vision acuity, 1 (20/40) to 7 (worse than 20/400); longitudinal | confirmed (public page) |
| Covariate | `age_death` | age at death | confirmed (public page) |
| Covariate | `msex` | sex, 1 = male | confirmed (public page) |
| Covariate | `educ` | years of education | confirmed (public page) |
| Covariate | `apoe_genotype` | APOE genotype 22–44; codes 24, 34, 44 = ε4 carrier | confirmed (public page) |
| Timing / S9 | `age_bl`, `died`, `dcfdx`, `cogn_global` | baseline age; death; diagnosis at each visit (1 = no cognitive impairment); global cognition | confirmed (public pages) |
| Keys | `projid`, `study`, `fu_year` | participant id, cohort, years since baseline | to be confirmed against the codebook on receipt (not on the public list) |
| Age at each visit | none listed | computed as `age_bl` + `fu_year` unless the delivery has a visit-age column | to be confirmed against the codebook on receipt |
| LC measure | none listed | the public variable index, the pathology and neuroimaging overviews and the substudy page list no LC variable | to be confirmed with RADC |
| Hearing | none listed | no hearing variable in the public variable index | not available |

**Exposure timing.** The four H3 exposures are measured at yearly visits. The primary value is each participant's earliest non-missing value, which in practice is the baseline value, the one furthest from death and so least affected by the years of decline before it.

## 4. Causal assumptions (fixed now)

```
X (common causes): age at death, sex, years of education, APOE ε4 carrier, cohort
A: amyloid load        T: social / sensory exposure        M: LC measure (often unobserved)
Y: Braak stage         U: unmeasured (prodromal decline, lifelong traits)    S: autopsy

X → A, X → T, X → M, X → Y
A → Y                         (H1, H2)
T → M → Y,  T → Y             (H3 is the total; H4 the part through M)
A → M (allowed)
U → T,  U → Y                 (reverse causation; addressed by S5, S6, E-values)
X, T, cognition → S           (selection to autopsy; addressed by S9)
```

Amyloid is not in the primary H3 model. It may lie on the path from exposure to tau, or reflect disease that also changes social life; S7 adds it.

The associations can be read causally only if, given X, there is no unmeasured confounding of the exposure and Braak stage, and selection into autopsy does not depend on both. Amyloid and Braak stage are measured at the same autopsy, so their order is assumed from biology, not observed. For H4 the usual sequential ignorability assumptions are added: no unmeasured exposure–LC and LC–Braak confounding given X and the exposure, and no LC–Braak confounder affected by the exposure. The last is the most likely to fail: general neurodegeneration damages the LC and raises tau. Mediation analysis makes strong causal assumptions from observational data, and the write-up will say so in its abstract.

## 5. Methods

### 5.1 Models

* **H1, H3:** proportional-odds (cumulative logit) regression of `braaksc` (7 levels) on the exposure and the adjustment set age at death, sex, education, APOE ε4 carrier and cohort. The cohort indicator is dropped, and the drop reported, where it is constant (MAP-only exposures).
* **H2:** logistic regression of Braak stage I or above versus 0 on amyloid and the same adjustment set.
* **Scale:** each exposure is standardised by its SD in the analysis sample, and oriented so that the predicted risk side is positive: smaller `soc_net`, larger `social_isolation`, smaller `late_life_soc_act`, larger `vision` code.
* **Intervals:** Wald intervals from the observed information, maximum likelihood with exact gradients (cross-checked against an independent implementation).

### 5.2 Secondary estimates (two-sided, no verdict)

* D1: the age at death association with Braak stage, given amyloid (how much age adds once amyloid is known; this informs Module 2's age slope).
* D2: the amyloid × age interaction.
* D3: the age association in amyloid-free brains (A score 0), if at least 30 brains are available.
* D4 (LC branch only): the LC × amyloid interaction, a direct look at gating.

### 5.3 Unmeasured confounding of the associations

For every confirmatory odds ratio I report the E-value (VanderWeele and Ding 2017) for the estimate and for the adjusted interval limit nearer the null, using the square-root conversion of the odds ratio because higher Braak stages are common.

### 5.4 Missing data

* **Primary:** complete cases per analysis. Every excluded participant is counted and logged with the variables that were missing.
* **Structural absence** (ROS has no `social_isolation` or `late_life_soc_act`) is not missing data. Those analyses are MAP only.
* **S8:** if more than 5% of the eligible participants are excluded for other reasons, the analysis is repeated with multiple imputation by chained equations (20 imputations, 10 cycles, predictive mean matching, over the outcome, the exposure and the adjustment set), pooled by Rubin's rules.

### 5.5 The conditional LC branch (H4)

**What decides the branch, and when.** The branch is decided by what RADC delivers, never by results. Before opening the data, I read the delivered codebook (names and descriptions only) and RADC's written answer about LC measures. If an LC measure is there, I enter it in `config/mediation.yaml` (`lc.variables`), commit the change, and only then open the data. If several are delivered I use the first available of: LC neuronal density, LC neuronal loss rating, LC tangle density. Each is oriented so that higher means more damage.

**If LC measures exist:**

* The LC measure is z-scored. Its model is linear, given the exposure and the adjustment set.
* Braak stage is modelled by ordered probit, given the exposure, the LC measure and the adjustment set, with no exposure × LC interaction.
* Natural indirect and direct effects on the expected Braak stage are computed exactly from the two models (the latent outcome is normal given the covariates), for an exposure contrast from its 25th to its 75th percentile, moved in the risk direction.
* Intervals: 2000 case-resampling bootstrap replicates, percentile intervals, Bonferroni over the 4 exposures.
* **Sensitivity to unmeasured LC–Braak confounding:** ρ, the correlation of the two models' errors, from −0.95 to 0.95 in steps of 0.01. For each ρ the outcome model is corrected exactly (the ordered-probit analogue of Imai, Keele and Yamamoto 2010). I report ρ\*, at which the NIE is 0, and the smallest |ρ| at which its 95% interval includes 0.
* S3 for H4: the same decomposition with Braak stage treated as a number 0–6 and linear models, using the Imai–Keele–Yamamoto ρ in closed form.
* The branch needs at least 100 autopsied participants with the LC measure and the exposure. With fewer, results are descriptive, with no verdict.

**If not:** H4 is not run. The report states that this dataset cannot test mediation through the LC, and Killer 2 stays untested.

### 5.6 Multiple testing

Bonferroni within each family, as in section 2: Family A (H1, H2), Family B (H3 × 4), Family C (H4 × 4). The families answer different questions and are not pooled. Secondary estimates and sensitivity analyses are not adjusted and carry no verdict of their own.

## 6. Operating characteristics, measured on simulated data before any real data

`sim/mediation_power.py` generates synthetic autopsy cohorts shaped like ROS/MAP (`sim/mediation_sim.py`). Every effect size and distribution in them is an assumption, not a RADC value. Results: `manuscript/notes/module3_power.md`, 200 data sets per cell for H1–H3, 100 per cell for the LC branch. Simulated sample sizes: 500, 1000 and 1500 autopsied participants, half MAP. The share of brains at Braak stage 0 is 3% or 10%.

* **No effects at all:** each test came out "supports" in at most 2% of data sets, and AGAINST in at most 4%.
* **H1:** "supports" in 200 of 200 data sets in every world where amyloid speeds cortical spread (odds ratio 2.0 per SD), including the world where amyloid does nothing at cortical entry. H1 does not discriminate.
* **H2** discriminates.
  * Where amyloid has no effect on entry, it said "supports" in at most 2%.
  * Where amyloid raises the odds of entry 1.5-fold per SD, it said "supports" in 23% of data sets at 500 brains with 3% at stage 0, rising to 97% at 1500 brains with 10% at stage 0.
  * Its power depends on how many brains are at stage 0, which I do not know.
* **H3:** with a true odds ratio of 1.15 per SD, "supports" in 7–72% of data sets, depending on the exposure and sample size; with 1.3 per SD, 24–99%. The MAP-only exposures have about half the sample and the least power.
* **The practically-null verdict was never reached**, in any cell, including the null world (0 of 200 in every cell). So in practice H1–H3 can end only as supports, AGAINST (opposite side), inconclusive or uninformative.
* **Identifiability.** Two simulated worlds had the same total effect of the exposures on Braak stage: in one the whole effect runs through the LC, in the other none of it does.
  * H3 gave the same result in both: mean log odds ratio 0.213 and "supports" in 93%.
  * Only the NIE told them apart: predicted-direction NIE in 100% versus 1%.
  * Without LC measures, the Standard dataset cannot say whether any association runs through the LC.
* **LC branch (exposure `vision`):**
  * With no exposure → LC path, a predicted-direction NIE appeared in 1%, and an opposite one in at most 2%.
  * With a modest path, a predicted-direction NIE in 46% at 300 participants and 99% at 800; with a strong path, in 95% and 100%.
  * The median ρ\* was 0.18–0.27, so even true simulated effects would be called robust only against weak confounding.

## 7. Validation that must pass first (`make validate-mediation`)

1. Exact identities of the linear decomposition and its ρ formula.
2. Recovery and 95% bootstrap coverage under sequential ignorability.
3. The ρ correction recovering the truth when unmeasured confounding is simulated.
4. A null mediator triggering Killer 2.
5. **Reproduction of Imai, Keele and Tingley (2010), Table 1** (JOBS II, from the authors' CC0 replication data, n = 1,285).
   * ACME −0.0157 against the published −.016 [−.03, −.002].
   * Direct effect −0.0316 against −.032.
   * One discrepancy, reported rather than hidden: the effect on the mediator is 0.1008 in the authors' own data and code, which rounds to .101, while the table note prints .100. It is checked to 0.001.
6. **Reproduction of Imai, Keele and Yamamoto (2010), Section 6.2** (Nelson et al. framing data, n = 136).
   * ACME −0.510 and total effect −0.540, as in Table 3.
   * The ACME changes sign at ρ = −0.675, against the published "unless ρ is less than −0.68" on a 0.01 grid.
7. The proportional-odds fit recovering a known log odds ratio with nominal coverage.
8. The cortical-entry test seeing amyloid only when amyloid acts on entry.
9. The ordinal mediation formulas, with and without ρ, recovering natural effects computed independently by Monte Carlo.
10. Bootstrap coverage of the ordinal NIE.

The published values were extracted by Claude and are marked unverified until I check them (`data/extracted/mediation_validation_published.csv`).

## 8. Sensitivity analyses (all pre-specified)

* S1. Thal-based A score (`a_score_st4`) instead of amyloid density (H1, H2).
* S2. Braak stage in 4 levels (`b_score`) instead of 7 (H1, H3).
* S3. Braak stage as a number 0–6, linear regression with robust errors (H1, H3); linear mediation for H4.
* S4. Exposure as the mean of all visits, and as the last value before death (H3).
* S5. Only exposures measured at a visit with no cognitive impairment (`dcfdx` = 1) (H3).
* S6. Only exposures measured at least 5 years before death (H3).
* S7. H3 additionally adjusted for amyloid: the association not running through amyloid.
* S8. Multiple imputation, if triggered (section 5.4).
* S9. Inverse-probability-of-autopsy weights among deceased participants, predicted from baseline age, sex, education, cohort, baseline global cognition and, for H3, the exposure itself. Weights stabilised, truncated at the 99th percentile, 500 bootstrap replicates with the weights re-estimated in each.
* S10. ROS and MAP separately, where the exposure exists in both.

**Robustness.** A "supports" or AGAINST verdict is called robust only if S1, S2, S3, S5 and S6, where they apply, put the adjusted interval on the same side of the null. A practically-null verdict is robust only if none of them excludes the null. The others are reported for transparency.

## 9. Limitations, stated in advance

* No LC or pre-tangle stages: the LC → cortex step itself is not observed (stated at the top). H2's cortical-entry split is the nearest these data come, and stage 0 mixes brains with and without LC tau.
* Silver-stain Braak staging misses pre-tangle tau.
* Amyloid and Braak stage are measured at the same autopsy. Their order is assumed, not observed.
* Only autopsied participants have a Braak stage. Selection is addressed by S9 but not removed.
* Social activity, loneliness and network size can fall because of disease that has already begun; S5 and S6 address this only partly.
* `social_isolation` measures loneliness, not objective isolation. `vision` is corrected near acuity, which reflects eye disease as much as sensory input. ROS asks different network questions from MAP.
* The MAP-only exposures have about half the sample.
* No hearing measure, the exposure my brief named first, is available.
* Every variable name is to be confirmed against the delivered codebook. The keys and the visit age are not on the public list.
* Published ROS/MAP papers have studied some of these variables, and I know of that literature in general terms (see section 10).

## 10. What has been seen

Nothing from the RADC data: no values, no summaries, no counts. Claude read only RADC's public documentation pages (variable names, definitions, codings, and which cohorts collected each variable). The hypotheses and their directions come from the LC-first framework in my project brief, which was written before this dataset was chosen. Results from published ROS/MAP analyses are part of the public literature, which both I and the AI model know in general terms. That is why H1 is labelled non-discriminating, and why every verdict rule above was fixed and simulated before the data arrive.
