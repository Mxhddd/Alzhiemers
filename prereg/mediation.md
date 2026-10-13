# Pre-registration — Module 3: does LC integrity mediate risk factors → tau? (Killer 2)

**Status: COMPLETE DRAFT with one open item that only the researcher can
close: the cohort.** Every methodological decision is made below. The
cohort (name, release, variable coding) cannot be chosen by Claude. It
needs data access, and the brief forbids substituting a public dataset
without asking. Freeze this document only once §4 names the cohort.
Claude must not run `make freeze-prereg`. The DAG is fixed here before
anything is fitted, as the brief requires. Configuration:
`config/mediation.yaml`.

> **DO NOT FREEZE until §4 names the cohort** (data access is the researcher's).

## Change log

| Date | Change |
|---|---|
| 2026-10-13 | Complete draft. Replaces `mediation.DRAFT.md`. Decisions: primary exposure, the operational definitions, amyloid's role in the DAG, the missing-data rule, the Killer 2 criterion and its robustness reading. |

## 1. Killer 2 (stated first)

LC integrity does **not** mediate the association between sensory/social
risk factors and tau pathology.

**Killer 2 is met** for the primary exposure (hearing loss) if the average
causal mediation effect (ACME) 95% bootstrap CI includes 0 or lies on the
opposite side to the prediction. The predicted ACME is positive: risk
factor → lower LC integrity → more tau. If Killer 2 is met, the hypothesis
is wrong in its current form, and the report says so in its first
paragraph.

**Robustness reading.** The ρ sensitivity analysis gives the correlation
of mediator and outcome errors at which the ACME becomes 0 (ρ\*). A
positive ACME with |ρ\*| < 0.1 is reported as **fragile**: a small amount
of unmeasured mediator–outcome confounding would erase it. Fragile is not
the same as "supports".

## 2. DAG (fixed now)

```
X (common causes): age, sex, years of education, APOE ε4 carrier status,
                   hypertension, diabetes, amyloid status (see below)
T (exposure):      hearing loss [primary]; social isolation, physical inactivity [secondary]
M (mediator):      LC integrity
Y (outcome):       early tau burden

X → T,  X → M,  X → Y
T → M → Y     (the mediated path under test)
T → Y         (direct path)
```

**Amyloid is a common cause of M and Y** (X → M, X → Y) in the primary
analysis. Amyloid pathology is associated with lower LC integrity, and
it accelerates tau spread. Two sensitivity analyses:

- (a) not adjusting for amyloid;
- (b) restricted to amyloid-positive participants, the population of
  Killer 1.

If hearing loss itself raises amyloid, then amyloid is on a T → Y path and
adjusting for it violates assumption 4 below. That is why (a) is
reported, and any difference between (a) and the primary analysis is
discussed.

## 3. Assumptions required for identification

1. **No unmeasured exposure–outcome confounding** given X.
2. **No unmeasured exposure–mediator confounding** given X.
3. **No unmeasured mediator–outcome confounding** given X and T.
   - This is the assumption most likely to fail: general
     neurodegeneration lowers LC contrast and raises tau.
   - The ρ analysis quantifies how much violation the result survives.
4. **No mediator–outcome confounder affected by the exposure.** Amyloid
   and depression are the candidates; see §2.
5. **Linear models with no T × M interaction** (primary). The interaction
   model is a secondary analysis.
6. **Positivity and consistency.**
7. **Measurement:** neuromelanin-sensitive LC contrast is a valid proxy
   of LC regulatory integrity. This is strong and partly untestable.
8. **Temporal order:** T is measured before M, and M before Y.
   - With cross-sectional data, the result is reported as a
     decomposition of an association, not a causal mediation estimate.
   - Killer 2 can still be met (no mediated association), but support
     cannot be claimed.

Mediation analysis makes strong causal assumptions from observational
data. The write-up says so in its abstract.

## 4. Data (cohort to be named by the researcher)

**Required:** in the same participants, measured in this order:

1. the exposures;
2. neuromelanin-sensitive LC MRI;
3. tau-PET.

**Candidates:** cohorts with LC MRI and tau-PET exist, e.g. the Harvard
Aging Brain Study. Choosing and obtaining one is the researcher's
decision and access. **Cohort:** [to be named before freezing].

Operational definitions, fixed now and adapted only in units to the
cohort's instruments (any adaptation is logged):

| Variable | Definition |
|---|---|
| Hearing loss (primary T) | Audiometric pure-tone average ≥ 25 dB HL in the better ear, if measured; otherwise the cohort's self-report item. The source is stated. |
| Social isolation | The cohort's standard isolation scale, dichotomised at its published cut-off. |
| Physical inactivity | Below 150 minutes of moderate activity per week. |
| LC integrity (M) | Rostral–middle LC contrast ratio from neuromelanin-sensitive MRI (higher = more intact), z-scored. |
| Early tau (Y) | Tau-PET SUVR in an entorhinal / Braak I–II meta-ROI (inferior cerebellar grey reference), z-scored. |

## 5. Analysis

- **Outcome and mediator models:** linear regressions of M on T and X,
  and of Y on T, M and X.
- **Effects estimated:** ACME, ADE, total effect and proportion mediated,
  with 5,000 nonparametric bootstrap replicates and percentile 95% CIs.
- **ρ sensitivity:** ρ from −0.95 to 0.95 in steps of 0.01. Report ρ\*
  and the smallest |ρ| at which the ACME CI includes 0.
- **Missing data:** complete-case analysis, with every excluded
  participant counted and logged by reason. Sensitivity: multiple
  imputation (20 imputations, chained equations including T, M, Y and X)
  if more than 5% of participants are excluded.
- **Secondary:**
  - the T × M interaction model;
  - the secondary exposures;
  - amyloid sensitivity analyses (a) and (b).
- **Multiplicity:** Killer 2 rests on the primary exposure. Secondary
  exposures are reported with Bonferroni-adjusted (3-test) CIs alongside.

## 6. Validation that must pass first (`make validate-mediation`)

- Simulation recovery of known ACME, ADE and ρ behaviour.
- Reproduction of a published mediation estimate. Its values are to be
  extracted from the primary paper with locators; that is not done in
  this round.
- See `tests/validation/test_mediation.py`.

## 7. Limitations, stated in advance

- Observational mediation rests on assumptions 1–4.
- LC MRI is a proxy.
- Cross-sectional designs cannot order T, M and Y.
- The cohort is not yet available.
