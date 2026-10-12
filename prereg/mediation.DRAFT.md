# Pre-registration DRAFT — Module 3: does LC integrity mediate risk factors → tau? (Killer 2)

> **DRAFT by Claude for researcher review; not frozen.** The brief requires
> the DAG to be pre-specified before anything is fitted. Edit, save as
> `prereg/mediation.md`, upload to OSF, `make freeze-prereg NAME=mediation OSF=<url>`.

## Killer 2 (stated first)

LC integrity does **not** mediate the association between sensory/social risk
factors and tau pathology. Operationally, for each pre-registered exposure:
Killer 2 is met if the ACME 95% CI includes 0 or lies opposite to the
predicted direction. If met for [all / the primary] exposure(s), the
hypothesis is wrong in its current form and we say so.

## DAG (nodes and directed edges)

```
Exposures T: hearing loss, social isolation, physical inactivity   [researcher: final list, coding]
Mediator  M: LC integrity (neuromelanin-sensitive MRI contrast)     [researcher: measure, region]
Outcome   Y: tau burden (tau-PET SUVR, early Braak regions)          [researcher: measure]
Covariates X (common causes): age, sex, education, APOE4, vascular risk,
                              socioeconomic position, [amyloid status?] [researcher]

X → T,  X → M,  X → Y
T → M → Y          (the mediated path under test)
T → Y              (direct path)
```

[Researcher to decide: is amyloid a confounder (X → M, X → Y), an effect
modifier (stratify), or itself on the path (T → amyloid → Y)? Each choice
changes what the ACME means.]

## Assumptions required for identification (each stated plainly)

1. **No unmeasured exposure–outcome confounding** given X.
2. **No unmeasured exposure–mediator confounding** given X.
3. **No unmeasured mediator–outcome confounding** given X and T. This is the
   assumption most likely to fail (e.g. general neurodegeneration lowers LC
   contrast and raises tau). The ρ sensitivity analysis quantifies how much
   violation the result survives; it is reported for every exposure.
4. **No mediator–outcome confounder affected by the exposure** (e.g. if
   hearing loss → depression → both LC and tau, depression must be modelled).
5. **Linear models, no exposure × mediator interaction** (tested as a
   secondary analysis: [researcher]).
6. **Positivity and consistency**: every exposure level occurs at every
   covariate pattern; the exposure is well defined.
7. **Measurement**: LC MRI contrast is a valid proxy of LC regulatory
   integrity — a strong and partly untestable assumption.
8. **Temporal order**: exposure measured before mediator before outcome
   [cohort design permitting]; cross-sectional data cannot establish it.

Mediation analysis makes strong causal assumptions from observational data;
the write-up will say so in the abstract.

## Analysis
Bootstrap (5000, percentile 95% CIs) for ACME, ADE, total, proportion
mediated; ρ ∈ [−0.95, 0.95]; report ρ at which ACME = 0 and the smallest |ρ|
at which its CI includes 0. Missing data: [raise / complete-case, justified].
