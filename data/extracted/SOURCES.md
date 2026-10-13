# Extracted values — sources and verification

Every value in `data/extracted/*.csv` was extracted by Claude from the primary
paper named in its row and is marked `verified: false` until the researcher
checks it against the PDF and flips it to `true`. Analyses may run on
unverified values, but every output built on them carries an UNVERIFIED
banner. Nothing here is estimated, interpolated or carried over from another
study; a value that could not be found is left out and listed below.

These files are committed (unlike `data/raw/`) so that verification flips are
tracked in git.

## Sources

| id | Citation | Access used |
|---|---|---|
| matthews2016 | Matthews FE, Stephan BC, Robinson L, et al. A two decade dementia incidence comparison from the Cognitive Function and Ageing Studies I and II. *Nat Commun* 2016;7:11398. doi:10.1038/ncomms11398. PMID 27092707, PMC4838896 | Europe PMC full-text XML (open access) |
| jansen2015 | Jansen WJ, Ossenkoppele R, Knol DL, et al. Prevalence of cerebral amyloid pathology in persons without dementia: a meta-analysis. *JAMA* 2015;313(19):1924–1938. doi:10.1001/jama.2015.4668. PMID 25988462, PMC4486209 | NCBI E-utilities PMC XML (author manuscript) |
| ossenkoppele2022 | Ossenkoppele R, Pichet Binette A, Groot C, et al. Amyloid and tau PET-positive cognitively unimpaired individuals are at high risk for future cognitive decline. *Nat Med* 2022;28(11):2381–2387. doi:10.1038/s41591-022-02049-x. PMID 36357681, PMC9671808 | Europe PMC full-text XML (open access) |
| braak2011 | Braak H, Thal DR, Ghebremedhin E, Del Tredici K. Stages of the pathologic process in Alzheimer disease: age categories from 1 to 100 years. *J Neuropathol Exp Neurol* 2011;70(11):960–969. doi:10.1097/NEN.0b013e318232a379. PMID 22002422 | **Abstract only** (PubMed). Full text is paywalled and the publisher site blocks automated access |

## Notes for the verifier

- **matthews2016 Table 1**: rates are by *age at baseline interview*,
  all-cause dementia (GMS-AGECAT algorithm), modelled with weighted Poisson
  regression and multiple imputation (Methods). The 85+ band is open-ended
  (`age_upper = open`); where the model needs a single age for it, the
  config states the choice.
- **jansen2015 Table 2**: modelled (GEE) prevalence at exact ages; the PMC
  copy is an author manuscript with visible transcription errors elsewhere in
  the tables (Table 3 "128/53C"; Table 2 MCI APOE-ε4+ age 85 CI printed
  "86.3 (32.9-89.2)"). None of those cells were extracted; the extracted
  column (normal cognition, total) looked clean but please check it.
- **ossenkoppele2022 Table 2**: HRs vs A−T− reference. Note: in the PMC copy
  of **Table 1** the "Progression to MCI" percentages for A−T− and A+T−
  appear swapped relative to the Results text (text: 26/781 = 3.3% for
  A−T−, 26/292 = 8.9% for A+T−). Table 1 was not used.
- **braak2011**: only abstract totals are extracted (used as consistency
  checks: the full stage-by-age table must reproduce them).

## Not found / blocked

- **Braak 2011 stage counts by age** (`braak2011_stage_counts.csv` is empty):
  the per-decade × per-stage counts are in the paywalled full text.
  Needed: the PDF, placed at `data/raw/papers/braak2011.pdf` (git-ignored).
- **Braak 2011 amyloid-plaque prevalence by decade** (same cohort; primary
  A(t) for model A): same PDF.
- **An amyloid-stratified hazard ratio for incident tau-PET positivity**
  (closer to the LC→cortex transition than clinical progression): searched,
  not found in open access. The γ derivation uses ossenkoppele2022 as a proxy
  and says so.

## Independent second-pass check (2026-10-13)

A second reader re-read every extracted value. It saw only the source,
the locator and a description of the quantity, never the stored numbers.
Its readings are in `recheck/`, and `python -m src.common.apply_recheck`
compares them with the stored values.

- **Result:** 86 of 87 rows agree exactly.
- **The one discrepancy** is the description of
  `braak2011_abstract_totals.csv` row 8. The abstract says "4% of all
  cases", not "percent of cases in that decade". The number agrees, and the
  unit text now quotes the source.
- **The `machine_rechecked` column** records the outcome per row. `verified`
  is untouched; only a human sets it.

Oddities the second reader noticed in cells NOT extracted (left as printed):

- **Jansen 2015 Table 2, age 85:**
  - Normal-cognition APOE-ε4− is printed "34.2 (27.7 41.4)", with the dash missing.
  - MCI APOE-ε4+ is printed "86.3 (32.9-89.2)". A lower bound of 32.9 is implausible; it is probably 82.9.
  - The footnote prints "Mo estimate" for "No estimate".
- **Ossenkoppele 2022 Table 2:** all-cause dementia, A+T− at 80+, is printed "0.90 (0.50–15.98)". That CI does not fit the point estimate on the log scale. It is not used.
- **Braak 2011 abstract:** two numbers are written in words ("Two thousand three hundred and thirty two"; "Fifty-eight").
