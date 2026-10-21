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
| braak2011 | Braak H, Thal DR, Ghebremedhin E, Del Tredici K. Stages of the pathologic process in Alzheimer disease: age categories from 1 to 100 years. *J Neuropathol Exp Neurol* 2011;70(11):960–969. doi:10.1097/NEN.0b013e318232a379. PMID 22002422 | **Full text PDF supplied by the researcher on 2026-10-13**, saved as `data/raw/papers/braak2011.pdf` (sha256 9bc53256…f2bd2). **Supplemental Digital Content 1 and 2** were downloaded on 2026-10-13 from the publisher's open permalinks https://links.lww.com/NEN/A273 and https://links.lww.com/NEN/A274 (HTTP 200, no login, no bot check). Saved as `braak2011_sdc1.doc` (sha256 fdccf6bc…2b256c) and `braak2011_sdc2.doc` (sha256 44b988b7…7466). The earlier abstract-only extraction is kept in `braak2011_abstract_totals.csv` |

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
- **braak2011** (full extraction 2026-10-13, `python -m src.extraction.braak2011`):
  - `braak2011_table1.csv`: Table 1, sex by age category (53 cells).
  - `braak2011_table2.csv`: Table 2, amyloid phase 0–4 by AT8 stage group
    (0, a–c, 1a–1b, I … VI), n and row %, plus the column headers (120 cells).
    There is no age dimension.
  - `braak2011_text_counts.csv`: every count printed in the Results and
    Discussion text, the Figure 1 panel titles and the totals above the
    Figure 2 bars (51 rows). Bar heights were not read.
  - `braak2011_stage_counts.csv`: **Supplementary Table 1**, the 12 AT8 stages
    (0, a, b, c, 1a, 1b, I–VI) by 10 age categories (120 counts). This is
    Module 2's input. SDC 1 prints the first age category as "00-10"; the
    paper's Table 1 prints "1–10" and the cohort is aged 1–100, so the band
    is stored as 1–10 (`age_label_printed` keeps "00-10").
  - `braak2011_sdc1_groups.csv`: SDC 1 grouped columns and printed row
    percentages, plus its total row (144 cells).
  - `braak2011_sdc2_amyloid_by_age.csv`: **Supplementary Table 2**, amyloid
    phase 0–4 by age category, n and row % (132 cells). Module 2's A(t) is
    derived from it in code.
  - `braak2011_sum_checks.md`: 51 consistency checks. **All 43 count checks
    pass**: rows and columns sum to the printed totals, and the totals agree
    across Table 1, Table 2, both supplements and the text. Two checks of
    printed percentages or prose numbers show the **source disagreeing with
    itself**:
    - SDC 1 prints 34.3% for 21/61 (34.4%), age 21–30, stages a–c.
    - The Discussion says "52 cases … (pretangle Stage 'a')". The Results
      and SDC 1 give 40 stage-a cases; 52 = stage a + stage b. The 52 is
      stored as printed and not used.
  - **Per-brain amyloid by tau stage by age is not in the paper or its
    supplements.** Table 2 crosses amyloid with tau stage without age. SDC 1
    crosses tau stage with age without amyloid. SDC 2 crosses amyloid with
    age without tau stage. Figure 3 plots tau stage against age and amyloid
    against age in separate panels, which cannot be linked per brain. Three
    two-way tables do not determine the three-way table.
  - **Blind second reading (2026-10-14):** a separate agent read rendered page
    images of the paper and both supplements, never the stored values or the
    text layer, and transcribed every table, text count and figure label.
    `python -m src.extraction.recheck_braak2011` compared all of them:
    **620 of 620 values agree**, so `machine_rechecked` is true on every row
    (log: `recheck/braak2011/recheck_log.md`). `verified` stays false.
    The second reader independently found the same two inconsistencies (34.3%
    and the Discussion's "52"). It also noted three wording points that change
    no number:
    - the Results say 1b cases are "designated in Table 2", but Table 2 has
      only a combined 1a–1b row (the 236 is in SDC 1);
    - the Discussion's 1,291 "pretangle or early NFT stages without plaques"
      is phase 0 minus the 10 no-change cases, which also includes phase-0
      brains at stages III–V;
    - the Results cite "Fig. 1A", but Figure 1 has no panel letters.

## Not found / blocked

- ~~Braak 2011 stage counts by age~~ and ~~amyloid by age~~: resolved
  2026-10-13 from the paper's open supplements (above).
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

## Module 3 validation: published mediation analyses (added 2026-10-21)

File: `data/extracted/mediation_validation_published.csv` (all rows `verified: false`).
Extracted by Claude from the open full texts; used only to validate the
estimator in `tests/validation/test_mediation.py`, never as Module 3 data.

- `imai_keele_tingley2010`: Imai K, Keele L, Tingley D. A general approach to
  causal mediation analysis. *Psychol Methods* 2010;15(4):309–334.
  doi:10.1037/a0020761. Read from the authors' copy
  (https://imai.fas.harvard.edu/research/files/BaronKenny.pdf), Table 1 and its note.
  Data and code: the authors' replication archive, Harvard Dataverse
  doi:10.7910/DVN/UMEYXD (CC0 1.0), file `Jobs-NoMiss-Cont.dta` (n = 1,285) and
  `Table-1.R` (model: `depress2 ~ treat + job_seek + depress1 + age + educ + sex +
  nonwhite + econ_hard + marital + occp + income`, factors as stored).
- `imai_keele_yamamoto2010`: Imai K, Keele L, Yamamoto T. Identification,
  inference and sensitivity analysis for causal mediation effects. *Stat Sci*
  2010;25(1):51–71. doi:10.1214/10-STS321; arXiv:1011.1079 (read there).
  Table 3 and Section 6.2. Data and code: Harvard Dataverse
  doi:10.7910/DVN/APNMP7 (CC0 1.0), `nelson_data.RData` (n = 136),
  `nelson-analysis.R`, `sem-sens.R` (models without covariates:
  `idisord ~ vidcdum`, `kspeech ~ vidcdum + idisord`).
- Download and conversion: `python -m src.mediation.fetch_validation_data`
  (sha256-pinned; writes `data/raw/validation/`, not committed).
- Note: the `jobs` data shipped in the R package `mediation` (n = 899) is a
  different extract and does NOT reproduce Table 1 (effect on the mediator
  0.077 against the published .100); only the Dataverse file is used.

## Module 6: locus coeruleus atlas (added 2026-10-21)

File: `data/extracted/network_lc_atlas.csv` (both rows `verified: false`).
Extracted by Claude through the NCBI E-utilities PMC XML (the PMC web page
showed a browser check and was not used). No atlas file was downloaded.

- `ye2021`: Ye R, Rua C, O'Callaghan C, et al. An in vivo probabilistic atlas
  of the human locus coeruleus at ultra-high field. *NeuroImage*
  2021;225:117487. doi:10.1016/j.neuroimage.2020.117487. PMC7779564. The two
  thresholds (5% liberal, 25% stringent) name the two distributed versions of
  the atlas; the "Data and code availability statement" says: "The thresholded
  atlases are freely available at https://www.nitrc.org/projects/lc_7t_prob/ ."
- `keren2009`: Keren NI, Lozar CT, Harris KC, Morgan PS, Eckert MA. In vivo
  mapping of the human locus coeruleus. *NeuroImage* 2009;47:1261–1267.
  doi:10.1016/j.neuroimage.2009.06.012. PMC3671394. Not used. Recorded because
  the 2026-10-13 Module 6 draft described this map as "thresholded at 25%";
  the paper (Discussion) distributes "two binary templates representing the
  1 SD and 2 SD extents of LC spatial distribution", so that description was
  wrong.
