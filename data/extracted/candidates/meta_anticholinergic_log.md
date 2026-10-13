# Anticholinergic exposure and incident dementia: extraction log

Run date: 2026-10-13. Output CSV: `meta_anticholinergic_rows.csv` (116 rows: 45 `primary_study`, 43 `ma_component`, 28 `ma_pooled`).
Saved source files (used for extraction, exactly as downloaded): `meta_anticholinergic/` (23 XML, 1 supplementary zip, plus the docx extracted unchanged from it). SHA-256 list: `scripts_anticholinergic/saved_files.sha256`.
Helper scripts: `scripts_anticholinergic/` (`epmc_search.py`, `xml2txt.py`, `docx2txt.py`, `screen_ma.py`, `build_rows.py` + `rows_spec.py`, `verify_csv.py`). The machine-readable query log is `scripts_anticholinergic/queries.log`, copied in full at the end of this file.

## 1. Sources, access and blocked hosts

- **Sources used:** the Europe PMC REST API (`search`, `{PMCID}/fullTextXML` and `{PMCID}/supplementaryFiles`) and NCBI E-utilities (`efetch db=pmc`). Requests went through curl or Python urllib. No API key and no email address were used anywhere. No other host was contacted.
- **Blocked hosts:** none. There was no Cloudflare, login, paywall or captcha.
- **Errors and limits:**
  - NCBI efetch returned **HTTP 429 (rate limit)** once, for PMC4358759. I paused and retried, and the retry succeeded.
  - Europe PMC `fullTextXML` returned **HTTP 500** for every non-open-access record I tried: PMC8169439 (Cochrane, Taylor-Rowan 2021), PMC4358759 (Gray 2015), PMC4541969 (Ruxton 2015), PMC8175463 (Lockery 2021), PMC9974875 (Gildengers 2023) and PMC2933398 (Carrière 2009).
  - NCBI efetch then gave the full text for two of these as **NIH author manuscripts**: PMC4358759 (Gray 2015) and PMC9974875 (Gildengers 2023). Europe PMC marks both `isOpenAccess=N`. The Gray and Gildengers rows are flagged this way in `data_source`.
  - For PMC8169439, PMC8175463 and PMC2933398, efetch returned only front matter with the comment "The publisher of this article does not allow downloading of the full text in XML form". These were not used.
  - Europe PMC `supplementaryFiles` returned **HTTP 500 twice** for PMC7583519 (Pieper 2020), so its appendix could not be retrieved.

## 2. Verbatim and locator conventions

- **Paragraph text:** the whitespace-normalised text content of the XML `<p>`. Every run of Unicode whitespace, including NBSP and thin space, becomes one ASCII space. All other characters are kept, including U+2010 hyphens, en and em dashes, § and ‡ markers, and ≥ and ≧.
- **Table rows:** the cell texts of one `<tr>`, whitespace-normalised and joined by ` | ` in row order. The `locator` names the column that the estimate comes from.
- **Multi-part verbatims** (used for pooled estimate plus method statements) join the snippets with ` [...] `.
- **Supplementary docx** (Dmochowski Supplemental Table 4/5): cell text is taken from `word/document.xml`, whitespace-normalised. Table 4 holds many estimates per cell, so each component row quotes the exact sub-string for that estimate.
- **Checks:** `build_rows.py` checks each verbatim against text regenerated from the *saved* files. `verify_csv.py` checks again against the raw XML or docx: every ` | ` piece must be present, and estimate, ci_lower and ci_upper must each appear inside the verbatim. Both report 0 problems.
- **pmcid/doi on `ma_component` rows:** these refer to the **meta-analysis**, which is where the number was copied from. `study_id` names the component study, and `parent_meta_analysis` names the meta-analysis.
- No number was computed, converted, rounded or inferred. The only parsing in the scripts splits the Table S5 CI string ("1.30 to 1.93") into its two printed endpoints.
- **Informal plausibility checks:** to rank the candidates I did a few rough mental re-pools (described in section 4). No computed number from them was recorded anywhere.

## 3. Goal 1: primary studies

### 3a. Included, with one main estimate each (`is_main_estimate` starting "yes")

| study_id | PMCID | access | main comparison (verbatim numbers in CSV) |
|---|---|---|---|
| Gray2015 | PMC4358759 | **NIH author manuscript (NCBI efetch); EPMC isOA=N** | >1095 TSDD vs 0, dementia, adjusted HR (Table 3) |
| Richardson2018 | PMC5915701 | OA | ACB-3 >1460 DDD vs 0, OR adjusted at end of DEP (Table 2) |
| Coupland2019 | PMC6593623 | OA | >1095 TSDD vs nonuse, 1–11 y window, adjusted OR (Table 4) |
| Grossi2019 | PMC6802337 | OA | recurrent ACB3 use vs none, adjusted IRR (Table 2). The paper's main models are "ever-use"; the ever-use estimate is recorded as non-main |
| Suh2020 | PMC7391507 | OA | KABS High (≥2) vs Minimal, 2–10 y, aOR (Table 4) |
| Joung2019 | PMC6494800 | OA | strong anticholinergics ≥120 vs 0–9 doses/yr, **AD outcome**, Model II HR (Table 3) |
| Yang2017 | PMC5383287 | OA | oxybutynin vs matched controls, after-match HR (Table 3). **Provisional choice:** only drug-specific HRs exist, and oxybutynin is the drug the paper's conclusion singles out |
| Liu2020 | PMC7033580 | OA | with vs without anticholinergics, aHR (Table 2). This is the paper's primary binary exposure; graded models are recorded as non-main. **CIs in this paper are markedly asymmetric (copied as printed)** |
| Poonawalla2023 | PMC10311860 | OA | ≥4 ACH meds vs no ACH exposure (time-varying), adjusted HR (Table 2). Outcome is dementia/AD combined |
| Zhu2025 | PMC12538848 | OA | strong anticholinergics ≥1095 DDD vs 0, 1-y lag, Model 3 AOR (Table 2) |
| Iyen2024 | PMC11580265 | OA | any OAB anticholinergic use vs none, adjusted OR (Table 4). This is the only class-wide estimate; per-drug >1095 TSDD rows are recorded as non-main |
| Chi2025 | PMC12695687 | OA | TSDD>0 vs TSDD=0, Model 2 HR (Table 3). **Discrepancy:** the abstract calls "HR 1.55, 1.08–2.24" adjusted, but in Table 3 1.5 (1.08, 2.24) is Model 1 (univariate) |
| PhamNguyen2025 | PMC12022375 | OA | time-varying ACB 3+ vs 0, Parkinson disease dementia (special population) |
| Amrein2026 | PMC13570850 | OA | long-term antidepressant users with ≥1 TCA vs PS-matched non-users, OR (Results text). Exposure is duration-defined |
| Li2025 | PMC12274944 | OA | OAB anticholinergic vs **beta-3 agonist (active comparator)**, age 18–64, Fine-Gray HR (Table 3). Flagged as not "highest vs none" |
| Gildengers2023 | PMC9974875 | **NIH author manuscript (NCBI efetch); EPMC isOA=N** | total ACR score >3 vs 0, dementia (CDR≥1), HR (Table 4) |

### 3b. Included with no main estimate (all rows `is_main_estimate = no`)

- **Asiimwe2025** (PMC12596250): only a continuous per-drug HR (ACMI count), for UK Biobank and for All of Us. There is no categorical comparison.
- **Mur2022** (PMC9005668): the main text gives per-scale HRs only as a range ("1.027 to 1.125"). The exact values are in supplementary Tables S7–S9, which were not retrieved or extracted. I recorded only the trajectory-slope HR.
- **Park2017** (PMC5215897): the anticholinergic OR is univariate (crude) under a matched design. It is recorded because it is a Dmochowski 2021 component, so the two sources can be cross-checked.

### 3c. Considered and rejected for Goal 1

- **Chuang 2017**, PMC5651434 (OA): the outcome is "mild cognitive impairment or Alzheimer's disease" (composite, Table 2), not dementia or AD alone. Exposure groups are possible vs definite users, not graded. Not extracted.
- **Moga 2017**, PMC5408467 (OA): the outcome is cognitive decline, not incident dementia.
- **Welk & McClure 2022**, PMC9732452: the outcome is cognitive change. Not downloaded.
- **Lockery 2021** (JGIM), PMC8175463: no full text through either allowed source (EPMC 500; NCBI says publisher does not allow).
- **Carrière 2009** (Arch Intern Med), PMC2933398: same as Lockery 2021.
- **Not in Europe PMC full text (isOA N, inEPMC N), so not retrievable:**
  - Heath 2018 (JAGS)
  - Hong 2019 (Parkinsonism Relat Disord)
  - Whalley 2012
  - Jessen 2010
  - Papenberg 2017
  - Welk 2020 (BJU Int)
  - Sheu 2019
  - Cai 2013
  - Campbell 2010 and 2018 (PMC5362375 is the 2016 Campbell paper on cognitive impairment: author manuscript, outcome not dementia)
- **Mur 2025 simulation study** (PMC12643546) and **Richardson 2024 conference abstract** (PMC11712707): skipped as not primary estimates, or abstract only.
- **Gray 2015 open-access status:** it is NOT open access in Europe PMC. Its full text is reachable only as an NIH author manuscript through NCBI efetch. The extraction was done but is flagged.

## 4. Goal 2: meta-analysis candidates

Searches for anticholinergic–dementia meta-analyses with open-access full text found only **Pieper 2020** (Age Ageing, PMC7583519) and **Dmochowski 2021** (Neurourol Urodyn, PMC7821204).

- **Zheng 2021** (Neurosci Biobehav Rev): not in Europe PMC full text.
- **Kotochinsky 2025** (Neurol Sci): not in Europe PMC full text.
- **Cochrane Taylor-Rowan 2021** (CD013540, PMC8169439): full text not downloadable as XML.
- **Ruxton 2015** (PMC4541969): EPMC 500. It is also about cognitive impairment, not dementia.

Neither anticholinergic meta-analysis tabulates exactly which component estimate entered its main pooled estimate, so I also screened other open-access dementia–drug meta-analyses (`screen_ma.py`).

### Ranked candidates (rows in CSV)

**1. Moraros 2017, Depress Anxiety, PMC5347943, doi 10.1002/da.22584.** Antidepressant use and dementia or cognitive impairment. **Most complete tabulation.**
- **Pooled estimate:** OR = 2.17 (95% CI: 1.41–3.33; I2 = 84.9%).
- **Method as stated:** "a random effects model defined by the equation Y i = μ + u i + ε i" and "STATA version 13.1 using the “Metan” command". The estimator (for example DerSimonian–Laird) is not named, and Knapp–Hartung is not mentioned.
- **Components:** all 5 crude ORs with 95% CI **and study weights** are in Table 3.
- **Extra checkpoints in the text:** 10 subgroup pooled ORs (with named component studies) and 2 leave-one-out pooled ORs, all entered as `ma_pooled` rows.
- **Caveats:**
  - The component ORs are **crude**.
  - The outcome mixes dementia, MCI and other cognitive disorders.
  - One component is cross-sectional (Sachdev).
  - Exposure is antidepressants, not anticholinergics specifically.
- **Internal inconsistencies in the source:**
  - Cherbuin's upper CI is "7.06" in Table 3 but "1.62–7.07" in the text.
  - The "no baseline depression" subgroup lists Goveas, who is also in the depression subgroup (it probably should be Kessing). The composition is flagged in the row.
  - The text calls Goveas 2012 an RCT, while Table 3 says "Prospective cohort".
- **Plausibility check:** an informal mental re-pool using the tabulated weights was consistent with the reported pooled OR.

**2. Dmochowski 2021, Neurourol Urodyn, PMC7821204, doi 10.1002/nau.24536.** Anticholinergics and incident dementia. **On-topic but only partially tabulated.**
- **Pooled estimate:** RR 1.46 (95% CI: 1.17–1.81; 95% PI: 0.70–3.04), from 6 studies.
- **Method as stated:** "random‐effects models", Stata v16. The estimator is not named, and Knapp–Hartung is not mentioned. The selection rule is "strongest anticholinergic exposure pattern" per study.
- **Where the components are:** Supplemental Table 4 (docx) gives every candidate estimate as text, 27 rows for the 6 studies. The single estimate actually pooled for each study appears only in the **Figure 2 forest-plot image**, which was not used.
- **What the text and tables pin down:**
  - The single-estimate studies are unambiguous: Ancelin, Hong and Park.
  - The selection rule points to >1095 TSDD for Coupland and Gray.
  - For **Richardson** it is ambiguous whether ACB-1, ACB-2 or ACB-3 >1460 DDD was used. The Results statement that study RRs "ranged from 1.05 to 2.63" matches only the ACB-1 value (1.05), but this is unverified. The rows mark it `uncertain`.
- **Validation checkpoints available:**
  - Supplemental Table 5 gives the **leave-one-out pooled RR and CI for each study** (6 rows). The pipeline can use these to identify the Richardson input without reading the image.
  - Table 2 gives the design subgroups (case-control and cohort).
  - The weakest-exposure sensitivity analysis gives 1.30 (1.14–1.48).

**3. Wang GH 2023, J Clin Med, PMC10573169, doi 10.3390/jcm12196342.** Antidepressants and dementia in older adults with depression.
- **Pooled estimate:** RR 1.21 (95% CI 1.12–1.29), I2 = 71%.
- **Method as stated:** a fixed-effect model within articles for multiple class or dose estimates, then a random-effects model between articles, in RevMan. The estimator is not named.
- **Components:** all component estimates and CIs are in Table 1 (11 rows). Two studies (Lee 2016 with 5 classes, Su 2020 with 2 doses) need the stated within-study fixed-effect step first. The intermediate within-study pooled values are only in the Figure 3 image.
- **Caveat:** Chatterjee 2015 contributes an "Anticholinergic ADs" estimate.

**Rejected as validation target, but pooled values kept as `ma_pooled` rows:**
- **Pieper 2020** (PMC7583519). Pooled ORs: any use 1.20 (1.09–1.32), ≥90 days 1.23 (1.17–1.29), ≥365 days 1.50 (1.22–1.85). These are random-effects inverse-variance pooled, in Stata 14.0. Study-level ORs are **only in the Figure 1 image**, and the appendix docx (aa-19-0341-file002_afaa090.docx) could not be retrieved (EPMC 500 ×2).

### Other meta-analyses screened and rejected (no rows)

- **Wang YC 2018** (Behav Neurol, PMC6079596; antidepressants):
  - The TCA-specific pooled RR 2.131 (1.427–3.184), with I2 and tau2 stated and computed in CMA v3, is described as coming from four studies.
  - Table 1 lists four TCA estimates as free text.
  - An informal mental re-pool of those four tabulated values did **not** appear to reproduce 2.131, so the tabulated values may not be the ones that were pooled.
  - Rejected for validation; usable only with caution.
- **Zhong 2015** (PLoS One, PMC4446315; benzodiazepines): Table 2 has subgroup and sensitivity pooled values only. Components are in the figure.
- **He 2019** (J Clin Neurol, PMC6325366; benzodiazepines, DerSimonian–Laird): components only in Figure 2.
- **Li 2019** (Medicine, PMC6408083; PPI): tables are images, with no text content in the XML.
- **Song 2019** (PLoS One, PMC6605652; PPI): Table 1 lists estimates. However, which estimate was pooled for Imfeld and Taipale (AD or vascular dementia) is ambiguous, and the pooled CI is symmetric on the linear scale.
- **Peng 2026** (Brain Sci, PMC12938945; PPI): Table 1 lists AD and dementia estimates for 18 studies. Which estimate was pooled is ambiguous, it includes RCTs, and some CIs look erroneous.
- **Zhang 2023** (Front Neurol, PMC10073491; anti-seizure drugs): Table 1 has no estimates.
- **Moraros-type check not possible:**
  - Xu 2017 (PMC5637833)
  - Dastidar (PMC13324608; no tables)
- **Du 2025** (statins, PMC11830700; DerSimonian–Laird): main components are only in figures; Table 2 covers dose and duration only.
- **Xie 2024** (DMARDs, PMC10900342): subgroup tables only.
- **Nerattini 2023** (PMC10625913) and **Xing 2022** (PMC9878281): no study-level estimate tables.

### Recommendation

For **validating the pooling code**, use **Moraros 2017**:
- Every component, its CI and its **weight** are tabulated.
- The software is stated (Stata `metan`, random effects).
- The text adds 12 extra pooled checkpoints with named subsets.

For an **on-topic anticholinergic check**, use **Dmochowski 2021**:
- Its supplementary Table S5 leave-one-out results let you test which candidate estimate per study reproduces the reported pooled RR.

Neither paper names the tau² estimator or says whether Knapp–Hartung was used. `metan`'s default random-effects method is DerSimonian–Laird, and Stata 16's `meta` default is REML, but this is software background knowledge, **not stated in either paper**.

## 5. Full query and download log (copied from `queries.log`)

```
2026-10-13T14:14:09Z  EPMC  TITLE:"anticholinergic" AND (TITLE:dementia OR TITLE:alzheimer) AND (OPEN_ACCESS:y)  hits=50
2026-10-13T14:14:27Z  EPMC  (anticholinergic) AND (dementia OR alzheimer) AND ("meta-analysis" OR "meta analysis") AND (OPEN_ACCESS:y) AND PUB_TYPE:"review"  hits=873 (too broad; not screened beyond first page)
2026-10-13T14:14:32Z  EPMC  (TITLE:anticholinergic OR TITLE:antimuscarinic OR TITLE:"anticholinergic burden") AND (TITLE:"meta-analysis" OR TITLE:"meta analysis" OR TITLE:"systematic review") AND HAS_FT:y  hits=43
2026-10-13T14:14:39Z  EPMC  (ABSTRACT:anticholinergic OR TITLE:anticholinergic) AND (TITLE:dementia OR TITLE:alzheimer OR TITLE:cognitive OR TITLE:cognition) AND (ABSTRACT:"meta-analysis" OR ABSTRACT:"pooled") AND HAS_FT:y  hits=14
2026-10-13T14:14:48Z  EPMC  AUTH:"Gray SL" AND TITLE:"anticholinergic" AND TITLE:"dementia"  hits=7
2026-10-13T14:14:49Z  EPMC  (TITLE:anticholinergic OR TITLE:antimuscarinic OR TITLE:"anticholinergic burden" OR TITLE:"anticholinergic drugs") AND (TITLE:dementia OR TITLE:alzheimer*) AND HAS_FT:y AND PUB_YEAR:[2010 TO 2026]  hits=88
2026-10-13T14:15:09Z  EPMC_fullTextXML  PMC7583519  http=200
2026-10-13T14:15:10Z  EPMC_fullTextXML  PMC7821204  http=200
2026-10-13T14:15:12Z  EPMC_fullTextXML  PMC8169439  http=500
2026-10-13T14:15:15Z  EPMC_fullTextXML  PMC4358759  http=500
2026-10-13T14:15:17Z  EPMC_fullTextXML  PMC5915701  http=200
2026-10-13T14:15:19Z  EPMC_fullTextXML  PMC6593623  http=200
2026-10-13T14:15:27Z  NCBI_efetch_pmc   PMC4358759  http=429 (rate limit; retried)
2026-10-13T14:15:29Z  NCBI_efetch_pmc   PMC8169439  http=200 (front matter only: publisher does not allow full-text XML)
2026-10-13T14:15:37Z  NCBI_efetch_pmc   PMC4358759  http=200 (author manuscript full text)
2026-10-13T14:17:24Z  EPMC_supplementaryFiles  PMC7583519  http=500
2026-10-13T14:18:56Z  EPMC_supplementaryFiles  PMC7821204  http=200
2026-10-13T14:19:46Z  EPMC  (ABSTRACT:"anticholinergic" OR ABSTRACT:"antimuscarinic") AND (ABSTRACT:"dementia") AND (ABSTRACT:"meta-analysis" OR ABSTRACT:"pooled") AND OPEN_ACCESS:y AND NOT TITLE:protocol  hits=6
2026-10-13T14:19:57Z  EPMC  (TITLE:anticholinergic OR TITLE:antimuscarinic OR TITLE:"overactive bladder" OR TITLE:"bladder") AND (TITLE:dementia OR TITLE:alzheimer OR TITLE:cognitive) AND (TITLE:meta OR TITLE:"systematic review" OR TITLE:"pooled")  hits=25
2026-10-13T14:21:24Z  EPMC_supplementaryFiles(retry)  PMC7583519  http=500
2026-10-13T14:21:25Z  EPMC_fullTextXML  PMC4541969  http=500
2026-10-13T14:21:51Z  EPMC  (TITLE:benzodiazepine* OR TITLE:"proton pump" OR TITLE:antidepressant* OR TITLE:antihistamine* OR TITLE:"sedative" OR TITLE:"z-drug*" OR TITLE:"hypnotic*" OR TITLE:"medication*" OR TITLE:"drug use") AND (TITLE:dementia OR TITLE:alzheimer*) AND (TITLE:"meta-analysis" OR TITLE:"meta analysis") AND OPEN_ACCESS:y  hits=20
2026-10-13T14:22:12-14:22:31Z  EPMC_fullTextXML  PMC4446315 PMC6325366 PMC6079596 PMC6408083 PMC6605652 PMC12938945 PMC10073491 PMC10573169 PMC5347943 PMC13324608 PMC5637833  all http=200
2026-10-13T14:23:25Z  EPMC  (anticholinergic OR antimuscarinic) AND (dementia OR "Alzheimer") AND ("pooled hazard ratio" OR "pooled odds ratio" OR "pooled risk ratio" OR "pooled relative risk" OR "pooled HR" OR "pooled OR" OR "summary relative risk") AND (TITLE:"meta-analysis" OR TITLE:"meta analysis" OR TITLE:"systematic review") AND OPEN_ACCESS:y  hits=108
2026-10-13T14:24:50Z  EPMC  (ABSTRACT:"anticholinergic" OR ABSTRACT:"antimuscarinic") AND (ABSTRACT:"incident dementia" OR ABSTRACT:"risk of dementia" OR ABSTRACT:"dementia risk" OR ABSTRACT:"Alzheimer") AND (ABSTRACT:"hazard ratio" OR ABSTRACT:"odds ratio" OR ABSTRACT:"HR" OR ABSTRACT:"OR") AND (ABSTRACT:cohort OR ABSTRACT:"case-control") AND OPEN_ACCESS:y AND NOT PUB_TYPE:review  hits=7
2026-10-13T14:25:02-14:25:47Z  EPMC_fullTextXML  PMC6802337 PMC7391507 PMC6494800 PMC5651434 PMC5383287 PMC7033580 PMC9005668 PMC10311860 PMC12538848 PMC12596250 PMC11580265 PMC12695687 PMC12022375 PMC13570850 PMC12274944  all http=200
2026-10-13T14:25:48Z  EPMC_fullTextXML  PMC8175463  http=500
2026-10-13T14:25:50Z  EPMC_fullTextXML  PMC9974875  http=500
2026-10-13T14:25:58Z  NCBI_efetch_pmc   PMC8175463  http=200 (front matter only: publisher does not allow full-text XML)
2026-10-13T14:26:00Z  NCBI_efetch_pmc   PMC9974875  http=200 (author manuscript full text)
2026-10-13T14:29:40Z  EPMC  TITLE:"polypharmacy and dementia" AND AUTH:"Park"  hits=1
2026-10-13T14:29:41Z  EPMC  TITLE:"Bladder antimuscarinics and cognitive decline"  hits=1
2026-10-13T14:29:41Z  EPMC  TITLE:"Cumulative antidepressant use and risk of dementia"  hits=1 (not in EPMC full text)
2026-10-13T14:29:42Z  EPMC  TITLE:"Antiparkinsonism anticholinergics increase dementia risk"  hits=1 (not in EPMC full text)
2026-10-13T14:29:43Z  EPMC  (TITLE:"anticholinergic") AND (TITLE:"dementia" OR TITLE:"Alzheimer" OR TITLE:"cognitive impairment") AND (AUTH:"Cai X" OR AUTH:"Carriere I" OR AUTH:"Jessen F" OR AUTH:"Papenberg G" OR AUTH:"Campbell N" OR AUTH:"Shah RC" OR AUTH:"Whalley LJ" OR AUTH:"Hsu WH" OR AUTH:"Welk B" OR AUTH:"Lockery JE")  hits=13
2026-10-13T14:29:59Z  EPMC_fullTextXML  PMC5215897  http=200
2026-10-13T14:30:02Z  EPMC_fullTextXML  PMC5408467  http=200
2026-10-13T14:30:04Z  EPMC_fullTextXML  PMC2933398  http=500
2026-10-13T14:30:07Z  NCBI_efetch_pmc   PMC2933398  http=200 (front matter only: publisher does not allow full-text XML)
2026-10-13T14:31:14Z  EPMC  (TITLE:dementia OR TITLE:alzheimer*) AND (TITLE:"meta-analysis") AND (TITLE:risk OR TITLE:association OR TITLE:incident) AND "DerSimonian" AND ("weight (%)" OR "Weight (%)" OR "% weight") AND OPEN_ACCESS:y AND (TITLE:drug* OR TITLE:medication* OR TITLE:use OR TITLE:users OR TITLE:exposure OR TITLE:therapy)  hits=3
2026-10-13T14:31:27Z  EPMC  (TITLE:dementia OR TITLE:alzheimer*) AND (TITLE:"meta-analysis" OR TITLE:"meta analysis") AND ("DerSimonian" OR "REML" OR "restricted maximum likelihood") AND OPEN_ACCESS:y AND (TITLE:drug* OR ... OR TITLE:anticholinergic*)  hits=14
2026-10-13T14:31:35-14:31:40Z  EPMC_fullTextXML  PMC11830700 PMC10625913 PMC10900342 PMC9878281  all http=200
```

Downloaded but not used for any CSV number (left in `staging_anticholinergic/`, not in `meta_anticholinergic/`):
- PMC5651434
- PMC5408467
- PMC4446315
- PMC6325366
- PMC6079596
- PMC6408083
- PMC6605652
- PMC12938945
- PMC10073491
- PMC13324608
- PMC5637833
- PMC11830700
- PMC10625913
- PMC10900342
- PMC9878281
- the NCBI front-matter stubs
