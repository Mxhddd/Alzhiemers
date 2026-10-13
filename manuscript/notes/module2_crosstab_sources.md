# Open sources for amyloid × tau stage × age (item 8, 2026-10-13)

> **Search report only. Nothing here has been fitted** (instruction:
> "Do not fit anything on cross-table data yet"). Numbers below are only
> descriptive counts of what each table contains, read from the saved open
> full texts, to show what is there. Searches used:
> - Europe PMC REST;
> - NCBI E-utilities;
> - the Semantic Scholar API without a key (one call succeeded after rate
>   limiting);
> - one institutional repository.
>
> No email address was used anywhere. Query log:
> `data/extracted/candidates/queries_2026-10-13.txt`.

## 1. Braak et al. 2011 (stage-by-age tables): still not obtained

| Route | Result |
|---|---|
| Europe PMC (DOI and PMID 22002422) | Record found, **no full text**, isOpenAccess = N, no PMCID |
| NCBI PMC | No PMCID, so no copy |
| Semantic Scholar | `isOpenAccess: true`, status **BRONZE** (free to read at the publisher, **no licence**). PDF link: academic.oup.com/jnen/article-pdf/70/11/960/17418183/70-11-960.pdf |
| That PDF link | **Blocked: HTTP 403, Cloudflare "Just a moment..." challenge.** One plain request was made and nothing was done to get around it |
| Ulm University repository (OPARU) | Search URL returned HTTP 404 |
| CORE API | Requires an API key obtained with an email address: **not used** (blocked by the no-email rule) |

The publisher's free-to-read PDF exists. A person with a normal browser can
download it without an institutional login. Saved as
`data/raw/papers/braak2011.pdf`, it would unblock the stage-by-age table.

## 2. Sources with amyloid and tau stage *per person* (individual-level), open access

These are the most useful finds. Each lists every brain with age, tau stage
and amyloid stage, so an amyloid × tau-stage × age cross-table can be
*counted* from them, with no modelling or interpolation.

| Source | Licence | What is tabulated | N, ages | Locator |
|---|---|---|---|---|
| Braak H, Del Tredici K. *J Neuropathol Exp Neurol* 2020 (PMC6970449) | CC BY-NC 4.0 | Per case: sex, **age**, **NFT stage** (I–VI), **Aβ (Thal) phase** (0–5), APOE, AGD, TDP-43, HS | 95 cases, ages 34–97; NFT I 8, II 22, III 23, IV 11, V 28, VI 3; Aβ phase 0 in 39 | the single TABLE ("Demographic and neuropathological data from the 95 cases") |
| Braak H, Del Tredici K. *Cereb Cortex* 2018;28:3372 (PMC6095209) | CC BY-NC 4.0 | Per case: **age**, sex, **NFT stage** (0–IV), **Aβ phase**, α-syn, layer-pattern codes | 54 rows (one with blanks), ages 12–86; mostly NFT 0–I, Aβ 0 in 47 | Table 1 |
| Elobeid A, Soininen H, Alafuzoff I. *Acta Neuropathol* 2012;123:97 (PMC3249177) | CC BY-NC 4.0 | Per case: **age**, sex, APOE, cause of death, **cortical β-amyloid (+/−)**, Braak stage, and **locus coeruleus tau (pretangles / neurites / tangles)** | 31 subjects with hyperphosphorylated tau, ages 38–50 (from the 95 examined) | Table 4 |

Caveat: these are selected series (consecutive autopsies at one centre,
or cases chosen for a staging study), not population samples. NFT staging
uses I–VI, plus 0 in the 2018 table. The 2011 paper's pretangle stages
a–c / 1a–1b do not appear as separate categories in the 2018 and 2020
tables. Elobeid 2012 reports LC pretangles directly.

## 3. Aggregate cross-tables (counts or percentages by group), open access

| Source | Licence | What is cross-tabulated | N | Locator |
|---|---|---|---|---|
| Braak H, Mayer B, Feldengut S, Schön M, Del Tredici K. *Acta Neuropathol* 2025;149 (PMC12102137) | CC BY 4.0 | **Age group (decades) × NFT stage I–III in cases WITHOUT Aβ deposits**: the amyloid-negative stratum only | 308 | Table 1 |
| Alafuzoff I, Libard S. *J Alzheimers Dis* 2026;110:397 (PMC12960740) | CC BY 4.0 | **Age group (50–59 … ≥90) × combinations of altered proteins** (HPτ, Aβ, α-syn, TDP-43), overall and by dementia status | 1,825 (1,453 non-demented) | Tables 2–3 |
| Alafuzoff I, Libard S. *Int J Mol Sci* 2024;25:4065 (PMC11012171) | CC BY 4.0 | **Age group × PART (tau without Aβ) vs ADNC (tau with Aβ)**, with lesion percentages | 1,184 | Table 5 |
| Elobeid A, Libard S, Leino M, Popova SN, Alafuzoff I. *J Neuropathol Exp Neurol* 2016;75:316 (PMC4793886) | CC BY-NC 4.0 | Age group × **Braak stage including a,b** (Table 5); age group × **Thal phase** (Table 7). Margins only, not crossed with each other; cognitively unimpaired | as in tables | Tables 5, 7 |
| Nelson PT et al. *J Neuropathol Exp Neurol* 2012 (PMC3560290) | NIH author manuscript (free, not CC) | **Braak stage × neuritic Aβ plaque density** (no age) | as in table | Table 4 |
| Neltner JH et al. *Neurobiol Aging* 2016 (PMC4688098) | NIH author manuscript | **Braak stage × CERAD** in centenarian / extreme-old cohorts (no age bands) | combined cohort | Table 2 |
| Pletnikova O et al. *Cogn Behav Neurol* 2015 (PMC4588069) | NIH author manuscript | Aβ status by APOE and other factors, ages 30–50 | 95 aged 40–50 (154 total) | Tables 4, 6 |
| Savola S et al. *Brain Pathol* 2022 (PMC9305229) and *Alzheimers Dement* 2025 (PMC12371459) | CC BY 4.0 | PART vs low/high ADNC with CERAD/Thal distributions, ages 85+ only (Vantaa 85+, population-based) | 301 / 210 | Tables 1–2 / 1 |
| Wharton SB et al. *Acta Neuropathol Commun* 2016 (PMC4746919) | CC BY 4.0 | Braak NFT stage distributions (CFAS, population-based); no amyloid cross | CFAS | Tables 1–5 |
| Suemoto CK et al. *PLoS Med* 2017 (PMC5369698) | CC BY 4.0 | Lesion frequencies by dementia status (Brazil); no age × stage × amyloid cross | 1,092 | Tables 2–3 |

## 4. What this would allow (not done)

- The three individual-level tables in §2 make it possible to count
  amyloid-positive vs amyloid-negative brains at each tau stage within
  age bands. This is the identifying design (option (b) in
  `module2_identifiability.md`): it compares brains *of the same age*.
- Together they hold 180 people (95 + 54 + 31). That is small and from selected
  series, so any test from them would be weak. They are not population
  samples.
- Braak 2025 gives the tau-stage-by-age distribution in amyloid-negative
  brains only.
- Alafuzoff 2024/2026 give amyloid ± by tau ± by age at large N (1,184 /
  1,825), but without tau *stage*.
- Per the instruction, nothing has been extracted into `data/extracted/`
  and nothing has been fitted. Extraction (with verified: false), a model
  and a pre-registration would come first.
