# data/raw/

Downloaded or hand-entered inputs. Files here are **never edited and never
committed** (see `.gitignore`); only this README is tracked.

Files that the researcher transcribes by hand — code must never create,
scrape or fill them in. If one is missing, the pipeline stops and says so.

| File | Used by | Provided by |
|---|---|---|
| `braak_counts.csv` (`age_lower, age_upper, stage, n, source`) | Module 2 | researcher, hand-transcribed |
| `extraction.csv` (study id, drug class, drug, comparator_type, comparator, effect measure, estimate, CI bounds, covariates, 7 ROBINS-I domains) | Module 5 | researcher, hand-extracted |

Downloaded files (GWAS summary statistics, QTLs, connectomes) will be listed
here, with the URL, access date and checksum, by the session that adds the
download script.
