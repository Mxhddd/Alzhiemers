# data/raw/

Downloaded or hand-entered inputs. Files here are **never edited and never
committed** (see `.gitignore`); only this README is tracked.

Files that the researcher transcribes by hand — code must never create,
scrape or fill them in. If one is missing, the pipeline stops and says so.

| File | Used by | Provided by |
|---|---|---|
| `braak_counts.csv` (`age_lower, age_upper, stage, n, source`) | Module 2 | researcher, hand-transcribed |
| `validation/mr_ldl_hmgcr_cad_exposure.csv`, `validation/mr_ldl_hmgcr_cad_outcome.csv` (`SNP, effect_allele, other_allele, beta, se, eaf, source`) | Session 1 MR validation | researcher, from the reproduced paper's supplement or the GWAS it names |
| `extraction.csv` (study id, drug class, drug, comparator_type, comparator, effect measure, estimate, CI bounds, covariates, 7 ROBINS-I domains) | Module 5 | researcher, hand-extracted |

Module 1 inputs (paths and column maps go in `config/mr_targets.yaml`): the
Bellenguez 2022 stage-1 file (GCST90027158, GRCh37), eQTLGen cis-eQTL and
allele-frequency files, MetaBrain / pQTL files where used, a resting-heart-rate
GWAS for the ADRB1 positive control, a 1000 Genomes EUR PLINK fileset for LD
clumping, optional proxy-free and FinnGen files, and
`validation/gene_desert_variants.csv` for the null test.

Downloaded files (GWAS summary statistics, QTLs, connectomes) will be listed
here, with the URL, access date and checksum, by the session that adds the
download script.
