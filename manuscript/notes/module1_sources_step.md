# Module 1 — added instrument sources, validation and gene deserts (2026-10-21)

Researcher decisions of 2026-10-21, applied. **No row of any Alzheimer's
outcome file was read** (header line of FinnGen R13 G6_ALZHEIMER only). No MR,
liftover, filter or summary was run on AD data; `make mr` / `run_mr.R` and the
gene-desert null test were not run.

## Access terms (checked 2026-10-21)

| Source | Full data need | Done |
|---|---|---|
| UKB-PPP (Sun 2023) | Synapse account (anonymous download refused: "Anonymous users have only READ access permission"); AWS bucket is controlled access; portal ukb-ppp.gwas.eu unreachable | open Supplementary Table 9 (CC BY) used |
| deCODE (Ferkingstad 2021) | form: name, e-mail, institution, accept terms (scientific use, no redistribution), reCAPTCHA | open supplement checked: among targets only LGMN has cis pQTLs (untestable); not used |
| MetaBrain (de Klein 2023) | Google form: name, e-mail, institute, intended use; link e-mailed | open browser table (MOLGENIS API, website terms: no redistribution, cite, no re-identification) used, Cortex-EUR rank 1 |
| FinnGen R13 | public bucket, anonymous download works (finngen.fi access page refused automated access; it describes a form) | G6_ALZHEIMER downloaded, quarantined |
| Not adopted (options for the researcher only) | GTEx v8 brain-tissue cis-eQTLs (GTEx portal); INTERVAL plasma pQTLs (Sun et al. 2018). Access terms not checked here | not downloaded |

## Derived inputs

- `src/mr/prepare_qtl_sources.py`: UKB-PPP ST9 rows for targets (26; DBH cis 1,
  LGMN cis 1, rest trans), lifted 26/26, lifted position equals the hg37 position
  in the ST9 id for 26/26. MetaBrain: 8 Cortex-EUR rank-1 rows for targets,
  lifted 8/8; 2 palindromic (ADRA1A rs7015698, CACNA1D rs150032595) get no EAF
  under the rule fixed beforehand and are dropped.
- `src/mr/gene_deserts.py`: 458 deserts (shortest 667,564 bp), 357 variants.
- `src/mr/build_hmgcr_validation_inputs.py`: 7 SNPs; GLGC values equal the
  paper's Table 3.

## Results

Positive controls (`results/tables/mr_positive_controls.csv`,
`..._source_selection.csv`; cross-check 154/154 agree): no target passes in the
±100 kb window with any source; ADRA2A passes in the ±1,000 kb window with
MetaBrain rs4558098 (Wald 0.0217, SE 0.0096, p 0.024). HMGCR validation: PASS
(log OR 0.3676 [0.2155, 0.5197] vs published 0.3676 [0.2151, 0.5188]). Full
tables in prereg/mr.md section 0.
