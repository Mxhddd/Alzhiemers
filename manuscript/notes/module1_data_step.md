# Module 1 — data step (2026-10-15)

This note covers the downloads, the derived inputs, the liftover check, the LD
panel, the validation run, the positive-control gate and what is still open.
Scripts: `src/mr/download_inputs.sh`, `src/mr/prepare_inputs.py`,
`src/mr/liftover.py`, `src/mr/run_positive_controls.R`,
`src/mr/pipeline_check_hmgcr.R`. Config: `config/mr_targets.yaml` (each new
value has a comment naming its source).

**No row of any Alzheimer's or dementia outcome file was read.** The three AD
files (Bellenguez 2022, Kunkle 2019, FinnGen R10 G6_ALZHEIMER) were
downloaded into `data/raw/mr/quarantine_outcome/`. They were checksummed, and
only their header line (`zcat … | head -1`) was looked at, to write the column
maps. No MR, liftover, filter or summary was run on them. `make mr` /
`src/mr/run_mr.R` was not run. `prereg/mr.md` is not frozen and was not
edited. The only other AD-related files read were study-level metadata (GWAS
Catalog `-meta.yaml`, `README.txt`, `md5sum.txt`, REST API study records, and
the FinnGen manifests). The Bellenguez file's md5 matches its `md5sum.txt`
(`9d23b9ba…`).

## 1. Downloads

All downloads are in `data/raw/mr/` (git-ignored). The full manifest is in
`data/raw/mr/MANIFEST.tsv`, with path, URL, bytes, sha256, download time and
purpose. `bash src/mr/download_inputs.sh` re-creates it; the script resumes
partial files. No login, token or e-mail address was used anywhere.

| File (under data/raw/mr/) | Source | Bytes | sha256 (first 12) | Purpose |
|---|---|---|---|---|
| eqtlgen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | molgenis26.gcc.rug.nl eQTLGen | 3,880,187,026 | 1133a89fe33b | exposure: full cis-eQTLs, blood, hg19 |
| eqtlgen/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz | eQTLGen | 240,045,342 | edffd4309459 | exposure allele frequencies |
| eqtlgen/README_cis, README_allele_frequency | eQTLGen | 1,645 / 723 | 9f61a616f627 / 029fd7769631 | column definitions |
| gwascatalog/GCST007609/ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.rhrmean.assoc.gz (+ readme.txt) | GWAS Catalog FTP | 149,492,529 | 40433f4dea1e | control: resting heart rate |
| gwascatalog/GCST006630/Evangelou_30224653_DBP.txt.gz (+ README) | GWAS Catalog FTP | 135,219,424 | c55cfc8b50ab | control: diastolic BP |
| gwascatalog/GCST006624/Evangelou_30224653_SBP.txt.gz (+ README) | GWAS Catalog FTP | 136,130,586 | 71d5c8af9687 | control: systolic BP |
| gwascatalog/GCST90002232/GCST90002232_buildGRCh37.tsv.gz (+ meta, README, md5) | GWAS Catalog FTP | 499,677,827 | 030182eb8d93 | control: fasting glucose (md5 matches `md5sum.txt`) |
| gwascatalog/GCST007561/Dashti_30846698_sleepdurationsumstats.txt.zip (+ README) | GWAS Catalog FTP | 378,315,754 | 77ebce3fcabb | control: sleep duration |
| gwascatalog/GCST002222/jointGwasMc_LDL.txt.gz | GWAS Catalog FTP | 58,851,380 | ddf75ed591ff | LDL (HMGCR check / validation inputs) |
| gwascatalog/GCST003116/cad.add.160614.website.txt (+ readme) | GWAS Catalog FTP | 791,997,790 | 8b6773782707 | CAD (HMGCR check / validation inputs) |
| liftover/hg38ToHg19.over.chain.gz | hgdownload.soe.ucsc.edu | 1,246,411 | 14a712e8e147 | GRCh38 → GRCh37 chain |
| ld_1kg_eur/g1000_eur.zip (unzipped to g1000_eur.bed/bim/fam) | MAGMA reference data (vu.data.surfsara.nl, linked from cncr.nl/research/magma) | 511,626,945 | 83a48fd9dcaa | LD panel |
| finngen_r10/R10_manifest.tsv | storage.googleapis.com/finngen-public-data-r10 | 711,459 | f5a97e7aad28 | FinnGen endpoint counts |
| quarantine_outcome/GCST90027158/GCST90027158_buildGRCh38.tsv.gz (+ meta, README, md5) | GWAS Catalog FTP | 755,201,909 | 013134d22a2a | **AD outcome, not read** |
| quarantine_outcome/GCST007511/Kunkle_etal_Stage1_results.txt | GWAS Catalog FTP | 569,380,789 | ae280a70eb82 | **AD proxy-free sensitivity, not read** |
| quarantine_outcome/finngen_R10_G6_ALZHEIMER/finngen_R10_G6_ALZHEIMER.gz | FinnGen R10 public bucket | 804,900,131 | 56b53bf83764 | **AD replication, not read** |

Total is about 12 GB on disk, including the unzipped 3.5 GB LD panel.

**Download problems:**
- The MRC-IEU LD panel at `fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz` could not
  be used. Over HTTPS the certificate has expired; plain HTTP is refused by
  the proxy (405). TLS verification was not turned off, and the MAGMA copy
  of the same 1000G Phase 3 EUR data was used instead.
- The sleep-duration zip was cut off once (curl exit 56). The resumable
  script finished it on the second run.

**Tools installed** (environment only, not in the lockfiles):
- PLINK v1.90b7.2 (11 Dec 2023), from the official
  `s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip` (zip
  sha256 `1162eca2…91bc3`), installed as `/usr/local/bin/plink`.
- `pyliftover==0.4.1`, in a separate venv
  `/root/.local/share/lcfirst-liftover-venv`. It was not added to
  `uv.lock`, to avoid touching lockfiles while another session works.

## 2. Derived inputs (`src/mr/prepare_inputs.py`)

All derived files are in `data/interim/mr/` (git-ignored). Original column
names are kept and new columns are appended. Each file is restricted to the
target genes' regions: Ensembl GRCh37 gene body ±1000 kb (the pre-registered
sensitivity window), plus the gene's eQTLGen cis range. This contains every
possible instrument; the pipeline then applies its own ±100 kb or ±1000 kb
window.

- **eQTLGen** → `eqtlgen_cis_targets_with_eaf.tsv`.
  - 127,341,798 rows were read; 50,058 belong to the 11 genes (10 targets +
    HMGCR).
  - The eaf was joined by rsID with an allele check: 42,338 rows with the
    AF file's AlleleB = assessed allele, 7,720 with AlleleA = assessed
    (eaf = 1 − AF), and 0 mismatches.
  - Build check: SNPPos equals the AF file's hg19 position for all 50,058
    rows.
  - **SLC6A2, ADRA1A, MAOA and HCRTR2 are not in eQTLGen.** This was checked
    by symbol and by Ensembl ID over the whole file: 19,250 genes, autosomes
    only.
- **Control GWAS.** Rows read → kept in regions → written:

| Trait | Read | In regions | Written | Notes |
|---|---|---|---|---|
| resting heart rate | 5,264,923 | 43,437 | 42,531 | rsID from 1000G EUR (906 unmapped dropped). The file gives MAF only: A1 frequency is oriented with the 1000G EUR frequency for 36,582 variants; 5,949 palindromic variants get eaf NA, and the reader drops them |
| diastolic BP | 7,160,657 | 58,940 | 58,334 | no rsIDs in the file: mapped from 1000G EUR by chr, pos and allele pair (606 unmapped dropped) |
| systolic BP | 7,088,121 | 58,387 | 57,786 | as for DBP (601 unmapped dropped) |
| fasting glucose | 34,064,006 | 276,661 | 148,119 | rsID mapping drops 128,542 variants not in 1000G EUR (many are non-EUR or rare). Eaf as given for 111,963; filled from 1000G EUR for 30,730 non-palindromic NA rows; 5,426 palindromic NA rows left NA |
| sleep duration | 14,661,601 | 117,627 | 117,627 | own rsIDs |
| LDL (GLGC 2013) | 2,437,751 | 21,754 | 21,754 | own rsIDs, hg19 from SNP_hg19 |
| CAD (CARDIoGRAMplusC4D) | 9,455,778 | 76,394 | 76,394 | own rsIDs (markername) |

The 1000G EUR panel has 174,337 variants in the regions.

**Gene coordinates:** `data/interim/gene_regions_GRCh37.csv`, from Ensembl
GRCh37 via biomaRt, as `src/mr/regions.R` does it, for all 10 targets plus
HMGCR. The archive reports its release as "Ensembl Genes 116". biomaRt gave
transient 504/405 errors, and retries succeeded.

## 3. Liftover

`src/mr/liftover.py` uses the UCSC hg38ToHg19 chain. It drops (and logs, by
reason) any position that is unmapped, maps to several places, maps to the
minus strand or maps to another chromosome. It refuses any file under
`quarantine_outcome/` unless `python -m src.common.prereg check mr` passes.

**Check on non-AD data.** I took the first 299,999 chr1 rows of FinnGen R10
`HEIGHT_IRN` (GRCh38; a byte-range slice in the scratchpad, not kept).
- Lifted: 299,231 written. Dropped: 458 unmapped, 306 opposite strand, 4 on
  a different chromosome.
- Of the 176,230 rsIDs also in 1000G EUR, **176,225 (99.997%) got exactly
  the panel's GRCh37 position**; 5 differ.

**No analysis file needed liftover.** Every exposure and control file chosen
is GRCh37. The only GRCh38 files are the quarantined AD outcomes, and lifting
them means reading their rows.

## 4. LD panel

`instruments.ld_reference_bfile = data/raw/mr/ld_1kg_eur/g1000_eur`:
- 1000 Genomes Phase 3 EUR, 503 samples, GRCh37;
- bi-allelic SNPs with rsIDs, no duplicate IDs;
- 22,665,064 variants, chromosomes 1–23;
- the MAGMA "European" reference.

Clumping uses PLINK 1.9 (`mr_plink_r2_fun`). Three candidate instruments were
dropped as absent from the panel (ADRB1 ×2, DBH ×1); they are logged.

## 5. Choices made (also in config comments)

- **Resting heart rate: GCST007609** (Zhu 2019, UKB BOLT-LMM, n = 458,969,
  GRCh37). It is the only European resting-heart-rate GWAS on the GWAS
  Catalog with full GRCh37 statistics. The larger MVP files (e.g.
  GCST90476338) are GRCh38.
- **DBP / SBP: Evangelou 2018** (GCST006630 / GCST006624, n = 757,601,
  GRCh37). The Keaton 2024 files named in the 2026-10-13 control decision are
  GRCh38. Evangelou is the largest open GRCh37 BP GWAS.
- **Fasting glucose: MAGIC Chen 2021**, European (GCST90002232,
  n = 200,622, GRCh37).
- **Sleep duration: Dashti 2019** (GCST007561, n = 446,118, GRCh37), the
  study named in the control decision.
- **LD panel:** MAGMA 1000G EUR, because the MRC-IEU file was unreachable
  (§1).
- **Bellenguez `outcome.n_cases = 85,934`**:
  - this is 39,106 clinically diagnosed + 46,828 proxy cases;
  - `n_controls = 401,577`;
  - source: GWAS Catalog API study record. The total, 487,511, equals the
    file's metadata.
- **FinnGen replication = R10 G6_ALZHEIMER** (10,520 cases / 401,661
  controls, from the R10 manifest), as the brief asked. See §8 for R13.
- **R yaml reads a bare key `n` as `FALSE` (YAML 1.1 boolean).**
  - This broke the `n` column map for eQTLGen, and the existing
    `positive_control.n` field was silently parsed as `FALSE` too.
  - Fix: every `n:` key in `config/mr_targets.yaml` is now quoted (`"n":`).
    No code change was needed.

## 6. Validation (`make validate-mr`, 2026-10-15 10:31 UTC) — FAILED

| Target | Result | Why |
|---|---|---|
| mr_ldl_hmgcr_cad | FAIL (error) | Target is incomplete: the published citation, estimate, CI, conversion factor and `inputs.genome_build` are RESEARCHER fields "Claude must not fill these in"; the variant-level input CSVs in `data/raw/validation/` are researcher-provided (data/raw/README.md). No comparison with a published value was possible. |
| mr_null_gene_desert | FAIL (error) | `inputs.genome_build` and `inputs.desert_variants` are RESEARCHER fields. In addition, by design this test reads the **AD outcome** rows (`mr_cfg$outcome$file`), which is quarantined until the freeze; `outcome.file` is null, so it cannot run before then. |

Nothing was tuned or filled in to make either target pass.

**Informal machinery check (not a validation):**
`src/mr/pipeline_check_hmgcr.R` runs the Module 1 pipeline on LDL (GLGC 2013,
HMGCR ±100 kb) → CAD (CARDIoGRAMplusC4D 1000G).
- One instrument, rs12916, with F = 372, harmonised by allele swap.
- Wald ratio b = 0.489 (SE 0.129), OR 1.63 per SD of LDL, 95% CI 1.27–2.10,
  p = 1.6 × 10⁻⁴.
- The positive sign is the one statin biology predicts.
- TwoSampleMR cross-check: 7/7 quantities agree.
- Output: `results/tables/mr_pipeline_check_hmgcr_ldl_cad.csv`.

No published number was compared.

## 7. Positive controls (`Rscript src/mr/run_positive_controls.R`)

Gate: p < 0.05 and the expected sign. Exposure: eQTLGen. Primary window ±100
kb. The ±1000 kb sensitivity run gave the same instruments and the same
estimates. Output: `results/tables/mr_positive_controls.csv`, the TwoSampleMR
cross-check in `results/tables/mr_positive_controls_crosscheck.csv` (112/112
quantities agree), and drops in `results/logs/mr_positive_controls_dropped.md`.

| Target | Control trait | Expected sign | Instruments (F min) | Method | b (SE) [95% CI] | p | Gate |
|---|---|---|---|---|---|---|---|
| ADRB1 | resting heart rate (bpm) | + | rs4917675, rs117648943 (42.9) | IVW | −0.0104 (0.260) [−0.519, 0.498] | 0.97 | **FAIL** (null) |
| ADRA2A | fasting glucose (mmol/L) | + | rs11195432 (39.7) | Wald ratio | 0.0401 (0.0365) [−0.031, 0.112] | 0.27 | **FAIL** (sign right, not significant) |
| DBH | systolic BP (mmHg) | + | rs2519131, rs2427998 (58.5) | IVW | −0.303 (0.173) [−0.642, 0.037] | 0.081 | **FAIL** (sign opposite, not significant) |
| CACNA1D | systolic BP (mmHg) | + | rs9830632 (94.9) | Wald ratio | −2.43 (0.396) [−3.20, −1.65] | 8.9 × 10⁻¹⁰ | **FAIL** (strongly significant, sign opposite) |
| SLC6A2 | resting heart rate | − | — | — | — | — | not run: gene not in eQTLGen |
| ADRA1A | diastolic BP | + | — | — | — | — | not run: gene not in eQTLGen |
| HCRTR2 | sleep duration | − | — | — | — | — | not run: gene not in eQTLGen |
| LGMN, SLC18A2, MAOA | — | — | — | — | — | — | untestable (no control pre-registered); MAOA also absent from eQTLGen |

Estimates are per SD of blood expression. The units come from each control
GWAS.

**Reading.** No target passes its gate with blood eQTL instruments, so under
the pre-registered rule every target would be excluded before the AD outcome
is read. Likely reasons (no instrument or threshold was changed in response):
- Blood expression may not proxy target activity in heart, vasculature,
  pancreas or brain. This is the brief's stated limitation.
- Instrument counts are small: 1–2 after clumping at r² < 0.001.
- **CACNA1D:** rs9830632 sits in a dense 3p21 BP locus. A strong
  opposite-sign result is more consistent with LD/pleiotropy with a
  neighbouring gene, or a tissue-discordant eQTL, than with Cav1.3 biology.
  Colocalisation was not run here.
- **ADRB1:** the result is null; the heart-rate signal at ADRB1 may not be
  expression-mediated in blood.

## 8. What still needs the researcher

1. **Validation target `mr_ldl_hmgcr_cad`.** Choose the paper and transcribe
   the published fields: citation, DOI, locator, method, scale, estimate, CI,
   conversion factor, build. Supply `data/raw/validation/mr_ldl_hmgcr_cad_
   {exposure,outcome}.csv`. The GLGC LDL and CAD GWAS the inputs could be
   extracted from are now downloaded (§1).
2. **Validation target `mr_null_gene_desert`.**
   - Fill `genome_build` and `desert_variants`, or approve a script to make
     the desert list.
   - **Design conflict:** the null test reads the AD outcome, but the AD
     outcome must stay unread until the freeze, while the validation must pass
     before the run. Either allow the null test after the freeze, or switch it
     to a non-AD outcome (e.g. CAD).
3. **AD outcome on GRCh38.** The GWAS Catalog Bellenguez file is GRCh38
   (file name and metadata), not GRCh37 as brief §5 says.
   - `outcome.file` stays null.
   - After the freeze, lift it with `src/mr/liftover.py` (the guard then
     allows it) and point `outcome.file` at the lifted copy.
   - The column map is filled from the header.
   - The same applies to the FinnGen replication file (`replication.file`
     null; columns and N filled).
4. **Proxy-free file (Kunkle 2019).** Its header has no allele-frequency
   column and different names from `outcome.columns`, which run_mr.R reuses
   for it. It needs a reformatted copy with an EAF after the freeze.
   `proxy_free_file` stays null.
5. **FinnGen release.** prereg §4 says "release current at run time". That is
   R13 (G6_ALZHEIMER 16,652 / 483,534 in the R13 manifest), not R10. Decide
   which one.
6. **Exposures for SLC6A2, ADRA1A, HCRTR2 and MAOA** (not in eQTLGen).
   - Only MetaBrain (cortex) or pQTL sources could instrument them.
   - The config `metabrain` and `pqtl` blocks are still null: file, columns
     and build. Not downloaded in this step.
7. **Positive controls all failed (§7).** Every target is excluded under the
   pre-registered gate. Decide whether that stands as the result (the gate
   working) or whether the controls or exposures need revisiting before the
   freeze. A change must be made before any AD data are read.
8. **`R.utils` is not in the renv library.** `data.table::fread` needs it to
   read `.gz` files such as the Bellenguez `.tsv.gz`. The derived files here
   are uncompressed, so nothing failed yet. Add it before the main run, or
   decompress the lifted outcome.
9. **`pyliftover`** is not in `uv.lock`; add it if the liftover is to run
   under `uv`.
10. **Possible code issue, not changed** (`src/mr/targets.R`,
    `mr_primary_row`). When base R and TwoSampleMR disagree, the status
    becomes "estimated (TwoSampleMR value used…)", and `mr_primary_row` then
    stops with "no primary MR estimate available". It did not trigger here:
    every comparison agreed.
11. **Alleles outside A/C/G/T.** `mr_harmonise` stops on any allele outside
    A/C/G/T anywhere in the outcome table, and `mr_read_sumstats` does not
    drop such rows. The CAD file has them; `pipeline_check_hmgcr.R` drops and
    logs them itself. Check the AD files' allele coding after the freeze.

Still null in `config/mr_targets.yaml`:
- `outcome.file`, `outcome.proxy_free_file`;
- `replication.file`;
- `exposures.metabrain.{file,columns,build}`, `exposures.pqtl.{file,columns,build}`;
- `positive_control: null` for LGMN, SLC18A2 and MAOA, which is the
  pre-registered "untestable".
