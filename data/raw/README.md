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

Downloaded files (GWAS summary statistics, QTLs, connectomes) are listed
here, with the URL, access date and checksum, by the session that adds the
download script.

## Module 1 downloads (`data/raw/mr/`, 2026-10-15)

Script: `bash src/mr/download_inputs.sh` (resumable; rebuilds
`data/raw/mr/MANIFEST.tsv` with path, URL, bytes, sha256, UTC download time and
purpose). Choices, derived files and checks: `manuscript/notes/module1_data_step.md`.
Files under `quarantine_outcome/` are Alzheimer's outcomes: downloaded and
checksummed only, **not read** (header line only) until `prereg/mr.md` is frozen.

| Path (under data/raw/mr/) | URL | Bytes | sha256 |
|---|---|---|---|
| eqtlgen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | 3880187026 | `1133a89fe33b795dde7e480e39841fbd0787bfdb5d18ef19058fa749a936f814` |
| eqtlgen/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz | https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz | 240045342 | `edffd43094598d08db6b423f7b3ce68b000a02f0b993ea4cfcf1f61a8175485d` |
| eqtlgen/README_cis | https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_cis | 1645 | `9f61a616f6276de6d45a7e35577b50596b95d3e76e2ad3c3932cb18f10cd6bdc` |
| eqtlgen/README_allele_frequency | https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_allele_frequency | 723 | `029fd77696317478f4ca0c232cd5ad3150f14d3883fd614e6490a0c60ed1144e` |
| gwascatalog/GCST007609/ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.rhrmean.assoc.gz | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007609/ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.rhrmean.assoc.gz | 149492529 | `40433f4dea1edaa1cfc764f47d91e318b530af76eb88fe5f3af9dfb6fc4bd993` |
| gwascatalog/GCST007609/readme.txt | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007609/readme.txt | 494 | `daa2d762d294f98d001fca010654881e0ddae5c10164a1f18a64932c4d6a6179` |
| gwascatalog/GCST006630/Evangelou_30224653_DBP.txt.gz | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006630/Evangelou_30224653_DBP.txt.gz | 135219424 | `c55cfc8b50abee8e13ef88d70ab52678b1f45b45cf5974bdab91f8c8b361fd59` |
| gwascatalog/GCST006630/README_Evangelou_30224653.txt | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006630/README_Evangelou_30224653.txt | 875 | `803144b2980ede246904d345491b339b0a609b53d0fba47de0c18b6196476e9d` |
| gwascatalog/GCST006624/Evangelou_30224653_SBP.txt.gz | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006624/Evangelou_30224653_SBP.txt.gz | 136130586 | `71d5c8af9687faaa1b941fbebb49b6907e3ae0513162dac80b05d452b34133f2` |
| gwascatalog/GCST006624/README_Evangelou_30224653.txt | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006624/README_Evangelou_30224653.txt | 875 | `803144b2980ede246904d345491b339b0a609b53d0fba47de0c18b6196476e9d` |
| gwascatalog/GCST90002232/GCST90002232_buildGRCh37.tsv.gz | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002232/GCST90002232_buildGRCh37.tsv.gz | 499677827 | `030182eb8d93d4e4c7c3f7026a3c50bdb539cfaceebd735673195091e476b047` |
| gwascatalog/GCST90002232/GCST90002232_buildGRCh37.tsv.gz-meta.yaml | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002232/GCST90002232_buildGRCh37.tsv.gz-meta.yaml | 1437 | `5165d2613b23ff1ac032d448c8b3bce19c27b55a94c474d5891a0a910a6f73a5` |
| gwascatalog/GCST90002232/README.txt | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002232/README.txt | 688 | `c7025cdb3f69f5fd525bae9a149657be62ee67e145abeacf584e845a7ed2fd89` |
| gwascatalog/GCST90002232/md5sum.txt | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002232/md5sum.txt | 228 | `3c15f2e587dad282921f34d2069775262b09dcdff9231edb6f77f43281d37b04` |
| gwascatalog/GCST007561/Dashti_30846698_sleepdurationsumstats.txt.zip | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007561/Dashti_30846698_sleepdurationsumstats.txt.zip | 378315754 | `77ebce3fcabb32a36e028390c7a34ce7f6aac11e2b663014132be1dc0f82d7c8` |
| gwascatalog/GCST007561/Dashti_30846698_fullUKBB_Sleepduration_summary_stats_README.txt | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007561/Dashti_30846698_fullUKBB_Sleepduration_summary_stats_README.txt | 489 | `cfce44b3e09f783da29a83d936cbd3e35504073a12bf4083d0be53e032dc32b3` |
| gwascatalog/GCST002222/jointGwasMc_LDL.txt.gz | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST002001-GCST003000/GCST002222/jointGwasMc_LDL%20(1).txt.gz | 58851380 | `ddf75ed591ffcce8c3e81ae2f4c77f4428cd7dbf4628b24467aa8cf6e85e3dce` |
| gwascatalog/GCST003116/cad.add.160614.website.txt | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003116/cad.add.160614.website.txt | 791997790 | `8b677378270772fbd49c140b75acddb0817b1f321d0cdaca48e91c88758b6f4b` |
| gwascatalog/GCST003116/cad.add.readme | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003116/cad.add.readme | 237 | `4addf882b7c0e01b405bfd52d39d575bd5ad933cf1e63553bea9d327178d5f3b` |
| liftover/hg38ToHg19.over.chain.gz | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz | 1246411 | `14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1` |
| ld_1kg_eur/g1000_eur.zip | https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download | 511626945 | `83a48fd9dcaa0b9a874b18c63143a4ede93f05505b215b0bd8790130a0d7a954` |
| finngen_r10/R10_manifest.tsv | https://storage.googleapis.com/finngen-public-data-r10/summary_stats/R10_manifest.tsv | 711459 | `f5a97e7aad281ec1df05308d58b8a23587b2bc6568992be991c0a2221111ff63` |
| quarantine_outcome/GCST90027158/GCST90027158_buildGRCh38.tsv.gz | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/GCST90027158_buildGRCh38.tsv.gz | 755201909 | `013134d22a2a0c2dd0cb683192d0c109da1aec83ba54d94a837a8b71bcbc4a1f` |
| quarantine_outcome/GCST90027158/GCST90027158_buildGRCh38.tsv.gz-meta.yaml | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/GCST90027158_buildGRCh38.tsv.gz-meta.yaml | 1238 | `f27c635833b4200a0790c247d5d6e44f4352461a8317ac168a4a829cd7152151` |
| quarantine_outcome/GCST90027158/README.txt | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/README.txt | 526 | `57cdd26f9b69bcfb17094c6882162fd72bdc43929e0b805b8ea8998296c34540` |
| quarantine_outcome/GCST90027158/md5sum.txt | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/md5sum.txt | 140 | `92d997145d497fdc217bed9a7c176a86c5b48c624e86e219a21fc6e7257eeaea` |
| quarantine_outcome/GCST007511/Kunkle_etal_Stage1_results.txt | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007511/Kunkle_etal_Stage1_results.txt | 569380789 | `ae280a70eb82d50c2f90853f1cce592a43b2fa6dc497c5f3cc5381be9c1872ee` |
| quarantine_outcome/finngen_R10_G6_ALZHEIMER/finngen_R10_G6_ALZHEIMER.gz | https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_G6_ALZHEIMER.gz | 804900131 | `56b53bf8376e88f05ffb7c7fcb637c4adb5dcfeb3a4c55302670aba1be36d6ec` |

Downloaded 2026-10-15 UTC (per-file times in MANIFEST.tsv). The MAGMA zip is
unpacked in place to `ld_1kg_eur/g1000_eur.{bed,bim,fam,synonyms}`.
