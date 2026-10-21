# Noradrenergic drug targets and Alzheimer's disease: a pre-registered drug-target Mendelian randomisation

Mohammad Hassan, independent researcher · October 2026

This is my plan for testing, with genetic proxies, whether lowering the activity of noradrenergic drug targets changes the risk of Alzheimer's disease (Module 1 of my project). I am registering it before any row of any Alzheimer's outcome file has been read, so that which targets are tested, how, and what counts as failure cannot move once the outcome is seen. Some work on other traits has already been done; section 0 says exactly what, and the change log below says when each part of the plan was set. The code and settings are in the project repository (`config/mr_targets.yaml` and `config/mr_methods.yaml`), at the commit where I register this file. The analysis code refuses to read the outcome until this document is frozen and the validation in section 6 has passed.

## Tools

I used Claude, Anthropic's AI model, through Claude Code to write the code, extract the published values and draft this text. Claude also downloaded the public data, checked the access terms of each data source, and ran the positive controls and the HMGCR validation described in section 0. Where I left a methodological setting to it, it chose one and wrote its reason next to the setting in the config files. I made the main decisions of 2026-10-21 myself (change log). By registering this plan I adopt all of these choices as mine.

## Change log

None of these changes was made after seeing any Alzheimer's outcome data. No outcome row has been read.

| Date | Change | Reason |
|---|---|---|
| 2026-10-13 | First complete version: targets, predictions, a positive control per target (three targets declared untestable), cis window, estimators, multiplicity, replication and proxy-free sensitivity, TwoSampleMR cross-check. | — |
| 2026-10-14 | **Cross-check made operational.** "Disagree" now has a fixed definition: b, SE, p, Q and the Egger intercept must match TwoSampleMR within absolute 10⁻¹⁰ or relative 10⁻⁸; the weighted-median SE and p within 15% relative (`config/mr_methods.yaml`, `crosscheck`). | A cross-check needs a stated tolerance. TwoSampleMR draws its bootstrap samples in a different order, so without a tolerance every bootstrap SE would be a "disagreement". |
| 2026-10-15 | Data step: public inputs downloaded; the Alzheimer's outcome files downloaded and quarantined (header lines only); positive controls run with blood eQTLs (eQTLGen). The plan was not changed. | No target passed its positive control with blood eQTLs (section 0), and four targets have no blood eQTL at all. |
| 2026-10-21 | **Brain eQTLs and plasma pQTLs added as instrument sources, alongside eQTLGen** (MetaBrain cortex; UKB-PPP and deCODE plasma), with a fixed rule for choosing a target's primary source (section 5.2). The positive controls, their traits and the gate are unchanged. The blood-only result is kept and reported. | Blood expression may not proxy target activity, and four targets were not in eQTLGen. This is my decision (option b), taken after seeing the blood-only positive controls and before any outcome data. Section 0 says what this means for interpretation. |
| 2026-10-21 | **HMGCR validation filled in** (section 6.1): Huang et al. 2021, eLife, with the comparison rule fixed before the run. | The validation target had no published estimate. I asked for one whose instruments and outcome are open data. |
| 2026-10-21 | **Gene-desert null test defined** (section 6.2): Ovcharenko et al. 2005, applied to GENCODE release 19 on GRCh37. The test runs after the freeze, before any target's MR. | The test reads the Alzheimer's outcome by design, so it cannot run before the freeze. |
| 2026-10-21 | **Replication release fixed to FinnGen R13**, endpoint G6_ALZHEIMER (16,652 cases, 483,534 controls). | The plan said "the release current at run time"; that is R13, and its files can be downloaded without an account. |

## 0. What has already been done, before registration

I am stating this first, because it bears on how much the plan could have been shaped by results.

**Downloaded.** All exposure, positive-control and validation inputs, and three Alzheimer's outcome files (Bellenguez 2022 stage 1, Kunkle 2019 stage 1, FinnGen R13 G6_ALZHEIMER). The outcome files are kept apart (`data/raw/mr/quarantine_outcome/`). For each, only the header line was read, to write the column map; their checksums are recorded. No MR, liftover, filter or summary has been run on them. A FinnGen R10 copy, downloaded on 2026-10-15, was deleted unread on 2026-10-21.

**Not read.** No row of any Alzheimer's or dementia outcome file.

**Run on other traits.**

* Positive controls with blood eQTLs (2026-10-15), and with the new sources (2026-10-21, after the source-selection rule in section 5.2 was written into the config). Results are below.
* The HMGCR validation (section 6.1), after its target estimate and comparison rule were written down.

What I knew when I made the 2026-10-21 decision: the blood-only positive controls had all failed. Adding sources can only add targets to the tested set; it cannot change any target's Alzheimer's result, which nobody has seen. The positive-control traits, the expected signs and the gate are exactly as registered on 2026-10-13.

### 0.1 Positive controls, blood eQTLs only (option a)

{{PC_BLOOD}}

### 0.2 Positive controls with every available source (option b)

{{PC_ALL}}

### 0.3 HMGCR validation

{{HMGCR_RESULT}}

## 1. Question and hypotheses

The locus coeruleus sets the *rate* of Alzheimer's progression through its regulatory function. If so, genetic proxies for lower activity of noradrenergic drug targets should shift Alzheimer's risk in the directions in section 3.

Each target is its own prediction. For a target with a predicted direction, the hypothesis is that its primary estimate lies on the predicted side of zero.

## 2. What would count against the hypothesis

I am stating these before the data section, on purpose.

* A target whose primary estimate is **discordant** counts against the hypothesis as stated. Discordant means opposite to the prediction, with a 95% confidence interval that excludes zero.
* A set of estimates that are all null, even if they point the predicted way, is reported as "no evidence", not as support.
* A target is read at the 95% level. I also report whether each result survives Bonferroni correction over the targets actually tested (those that pass their positive control). A conclusion that changes between the two levels is called fragile.
* A target that cannot be tested (no defensible positive control, no instrument, or a failed positive control with every source) is reported as **untestable with current instruments**. It counts neither for nor against.

## 3. Targets, predictions and positive controls

Instruments are oriented to the allele that increases the measured level (expression for eQTLs, protein for pQTLs). "Lower activity protective" therefore means a **positive** estimate per higher level.

| Gene | Drug analogue | Lower activity predicted | Expected sign per higher level | Positive control (sign per higher level) |
|---|---|---|---|---|
| SLC6A2 | atomoxetine | protective | + | resting heart rate (−) |
| ADRB1 | lipophilic β-blockers | protective | + | resting heart rate (+) |
| ADRA1A | prazosin, doxazosin | protective | + | diastolic blood pressure (+) |
| ADRA2A | guanfacine (agonist) | **harmful** | − | fasting glucose (+) |
| LGMN | — | protective | + | **none defensible: untestable** |
| SLC18A2 | — | **harmful** | − | **none defensible: untestable** |
| MAOA | MAO-A inhibitors | protective | + | **none defensible (X-linked): untestable** |
| DBH | disulfiram, nepicastat | exploratory | none | systolic blood pressure (+) |
| CACNA1D | isradipine | protective | + | systolic blood pressure (+) |
| HCRTR2 | suvorexant, lemborexant | exploratory | none | sleep duration (−) |

The evidence for each control is in `manuscript/notes/mr_positive_controls.md`.

**The positive control is a hard gate.** A target's instruments from a given source pass if they move the control trait with p < 0.05 and the expected sign. The gate is evaluated before the outcome is read. A target whose instruments fail the gate with every available source is excluded and reported as untestable with current instruments, never as a null result.

## 4. Data

**Outcome.** Bellenguez et al. 2022, stage 1, GWAS Catalog GCST90027158: 85,934 cases (39,106 clinically diagnosed and 46,828 UK Biobank proxy cases) and 401,577 controls. The GWAS Catalog file is on GRCh38; it is lifted to GRCh37 after the freeze (section 5.1).

**Proxy-free sensitivity.** Kunkle et al. 2019 stage 1 (GCST007511), clinically diagnosed cases only. Its file has no allele-frequency column; after the freeze it is given one from 1000 Genomes EUR (non-palindromic variants only).

**Replication.** FinnGen release R13, endpoint G6_ALZHEIMER: 16,652 cases and 483,534 controls (official R13 manifest). GRCh38, lifted to GRCh37 after the freeze.

**Exposures.** Four sources. I checked the access terms of each on 2026-10-21; only what can be downloaded without an account, form or agreement was downloaded.

| Source | Measure | What the full data need | What is used here |
|---|---|---|---|
| eQTLGen (Võsa et al. 2021) | blood cis-eQTLs, 31,684 people | open download | full cis-eQTL file |
| MetaBrain (de Klein et al. 2023), Cortex-EUR | brain cis-eQTLs, 2,683 people | a form with name, e-mail and institute | the open table behind the MetaBrain browser: the top variant of each gene (conditional rank 1), GRCh38 lifted to GRCh37 |
| UKB-PPP (Sun et al. 2023) | plasma cis-pQTLs (Olink), 34,557 people of European ancestry | a Synapse account | Supplementary Table 9 of the paper (CC BY): the primary pQTL of each locus, GRCh38 lifted to GRCh37 |
| deCODE (Ferkingstad et al. 2021) | plasma pQTLs (SomaScan), 35,559 Icelanders | a form with name, e-mail, institution and terms | not used: the open supplement has cis pQTLs among the targets only for LGMN, which is untestable |

The open MetaBrain and UKB-PPP tables contain only lead variants. With them each target gets at most one instrument per source, and colocalisation over the cis window cannot be run for these two sources.

**Positive-control traits** (all GRCh37): resting heart rate (Zhu et al. 2019, GCST007609); diastolic and systolic blood pressure (Evangelou et al. 2018, GCST006630 and GCST006624); fasting glucose (Chen et al. 2021, GCST90002232); sleep duration (Dashti et al. 2019, GCST007561).

**LD reference.** 1000 Genomes Phase 3 EUR (503 people, GRCh37).

## 5. Methods (all settings in the config files)

### 5.1 Instruments

* Cis window **±100 kb** of the gene body (Ensembl GRCh37 coordinates; release logged). Sensitivity: **±1,000 kb**.
* Filters: p < 5×10⁻⁸, minor-allele frequency ≥ 0.01, F > 10, LD clumping at r² < 0.001 within 10 Mb.
* Genome build GRCh37 everywhere. GRCh38 files are lifted with the UCSC hg38ToHg19 chain; a position that does not map, maps to several places, maps to the other strand or to another chromosome is dropped and logged. A build mismatch stops the run.
* eQTLGen gives Z-scores only; beta and SE are derived from Z, allele frequency and N (Zhu et al. 2016).

### 5.2 Choosing a target's instrument source

This rule was written into the config on 2026-10-21, before any positive control with the new sources was run.

1. A source is eligible for a target only if its instruments pass that target's positive control in the same window.
2. Among eligible sources, the first in this order is the **primary source**: **UKB-PPP, deCODE, MetaBrain, eQTLGen**. The others are reported as secondary analyses. Sources are never pooled.
3. Reason for the order: a drug acts on a protein, so protein level is one step closer to drug action than mRNA; among transcript sources, cortex is the tissue of the disease and of the locus coeruleus's projections, while blood is a convenience tissue. UKB-PPP comes before deCODE because its participants are of European ancestry like the outcome GWAS, its open table has full-precision effects, and antibody assays are less exposed to binding artefacts than aptamer assays.
4. The **blood-only** result (eQTLGen alone, the design as it stood on 2026-10-15) is always reported alongside, labelled "option (a), blood only".

**Orientation.** eQTL instruments are coded to the expression-increasing allele and pQTL instruments to the protein-increasing allele. The Wald ratio and IVW do not depend on which allele is coded, so this only fixes the reporting: every estimate is per 1 SD higher genetically predicted level, and the signs in section 3 apply unchanged to protein instruments.

### 5.3 Analysis of each target

1. **Positive control** (section 3), before the outcome is read.
2. **Primary estimator:** IVW with multiplicative random effects for 2 or more instruments; Wald ratio for 1. Secondary: MR-Egger and weighted median (each needs at least 3 instruments).
3. **Cross-check:** every estimate is re-run with TwoSampleMR and compared at the tolerances in the change log (2026-10-14). If they disagree, the TwoSampleMR result is used and the disagreement is reported.
4. **Colocalisation:** ABF coloc over the cis window, with p1 = p2 = 10⁻⁴ and p12 = 5×10⁻⁶; supportive if PP4 ≥ 0.8. Only possible for eQTLGen (sections 4 and 9).
5. **Steiger directionality:** supportive evidence only.
6. **Concordance label** of each primary estimate: concordant; direction as predicted, CI includes null; opposite, CI includes null; DISCORDANT.

## 6. Checks that must pass first (`make validate-mr`)

### 6.1 Reproducing a published HMGCR estimate

**Target, fixed before the run.** Huang et al. 2021 (eLife 10:e73873), Supplementary file 1, Table 6: genetically proxied HMGCR (LDL-C) on coronary heart disease, IVW, OR **1.44** (95% CI **1.24–1.68**) per unit of the LDL-C GWAS effect (1 SD). The 7 instruments are listed in their Table 3 (variants within ±100 kb of HMGCR, p < 5×10⁻⁸, r² < 0.30). Their exposure GWAS is GLGC 2013 and their outcome CARDIoGRAMplusC4D 1000 Genomes; both are open, and the 7 SNPs' LDL-C statistics in our GLGC file equal their Table 3 to every printed digit.

**Rule.** Our IVW estimate for the same 7 SNPs must fall inside the published 95% CI on the log scale. Whether the published estimate falls inside our CI is reported but does not decide. If it fails, that is reported and nothing is adjusted.

Ference et al. 2016 (NEJM) was considered. Its website refused automated access, so its confidence intervals were not read. It also used individual-level data and its own composite outcome, so it cannot be reproduced from open data. The candidates are listed in `data/extracted/mr_sources.md`.

### 6.2 Gene-desert null test

**Definition.** Ovcharenko et al. 2005 (Genome Res 15:137) define gene deserts as "the 3% longest intergenic intervals" of their gene annotation; in their annotation the shortest was 640 kb. I apply the same rule:

* genes: GENCODE release 19 (GRCh37) protein-coding gene records; overlapping genes merged;
* autosomes; segments bounded by chromosome ends, centromeres, telomeres, heterochromatin and short-arm gaps, and assembly gaps longer than 250 kb (NCBI GRCh37.p13 gap table), as in the paper's Methods; an interval touching a segment boundary is excluded;
* deserts: the 3% longest of the remaining intervals.

{{DESERT_RESULT}}

**Variants.** From the 1000 Genomes EUR panel: bi-allelic SNVs with an rsID and minor-allele frequency ≥ 0.01, at least 100 kb inside a desert from both flanking genes, thinned to at least 50,000 bp apart, then LD-pruned to r² < 0.001 within 10 Mb.

**Test** (fixed in `tests/validation/targets/mr_null_gene_desert.yaml`): 1,000 random sets of 5 desert variants present in the outcome, each given a synthetic exposure, analysed with IVW. It passes if the share of sets with p < 0.05 does not exceed the 99% binomial quantile under the null.

**When.** This test reads the Alzheimer's outcome, so it runs after the freeze and before any target's MR. If it fails, no target's MR is run and the failure is reported.

## 7. Inference criteria

* Each target with a prediction is read at the 95% level, with Bonferroni correction over the targets that pass their positive control reported alongside (section 2).
* DBH and HCRTR2 are exploratory: their estimates are reported without a verdict.
* The primary result of a target comes from its primary source (section 5.2) in the ±100 kb window.

## 8. Sensitivity and secondary analyses

1. The ±1,000 kb window, with the source rule re-applied in that window.
2. Option (a), blood only: eQTLGen instruments for every target, whatever the hierarchy says.
3. Every eligible secondary source, reported per target.
4. Proxy-free outcome (Kunkle 2019).
5. Replication in FinnGen R13.
6. MR-Egger and weighted median where there are at least 3 instruments.

A verdict is called robust only if the primary estimate's concordance label is the same under analyses 1, 4 and 5.

## 9. Limitations, stated in advance

* Expression in blood or cortex, or protein in plasma, may not reflect target activity in the locus coeruleus.
* The open MetaBrain and UKB-PPP tables give one instrument per target, so those estimates are Wald ratios, colocalisation cannot be run, and pleiotropy cannot be examined with multi-instrument methods. The full files need a form or an account that I have not yet completed.
* UKB-PPP, the outcome (UK Biobank proxy cases) and several positive-control GWAS share UK Biobank participants. With instruments this strong (F far above 10) the bias from sample overlap should be small, but it is not zero.
* The main outcome includes proxy cases.
* Positive controls test that the instruments capture the target's biology, not the Alzheimer's pathway.
* LGMN, SLC18A2 and MAOA cannot be tested.
* All published values were extracted with Claude's help and are unverified until I check them against the sources.
