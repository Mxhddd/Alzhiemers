# Pre-registration — Module 1: drug-target Mendelian randomization

**Status: COMPLETE DRAFT, ready for the researcher to review, upload to OSF
and freeze — not yet frozen.** Written by Claude under the standing policy.
Claude must not run `make freeze-prereg`. Configuration:
`config/mr_targets.yaml` and `config/mr_methods.yaml` at the commit that
freezes this file. No AD outcome data have been read. The module refuses
to run until this file is frozen and its validation passes.

## Change log

| Date | Change |
|---|---|
| 2026-10-13 | Complete draft. Replaces `mr.DRAFT.md`. Decisions filled in: cis window, positive control per target (and three targets declared untestable), multiple testing, replication and proxy-free sensitivity, TwoSampleMR cross-check. |
| 2026-10-14 | **Cross-check (step 4) made operational.** "Disagree" now has a fixed definition (`config/mr_methods.yaml` `crosscheck`). b, SE, p, Q and the Egger intercept must match TwoSampleMR within floating-point noise (absolute 10⁻¹⁰ or relative 10⁻⁸). The weighted-median SE and p must match within 15% relative, because TwoSampleMR draws its bootstrap samples in a different order. Every comparison is written to `results/tables/mr_crosscheck_twosamplemr.csv`, whether it agrees or not. Before any real data: on simulated data, the base-R ports agree with TwoSampleMR 0.7.12 for all four estimators (analytic quantities to < 10⁻¹²). Harmonisation of 400 simulated variants keeps the same variants with identical aligned effects (`tests/unit/test_mr_crosscheck.R`). | TwoSampleMR became installable (Posit Package Manager binaries, r-universe). A cross-check needs a stated tolerance, or "however small" would make every bootstrap SE a disagreement. |

## 1. Hypothesis under test

The locus coeruleus sets the *rate* of Alzheimer's progression through
its regulatory function. Genetic proxies for lower activity of
noradrenergic drug targets should therefore shift AD risk in the
directions below.

## 2. Results that would count against the hypothesis (stated first)

- A pre-specified target whose primary estimate is **discordant** counts
  against the hypothesis as stated. Discordant means opposite to the
  prediction, with a 95% CI that excludes zero.
- A uniformly null set of concordant-direction estimates is reported as
  "no evidence", not as support.
- Each target is its own falsifiable prediction.
  - It is read at 95%.
  - Whether the result survives Bonferroni correction over the targets
    actually tested (those that pass their positive control) is reported
    alongside.
  - A conclusion that changes between the two levels is called fragile.

## 3. Targets, predictions and positive controls

Expression instruments are oriented to the expression-*increasing* allele.
"Lower activity protective" therefore means a **positive** estimate per
higher expression.

| Gene | Drug analogue | Lower activity predicted | Expected sign per higher expression | Positive control (sign per higher expression) |
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

Evidence and reasoning for each control are in
`manuscript/notes/mr_positive_controls.md` (GWAS Catalog, queried
2026-10-13).

**The control is a hard gate.** A target passes if its instruments move
the control trait with p < 0.05 and the expected sign. The gate is
evaluated before the AD outcome is read.

- A target that fails its control is excluded and logged.
- A target with no defensible control (LGMN, SLC18A2, MAOA) is reported
  as **untestable with current instruments**, never as a null result.

## 4. Data

- **Outcome:** Bellenguez et al. 2022, GWAS Catalog GCST90027158 (stage 1,
  GRCh37). This includes UK Biobank proxy cases.
- **Proxy-free sensitivity:** Kunkle et al. 2019 stage 1 (clinically
  diagnosed AD, GCST007511), where the instruments are available in it.
- **Replication:** FinnGen public release current at run time, Alzheimer's
  disease endpoint. The release and endpoint code are logged when
  downloaded; if the code differs from G6_ALZHEIMER, the substitution is
  logged and justified.
- **Exposures:**
  - cis-eQTLs: eQTLGen (blood) and MetaBrain (cortex);
  - cis-pQTLs: UKB-PPP / deCODE, where available, lifted over to GRCh37.
  - The build is asserted everywhere; a mismatch stops the run.
- **Control traits:** published GWAS summary statistics on GRCh37 for
  resting heart rate, diastolic and systolic blood pressure, fasting
  glucose and sleep duration. The specific files are logged at download.

## 5. Analysis plan (all settings in the config files)

1. **Instruments:** the cis window is **±100 kb** of the gene body (Ensembl
   GRCh37 coordinates via biomaRt; release logged). Filters: p < 5×10⁻⁸,
   MAF ≥ 0.01, F > 10, LD clumping at r² < 0.001 within 10 Mb (1000
   Genomes EUR).
   - Why ±100 kb: most cis-QTL effects lie within 100 kb, and a narrow
     window limits pleiotropy from neighbouring genes.
   - **Sensitivity:** ±1,000 kb (eQTLGen's own cis definition).
2. **Positive control** (§3), before the outcome is read.
3. **Primary estimator:**
   - IVW with multiplicative random effects when there are ≥ 2 instruments;
   - Wald ratio with 1 instrument.
   - Secondary: MR-Egger and weighted median.
4. **Cross-check:** every estimate is re-run with `TwoSampleMR` and compared
   with the base-R implementation. Every discrepancy, however small, is
   reported. If they disagree, the `TwoSampleMR` result is the one used.
5. **Colocalisation:** ABF coloc over the cis window, with p1 = p2 = 10⁻⁴
   and p12 = 5×10⁻⁶. Supportive if PP4 ≥ 0.8. SuSiE-coloc for
   multi-signal regions if susieR is available; otherwise it is reported
   as not run.
6. **Steiger directionality:** supportive evidence only.
7. **Concordance labels:**
   - concordant;
   - direction as predicted, CI includes null;
   - opposite, CI includes null;
   - DISCORDANT.
8. **Multiplicity:** §2.

## 6. Validation that must pass first (`make validate-mr`)

- **Reproduction** of a published two-sample MR of LDL-C (HMGCR variants)
  on coronary artery disease. The estimate must fall within the published
  CI.
  - The published estimate and the variant-level inputs are to be
    extracted from the primary paper with locators and `verified: false`.
  - That extraction is Module 1 work, not done in this round.
- **Gene-desert null test:** the IVW false-positive rate is not above
  nominal.
- If the reproduction misses the published CI, that is reported and
  nothing is adjusted to make it match.

## 7. Limitations, stated in advance

- Expression in blood or cortex may not reflect target activity in the LC.
- Few instruments per target means low power and wide CIs.
- The outcome GWAS includes proxy cases.
- Positive controls test the instrument's biology, not the AD pathway.
- Three targets are untestable.
