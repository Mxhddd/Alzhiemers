# Pre-registration DRAFT — Module 1: drug-target Mendelian randomization

> **Status: DRAFT written by Claude from PROJECT_BRIEF.md for the researcher to
> review.** It is not frozen. To freeze: edit as needed, save as `prereg/mr.md`,
> upload to OSF, then `make freeze-prereg NAME=mr OSF=https://osf.io/<id>/` and
> commit both files. `make mr` refuses to run until that is done.

## Hypothesis under test

The locus coeruleus sets the *rate* of Alzheimer's progression through its
regulatory function. Genetic proxies for lower activity of noradrenergic drug
targets should therefore shift Alzheimer's risk in the directions below.

## Results that would count against the hypothesis

Any pre-specified target whose primary estimate is **discordant** (opposite to
the prediction with a CI that excludes zero) counts against the hypothesis as
stated. A uniformly null set of concordant-direction estimates is reported as
"no evidence", not as support.

## Targets and predictions (copied from the brief)

| Gene | Drug analogue | Predicted effect of lower activity | Expected sign of estimate per higher expression |
|---|---|---|---|
| SLC6A2 | Atomoxetine | Protective | + |
| ADRB1 | Lipophilic β-blockers | Protective | + |
| ADRA1A | Prazosin, doxazosin | Protective | + |
| ADRA2A | Guanfacine (agonist) | Harmful | − |
| LGMN | — | Protective | + |
| SLC18A2 | — | Harmful | − |
| MAOA | MAO-A inhibitors | Protective | + |
| DBH | Disulfiram, nepicastat | Exploratory | none |
| CACNA1D | Isradipine | Protective | + |
| HCRTR2 | Suvorexant, lemborexant | Exploratory | none |

## Data

- Outcome: Bellenguez et al. 2022, GWAS Catalog GCST90027158 (stage 1, GRCh37,
  includes UK Biobank proxy cases). Sensitivity: proxy-free statistics where
  available. Replication: FinnGen public AD endpoint [researcher to name].
- Exposures: cis-eQTLs (eQTLGen blood; MetaBrain cortex); cis-pQTLs (UKB-PPP /
  deCODE) where available.

## Analysis plan (all settings in config/mr_targets.yaml and config/mr_methods.yaml)

1. Instruments: cis window ±100 kb of the gene body (Ensembl GRCh37), p < 5e-8,
   MAF ≥ 0.01, F > 10, LD clumping r² < 0.001 within 10 Mb.
2. **Positive control (hard gate), evaluated before the AD outcome is read.**
   ADRB1 instruments must move resting heart rate in the expected direction
   (p < 0.05). Targets without a positive control are excluded. [Researcher:
   specify controls for the other targets or accept their exclusion.]
3. Primary estimate: IVW (multiplicative random effects) with ≥ 2 instruments,
   Wald ratio with 1. Secondary: MR-Egger, weighted median.
4. Colocalisation: ABF coloc (p1 = p2 = 1e-4, p12 = 5e-6); supported if PP4 ≥ 0.8.
   SuSiE-coloc for multi-signal regions [currently not runnable: susieR].
5. Steiger directionality reported as supportive evidence.
6. Concordance labels: concordant / direction as predicted, CI includes null /
   opposite, CI includes null / DISCORDANT.
7. Multiple testing: [Researcher to decide — e.g. Bonferroni over the 8
   pre-specified targets, or report unadjusted with the count of tests.]

## Validation that must pass first

- Reproduction of a published LDL-C (HMGCR) → CAD two-sample MR estimate
  [paper to be named by the researcher].
- Gene-desert null test: IVW false-positive rate not above nominal.
