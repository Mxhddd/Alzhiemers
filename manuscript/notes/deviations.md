# Departures from the registered plans

Registered plans: Module 2 base-rate model, https://osf.io/xb7fy/ (frozen
2026-10-15, `prereg/baserate.freeze.yaml`); Module 5 meta-analysis,
https://osf.io/9kb4w/ (frozen 2026-10-15, `prereg/metaanalysis.freeze.yaml`).
Every result below is built on extracted values marked `verified: false`. It is
provisional until the researcher has checked the inputs.

A *departure* is anything run differently from the registered text. An
*implementation note* records a place where the registered text leaves a
detail open, or where the code's exact behaviour is worth stating. Notes are
listed so nothing is hidden, but they are not departures.

## Procedure (both modules)

* **Not a departure.** The registered OSF documents are PDFs with no text
  layer, so they cannot be compared with the repository Markdown byte for
  byte. The PDFs are kept byte-identical in `prereg/registered/` (sha256 equal
  to OSF's own hash, checked by `tests/unit/test_prereg_registered.py`).
  `prereg/baserate.md` and `prereg/metaanalysis.md` were replaced with
  faithful transcriptions of the registered text, and those transcriptions
  are what the freeze hashes. Every registered number and rule matched
  `config/*.yaml` before freezing; no registered text and no setting was
  changed.

## Module 5 — meta-analysis (osf.io/9kb4w)

**Departures: none.** The run (`make metaanalysis`, 2026-10-15) repeated the
positive control first: 1.62 (1.26–2.08), k = 7, recovered. It then applied
§2, §7 and §8 as written.

Implementation notes:

1. *Single-study intervals.* §4 says a single study "is shown with its own
   transcribed interval". The forest plot and `meta_pooled_by_class.csv`
   show the transcribed bounds. The verdict table instead uses estimate ±
   1.96·SE, with SE backed out of the transcribed interval (the same SE that
   pooling uses). The only predicted contrast that rests on a single study is
   DHP (Colbourne 2022). There, both intervals are 0.82 (0.72–0.93) to two
   decimals and give the same verdict. Checked after the run; nothing
   changes.
2. *Orientation of the exploratory contrast.* The config key is named
   `z_drug_vs_dual_orexin_antagonist`, but the estimate is oriented as the
   registered table says: dual orexin antagonist versus Z-drug. The paper's
   Z-drug-versus-DORA ratio is inverted through `invert_ratio`. No verdict is
   attached (registered as exploratory).
3. *DerSimonian–Laird sensitivity.* Only the τ² estimator changes; the
   Knapp–Hartung interval is kept, as §8 item 1 says ("DerSimonian–Laird
   instead of REML").
