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

## Module 2 — base-rate model (osf.io/xb7fy)

Run completed 2026-10-22 21:38 UTC ("ALL STAGES COMPLETE", log
results/logs/real/baserate_run8.log). All results rest on extracted values
marked `verified: false`; they are provisional.

Notes so far (not departures):

1. *The step-size rule binds on the real data.* §4 requires the derived error
   bound h²Φ(θ̂) ≤ 10⁻³ at every fitted estimate, halving the step until it
   holds. For variant A on the real counts:

   | step | bound | Radau max \|Δp\| |
   |---|---|---|
   | 0.25 y | 2.3×10⁻² | 7.2×10⁻⁶ |
   | 0.125 y | 5.7×10⁻³ | 1.8×10⁻⁶ |
   | 0.0625 y | 1.4×10⁻³ | 4.5×10⁻⁷ |

   So the fit needs h = 0.03125 y. The independent Radau audit passes easily
   at every step: the bound is conservative. The rule is applied as written,
   and nothing about it was changed.
2. *Run time.* A likelihood-and-gradient evaluation costs 34 ms at h = 0.25
   and 271 ms at h = 0.03125. Cross-validation folds use the step the primary
   fit settled on (§4–5). The whole registered run is therefore expected to
   take on the order of two days of compute. The first attempt was stopped by
   the tool's 2-hour limit for background jobs while the first stage was
   still running. That stage held both the maximum-likelihood fit and NUTS,
   so no checkpoint had been written. The run was restarted unchanged, as a
   detached process. Stage checkpoints protect every finished stage.
3. *Container restart (2026-10-15 → 2026-10-21).* The second attempt ran as
   a detached process. It was lost when the cloud container was reclaimed:
   its last log line is 15:57 UTC on 15 October, while it was fitting
   variant A at h = 0.03125 y. Again, no stage had finished. Checkpoints were
   then added inside stages (each step-size attempt, each NUTS chain, each
   cross-validation fold; commit 498a664). Tests show this changes no
   computation: the per-chain NUTS draws equal those of the multi-chain run.
   Module 2 was re-validated on the new code before the run restarted.
   Not a departure: no setting, rule or input changed.
4. *Second container restart (2026-10-21, between 08:26 and 10:36 UTC).* The
   run was lost again while fitting variant A at h = 0.03125 y. The three
   finished step-size attempts (0.25, 0.125 and 0.0625 y) were saved as
   sub-checkpoints and reused on resume, so only the in-progress fit was
   lost. `scripts/ensure_baserate_running.sh` now restarts the run unchanged
   whenever it has died, and an hourly scheduled check calls it. Not a
   departure.
5. *Third container restart (2026-10-21, after 12:13 UTC; resumed 15:38).*
   Variant A's MLE at h = 0.03125 y passed the step rule (bound 3.6×10⁻⁴,
   Radau max |Δp| 1.1×10⁻⁷) and was saved at 12:13. The container was then
   reclaimed while the first NUTS chain was sampling, during a gap in which
   this session was paused by the account's usage limit. The resume reused
   all four saved MLE attempts. Not a departure.
6. *Robustness rule could not be applied in full (implementation gap, reported,
   not fixed after the fact).* §8 says a §2 conclusion is robust only if it
   holds in sensitivity analyses 1–5. The registered code re-evaluates the
   falsification test under analyses 1–5, but not criterion 2: the
   cross-validation is run only for the primary design and for analysis 7.
   So it cannot be shown that the criterion-2 verdict holds in 1–5, and the
   conclusion is therefore **not** called robust. Re-running the
   cross-validation under analyses 3–5 would be an addition made after
   seeing the result. It is left to the researcher; if it is done, it must
   be labelled as such.
7. *Whole-band cross-validation (analysis 7): unverified fold fits.* One of 5
   fold fits had no start passing the convergence test, for B and for C
   (A: 0 of 5). As §5 prescribes, each used its best start and is reported as
   unverified. Analysis 7 is shown for transparency only.

**Departures: none.** Every setting, rule and input is as registered.
