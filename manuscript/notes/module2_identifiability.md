# Module 2 — what stage-by-age counts can and cannot test

> **Status: open issue raised by Claude on 2026-10-13, awaiting the
> researcher's decision.** All numbers below come from SIMULATED data (no
> real data have been fitted). `prereg/baserate.md` is marked DO NOT FREEZE
> until this is resolved.

## 1. The problem: criterion 1 reads ageing as amyloid gating

Variant A (as specified in PROJECT_BRIEF.md §5) gives the LC → cortex
transition (Braak stage c → 1a) the rate

    h(t) = μ0 (1 + (HR_A − 1) A(t)),

where A(t) is the *population* amyloid prevalence at age t. Two facts combine:

1. The data are stage counts by age band. They carry no per-brain amyloid
   status, so A(t) enters only as a known function of age.
2. A(t) rises with age, and h(t) has no other way to change with age.

So **any** increase of the c → 1a rate with age is fitted as HR_A > 1,
whether its cause is amyloid, age itself, or anything else that rises with
age. Criterion 1 ("supports gating" if the 95% interval for log HR_A is
above 0) and the γ consistency check therefore test "does this transition
speed up with age?", not "does amyloid speed it up?". Nearly every
age-related neuropathological transition speeds up with age, so in practice
the criterion could essentially only come out in favour of the hypothesis.
This is the situation the standing guardrail names ("an analysis that
cannot come out against the hypothesis → stop and tell").

## 2. Simulated demonstration (`python -m sim.baserate_age_confound 0.0 0.1 0.3`)

Truth: the dry run's synthetic world (12 Braak-2011 stages, 10 decade bands,
230 brains each), but with **no amyloid effect**: the c → 1a rate is
exp(α + β (t − 50)/10), with its age-50 value equal to the A world's.
Variant A is fitted to the *expected* counts (230 × true probabilities),
which gives the large-sample value it converges to, free of Monte Carlo
noise. SEs are at n = 2,300.

| true age slope β of c → 1a (per decade), no amyloid effect | fitted HR_A | log HR_A, 95% CI | criterion 1 verdict | γ check vs published 2.43 |
|---|---|---|---|---|
| 0 | 1.00 | −0.76 to 0.76 | inconclusive | — |
| 0.1 (same slope as the neighbouring edges in the simulation) | **2.22** | **0.28 to 1.31** | **supports gating** | **consistent** |
| 0.3 | **7.58** | **1.61 to 2.44** | **supports gating** | disagrees (too large) |

A modest, amyloid-free age slope is enough for the pre-registered analysis
to report amyloid gating at almost exactly the published hazard ratio.

Even with no age slope at all (first row), the 95% interval is ±0.76, so the
"against" outcome (upper end < 0.322) needs the estimate to fall below
−0.44 by chance. Under the Wald approximation with SE 0.388, that happens
about 13% of the time.

## 3. Can cross-validation (criterion 2) discriminate? (`python -m sim.baserate_discrimination`)

> Superseded by §6 (the power study, 2026-10-14). Kept for the record.

The truth is each variant in turn. The table gives the in-sample
log-likelihood gap between the true variant and each wrong one, fitted to
expected counts. That gap is the *average* evidence a study of 2,300 brains
in 10 decade bands would have (n × KL divergence).

| truth | gap to A_gated | gap to B_independent | gap to C_multihit |
|---|---|---|---|
| A_gated (HR 2; the dry-run world) | — | 8.27 | 35.72 |
| B_independent | 116.19 | — | 17.06 |
| C_multihit (2 hits) | 32.90 | 6.08 | — |

In the full-settings dry run (one simulated data set from the A world), the
in-sample gap A vs B was 9.1 (log-likelihoods −166.55 vs −175.63). But
5-fold cross-validation by age band put **B ahead of A by 0.43**: "not
distinguishable" under the margin of 2. Held-out bands are dealt
round-robin, so each fold must interpolate two whole decades. That loses
most of A's in-sample advantage.

Reading: criterion 2 **can** come out against A. If the truth is B or C, A
loses by about 116 or 33 units in-sample, far beyond the margin of 2. If the
truth is A, it struggles to separate A from B (8.3 in-sample, −0.43 held-out
in the one dry-run data set). B vs C is also close when C is true (6.1). So
the criterion is capable of "against" but weak for "supports". That
asymmetry is acceptable under the guardrail but should be stated in the
pre-registration. These are in-sample gaps on expected counts; held-out gaps
are smaller, and a proper CV power study (many simulated data sets per truth)
has not been run (about 30 min per data set at full settings).

## 4. Other dry-run observations (simulated; not caused by the issue above)

* The falsification test passed for all three variants in the primary
  scenario. That is expected: the synthetic incidence envelope was built
  around the true stage-V incidence. It therefore does not discriminate
  between variants here, and was never meant to.
* The threshold-IV sensitivity scenario failed falsification for all three
  variants. This is an artefact of the same construction: the envelope is
  built around stage-V incidence, and stage-IV incidence is higher.
* Band-average evaluation (sensitivity) failed falsification for A. The
  synthetic counts were generated at band midpoints, so band-averaging is a
  misspecified likelihood in this world. The run stopped (2-hour background
  limit) before B and C were refitted, so whether they fail too is unknown.

## 5. Options (decided 2026-10-13: option (c))

**(a) Age-adjusted gating.** h(t) = exp(α + β (t − 50)/10) · (1 + (HR_A − 1) A(t)).
HR_A is then identified only by the part of A(t)'s shape that a log-linear
age trend cannot mimic. This removes the bias, but the interval for HR_A
will be wide. Criterion 1 would then almost always be "inconclusive"
in *both* directions: honest, but not a test.

**(b) Per-brain amyloid × tau-stage counts by age.** This is the only design
here that can actually test gating. It uses a bivariate chain (amyloid
status × tau stage) with a shared age slope on c → 1a and the rate
multiplied by HR_A in amyloid-positive brains. HR_A is then identified from
amyloid-positive vs amyloid-negative brains *of the same age*. Possible
sources:

* the Braak 2011 full text, if it cross-tabulates stage by plaque status
  (cannot be checked: the paper is not reachable from here);
* the NACC Neuropathology Data Set, which requires the researcher's data
  application;
* other autopsy series that report stage × amyloid phase by age.

**(c) Drop criterion 1 from Module 2.** Module 2 keeps the base rate, the
variant comparison (criterion 2) and the falsification test. The amyloid-
gating question is left to modules with per-person data (Killer 1:
LC function → tau accumulation in amyloid-positive people).

**Recommendation:** (c) now, with (a) reported as a descriptive estimate
labelled "not a test". Pursue (b) if per-brain data can be obtained. Keep
the brief's literal variant A only as a sensitivity analysis labelled
"confounded with age; not evidence either way". The γ check becomes
descriptive under (a) and (c).

Nothing has been changed in the model, configuration or pre-registration
criteria pending this decision. To keep the gate shut meanwhile, a
validation target was added:
`tests/validation/test_baserate.py::test_age_acceleration_is_not_read_as_amyloid_gating`.
It uses a 3-stage truth with no amyloid effect, an LC → cortex age slope of
0.3 per decade, and 3,000 brains per band. Variant A returns log HR_A = 2.12
(SE 0.08), so `make validate-baserate` now FAILS, and
`make baserate` refuses to run on real data.

## 6. Power study of criterion 2 (2026-10-14, `sim/baserate_power.py`)

**Design.** 100 simulated data sets of the real size: 12 stages, 10 decade
bands, 230 brains each. Twenty per world:

* LC-first true, no amyloid effect, LC → cortex age slope 0, 0.1 or 0.3 per
  decade (the null scenarios required by the revision);
* B true;
* C true (2 hits).

Each data set is analysed as criterion 2 prescribes, with two shortcuts for
speed:

* single warm-started fits;
* step 0.5 years.

**First design (cross-validation with whole age bands held out)** — results
archived in `results/dryrun/power/whole_band_design/`:

| truth | supports | against | inconclusive |
|---|---|---|---|
| LC-first, slope 0 | 17 | 1 | 2 |
| LC-first, slope 0.1 | 17 | 2 | 1 |
| LC-first, slope 0.3 | 19 | 1 | 0 |
| B | 0 | 20 | 0 |
| C | 0 | 18 | 2 |

Two of the four false "against" verdicts had implausibly large gaps (A
behind by 62 and 24 units in an A world). `sim/baserate_power_refit.py`
refitted them with the analysis's own multi-start settings. The diagnosis:

* In the fold that holds out the 91–100 band, stage VI is absent from the
  training data. Every model puts its V → VI rate on the search-box bound
  and predicts the 27 stage-VI brains at about 10⁻¹². That band costs about
  −730 log-likelihood units for every model.
* The verdict then turns on numerical detail. The multi-start refit of one
  replicate scored −856 / −883 / −876 for A / B / C, which reads "supports";
  the power study's single-start fit had read "against".
* In 2 of 4 refitted replicates, one fold fit had no start passing the
  convergence test, and the code stopped. The pre-registration says such
  fits are "reported as such".
* The fourth replicate (slope 0.3, rep 8) was a genuine false "against": B
  fits better than the true A even with multi-start fits (−229.8 vs −238.8).

**Revision.** Criterion 2 now uses cross-validation over brains within age
bands. Whole-band CV is kept as sensitivity analysis 7. Unverifiable fold
fits are flagged, not fatal. See the prereg change log.

**Revised design, same 100 simulated worlds:**

| truth | supports | against | inconclusive |
|---|---|---|---|
| LC-first, slope 0 | 19 | 0 | 1 |
| LC-first, slope 0.1 | 19 | 0 | 1 |
| LC-first, slope 0.3 | 19 | 1 | 0 |
| B | 0 | 20 | 0 |
| C | 0 | 16 | 4 |

**Answer to "can the variant comparison discriminate?"** Yes, in these
simulated worlds, with known error rates.

* LC-first true: "supports" in 57 of 60 (95%, Wilson 95% CI 86–98%);
  "against" in 1 of 60 (1.7%, 0.3–8.9%). The single false "against" is the
  genuine case above.
* B true: "against" in 20 of 20.
* C true: "against" in 16 of 20, inconclusive in 4.
* "Supports" never appeared when B or C was true: 0 of 40, upper 95% bound
  8.8%.

The age slope on the LC → cortex step (0 to 0.3 per decade) does not erode
"supports", because variant A carries its own age slope.

**What Module 2 can and cannot claim.**

* It **can** say which of three model classes best explains the stage-by-age
  distribution, with the error rates above, plus the incidence screen.
* It **cannot** speak to amyloid gating (section 1).
* It **cannot** distinguish LC-first from other models outside these three
  classes.
* The worlds are built around the dry-run parameter values; real data could
  be harder.
