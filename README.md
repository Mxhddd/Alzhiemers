# LC-First

A computational test of the hypothesis that the **locus coeruleus (LC) sets the
rate of Alzheimer's disease progression** rather than being its origin. The
full scientific brief, ground rules and build order are in
[`PROJECT_BRIEF.md`](PROJECT_BRIEF.md). Read it first.

**The two results that would kill the hypothesis** (they lead every write-up):

1. LC function does not predict the rate of tau accumulation in amyloid-positive
   people, after adjusting for baseline tau, age and APOE4.
2. LC integrity does not mediate the association between sensory/social risk
   factors and tau pathology.

## Status

| Session | Module | Status |
|---|---|---|
| 0 | Scaffold and provenance | **built** |
| 1 | Validation harness + MR validation (LDL–CAD via HMGCR) | **built, gate not yet passed**: the published LDL–CAD estimate and variant tables are still to be extracted (Module 1 work, not yet done) |
| 2 | Module 1 — drug-target Mendelian randomization (R) | **built; pre-registration complete** (`prereg/mr.md`: cis window ±100 kb, sensitivity ±1,000 kb; positive control per target from GWAS Catalog evidence; LGMN/SLC18A2/MAOA untestable). Cannot run: validation gate, downloads, TwoSampleMR cross-check and the ±1,000 kb sensitivity code are still to do |
| 3 | Module 2 — Markov base-rate model (Python) | **revised (option c) and ready to freeze.** The amyloid-gating criterion was dropped, with a change log in `prereg/baserate.md`. A is age-adjusted, and its HR is "descriptive, not a test of gating". The literal A is a sensitivity analysis "confounded with age". `make validate-baserate` 13/13, with guards that fail loudly if gating is re-enabled. Staged, checkpointed run; uncertainty-aware falsification test. Blocked on the Braak 2011 stage-by-age table (publisher Cloudflare) |
| 4 | Module 5 — meta-analysis and triangulation (R) | **extraction built from open sources; validation 4/4 PASS** (Moraros 2017 reproduced; anticholinergic control recovered); pre-registration complete. The real run is refused at the freeze gate, by design |
| 5 | Module 4 — ODE dynamical model (Python) | **machinery built; pre-registration complete** (`prereg/dynamics.md`). The analytic fixed points are derived and machine-checked. The derivation showed the proposal's damage term made bistability impossible, so it was generalised. The LC model itself is not yet coded |
| 6 | Module 3 — mediation analysis (Killer 2) | **built; pre-registration complete except the cohort** (DO NOT FREEZE until the researcher names it); 4/5 validation checks pass |
| 7 | Module 6 — network diffusion | **built; spin test validated on a known null** (6% vs 57% false positives); pre-registration complete except the datasets (DO NOT FREEZE until named) |

Pre-registrations: `make todo` lists which are ready to freeze. `make prereg-package`
builds the OSF upload package. **Freezing is the researcher's step only.**

## Extracted published values (`data/extracted/`)

Numbers taken from published papers by Claude live in `data/extracted/`,
committed, one source per row (paper, table/figure, exact quantity) and a
`verified` flag that starts `false`. `data/extracted/SOURCES.md` lists the
papers, how each was accessed, oddities spotted in the source, and what
could not be found. Any figure or table built on unverified values carries a
red UNVERIFIED banner (figures) or a leading `UNVERIFIED_INPUTS` column
(tables); the run-log footer lists them. Derivations Claude made for the
researcher to check are in `manuscript/derivations/`, each marked as not
independently checked.

## What the researcher must supply

`make todo` prints it, computed from the repository (every `null` field in
`config/` and `tests/validation/targets/`, every referenced data file that is
absent, and every pre-registration still in draft). Nothing on that list is
ever filled in by code.

## Quick start

Requirements: Python 3.11 with [uv](https://docs.astral.sh/uv/), R ≥ 4.3, git, make.

```bash
make setup              # install pinned Python + R environments
make test               # unit tests (Python and R)
make provenance-check   # write a run log with a full provenance header
make validate           # reproductions of published results (fails until targets are supplied)
make all                # everything built so far
```

`make all` currently **fails at `make validate`, on purpose**: the LDL–CAD
target has not been transcribed yet, and an unvalidated MR module must not run.

## How the shared infrastructure works (`src/common/`)

| File | What it does |
|---|---|
| `paths.py` | Canonical locations of `config/`, `data/`, `sim/`, `results/`, … |
| `provenance.py` | Git commit + dirty-tree status, SHA-256 of files, package versions |
| `runlog.py` | `start_run(...)`: opens `results/logs/<module>/<UTC time>_<commit>.log` with a YAML provenance header; a drop ledger (`run.drop(stage=, item=, reason=)`); output hashes; a footer with status and drop counts. Errors are logged as `FAILED` and re-raised, never swallowed. |
| `provenance.R` | The same log format for the R modules (`lc_with_run`, `lc_drop`, …). |
| `config.py` | Strict YAML: no default values, missing keys raise, **duplicate keys raise** (plain YAML silently keeps the last one). A unit test fails if any setting in `config/` lacks a comment saying why. |
| `io.py` | `analysis_input(path, provided_by=...)`: refuses anything under `sim/` or `tests/` (simulated data never reaches an analysis), refuses paths outside `data/`, and when a hand-supplied file is missing raises an error naming who must provide it. |
| `prereg.py` | The pre-registration gate — see [`prereg/README.md`](prereg/README.md). Called from R via `lc_require_preregistration()` so there is one implementation. |
| `validation.py` | The validation gate. `make validate` records PASSED/FAILED plus a **code fingerprint** (hash of the shared code, the module's code, the validation harness and its method config) in `results/validation/<module>.status.yaml`. `require_validated(module)` refuses to run a module whose last validation failed or whose code has changed since. |
| `io.R` | R version of the input guard. |
| `plotting.py` | Shared figure theme: bundled font, colour-blind-checked fixed palette, figures stamped with the commit hash and byte-reproducible. |

## Module 1 building blocks (`src/mr/`, Session 1)

| File | What it does |
|---|---|
| `harmonise.R` | `mr_harmonise()` puts exposure and outcome on the same effect allele: direct, swapped, strand-flipped, palindromes resolved by allele frequency (only if MAF ≤ the configured threshold), everything else dropped **with a reason**. Stops on genome-build mismatch, duplicate IDs or missing values. |
| `estimators.R` | Wald ratio, IVW (multiplicative random effects), MR-Egger, weighted median (Bowden 2016, parametric bootstrap). Base R, following TwoSampleMR conventions; settings from `config/mr_methods.yaml`. Too few variants gives a "not estimable" row, never a number. |

Why base R and not TwoSampleMR: the package could not be installed in the build
environment (CRAN/GitHub R repositories blocked). The estimators are tested
against hand arithmetic, independent `lm()` fits, allele-recoding invariance,
CI coverage on simulated data, and Bowden 2016's robustness property. A
cross-check against TwoSampleMR on identical inputs is pending network access.

## Module 1 pipeline (`make mr`, Session 2)

`src/mr/run_mr.R` refuses to start unless `make validate` passed on the
current code **and** `prereg/mr.md` is frozen (the complete draft is
`prereg/mr.md`, awaiting the researcher's freeze). Then, per target and exposure source:

1. `mr_gene_region()` — Ensembl GRCh37 coordinates via biomaRt, release logged,
   cached in `data/interim/`.
2. `mr_read_sumstats()` — column-mapped reading; eQTLGen Z-scores converted to
   beta/SE (Zhu 2016); build asserted; every removed row logged.
3. `mr_build_instruments()` — p, MAF, F filters, PLINK-based LD clumping.
4. **Positive control (hard gate) before the AD outcome is read** — a target
   with no control, or a failing one, is excluded and logged.
5. Harmonise → IVW / Egger / weighted median (or Wald ratio) → ABF coloc over
   the cis window → Steiger → proxy-free sensitivity and FinnGen replication
   when supplied.
6. `results/tables/mr_summary.csv` (one row per gene × source × method, with
   F, Egger intercept, PP4, Steiger, positive-control status, concordance with
   the prediction), `results/logs/mr_dropped.md`, and
   `results/figures/mr_forest.{pdf,png}` (pre-registered direction shaded;
   excluded targets kept as labelled gaps).

Not implemented, reported as such in the output: SuSiE-coloc (needs susieR).

## The validation harness (`tests/validation/`)

A validation target is a YAML file in `tests/validation/targets/` holding a
published result **transcribed by the researcher**. Fields left `null` make
`make validate` fail with a list of what is missing — never skip.

To unblock Module 1, fill in `tests/validation/targets/mr_ldl_hmgcr_cad.yaml`
(citation, method, estimate and CI as printed, scale, conversion factor, genome
build) and put the per-variant tables at
`data/raw/validation/mr_ldl_hmgcr_cad_{exposure,outcome}.csv`
(`SNP, effect_allele, other_allele, beta, se, eaf, source`). The test passes
only if our estimate, by the paper's method, falls inside the published CI.
Results go to `results/tables/validation_mr.csv`; both directions (ours in
theirs, theirs in ours) are reported.

### What a run log header looks like

```yaml
---
lc_first_run_header: 1
module: provenance_check
language: python
started_utc: 2026-10-12T22:30:00.123456Z
git:
  commit: 3f2c…(40 hex)
  dirty: false
  uncommitted_files: []
config:
  path: config/provenance_check.yaml
  sha256: 9a1e…
lockfiles:
  uv.lock: …
  renv.lock: …
seed: 20261012
packages:
  python: 3.11.17
  pyyaml: 6.0.3
  matplotlib: 3.11.2
platform: Linux-…
command: […]
preregistration: null
---
```

A dirty working tree is recorded (with the list of changed files) and flagged
with a `WARNING` line: such a result cannot be reproduced from the commit alone.

## Repository layout

```
config/        one YAML per analysis; every value commented
prereg/        frozen pre-registrations (OSF link + local copy)
data/raw/      hand-entered or downloaded inputs — never committed, never edited
data/interim/, data/processed/
src/common/    shared infrastructure (above)
src/<module>/  one directory per module
tests/unit/, tests/validation/
sim/           simulated-data generators (never read by an analysis)
results/       figures/, tables/, logs/ — generated by `make all`, not committed
manuscript/
env/           environment notes
```
