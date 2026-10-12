# PROJECT_BRIEF.md — LC-First: a computational test of the locus coeruleus rate-controller hypothesis

> **How to use this file.** Put it at the root of a new repository and open Claude Code there. Start with: *"Read PROJECT_BRIEF.md. Confirm you understand the ground rules and the build order, then do Session 0 only. Stop and report before starting Session 1."* Work one session at a time. Do not ask for the whole repo in one go.

## 1. What this project is

I am testing a specific hypothesis about Alzheimer's disease using only public data, mathematics, and code. I have no laboratory, no institutional data access, and no clinical collaborator yet. The goal of this repository is to produce three or four pre-registered, reproducible analyses that are good enough to publish as preprints and to earn a collaborator.

You are the engineer. I am the researcher. I do the mathematics by hand and check your output against it. You build pipelines that are auditable, reproducible, and honest about failure.

I am not a professional programmer or statistician. Explain what you build. When a statistical choice is available, tell me what you picked and why, rather than silently choosing. If I ask for something that is methodologically wrong, say so.

## 2. The scientific hypothesis

Read this section carefully. Several design decisions below only make sense in its light.

### The claim

**The locus coeruleus (LC) is not the origin of Alzheimer's disease. It is the rate controller of its progression.**

The LC is a small brainstem nucleus, the brain's main source of noradrenaline. Tau pathology appears there earlier than anywhere else — in most people by early adulthood. But most people never develop Alzheimer's. So LC tau cannot be what distinguishes cases from non-cases.

What varies between people is LC *regulatory function*: the ratio of phasic (task-locked) to tonic (background) firing, the integrity of the α2A autoreceptors that restrain firing, and the amplitude of noradrenaline oscillations during non-REM sleep that drive glymphatic clearance.

The claim is that in someone who already has amyloid pathology, the *rate* at which tau spreads to cortex is set by how well the LC is still regulating itself. Good regulation, slow spread, death from something else first. Poor regulation, fast spread, dementia.

### Why this framing rather than the usual one

| Standard LC hypothesis | This version |
|---|---|
| Claims the LC is where AD starts | Claims the LC sets how fast it progresses |
| Contradicted by near-universal LC tau | Predicts near-universal LC tau; it is the substrate, not the discriminator |
| Needs amount of noradrenaline to be high or low | Needs regulation to be intact; permits opposite drugs to help different patients |
| Origin is near-unprovable in humans | Rate is directly measurable in longitudinal data |

### Six dimensions the model should unify

1. **Modifiable risk converges here.** Hearing loss, vision loss, social isolation, physical inactivity, depression and poor sleep all plausibly act on noradrenergic drive or its regulation. The LC is a candidate common path for a list that currently has no theory.
2. **Cognitive reserve becomes physical.** Reserve is usually measured by proxy (education, occupation). Under this model, reserve *is* LC regulatory capacity, measurable by neuromelanin MRI and pupillometry.
3. **Sex and menopause.** Women are roughly two-thirds of cases. Estrogen loss as a sex-specific regulatory shock at a specific life stage predicts both the excess risk and its timing.
4. **APOE4 as a noradrenergic handling defect.** APOE4 inhibits VMAT2 in LC neurons, pushing noradrenaline into the cytosol, where MAO-A converts it to DOPEGAL, which activates asparagine endopeptidase (AEP/legumain), which cleaves tau at N368 into propagation-prone forms.
5. **Three timescales.** Clearance is set nightly; regulatory state shifts over months; tau accumulates over decades. A slow variable driven by a fast one.
6. **The escapees.** People with heavy pathology and intact cognition should show preserved LC function, not less pathology.

### The two results that would kill it

These go at the front of every write-up, not the back:

- **Killer 1.** LC function does not predict the rate of tau accumulation in amyloid-positive people, after adjusting for baseline tau, age and APOE4.
- **Killer 2.** LC integrity does not mediate the association between sensory/social risk factors and tau pathology.

If either comes back null, the hypothesis is wrong in its current form and we say so in print.

### The honest risk

A theory that explains sex differences, reserve, genetics, lifestyle risk and resilience is in danger of explaining everything and forbidding nothing. The two killers above are the guard against that. If you ever find yourself building an analysis that cannot come out against the hypothesis, stop and tell me.

## 3. Ground rules — these override everything else

1. **Validation before science.** Every analysis module ships with a test that reproduces a published result. A module that cannot recover a known answer must not be used on a new question. Build the validation test first, not last.
2. **Never invent data.** Do not generate, guess, interpolate or "reconstruct" any empirical number — not autopsy counts, not effect sizes, not prevalences. If an analysis needs a number I have not supplied in a data file, stop and ask me for it. Simulated data is allowed only in files under `tests/` or `sim/`, must be generated by a visible script, and must never be read by an analysis pipeline.
3. **No silent defaults.** Every threshold, prior, weighting scheme and model choice lives in a YAML config with a comment saying why. At the end of each session, list every statistical choice you made that this brief left open.
4. **Fail loudly.** No `try/except` or `na.rm = TRUE` that swallows a problem and returns a plausible number. Missing data raises. Dropped records are counted and logged with reasons.
5. **Reproducible end to end.** `make all` regenerates every figure and table from raw inputs. No number is ever typed into a manuscript by hand. Seeds are set and recorded. Package versions are pinned (`requirements.txt`, `renv.lock`).
6. **Provenance in every log.** Each run writes a log header with the git commit hash, the config file hash, the timestamp, and the versions of key packages.
7. **Pre-registration is a gate.** Modules that test the hypothesis refuse to run on real target data unless a pre-registration file exists in `prereg/` with a frozen timestamp. Enforce this in code, not by good intentions.
8. **Honest output.** Dropped targets, failed instruments and wide confidence intervals are results and are reported. If a result looks suspiciously clean, say so.

## 4. Repository structure

```
lc-first/
  PROJECT_BRIEF.md       # this file
  README.md              # what each module does and how to run it
  Makefile               # per-module targets plus `make all`
  config/                # one YAML per analysis
  prereg/                # frozen pre-registrations (OSF links + local copies)
  data/
    raw/                 # downloaded or hand-entered; never edited, never committed
    interim/
    processed/
  src/
    common/              # IO, logging, provenance, plotting theme
    mr/                  # Module 1 (R)
    baserate/            # Module 2 (Python)
    mediation/           # Module 3 (R or Python)
    dynamics/            # Module 4 (Python)
    metaanalysis/        # Module 5 (R)
    network/             # Module 6 (Python)
  tests/
    validation/          # reproductions of published results — the important ones
    unit/
  sim/                   # simulated-data generators, clearly separated
  results/
    figures/
    tables/
    logs/
  manuscript/
  env/
```

Python 3.11+ with `uv`. R 4.3+ with `renv`. Git from the first commit.

## 5. Modules

Build in the order given in §6. Each module spec below is the contract.

### Module 0 — scaffold and provenance

Repository structure above; logging utility that records commit hash, config hash and package versions; shared plotting theme; Makefile; README skeleton; `.gitignore` that excludes `data/raw/`. No analysis code.

### Module 1 — Drug-target Mendelian randomization (R)

**Question.** Do genetic variants that mimic noradrenergic drugs change Alzheimer's risk?

**Targets and pre-registered predictions.** Under the rate-controller framing, the predictions are about *regulation*, which matters for how results are interpreted:

| Gene | Protein | Drug analogue | Predicted direction of lower activity |
|---|---|---|---|
| SLC6A2 | Noradrenaline transporter | Atomoxetine | Protective |
| ADRB1 | β1 adrenergic receptor | Lipophilic β-blockers | Protective |
| ADRA1A | α1A receptor | Prazosin, doxazosin | Protective |
| ADRA2A | α2A autoreceptor | Guanfacine (agonist) | **Harmful** (loss of restraint) |
| LGMN | Asparagine endopeptidase | none in humans yet | Protective |
| SLC18A2 | VMAT2 | none | **Harmful** |
| MAOA | Monoamine oxidase A | MAO-A inhibitors | Protective |
| DBH | Dopamine β-hydroxylase | Disulfiram, nepicastat | Exploratory, direction not pre-specified |
| CACNA1D | Cav1.3 calcium channel | Isradipine | Protective |
| HCRTR2 | Orexin receptor 2 | Suvorexant, lemborexant | Exploratory |

Note: MAOA is X-linked and may be absent from some summary files — check and report rather than working around it.

**Data.**

- Exposure instruments: cis-eQTLs from eQTLGen (blood) and MetaBrain (cortex); cis-pQTLs from UK Biobank Pharma Proteomics Project / deCODE where available.
- Outcome: Bellenguez et al. 2022, GWAS Catalog GCST90027158. **Important:** the Catalog hosts the stage 1 discovery results only, on build GRCh37, and they include UK Biobank proxy cases. Run a sensitivity analysis excluding proxy-case contributions where possible.
- Replication: FinnGen public Alzheimer's endpoints.

**Config shape** (`config/mr_targets.yaml`): outcome id and build; per-target window, predicted direction, positive control; global settings for F-statistic minimum (10), LD r² (0.001), LD window (10000 kb), MAF minimum (0.01), coloc PP4 threshold (0.8).

Gene coordinates are not hardcoded. A `gene_region()` function looks them up from the Ensembl GRCh37 archive via biomaRt and logs the release used. A build mismatch between exposure and outcome produces wrong instruments with no warning, so assert the build everywhere.

**Functions.** `fetch_exposure()`, `build_instruments()` (F > 10, report drops), `harmonise()` (report palindromes resolved and dropped), `run_mr()` (IVW, MR-Egger, weighted median, Wald ratio), `run_coloc()` (coloc, SuSiE-coloc for multiple signals), `steiger_test()`, `positive_control()` (a hard gate — an instrument that fails is excluded automatically), `report_target()`.

**Validation test.** Reproduce a published two-sample MR: LDL cholesterol via HMGCR variants on coronary artery disease, asserting the estimate falls within the published confidence interval. Also: a harmonisation unit test with hand-built strand flips and palindromes; a null test using variants in no gene region.

**Outputs.** `results/tables/mr_summary.csv` (one row per gene × method, with F-statistics, Egger intercept, PP4, Steiger direction, positive-control status); a forest plot marking the pre-registered predicted directions so a reader can see which held; `results/logs/mr_dropped.md` listing every excluded target with its reason.

### Module 2 — Markov base-rate model (Python)

**Question.** If nearly everyone has LC tau by 40, what has to be true for only a minority to develop Alzheimer's?

**Data.** I will transcribe the Braak staging counts by hand into `data/raw/braak_counts.csv` with columns `age_lower, age_upper, stage, n, source`. **Do not create this file, scrape it, or fill in any values.** If it is missing, stop and tell me.

**Model.** Continuous-time Markov chain over tau stages with age-dependent rates, dp/dt = p(t) Q(t), q_{k,k+1}(t) = exp(α_k + β_k t). Three variants sharing one interface so switching is trivial and none is favoured by construction:

- **A, gated:** LC→cortex rate depends on amyloid status, μ(t) = μ₀ + γ·A(t)
- **B, independent:** cortical tau arises unrelated to LC stage
- **C, multi-hit:** k independent hits; incidence scales as t^(k−1)

**Functions.** `build_generator()` (assert rows sum to zero), `stage_probabilities()` (matrix exponential), `log_likelihood()`, `fit_mle()`, `fit_bayesian()` (PyMC; report R-hat and divergences), `predict_incidence()`, `compare_models()` (k-fold cross-validation by age band).

**Age binning.** I will supply the bin width from an error analysis I do by hand — Q(t₁) and Q(t₂) do not commute, so the matrix exponential is only an approximation over a bin. Ask me for the number; do not pick one.

**Validation tests.** (a) Simulation recovery at three sample sizes. (b) Analytic check: a two-state chain with constant rate must match 1 − exp(−λt) to 8 decimal places. (c) Stage probabilities sum to 1 at every age. (d) **The falsification test:** fit only on autopsy data, then compare the predicted dementia incidence curve to published incidence the model never saw.

**Outputs.** Parameter table with credible intervals; predicted vs observed stage prevalence; predicted vs observed incidence; the estimated γ with uncertainty. I will have derived an expected γ by hand from published amyloid hazard ratios — flag loudly if the fit disagrees.

### Module 3 — Mediation analysis (NEW, and the sharpest test)

**Question — this is Killer 2.** Does LC integrity mediate the association between modifiable risk factors (hearing loss, social isolation, physical inactivity) and tau pathology?

This module is the most novel thing in the repository. It tests dimension 1 of the hypothesis directly, and it is the analysis most likely to be genuinely new.

**Data.** Requires a cohort with LC measures, risk-factor measures and tau outcomes. I do not yet have access. Build the module against simulated data in `sim/`, with the real-data path clearly marked and failing loudly if the file is absent. Do not substitute any public dataset without asking me first.

**Method.** Causal mediation (`mediation` in R or DoWhy/statsmodels in Python): estimate the average causal mediation effect, the direct effect, and the proportion mediated, with bootstrapped confidence intervals. **Include sensitivity analysis for unmeasured confounding** (the ρ parameter) — mediation results are notoriously fragile to it, and a result without this analysis is not publishable.

**Pre-specify the directed acyclic graph** before fitting anything, as a file in `prereg/`. State every assumption that identification requires. Mediation analysis makes strong causal assumptions from observational data, and the write-up must say so plainly.

### Module 4 — ODE dynamical model (Python)

**Question.** Can one parameter set make both atomoxetine (raises noradrenaline) and prazosin (blocks α1 signalling) protective? Is the system bistable?

**Structure.** `model.py` (pure `rhs(t, y, params)`), `params.py` (dataclass carrying units, source and confidence per parameter), `simulate.py`, `fixedpoints.py` (roots plus Jacobian eigenvalues), `bifurcation.py`, `sensitivity.py` (SALib Sobol), `figures.py`.

**Parameter discipline.** Every parameter carries `source` and `confidence`. Parameters marked `weak` are never reported as fitted values — they are swept.

**Analyses.** Baseline trajectory; fixed points and stability; bifurcation over the second-hit parameter; drug simulation; dose–response sweep (the U-shape test); Sobol sensitivity; and critically, a **Latin hypercube robustness run** reporting the fraction of plausible parameter space in which the U-shape appears. If it is 3%, the figure caption says 3%.

**Validation.** Non-negativity and conservation; an analytic fixed point I derive by hand; solver independence (two stiff solvers, tolerances tightened 10×); degenerate case with zero coupling reducing to known exponentials.

### Module 5 — Meta-analysis and triangulation (R)

**Question.** Across published cohorts, do noradrenergic drug classes show a consistent, brain-penetrance-graded association with dementia incidence?

**Data.** I extract studies by hand into `data/raw/extraction.csv`. You build the analysis around it; you never populate it. Columns include study id, drug class, drug, `comparator_type` (`active` or `non_user`), comparator, effect measure, estimate, CI bounds, covariates, and the seven ROBINS-I domain ratings.

**The informative contrasts** (active comparators, so confounding by indication largely cancels): lipophilic vs hydrophilic β-blockers; terazosin/doxazosin/alfuzosin vs tamsulosin; Z-drugs vs dual orexin antagonists; atomoxetine vs stimulants. Plus **anticholinergic burden as a positive control** — if the pipeline cannot recover that known association, it is broken.

**Functions.** `load_extraction()` (validate hard, fail on any gap), `to_log_scale()`, `pool()` (`metafor::rma` random effects, with prediction intervals), `subgroup()`, `publication_bias()` (funnel, Egger, trim-and-fill), `leave_one_out()`, `forest_plot()`, and `triangulate()` — join pooled epidemiological estimates to Module 1 MR estimates by target.

**The triangulation figure is the centrepiece** and does not exist in the literature: per target, the pooled observational estimate, the MR estimate, the mechanistic prediction, and a concordance flag, on a common log scale. Missing estimates are shown as gaps, not omitted.

**Validation.** Reproduce a published pooled estimate. Hand-computable three-study example (I will supply the arithmetic). Simulated homogeneous data returns I² near zero.

### Module 6 — Network diffusion (Python)

**Question — this is Killer 1's geometric cousin.** Does tau spread better fit an LC seed than an entorhinal seed?

**Model.** ẋ = −βLx plus growth, on a structural connectome, solved by eigendecomposition of the graph Laplacian. Compare seeds (LC, entorhinal, controls) on correlation with observed regional tau, maximised over diffusion time.

**The statistical trap, non-negotiable.** Brain maps are spatially autocorrelated, so unrelated maps correlate far above chance. Null distributions come from spatially-constrained surrogates (spin tests), never from shuffled regions. Build the spin test before the comparison.

**Expect limited power at late stages:** high-eigenvalue modes decay fast, so late-stage patterns may be nearly seed-independent. Run against the earliest-stage patterns available and report the crossover time.

**Data.** Connectomes from the Human Connectome Project; tau-PET from ADNI (requires an application with verified institutional affiliation — I may not have this; build against simulated data and mark the real path clearly).

## 6. Build order and stop points

Stop and report after each. Do not continue to the next without my go-ahead.

| Session | Build | Gate before proceeding |
|---|---|---|
| 0 | Module 0 scaffold | Repo runs, logging records provenance |
| 1 | `tests/validation/` harness + MR validation test | LDL–CAD estimate falls in the published CI |
| 2 | Module 1 | ADRB1 positive control passes (instrument moves resting heart rate) |
| 3 | Module 2 | Simulation recovery passes; `braak_counts.csv` supplied by me |
| 4 | Module 5 | Reproduces a published pooled estimate; anticholinergic control recovered |
| 5 | Module 4 | Analytic fixed point matches my hand derivation |
| 6 | Module 3 | DAG pre-registered; sensitivity analysis implemented |
| 7 | Module 6 | Spin test validated against a known null |

Session 1 produces no new science and must not be skipped. It is the session that makes everything after it believable.

## 7. What you must not do

- Do not invent, impute or reconstruct empirical numbers. Ask me.
- Do not substitute a different dataset for one I named without asking.
- Do not run a hypothesis-testing module on real target data before `prereg/` holds a frozen pre-registration. Enforce this in code.
- Do not silently drop records, instruments, studies or trials. Count, log, report.
- Do not tune a model until it agrees with the hypothesis. If a result contradicts the hypothesis, that is the finding.
- Do not write analysis code that cannot produce a result against the hypothesis.
- Do not claim a test passes without running it.

## 8. Reporting back

After every session, give me:

1. What you built and how to run it.
2. Every statistical choice you made that this brief left open.
3. Any place where you were tempted to guess a value, and what you need from me instead.
4. What the next session should do, and what could go wrong in it.

## 9. Start
