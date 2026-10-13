# Pre-registration — Module 6: network diffusion (LC vs entorhinal seed)

**Status: COMPLETE DRAFT with one open item that only the researcher can
close: data access.** The connectome needs HCP data-use registration, and
tau-PET needs an ADNI application with verified institutional
affiliation. All methods are fixed below. Freeze only once §3 names the
exact datasets and releases. Claude must not run `make freeze-prereg`.
Configuration: `config/network.yaml`.

> **DO NOT FREEZE until §3 names the datasets and releases** (data access is the researcher's).

## Change log

| Date | Change |
|---|---|
| 2026-10-13 | Complete draft. Replaces `network.DRAFT.md`. Decisions: parcellation including an LC node, the tau group, control seeds, the significance level, and the crossover rule. |

## 1. Prediction and what counts against it (stated first)

**Prediction:** the LC seed fits the earliest-stage regional tau pattern
better than the entorhinal (EC) seed. Let Δr = r_LC − r_EC, where each
seed's Pearson r with observed regional tau is maximised over diffusion
time. The prediction is Δr > 0 with spin-test p < 0.05.

**Counts against:**

- Δr ≤ 0; or
- Δr > 0 but not distinguishable from the spin null (p ≥ 0.05).

**Uninformative, not support:** both seeds' best fits lie beyond the
crossover time, where the diffusion patterns from different seeds have
become nearly identical. The data then cannot discriminate.

## 2. Methods fixed in advance

- **Model:** ẋ = −βLx plus growth, solved by eigendecomposition of the
  normalised graph Laplacian.
- **Diffusion time τ:** 400 log-spaced values from 0.01 to 50. Each
  seed's r is maximised over τ.
- **Significance:** only from spin nulls.
  - 10,000 rotations with Hungarian reassignment. The right hemisphere is
    mirrored, and the τ maximisation is repeated inside every null.
  - Validated on a known null: 6% false positives vs 57% for naive
    shuffling (`tests/validation/test_network.py`).
- **Parcellation:** Desikan–Killiany cortex (68), plus subcortical nuclei
  and an **LC node** from a published brainstem atlas: the LC probability
  map of Keren et al. 2009, thresholded at 25%, or the Harvard Ascending
  Arousal Network atlas. The choice is logged.
- **Structural connectivity:** group-average HCP tractography to all
  nodes.
  - *Limitation:* tractography to the LC is unreliable.
  - *Sensitivity:* LC connectivity replaced by the published LC efferent
    projection targets (binary).
- **Control seeds:** pericalcarine cortex (V1, a late-tau region) and
  the precuneus (default-mode hub).
- **Crossover time:** the first τ at which the LC and EC predicted
  patterns correlate at ≥ 0.95. Report it.

## 3. Data (to be named before freezing)

- **Tau:** ADNI flortaucipir, regional SUVR (inferior cerebellar grey
  reference, partial-volume corrected). Use the earliest-stage group:
  amyloid-positive, cognitively unimpaired participants in PET Braak
  stage I–II. N is reported.
- **Connectome:** HCP (release and pipeline named at download).

## 4. Limitations, stated in advance

- Group-level tau maps.
- LC tractography is unreliable.
- Late stages carry limited seed information (high-eigenvalue modes
  decay fast).
- Diffusion may not be the mechanism of spread.
