# Does early cortical tau spread from the locus coeruleus or from the entorhinal cortex? A pre-registered network diffusion test

Mohammad Hassan, independent researcher · October 2026

This is my plan for Module 6 of my project: asking whether the earliest cortical tau pattern seen on tau PET is better explained by diffusion along the brain's structural connections starting from the locus coeruleus (LC) or starting from the entorhinal cortex. I'm writing it down before I have any of the data, so that what counts as support and what counts as failure is fixed in advance. **I have not been granted access to ADNI or to the Human Connectome Project (HCP) yet. Nothing from either has been requested, downloaded, opened or seen.** The code and settings are in the project repository (`config/network.yaml`, `src/network/`), at the commit where I register this file.

## Tools

I used Claude, Anthropic's AI model, through Claude Code to write the code, extract the published values, run the simulations and write the first version of this text. I told it to decide methodological questions for me and to explain each decision, which it did in the config file (every setting carries its reason) and in the change log below. By registering this plan I adopt those decisions as mine. Claude read only public documentation pages about ADNI and HCP (listed in section 3), never any data. Because the ADNI Data Use Agreement forbids putting participant-level ADNI data into third-party AI tools, Claude will never see ADNI participant-level data: I will run that step myself (section 3). The analysis code refuses to run on real data until this plan is frozen.

## Change log

This document has never been frozen. Every change below was made before any ADNI or HCP data existed on my side; none is a change made after seeing a result.

| Date | Change | Reason |
|---|---|---|
| 2026-10-13 | First complete draft: LC versus entorhinal seed, spin test, control seeds, crossover time. Data access left open. | — |
| 2026-10-21 | **The entorhinal regions are left out of every correlation** (and every seed's own regions). The first draft's design is kept as sensitivity analysis 4. | Early tau sits in the entorhinal cortex. Scored on its own value, the entorhinal seed starts with a lead no LC-seeded prediction can match, because the LC is not a cortical region and has no value in the map. That design leaned against the hypothesis by construction. |
| 2026-10-21 | **Primary group changed from "amyloid-positive, cognitively normal, PET Braak stage I–II" to all amyloid-positive cognitively normal people**, as a W-score map against amyloid-negative cognitively normal people. Amyloid-positive MCI is sensitivity analysis 9. | Selecting people by tau in the entorhinal region selects the map in favour of the entorhinal seed. PET Braak staging would also need regional cut-offs that I cannot source before access. The W-score against amyloid-negative people removes the off-target binding pattern shared by both groups. |
| 2026-10-21 | **Verdict rule rewritten (section 2).** The LC − entorhinal difference is now tested one-sided against the spun maps' own differences, in each direction. "Supports" now also needs the LC seed to beat the spin null by itself and to fit before the crossover time. A positive but non-significant difference is now "inconclusive", not "against". | The first draft's code compared the observed difference with a null that had been re-centred to zero, which mixes two different questions. A positive difference that the spin null cannot distinguish from chance is lack of evidence, not evidence against. A difference of zero or below still counts against. |
| 2026-10-21 | **Operating characteristics measured on simulated data** (section 2): 1 "supports" in 200 null worlds. | To know what the rule does before using it. Simulated data only. |
| 2026-10-21 | **Data named (section 3):** the ADNI tables and variables, the HCP release, subjects and files, the LC atlas, and the tractography settings. The LC atlas is now Ye et al. 2021 (7T); the first draft's "Keren 2009 at 25%" was wrong, since that paper distributes binary 1 SD and 2 SD templates. Uniform LC projection replaces "published LC efferent targets" as the answer to unreliable LC tractography. | Named from public documentation only. A target list for LC efferents would have needed published values I have not extracted; the uniform projection needs none and is labelled an assumption. |

## 1. Question

In Braak's account, tau appears in the LC before it appears in the cortex. If tau then spreads along neural connections, the earliest cortical tau pattern should look like diffusion along the LC's own connections. If cortical tau instead starts in the entorhinal cortex, it should look like diffusion from there. Which seed explains the earliest cortical tau pattern better: the LC or the entorhinal cortex?

## 2. Prediction, and what would count against it

I am stating these first, on purpose.

**Definitions.** For each seed, r is the correlation (Pearson) between the predicted and the observed regional tau map, maximised over diffusion time (section 4). Δr = r_LC − r_EC, where EC is the entorhinal seed. Both are computed on the same cortical regions: the 68 Desikan–Killiany regions minus the two entorhinal regions. Every spin p-value is one-sided, uses 10,000 spins, and adds 1 to numerator and denominator.

* p_LC: how often a spun map fits the LC seed at least as well as the real map does.
* p(LC better): how often a spun map gives a difference at least as large as Δr.
* p(EC better): how often a spun map gives a difference at least as small as Δr.
* Crossover time: the first diffusion time at which the LC and entorhinal predictions correlate at 0.95 or more. After it, the two predictions are nearly the same pattern.

**Prediction:** Δr > 0, and it is larger than spatial autocorrelation alone produces.

**Verdict, with the rules applied in this order** (`src/network/verdict.py`):

1. **Uninformative:** both seeds' best fits lie at or beyond the crossover time. The data then cannot tell the seeds apart.
2. **Against:** Δr ≤ 0. The prediction failed. It is labelled *against (beyond the spin null)* if, in addition, p(EC better) < 0.05.
3. **Supports:** Δr > 0, p(LC better) < 0.05, p_LC < 0.05, and the LC seed's best fit lies before the crossover time.
4. **Inconclusive:** anything else. This includes a positive Δr that the spin null cannot distinguish from chance.

**The result counts against the hypothesis** if the primary verdict is *against*. **It counts in its favour** only if the primary verdict is *supports*. I will call either verdict **robust** only if sensitivity analyses 1, 2, 3, 5, 6, 7 and 8 (section 6) all point the same way (supports / against / neither); an analysis that cannot be computed makes the verdict not robust.

**Operating characteristics, measured on simulated data before any real data** (`sim/network_calibration.py`; full table in `manuscript/notes/module6_calibration.md`). The synthetic worlds have 34 parcels per hemisphere, a distance-dependent cortical connectome, and one LC node per hemisphere wired to its whole hemisphere with random weights. The analysis code is the real one, with 1,000 spins and 200 data sets per world. "Noise" is a smooth random map with 0.5, 1 or 2 times the signal's spread.

| Simulated world | Supports | Against (both labels) |
|---|---|---|
| Map unrelated to the brain graph | 1/200 | 98/200 |
| Map spread from the LC, noise 0.5 / 1 / 2 | 199/200 · 180/200 · 106/200 | 0/200 · 2/200 · 21/200 |
| Map spread from the entorhinal cortex, noise 0.5 / 1 / 2 | 1/200 · 2/200 · 3/200 | 159/200 · 127/200 · 117/200 |

So in these worlds a false "supports" was rare (at most 3 in 200), and an unrelated map gives "against" about half the time, because a difference of zero or below counts against whatever its p-value. These numbers describe synthetic worlds. They are not guarantees for the real data, where the LC's connections come from tractography that is known to be unreliable.

**What a result against would and would not mean.** In Braak's sequence, tau that starts in the LC reaches the cortex first in the transentorhinal region. A cortical map produced that way looks entorhinal-seeded. *Against* therefore says that early cortical tau is not explained by diffusion along the LC's own structural connections. It does not rule out an LC origin with entry into the cortex through the transentorhinal region alone. *Supports* would say that the LC's connection profile explains the early cortical pattern better than the entorhinal seed does.

**Control seeds and ranking (descriptive, no verdict).** Two control seeds are fitted the same way: the pericalcarine cortex (V1, late in Braak staging) and the precuneus (a highly connected hub). If a control seed fits better than both the LC and the entorhinal seed, I will say that the diffusion model with these seeds does not capture the early pattern. I will also rank every bilateral pair of nodes as a seed (each evaluated without the entorhinal regions and without its own regions) and report where the LC falls.

## 3. Data

**Access status.** I have not applied to ADNI and have not registered with HCP. Nothing has been downloaded. Every file and variable name below comes from public documentation; anything not confirmed there is marked *to be confirmed on access*.

**Tau PET and clinical data (ADNI)**, through the LONI Image and Data Archive after approval of an application under the ADNI Data Use Agreement:

| Table | Use | Variables |
|---|---|---|
| UC Berkeley – Tau PET PVC 6mm Res analysis (`UCBERKELEY_TAUPVC_6MM`) | Primary regional tau | `CTX_LH_<REGION>_SUVR`, `CTX_RH_<REGION>_SUVR` for the 34 Desikan–Killiany regions per side; `SCANDATE`; `TRACER`; `INFERIORCEREBELLUM_SUVR` (reference region: inferior cerebellar grey matter) |
| UC Berkeley – Tau PET 6mm Res analysis (`UCBERKELEY_TAU_6MM`) | Sensitivity 8 (no partial-volume correction) | as above; `qc_flag` (2 pass, 1 partial, 0 fail) |
| UC Berkeley – Amyloid PET 6mm Res analysis (`UCBERKELEY_AMY_6MM`) | Amyloid status | `AMYLOID_STATUS` (the lab's own threshold, whole-cerebellum reference), `SCANDATE` |
| Diagnostic Summary (`DXSUM`) | Cognitive status | `DIAGNOSIS` (1 CN, 2 MCI, 3 dementia), `EXAMDATE` |
| Participant Demographics (`PTDEMOG`) | Covariates | `PTDOBYY` (year of birth), `PTGENDER` (1 male, 2 female) |

*To be confirmed on access:* the exact file names and download dates; that `RID` identifies participants in every table; that the PVC table carries `qc_flag`; the value of `TRACER` for flortaucipir (assumed `FTP`) and of `AMYLOID_STATUS` for positive (assumed 1); the partial-volume-correction method and reference region (in the Berkeley Tau PET Processing Methods PDF, which is only on LONI); whether `PTDOBYY` and `DIAGNOSIS` are filled for every phase. If any of these differs, I will correct the name in the config before running and record it as a deviation; the selection rules below do not change.

**Selection** (`src/network/adni_tau.py`, settings in `config/network.yaml`, key `adni`):

* Tau scans with tracer flortaucipir (other tau tracers are not comparable) and visual QC pass (`qc_flag` = 2). The first such scan per participant.
* Amyloid status from the amyloid scan nearest in time to the tau scan, within 365 days. Diagnosis from the diagnostic summary nearest in time, within 180 days. Age = scan date minus year of birth, assuming birth at mid-year.
* Reference group: amyloid-negative, cognitively normal (CN). Primary target group: amyloid-positive CN. Sensitivity group: amyloid-positive MCI.
* For each region, a linear regression of SUVR on age and sex in the reference group. Each target participant's W-score is (SUVR − predicted) / residual SD. The group map is the mean W-score per region.
* A group map, and the reference group, need at least 30 people (an assumption: smaller group means are dominated by individual noise). Below that, the map is not built and the analysis it feeds is reported as not computable.
* Every exclusion is counted and reported (`results/tables/network_adni_counts.csv`). Nothing is imputed.

**How the ADNI data will be handled.** The ADNI Data Use Agreement (PDF dated April 2026) forbids putting ADNI data into third-party AI tools that may retain it, and forbids redistributing participant-level data. So I will download the ADNI tables and run `src/network/adni_tau.py` on my own computer, with no AI tool involved. Only its group-level outputs (68 regional means per map, and counts) leave that computer. The agreement allows summary data to be reported.

**Structural connectome (HCP).** HCP Young Adult, 1200 Subjects Release (S1200), 3T, open-access tier only, through ConnectomeDB after accepting the Open Access Data Use Terms:

* Subjects: the ConnectomeDB preset group "100 Unrelated Subjects" (100 people). Because the HCP itself selected them as unrelated, no restricted family-structure data are needed.
* Packages: *Diffusion Preprocessed* (`T1w/Diffusion/data.nii.gz`, `bvals`, `bvecs`, `nodif_brain_mask.nii.gz`) and *Structural Preprocessed* (`T1w/aparc+aseg.nii.gz`, `MNINonLinear/xfms/standard2acpc_dc.nii.gz`). File names from the S1200 Reference Manual, Appendix III, Section B.
* I will run the tractography myself, after accepting the HCP terms.

**LC atlas.** Ye et al. 2021 (*NeuroImage* 225:117487), a 7T probabilistic LC atlas distributed in thresholded versions at nitrc.org/projects/lc_7t_prob. Primary: the 5% (liberal) version; sensitivity 6: the 25% (stringent) version (`data/extracted/network_lc_atlas.csv`, unverified). Voxels left of the midline form the left LC node, the rest the right. *To be confirmed on access:* the file names and that the atlas is in MNI152 space as described.

**Sphere coordinates for the spin test.** The centroid of each Desikan–Killiany parcel on the FreeSurfer fsaverage sphere (`lh.sphere`, `rh.sphere`, `lh.aparc.annot`, `rh.aparc.annot`): the mean of the parcel's vertex coordinates, projected back to the sphere.

**Public pages read (2026-10-21), and nothing else:** adni.loni.usc.edu/data-samples/adni-data/ (access procedure); adni.loni.usc.edu/data-samples/adni-data/neuroimaging/pet/ (PET processing overview); the ADNI Data Use Agreement PDF (adni.loni.usc.edu/wp-content/themes/adni_2023/documents/ADNI_Data_Use_Agreement.pdf); the ADNI public data dictionary search (adni.loni.usc.edu/data-samples/data-dictionary-search/); humanconnectome.org/study/hcp-young-adult/data-use-terms, the WU-Minn HCP Open Access Data Use Terms page, and the 1200 Subjects Data Release page; the S1200 Reference Manual and its Appendix III (PDFs on humanconnectome.org). The LC atlas papers were read through NCBI E-utilities.

## 4. Methods

**Connectome construction** (outside the repository, settings fixed in `config/network.yaml`, key `connectome`). For each of the 100 subjects, with MRtrix3 3.0.4:

* Three-tissue response functions and multi-shell multi-tissue constrained spherical deconvolution.
* A five-tissue-type image from `aparc+aseg`. The LC atlas is warped from MNI space to the subject with `standard2acpc_dc` (nearest neighbour) and set as sub-cortical grey matter, so streamlines can end in it.
* 10 million streamlines, iFOD2 with anatomically constrained tractography, backtracking and dynamic seeding; maximum length 250 mm; SIFT2 weights.
* 84 nodes: 68 Desikan–Killiany cortical regions, 14 subcortical regions (thalamus, caudate, putamen, pallidum, hippocampus, amygdala, accumbens; both sides), and 2 LC nodes. Endpoints are assigned to the nearest node within 2 mm.
* Edge weight: the SIFT2-weighted streamline count times the subject's SIFT2 proportionality coefficient.

The group connectome is the mean across the 100 subjects, absent edges counting as zero (`python -m src.network.connectome`). Sensitivity 7 keeps only edges present in at least half the subjects.

**Model.** dx/dt = −βLx + αx on the group connectome, where L is the degree-normalised graph Laplacian (Raj et al. 2012). It is solved through the eigendecomposition of L. β and t enter only as τ = βt. Uniform growth (α) multiplies every region by the same factor and cannot change a correlation, so it is not fitted. A seed puts unit mass on its nodes, split equally (both LC nodes; both entorhinal regions).

**Fit.** For each seed, the predicted cortical pattern is computed at 400 log-spaced values of τ from 0.01 to 50. The correlation with the observed map is computed over the evaluation regions at each τ, and the maximum is kept.

**Spin null** (Alexander-Bloch et al. 2018; Váša et al. 2018). Each hemisphere's parcel centroids are rotated by a uniformly random rotation, the right hemisphere by the mirror image of the left's rotation, and parcels are reassigned to the rotated positions without duplicates (Hungarian assignment). The whole 68-region map is spun and the evaluation regions are read from the spun map. Every spun map is fitted exactly like the real one, including the maximisation over τ, so the seed with more freedom gets no free advantage. 10,000 spins, seed 20261021. Shuffling regions is never used.

**Outputs.** `results/tables/network_seed_comparison.csv` (every seed in every analysis, with p-values, crossover time and verdict), `results/tables/network_seed_ranking.csv`, and a figure of correlation against τ for each seed with the crossover marked.

## 5. Inference criteria

* One primary comparison (LC versus entorhinal), one-sided tests at 0.05, so there is no multiplicity correction. The verdict rules are in section 2.
* All p-values come from the spin null. A p-value is never below 1/10,001.
* Control seeds, the seed ranking and sensitivity analysis 4 and 9 are descriptive.

## 6. Sensitivity analyses (all pre-specified; `config/network.yaml`, key `sensitivity`)

Each changes one thing from the primary analysis and is run with its own 10,000 spins.

1. Unnormalised Laplacian.
2. Spearman instead of Pearson correlation.
3. Spins with nearest-neighbour reassignment (the original method) instead of Hungarian assignment.
4. Entorhinal regions kept in the evaluation (the first draft's design). It favours the entorhinal seed (section 2), so it is shown for transparency only.
5. Uniform LC projection: each LC node's tractography edges are replaced by equal weights to every cortical region of its own hemisphere, keeping its total weight. This is an assumption (a diffuse projection with no regional preference). It answers the unreliability of tractography to the LC.
6. LC mask from the 25% version of the atlas.
7. Consistency-thresholded connectome (edges present in at least half the subjects).
8. Tau SUVR without partial-volume correction.
9. Amyloid-positive MCI instead of amyloid-positive CN. A later stage, expected near or past the crossover time, so shown for transparency only.

The verdict is called robust only if analyses 1, 2, 3, 5, 6, 7 and 8 point the same way as the primary analysis.

## 7. Checks that must pass first

* **Validation of the spin test on a known null** (`tests/validation/test_network.py`, 6 tests, all passed on 2026-10-21 on the code registered here). Pairs of independent but spatially smooth maps on a synthetic sphere: the spin test gave false positives in 6.0% of 200 pairs at the 5% level (allowed up to 9.0%), naive region shuffling in 56.5%. The other tests check that the spin test still detects a real association, that spins are proper rotations and permutations, the Laplacian spectrum of a ring, the eigen-solution against the matrix exponential, and recovery of a known seed and diffusion time.
* The analysis refuses to run unless this validation passed on the current code and this document is frozen. The ADNI step refuses to run unless this document is frozen.

## 8. Limitations, stated in advance

* **Tractography to the LC is unreliable.** The LC is a few millimetres wide and the diffusion data have 1.25 mm voxels. If the LC's connection profile is mostly noise, the LC seed predicts a near-random pattern and the verdict will lean towards *against* or *inconclusive*. Sensitivity analysis 5 replaces it with a uniform projection, which is itself an assumption.
* **An entorhinal-like pattern is also what Braak's LC-first sequence predicts in the cortex** (section 2). The test can show that diffusion along the LC's connections does not explain the cortical map; it cannot rule out an LC origin with a single cortical entry point.
* Tau PET cannot measure tau in the LC itself, so the LC enters only as a seed.
* Group-level maps: individual patterns differ, and early-stage individual maps are too noisy to fit one by one.
* The connectome comes from young adults (HCP) and the tau maps from older adults (ADNI).
* Late stages carry little seed information, because high-eigenvalue modes decay fast. That is why the primary group is the earliest available, and why the crossover time is reported.
* Diffusion along connections may not be the mechanism of spread.
* Every name marked *to be confirmed on access* may need correcting, which I will record as a deviation.
* Every extracted value is unverified until I have checked it against its source.
