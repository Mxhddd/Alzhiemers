"""prereg/network.md must state exactly what config/network.yaml does.

Each case pairs a config value with the phrase in the pre-registration that
fixes it. A config edit not mirrored in the document (or the reverse) fails
here, so the document handed over for freezing cannot drift from the code.
"""
import re

import pytest

from src.common.config import load_config
from src.common.paths import CONFIG_DIR, REPO_ROOT

DOC = re.sub(r"\s+", " ", (REPO_ROOT / "prereg" / "network.md").read_text(encoding="utf-8"))
CFG = load_config(CONFIG_DIR / "network.yaml")

CASES = [
    # (config key, expected config value, phrase that must appear in the prereg)
    ("seed", 20261021, "10,000 spins, seed 20261021"),
    ("alpha", 0.05, "one-sided tests at 0.05"),
    ("alpha", 0.05, "p(LC better) < 0.05, p_LC < 0.05"),
    ("model.laplacian", "normalized", "L is the degree-normalised graph Laplacian"),
    ("model.tau_min", 0.01, "400 log-spaced values of τ from 0.01 to 50"),
    ("model.tau_max", 50.0, "400 log-spaced values of τ from 0.01 to 50"),
    ("model.n_tau", 400, "400 log-spaced values of τ from 0.01 to 50"),
    ("model.crossover_threshold", 0.95, "correlate at 0.95 or more"),
    ("spin.n_perm", 10000, "Every spin p-value is one-sided, uses 10,000 spins"),
    ("spin.method", "hungarian", "without duplicates (Hungarian assignment)"),
    ("evaluation.correlation", "pearson", "the correlation (Pearson)"),
    ("evaluation.exclude", ["ctx-lh-entorhinal", "ctx-rh-entorhinal"],
     "the 68 Desikan–Killiany regions minus the two entorhinal regions"),
    ("seeds.lc", ["LC-lh", "LC-rh"], "(both LC nodes; both entorhinal regions)"),
    ("seeds.entorhinal", ["ctx-lh-entorhinal", "ctx-rh-entorhinal"], "(both LC nodes; both entorhinal regions)"),
    ("seeds.controls", {"V1": ["ctx-lh-pericalcarine", "ctx-rh-pericalcarine"],
                        "precuneus": ["ctx-lh-precuneus", "ctx-rh-precuneus"]},
     "the pericalcarine cortex (V1, late in Braak staging) and the precuneus"),
    ("robustness_set", ["S1_unnormalized_laplacian", "S2_spearman", "S3_nearest_spin", "S5_lc_uniform",
                        "S6_lc_mask_25", "S7_consistency_threshold", "S8_no_pvc"],
     "sensitivity analyses 1, 2, 3, 5, 6, 7 and 8 (section 6) all point the same way"),
    ("connectome.release", "HCP-YA S1200 3T", "1200 Subjects Release (S1200), 3T"),
    ("connectome.subject_group", "100 Unrelated Subjects", "\"100 Unrelated Subjects\" (100 people)"),
    ("connectome.n_subjects", 100, "the mean across the 100 subjects"),
    ("connectome.packages", ["Diffusion Preprocessed", "Structural Preprocessed"],
     "*Diffusion Preprocessed* (`T1w/Diffusion/data.nii.gz`"),
    ("connectome.software", "MRtrix3 3.0.4", "with MRtrix3 3.0.4"),
    ("connectome.fod", "msmt_csd", "multi-shell multi-tissue constrained spherical deconvolution"),
    ("connectome.algorithm", "iFOD2_ACT", "iFOD2 with anatomically constrained tractography"),
    ("connectome.n_streamlines", 10_000_000, "10 million streamlines"),
    ("connectome.max_length_mm", 250, "maximum length 250 mm"),
    ("connectome.radial_search_mm", 2, "nearest node within 2 mm"),
    ("connectome.weight", "sift2_mu", "SIFT2-weighted streamline count times the subject's SIFT2 proportionality"),
    ("connectome.n_nodes", 84, "84 nodes: 68 Desikan–Killiany cortical regions, 14 subcortical regions"),
    ("connectome.lc_atlas", "ye2021_7t", "Ye et al. 2021 (*NeuroImage* 225:117487), a 7T probabilistic LC atlas"),
    ("connectome.lc_threshold_primary", 5, "Primary: the 5% (liberal) version"),
    ("connectome.lc_threshold_sensitivity", 25, "sensitivity 6: the 25% (stringent) version"),
    ("connectome.lc_in_5tt", "sgm", "set as sub-cortical grey matter, so streamlines can end in it"),
    ("connectome.consistency_min_fraction", 0.5, "only edges present in at least half the subjects"),
    ("adni.tau_tracer", "FTP", "flortaucipir (assumed `FTP`)"),
    ("adni.qc_pass", [2], "visual QC pass (`qc_flag` = 2)"),
    ("adni.amyloid_positive_value", 1, "`AMYLOID_STATUS` for positive (assumed 1)"),
    ("adni.amyloid_window_days", 365, "nearest in time to the tau scan, within 365 days"),
    ("adni.diagnosis_window_days", 180, "nearest in time, within 180 days"),
    ("adni.diagnosis_codes", {"cu": 1, "mci": 2}, "`DIAGNOSIS` (1 CN, 2 MCI, 3 dementia)"),
    ("adni.female_code", 2, "`PTGENDER` (1 male, 2 female)"),
    ("adni.w_score_covariates", ["age", "sex"], "a linear regression of SUVR on age and sex in the reference group"),
    ("adni.min_group_n", 30, "need at least 30 people"),
    ("adni.suvr_column", "CTX_{hemi}_{region}_SUVR", "`CTX_LH_<REGION>_SUVR`, `CTX_RH_<REGION>_SUVR`"),
    ("adni.columns.amyloid_status", "AMYLOID_STATUS", "`AMYLOID_STATUS` (the lab's own threshold"),
    ("adni.columns.diagnosis", "DIAGNOSIS", "`DIAGNOSIS` (1 CN, 2 MCI, 3 dementia), `EXAMDATE`"),
    ("adni.columns.dx_date", "EXAMDATE", "`DIAGNOSIS` (1 CN, 2 MCI, 3 dementia), `EXAMDATE`"),
    ("adni.columns.birth_year", "PTDOBYY", "`PTDOBYY` (year of birth)"),
    ("adni.columns.qc", "qc_flag", "`qc_flag` (2 pass, 1 partial, 0 fail)"),
    ("adni.columns.tau_date", "SCANDATE", "`SCANDATE`; `TRACER`"),
    ("adni.columns.id", "RID", "that `RID` identifies participants in every table"),
    ("adni.files.tau_pvc", "data/raw/adni/UCBERKELEY_TAUPVC_6MM.csv", "(`UCBERKELEY_TAUPVC_6MM`)"),
    ("adni.files.tau_nopvc", "data/raw/adni/UCBERKELEY_TAU_6MM.csv", "(`UCBERKELEY_TAU_6MM`)"),
    ("adni.files.amyloid", "data/raw/adni/UCBERKELEY_AMY_6MM.csv", "(`UCBERKELEY_AMY_6MM`)"),
    ("adni.files.diagnosis", "data/raw/adni/DXSUM.csv", "Diagnostic Summary (`DXSUM`)"),
    ("adni.files.demographics", "data/raw/adni/PTDEMOG.csv", "Participant Demographics (`PTDEMOG`)"),
]

SENSITIVITY_PHRASES = {
    "S1_unnormalized_laplacian": ({"laplacian": "unnormalized"}, "1. Unnormalised Laplacian."),
    "S2_spearman": ({"correlation": "spearman"}, "2. Spearman instead of Pearson correlation."),
    "S3_nearest_spin": ({"spin_method": "nearest"}, "3. Spins with nearest-neighbour reassignment"),
    "S4_entorhinal_evaluated": ({"exclude": []}, "4. Entorhinal regions kept in the evaluation"),
    "S5_lc_uniform": ({"lc_projection": "uniform_ipsilateral_cortex"},
                      "equal weights to every cortical region of its own hemisphere, keeping its total weight"),
    "S6_lc_mask_25": ({"connectome": "connectome_lc25"}, "6. LC mask from the 25% version of the atlas."),
    "S7_consistency_threshold": ({"connectome": "connectome_consistent"}, "7. Consistency-thresholded connectome"),
    "S8_no_pvc": ({"tau": "tau_nopvc"}, "8. Tau SUVR without partial-volume correction."),
    "S9_apos_mci": ({"tau": "tau_mci"}, "9. Amyloid-positive MCI instead of amyloid-positive CN."),
}


@pytest.mark.parametrize("key, value, phrase", CASES, ids=[f"{c[0]}-{i}" for i, c in enumerate(CASES)])
def test_prereg_states_what_the_config_does(key, value, phrase):
    assert CFG.get(key) == value, f"config {key} = {CFG.get(key)!r}, but the pre-registration fixes {value!r}"
    assert re.sub(r"\s+", " ", phrase) in DOC, f"prereg/network.md no longer says: {phrase!r}"


def test_sensitivity_analyses_are_exactly_the_nine_in_the_document():
    sens = CFG.get("sensitivity")
    assert list(sens) == list(SENSITIVITY_PHRASES)
    for name, (change, phrase) in SENSITIVITY_PHRASES.items():
        assert sens[name] == change, name
        assert phrase in DOC, phrase


def test_data_paths_point_at_the_group_level_files():
    for key, part in (("data.tau", "apos_cu_pvc"), ("data.tau_nopvc", "apos_cu_nopvc"), ("data.tau_mci", "apos_mci_pvc"),
                      ("data.connectome", "lc5"), ("data.connectome_lc25", "lc25"),
                      ("data.connectome_consistent", "consistency50")):
        assert part in CFG.get(key)


def test_document_has_no_freeze_block_placeholders_or_draft_status():
    assert "DO NOT FREEZE" not in DOC
    assert "{{" not in DOC and "}}" not in DOC
    assert "Status:" not in DOC
    assert "Nothing from either has been requested, downloaded, opened or seen" in DOC
