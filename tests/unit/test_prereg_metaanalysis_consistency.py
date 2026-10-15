"""prereg/metaanalysis.md (the registered text) must state what config/metaanalysis.yaml does.

Same idea as test_prereg_baserate_consistency.py: each case pairs a config value
with the registered sentence that fixes it.
"""
import re

import pytest

from src.common.config import load_config
from src.common.paths import CONFIG_DIR, REPO_ROOT

DOC = re.sub(r"\s+", " ", (REPO_ROOT / "prereg" / "metaanalysis.md").read_text())
CFG = load_config(CONFIG_DIR / "metaanalysis.yaml")

CASES = [
    ("pooling.method", "REML", "Random effects: REML, with Hartung–Knapp–Sidik–Jonkman confidence intervals"),
    ("pooling.knapp_hartung", True, "I am keeping Knapp–Hartung as originally planned."),
    ("pooling.rma_control.stepadj", 0.5, "Fisher-scoring step halving is used to help convergence"),
    ("pooling.min_studies", 2, "Pooling needs at least 2 studies."),
    ("pooling.combinable_measures", ["HR", "RR", "IRR"],
     "Hazard ratios, risk ratios and incidence rate ratios are pooled together. Odds ratios are pooled separately."),
    ("publication_bias.min_studies", 10, "only when there are at least 10 studies"),
    ("extraction.file", "data/extracted/meta_extraction.csv", "The file is `data/extracted/meta_extraction.csv`."),
    ("extraction.log_symmetry_tolerance", 0.15, "asymmetric on the log scale by more than 15% of the half-width"),
    ("triangulation.comparator_type", "active", "The confirmatory evidence is the set of **active-comparator contrasts**."),
    ("positive_control.drug_class", "anticholinergic", "A positive control:** anticholinergic use versus non-use"),
    ("positive_control.expected_direction", "harmful", "from k = 7, harmful, with an interval that excludes 1"),
    ("sensitivity.dl_estimator", True, "DerSimonian–Laird instead of REML."),
    ("sensitivity.no_knapp_hartung", True, "Plain z-based confidence intervals instead of Knapp–Hartung."),
    ("sensitivity.combine_or", True, "Odds ratios pooled with the hazard, risk and rate ratios"),
    ("sensitivity.include_e1_rows", True, "Rows excluded under E1 (Parkinson's populations) put back in."),
]


@pytest.mark.parametrize("key, value, phrase", CASES, ids=[c[0] for c in CASES])
def test_prereg_states_what_the_config_does(key, value, phrase):
    assert CFG.get(key) == value, f"config {key} = {CFG.get(key)!r}, but the pre-registration fixes {value!r}"
    assert re.sub(r"\s+", " ", phrase) in DOC, f"prereg/metaanalysis.md no longer says: {phrase!r}"


def test_contrast_predictions_match_the_registered_table():
    rows = {
        "lipophilic_vs_hydrophilic_beta_blocker": ("| lipophilic vs hydrophilic β-blockers | ADRB1 | < 1", "below_1"),
        "alpha1_blocker_vs_tamsulosin": ("| terazosin/doxazosin/alfuzosin vs tamsulosin | ADRA1A | < 1 |", "below_1"),
        "atomoxetine_vs_stimulant": ("| atomoxetine vs stimulants | SLC6A2 | < 1 |", "below_1"),
        "brain_penetrant_vs_nonpenetrant_dhp_ccb":
            ("| brain-penetrant DHP calcium-channel blockers vs amlodipine | CACNA1D | < 1 |", "below_1"),
        "z_drug_vs_dual_orexin_antagonist": ("| dual orexin antagonists vs Z-drugs | HCRTR2 | none;", "none"),
    }
    preds = CFG.get("contrast_predictions")
    assert sorted(preds) == sorted(rows)
    assert sorted(CFG.get("contrasts")) == sorted(rows)
    for contrast, (line, pred) in rows.items():
        assert line in DOC, line
        assert preds[contrast] == pred, contrast
    assert sum(p != "none" for p in preds.values()) == 4
    assert '"any of the 4 predicted contrasts AGAINST"' in DOC
