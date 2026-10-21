"""prereg/mediation.md must state exactly what config/mediation.yaml does.

Each case pairs a config value with the phrase in the pre-registration that
fixes it, so the document handed over for freezing cannot drift from the code
that will run (model: test_prereg_baserate_consistency.py).
"""
import re

import pytest

from src.common.config import load_config
from src.common.paths import CONFIG_DIR, REPO_ROOT

DOC = re.sub(r"\s+", " ", (REPO_ROOT / "prereg" / "mediation.md").read_text())
CFG = load_config(CONFIG_DIR / "mediation.yaml")

CASES = [
    ("data.cohorts", ["ROS", "MAP"], "MARS, Clinical Core and Latino Core participants are excluded"),
    ("data.columns.id", "projid", "`projid`, `study`, `fu_year`"),
    ("data.columns.visit_age", "age_bl_plus_fu_year", "computed as `age_bl` + `fu_year`"),
    ("data.codebook_confirmed", False, "(`data.codebook_confirmed`)"),
    ("variables.outcome", "braaksc", "| Outcome | `braaksc` |"),
    ("variables.outcome_sensitivity", "b_score", "Braak stage in 4 levels (`b_score`) instead of 7"),
    ("variables.amyloid", "amylsqrt_est_8reg", "| Exposure H1/H2 | `amylsqrt_est_8reg` |"),
    ("variables.amyloid_sensitivity", "a_score_st4", "Thal-based A score (`a_score_st4`) instead of amyloid density"),
    ("variables.apoe4_carrier_codes", [24, 34, 44], "codes 24, 34, 44 = ε4 carrier"),
    ("variables.normal_cognition_code", 1, "(`dcfdx` = 1)"),
    ("exposure_timing.primary", "earliest", "earliest non-missing value"),
    ("exposure_timing.sensitivity", ["mean", "last"], "the mean of all visits, and as the last value before death"),
    ("exposure_timing.min_lag_years", 5, "at least 5 years before death"),
    ("model.link", "logit", "proportional-odds (cumulative logit) regression"),
    ("model.covariates", ["age_death", "msex", "educ", "apoe4", "map"],
     "age at death, sex, education, APOE ε4 carrier and cohort"),
    ("model.entry_split", 1, "Braak stage I or above versus 0"),
    ("model.ci_level", 0.95, "Unadjusted 95% intervals"),
    ("model.equivalence_log_or", 0.10536, "|log OR| < 0.10536"),
    ("model.uninformative_log_width", 1.386, "log width above 1.386"),
    ("families.amyloid", ["H1", "H2"], "Family A is H1 and H2 (2 tests, 97.5% intervals)"),
    ("families.social_sensory", ["soc_net", "social_isolation", "late_life_soc_act", "vision"],
     "Family B is H3 for the four exposures (4 tests, 98.75% intervals)"),
    ("missing.primary", "complete_case_logged", "complete cases per analysis"),
    ("missing.mi_threshold_share", 0.05, "more than 5% of the eligible participants"),
    ("missing.mi_imputations", 20, "20 imputations"),
    ("missing.mi_cycles", 10, "10 cycles"),
    ("autopsy_weights.truncate_quantile", 0.99, "truncated at the 99th percentile"),
    ("autopsy_weights.n_boot", 500, "500 bootstrap replicates"),
    ("lc.variables", [], "(`lc.variables`)"),
    ("lc.priority", ["lc_neuronal_density", "lc_neuronal_loss_rating", "lc_tangle_density"],
     "LC neuronal density, LC neuronal loss rating, LC tangle density"),
    ("lc.min_n", 100, "at least 100 autopsied participants with the LC measure"),
    ("mediation.contrast_quantiles", [0.25, 0.75], "from its 25th to its 75th percentile"),
    ("mediation.n_boot", 2000, "2000 case-resampling bootstrap replicates"),
    ("mediation.family_size", 4, "Family C: 4 tests"),
    ("mediation.rho_min", -0.95, "from −0.95 to 0.95 in steps of 0.01"),
    ("mediation.rho_max", 0.95, "from −0.95 to 0.95 in steps of 0.01"),
    ("mediation.rho_step", 0.01, "from −0.95 to 0.95 in steps of 0.01"),
    ("mediation.fragile_abs_rho", 0.1, "below 0.1 in absolute value is reported as **fragile**"),
    ("sensitivity.robustness_set", ["S1", "S2", "S3", "S5", "S6"], "S1, S2, S3, S5 and S6"),
    ("killer2.predicted_sign", 1, "the natural indirect effect (NIE) through the LC measure on expected Braak stage is positive"),
]


@pytest.mark.parametrize("key, value, phrase", CASES, ids=[c[0] for c in CASES])
def test_prereg_states_what_the_config_does(key, value, phrase):
    assert CFG.get(key) == value, f"config {key} = {CFG.get(key)!r}, but the pre-registration fixes {value!r}"
    assert re.sub(r"\s+", " ", phrase) in DOC, f"prereg/mediation.md no longer says: {phrase!r}"


def test_exposures_and_directions_match_the_document():
    ex = {e["name"]: (e["risk_direction"], e["cohorts"]) for e in CFG.get("exposures")}
    assert ex == {"soc_net": (-1, ["ROS", "MAP"]), "social_isolation": (1, ["MAP"]),
                  "late_life_soc_act": (-1, ["MAP"]), "vision": (1, ["ROS", "MAP"])}
    assert ("smaller `soc_net`, larger `social_isolation`, smaller `late_life_soc_act`, larger `vision` code"
            in DOC)
    assert "MAP only" in DOC


def test_seed_and_predictors_of_autopsy():
    assert CFG.get("seed") == 20261021
    assert CFG.get("autopsy_weights.predictors") == ["age_bl", "msex", "educ", "map", "cogn_baseline"]
    assert "baseline age, sex, education, cohort, baseline global cognition" in DOC


def test_document_states_the_limits_of_the_dataset_and_carries_no_freeze_block():
    assert "these data cannot test the LC → cortex step itself" in DOC
    assert "This test can hardly fail, and it does not discriminate" in DOC
    assert "DO NOT FREEZE" not in DOC and "{{" not in DOC
