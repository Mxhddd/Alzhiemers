"""prereg/baserate.md must state exactly what config/baserate.yaml does.

Each case pairs a config value with the sentence in the pre-registration that
fixes it. A config edit that is not mirrored in the document (or the reverse)
fails here, so the document handed over for freezing cannot drift from the
code that will run.
"""
import math
import re

import pytest

from src.common.config import load_config
from src.common.paths import CONFIG_DIR, REPO_ROOT

DOC = re.sub(r"\s+", " ", (REPO_ROOT / "prereg" / "baserate.md").read_text())
CFG = load_config(CONFIG_DIR / "baserate.yaml")

CASES = [
    # (config key, expected config value, phrase that must appear in the prereg)
    ("amyloid.gating_is_a_test", False, "Module 2 does not test amyloid gating"),
    ("variants.lc_first", "A_gated_age", "A, LC first (age-adjusted)"),
    ("variants.c_hits", 2, "after m = 2 independent hits"),
    ("variants.c_hits_sensitivity", 3, "Variant C with 3 hits"),
    ("cross_validation.k", 5, "5-fold cross-validation by whole age band"),
    ("cross_validation.margin", 2.0, "by more than 2 held-out log-likelihood units"),
    ("cross_validation.n_random_starts", 4, "full-data MLE plus 4 random starts"),
    ("mle.n_starts", 20, "20 random starts"),
    ("mle.max_newton_decrement", 1.0e-3, "Newton decrement is < 10⁻³"),
    ("mle.bounds.log_rate", [-20.0, 2.0], "per-year rates e⁻²⁰ to e²"),
    ("mle.bounds.slope", [-2.0, 2.0], "at most an e²-fold change per decade"),
    ("mle.bounds.log_hr", [-5.0, 5.0], "amyloid HR e⁻⁵ to e⁵"),
    ("bayes.chains", 4, "4 chains × 1,000 draws after 1,000 warm-up"),
    ("bayes.draws", 1000, "4 chains × 1,000 draws"),
    ("bayes.tune", 1000, "after 1,000 warm-up"),
    ("bayes.init_overdispersion", 2.0, "N(MLE, 4 × Laplace covariance)"),
    ("bayes.log_rate_prior", {"mu": -5.0, "sigma": 3.0}, "log-rate N(−5, 3)"),
    ("bayes.slope_prior", {"mu": 0.0, "sigma": 1.0}, "log-slope per decade N(0, 1)"),
    ("bayes.log_hr_prior", {"mu": 0.0, "sigma": 1.0}, "log HR_A N(0, 1)"),
    ("bayes.target_accept", 0.9, "Target acceptance probability 0.9"),
    ("bayes.max_rhat", 1.01, "max R-hat ≤ 1.01"),
    ("bayes.max_divergences", 0, "zero divergences"),
    ("solver.step_years", 0.25, "step 0.25 years"),
    ("solver.max_prob_error", 1.0e-3, "h²Φ(θ̂) must be ≤ 10⁻³"),
    ("solver.max_loglik_error", 0.1, "|ΔLL| < 0.1"),
    ("solver.start_age", 1, "Everyone is at stage 0 at age 1"),
    ("solver.age_centre", 50, "exp(α_k + β_k (t − 50)/10)"),
    ("solver.age_scale", 10, "exp(α_k + β_k (t − 50)/10)"),
    ("solver.evaluate_at", "midpoint", "stage probabilities at the band midpoint"),
    ("solver.evaluate_at_sensitivity", "band_average", "Band-average instead of midpoint evaluation"),
    ("data.dementia_threshold_stage", "V", "entering Braak stage **V**"),
    ("data.dementia_threshold_stage_sensitivity", "IV", "Dementia threshold stage IV instead of V"),
    ("data.open_band_evaluation_age", 90, "85+ is evaluated at age 90"),
    ("data.open_band_evaluation_age_sensitivity", 87.5, "evaluated at 87.5 instead of 90"),
    ("data.braak_counts", "data/extracted/braak2011_stage_counts.csv", "`data/extracted/braak2011_stage_counts.csv`"),
    ("data.amyloid_prevalence", "data/extracted/braak2011_sdc2_amyloid_by_age.csv",
     "`data/extracted/braak2011_sdc2_amyloid_by_age.csv`"),
    ("data.amyloid_positive_min_phase", 1, "brains in phase 1–4 / all brains in the band"),
    ("data.amyloid_prevalence_sensitivity", "data/extracted/amyloid_prevalence_jansen2015.csv",
     "`data/extracted/amyloid_prevalence_jansen2015.csv`"),
    ("data.dementia_incidence", "data/extracted/dementia_incidence_cfas.csv", "`data/extracted/dementia_incidence_cfas.csv`"),
    ("falsification.n_draws", 200, "Predictions are computed for 200 parameter draws"),
    ("gamma_check.expected_log_hr", 0.888, "log HR 0.888 (0.322–1.449)"),
    ("gamma_check.expected_log_hr_ci_lower", 0.322, "log HR 0.888 (0.322–1.449)"),
    ("gamma_check.expected_log_hr_ci_upper", 1.449, "log HR 0.888 (0.322–1.449)"),
    ("gamma_check.uninformative_log_width", 4.605, "log width > 4.605"),
]


@pytest.mark.parametrize("key, value, phrase", CASES, ids=[c[0] for c in CASES])
def test_prereg_states_what_the_config_does(key, value, phrase):
    assert CFG.get(key) == value, f"config {key} = {CFG.get(key)!r}, but the pre-registration fixes {value!r}"
    assert re.sub(r"\s+", " ", phrase) in DOC, f"prereg/baserate.md no longer says: {phrase!r}"


def test_falsification_tolerance_is_the_stated_factor_of_two():
    assert math.isclose(CFG.get("falsification.max_abs_log_ratio"), math.log(2), abs_tol=1e-4)
    assert "widened by a factor of 2 on both sides" in DOC


def test_fitted_variants_are_the_three_in_the_document():
    assert list(CFG.get("variants.fit")) == ["A_gated_age", "B_independent", "C_multihit"]
    for v in ("**A, LC first (age-adjusted):**", "**B, independent:**", "**C, multi-hit:**"):
        assert v in DOC


def test_document_carries_no_freeze_block_and_no_pending_inputs():
    assert "DO NOT FREEZE" not in DOC
    assert "pending" not in DOC.lower()
