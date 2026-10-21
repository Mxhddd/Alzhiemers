"""prereg/dynamics.md must state exactly what config/dynamics.yaml and
config/dynamics_params.yaml do.

Each case pairs a config value with the phrase in the pre-registration that
fixes it. A config edit not mirrored in the document (or the reverse) fails here.
"""
import re

import pytest

from src.common.config import load_config
from src.common.paths import CONFIG_DIR, REPO_ROOT
from src.dynamics.params import ParameterSet

DOC = re.sub(r"\s+", " ", (REPO_ROOT / "prereg" / "dynamics.md").read_text(encoding="utf-8"))
CFG = load_config(CONFIG_DIR / "dynamics.yaml")
PS = ParameterSet.from_yaml(CONFIG_DIR / "dynamics_params.yaml")

CASES = [
    ("seed", 20261019, "seed 20261019"),
    ("pilot.seed", 20261021, "its own seed 20261021"),
    ("model.readout", "steady_state_from_healthy_start", "The steady-state tau burden T\\* reached from the healthy start"),
    ("model.healthy_start", {"E": 0.0, "Y": 0.0, "R": 1.0, "T": 0.0}, "The healthy start is (E, Y, R, T) = (0, 0, 1, 0)"),
    ("simulation.max_solver_disagreement", 1.0e-3, "Radau and BDF must agree within a relative 10⁻³ with tolerances tightened 10×"),
    ("simulation.method", "Radau", "Radau and BDF"),
    ("simulation.readout_check_n", 50, "For the first 50 Latin-hypercube sets"),
    ("simulation.readout_check_max_rel_diff", 1.0e-3, "reproduce the analytic T\\* within a relative 10⁻³"),
    ("fixed_points.r_grid", {"log10_min": -9.0, "log10_switch": -2.0, "n_log": 400, "n_lin": 3600},
     "400 log-spaced points on [10⁻⁹, 10⁻²] plus 3,600 linear points on [10⁻², 1]"),
    ("bifurcation.parameter", "hA", "in the second-hit parameter hA, that is, is there a tipping point"),
    ("bifurcation.range", [0.0, 1.0], "the set of hA in [0, 1] with two stable fixed points"),
    ("bifurcation.method", "exact_window", "computed exactly from Φ(R)"),
    ("bifurcation.n_values", 201, "The 201-point scan is kept as a cross-check"),
    ("dose_response.reference_dose", 1.0, "T\\*(drug at 1 IC50) / T\\*(no drug)"),
    ("dose_response.protective_ratio", 0.95, "*protective* below 0.95"),
    ("dose_response.harmful_ratio", 1.05, "*harmful* above 1.05"),
    ("dose_response.dose_grid", {"include_zero": True, "min": 0.01, "max": 10.0, "n_log": 41},
     "over 0 plus 41 log-spaced doses from 0.01 to 10 IC50"),
    ("dose_response.min_relative_rise", 0.05, "both ends exceed it by more than 5%"),
    ("sensitivity.sobol_n_base", 1024, "Saltelli sampling, N = 1,024"),
    ("sensitivity.lhs_n_samples", 1000, "a Latin hypercube of 1,000 parameter sets"),
    ("robustness.fragile_upper", 0.10, "**fragile** if the upper bound is below 10%"),
    ("robustness.robust_lower", 0.50, "**robust** if the lower bound is at least 50%"),
    ("robustness.ci_level", 0.95, "Wilson 95% CI"),
    ("sensitivity_analyses.SA1_central_half", {"shrink_to": 0.5}, "every range shrunk to its central half"),
    ("sensitivity_analyses.SA2_no_alpha1_clearance", {"fix": {"eps": 0.0}}, "**SA2:** ε = 0"),
    ("sensitivity_analyses.SA3_eta_linear", {"range": {"eta": [0.0, 10.0]}, "scale": {"eta": "linear"}},
     "**SA3:** η uniform on [0, 10]"),
    ("sensitivity_analyses.SA4_K_draft_range", {"range": {"K": [0.1, 10.0]}}, "**SA4:** K log-uniform on [0.1, 10]"),
    ("sensitivity_analyses.SA5_thresholds", [[0.90, 1.10], [0.99, 1.01]], "0.90/1.10 and 0.99/1.01"),
    ("pilot.morris_trajectories", 50, "750 model evaluations"),
]


@pytest.mark.parametrize("key, value, phrase", CASES, ids=[c[0] for c in CASES])
def test_prereg_states_what_the_config_does(key, value, phrase):
    assert CFG.get(key) == value, f"config {key} = {CFG.get(key)!r}, but the pre-registration fixes {value!r}"
    assert re.sub(r"\s+", " ", phrase) in DOC, f"prereg/dynamics.md no longer says: {phrase!r}"


def test_morris_evaluation_count_matches_the_document():
    n = CFG.get("pilot.morris_trajectories") * (len(PS.params) + 1)
    assert n == 750 and "750 model evaluations" in DOC


SYMBOL = {"kappa": "κ", "u0": "u0", "phi": "φ", "v": "v", "m": "m", "r_rep": "r_rep", "r_dam": "r_dam", "K": "K",
          "q": "q", "sigma": "σ", "eta": "η", "eps": "ε", "hA": "hA", "c0": "c0"}


def _fmt(x):
    return f"{x:g}"


@pytest.mark.parametrize("name", sorted(SYMBOL))
def test_parameter_table_matches_the_parameter_file(name):
    p = PS.params[name]
    lo, hi = p.range
    rng = f"{_fmt(lo)}–{_fmt(hi)}" + (", linear" if p.scale == "linear" else "")
    row = re.search(rf"\| {re.escape(SYMBOL[name])} \| [^|]+ \| ([^|]+) \|", DOC)
    assert row, f"no row for {name} in the parameter table"
    assert row.group(1).strip().startswith(rng), f"{name}: document says {row.group(1)!r}, config {rng!r}"
    assert p.confidence == "weak" and p.source.startswith("assumption:")


def test_all_parameters_are_assumptions_as_the_document_says():
    assert len(PS.params) == 14 and "All 14 parameters are therefore assumptions" in DOC
    assert "No published value sets any parameter" in DOC


def test_validation_fixed_point_is_the_one_in_the_document():
    t = load_config(REPO_ROOT / "tests" / "validation" / "targets" / "dynamics_analytic_fixed_point.yaml")
    assert t.get("fixed_point") == {"E": 0.5, "Y": 0.125, "R": 0.5, "T": 12.0}
    assert "(E, Y, R, T) = (1/2, 1/8, 1/2, 12)" in DOC and t.get("stability") == "stable"
    assert "derived by Claude, machine-checked" in DOC


def test_document_carries_no_freeze_block_no_pending_inputs_and_no_unfilled_placeholders():
    assert "DO NOT FREEZE" not in DOC
    assert "pending" not in DOC.lower()
    assert "{{" not in DOC and "}}" not in DOC
