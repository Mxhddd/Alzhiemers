"""Module 4 validation on the RESEARCHER'S LC model (brief §5): analytic
fixed point vs hand derivation, non-negativity and conservation, solver
independence, and the zero-coupling degenerate case. These fail until the
model (src/dynamics/lc_model.py), its parameters and the hand derivation exist."""

from pathlib import Path

import numpy as np
import pytest

from src.common.config import load_config
from src.common.paths import CONFIG_DIR, REPO_ROOT
from src.dynamics import lc_model
from src.dynamics.fixedpoints import find_fixed_points, grid_guesses
from src.dynamics.params import ParameterSet
from src.dynamics.simulate import simulate, solver_independence

pytestmark = pytest.mark.validation
TARGET = REPO_ROOT / "tests" / "validation" / "targets" / "dynamics_analytic_fixed_point.yaml"


def _setup():
    params = ParameterSet.from_yaml(CONFIG_DIR / "dynamics_params.yaml")
    if not params.params:
        raise AssertionError("config/dynamics_params.yaml has no parameters (researcher to write)")
    cfg = load_config(CONFIG_DIR / "dynamics.yaml")
    return lc_model.build_model(), params, cfg


def test_fixed_point_matches_hand_derivation():
    model, params, cfg = _setup()
    t = load_config(Path(TARGET))
    p = params.values(t.get("parameters"))
    fps = find_fixed_points(model, p, grid_guesses(cfg.get("fixed_points.bounds"), model.state_names,
                                                   cfg.get("fixed_points.grid_per_dim")))
    want = np.array([t.get("fixed_point")[s] for s in model.state_names], float)
    tol = t.get("tolerance")
    match = [f for f in fps if np.max(np.abs(f.y - want)) <= tol]
    assert len(match) == 1, f"hand fixed point {want} not found among {[f.y for f in fps]}"
    stab = t.get("stability")
    got = "stable" if match[0].stable else ("saddle" if match[0].kind == "saddle" else "unstable")
    assert got == stab, f"stability: numerical {match[0].kind} vs hand {stab}"


def test_nonnegativity_conservation_and_solver_independence():
    model, params, cfg = _setup()
    p = params.values()
    t0, t1 = cfg.get("simulation.t_span")
    t_eval = np.linspace(t0, t1, cfg.get("simulation.n_out"))
    simulate(model, p, lc_model.initial_state(p), (t0, t1), t_eval, method=cfg.get("simulation.method"),
             rtol=cfg.get("simulation.rtol"), atol=cfg.get("simulation.atol"))
    rep = solver_independence(model, p, lc_model.initial_state(p), (t0, t1), t_eval,
                              rtol=cfg.get("simulation.rtol"), atol=cfg.get("simulation.atol"))
    assert rep["converging"] and rep["cross_solver_tight"] <= cfg.get("simulation.max_solver_disagreement"), rep


def test_zero_coupling_reduces_to_known_exponentials():
    model, params, cfg = _setup()
    p = params.values({k: 0.0 for k in lc_model.coupling_params()})
    t0, t1 = cfg.get("simulation.t_span")
    t_eval = np.linspace(t0, t1, cfg.get("simulation.n_out"))
    y0 = lc_model.initial_state(p)
    tr = simulate(model, p, y0, (t0, t1), t_eval, method="Radau", rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(tr.y, lc_model.uncoupled_solution(t_eval, p, y0), rtol=1e-6, atol=1e-9)
