"""Module 4 entry point: python -m src.dynamics.run_dynamics  (or `make dynamics`)

Refuses to start unless Module 4's validation passed on the current code
(including the researcher's analytic fixed point) and prereg/dynamics.md is
frozen. Analyses (brief §5): baseline trajectory with solver-independence
check; fixed points and stability; bifurcation over the second-hit parameter;
drug simulation and dose-response (U-shape test); Sobol sensitivity; Latin
hypercube robustness — the fraction of plausible parameter space in which
the U-shape appears, and in which BOTH drugs are protective.
Weak parameters are swept, never reported as point values.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from src.common.config import load_config
from src.common.paths import CONFIG_DIR, TABLES_DIR
from src.common.plotting import save_figure
from src.common.prereg import require_preregistration
from src.common.runlog import start_run
from src.common.validation import require_validated
from src.dynamics import figures, lc_model
from src.dynamics.bifurcation import bistable_range, folds, scan
from src.dynamics.fixedpoints import find_fixed_points, grid_guesses
from src.dynamics.params import ParameterSet
from src.dynamics.sensitivity import is_u_shaped, lhs_robustness, sobol_indices, wilson_interval
from src.dynamics.simulate import simulate, solver_independence

CONFIG = CONFIG_DIR / "dynamics.yaml"
PARAMS = CONFIG_DIR / "dynamics_params.yaml"


def main() -> int:
    require_validated("dynamics")
    prereg = require_preregistration("dynamics")
    cfg = load_config(CONFIG)
    params = ParameterSet.from_yaml(PARAMS)
    if not params.params:
        raise ValueError(f"{PARAMS} is empty: the researcher writes the parameters")
    model = lc_model.build_model()
    drugs = lc_model.drugs()
    for need in ("atomoxetine", "prazosin"):
        if need not in drugs:
            raise KeyError(f"lc_model.drugs() must define '{need}' (brief §5 question)")
    seed = cfg.get("seed")
    t0, t1 = cfg.get("simulation.t_span")
    t_eval = np.linspace(t0, t1, cfg.get("simulation.n_out"))
    sim_kw = dict(method=cfg.get("simulation.method"), rtol=cfg.get("simulation.rtol"), atol=cfg.get("simulation.atol"))
    doses = cfg.get("dose_response.doses")
    rise = cfg.get("dose_response.min_relative_rise")

    def outcome_for(p):
        return lc_model.outcome(simulate(model, p, lc_model.initial_state(p), (t0, t1), t_eval, **sim_kw))

    def dose_curve(p, drug):
        d = np.asarray(doses[drug], float)
        if d[0] != 0:
            raise ValueError(f"dose grid for {drug} must start at 0 (no drug)")
        return d, np.array([outcome_for(drugs[drug].apply(p, x)) for x in d])

    with start_run(module="dynamics", config_path=CONFIG, seed=seed,
                   packages=["numpy", "scipy", "SALib"], preregistration=prereg) as run:
        TABLES_DIR.mkdir(parents=True, exist_ok=True)
        run.info(f"parameters: {len(params.params)} ({len(params.weak())} weak, swept not reported)")
        rows = [{"parameter": k, "value": (p.value if p.confidence != "weak" else np.nan), "units": p.units,
                 "source": p.source, "confidence": p.confidence,
                 "range_low": p.range[0] if p.range else np.nan, "range_high": p.range[1] if p.range else np.nan}
                for k, p in params.params.items()]
        pd.DataFrame(rows).to_csv(TABLES_DIR / "dynamics_parameters.csv", index=False)
        run.record_output(TABLES_DIR / "dynamics_parameters.csv")
        p0 = params.values()

        y0 = lc_model.initial_state(p0)
        traj = simulate(model, p0, y0, (t0, t1), t_eval, **sim_kw)
        rep = solver_independence(model, p0, y0, (t0, t1), t_eval, rtol=sim_kw["rtol"], atol=sim_kw["atol"])
        run.info(f"solver independence: {rep}")
        if not (rep["converging"] and rep["cross_solver_tight"] <= cfg.get("simulation.max_solver_disagreement")):
            raise RuntimeError(f"results depend on the numerics: {rep}")
        save_figure(figures.trajectory(traj, model.state_names), "dynamics_baseline", run)

        guesses = grid_guesses(cfg.get("fixed_points.bounds"), model.state_names, cfg.get("fixed_points.grid_per_dim"))
        fps = find_fixed_points(model, p0, guesses)
        fp_df = pd.DataFrame([{**dict(zip(model.state_names, f.y)), "stable": f.stable, "kind": f.kind,
                               "eigenvalues": ";".join(f"{e:.4g}" for e in f.eigenvalues)} for f in fps])
        fp_df.to_csv(TABLES_DIR / "dynamics_fixed_points.csv", index=False)
        run.record_output(TABLES_DIR / "dynamics_fixed_points.csv")
        run.info(f"{len(fps)} fixed points, {sum(f.stable for f in fps)} stable")

        bp = cfg.get("bifurcation.parameter")
        lo, hi = cfg.get("bifurcation.range")
        rows_b = scan(model, p0, bp, np.linspace(lo, hi, cfg.get("bifurcation.n_values")), guesses)
        br = bistable_range(rows_b)
        run.info(f"bifurcation over {bp}: folds {folds(rows_b)}; bistable range {br if br else 'NONE (not bistable)'}")
        pd.DataFrame([{k: r[k] for k in ("value", "n_fixed", "n_stable")} for r in rows_b]).to_csv(
            TABLES_DIR / "dynamics_bifurcation.csv", index=False)
        run.record_output(TABLES_DIR / "dynamics_bifurcation.csv")
        save_figure(figures.bifurcation(rows_b, bp, len(model.state_names) - 1, model.state_names[-1]),
                    "dynamics_bifurcation", run)

        curves = {}
        for name in drugs:
            d, y = dose_curve(p0, name)
            curves[name] = (d, y, is_u_shaped(d, y, min_relative_rise=rise))
        protective = {n: bool(np.min(c[1][1:]) < c[1][0]) for n, c in curves.items()}
        run.info(f"at baseline parameters: protective at some dose = {protective}; "
                 f"both atomoxetine and prazosin protective = {protective['atomoxetine'] and protective['prazosin']}")
        save_figure(figures.dose_response(curves), "dynamics_dose_response", run)

        swept = {k: p.range for k, p in params.params.items() if p.range is not None}
        sob = sobol_indices(lambda q: outcome_for(params.values(q)), swept,
                            n_base=cfg.get("sensitivity.sobol_n_base"), seed=seed)
        pd.DataFrame({"parameter": sob["names"], "S1": sob["S1"], "S1_conf": sob["S1_conf"],
                      "ST": sob["ST"], "ST_conf": sob["ST_conf"]}).to_csv(TABLES_DIR / "dynamics_sobol.csv", index=False)
        run.record_output(TABLES_DIR / "dynamics_sobol.csv")
        save_figure(figures.sobol_bars(sob), "dynamics_sobol", run)

        fractions = {}
        for name in ("atomoxetine", "prazosin"):
            fractions[f"U-shape: {name}"] = lhs_robustness(lambda q, n=name: dose_curve(params.values(q), n), swept,
                                                           n_samples=cfg.get("sensitivity.lhs_n_samples"), seed=seed,
                                                           min_relative_rise=rise)
        both = fractions["U-shape: atomoxetine"]["samples"]
        k = 0
        for row in both:
            q = params.values(dict(zip(fractions["U-shape: atomoxetine"]["names"], row)))
            ys = [dose_curve(q, n)[1] for n in ("atomoxetine", "prazosin")]
            k += all(np.min(y[1:]) < y[0] for y in ys)
        n = len(both)
        fractions["both drugs protective"] = {"fraction": k / n, "ci95": wilson_interval(k, n), "n_samples": n}
        for name, f in fractions.items():
            run.info(f"robustness — {name}: {100 * f['fraction']:.1f}% of {f['n_samples']} samples "
                     f"(95% CI {100 * f['ci95'][0]:.1f}-{100 * f['ci95'][1]:.1f}%)")
        pd.DataFrame([{"quantity": n_, "fraction": f["fraction"], "ci95_low": f["ci95"][0], "ci95_high": f["ci95"][1],
                       "n_samples": f["n_samples"]} for n_, f in fractions.items()]).to_csv(
            TABLES_DIR / "dynamics_robustness.csv", index=False)
        run.record_output(TABLES_DIR / "dynamics_robustness.csv")
        save_figure(figures.robustness(fractions), "dynamics_robustness", run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
