"""Module 4 entry point: python -m src.dynamics.run_dynamics  (or `make dynamics`)

The registered analyses of prereg/dynamics.md. Refuses to start unless
Module 4's validation passed on the current code and prereg/dynamics.md is
frozen. Weak parameters are swept, never reported as point values.

1. Reference point (illustration only): trajectory from the healthy start,
   solver independence, fixed points (analytic and 4-D numerical), hA
   bifurcation diagram, dose-responses.
2. Latin hypercube (n = sensitivity.lhs_n_samples): every registered readout
   per sample; fractions with Wilson 95% CIs and the registered labels.
3. ODE re-check of the analytic readout on the first readout_check_n samples
   (stops the run if it fails: nothing is interpreted).
4. Sobol indices (Saltelli, N = sensitivity.sobol_n_base) for the registered outputs.
5. Sensitivity analyses SA1-SA5 (fractions and labels recomputed).
6. The verdict table, by the rules of prereg §2.
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
from src.dynamics.analysis import Space, dose_grid, evaluate
from src.dynamics.bifurcation import bistable_range, folds, scan
from src.dynamics.fixedpoints import find_fixed_points, grid_guesses
from src.dynamics.params import ParameterSet
from src.dynamics.sensitivity import label_fraction, lhs_sample, sobol_multi
from src.dynamics.simulate import simulate, solver_independence

CONFIG = CONFIG_DIR / "dynamics.yaml"
PARAMS = CONFIG_DIR / "dynamics_params.yaml"

# What each fraction would mean (prereg §2). "against" is the reading when the
# label is the one in the second slot.
READING = {
    "both_protective": ("fragile", "one parameter set making both atomoxetine and prazosin protective is not a "
                        "robust property of this model class"),
    "u_shape_atomoxetine": ("fragile", "a U-shaped atomoxetine dose-response is not a robust prediction"),
    "bistable": ("fragile", "a tipping point in LC regulation is not a robust property of this model class"),
    "harmful_atomoxetine": ("robust", "the model class predicts NET inhibition to be harmful, against the "
                            "pre-registered Module 1 direction"),
    "harmful_prazosin": ("robust", "the model class predicts alpha1 blockade to be harmful, against the "
                         "pre-registered Module 1 direction"),
}
SUPPORT = {"both_protective": "robust", "u_shape_atomoxetine": "robust", "bistable": "robust"}


def lhs_fractions(space: Space, cfg, n: int, seed: int, thresholds=None) -> tuple[dict, pd.DataFrame]:
    names, U = lhs_sample(space.sampling_ranges(), n_samples=n, seed=seed)
    rows = []
    for row in U:
        p = space.to_params(dict(zip(names, row)))
        r = evaluate(p, cfg, thresholds=thresholds)
        rows.append({**{k: p[k] for k in lc_model.PARAMETERS}, **r})
    df = pd.DataFrame(rows)
    lab = {k: label_fraction(int(df[k].sum()), len(df), fragile_upper=cfg.get("robustness.fragile_upper"),
                             robust_lower=cfg.get("robustness.robust_lower"), level=cfg.get("robustness.ci_level"))
           for k in cfg.get("robustness.fractions")}
    return lab, df


def main() -> int:
    require_validated("dynamics")
    prereg = require_preregistration("dynamics")
    cfg = load_config(CONFIG)
    params = ParameterSet.from_yaml(PARAMS)
    model = lc_model.build_model()
    seed = cfg.get("seed")
    t0, t1 = cfg.get("simulation.t_span")
    t_eval = np.linspace(t0, t1, cfg.get("simulation.n_out"))
    sim_kw = dict(method=cfg.get("simulation.method"), rtol=cfg.get("simulation.rtol"), atol=cfg.get("simulation.atol"))
    space = Space.from_params(params)
    doses = dose_grid(cfg.get("dose_response.dose_grid"))

    with start_run(module="dynamics", config_path=CONFIG, seed=seed,
                   packages=["numpy", "scipy", "SALib"], preregistration=prereg) as run:
        TABLES_DIR.mkdir(parents=True, exist_ok=True)

        def out_csv(df, stem):
            df.to_csv(TABLES_DIR / f"{stem}.csv", index=False)
            run.record_output(TABLES_DIR / f"{stem}.csv")

        run.info(f"parameters: {len(params.params)} ({len(params.weak())} weak, swept not reported)")
        out_csv(pd.DataFrame([{"parameter": k, "units": p.units, "source": p.source, "confidence": p.confidence,
                               "range_low": p.range[0], "range_high": p.range[1], "scale": p.scale}
                              for k, p in params.params.items()]), "dynamics_parameters")

        # 1. reference point (illustration; its values are not results)
        p0 = params.values()
        y0 = lc_model.initial_state(p0)
        traj = simulate(model, p0, y0, (t0, t1), t_eval, **sim_kw)
        rep = solver_independence(model, p0, y0, (t0, t1), t_eval, rtol=sim_kw["rtol"], atol=sim_kw["atol"])
        run.info(f"solver independence at the reference point: {rep}")
        if not (rep["converging"] and rep["cross_solver_tight"] <= cfg.get("simulation.max_solver_disagreement")):
            raise RuntimeError(f"results depend on the numerics: {rep}")
        save_figure(figures.trajectory(traj, model.state_names), "dynamics_reference_trajectory", run)
        ana = lc_model.fixed_points(p0)
        guesses = grid_guesses(cfg.get("fixed_points.bounds"), model.state_names, cfg.get("fixed_points.grid_per_dim"))
        num = find_fixed_points(model, p0, np.vstack([guesses] + [f["y"] for f in ana]))
        if len(num) != len(ana) or any(not np.allclose(n.y, a["y"], rtol=1e-6) or n.stable != a["stable"]
                                       for n, a in zip(sorted(num, key=lambda f: f.y[2]), ana)):
            raise RuntimeError("analytic and numerical fixed points disagree at the reference point")
        bp = cfg.get("bifurcation.parameter")
        lo, hi = cfg.get("bifurcation.range")
        rows_b = scan(model, p0, bp, np.linspace(lo, hi, cfg.get("bifurcation.n_values")), guesses)
        run.info(f"reference point, scan over {bp}: folds {folds(rows_b)}; grid bistable range "
                 f"{bistable_range(rows_b)}; exact window {lc_model.bistable_window(p0, lo, hi)}")
        save_figure(figures.bifurcation(rows_b, bp, 2, "R"), "dynamics_reference_bifurcation", run)
        curves = {}
        for name in ("atomoxetine", "prazosin"):
            y = lc_model.dose_response(p0, name, doses)
            curves[name] = (doses, y / y[0], False)
        save_figure(figures.dose_response(curves), "dynamics_reference_dose_response", run)

        # 2. Latin hypercube: the registered fractions
        n = cfg.get("sensitivity.lhs_n_samples")
        lab, df = lhs_fractions(space, cfg, n, seed)
        out_csv(df, "dynamics_lhs_samples")
        for k, f in lab.items():
            run.info(f"robustness - {k}: {f['k']}/{f['n']} = {100 * f['fraction']:.1f}% "
                     f"(Wilson 95% CI {100 * f['ci95'][0]:.1f}-{100 * f['ci95'][1]:.1f}%): {f['label']}")
        save_figure(figures.robustness({k: lab[k] for k in ("both_protective", "u_shape_atomoxetine", "bistable")}),
                    "dynamics_robustness", run)

        # 3. ODE re-check of the analytic readout
        n_chk = cfg.get("simulation.readout_check_n")
        worst = 0.0
        for _, r in df.head(n_chk).iterrows():
            p = {k: float(r[k]) for k in lc_model.PARAMETERS}
            fp = lc_model.healthy_steady_state(p)
            yfp = np.array([fp[s] for s in lc_model.STATES])
            slow = np.min(np.abs(np.linalg.eigvals(model.jacobian(yfp, p)).real))
            tr = simulate(model, p, y0, (0.0, 50.0 / slow), np.array([0.0, 50.0 / slow]), **sim_kw)
            worst = max(worst, abs(tr.y[3, -1] - fp["T"]) / fp["T"])
        run.info(f"ODE re-check of the readout on {n_chk} samples: max relative difference {worst:.2e}")
        if worst > cfg.get("simulation.readout_check_max_rel_diff"):
            raise RuntimeError(f"analytic readout disagrees with ODE integration ({worst:.2e}); nothing is interpreted")

        # 4. Sobol
        outs = cfg.get("sensitivity.sobol_outputs")
        sob = sobol_multi(lambda q: {k: float(v) for k, v in evaluate(space.to_params(q), cfg).items() if k in outs},
                          space.sampling_ranges(), outs, n_base=cfg.get("sensitivity.sobol_n_base"), seed=seed)
        rows = []
        for o, res in sob.items():
            if "constant" in res:
                run.info(f"Sobol {o}: constant ({res['constant']}) over the sample; no indices")
                continue
            rows += [{"output": o, "parameter": nm, "S1": a, "S1_conf": b, "ST": c, "ST_conf": d}
                     for nm, a, b, c, d in zip(res["names"], res["S1"], res["S1_conf"], res["ST"], res["ST_conf"])]
            save_figure(figures.sobol_bars(res), f"dynamics_sobol_{o}", run)
        out_csv(pd.DataFrame(rows), "dynamics_sobol")

        # 5. sensitivity analyses
        sa_rows = []
        sa = cfg.get("sensitivity_analyses")
        variant_labels = {}
        for name, spec in sa.items():
            if name.startswith("SA5"):
                for prot, harm in spec:
                    l5, _ = lhs_fractions(space, cfg, n, seed, thresholds=(prot, harm))
                    sa_rows += [{"analysis": f"{name} ({prot}, {harm})", "quantity": k, **_flat(f)} for k, f in l5.items()]
                continue
            lv, _ = lhs_fractions(space.variant(spec), cfg, n, seed)
            variant_labels[name] = lv
            sa_rows += [{"analysis": name, "quantity": k, **_flat(f)} for k, f in lv.items()]
        out_csv(pd.DataFrame([{"analysis": "primary", "quantity": k, **_flat(f)} for k, f in lab.items()] + sa_rows),
                "dynamics_robustness")

        # 6. verdicts (prereg §2)
        vrows = []
        for k, (against_label, meaning) in READING.items():
            f = lab[k]
            same = all(variant_labels[v][k]["label"] == f["label"] for v in variant_labels)
            verdict = ("AGAINST: " + meaning if f["label"] == against_label
                       else "SUPPORTS" if SUPPORT.get(k) == f["label"] else "no verdict (depends on parameters)")
            vrows.append({"quantity": k, "fraction": f["fraction"], "ci95_low": f["ci95"][0], "ci95_high": f["ci95"][1],
                          "label": f["label"], "verdict": verdict, "label_same_in_SA1_to_SA4": same})
            run.info(f"VERDICT {k}: {verdict} (label {f['label']}; same in SA1-SA4: {same})")
        out_csv(pd.DataFrame(vrows), "dynamics_verdicts")
    return 0


def _flat(f: dict) -> dict:
    return {"k": f["k"], "n": f["n"], "fraction": f["fraction"], "ci95_low": f["ci95"][0], "ci95_high": f["ci95"][1],
            "label": f["label"]}


if __name__ == "__main__":
    sys.exit(main())
