"""Module 2 entry point: fit variants A, B, C to the autopsy stage counts.

    python -m src.baserate.run_baserate        (or `make baserate`)

Refuses to start unless the Module 2 validation passed on the current code
and prereg/baserate.md is frozen. Then:
  1. loads the researcher's files strictly (nothing imputed);
  2. per variant: multi-start MLE, step-error audit, NUTS posterior with
     R-hat / ESS / divergences;
  3. k-fold cross-validation by age band;
  4. predicted vs observed stage prevalence;
  5. falsification test: predicted dementia incidence vs published incidence
     the model never saw, judged by the researcher's pre-set tolerance;
  6. gamma vs the researcher's hand-derived value — flagged loudly if outside
     the credible interval.
"""

from __future__ import annotations

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.baserate.data import load_age_bands, load_braak_counts
from src.baserate.fit import fit_bayesian, fit_mle
from src.baserate.model import AGated, BIndependent, ChainSpec, CMultiHit, piecewise_constant
from src.baserate.predict import compare_models, predict_incidence, stage_prevalence_table
from src.baserate.solve import max_step_error
from src.common.config import load_config
from src.common.io import analysis_input
from src.common.paths import CONFIG_DIR, REPO_ROOT, TABLES_DIR
from src.common.plotting import INK, apply_theme, save_figure, series_colors
from src.common.prereg import require_preregistration
from src.common.runlog import start_run
from src.common.validation import require_validated

CONFIG = CONFIG_DIR / "baserate.yaml"
RESEARCHER_FIELDS = ["data.stage_order", "data.first_cortical_stage", "data.dementia_threshold_stage",
                     "solver.step_years", "falsification.max_abs_log_ratio"]


def main() -> int:
    require_validated("baserate")
    prereg = require_preregistration("baserate")
    cfg = load_config(CONFIG)
    missing = []
    for k in RESEARCHER_FIELDS:
        try:
            cfg.get(k)
        except KeyError:
            missing.append(k)
    if missing:
        raise KeyError(f"{CONFIG}: the researcher must supply {missing}")

    with start_run(module="baserate", config_path=CONFIG, seed=cfg.get("seed"),
                   packages=["numpy", "scipy", "pandas", "pymc", "arviz"], preregistration=prereg) as run:
        TABLES_DIR.mkdir(parents=True, exist_ok=True)
        stages = [str(s) for s in cfg.get("data.stage_order")]
        spec = ChainSpec.from_labels(stages, cfg.get("data.first_cortical_stage"))
        step, start = cfg.get("solver.step_years"), cfg.get("solver.start_age")
        centre, scale = cfg.get("solver.age_centre"), cfg.get("solver.age_scale")
        if cfg.get("solver.evaluate_at") != "midpoint":
            raise ValueError("only solver.evaluate_at: midpoint is implemented")
        data = load_braak_counts(analysis_input(REPO_ROOT / cfg.get("data.braak_counts"),
                                                provided_by="researcher, hand-transcribed Braak counts"), stages)
        run.info(f"{len(data)} age-band groups, {sum(g.n for g in data)} brains, {len(stages)} stages")
        amy = load_age_bands(analysis_input(REPO_ROOT / cfg.get("data.amyloid_prevalence"),
                                            provided_by="researcher, amyloid prevalence by age"), "prevalence")
        if (amy["prevalence"] > 1).any():
            raise ValueError("amyloid prevalence must be a probability")
        amyloid = piecewise_constant(list(amy[["age_lower", "age_upper", "prevalence"]].itertuples(index=False, name=None)))
        builders = {"A_gated": lambda: AGated(spec, amyloid, centre, scale),
                    "B_independent": lambda: BIndependent(spec, centre, scale),
                    "C_multihit": lambda: CMultiHit(spec, cfg.get("variants.c_hits"), centre, scale)}
        models = [builders[v]() for v in cfg.get("variants.fit")]
        seed = cfg.get("seed")

        param_rows, diag_rows = [], []
        fits = {}
        for m in models:
            fit = fit_mle(m, data, step=step, start_age=start, cfg_mle=cfg.get("mle"), seed=seed)
            fits[m.name] = fit
            for msg in fit.messages:
                if msg.startswith("WARNING") or msg.startswith("NOT"):
                    run.warning(f"{m.name}: {msg}")
            err = max_step_error(m, fit.theta, [g.midpoint for g in data], step=step, start_age=start)
            run.info(f"{m.name}: logLik {fit.loglik:.2f}, AIC {fit.aic:.2f}, {fit.n_at_best}/{fit.converged_starts} "
                     f"starts at best, Hessian {'ok' if fit.hessian_ok else 'NOT positive definite'}, "
                     f"max |p_step - p_ode| = {err:.2e} at step {step}")
            idata, summ, diag = fit_bayesian(m, data, step=step, start_age=start, cfg_bayes=cfg.get("bayes"), seed=seed)
            bad = diag["max_rhat"] > 1.01 or diag["divergences"] > 0
            (run.warning if bad else run.info)(f"{m.name} NUTS: max R-hat {diag['max_rhat']:.3f}, "
                                               f"min bulk ESS {diag['min_ess_bulk']:.0f}, divergences {diag['divergences']}")
            diag_rows.append({"model": m.name, "loglik": fit.loglik, "aic": fit.aic, "step_error": err,
                              "hessian_ok": fit.hessian_ok, "starts_at_best": fit.n_at_best, **diag})
            for i, name in enumerate(m.param_names):
                key = name.replace("[", "_").replace("]", "")
                r = summ.loc[key]
                lo_col, hi_col = [c for c in summ.columns if c.startswith("hdi_")]
                param_rows.append({"model": m.name, "parameter": name, "mle": fit.theta[i], "mle_se": fit.se[i],
                                   "posterior_mean": r["mean"], "ci_lower": r[lo_col], "ci_upper": r[hi_col],
                                   "r_hat": r["r_hat"], "ess_bulk": r["ess_bulk"]})
        params = pd.DataFrame(param_rows)
        for name, df in (("baserate_parameters.csv", params), ("baserate_diagnostics.csv", pd.DataFrame(diag_rows))):
            df.to_csv(TABLES_DIR / name, index=False)
            run.record_output(TABLES_DIR / name)

        # gamma vs researcher's hand derivation
        g = params[(params["model"] == "A_gated") & (params["parameter"] == "gamma")]
        expected = cfg.as_dict()["gamma_check"]["expected"]
        if len(g) and expected is None:
            run.warning("gamma_check.expected not supplied: the fitted gamma is NOT compared with the hand derivation")
        elif len(g):
            inside = g["ci_lower"].iloc[0] <= expected <= g["ci_upper"].iloc[0]
            msg = (f"gamma posterior {g['posterior_mean'].iloc[0]:.4g} [{g['ci_lower'].iloc[0]:.4g}, "
                   f"{g['ci_upper'].iloc[0]:.4g}] vs researcher's {expected}")
            if inside:
                run.info(msg + ": consistent")
            else:
                run.warning("!!! GAMMA DISAGREES WITH THE HAND DERIVATION !!! " + msg)

        cv = compare_models(models, data, k=cfg.get("cross_validation.k"), step=step, start_age=start,
                            cfg_mle=cfg.get("mle"), seed=seed)
        cv_df = pd.DataFrame([{k: v for k, v in r.items() if k != "per_fold"} for r in cv])
        cv_df.to_csv(TABLES_DIR / "baserate_cv.csv", index=False)
        run.record_output(TABLES_DIR / "baserate_cv.csv")

        prev = pd.DataFrame([dict(row, model=m.name) for m in models
                             for row in stage_prevalence_table(m, fits[m.name].theta, data, step=step, start_age=start)])
        prev.to_csv(TABLES_DIR / "baserate_prevalence_fit.csv", index=False)
        run.record_output(TABLES_DIR / "baserate_prevalence_fit.csv")

        # Falsification: incidence the model never saw
        inc = load_age_bands(analysis_input(REPO_ROOT / cfg.get("data.dementia_incidence"),
                                            provided_by="researcher, published dementia incidence"), "rate_per_1000_py")
        thr = stages.index(str(cfg.get("data.dementia_threshold_stage")))
        mids = ((inc["age_lower"] + inc["age_upper"]) / 2).to_numpy()
        rows = []
        tol = cfg.get("falsification.max_abs_log_ratio")
        for m in models:
            pred = 1000 * predict_incidence(m, fits[m.name].theta, mids, threshold=thr, step=step, start_age=start)
            lr = np.log(pred / inc["rate_per_1000_py"].to_numpy())
            verdict = "PASS" if np.all(np.abs(lr) <= tol) else "FAIL"
            (run.info if verdict == "PASS" else run.warning)(
                f"falsification {m.name}: {verdict} (max |log ratio| {np.abs(lr).max():.2f}, tolerance {tol})")
            for i in range(len(inc)):
                rows.append({"model": m.name, "age_lower": inc["age_lower"][i], "age_upper": inc["age_upper"][i],
                             "observed_per_1000py": inc["rate_per_1000_py"][i], "predicted_per_1000py": pred[i],
                             "log_ratio": lr[i], "model_verdict": verdict, "source": inc["source"][i]})
        inc_df = pd.DataFrame(rows)
        inc_df.to_csv(TABLES_DIR / "baserate_incidence_falsification.csv", index=False)
        run.record_output(TABLES_DIR / "baserate_incidence_falsification.csv")

        apply_theme()
        colors = series_colors(len(models))
        fig, ax = plt.subplots(figsize=(6, 4))
        obs = inc_df.drop_duplicates(["age_lower", "age_upper"])
        ax.plot((obs["age_lower"] + obs["age_upper"]) / 2, obs["observed_per_1000py"], "o", color=INK["primary"],
                label="published (never fitted)")
        for m, c in zip(models, colors):
            d = inc_df[inc_df["model"] == m.name]
            ax.plot((d["age_lower"] + d["age_upper"]) / 2, d["predicted_per_1000py"], "-", color=c,
                    label=f"{m.name} ({d['model_verdict'].iloc[0]})")
        ax.set_yscale("log")
        ax.set_xlabel("Age (years)")
        ax.set_ylabel("Dementia incidence per 1000 person-years")
        ax.legend()
        save_figure(fig, "baserate_incidence_falsification", run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
