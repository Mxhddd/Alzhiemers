"""Module 3 entry point (Killer 2): python -m src.mediation.run_mediation

Real-data path ONLY: reads config data.file from data/raw/ and fails loudly
while it is absent. Simulated cohorts are used only in tests/validation.
Refuses to start unless Module 3's validation passed on the current code and
the DAG pre-registration (prereg/mediation.md) is frozen.
"""

from __future__ import annotations

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.common.config import load_config
from src.common.io import analysis_input
from src.common.paths import CONFIG_DIR, REPO_ROOT, TABLES_DIR
from src.common.plotting import INK, apply_theme, save_figure, series_colors
from src.common.prereg import require_preregistration
from src.common.runlog import start_run
from src.common.validation import require_validated
from src.mediation.mediation import bootstrap, estimate, killer2_verdict, summarise

CONFIG = CONFIG_DIR / "mediation.yaml"


def load_cohort(cfg, run) -> pd.DataFrame:
    path = analysis_input(REPO_ROOT / cfg.get("data.file"),
                          provided_by="researcher: cohort with LC measures, risk factors and tau (no access yet)")
    df = pd.read_csv(path)
    cols = [*cfg.get("data.exposures"), cfg.get("data.mediator"), cfg.get("data.outcome"), *cfg.get("data.covariates")]
    absent = [c for c in cols if c not in df.columns]
    if absent:
        raise KeyError(f"{path}: columns {absent} not found")
    df = df[cols]
    bad = df.isna().any(axis=1)
    if bad.any():
        if cfg.get("data.missing") == "raise":
            raise ValueError(f"{int(bad.sum())} rows have missing values in {cols} (config data.missing = raise)")
        if cfg.get("data.missing") != "complete_case_logged":
            raise ValueError("data.missing must be 'raise' or 'complete_case_logged'")
        for i in df.index[bad]:
            run.drop(stage="load", item=f"row {i}", reason="missing: " + ",".join(df.columns[df.loc[i].isna()]))
        df = df[~bad]
    return df


def main() -> int:
    require_validated("mediation")
    prereg = require_preregistration("mediation")
    cfg = load_config(CONFIG)
    if cfg.get("model.type") != "linear_no_interaction":
        raise NotImplementedError("only model.type = linear_no_interaction is implemented")
    with start_run(module="mediation", config_path=CONFIG, seed=cfg.get("seed"),
                   packages=["numpy", "pandas"], preregistration=prereg) as run:
        TABLES_DIR.mkdir(parents=True, exist_ok=True)
        df = load_cohort(cfg, run)
        run.info(f"cohort: {len(df)} complete rows")
        covs = cfg.get("data.covariates")
        X = df[covs].to_numpy(float) if covs else None
        M = df[cfg.get("data.mediator")].to_numpy(float)
        Y = df[cfg.get("data.outcome")].to_numpy(float)
        grid = np.round(np.arange(cfg.get("sensitivity.rho_min"), cfg.get("sensitivity.rho_max") + 1e-9,
                                  cfg.get("sensitivity.rho_step")), 10)
        level = cfg.get("bootstrap.ci_level")
        rows, curves = [], {}
        for k, expo in enumerate(cfg.get("data.exposures")):
            T = df[expo].to_numpy(float)
            pe = estimate(T, M, Y, X)
            s = summarise(pe, bootstrap(T, M, Y, X, n_boot=cfg.get("bootstrap.n_boot"), seed=cfg.get("seed") + k),
                          level=level, rho_grid=grid)
            verdict = killer2_verdict(s, cfg.get("killer2.predicted_acme_sign"))
            run.info(f"{expo}: {verdict}")
            sens = s["sensitivity"]
            rows.append({"exposure": expo, **{f"{q}_{f}": v for q, d in s["estimates"].items() for f, v in d.items()},
                         "rho_at_zero_acme": sens["rho_at_zero_acme"],
                         "smallest_abs_rho_with_ci_including_zero": sens["rho_ci_includes_zero"],
                         "killer2": verdict, "n": len(T)})
            curves[expo] = sens
        out = pd.DataFrame(rows)
        out.to_csv(TABLES_DIR / "mediation_results.csv", index=False)
        run.record_output(TABLES_DIR / "mediation_results.csv")

        apply_theme()
        fig, ax = plt.subplots(figsize=(6, 4))
        for (expo, sens), c in zip(curves.items(), series_colors(len(curves))):
            ax.fill_between(sens["rho"], sens["ci_lower"], sens["ci_upper"], color=c, alpha=0.15, lw=0)
            ax.plot(sens["rho"], sens["acme"], color=c, label=expo)
        ax.axhline(0, color=INK["axis"], lw=0.8)
        ax.axvline(0, color=INK["axis"], lw=0.8, ls=":")
        ax.set_xlabel("rho: correlation of mediator- and outcome-model errors (unmeasured confounding)")
        ax.set_ylabel("ACME (indirect effect through LC integrity)")
        ax.set_title("Sensitivity of the mediated effect to unmeasured confounding (bands: 95% bootstrap)",
                     fontsize=8, loc="left")
        ax.legend()
        save_figure(fig, "mediation_sensitivity", run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
