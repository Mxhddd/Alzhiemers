"""Module 5 figures: pooled forest plot and the triangulation figure.

The triangulation figure (brief §5, "the centrepiece"): per MR target, the
pooled observational estimate (drug use), the MR estimate and the
mechanistic prediction, all as the effect of LOWER target activity on log
risk, with a concordance flag. Missing estimates are drawn as gaps with a
label, never omitted. The two estimates are on different scales (drug use vs
per-SD expression), so only their signs and intervals are comparable — the
caption says so.

Usage: python -m src.metaanalysis.plots
"""

from __future__ import annotations

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.common.paths import CONFIG_DIR, TABLES_DIR
from src.common.plotting import DIVERGING, INK, apply_theme, save_figure, series_colors
from src.common.runlog import start_run

PRED_SIGN = {"protective": -1, "harmful": 1, "exploratory": 0}


def triangulation_figure(tri: pd.DataFrame):
    need = {"target", "prediction_lower_activity", "obs_est", "obs_ci_lower", "obs_ci_upper", "mr_est",
            "mr_ci_lower", "mr_ci_upper", "obs_vs_prediction", "mr_vs_prediction", "obs_vs_mr"}
    if need - set(tri.columns):
        raise ValueError(f"triangulation.csv lacks {sorted(need - set(tri.columns))}")
    apply_theme()
    obs_c, mr_c = series_colors(2)
    n = len(tri)
    vals = tri[["obs_ci_lower", "obs_ci_upper", "mr_ci_lower", "mr_ci_upper"]].to_numpy(float)
    lim = 1.15 * (np.nanmax(np.abs(vals)) if np.isfinite(vals).any() else 1.0)
    fig, ax = plt.subplots(figsize=(8, 0.55 * n + 1.6))
    for i, r in enumerate(tri.itertuples()):
        y = n - 1 - i
        s = PRED_SIGN[r.prediction_lower_activity]
        if s:
            lo, hi = (0, lim) if s > 0 else (-lim, 0)
            ax.axhspan(y - 0.42, y + 0.42, xmin=(lo + lim) / (2 * lim), xmax=(hi + lim) / (2 * lim),
                       color=DIVERGING["high"] if s > 0 else DIVERGING["low"], alpha=0.12, lw=0)
        for off, est, lo_, hi_, col, lab in ((0.15, r.obs_est, r.obs_ci_lower, r.obs_ci_upper, obs_c, "observational"),
                                             (-0.15, r.mr_est, r.mr_ci_lower, r.mr_ci_upper, mr_c, "MR")):
            if np.isfinite(est):
                ax.plot([lo_, hi_], [y + off] * 2, color=col, lw=1.5)
                ax.plot(est, y + off, "o", color=col, ms=5, label=lab if i == 0 else None)
            else:
                ax.text(-lim * 0.98, y + off, f"no {lab} estimate", fontsize=6.5, va="center", color=INK["muted"])
        flag = f"obs: {r.obs_vs_prediction} · MR: {r.mr_vs_prediction} · obs vs MR: {r.obs_vs_mr}"
        ax.text(lim * 1.03, y, flag, fontsize=6.5, va="center", color=INK["secondary"])
    ax.axvline(0, color=INK["axis"], lw=0.8)
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(tri["target"])[::-1])
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-0.7, n - 0.3)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Effect of LOWER target activity on log risk of dementia / AD\n"
                  "(observational: log ratio for drug use; MR: log-OR per SD lower expression — different scales; "
                  "compare signs, not sizes)", fontsize=7.5)
    ax.set_title("Triangulation — shaded half-axis = mechanistic prediction (red: lower activity harmful; "
                 "blue: protective; none: exploratory)", fontsize=8, loc="left")
    ax.legend(loc="lower left", fontsize=7)
    return fig


def forest_figure(pooled: pd.DataFrame):
    apply_theme()
    d = pooled.reset_index(drop=True)
    n = len(d)
    fig, ax = plt.subplots(figsize=(7, 0.4 * n + 1.2))
    for i, r in enumerate(d.itertuples()):
        y = n - 1 - i
        if np.isfinite(r.est):
            ax.plot(np.exp([r.ci_lower, r.ci_upper]), [y, y], color=INK["primary"], lw=1.5)
            ax.plot(np.exp(r.est), y, "s", color=INK["primary"], ms=5)
            if np.isfinite(r.pi_lower):
                ax.plot(np.exp([r.pi_lower, r.pi_upper]), [y, y], color=INK["muted"], lw=0.8, ls=":")
        ax.text(1.02, y, f"k={r.k} · {r.status}", transform=ax.get_yaxis_transform(), fontsize=6.5,
                va="center", color=INK["secondary"])
    ax.axvline(1, color=INK["axis"], lw=0.8)
    ax.set_xscale("log")
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(d["group"])[::-1])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Pooled ratio (random effects; dotted = 95% prediction interval)")
    return fig


def main() -> int:
    with start_run(module="metaanalysis_plots", config_path=CONFIG_DIR / "metaanalysis.yaml", seed=None,
                   packages=["pandas", "matplotlib"]) as run:
        for name in ("triangulation.csv", "meta_pooled_by_class.csv"):
            if not (TABLES_DIR / name).is_file():
                raise FileNotFoundError(f"{TABLES_DIR / name} missing: run src/metaanalysis/run_meta.R first")
        save_figure(triangulation_figure(pd.read_csv(TABLES_DIR / "triangulation.csv")), "triangulation", run)
        save_figure(forest_figure(pd.read_csv(TABLES_DIR / "meta_pooled_by_class.csv")), "meta_forest", run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
