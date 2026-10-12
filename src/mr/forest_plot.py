"""Module 1 forest plot: primary MR estimate per target, with the
pre-registered prediction marked, so a reader sees at a glance which held.

Reads results/tables/mr_summary.csv (written by src/mr/run_mr.R). Excluded
targets stay in the figure as labelled gaps — they are results too.

Usage: python -m src.mr.forest_plot --summary results/tables/mr_summary.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.common.paths import CONFIG_DIR, TABLES_DIR
from src.common.plotting import DIVERGING, INK, apply_theme, save_figure
from src.common.runlog import start_run

PREDICTED_SIGN = {"protective": 1, "harmful": -1, "exploratory": 0}


def forest_rows(summary: pd.DataFrame) -> pd.DataFrame:
    """One row per gene x exposure source: the primary main-analysis estimate,
    or the exclusion reason. Fails if a gene/source has neither."""
    need = {"gene", "exposure_source", "predicted_lower_activity", "status"}
    missing = need - set(summary.columns)
    if missing:
        raise ValueError(f"mr_summary.csv lacks columns {sorted(missing)}")
    rows = []
    for (gene, src), g in summary.groupby(["gene", "exposure_source"], sort=False):
        if "is_primary" in g and "analysis" in g:
            prim = g[(g["is_primary"] == True) & (g["analysis"] == "main")]  # noqa: E712
        else:
            prim = g.iloc[0:0]
        if len(prim) == 1:
            r = prim.iloc[0]
            rows.append({"label": f"{gene} [{src}]", "pred": r["predicted_lower_activity"],
                         "b": r["b"], "lo": r["ci_lower"], "hi": r["ci_upper"],
                         "note": f"{r['method']}, {int(r['nsnp'])} SNP · {r['concordance']}"})
        elif len(prim) == 0 and g["status"].astype(str).str.startswith("excluded").all():
            rows.append({"label": f"{gene} [{src}]", "pred": g["predicted_lower_activity"].iloc[0],
                         "b": np.nan, "lo": np.nan, "hi": np.nan, "note": g["status"].iloc[0]})
        else:
            raise ValueError(f"{gene} [{src}]: expected one primary estimate or an exclusion, got {len(prim)}")
    return pd.DataFrame(rows)


def draw(rows: pd.DataFrame):
    apply_theme()
    n = len(rows)
    fig, ax = plt.subplots(figsize=(7.5, 0.42 * n + 1.2))
    finite = rows[["lo", "hi"]].to_numpy(dtype=float)
    span = np.nanmax(np.abs(finite)) if np.isfinite(finite).any() else 1.0
    lim = 1.15 * span
    for i, r in enumerate(rows.itertuples()):
        y = n - 1 - i
        sign = PREDICTED_SIGN[r.pred]
        if sign:  # shade the half-axis the prediction points to
            x0, x1 = (0, lim) if sign > 0 else (-lim, 0)
            ax.axhspan(y - 0.4, y + 0.4, xmin=(x0 + lim) / (2 * lim), xmax=(x1 + lim) / (2 * lim),
                       color=DIVERGING["high"] if sign > 0 else DIVERGING["low"], alpha=0.12, lw=0)
        if np.isfinite(r.b):
            ax.plot([r.lo, r.hi], [y, y], color=INK["primary"], lw=1.5, solid_capstyle="round")
            ax.plot(r.b, y, "o", color=INK["primary"], ms=5)
        else:
            ax.text(0, y, "no estimate", ha="center", va="center", fontsize=7, color=INK["muted"])
        ax.text(lim * 1.03, y, r.note, va="center", fontsize=7, color=INK["secondary"])
    ax.axvline(0, color=INK["axis"], lw=0.8)
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(rows["label"])[::-1])
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-0.7, n - 0.3)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("MR estimate: log-OR of Alzheimer's per SD higher expression")
    ax.set_title("Drug-target MR — shaded half-axis = pre-registered predicted direction\n"
                 "(red: higher expression harmful, i.e. lower activity protective; blue: the reverse; "
                 "none: exploratory)", fontsize=8, loc="left")
    return fig


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--summary", type=Path, default=TABLES_DIR / "mr_summary.csv")
    args = p.parse_args(argv)
    if not args.summary.is_file():
        raise FileNotFoundError(f"{args.summary} not found: run `make mr` first")
    with start_run(module="mr_forest", config_path=CONFIG_DIR / "mr_targets.yaml", seed=None,
                   packages=["pandas", "matplotlib"]) as run:
        rows = forest_rows(pd.read_csv(args.summary))
        run.info(f"{len(rows)} target x source rows; {rows['b'].isna().sum()} without an estimate")
        save_figure(draw(rows), "mr_forest", run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
