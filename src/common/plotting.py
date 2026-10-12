"""Shared plotting theme for every figure in the project.

Choices (each explained so they can be challenged):

* **Font: DejaVu Sans.** It ships inside matplotlib, so a figure renders
  identically on any machine; a system font would not.
* **Categorical colours: a fixed 8-slot order** validated for colour-vision
  deficiency (worst adjacent-pair CVD ΔE 9.1, normal-vision ΔE 19.6 on a
  near-white surface). Slots are assigned in order and never cycled; a 9th
  series must be folded into "Other" or split into panels. In charts where
  every pair can sit next to every other (scatter plots, small multiples),
  only the first three slots are safe. Slots 3, 4, 5 are below 3:1 contrast
  on white, so series in those colours need direct labels.
* **Ink:** text is near-black/grey, never the series colour.
* **Byte-reproducible output:** PDF creation dates and SVG hash salts are
  fixed, so re-running ``make all`` on the same commit yields identical files
  and their hashes in the run log can be compared across runs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless, deterministic backend
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from src.common.paths import FIGURES_DIR  # noqa: E402

# Categorical slots, in the order they must be used.
PALETTE: tuple[str, ...] = (
    "#2a78d6",  # 1 blue
    "#eb6834",  # 2 orange
    "#1baf7a",  # 3 aqua      (<3:1 on white: label directly)
    "#eda100",  # 4 yellow    (<3:1 on white: label directly)
    "#e87ba4",  # 5 magenta   (<3:1 on white: label directly)
    "#008300",  # 6 green
    "#6250d6",  # 7 violet
    "#e34948",  # 8 red
)
# Diverging pair (e.g. protective vs harmful) with a neutral grey midpoint.
DIVERGING = {"low": "#2a78d6", "mid": "#f0efec", "high": "#e34948"}

INK = {
    "primary": "#0b0b0b",
    "secondary": "#52514e",
    "muted": "#898781",
    "grid": "#e1e0d9",
    "axis": "#c3c2b7",
}

THEME_RC = {
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "text.color": INK["primary"],
    "axes.labelcolor": INK["secondary"],
    "axes.edgecolor": INK["axis"],
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": INK["grid"],
    "grid.linewidth": 0.5,
    "xtick.color": INK["muted"],
    "ytick.color": INK["muted"],
    "xtick.labelcolor": INK["secondary"],
    "ytick.labelcolor": INK["secondary"],
    "axes.prop_cycle": matplotlib.cycler(color=list(PALETTE)),
    "lines.linewidth": 1.5,
    "lines.markersize": 5,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "figure.dpi": 100,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,  # embed TrueType so journals can edit text
    "svg.hashsalt": "lc-first",
}


def apply_theme() -> None:
    matplotlib.rcParams.update(THEME_RC)


def series_colors(n: int) -> list[str]:
    """First ``n`` categorical slots. Raises rather than cycling past 8."""
    if n > len(PALETTE):
        raise ValueError(
            f"{n} series requested but only {len(PALETTE)} distinguishable colours "
            "exist; fold the tail into 'Other' or use small multiples"
        )
    return list(PALETTE[:n])


def save_figure(fig: Figure, stem: str, run, *, figures_dir: Path = FIGURES_DIR,
                formats: tuple[str, ...] = ("pdf", "png")) -> list[Path]:
    """Save ``fig`` as results/figures/<stem>.<fmt>, stamping the commit hash
    into the file metadata and recording each file's hash in the run log."""
    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)
    commit = run.header["git"]["commit"]
    paths = []
    for fmt in formats:
        out = figures_dir / f"{stem}.{fmt}"
        if fmt == "pdf":
            meta = {"Creator": f"lc-first {commit}", "Producer": None, "CreationDate": None}
        elif fmt == "png":
            meta = {"Software": f"lc-first {commit}"}
        else:
            meta = {}
        fig.savefig(out, format=fmt, metadata=meta)
        run.record_output(out)
        paths.append(out)
    plt.close(fig)
    return paths
