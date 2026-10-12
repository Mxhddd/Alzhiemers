"""Module 4 figures (shared theme). Captions state numbers honestly: the
robustness figure's title carries the fraction of parameter space showing
the U-shape, whatever it is."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from src.common.plotting import INK, apply_theme, series_colors


def trajectory(traj, state_names):
    apply_theme()
    cols = series_colors(len(state_names))
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for i, (s, c) in enumerate(zip(state_names, cols)):
        ax.plot(traj.t, traj.y[i], color=c, label=s)
    ax.set_xlabel("time (model units)")
    ax.legend()
    return fig


def bifurcation(rows, param, state_index, state_name):
    apply_theme()
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for r in rows:
        for fp in r["fixed_points"]:
            ax.plot(r["value"], fp.y[state_index], "o", ms=2.2,
                    color=INK["primary"] if fp.stable else INK["muted"], mfc=None if fp.stable else "none")
    ax.set_xlabel(param)
    ax.set_ylabel(f"{state_name} at fixed point")
    ax.set_title("filled: stable · open: unstable", fontsize=8, loc="left")
    return fig


def dose_response(curves: dict[str, tuple[np.ndarray, np.ndarray, bool]]):
    apply_theme()
    cols = series_colors(len(curves))
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for (name, (d, y, u)), c in zip(curves.items(), cols):
        ax.plot(d, y, "-o", ms=3, color=c, label=f"{name} ({'U-shaped' if u else 'not U-shaped'})")
    ax.set_xlabel("dose (model units)")
    ax.set_ylabel("outcome (lower = better)")
    ax.legend()
    return fig


def sobol_bars(res):
    apply_theme()
    fig, ax = plt.subplots(figsize=(6, 0.35 * len(res["names"]) + 1.2))
    y = np.arange(len(res["names"]))
    ax.barh(y + 0.18, res["ST"], height=0.34, xerr=res["ST_conf"], color=series_colors(2)[0], label="total-order")
    ax.barh(y - 0.18, res["S1"], height=0.34, xerr=res["S1_conf"], color=series_colors(2)[1], label="first-order")
    ax.set_yticks(y)
    ax.set_yticklabels(res["names"])
    ax.set_xlabel("Sobol index")
    ax.legend()
    return fig


def robustness(fractions: dict[str, dict]):
    apply_theme()
    fig, ax = plt.subplots(figsize=(6, 0.5 * len(fractions) + 1.4))
    names = list(fractions)
    for i, n in enumerate(names):
        f = fractions[n]
        ax.plot([100 * f["ci95"][0], 100 * f["ci95"][1]], [i, i], color=INK["primary"], lw=1.5)
        ax.plot(100 * f["fraction"], i, "o", color=INK["primary"])
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names)
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of plausible parameter space (Latin hypercube; Wilson 95% CI)")
    ax.set_title("; ".join(f"{n}: {100 * fractions[n]['fraction']:.1f}%" for n in names), fontsize=8, loc="left")
    return fig
