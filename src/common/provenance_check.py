"""Session 0 gate: run end to end and record provenance. No analysis.

Writes a run log under results/logs/provenance_check/ and a palette swatch to
results/figures/theme_check.{pdf,png} so the theme can be inspected by eye.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt

from src.common.config import load_config
from src.common.plotting import PALETTE, apply_theme, save_figure
from src.common.runlog import read_header, start_run


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args(argv)
    cfg = load_config(args.config)

    with start_run(module="provenance_check", config_path=args.config,
                   seed=cfg.get("seed"), packages=cfg.get("packages")) as run:
        apply_theme()
        fig, ax = plt.subplots(figsize=(6, 1.6))
        for i, colour in enumerate(PALETTE):
            ax.bar(i, 1, color=colour, width=0.9)
            ax.text(i, -0.12, f"{i + 1}\n{colour}", ha="center", va="top", fontsize=7)
        ax.set_axis_off()
        ax.set_title("LC-First categorical palette (slot order is fixed)")
        save_figure(fig, "theme_check", run)

        header = read_header(run.log_path)
        run.info(f"header commit={header['git']['commit']} "
                 f"config_sha256={header['config']['sha256']} "
                 f"dirty={header['git']['dirty']}")
    print(run.log_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
