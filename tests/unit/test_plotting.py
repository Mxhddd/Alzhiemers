import matplotlib.pyplot as plt
import pytest

from src.common.config import load_config
from src.common.paths import CONFIG_DIR
from src.common.plotting import PALETTE, apply_theme, save_figure, series_colors
from src.common.runlog import start_run


def test_palette_is_fixed_and_unique():
    assert len(PALETTE) == 8 and len(set(PALETTE)) == 8


def test_no_cycling_past_eight():
    assert series_colors(3) == list(PALETTE[:3])
    with pytest.raises(ValueError):
        series_colors(9)


def test_save_figure_is_byte_reproducible_and_logged(tmp_path):
    cfg = CONFIG_DIR / "provenance_check.yaml"
    apply_theme()
    digests = []
    for i in range(2):
        with start_run(module="unit", config_path=cfg, seed=None, packages=[],
                       logs_dir=tmp_path / "logs") as run:
            fig, ax = plt.subplots()
            ax.plot([0, 1], [0, 1])
            paths = save_figure(fig, "t", run, figures_dir=tmp_path / f"fig{i}")
            digests.append({p.suffix: run.outputs[str(p.resolve())] for p in paths})
    assert digests[0] == digests[1]
    assert load_config(cfg).get("seed")  # config used above is intact
