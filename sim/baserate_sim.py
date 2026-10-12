"""SIMULATED autopsy stage counts for Module 2 validation (simulation recovery).

Ground rule 2: simulated data live only under sim/ or tests/, are generated
by a visible script (this one), and are never read by an analysis pipeline —
src/common/io.py refuses any path under sim/. Nothing here is written to
disk by default; tests call ``simulate_counts`` directly.
"""

from __future__ import annotations

import numpy as np

from src.baserate.data import StageCounts
from src.baserate.model import StageModel
from src.baserate.solve import stage_probabilities


def simulate_counts(model: StageModel, theta: np.ndarray, bands: list[tuple[float, float]],
                    n_per_band: int, *, step: float, start_age: float, seed: int) -> list[StageCounts]:
    rng = np.random.default_rng(seed)
    P = stage_probabilities(model, theta, [(lo + hi) / 2 for lo, hi in bands], step=step, start_age=start_age)
    return [StageCounts(float(lo), float(hi), "SIMULATED", rng.multinomial(n_per_band, p / p.sum()))
            for (lo, hi), p in zip(bands, P)]
