"""Strict loaders for the researcher's hand-transcribed files.

Nothing is imputed. In particular, a stage that is absent from an age band is
NOT taken to be zero: every (age band, source) group must list every stage
explicitly, with n = 0 where the paper reports none. Otherwise the loader
stops and says which group is incomplete.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

BRAAK_COLUMNS = ["age_lower", "age_upper", "stage", "n", "source"]
BAND_COLUMNS = ["age_lower", "age_upper", "source"]


@dataclass(frozen=True)
class StageCounts:
    """One observation group: counts per stage in one age band from one source."""

    age_lower: float
    age_upper: float
    source: str
    counts: np.ndarray  # length K, in stage_order

    @property
    def midpoint(self) -> float:
        return 0.5 * (self.age_lower + self.age_upper)

    @property
    def n(self) -> int:
        return int(self.counts.sum())


def _read_exact(path: Path, columns: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"stage": str, "source": str}, keep_default_na=False, na_values=[""])
    if list(df.columns) != columns:
        raise ValueError(f"{path}: columns must be exactly {columns}, found {list(df.columns)}")
    if df.isna().any().any():
        rows = df.index[df.isna().any(axis=1)].tolist()
        raise ValueError(f"{path}: empty cells in rows {[r + 2 for r in rows]} (line numbers)")
    bad = df["age_upper"] <= df["age_lower"]
    if bad.any():
        raise ValueError(f"{path}: age_upper <= age_lower in lines {[i + 2 for i in df.index[bad]]}")
    return df


def load_braak_counts(path: Path, stage_order: list[str]) -> list[StageCounts]:
    df = _read_exact(Path(path), BRAAK_COLUMNS)
    stage_order = [str(s) for s in stage_order]
    unknown = sorted(set(df["stage"]) - set(stage_order))
    if unknown:
        raise ValueError(f"{path}: stage labels {unknown} are not in stage_order {stage_order}")
    n = df["n"]
    if not (np.issubdtype(n.dtype, np.integer) and (n >= 0).all()):
        raise ValueError(f"{path}: n must be non-negative integers")
    dup = df.duplicated(["age_lower", "age_upper", "source", "stage"])
    if dup.any():
        raise ValueError(f"{path}: duplicated (band, source, stage) rows at lines {[i + 2 for i in df.index[dup]]}")
    groups = []
    for (lo, hi, src), g in df.groupby(["age_lower", "age_upper", "source"], sort=True):
        missing = [s for s in stage_order if s not in set(g["stage"])]
        if missing:
            raise ValueError(f"{path}: band {lo}-{hi} ({src}) does not list stages {missing}; "
                             "enter n = 0 explicitly if the source reports none")
        counts = np.array([int(g.loc[g["stage"] == s, "n"].iloc[0]) for s in stage_order])
        if counts.sum() == 0:
            raise ValueError(f"{path}: band {lo}-{hi} ({src}) has zero total count")
        groups.append(StageCounts(float(lo), float(hi), str(src), counts))
    return groups


def load_age_bands(path: Path, value_column: str) -> pd.DataFrame:
    """Amyloid prevalence or dementia incidence by age band."""
    cols = ["age_lower", "age_upper", value_column, "source"]
    df = _read_exact(Path(path), cols)
    if (df[value_column] < 0).any():
        raise ValueError(f"{path}: negative {value_column}")
    return df
