"""Strict loaders for Module 2's extracted inputs (data/extracted/).

Nothing is imputed. A stage that is absent from an age band is NOT taken to
be zero: every (age band, source) group must list every stage explicitly,
with n = 0 where the paper reports none; otherwise the loader stops and says
which group is incomplete. Every file must carry provenance and a
`verified` flag (src/common/extracted.py); unverified rows stamp the run.

Age bands are written as in the papers, usually with inclusive integer upper
bounds ("61-70"). With upper_inclusive=True a band [lo, hi] becomes the
interval [lo, hi + 1).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.extracted import load_extracted


@dataclass(frozen=True)
class StageCounts:
    """One observation group: counts per stage in one age band from one source."""

    age_lower: float
    age_upper: float   # exclusive
    source: str
    counts: np.ndarray  # length K, in stage_order

    @property
    def midpoint(self) -> float:
        return 0.5 * (self.age_lower + self.age_upper)

    @property
    def n(self) -> int:
        return int(self.counts.sum())

    def eval_points(self, mode: str, n_points: int) -> tuple[np.ndarray, np.ndarray]:
        """Ages and weights at which the band's stage probabilities are taken.
        midpoint: the band centre. band_average: n equally spaced ages (cell
        centres), equal weights — i.e. ages uniform within the band."""
        if mode == "midpoint":
            return np.array([self.midpoint]), np.array([1.0])
        if mode == "band_average":
            w = (self.age_upper - self.age_lower) / n_points
            return self.age_lower + w * (np.arange(n_points) + 0.5), np.full(n_points, 1.0 / n_points)
        raise ValueError(f"unknown evaluation mode {mode!r}")


def _numeric(df: pd.DataFrame, cols: list[str], path) -> pd.DataFrame:
    for c in cols:
        try:
            df[c] = pd.to_numeric(df[c])
        except ValueError as e:
            raise ValueError(f"{path}: column {c} must be numeric") from e
    return df


def load_braak_counts(path: Path, stage_order: list[str], *, upper_inclusive: bool, run=None) -> list[StageCounts]:
    need = ["age_lower", "age_upper", "stage", "n", "source_id", "locator", "verified"]
    df = load_extracted(Path(path), provided_by="extracted from Braak et al. 2011 (see data/extracted/SOURCES.md)",
                        run=run)
    if list(df.columns) != need:
        raise ValueError(f"{path}: columns must be exactly {need}, found {list(df.columns)}")
    df = _numeric(df, ["age_lower", "age_upper", "n"], path)
    if ((df["n"] < 0) | (df["n"] != np.floor(df["n"]))).any():
        raise ValueError(f"{path}: n must be non-negative integers")
    if upper_inclusive:
        df["age_upper"] = df["age_upper"] + 1
    if (df["age_upper"] <= df["age_lower"]).any():
        raise ValueError(f"{path}: age_upper <= age_lower")
    stage_order = [str(s) for s in stage_order]
    unknown = sorted(set(df["stage"]) - set(stage_order))
    if unknown:
        raise ValueError(f"{path}: stage labels {unknown} are not in stage_order {stage_order}")
    if df.duplicated(["age_lower", "age_upper", "source_id", "stage"]).any():
        raise ValueError(f"{path}: duplicated (band, source, stage) rows")
    groups = []
    for (lo, hi, src), g in df.groupby(["age_lower", "age_upper", "source_id"], sort=True):
        missing = [s for s in stage_order if s not in set(g["stage"])]
        if missing:
            raise ValueError(f"{path}: band {lo}-{hi} ({src}) does not list stages {missing}; "
                             "enter n = 0 explicitly if the source reports none")
        counts = np.array([int(g.loc[g["stage"] == s, "n"].iloc[0]) for s in stage_order])
        if counts.sum() == 0:
            raise ValueError(f"{path}: band {lo}-{hi} ({src}) has zero total count")
        groups.append(StageCounts(float(lo), float(hi), str(src), counts))
    return groups


def load_amyloid_bands(path: Path, *, upper_inclusive: bool, run=None) -> list[tuple[float, float, float]]:
    """(lower, upper_exclusive, prevalence in [0, 1]) bands from Braak-style files."""
    df = load_extracted(Path(path), provided_by="extracted amyloid prevalence (see data/extracted/SOURCES.md)", run=run)
    df = _numeric(df, ["age_lower", "age_upper", "prevalence"], path)
    if ((df["prevalence"] < 0) | (df["prevalence"] > 1)).any():
        raise ValueError(f"{path}: prevalence must be a proportion in [0, 1]")
    hi = df["age_upper"] + (1 if upper_inclusive else 0)
    return list(zip(df["age_lower"].astype(float), hi.astype(float), df["prevalence"].astype(float)))


def jansen_bands(path: Path, run=None) -> list[tuple[float, float, float]]:
    """Jansen 2015 point estimates at ages 50, 55, ..., 90 -> the 5-year bands
    centred on them (47.5-52.5, ...), which is the paper's own binning (Table 3)."""
    df = load_extracted(Path(path), provided_by="extracted from Jansen et al. 2015", run=run)
    df = _numeric(df, ["age", "prevalence_pct"], path)
    return [(a - 2.5, a + 2.5, p / 100.0) for a, p in zip(df["age"], df["prevalence_pct"])]


def load_incidence(path: Path, *, open_band_age: float, run=None) -> pd.DataFrame:
    """CFAS incidence rows -> per band: evaluation age and the [min, max]
    envelope of the published rates across sexes and study waves."""
    df = load_extracted(Path(path), provided_by="extracted from Matthews et al. 2016", run=run)
    df = _numeric(df, ["age_lower", "rate_per_1000_py"], path)
    rows = []
    for (lo, hi), g in df.groupby(["age_lower", "age_upper"], sort=True):
        if hi == "open":
            age = float(open_band_age)
            if age <= lo:
                raise ValueError("open-band evaluation age must exceed the band's lower bound")
        else:
            age = 0.5 * (float(lo) + float(hi) + 1)   # integer band [lo, hi] = [lo, hi + 1)
        rows.append({"band": f"{int(lo)}-{hi}", "eval_age": age, "n_estimates": len(g),
                     "observed_min": float(g["rate_per_1000_py"].min()),
                     "observed_max": float(g["rate_per_1000_py"].max())})
    return pd.DataFrame(rows).sort_values("eval_age").reset_index(drop=True)
