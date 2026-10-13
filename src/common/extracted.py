"""Loading researcher-checkable extracted values (data/extracted/).

Every extracted table must say where each value came from and whether the
researcher has verified it: a ``source_id`` (or ``source``) column, a
``locator`` column (table/figure) unless the file uses a single ``source``
string, and a ``verified`` column holding true/false. ``load_extracted``
refuses files without them and, when a run is given, marks the run
unverified if any row is still false — which stamps every output.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.common.io import analysis_input


def load_extracted(path: Path, *, provided_by: str, run=None) -> pd.DataFrame:
    path = analysis_input(path, provided_by=provided_by)
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    if df.empty:
        raise ValueError(f"{path}: no rows (values not yet extracted — see data/extracted/SOURCES.md)")
    if "verified" not in df.columns:
        raise ValueError(f"{path}: needs a 'verified' column (true/false per row)")
    if not ({"source_id", "source"} & set(df.columns)):
        raise ValueError(f"{path}: needs a 'source_id' or 'source' column")
    if "source_id" in df.columns and "locator" not in df.columns:
        raise ValueError(f"{path}: rows citing a source_id need a 'locator' (table/figure) column")
    v = df["verified"].str.strip().str.lower()
    bad = ~v.isin(["true", "false"])
    if bad.any():
        raise ValueError(f"{path}: 'verified' must be true or false (lines {[i + 2 for i in df.index[bad]]})")
    blank = (df == "").any(axis=1)
    if blank.any():
        raise ValueError(f"{path}: empty cells on lines {[i + 2 for i in df.index[blank]]} — leave a value "
                         "out entirely rather than blank, or fill it")
    df["verified"] = v == "true"
    n_unv = int((~df["verified"]).sum())
    if run is not None and n_unv:
        run.mark_unverified(f"{path.name}: {n_unv} of {len(df)} values unverified")
    return df
