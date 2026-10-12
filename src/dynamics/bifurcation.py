"""One-parameter bifurcation scan (natural-parameter continuation).

For each value of the swept parameter, fixed points are found from a grid of
guesses plus the previous value's fixed points; the number of stable fixed
points is recorded. Bistability = >= 2 stable fixed points. Fold (saddle-node)
locations are reported as the midpoints between adjacent parameter values
where the number of fixed points changes — so their precision is the scan
resolution, which is reported with them.
"""

from __future__ import annotations

import numpy as np

from src.dynamics.fixedpoints import find_fixed_points
from src.dynamics.model import ODEModel


def scan(model: ODEModel, p: dict, param: str, values: np.ndarray, guesses: np.ndarray) -> list[dict]:
    rows = []
    prev: list[np.ndarray] = []
    for v in np.asarray(values, float):
        q = dict(p, **{param: float(v)})
        g = np.vstack([guesses] + [f[None, :] for f in prev]) if prev else guesses
        fps = find_fixed_points(model, q, g)
        prev = [f.y for f in fps]
        rows.append({"value": float(v), "n_fixed": len(fps), "n_stable": sum(f.stable for f in fps),
                     "fixed_points": fps})
    return rows


def folds(rows: list[dict]) -> list[dict]:
    out = []
    for a, b in zip(rows, rows[1:]):
        if a["n_fixed"] != b["n_fixed"]:
            out.append({"location": 0.5 * (a["value"] + b["value"]), "resolution": b["value"] - a["value"],
                        "n_fixed_before": a["n_fixed"], "n_fixed_after": b["n_fixed"]})
    return out


def bistable_range(rows: list[dict]) -> tuple[float, float] | None:
    vals = [r["value"] for r in rows if r["n_stable"] >= 2]
    return (min(vals), max(vals)) if vals else None
