"""Shared pieces of the Module 4 analyses (registered run and pilot).

* ``Space``: the swept parameter box on each parameter's sampling scale
  (log10 for log-scale parameters), with the pre-registered sensitivity
  variants (central half; fixed values; changed range or scale).
* ``evaluate``: every pre-registered readout for one parameter set.
* ``dose_grid``: the registered dose grid.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from src.dynamics import lc_model
from src.dynamics.params import ParameterSet


@dataclass
class Space:
    base: dict[str, float]                      # every parameter's reference value
    ranges: dict[str, tuple[float, float]]      # natural-scale ranges of the swept parameters
    scales: dict[str, str]
    fixed: dict[str, float] = field(default_factory=dict)

    @classmethod
    def from_params(cls, params: ParameterSet) -> "Space":
        swept = [k for k, p in params.params.items() if p.range is not None]
        return cls(base=params.values(), ranges=params.ranges(swept), scales=params.scales(swept))

    def variant(self, spec: dict) -> "Space":
        """A pre-registered sensitivity variant: {'shrink_to': f} keeps the
        central fraction f of every range on its sampling scale; {'fix': {...}}
        holds parameters fixed; {'range': {...}, 'scale': {...}} replaces them."""
        ranges, scales, fixed = dict(self.ranges), dict(self.scales), dict(self.fixed)
        unknown = set(spec) - {"shrink_to", "fix", "range", "scale"}
        if unknown:
            raise ValueError(f"unknown sensitivity-variant keys {unknown}")
        if "shrink_to" in spec:
            f = float(spec["shrink_to"])
            for k, (lo, hi) in self.sampling_ranges().items():
                mid, half = 0.5 * (lo + hi), 0.5 * f * (hi - lo)
                a, b = mid - half, mid + half
                ranges[k] = (10 ** a, 10 ** b) if scales[k] == "log" else (a, b)
        for k, v in spec.get("fix", {}).items():
            ranges.pop(k)
            scales.pop(k)
            fixed[k] = float(v)
        for k, r in spec.get("range", {}).items():
            ranges[k] = (float(r[0]), float(r[1]))
        for k, s in spec.get("scale", {}).items():
            scales[k] = s
        for k, (lo, hi) in ranges.items():
            if scales[k] == "log" and lo <= 0:
                raise ValueError(f"{k}: log scale needs a positive range")
        return Space(base=self.base, ranges=ranges, scales=scales, fixed=fixed)

    def sampling_ranges(self) -> dict[str, tuple[float, float]]:
        return {k: ((np.log10(lo), np.log10(hi)) if self.scales[k] == "log" else (lo, hi))
                for k, (lo, hi) in self.ranges.items()}

    def to_params(self, sample: dict[str, float]) -> dict[str, float]:
        """A point of the sampling box -> a full parameter dict."""
        p = dict(self.base, **self.fixed)
        for k, x in sample.items():
            p[k] = float(10 ** x) if self.scales[k] == "log" else float(x)
        return p


def dose_grid(spec: dict) -> np.ndarray:
    d = np.logspace(np.log10(spec["min"]), np.log10(spec["max"]), spec["n_log"])
    return np.concatenate([[0.0], d]) if spec["include_zero"] else d


def evaluate(p: dict, cfg, *, with_dose_curves: bool = True, thresholds: tuple[float, float] | None = None) -> dict:
    """All pre-registered readouts (prereg §4) for one parameter set."""
    prot, harm = thresholds if thresholds else (cfg.get("dose_response.protective_ratio"),
                                                cfg.get("dose_response.harmful_ratio"))
    r = lc_model.readouts(p, reference_dose=cfg.get("dose_response.reference_dose"),
                          protective_ratio=prot, harmful_ratio=harm,
                          dose_grid=dose_grid(cfg.get("dose_response.dose_grid")) if with_dose_curves else None,
                          min_relative_rise=cfg.get("dose_response.min_relative_rise"),
                          bifurcation_range=tuple(cfg.get("bifurcation.range")))
    r["log_T0"] = float(np.log(r["T0"]))
    r["log_ratio_atomoxetine"] = float(np.log(r["ratio_atomoxetine"]))
    r["log_ratio_prazosin"] = float(np.log(r["ratio_prazosin"]))
    return r
