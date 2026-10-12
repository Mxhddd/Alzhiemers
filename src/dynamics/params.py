"""Parameters with provenance (brief §5, Module 4: "every parameter carries
source and confidence").

* ``confidence`` is one of strong / moderate / weak.
* A **weak** parameter must carry a plausible ``range`` and is never reported
  as a fitted or point value — it is swept (``reportable_values`` refuses it).
* Values are loaded from config/dynamics_params.yaml, written by the
  researcher; there are no defaults in code.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from src.common.config import load_config

CONFIDENCE = ("strong", "moderate", "weak")


@dataclass(frozen=True)
class Parameter:
    name: str
    value: float
    units: str
    source: str
    confidence: str
    range: tuple[float, float] | None = None

    def __post_init__(self):
        if self.confidence not in CONFIDENCE:
            raise ValueError(f"{self.name}: confidence must be one of {CONFIDENCE}")
        for attr in ("units", "source"):
            if not str(getattr(self, attr)).strip():
                raise ValueError(f"{self.name}: {attr} must be stated (use 'dimensionless' / a citation)")
        if self.range is not None:
            lo, hi = self.range
            if not lo < hi:
                raise ValueError(f"{self.name}: range must be (low, high) with low < high")
            if not lo <= self.value <= hi:
                raise ValueError(f"{self.name}: value {self.value} outside its plausible range {self.range}")
        if self.confidence == "weak" and self.range is None:
            raise ValueError(f"{self.name}: a weak parameter needs a plausible range to be swept over")


@dataclass
class ParameterSet:
    params: dict[str, Parameter] = field(default_factory=dict)

    @classmethod
    def from_list(cls, items: list[Parameter]) -> "ParameterSet":
        out = cls()
        for p in items:
            if p.name in out.params:
                raise ValueError(f"duplicate parameter {p.name}")
            out.params[p.name] = p
        return out

    @classmethod
    def from_yaml(cls, path: Path) -> "ParameterSet":
        cfg = load_config(path)
        items = []
        for name, spec in cfg.get("parameters").items():
            rng = spec.get("range")
            items.append(Parameter(name=name, value=float(spec["value"]), units=str(spec["units"]),
                                   source=str(spec["source"]), confidence=str(spec["confidence"]),
                                   range=tuple(rng) if rng is not None else None))
        return cls.from_list(items)

    def values(self, overrides: dict[str, float] | None = None) -> dict[str, float]:
        v = {k: p.value for k, p in self.params.items()}
        for k, x in (overrides or {}).items():
            if k not in v:
                raise KeyError(f"unknown parameter {k}")
            v[k] = float(x)
        return v

    def weak(self) -> list[str]:
        return [k for k, p in self.params.items() if p.confidence == "weak"]

    def ranges(self, names: list[str]) -> dict[str, tuple[float, float]]:
        out = {}
        for n in names:
            r = self.params[n].range
            if r is None:
                raise ValueError(f"{n} has no plausible range; it cannot be swept")
            out[n] = r
        return out

    def reportable_values(self) -> dict[str, float]:
        """Point values that may appear in a results table: never weak ones."""
        return {k: p.value for k, p in self.params.items() if p.confidence != "weak"}
