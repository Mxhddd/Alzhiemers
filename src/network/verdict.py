"""The pre-registered verdict for Module 6 (prereg/network.md, section 2).

Kept as a pure function so that the rule can be unit-tested and run on
simulated data (sim/network_calibration.py) exactly as it will run on the
real map.

Inputs, all computed on the same evaluation regions and the same spins:
  r_lc, r_ec       best Pearson (or Spearman) r over the diffusion-time grid
  tau_lc, tau_ec   the diffusion times of those best fits
  p_lc             spin p that the LC seed fits better than the spun maps do
  p_lc_better      spin p that r_lc - r_ec is larger than for the spun maps
  p_ec_better      spin p that r_lc - r_ec is smaller than for the spun maps
  crossover        first tau at which the LC and EC predictions correlate at
                   or above the crossover threshold (None: never in the grid)
"""

from __future__ import annotations

from dataclasses import dataclass

SUPPORTS = "SUPPORTS"
AGAINST = "AGAINST"
AGAINST_BEYOND_NULL = "AGAINST (beyond the spin null)"
INCONCLUSIVE = "INCONCLUSIVE"
UNINFORMATIVE = "UNINFORMATIVE"
LABELS = (SUPPORTS, AGAINST, AGAINST_BEYOND_NULL, INCONCLUSIVE, UNINFORMATIVE)


@dataclass(frozen=True)
class SeedComparison:
    r_lc: float
    r_ec: float
    tau_lc: float
    tau_ec: float
    p_lc: float
    p_lc_better: float
    p_ec_better: float
    crossover: float | None


def spin_p_upper(observed: float, null) -> float:
    """One-sided spin p for 'observed is larger than the null', +1 corrected."""
    import numpy as np
    null = np.asarray(null, float)
    return float((1 + np.sum(null >= observed)) / (1 + len(null)))


def spin_p_lower(observed: float, null) -> float:
    """One-sided spin p for 'observed is smaller than the null', +1 corrected."""
    import numpy as np
    null = np.asarray(null, float)
    return float((1 + np.sum(null <= observed)) / (1 + len(null)))


def verdict(c: SeedComparison, alpha: float) -> str:
    """Rules, in this order:
    1. UNINFORMATIVE if both seeds' best fits lie at or beyond the crossover
       time: the two predictions are then nearly the same pattern.
    2. AGAINST if r_lc <= r_ec (the prediction's direction failed);
       AGAINST (beyond the spin null) if, in addition, p_ec_better < alpha.
    3. SUPPORTS if r_lc > r_ec, p_lc_better < alpha, p_lc < alpha, and the LC
       best fit lies before the crossover time (or the crossover is never
       reached).
    4. INCONCLUSIVE otherwise."""
    beyond = c.crossover is not None
    if beyond and c.tau_lc >= c.crossover and c.tau_ec >= c.crossover:
        return UNINFORMATIVE
    if c.r_lc - c.r_ec <= 0:
        return AGAINST_BEYOND_NULL if c.p_ec_better < alpha else AGAINST
    lc_early = (not beyond) or c.tau_lc < c.crossover
    if c.p_lc_better < alpha and c.p_lc < alpha and lc_early:
        return SUPPORTS
    return INCONCLUSIVE


def direction(label: str) -> str:
    """Collapse a verdict to its direction for the robustness rule."""
    if label == SUPPORTS:
        return "supports"
    if label in (AGAINST, AGAINST_BEYOND_NULL):
        return "against"
    if label in (INCONCLUSIVE, UNINFORMATIVE):
        return "neither"
    raise ValueError(f"unknown verdict {label!r}")


def robust(primary: str, sensitivity: dict[str, str], robustness_set: list[str]) -> bool:
    """The primary verdict is robust only if every analysis in the robustness
    set points the same direction (supports / against / neither). A missing
    analysis (could not be computed) makes the verdict not robust."""
    want = direction(primary)
    return all(name in sensitivity and direction(sensitivity[name]) == want for name in robustness_set)
