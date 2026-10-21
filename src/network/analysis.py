"""One Module 6 analysis, without file input or output, so that the real run
(src/network/run_network.py) and the calibration on simulated data
(sim/network_calibration.py) execute the same code.

For each seed s, the evaluation regions are the cortical regions minus those
listed in `exclude` (the entorhinal regions, primary analysis) minus the
seed's own cortical regions; for the LC and the entorhinal seeds the two sets
are the same, so their correlations are computed on identical regions. Each
seed's correlation with the observed map is maximised over the diffusion-time
grid; every spun map is maximised over the same grid, so no seed gains from the
maximisation. Spins permute the full cortical map; the evaluation regions are
then read from the spun map.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from src.network.diffusion import (Diffuser, corr_curves, crossover_from_curves, laplacian, max_over_tau,
                                   prediction_curves, seed_vector, standardise_rows)
from src.network.verdict import SeedComparison, spin_p_lower, spin_p_upper, verdict

CHUNK = 2000  # spun maps per block (memory: CHUNK x n_tau floats per seed)


@dataclass
class SeedFit:
    label: str
    r: float
    tau: float
    p_spin: float
    curve: np.ndarray
    null_max: np.ndarray = field(repr=False)
    n_eval: int = 0


@dataclass
class Result:
    fits: dict[str, SeedFit]
    comparison: SeedComparison
    verdict: str
    delta: float


def _seed_fit(diff, x0, cort_nodes, mask, observed, spun, taus, method, label):
    nodes = cort_nodes[mask]
    P = prediction_curves(diff, x0, nodes, taus)
    Zp = standardise_rows(P, method)
    r_curve = corr_curves(Zp, standardise_rows(observed[mask], method))[0]
    r, t = max_over_tau(r_curve[None, :], taus)
    null = np.empty(len(spun))
    for s in range(0, len(spun), CHUNK):
        R = corr_curves(Zp, standardise_rows(spun[s:s + CHUNK][:, mask], method))
        null[s:s + CHUNK] = max_over_tau(R, taus)[0]
    return SeedFit(label, float(r[0]), float(t[0]), spin_p_upper(float(r[0]), null), r_curve, null, int(mask.sum())), P


def analyse(W: np.ndarray, names: list[str], cort_nodes: np.ndarray, observed: np.ndarray, perms: np.ndarray,
            seeds: dict[str, list[str]], exclude: list[str], *, taus: np.ndarray, laplacian_kind: str,
            method: str, crossover_threshold: float, alpha: float, lc: str = "LC",
            ec: str = "entorhinal") -> Result:
    """cort_nodes: indices (into names) of the cortical regions, left then right,
    in the order of `observed` and of the spin permutations."""
    cort_nodes = np.asarray(cort_nodes)
    observed = np.asarray(observed, float)
    if observed.shape != (len(cort_nodes),) or perms.shape[1] != len(cort_nodes):
        raise ValueError("observed map, cortical nodes and spins must have the same length")
    cort_names = [names[i] for i in cort_nodes]
    base = np.array([n not in set(exclude) for n in cort_names])
    diff = Diffuser.from_laplacian(laplacian(W, laplacian_kind))
    spun = observed[perms]
    fits, curves = {}, {}
    for label, regions in seeds.items():
        absent = [r for r in regions if r not in names]
        if absent:
            raise KeyError(f"seed {label}: regions {absent} not in the connectome")
        own = np.array([n in set(regions) for n in cort_names])
        x0 = seed_vector(len(names), [names.index(r) for r in regions])
        fits[label], curves[label] = _seed_fit(diff, x0, cort_nodes, base & ~own, observed, spun, taus,
                                               method, label)
    if fits[lc].n_eval != fits[ec].n_eval:
        raise ValueError("LC and entorhinal seeds must be evaluated on the same regions")
    d_obs = fits[lc].r - fits[ec].r
    d_null = fits[lc].null_max - fits[ec].null_max
    cross = crossover_from_curves(curves[lc], curves[ec], taus, crossover_threshold)
    comp = SeedComparison(r_lc=fits[lc].r, r_ec=fits[ec].r, tau_lc=fits[lc].tau, tau_ec=fits[ec].tau,
                          p_lc=fits[lc].p_spin, p_lc_better=spin_p_upper(d_obs, d_null),
                          p_ec_better=spin_p_lower(d_obs, d_null), crossover=cross)
    return Result(fits, comp, verdict(comp, alpha), d_obs)


def seed_ranking(W: np.ndarray, names: list[str], cort_nodes: np.ndarray, observed: np.ndarray,
                 seeds: dict[str, list[str]], exclude: list[str], *, taus: np.ndarray, laplacian_kind: str,
                 method: str) -> list[tuple[str, float]]:
    """Descriptive: every candidate seed's best r (same evaluation rule as
    analyse, no spins), sorted from best to worst."""
    cort_nodes = np.asarray(cort_nodes)
    cort_names = [names[i] for i in cort_nodes]
    base = np.array([n not in set(exclude) for n in cort_names])
    diff = Diffuser.from_laplacian(laplacian(W, laplacian_kind))
    out = []
    for label, regions in seeds.items():
        own = np.array([n in set(regions) for n in cort_names])
        mask = base & ~own
        x0 = seed_vector(len(names), [names.index(r) for r in regions])
        Zp = standardise_rows(prediction_curves(diff, x0, cort_nodes[mask], taus), method)
        R = corr_curves(Zp, standardise_rows(np.asarray(observed, float)[mask], method))
        out.append((label, float(max_over_tau(R, taus)[0][0])))
    return sorted(out, key=lambda t: -t[1])
