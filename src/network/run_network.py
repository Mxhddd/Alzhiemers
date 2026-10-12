"""Module 6 entry point: python -m src.network.run_network  (or `make network`)

Does observed regional tau fit an LC seed better than an entorhinal seed?
Every seed's correlation is maximised over diffusion time; significance comes
ONLY from spin nulls, and the null maximises over diffusion time too, so the
seed with more freedom gets no free advantage. The LC − entorhinal difference
gets its own spin null. The crossover time (after which LC and entorhinal
predictions are indistinguishable) is reported.

Refuses to start unless Module 6's validation passed on the current code and
prereg/network.md is frozen. Real data only; fails loudly while absent.
"""

from __future__ import annotations

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.common.config import load_config
from src.common.io import analysis_input
from src.common.paths import CONFIG_DIR, REPO_ROOT, TABLES_DIR
from src.common.plotting import INK, apply_theme, save_figure, series_colors
from src.common.prereg import require_preregistration
from src.common.runlog import start_run
from src.common.validation import require_validated
from src.network.diffusion import Diffuser, best_fit, crossover_tau, laplacian, seed_vector
from src.network.spin import perm_pvalue, spin_permutations

CONFIG = CONFIG_DIR / "network.yaml"


def load_inputs(cfg):
    W = pd.read_csv(analysis_input(REPO_ROOT / cfg.get("data.connectome"), provided_by="researcher (HCP connectome)"),
                    index_col=0)
    if list(W.index) != list(W.columns):
        raise ValueError("connectome row and column region names differ")
    reg = pd.read_csv(analysis_input(REPO_ROOT / cfg.get("data.regions"), provided_by="researcher (parcellation table)"))
    need = ["region", "hemisphere", "cortical", "sphere_x", "sphere_y", "sphere_z"]
    if list(reg.columns) != need:
        raise ValueError(f"regions.csv columns must be exactly {need}")
    reg = reg.set_index("region").loc[W.index]  # KeyError if a connectome region is missing
    tau = pd.read_csv(analysis_input(REPO_ROOT / cfg.get("data.tau"), provided_by="researcher (ADNI tau-PET)"))
    if list(tau.columns) != ["region", "value"] or tau["value"].isna().any():
        raise ValueError("tau_regional.csv must have columns region,value and no missing values")
    return W, reg, tau.set_index("region")["value"]


def main() -> int:
    require_validated("network")
    prereg = require_preregistration("network")
    cfg = load_config(CONFIG)
    with start_run(module="network", config_path=CONFIG, seed=cfg.get("seed"),
                   packages=["numpy", "scipy", "pandas"], preregistration=prereg) as run:
        TABLES_DIR.mkdir(parents=True, exist_ok=True)
        W, reg, tau = load_inputs(cfg)
        names = list(W.index)
        cort = reg["cortical"].astype(bool).to_numpy()
        lh = np.where(cort & (reg["hemisphere"] == "L").to_numpy())[0]
        rh = np.where(cort & (reg["hemisphere"] == "R").to_numpy())[0]
        ev = np.concatenate([lh, rh])
        missing = [names[i] for i in ev if names[i] not in tau.index]
        if missing:
            raise KeyError(f"tau values missing for cortical regions {missing[:10]} (no imputation)")
        observed = tau.loc[[names[i] for i in ev]].to_numpy(float)
        xyz = reg[["sphere_x", "sphere_y", "sphere_z"]].to_numpy(float)
        if np.isnan(xyz[ev]).any():
            raise ValueError("cortical regions need sphere coordinates for the spin test")

        diff = Diffuser.from_laplacian(laplacian(W.to_numpy(float), cfg.get("model.laplacian")))
        taus = np.geomspace(cfg.get("model.tau_min"), cfg.get("model.tau_max"), cfg.get("model.n_tau"))
        seeds = {"LC": cfg.get("seeds.lc"), "entorhinal": cfg.get("seeds.entorhinal"), **cfg.get("seeds.controls")}
        x0 = {}
        for label, regions in seeds.items():
            absent = [r for r in regions if r not in names]
            if absent:
                raise KeyError(f"seed {label}: regions {absent} not in the connectome")
            x0[label] = seed_vector(len(names), [names.index(r) for r in regions])
        run.info("uniform growth (alpha) cannot change spatial correlations and is not fitted")

        perms = spin_permutations(xyz[lh], xyz[rh], cfg.get("spin.n_perm"), seed=cfg.get("seed"),
                                  method=cfg.get("spin.method"))
        fits = {s: best_fit(diff, x, observed, ev, taus) for s, x in x0.items()}
        null_max = {s: np.empty(len(perms)) for s in x0}
        for k, p in enumerate(perms):
            spun = observed[p]
            for s, x in x0.items():
                null_max[s][k] = best_fit(diff, x, spun, ev, taus)[0]
        rows = []
        for s, (r, t_best, curve) in fits.items():
            p_spin = float((1 + np.sum(null_max[s] >= r)) / (1 + len(perms)))  # one-sided: better fit than null maxima
            rows.append({"seed": s, "best_r": r, "best_tau": t_best, "spin_p": p_spin, "n_perm": len(perms)})
            run.info(f"seed {s}: max r = {r:.3f} at tau {t_best:.3g}, spin p = {p_spin:.4f}")
        d_obs = fits["LC"][0] - fits["entorhinal"][0]
        d_null = null_max["LC"] - null_max["entorhinal"]
        p_diff = perm_pvalue(d_obs, d_null - d_null.mean())
        cross = crossover_tau(diff, x0["LC"], x0["entorhinal"], ev, taus, cfg.get("model.crossover_threshold"))
        run.info(f"LC - entorhinal: delta r = {d_obs:.3f}, spin p (two-sided, centred null) = {p_diff:.4f}; "
                 f"crossover tau = {cross if cross is not None else 'not reached in grid'}")
        for s in ("LC", "entorhinal"):
            if cross is not None and fits[s][1] >= cross:
                run.warning(f"{s} best fit lies beyond the crossover time: at this stage the seeds are not distinguishable")
        res = pd.DataFrame(rows)
        res["lc_minus_entorhinal_r"] = d_obs
        res["lc_minus_entorhinal_spin_p"] = p_diff
        res["crossover_tau"] = cross
        res.to_csv(TABLES_DIR / "network_seed_comparison.csv", index=False)
        run.record_output(TABLES_DIR / "network_seed_comparison.csv")

        apply_theme()
        fig, ax = plt.subplots(figsize=(6, 4))
        for (s, (_, _, curve)), c in zip(fits.items(), series_colors(len(fits))):
            ax.plot(taus, curve, color=c, label=s)
        if cross is not None:
            ax.axvline(cross, color=INK["muted"], ls=":", lw=1)
            ax.text(cross, ax.get_ylim()[0], " crossover", fontsize=7, color=INK["secondary"], va="bottom")
        ax.set_xscale("log")
        ax.set_xlabel("diffusion time tau = beta t")
        ax.set_ylabel("Pearson r with observed regional tau")
        ax.legend()
        save_figure(fig, "network_seed_fits", run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
