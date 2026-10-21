"""Module 6 entry point: python -m src.network.run_network  (or `make network`)

Does the earliest-stage regional tau map fit diffusion from an LC seed better
than diffusion from an entorhinal seed? (prereg/network.md)

* Every seed's correlation is maximised over diffusion time; significance comes
  ONLY from spin nulls, and every spun map is maximised over the same grid.
* The entorhinal regions (and every seed's own regions) are left out of the
  correlation, so the entorhinal seed is not scored on its own value.
* The pre-registered verdict (src/network/verdict.py), the crossover time, the
  control seeds, the descriptive ranking of every bilateral seed, and the nine
  pre-specified sensitivity analyses with the robustness rule.

Refuses to start unless Module 6's validation passed on the current code and
prereg/network.md is frozen. Inputs are group-level files (section 3 of the
pre-registration); fails loudly while they are absent.
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
from src.network.analysis import analyse, seed_ranking
from src.network.connectome import uniform_projection
from src.network.spin import spin_permutations
from src.network.verdict import robust

CONFIG = CONFIG_DIR / "network.yaml"
REGION_COLUMNS = ["region", "hemisphere", "cortical", "sphere_x", "sphere_y", "sphere_z"]


def load_connectome(path) -> pd.DataFrame:
    W = pd.read_csv(analysis_input(REPO_ROOT / path, provided_by="researcher (HCP group connectome)"), index_col=0)
    if list(W.index) != list(W.columns):
        raise ValueError(f"{path}: row and column region names differ")
    return W


def load_regions(path, names) -> pd.DataFrame:
    reg = pd.read_csv(analysis_input(REPO_ROOT / path, provided_by="researcher (node table)"))
    if list(reg.columns) != REGION_COLUMNS:
        raise ValueError(f"regions.csv columns must be exactly {REGION_COLUMNS}")
    return reg.set_index("region").loc[names]  # KeyError if a connectome region is missing


def load_tau(path) -> pd.Series:
    tau = pd.read_csv(analysis_input(REPO_ROOT / path, provided_by="researcher (ADNI group tau map, adni_tau.py)"))
    if list(tau.columns) != ["region", "value"] or tau["value"].isna().any():
        raise ValueError(f"{path}: columns must be region,value with no missing values")
    return tau.set_index("region")["value"]


def bilateral_seeds(names: list[str]) -> dict[str, list[str]]:
    """Every left/right pair of nodes (ctx-lh-x / ctx-rh-x, Left-x / Right-x, LC-lh / LC-rh)."""
    out = {}
    for n in names:
        for a, b in (("ctx-lh-", "ctx-rh-"), ("Left-", "Right-")):
            if n.startswith(a) and n.replace(a, b, 1) in names:
                out[n[len(a):]] = [n, n.replace(a, b, 1)]
        if n.endswith("-lh") and n[:-3] + "-rh" in names:
            out[n[:-3]] = [n, n[:-3] + "-rh"]
    return out


def main() -> int:
    require_validated("network")
    prereg = require_preregistration("network")
    cfg = load_config(CONFIG)
    with start_run(module="network", config_path=CONFIG, seed=cfg.get("seed"),
                   packages=["numpy", "scipy", "pandas"], preregistration=prereg) as run:
        TABLES_DIR.mkdir(parents=True, exist_ok=True)
        taus = np.geomspace(cfg.get("model.tau_min"), cfg.get("model.tau_max"), cfg.get("model.n_tau"))
        seeds = {"LC": cfg.get("seeds.lc"), "entorhinal": cfg.get("seeds.entorhinal"), **cfg.get("seeds.controls")}
        run.info("uniform growth (alpha) cannot change spatial correlations and is not fitted")

        primary = {"connectome": "connectome", "tau": "tau", "laplacian": cfg.get("model.laplacian"),
                   "correlation": cfg.get("evaluation.correlation"), "spin_method": cfg.get("spin.method"),
                   "exclude": list(cfg.get("evaluation.exclude")), "lc_projection": "tractography"}
        analyses = {"primary": primary}
        for name, change in cfg.get("sensitivity").items():
            unknown = set(change) - set(primary)
            if unknown:
                raise KeyError(f"sensitivity {name}: unknown settings {unknown}")
            analyses[name] = {**primary, **change}

        perms_cache, rows, verdicts, results = {}, [], {}, {}
        for name, spec in analyses.items():
            Wdf = load_connectome(cfg.get(f"data.{spec['connectome']}"))
            names = list(Wdf.index)
            reg = load_regions(cfg.get("data.regions"), names)
            cort = reg["cortical"].astype(bool).to_numpy()
            lh = np.where(cort & (reg["hemisphere"] == "L").to_numpy())[0]
            rh = np.where(cort & (reg["hemisphere"] == "R").to_numpy())[0]
            cort_nodes = np.concatenate([lh, rh])
            W = Wdf.to_numpy(float)
            if spec["lc_projection"] == "uniform_ipsilateral_cortex":
                for lc_name in cfg.get("seeds.lc"):
                    side = reg.loc[lc_name, "hemisphere"]
                    W = uniform_projection(W, names.index(lc_name), list(lh if side == "L" else rh))
            elif spec["lc_projection"] != "tractography":
                raise ValueError(f"unknown lc_projection {spec['lc_projection']!r}")
            tau = load_tau(cfg.get(f"data.{spec['tau']}"))
            missing = [names[i] for i in cort_nodes if names[i] not in tau.index]
            if missing:
                raise KeyError(f"tau values missing for cortical regions {missing[:10]} (no imputation)")
            observed = tau.loc[[names[i] for i in cort_nodes]].to_numpy(float)
            xyz = reg[["sphere_x", "sphere_y", "sphere_z"]].to_numpy(float)
            if np.isnan(xyz[cort_nodes]).any():
                raise ValueError("cortical regions need sphere coordinates for the spin test")
            key = (spec["spin_method"], tuple(cort_nodes))
            if key not in perms_cache:
                perms_cache[key] = spin_permutations(xyz[lh], xyz[rh], cfg.get("spin.n_perm"), seed=cfg.get("seed"),
                                                     method=spec["spin_method"])
            res = analyse(W, names, cort_nodes, observed, perms_cache[key], seeds, spec["exclude"], taus=taus,
                          laplacian_kind=spec["laplacian"], method=spec["correlation"],
                          crossover_threshold=cfg.get("model.crossover_threshold"), alpha=cfg.get("alpha"))
            results[name], verdicts[name] = res, res.verdict
            c = res.comparison
            for s, f in res.fits.items():
                rows.append({"analysis": name, "seed": s, "best_r": f.r, "best_tau": f.tau, "spin_p": f.p_spin,
                             "n_eval_regions": f.n_eval, "n_perm": len(perms_cache[key])})
            rows.append({"analysis": name, "seed": "LC minus entorhinal", "best_r": res.delta,
                         "spin_p_lc_better": c.p_lc_better, "spin_p_ec_better": c.p_ec_better,
                         "crossover_tau": c.crossover, "verdict": res.verdict})
            run.info(f"[{name}] r_LC {c.r_lc:.3f} (tau {c.tau_lc:.3g}, p {c.p_lc:.4f}); r_EC {c.r_ec:.3f} "
                     f"(tau {c.tau_ec:.3g}); delta {res.delta:.3f}, p(LC better) {c.p_lc_better:.4f}, "
                     f"p(EC better) {c.p_ec_better:.4f}; crossover {c.crossover}; verdict {res.verdict}")
            if name == "primary":
                ranking = seed_ranking(W, names, cort_nodes, observed, bilateral_seeds(names), spec["exclude"],
                                       taus=taus, laplacian_kind=spec["laplacian"], method=spec["correlation"])
                pd.DataFrame(ranking, columns=["seed", "best_r"]).assign(rank=lambda d: np.arange(1, len(d) + 1)) \
                    .to_csv(TABLES_DIR / "network_seed_ranking.csv", index=False)
                run.record_output(TABLES_DIR / "network_seed_ranking.csv")
                lc_rank = [s for s, _ in ranking].index("LC") + 1 if "LC" in dict(ranking) else None
                run.info(f"descriptive: LC ranks {lc_rank} of {len(ranking)} bilateral seeds")

        is_robust = robust(verdicts["primary"], verdicts, list(cfg.get("robustness_set")))
        run.info(f"PRIMARY VERDICT: {verdicts['primary']}; robust across {list(cfg.get('robustness_set'))}: {is_robust}")
        out = pd.DataFrame(rows)
        out["primary_verdict_robust"] = is_robust
        out.to_csv(TABLES_DIR / "network_seed_comparison.csv", index=False)
        run.record_output(TABLES_DIR / "network_seed_comparison.csv")

        apply_theme()
        res = results["primary"]
        fig, ax = plt.subplots(figsize=(6, 4))
        for (s, f), col in zip(res.fits.items(), series_colors(len(res.fits))):
            ax.plot(taus, f.curve, color=col, label=s)
        if res.comparison.crossover is not None:
            ax.axvline(res.comparison.crossover, color=INK["muted"], ls=":", lw=1)
            ax.text(res.comparison.crossover, ax.get_ylim()[0], " crossover", fontsize=7, color=INK["secondary"],
                    va="bottom")
        ax.set_xscale("log")
        ax.set_xlabel("diffusion time tau = beta t")
        ax.set_ylabel("correlation with observed regional tau")
        ax.legend()
        save_figure(fig, "network_seed_fits", run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
