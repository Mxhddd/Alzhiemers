"""Module 3 entry point: python -m src.mediation.run_mediation

Real-data path ONLY (the RADC ROS/MAP Standard dataset, under data/raw/radc/,
requested and not yet received). Refuses to start unless Module 3's
validation passed on the current code, prereg/mediation.md is frozen, and the
researcher has confirmed every variable name against the delivered codebook
(config data.codebook_confirmed). Simulated data are used only in tests and sim/.

What runs is fixed by config/mediation.yaml and prereg/mediation.md:
H1, H2 (amyloid), H3 (four social / sensory exposures), the pre-specified
sensitivity analyses S1-S10, and -- only if config lc.variables names an LC
measure -- H4, the mediation through the LC (Killer 2).
"""

from __future__ import annotations

import json
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

from src.common.config import load_config
from src.common.paths import CONFIG_DIR, TABLES_DIR
from src.common.plotting import INK, apply_theme, save_figure
from src.common.prereg import require_preregistration
from src.common.runlog import start_run
from src.common.validation import require_validated
from src.mediation.hypotheses import (association, autopsy_weights, mediation_verdict,
                                      multiply_imputed_association, verdict)
from src.mediation.mediation import bootstrap as linear_bootstrap
from src.mediation.mediation import estimate as linear_estimate
from src.mediation.ordinal import bootstrap_ordinal_mediation, fit_ordinal, fit_ordinal_mediation
from src.mediation.radc import load_radc

CONFIG = CONFIG_DIR / "mediation.yaml"


class Analysis:
    def __init__(self, cfg, run, people: pd.DataFrame):
        self.cfg, self.run, self.people = cfg, run, people
        v = cfg.get("variables")
        self.Y, self.AMY, self.AMY2, self.Y2 = v["outcome"], v["amyloid"], v["amyloid_sensitivity"], v["outcome_sensitivity"]
        self.covs = list(cfg.get("model.covariates"))
        self.aut = people[people["has_pathology"]].copy()
        self.rows: list[dict] = []
        self.kw = dict(equivalence=cfg.get("model.equivalence_log_or"), max_width=cfg.get("model.uninformative_log_width"),
                       level=cfg.get("model.ci_level"))

    # -- helpers -------------------------------------------------------------
    def complete(self, df: pd.DataFrame, cols: list[str], label: str) -> pd.DataFrame:
        bad = df[cols].isna().any(axis=1)
        for i in df.index[bad]:
            self.run.drop(stage=label, item=f"participant {i}",
                          reason="missing: " + ",".join(c for c in cols if pd.isna(df.at[i, c])))
        return df[~bad]

    def record(self, test: str, analysis: str, res: dict, m: int, *, confirmatory: bool, note: str = "") -> dict:
        ver = verdict(res["log_or"], res["se"], m=m, **self.kw)
        if analysis == "descriptive":
            ver["verdict"] = "descriptive: two-sided estimate, no verdict"
        row = {"test": test, "analysis": analysis, "confirmatory": confirmatory, "n": res["n"],
               "scale": res.get("scale", "log OR per SD (exp = odds ratio)"),
               "log_or": res["log_or"], "se": res["se"], "or": ver["or"], "ci95_lower": ver["ci95"][0],
               "ci95_upper": ver["ci95"][1], "adj_level": ver["adj_level"], "ciadj_lower": ver["ci_adj"][0],
               "ciadj_upper": ver["ci_adj"][1], "p": ver["p"], "verdict": ver["verdict"], "e_value": ver["e_value"],
               "e_value_ci": ver["e_value_ci"], "levels": json.dumps(res.get("levels", {})),
               "dropped_constant": ",".join(res.get("dropped_constant", [])), "note": note}
        self.rows.append(row)
        self.run.info(f"{test} [{analysis}] n={res['n']} OR={ver['or']:.3f} adj CI {ver['ci_adj'][0]:.3f}-"
                      f"{ver['ci_adj'][1]:.3f}: {ver['verdict']}")
        return row

    def expo(self, name: str) -> dict:
        return next(e for e in self.cfg.get("exposures") if e["name"] == name)

    def expo_sample(self, name: str, timing: str = "earliest") -> pd.DataFrame:
        e = self.expo(name)
        d = self.aut[self.aut[self.cfg.get("data.columns.study")].astype(str).isin(e["cohorts"])].copy()
        d["_x"] = d[f"{name}__{timing}"]
        return d

    # -- Family A ------------------------------------------------------------
    def amyloid_tests(self):
        cols = [self.Y, self.AMY] + self.covs
        d = self.complete(self.aut, cols, "H1/H2")
        entry = self.cfg.get("model.entry_split")
        for test, split in (("H1", None), ("H2", entry)):
            r = association(d, outcome=self.Y, exposure=self.AMY, direction=1, covariates=self.covs, split=split)
            self.record(test, "primary", r, 2, confirmatory=True)
            # S1 Thal-based A score instead of amyloid density
            d1 = self.complete(self.aut, [self.Y, self.AMY2] + self.covs, f"{test}-S1")
            self.record(test, "S1", association(d1, outcome=self.Y, exposure=self.AMY2, direction=1,
                                                covariates=self.covs, split=split), 2, confirmatory=False)
            if split is None:   # S2 collapsed Braak score (H1 only; H2's split is the same in both codings)
                d2 = self.complete(self.aut, [self.Y2, self.AMY] + self.covs, "H1-S2")
                self.record(test, "S2", association(d2, outcome=self.Y2, exposure=self.AMY, direction=1,
                                                    covariates=self.covs), 2, confirmatory=False)
                self.record(test, "S3", self.linear_assoc(d, self.AMY, 1), 2, confirmatory=False)
            self.cohort_split(test, d, self.AMY, 1, split)
            self.autopsy_weighted(test, self.AMY, 1, split)
        self.amyloid_descriptives(d)

    def amyloid_descriptives(self, d: pd.DataFrame):
        """Secondary, two-sided, no verdict: the age association given amyloid,
        the amyloid x age interaction, and age in amyloid-free brains."""
        age = self.cfg.get("variables.age_death")
        others = [c for c in self.covs if c != age] + [self.AMY]
        self.record("D1 age | amyloid", "descriptive",
                    association(d, outcome=self.Y, exposure=age, direction=1, covariates=others), 1, confirmatory=False,
                    note="log OR per SD of age at death")
        dz = d.copy()
        dz["_amy_x_age"] = ((dz[self.AMY] - dz[self.AMY].mean()) / dz[self.AMY].std()) * \
                           ((dz[age] - dz[age].mean()) / dz[age].std())
        self.record("D2 amyloid x age", "descriptive",
                    association(dz, outcome=self.Y, exposure="_amy_x_age", direction=1,
                                covariates=self.covs + [self.AMY]), 1, confirmatory=False,
                    note="interaction term (product of z-scores), scaled per SD of the product")
        neg = self.aut[pd.to_numeric(self.aut[self.AMY2], errors="coerce") == 0]
        neg = self.complete(neg, [self.Y, age] + [c for c in self.covs if c != age], "D3")
        if len(neg) >= 30 and neg[self.Y].nunique() > 1:
            self.record("D3 age | amyloid-free (A score 0)", "descriptive",
                        association(neg, outcome=self.Y, exposure=age, direction=1,
                                    covariates=[c for c in self.covs if c != age]), 1, confirmatory=False)
        else:
            self.run.info(f"D3 not estimated: {len(neg)} amyloid-free brains with complete data (needs 30)")

    # -- Family B ------------------------------------------------------------
    def social_tests(self):
        fam = self.cfg.get("families.social_sensory")
        m = len(fam)
        normal = self.cfg.get("variables.normal_cognition_code")
        lag = self.cfg.get("exposure_timing.min_lag_years")
        age = self.cfg.get("variables.age_death")
        for name in fam:
            e = self.expo(name)
            dirn = e["risk_direction"]
            base = self.expo_sample(name)
            d = self.complete(base, [self.Y, "_x"] + self.covs, f"H3 {name}")
            self.record(f"H3 {name}", "primary", association(d, outcome=self.Y, exposure="_x", direction=dirn,
                                                             covariates=self.covs), m, confirmatory=True)
            d2 = self.complete(base, [self.Y2, "_x"] + self.covs, f"H3 {name} S2")
            self.record(f"H3 {name}", "S2", association(d2, outcome=self.Y2, exposure="_x", direction=dirn,
                                                        covariates=self.covs), m, confirmatory=False)
            self.record(f"H3 {name}", "S3", self.linear_assoc(d, "_x", dirn), m, confirmatory=False)
            for timing in self.cfg.get("exposure_timing.sensitivity"):
                dt = self.complete(self.expo_sample(name, timing), [self.Y, "_x"] + self.covs, f"H3 {name} S4-{timing}")
                self.record(f"H3 {name}", f"S4-{timing}", association(dt, outcome=self.Y, exposure="_x", direction=dirn,
                                                                      covariates=self.covs), m, confirmatory=False)
            d5 = d[pd.to_numeric(d[f"{name}__dx"], errors="coerce") == normal]
            self.record(f"H3 {name}", "S5", association(d5, outcome=self.Y, exposure="_x", direction=dirn,
                                                        covariates=self.covs), m, confirmatory=False,
                        note="exposure measured at a visit with no cognitive impairment")
            d6 = d[(d[age] - d[f"{name}__age"]) >= lag]
            self.record(f"H3 {name}", "S6", association(d6, outcome=self.Y, exposure="_x", direction=dirn,
                                                        covariates=self.covs), m, confirmatory=False,
                        note=f"exposure measured >= {lag} years before death")
            d7 = self.complete(base, [self.Y, "_x", self.AMY] + self.covs, f"H3 {name} S7")
            self.record(f"H3 {name}", "S7", association(d7, outcome=self.Y, exposure="_x", direction=dirn,
                                                        covariates=self.covs + [self.AMY]), m, confirmatory=False,
                        note="additionally adjusted for amyloid (association not through amyloid)")
            self.maybe_impute(f"H3 {name}", base, "_x", dirn, m)
            self.cohort_split(f"H3 {name}", d, "_x", dirn, None, m=m)
            self.autopsy_weighted(f"H3 {name}", "_x", dirn, None, m=m, exposure=name)

    # -- shared sensitivity analyses ------------------------------------------
    def linear_assoc(self, d: pd.DataFrame, x: str, dirn: int) -> dict:
        """S3: Braak stage as a number 0-6, OLS with HC3 errors (coefficient per SD)."""
        z = dirn * (d[x] - d[x].mean()) / d[x].std()
        keep = [c for c in self.covs if d[c].nunique() > 1]
        X = sm.add_constant(np.column_stack([z, d[keep].to_numpy(float)]), has_constant="add")
        f = sm.OLS(d[self.Y].to_numpy(float), X).fit(cov_type="HC3")
        return {"log_or": float(f.params[1]), "se": float(f.bse[1]), "n": int(len(d)),
                "scale": "Braak stages per SD (linear; 'or' columns hold exp(coefficient))",
                "levels": {}, "dropped_constant": [c for c in self.covs if c not in keep]}

    def cohort_split(self, test, d, x, dirn, split, m=2):
        study = self.cfg.get("data.columns.study")
        for c in self.cfg.get("data.cohorts"):
            dc = d[d[study].astype(str) == c]
            if len(dc) < 30 or dc[x].nunique() < 2:
                continue
            self.record(test, f"S10-{c}", association(dc, outcome=self.Y, exposure=x, direction=dirn,
                                                      covariates=self.covs, split=split), m, confirmatory=False)

    def maybe_impute(self, test, base, x, dirn, m):
        cols = [self.Y, x] + self.covs
        share = float(base[cols].isna().any(axis=1).mean())
        if share <= self.cfg.get("missing.mi_threshold_share"):
            self.run.info(f"{test}: S8 not triggered ({share:.1%} incomplete)")
            return
        r = multiply_imputed_association(base, columns=cols, n_imputations=self.cfg.get("missing.mi_imputations"),
                                         cycles=self.cfg.get("missing.mi_cycles"),
                                         seed=self.cfg.get("seed"), outcome=self.Y, exposure=x, direction=dirn,
                                         covariates=self.covs)
        self.record(test, "S8", r, m, confirmatory=False, note=f"{share:.1%} incomplete; multiple imputation")

    def autopsy_weighted(self, test, x, dirn, split, m=2, exposure=None):
        """S9: inverse-probability-of-autopsy weights among deceased participants
        (baseline predictors; for a social / sensory exposure the exposure itself
        is a predictor too, since it is measured in life). Bootstrap interval."""
        cfgw = self.cfg.get("autopsy_weights")
        died = self.cfg.get("variables.died")
        study = self.cfg.get("data.columns.study")
        dead = self.people[self.people[died] == 1].copy()
        preds = list(cfgw["predictors"])
        if exposure is not None:
            dead = dead[dead[study].astype(str).isin(self.expo(exposure)["cohorts"])]
            dead["_x"] = dead[f"{exposure}__earliest"]
            preds.append("_x")
        preds = [c for c in preds if dead[c].nunique() > 1]   # e.g. the cohort indicator in MAP-only analyses
        dead = self.complete(dead, preds, f"{test} S9 weights")
        rng = np.random.default_rng(self.cfg.get("seed"))
        est = []
        point = None
        for b in range(cfgw["n_boot"] + 1):
            db = dead if b == 0 else dead.iloc[rng.integers(0, len(dead), len(dead))]
            try:
                w = autopsy_weights(db.reset_index(drop=True), predictors=preds, observed="has_pathology",
                                    truncate_q=cfgw["truncate_quantile"])
                dd = db.reset_index(drop=True).loc[w.index]
                ok = dd[[self.Y, x] + self.covs].notna().all(axis=1)
                r = association(dd[ok], outcome=self.Y, exposure=x, direction=dirn, covariates=self.covs,
                                split=split, weights=w[ok].to_numpy())
            except Exception as exc:   # a degenerate resample is counted, never hidden
                if b == 0:
                    raise
                self.run.info(f"{test} S9 bootstrap replicate {b} failed: {exc}")
                continue
            if b == 0:
                point = r
            else:
                est.append(r["log_or"])
        se = float(np.std(est, ddof=1))
        self.record(test, "S9", {"log_or": point["log_or"], "se": se, "n": point["n"], "levels": point["levels"],
                                 "dropped_constant": point["dropped_constant"]}, m, confirmatory=False,
                    note=f"autopsy-weighted; bootstrap SE from {len(est)} replicates")

    # -- robustness ----------------------------------------------------------
    def robustness(self) -> pd.DataFrame:
        out = pd.DataFrame(self.rows)
        rset = self.cfg.get("sensitivity.robustness_set")
        prim = out[out.confirmatory]
        rob = []
        for _, p in prim.iterrows():
            s = out[(out.test == p.test) & out.analysis.isin(rset)]
            if p.verdict == "supports":
                ok = all(_excludes_null(r, +1) for _, r in s.iterrows())
            elif p.verdict == "AGAINST":
                ok = all(_excludes_null(r, -1) for _, r in s.iterrows())
            elif p.verdict.startswith("AGAINST (practically"):
                ok = all(not _excludes_null(r, +1) and not _excludes_null(r, -1) for _, r in s.iterrows())
            else:
                ok = None
            rob.append({"test": p.test, "verdict": p.verdict, "robust": ok, "analyses_checked": ",".join(s.analysis)})
        return pd.DataFrame(rob)

    # -- H4: the conditional LC branch ----------------------------------------
    def lc_branch(self) -> pd.DataFrame | None:
        lcs = self.cfg.get("lc.variables")
        if not lcs:
            self.run.info("H4 / Killer 2 NOT RUN: config lc.variables is empty (no LC measure in the delivered "
                          "data). This dataset cannot test mediation through the LC.")
            return None
        lc = lcs[0]   # the researcher lists them in the pre-registered priority order
        name, sign = lc["name"], lc["damage_direction"]
        mcfg = self.cfg.get("mediation")
        grid = np.round(np.arange(mcfg["rho_min"], mcfg["rho_max"] + 1e-9, mcfg["rho_step"]), 10)
        rows = []
        for k, ex in enumerate(self.cfg.get("families.social_sensory")):
            e = self.expo(ex)
            base = self.expo_sample(ex)
            d = self.complete(base, [self.Y, "_x", name] + self.covs, f"H4 {ex}")
            if len(d) < self.cfg.get("lc.min_n"):
                self.run.info(f"H4 {ex}: {len(d)} participants with the LC measure (< {self.cfg.get('lc.min_n')}); "
                              "descriptive only, no verdict")
                rows.append({"exposure": ex, "n": len(d), "verdict": "not estimated (too few)"})
                continue
            T = e["risk_direction"] * (d["_x"] - d["_x"].mean()) / d["_x"].std()
            M = sign * (d[name] - d[name].mean()) / d[name].std()
            keep = [c for c in self.covs if d[c].nunique() > 1]
            X = d[keep].to_numpy(float)
            q0, q1 = np.quantile(T, mcfg["contrast_quantiles"])
            fit = fit_ordinal_mediation(T, M, d[self.Y], X)
            eff = fit.effects(q1, q0)
            boot = bootstrap_ordinal_mediation(T, M, d[self.Y], X, t1=q1, t0=q0, n_boot=mcfg["n_boot"],
                                               seed=self.cfg.get("seed") + k, rho_grid=grid)
            mv = mediation_verdict(eff["nie"], boot["nie"], m=mcfg["family_size"], level=self.kw["level"],
                                   predicted_sign=self.cfg.get("killer2.predicted_sign"))
            rho_star = fit.rho_star
            fragile = mv["verdict"].startswith("NIE in the predicted") and abs(rho_star) < mcfg["fragile_abs_rho"]
            curve_lo, curve_hi = np.quantile(boot["nie_curve"], [0.025, 0.975], axis=0)
            inc0 = (curve_lo <= 0) & (curve_hi >= 0)
            lin = linear_estimate(T.to_numpy(), M.to_numpy(), d[self.Y].to_numpy(float), X)
            lb = linear_bootstrap(T.to_numpy(), M.to_numpy(), d[self.Y].to_numpy(float), X,
                                  n_boot=self.cfg.get("linear.n_boot"), seed=self.cfg.get("seed") + 100 + k)
            la = np.quantile([b.acme for b in lb], [0.025, 0.975])
            rows.append({"exposure": ex, "lc_variable": name, "n": len(d), **{k2: eff[k2] for k2 in eff},
                         "nie_ci95_lower": mv["ci95"][0], "nie_ci95_upper": mv["ci95"][1],
                         "nie_ciadj_lower": mv["ci_adj"][0], "nie_ciadj_upper": mv["ci_adj"][1],
                         "verdict": mv["verdict"] + (" (FRAGILE: |rho*| < %.2f)" % mcfg["fragile_abs_rho"] if fragile else ""),
                         "rho_star": rho_star,
                         "smallest_abs_rho_ci_includes_0": float(grid[inc0][np.argmin(np.abs(grid[inc0]))]) if inc0.any() else None,
                         "S3_linear_acme": lin.acme, "S3_linear_acme_ci95": f"{la[0]:.4f} to {la[1]:.4f}",
                         "S3_linear_rho_tilde": lin.rho_tilde})
            self.run.info(f"H4 {ex}: NIE {eff['nie']:.4f} adj CI {mv['ci_adj'][0]:.4f}-{mv['ci_adj'][1]:.4f}; "
                          f"rho* {rho_star:.3f}: {rows[-1]['verdict']}")
        # secondary: does the LC measure's association with Braak depend on amyloid? (gating; two-sided)
        d = self.complete(self.aut, [self.Y, self.AMY, name] + self.covs, "H4 gating")
        if len(d) >= self.cfg.get("lc.min_n"):
            z = lambda s: (s - s.mean()) / s.std()   # noqa: E731
            dz = d.assign(_lc=sign * z(d[name]), _amy=z(d[self.AMY]))
            dz["_lc_x_amy"] = dz["_lc"] * dz["_amy"]
            self.record("D4 LC x amyloid", "descriptive",
                        association(dz, outcome=self.Y, exposure="_lc_x_amy", direction=1,
                                    covariates=self.covs + ["_lc", "_amy"]), 1, confirmatory=False)
        res = pd.DataFrame(rows)
        estimable = res[~res.verdict.str.startswith("not estimated")]
        overall = ("Killer 2 MET: no exposure shows mediation through the LC" if len(estimable) and
                   estimable.verdict.str.contains("Killer 2 MET").all() else
                   "Killer 2 NOT met: at least one exposure shows a predicted-direction NIE" if len(estimable) else
                   "Killer 2 not evaluable: too few participants with the LC measure")
        self.run.info(overall)
        res.attrs["overall"] = overall
        return res


def _excludes_null(row, side: int) -> bool:
    """Adjusted interval on the predicted (+1) or opposite (-1) side of the null.
    Every row stores exp(estimate) and exp(interval), so the null is 1 for odds
    ratios and for the S3 linear coefficient alike (exp(0) = 1)."""
    return bool(row.ciadj_lower > 1) if side > 0 else bool(row.ciadj_upper < 1)


def forest(rows: pd.DataFrame, run):
    apply_theme()
    p = rows[rows.analysis == "primary"].iloc[::-1]
    fig, ax = plt.subplots(figsize=(6.5, 0.45 * len(p) + 1.2))
    y = np.arange(len(p))
    ax.errorbar(p["or"], y, xerr=[p["or"] - p.ciadj_lower, p.ciadj_upper - p["or"]], fmt="o",
                color=INK["primary"], capsize=3)
    ax.axvline(1, color=INK["axis"], lw=0.8)
    ax.set_xscale("log")
    ax.set_yticks(y, [f"{t}  [{v}]" for t, v in zip(p.test, p.verdict)], fontsize=7)
    ax.set_xlabel("odds ratio per SD, oriented to the predicted risk side (Bonferroni-adjusted intervals)")
    ax.set_title("Module 3 confirmatory tests (right of 1 = predicted direction)", fontsize=8, loc="left")
    save_figure(fig, "mediation_confirmatory_forest", run)


def main() -> int:
    require_validated("mediation")
    prereg = require_preregistration("mediation")
    cfg = load_config(CONFIG)
    if cfg.get("data.codebook_confirmed") is not True:
        raise RuntimeError("config data.codebook_confirmed is not true: check every variable name in "
                           "config/mediation.yaml against the delivered RADC codebook (names only) first")
    with start_run(module="mediation", config_path=CONFIG, seed=cfg.get("seed"),
                   packages=["numpy", "pandas", "scipy", "statsmodels"], preregistration=prereg) as run:
        TABLES_DIR.mkdir(parents=True, exist_ok=True)
        people, info = load_radc(cfg)
        run.info(f"RADC participants: {info}")
        a = Analysis(cfg, run, people)
        a.amyloid_tests()
        a.social_tests()
        rows = pd.DataFrame(a.rows)
        rows.to_csv(TABLES_DIR / "mediation_tests.csv", index=False)
        run.record_output(TABLES_DIR / "mediation_tests.csv")
        rob = a.robustness()
        rob.to_csv(TABLES_DIR / "mediation_robustness.csv", index=False)
        run.record_output(TABLES_DIR / "mediation_robustness.csv")
        lc = a.lc_branch()
        if lc is not None:
            lc.to_csv(TABLES_DIR / "mediation_lc_branch.csv", index=False)
            run.record_output(TABLES_DIR / "mediation_lc_branch.csv")
        forest(rows, run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
