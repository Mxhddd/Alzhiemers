"""Load the RADC ROS/MAP Standard dataset into one row per participant.

Variable names come from RADC's PUBLIC documentation (Research Resource
Sharing Hub, https://www.radc.rush.edu/docs/var/standardDatasetVariables.htm
and the per-variable pages /docs/var/detail.htm?variable=<name>, read
2026-10-21). The LAYOUT of the delivered files (one cross-sectional file plus
one longitudinal file keyed by projid and fu_year, CSV) is an assumption to be
confirmed against the codebook on receipt; every name is in
config/mediation.yaml so it can be corrected there, and any absent column
stops the run.

Nothing here decides anything after seeing the data: exposure timing, codes
and the LC branch are fixed in the config before the data are opened.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.common.io import analysis_input
from src.common.paths import REPO_ROOT


class CodebookError(ValueError):
    pass


def _read(path_rel: str, provided_by: str) -> pd.DataFrame:
    path = analysis_input(REPO_ROOT / path_rel, provided_by=provided_by)
    return pd.read_csv(path)


def _need(df: pd.DataFrame, cols: list[str], where: str) -> None:
    absent = [c for c in cols if c not in df.columns]
    if absent:
        raise CodebookError(f"{where}: columns {absent} not found. Names come from RADC's public "
                            "documentation and must be confirmed against the delivered codebook; "
                            "correct them in config/mediation.yaml (before any analysis), never in code.")


def apoe4_carrier(genotype: pd.Series, carrier_codes: list[int]) -> pd.Series:
    g = pd.to_numeric(genotype, errors="coerce")
    out = g.isin(carrier_codes).astype(float)
    return out.where(g.notna())


def exposure_summaries(long: pd.DataFrame, *, id_col: str, visit_col: str, var: str, age_col: str | None,
                       age_bl: pd.Series, dx_col: str) -> pd.DataFrame:
    """Per participant: earliest non-missing value (primary), its visit, age and
    diagnosis; the mean of all non-missing values; the last non-missing value."""
    d = long[[id_col, visit_col, var, dx_col] + ([age_col] if age_col else [])].dropna(subset=[var])
    d = d.sort_values([id_col, visit_col])
    g = d.groupby(id_col, sort=False)
    first = g.head(1).set_index(id_col)
    last = g.tail(1).set_index(id_col)
    out = pd.DataFrame({
        f"{var}__earliest": first[var],
        f"{var}__visit": first[visit_col],
        f"{var}__dx": first[dx_col],
        f"{var}__mean": g[var].mean(),
        f"{var}__last": last[var],
    })
    if age_col:
        out[f"{var}__age"] = first[age_col]
    else:   # age at the visit = baseline age + years since baseline (fu_year)
        out[f"{var}__age"] = age_bl.reindex(out.index).to_numpy() + out[f"{var}__visit"].to_numpy()
    return out


def load_radc(cfg) -> tuple[pd.DataFrame, dict]:
    """Return (participants, info). One row per participant in the configured
    cohorts, all participants (alive or dead, autopsied or not), so that the
    inverse-probability-of-autopsy weights can be estimated."""
    v = cfg.get("variables")
    ids = cfg.get("data.columns")
    basic = _read(cfg.get("data.basic_file"), "researcher: RADC Standard dataset, cross-sectional file (requested, not received)")
    long = _read(cfg.get("data.long_file"), "researcher: RADC Standard dataset, longitudinal file (requested, not received)")
    expos = [e["name"] for e in cfg.get("exposures")]
    lc_vars = [x["name"] for x in cfg.get("lc.variables")]
    cross = [ids["id"], ids["study"], v["outcome"], v["outcome_sensitivity"], v["amyloid"], v["amyloid_sensitivity"],
             v["age_death"], v["age_baseline"], v["sex"], v["education"], v["apoe"], v["died"]] + lc_vars
    _need(basic, cross, cfg.get("data.basic_file"))
    age_col = ids["visit_age"] if ids["visit_age"] != "age_bl_plus_fu_year" else None
    lcols = [ids["id"], ids["visit"], v["diagnosis_visit"], v["cognition_visit"]] + expos + ([age_col] if age_col else [])
    _need(long, lcols, cfg.get("data.long_file"))
    if basic[ids["id"]].duplicated().any():
        raise CodebookError("cross-sectional file has duplicate participant ids")
    p = basic[cross].copy().set_index(ids["id"])
    cohorts = cfg.get("data.cohorts")
    p = p[p[ids["study"]].astype(str).str.strip().isin(cohorts)]
    p["apoe4"] = apoe4_carrier(p[v["apoe"]], cfg.get("variables.apoe4_carrier_codes"))
    p["map"] = (p[ids["study"]].astype(str).str.strip() == "MAP").astype(float)
    long = long[long[ids["id"]].isin(p.index)]
    for e in expos:
        s = exposure_summaries(long, id_col=ids["id"], visit_col=ids["visit"], var=e, age_col=age_col,
                               age_bl=p[v["age_baseline"]], dx_col=v["diagnosis_visit"])
        p = p.join(s, how="left")
    bl = long.sort_values([ids["id"], ids["visit"]]).dropna(subset=[v["cognition_visit"]]).groupby(ids["id"]).head(1)
    p["cogn_baseline"] = bl.set_index(ids["id"])[v["cognition_visit"]].reindex(p.index)
    p["has_pathology"] = p[v["outcome"]].notna()
    info = {"n_participants": int(len(p)), "n_died": int((p[v["died"]] == 1).sum()),
            "n_with_braak": int(p["has_pathology"].sum()),
            "n_by_cohort": p[ids["study"]].astype(str).value_counts().to_dict()}
    return p, info
