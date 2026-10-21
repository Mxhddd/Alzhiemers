"""ADNI tau-PET -> group-level regional tau maps for Module 6
(prereg/network.md, section 3). Real-data step; written before any ADNI data
were seen.

    python -m src.network.adni_tau        (refuses until prereg/network.md is frozen)

The ADNI Data Use Agreement forbids putting participant-level ADNI data into
third-party AI tools, so this step is run by the researcher on his own
computer, never inside an AI session. Its outputs are group-level summaries
(one value per cortical region, and counts), which the agreement allows to be
reported.

Steps (all settings in config/network.yaml, `adni`):
1. Tau scans: tracer = flortaucipir, visual QC pass, first such scan per
   participant.
2. Amyloid status: the amyloid PET scan nearest in time to the tau scan,
   within a window; the Berkeley lab's own positivity flag.
3. Diagnosis: the diagnostic summary nearest in time to the tau scan, within
   a window (CN or MCI).
4. Age at the tau scan from year of birth (mid-year assumed); sex.
5. Reference group A- CN; target groups A+ CN (primary) and A+ MCI.
6. Per region: linear regression of SUVR on age and sex in the reference
   group; each target participant's W-score = (SUVR - predicted) / residual
   SD; the group map is the mean W-score per region.
Every exclusion is counted and written out; nothing is imputed.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

DK_REGIONS = (
    "bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", "cuneus", "entorhinal", "fusiform",
    "inferiorparietal", "inferiortemporal", "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal",
    "lingual", "medialorbitofrontal", "middletemporal", "parahippocampal", "paracentral", "parsopercularis",
    "parsorbitalis", "parstriangularis", "pericalcarine", "postcentral", "posteriorcingulate", "precentral",
    "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal", "superiorfrontal", "superiorparietal",
    "superiortemporal", "supramarginal", "frontalpole", "temporalpole", "transversetemporal", "insula",
)
assert len(DK_REGIONS) == 34


def cortical_regions() -> list[str]:
    """Node names of the 68 Desikan-Killiany cortical regions (FreeSurfer LUT style)."""
    return [f"ctx-{h}-{r}" for h in ("lh", "rh") for r in DK_REGIONS]


def suvr_columns(pattern: str) -> dict[str, str]:
    """Map node name -> ADNI column name, e.g. ctx-lh-entorhinal -> CTX_LH_ENTORHINAL_SUVR."""
    return {f"ctx-{h}-{r}": pattern.format(hemi=h.upper(), region=r.upper())
            for h in ("lh", "rh") for r in DK_REGIONS}


def first_scans(tau: pd.DataFrame, *, id_col: str, date_col: str, tracer_col: str, tracer: str,
                qc_col: str, qc_pass: list, counts: dict) -> pd.DataFrame:
    counts["tau_rows"] = len(tau)
    t = tau[tau[tracer_col] == tracer]
    counts["tau_rows_other_tracer"] = len(tau) - len(t)
    ok = t[qc_col].isin(qc_pass)
    counts["tau_rows_failed_qc"] = int((~ok).sum())
    t = t[ok].copy()
    t[date_col] = pd.to_datetime(t[date_col])
    t = t.sort_values([id_col, date_col]).drop_duplicates(id_col, keep="first")
    counts["participants_with_tau_scan"] = len(t)
    return t


def attach_nearest(base: pd.DataFrame, other: pd.DataFrame, *, id_col: str, base_date: str, other_date: str,
                   value_col: str, out_col: str, max_days: int) -> pd.DataFrame:
    """For each base row, the value from the `other` row of the same participant
    nearest in time, if within max_days (ties: the earlier one). Missing
    values in `other` are dropped first. Rows without a match get NaN."""
    o = other[[id_col, other_date, value_col]].dropna().copy()
    o[other_date] = pd.to_datetime(o[other_date])
    m = base[[id_col, base_date]].reset_index().merge(o, on=id_col, how="left")
    m["_gap"] = (m[other_date] - m[base_date]).dt.days.abs()
    m = m[m["_gap"] <= max_days].sort_values(["index", "_gap", other_date])
    best = m.drop_duplicates("index", keep="first").set_index("index")[value_col]
    out = base.copy()
    out[out_col] = best.reindex(base.index)
    return out


def w_scores(ref: pd.DataFrame, tgt: pd.DataFrame, regions: list[str], covariates: list[str]) -> pd.DataFrame:
    """Per-region W-scores of the target rows against a regression on
    `covariates` fitted in the reference rows."""
    def design(df):
        return np.column_stack([np.ones(len(df))] + [df[c].to_numpy(float) for c in covariates])
    Xr, Xt = design(ref), design(tgt)
    if len(ref) <= Xr.shape[1]:
        raise ValueError("reference group too small for the W-score regression")
    Y = ref[regions].to_numpy(float)
    beta, *_ = np.linalg.lstsq(Xr, Y, rcond=None)
    resid = Y - Xr @ beta
    sd = np.sqrt((resid ** 2).sum(axis=0) / (len(ref) - Xr.shape[1]))
    if np.any(sd <= 0):
        raise ValueError("zero residual SD in the reference group")
    return pd.DataFrame((tgt[regions].to_numpy(float) - Xt @ beta) / sd, columns=regions, index=tgt.index)


def build_group_maps(tau: pd.DataFrame, amy: pd.DataFrame, dx: pd.DataFrame, demog: pd.DataFrame,
                     a: dict) -> tuple[dict[str, pd.Series], dict]:
    """a = the `adni` section of config/network.yaml (as a dict)."""
    c = a["columns"]
    counts: dict = {}
    cols = suvr_columns(a["suvr_column"])
    missing = [v for v in cols.values() if v not in tau.columns]
    if missing:
        raise KeyError(f"tau file lacks regional columns, e.g. {missing[:3]} (no imputation)")
    t = first_scans(tau, id_col=c["id"], date_col=c["tau_date"], tracer_col=c["tracer"], tracer=a["tau_tracer"],
                    qc_col=c["qc"], qc_pass=list(a["qc_pass"]), counts=counts)
    t = t.rename(columns={v: k for k, v in cols.items()})
    regions = list(cols)
    t = attach_nearest(t, amy, id_col=c["id"], base_date=c["tau_date"], other_date=c["amyloid_date"],
                       value_col=c["amyloid_status"], out_col="amyloid", max_days=int(a["amyloid_window_days"]))
    t = attach_nearest(t, dx, id_col=c["id"], base_date=c["tau_date"], other_date=c["dx_date"],
                       value_col=c["diagnosis"], out_col="dx", max_days=int(a["diagnosis_window_days"]))
    d = demog[[c["id"], c["birth_year"], c["sex"]]].dropna().drop_duplicates(c["id"], keep="first")
    t = t.merge(d, on=c["id"], how="left")
    date = t[c["tau_date"]]
    t["age"] = date.dt.year + (date.dt.dayofyear - 0.5) / 365.25 - (t[c["birth_year"]] + 0.5)
    t["sex"] = (t[c["sex"]] == a["female_code"]).astype(float).where(t[c["sex"]].notna())
    for col, label in (("amyloid", "no_amyloid_scan_in_window"), ("dx", "no_diagnosis_in_window"),
                       ("age", "no_birth_year"), ("sex", "no_sex")):
        n = int(t[col].isna().sum())
        counts[f"excluded_{label}"] = n
        t = t[t[col].notna()]
    n_reg = int(t[regions].isna().any(axis=1).sum())
    counts["excluded_missing_regional_suvr"] = n_reg
    t = t[~t[regions].isna().any(axis=1)]
    pos = t["amyloid"] == a["amyloid_positive_value"]
    cu = t["dx"] == a["diagnosis_codes"]["cu"]
    mci = t["dx"] == a["diagnosis_codes"]["mci"]
    groups = {"reference_aneg_cu": t[~pos & cu], "apos_cu": t[pos & cu], "apos_mci": t[pos & mci]}
    for g, df in groups.items():
        counts[f"n_{g}"] = len(df)
    counts["excluded_other_diagnosis"] = int((~cu & ~mci).sum())
    n_min = int(a["min_group_n"])
    if counts["n_reference_aneg_cu"] < n_min:
        raise ValueError(f"reference group has {counts['n_reference_aneg_cu']} < {n_min} participants")
    maps = {}
    for g in ("apos_cu", "apos_mci"):
        if counts[f"n_{g}"] < n_min:
            counts[f"{g}_status"] = f"not built: fewer than {n_min} participants"
            continue
        W = w_scores(groups["reference_aneg_cu"], groups[g], regions, list(a["w_score_covariates"]))
        maps[g] = W.mean(axis=0).rename("value")
        counts[f"{g}_status"] = "built"
    return maps, counts


def main() -> int:  # pragma: no cover - real-data path
    from src.common.config import load_config
    from src.common.io import analysis_input
    from src.common.paths import CONFIG_DIR, REPO_ROOT, TABLES_DIR
    from src.common.prereg import require_preregistration
    from src.common.runlog import start_run

    prereg = require_preregistration("network")
    cfg_path = CONFIG_DIR / "network.yaml"
    cfg = load_config(cfg_path)
    a = cfg.as_dict()["adni"]
    with start_run(module="network_adni_tau", config_path=cfg_path, seed=None,
                   packages=["numpy", "pandas"], preregistration=prereg) as run:
        def read(key):
            return pd.read_csv(analysis_input(REPO_ROOT / a["files"][key], provided_by="researcher (ADNI download)"))
        amy, dx, demog = read("amyloid"), read("diagnosis"), read("demographics")
        out_dir = REPO_ROOT / a["output_dir"]
        out_dir.mkdir(parents=True, exist_ok=True)
        rows = []
        for version in ("pvc", "nopvc"):
            maps, counts = build_group_maps(read(f"tau_{version}"), amy, dx, demog, a)
            for g, m in maps.items():
                p = out_dir / f"tau_map_{g}_{version}.csv"
                m.rename_axis("region").reset_index().to_csv(p, index=False)
                run.record_output(p)
            rows += [{"tau_file": version, "item": k, "value": v} for k, v in counts.items()]
            run.info(f"{version}: {counts}")
        TABLES_DIR.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(TABLES_DIR / "network_adni_counts.csv", index=False)
        run.record_output(TABLES_DIR / "network_adni_counts.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
