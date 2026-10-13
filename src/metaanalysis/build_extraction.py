"""Build data/extracted/meta_extraction.csv (Module 5) from the candidate rows.

    uv run python -m src.metaanalysis.build_extraction

Inputs (committed): data/extracted/candidates/*.csv. These are every
estimate Claude's extraction pass copied from open-access full texts, each
with a verbatim snippet and an exact locator. The source files are not
redistributed; their SHA-256 values are in candidates/source_files.sha256.

This script applies the selection rules below. They are pre-specified in
prereg/metaanalysis.md §3 and repeated here so the code is the record.
It never edits a number. It only chooses rows, assigns labels and records
inclusion decisions. Every estimate, CI bound and verbatim snippet is
copied unchanged. Rows not used by the analysis stay in the file with
include=false and the reason, so the pipeline counts and logs them.

Selection rules (one estimate per study per drug class and comparator type):
  S1  the paper's primary, fully adjusted estimate;
  S2  several analyses with equal stated priority -> the first listed;
  S3  only per-drug estimates against a shared reference -> the drug with the
      narrowest CI on the log scale (most information). Pooling them would
      count the shared reference group repeatedly;
  S4  only exposure-duration or dose categories -> the highest category vs
      none (the positive-control convention);
  S5  several follow-up windows -> the longest;
  S6  outcome: all-cause dementia when reported, otherwise Alzheimer's
      disease.
Exclusion rules (row kept, include=false):
  E1  dementia within another neurodegenerative disease (Parkinson's
      disease): a different disease process. Populations restricted by the
      drug's indication (hypertension, atrial fibrillation) are included;
  E2  the printed CI cannot yield a standard error: zero width, or asymmetric
      on the log scale beyond config extraction.log_symmetry_tolerance;
  E3  exposure contrast confounded by indication in a way the other studies'
      contrasts are not (antidepressant users vs non-users of antidepressants);
  E4  anticholinergic rows with an active comparator are kept as an
      active-comparator row and are not part of the positive control (which
      is use vs non-use).
"""

from __future__ import annotations

import math
import sys

import pandas as pd

from src.common.config import load_config
from src.common.paths import CONFIG_DIR, REPO_ROOT

CAND = REPO_ROOT / "data" / "extracted" / "candidates"
OUT = REPO_ROOT / "data" / "extracted" / "meta_extraction.csv"
OUT_VALID = REPO_ROOT / "data" / "extracted" / "meta_validation_moraros2017.csv"
NOT_ASSESSED = {f"robins_d{i}": "not_assessed" for i in range(1, 8)}

# noradrenergic candidates: (file, row index) -> labels. Indices refer to the
# concatenation [main file rows..., supplementary file rows...].
NORAD = {
    0: dict(drug_class="lipophilic_beta_blocker", contrast="lipophilic_vs_hydrophilic_beta_blocker", pen="high", rule="S1"),
    1: dict(drug_class="lipophilic_beta_blocker", contrast="lipophilic_vs_hydrophilic_beta_blocker", pen="high", rule="S2 (analysis 1 of 4, as-treated)"),
    5: dict(drug_class="alpha1_blocker", contrast="alpha1_blocker_vs_tamsulosin", pen="not_applicable", rule="S5 (10-year window)"),
    7: dict(drug_class="alpha1_blocker", contrast="alpha1_blocker_vs_tamsulosin", pen="not_applicable", rule="S3 (doxazosin: narrowest CI of terazosin/doxazosin/alfuzosin)"),
    9: dict(drug_class="alpha1_blocker", contrast="alpha1_blocker_vs_tamsulosin", pen="not_applicable", rule="S1",
            exclude="E2: CI printed as [0.84-0.84] (t = -867.9) on an IPTW pseudo-population; no usable standard error"),
    10: dict(drug_class="alpha1_blocker", contrast="alpha1_blocker_vs_tamsulosin", pen="not_applicable", rule="S1",
             exclude="E1: men with Parkinson's disease; outcome is dementia within PD"),
    11: dict(drug_class="alpha1_blocker", contrast="none", pen="not_applicable", rule="S1 (class-level, fully adjusted model b)"),
    17: dict(drug_class="alpha1_blocker", contrast="none", pen="not_applicable", rule="S6 (all-cause dementia)"),
    19: dict(drug_class="dhp_calcium_channel_blocker", contrast="brain_penetrant_vs_nonpenetrant_dhp_ccb", pen="high", rule="S1"),
    20: dict(drug_class="lipophilic_beta_blocker", contrast="other_active", pen="high", rule="S1 (active comparator: supplement initiators)"),
    21: dict(drug_class="hydrophilic_beta_blocker", contrast="none", pen="low", rule="S1"),
    22: dict(drug_class="lipophilic_beta_blocker", contrast="none", pen="high", rule="S1"),
    23: dict(drug_class="atomoxetine", contrast="none", pen="not_applicable", rule="S1 (self-controlled; within-person unexposed periods)"),
    24: dict(drug_class="mao_a_inhibitor", contrast="none", pen="not_applicable", rule="S1"),
    27: dict(drug_class="beta_blocker_any", contrast="none", pen="not_applicable", rule="S4 (>= 10 years of use vs non-use)"),
    28: dict(drug_class="beta_blocker_any", contrast="none", pen="not_applicable", rule="S1 (population restricted to AF: indication)"),
    29: dict(drug_class="dhp_calcium_channel_blocker", contrast="none", pen="not_applicable", rule="S1 (vessel-selective CCB, mostly DHP; AF population)"),
    30: dict(drug_class="dual_orexin_antagonist", contrast="z_drug_vs_dual_orexin_antagonist", pen="not_applicable", rule="S1",
             invert=True),
    31: dict(drug_class="beta_blocker_any", contrast="none", pen="not_applicable", rule="S1 (vs untreated hypertension)"),
    32: dict(drug_class="alpha_beta_blocker", contrast="none", pen="not_applicable", rule="S1 (vs untreated hypertension)"),
    33: dict(drug_class="central_alpha2_agonist", contrast="none", pen="not_applicable", rule="S1 (vs untreated hypertension)"),
}
# anticholinergic main estimates (study_id -> labels); rows taken where role = primary_study and is_main_estimate = yes...
ANTICHOL = {
    "Gray2015": {}, "Richardson2018": {}, "Coupland2019": {}, "Grossi2019": {}, "Suh2020": {}, "Joung2019": {},
    "Yang2017": {}, "Poonawalla2023": {}, "Zhu2025": {}, "Iyen2024": {}, "Chi2025": {}, "Gildengers2023": {},
    "Liu2020": {"exclude": "E2: CI asymmetric on the log scale as printed (e.g. 1.043 (0.958-1.212)); SE cannot be backed out"},
    "PhamNguyen2025": {"exclude": "E1: Parkinson's disease population; outcome is Parkinson disease dementia"},
    "Amrein2026": {"exclude": "E3: antidepressant (TCA) users vs non-users of antidepressants: confounded by indication (depression)"},
    "Li2025": {"comparator_type": "active", "contrast": "other_active",
               "note": "E4: active comparator (beta-3 agonist): not part of the use-vs-non-use positive control"},
}


def measure(label: str) -> str:
    s = label.strip().upper()
    for m in ("IRR", "HR", "OR", "RR"):
        if s.startswith(m):
            return m
    if s.startswith("RISK RATIO") or "RISK RATIO" in s:
        return "RR"
    raise ValueError(f"cannot map effect measure {label!r}")


def ndec(v: str) -> int:
    return len(v.split(".")[1]) if "." in v else 0


def asym_excess(est: str, lo: str, hi: str, tol: float) -> float:
    """Offset of the estimate from the log-CI midpoint minus what the rule
    allows (tol x half-width + the rounding shift for the printed decimals);
    > 0 means the CI fails the symmetry check (same rule as meta.R)."""
    e, l, h = float(est), float(lo), float(hi)
    a, b, c = math.log(e), math.log(l), math.log(h)
    shift = 0.5 * 10 ** -ndec(est) / e + (0.5 * 10 ** -ndec(lo) / l + 0.5 * 10 ** -ndec(hi) / h) / 2
    return abs(a - (b + c) / 2) - (tol * (c - b) / 2 + shift)


def main() -> int:
    cfg = load_config(CONFIG_DIR / "metaanalysis.yaml")
    tol = cfg.get("extraction.log_symmetry_tolerance")
    rd = lambda f: pd.read_csv(CAND / f, dtype=str, keep_default_na=False)  # noqa: E731  numbers kept as printed
    nor = pd.concat([rd("meta_noradrenergic_candidates.csv"), rd("meta_noradrenergic_candidates_supplementary.csv")],
                    ignore_index=True)
    ach = rd("meta_anticholinergic_candidates.csv")
    rows = []
    for i, lab in NORAD.items():
        r = nor.loc[i]
        rows.append(_row(r.study_id, r.pmcid, lab["drug_class"], r.drug, r.comparator_type, r.comparator, lab["contrast"],
                         r.outcome, r.effect_measure, r.estimate, r.ci_lower, r.ci_upper, r.ci_level, r.adjusted_for,
                         lab["pen"], r.locator, r.verbatim, lab["rule"], lab.get("exclude"), lab.get("invert", False), tol))
    main_ach = ach[(ach.role == "primary_study") & ach.is_main_estimate.astype(str).str.startswith("yes")]
    for sid, lab in ANTICHOL.items():
        sel = main_ach[main_ach.study_id == sid]
        if len(sel) != 1:
            raise ValueError(f"{sid}: expected exactly one main anticholinergic row, found {len(sel)}")
        r = sel.iloc[0]
        rows.append(_row(r.study_id, r.pmcid, "anticholinergic", r.exposure, lab.get("comparator_type", "non_user"),
                         r.comparator, lab.get("contrast", "none"), r.outcome, r.effect_measure, r.estimate, r.ci_lower,
                         r.ci_upper, r.ci_level, r.adjusted_for, "not_applicable", r.locator, r.verbatim,
                         "S1/S4 (highest exposure vs none)" + (f"; {lab['note']}" if "note" in lab else ""),
                         lab.get("exclude"), False, tol))
    df = pd.DataFrame(rows)
    unknown = set(df.drug_class) - set(cfg.get("drug_classes"))
    if unknown:
        raise ValueError(f"drug classes not in config/metaanalysis.yaml: {sorted(unknown)}")
    df.to_csv(OUT, index=False)
    print(f"{OUT.relative_to(REPO_ROOT)}: {len(df)} rows, {int((df.include == 'true').sum())} included")

    mor = ach[ach.parent_meta_analysis.astype(str).eq("Moraros2017") | (ach.study_id.eq("Moraros2017") & ach.role.eq("ma_pooled"))]
    pooled = mor[(mor.role == "ma_pooled")].iloc[[0]]
    comps = mor[mor.role == "ma_component"]
    v = pd.concat([pooled, comps])
    out = pd.DataFrame({
        "role": v.role.values, "study_id": v.study_id.values, "estimate": v.estimate.values,
        "ci_lower": v.ci_lower.values, "ci_upper": v.ci_upper.values,
        "ci_level": [float(str(c).rstrip("%")) / 100 if str(c).endswith("%") else float(c) for c in v.ci_level],
        "effect_measure": v.effect_measure.values, "source_id": "moraros2017", "locator": v.locator.values,
        "quantity": [f"{m}; {e}" for m, e in zip(v.effect_measure, v.exposure)], "verbatim": v.verbatim.values,
        "verified": "false", "machine_rechecked": "false"})
    out.to_csv(OUT_VALID, index=False)
    print(f"{OUT_VALID.relative_to(REPO_ROOT)}: pooled + {len(comps)} components")
    return 0


def _row(study_id, pmcid, drug_class, drug, comparator_type, comparator, contrast, outcome, eff, est, lo, hi, level,
         covariates, pen, locator, verbatim, rule, exclude, invert, tol):
    est, lo, hi = str(est).strip(), str(lo).strip(), str(hi).strip()      # numbers exactly as printed
    level = float(str(level).rstrip("%")) / 100 if str(level).endswith("%") else float(level)   # "95%" -> 0.95 (notation only)
    if exclude is None and float(hi) <= float(lo):
        exclude = "E2: zero-width CI"
    if exclude is None and asym_excess(est, lo, hi, tol) > 0:
        exclude = "E2: CI asymmetric on the log scale beyond the rounding-aware tolerance"
    return {"study_id": study_id, "drug_class": drug_class, "drug": str(drug), "comparator_type": comparator_type,
            "comparator": str(comparator), "contrast": contrast, "outcome": str(outcome), "effect_measure": measure(eff),
            "estimate": est, "ci_lower": lo, "ci_upper": hi, "ci_level": repr(level),
            "covariates": str(covariates), "brain_penetrance": pen, **NOT_ASSESSED,
            "source_id": str(pmcid), "locator": str(locator),
            "quantity": f"{eff} (as printed); selection rule {rule}", "verbatim": str(verbatim),
            "verified": "false", "machine_rechecked": "false",
            "invert_ratio": "true" if invert else "false",
            "include": "false" if exclude else "true", "exclusion_reason": exclude or "none"}


if __name__ == "__main__":
    sys.exit(main())
