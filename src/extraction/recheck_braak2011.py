"""Compare the blind second reading of Braak 2011 with the stored extraction.

    uv run python -m src.extraction.recheck_braak2011

The second reader (a separate agent) read the rendered page images of the
paper and of both supplements, never the stored values or the first
reader's text-layer route, and wrote data/extracted/recheck/braak2011/reread_*.csv.
This script compares every stored number with the re-read one and sets
`machine_rechecked` in the stored files:
  true          the numbers agree;
  discrepancy   they differ, or the second reader did not find the value.
`verified` is never touched: only a human sets it. Every comparison is
listed in data/extracted/recheck/braak2011/recheck_log.md.
"""
from __future__ import annotations

import re
import sys

import pandas as pd

from src.common.paths import REPO_ROOT

X = REPO_ROOT / "data" / "extracted"
RC = X / "recheck" / "braak2011"
SDC1_COLS = {"no_changes": "No changes", "a": "a", "b": "b", "c": "c", "a_c": "a-c", "s1a": "1a", "s1b": "1b",
             "s1a_1b": "1a-1b", "I": "I", "II": "II", "I_II": "I-II", "III": "III", "IV": "IV", "III_IV": "III-IV",
             "V": "V", "VI": "VI", "V_VI": "V-VI", "total": "Total"}
STAGE_OF = {"No changes": "0", "a": "a", "b": "b", "c": "c", "1a": "1a", "1b": "1b", "I": "I", "II": "II",
            "III": "III", "IV": "IV", "V": "V", "VI": "VI"}


def num(x) -> float | None:
    s = str(x).strip().replace(",", "").replace("%", "")
    if s in ("", "nan", "NOT_FOUND", "NOT_PRINTED", "None"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def norm_label(s: str) -> str:
    s = re.sub(r"[–—Y]", "-", str(s).strip())
    s = re.sub(r"\s+", " ", s)
    m = re.fullmatch(r"0*(\d+)\s*-\s*0*(\d+)", s)
    return f"{int(m.group(1))}-{int(m.group(2))}" if m else s.lower()


class Book:
    def __init__(self):
        self.log = []

    def compare(self, df, idx, stored, reread, what):
        a, b = num(stored), num(reread)
        ok = a is not None and b is not None and abs(a - b) < 1e-9
        prev = df.at[idx, "machine_rechecked"]
        df.at[idx, "machine_rechecked"] = "true" if ok and prev != "discrepancy" else "discrepancy"
        self.log.append((what, stored, reread, "agrees" if ok else "DISCREPANCY"))


def main() -> int:
    need = ["reread_table1.csv", "reread_table2.csv", "reread_sdc1.csv", "reread_sdc2.csv", "reread_text_counts.csv"]
    missing = [f for f in need if not (RC / f).exists()]
    if missing:
        raise SystemExit(f"second reading incomplete, missing: {missing}")
    bk = Book()
    files = {n: pd.read_csv(X / n, dtype=str, keep_default_na=False) for n in
             ("braak2011_table1.csv", "braak2011_table2.csv", "braak2011_text_counts.csv",
              "braak2011_stage_counts.csv", "braak2011_sdc1_groups.csv", "braak2011_sdc2_amyloid_by_age.csv")}
    for df in files.values():
        df["machine_rechecked"] = "pending"

    # Table 1
    t1, r1 = files["braak2011_table1.csv"], pd.read_csv(RC / "reread_table1.csv", dtype=str, keep_default_na=False)
    r1["key"] = r1["age_label"].map(norm_label)
    for i, row in t1.iterrows():
        key = "column header" if row["age_label_printed"] == "column header" else norm_label(row["age_label_printed"])
        rr = r1[r1["key"] == key]
        col = f"{row['sex']}_{row['measure']}" if row["sex"] != "total" else "total_n"
        val = rr[col].iloc[0] if len(rr) and col in rr else None
        bk.compare(t1, i, row["value"], val, f"Table 1 | {key} | {col}")

    # Table 2
    t2, r2 = files["braak2011_table2.csv"], pd.read_csv(RC / "reread_table2.csv", dtype=str, keep_default_na=False)
    r2["key"] = r2["at8_stage"].map(lambda s: "column total" if "header" in str(s).lower() else norm_label(s))
    for i, row in t2.iterrows():
        key = row["at8_stage"] if row["at8_stage"] == "column total" else norm_label(row["at8_stage"])
        rr = r2[r2["key"] == key]
        ph = row["amyloid_phase"]
        col = ("total" if ph == "Total" else "phase" + ph.split()[-1]) + "_" + row["measure"]
        val = rr[col].iloc[0] if len(rr) and col in rr else None
        bk.compare(t2, i, row["value"], val, f"Table 2 | {key} | {col}")

    # SDC 1: stage counts (band rows) and grouped cells (incl. total row)
    st, g1 = files["braak2011_stage_counts.csv"], files["braak2011_sdc1_groups.csv"]
    r3 = pd.read_csv(RC / "reread_sdc1.csv", dtype=str, keep_default_na=False)
    r3["key"] = r3["age_label"].map(norm_label)
    inv = {v: k for k, v in SDC1_COLS.items()}
    for i, row in st.iterrows():
        lab = next(k for k, v in STAGE_OF.items() if v == row["stage"])
        key = norm_label(row["age_label_printed"])
        rr = r3[r3["key"] == key]
        col = inv[lab] + "_n"
        val = rr[col].iloc[0] if len(rr) and col in rr else None
        bk.compare(st, i, row["n"], val, f"SDC 1 | {key} | {col}")
    for i, row in g1.iterrows():
        key = norm_label(row["age_label_printed"])
        rr = r3[r3["key"] == key]
        col = inv[row["column"]] + "_" + row["measure"]
        val = rr[col].iloc[0] if len(rr) and col in rr else None
        bk.compare(g1, i, row["value"], val, f"SDC 1 | {key} | {col}")

    # SDC 2
    s2, r4 = files["braak2011_sdc2_amyloid_by_age.csv"], pd.read_csv(RC / "reread_sdc2.csv", dtype=str,
                                                                       keep_default_na=False)
    r4["key"] = r4["age_label"].map(lambda s: "column header" if "header" in str(s).lower() else norm_label(s))
    for i, row in s2.iterrows():
        key = "column header" if row["age_label_printed"] == "column header" else norm_label(row["age_label_printed"])
        rr = r4[r4["key"] == key]
        ph = row["amyloid_phase"]
        col = ("total" if ph == "Total" else "phase" + ph.split()[-1]) + "_" + row["measure"]
        val = rr[col].iloc[0] if len(rr) and col in rr else None
        bk.compare(s2, i, row["value"], val, f"SDC 2 | {key} | {col}")

    # Text and figure counts (task text_i = stored row i)
    tc, r5 = files["braak2011_text_counts.csv"], pd.read_csv(RC / "reread_text_counts.csv", dtype=str,
                                                              keep_default_na=False)
    for i, row in tc.iterrows():
        rr = r5[r5["task"] == f"text_{i}"]
        val = rr["value"].iloc[0] if len(rr) else None
        bk.compare(tc, i, row["value"], val, f"Text | {row['quantity']}")

    for name, df in files.items():
        if (df["machine_rechecked"] == "pending").any():
            raise SystemExit(f"{name}: some rows were not compared")
        df.to_csv(X / name, index=False)
    n_ok = sum(o == "agrees" for *_, o in bk.log)
    lines = ["# Braak 2011: blind second reading vs stored extraction", "",
             "Generated by `python -m src.extraction.recheck_braak2011`. The second reader read rendered page "
             "images and never saw the stored values. `verified` is untouched.", "",
             f"**{n_ok} of {len(bk.log)} values agree.**", "",
             "## Discrepancies", "", "| value | stored | second reading |", "|---|---|---|"]
    lines += [f"| {w} | {a} | {b} |" for w, a, b, o in bk.log if o != "agrees"] or ["| (none) | | |"]
    lines += ["", "## Every comparison", "", "| value | stored | second reading | outcome |", "|---|---|---|---|"]
    lines += [f"| {w} | {a} | {b} | {o} |" for w, a, b, o in bk.log]
    (RC / "recheck_log.md").write_text("\n".join(lines) + "\n")
    print(f"{n_ok}/{len(bk.log)} values agree; log: {RC / 'recheck_log.md'}")
    for w, a, b, o in bk.log:
        if o != "agrees":
            print(f"DISCREPANCY {w}: stored {a!r} vs re-read {b!r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
