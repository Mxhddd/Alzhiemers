"""Extract every printed count from Braak et al. 2011 and its two supplements.

    uv run python -m src.extraction.braak2011

Source: Braak H, Thal DR, Ghebremedhin E, Del Tredici K. Stages of the
pathologic process in Alzheimer disease: age categories from 1 to 100 years.
J Neuropathol Exp Neurol 2011;70(11):960-969. doi:10.1097/NEN.0b013e318232a379
Inputs (git-ignored, data/raw/papers/):
  braak2011.pdf       full text, supplied by the researcher 2026-10-13
  braak2011_sdc1.doc  Supplemental Digital Content 1 (Supplementary Table 1),
                      https://links.lww.com/NEN/A273, the publisher's open permalink
  braak2011_sdc2.doc  Supplemental Digital Content 2 (Supplementary Table 2),
                      https://links.lww.com/NEN/A274
Routes: `pdftotext -layout` for the PDF text layer; LibreOffice (headless)
converts the .doc files to HTML, read with a stdlib table parser. Numbers
are taken only where they are printed in text or tables (and the totals
printed above the Figure 2 bars and in the Figure 1 panel titles); nothing
is read off bar heights.

Writes (all rows: source_id braak2011, verified false, machine_rechecked
false until the blind second pass, src/common/apply_recheck.py):
  data/extracted/braak2011_table1.csv                Table 1 cells (sex x age)
  data/extracted/braak2011_table2.csv                Table 2 cells (amyloid phase x AT8 stage)
  data/extracted/braak2011_text_counts.csv           counts printed in Results, Discussion, figures
  data/extracted/braak2011_stage_counts.csv          SDC 1, stage x age counts (Module 2 input)
  data/extracted/braak2011_sdc1_groups.csv           SDC 1 grouped columns and printed percentages
  data/extracted/braak2011_sdc2_amyloid_by_age.csv   SDC 2, amyloid phase x age (Module 2 A(t) input)
  data/extracted/braak2011_sum_checks.md             every internal consistency check, pass or fail
Exit code 1 if any check fails (the files are still written, so the failure
can be inspected).
"""
from __future__ import annotations

import re
import subprocess
import sys
import tempfile
from html.parser import HTMLParser
from pathlib import Path

import pandas as pd

from src.common.paths import REPO_ROOT

RAW = REPO_ROOT / "data" / "raw" / "papers"
OUT = REPO_ROOT / "data" / "extracted"
SRC = "braak2011"
DASH = "[Y–—-]"          # the PDF text layer encodes the en dash as 'Y'
BANDS = [(1, 10), (11, 20), (21, 30), (31, 40), (41, 50), (51, 60), (61, 70), (71, 80), (81, 90), (91, 100)]
# Supplementary Table 1 column order (AT8 stage); "No changes" is the model's stage "0".
SDC1_STAGES = {"No changes": "0", "a": "a", "b": "b", "c": "c", "1a": "1a", "1b": "1b",
               "I": "I", "II": "II", "III": "III", "IV": "IV", "V": "V", "VI": "VI"}
SDC1_GROUPS = ["No changes", "a-c", "1a-1b", "I-II", "III-IV", "V-VI"]


def meta(locator: str, quantity: str) -> dict:
    return {"source_id": SRC, "locator": locator, "quantity": quantity, "verified": "false",
            "machine_rechecked": "false"}


def pdf_text() -> str:
    pdf = RAW / "braak2011.pdf"
    return subprocess.run(["pdftotext", "-layout", str(pdf), "-"], check=True, capture_output=True,
                          text=True).stdout


class _Tables(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tables, self.row, self.cell = [], None, None

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            self.tables.append([])
        elif tag == "tr":
            self.row = []
        elif tag in ("td", "th"):
            self.cell = ""
        elif tag == "br" and self.cell is not None:
            self.cell += " "

    def handle_endtag(self, tag):
        if tag in ("td", "th") and self.cell is not None:
            self.row.append(" ".join(self.cell.split()))
            self.cell = None
        elif tag == "tr" and self.row is not None:
            if any(self.row):
                self.tables[-1].append(self.row)
            self.row = None

    def handle_data(self, data):
        if self.cell is not None:
            self.cell += data


def doc_tables(name: str) -> list[list[list[str]]]:
    with tempfile.TemporaryDirectory() as d:
        subprocess.run(["soffice", "--headless", "--norestore", "--convert-to", "html", "--outdir", d,
                        str(RAW / name)], check=True, capture_output=True, timeout=180)
        html = next(Path(d).glob("*.html")).read_text(encoding="utf-8", errors="replace")
    p = _Tables()
    p.feed(html)
    return p.tables


def band_of(label: str) -> tuple[int, int]:
    a, b = (int(x) for x in re.split(DASH, label))
    return (1 if a == 0 else a, b)      # SDC 1 prints the first band as "00-10"; the cohort is aged 1-100


def n_pct(cell: str) -> tuple[int, str | None]:
    m = re.fullmatch(r"(\d[\d,]*)\s*(?:\(([\d.]+)\))?", cell.strip())
    if not m:
        raise ValueError(f"cannot read 'n (pct)' from {cell!r}")
    return int(m.group(1).replace(",", "")), m.group(2)


# ------------------------------------------------------------------ PDF ----
def table1(txt: str) -> pd.DataFrame:
    rows = []
    for line in txt.splitlines():
        m = re.match(rf"^(\d+){DASH}(\d+)\s+(\d+)\s+([\d.]+)\s+(\d+)\s+([\d.]+)\s+(\d+)\s*$", line[:74])
        if not m:
            continue
        lo, hi = int(m.group(1)), int(m.group(2))
        label = f"{lo}-{hi}"
        for sex, n, pct in (("female", m.group(3), m.group(4)), ("male", m.group(5), m.group(6)),
                            ("total", m.group(7), None)):
            rows.append({"age_label_printed": label, "age_lower": lo, "age_upper": hi, "sex": sex,
                         "measure": "n", "value": n,
                         **meta("Table 1", f"Number of brains, {sex}, age category {label} y")})
            if pct is not None:
                rows.append({"age_label_printed": label, "age_lower": lo, "age_upper": hi, "sex": sex,
                             "measure": "pct", "value": pct,
                             **meta("Table 1", f"Percent of the age category's brains that are {sex}, "
                                               f"age category {label} y")})
    head = re.search(r"\(n = ([\d,]+)\)\s+\(n =\s*([\d,]+)\)\s+Total", txt)
    tot = re.search(r"Category, y.*\(n = ([\d,]+)\)", txt)
    for sex, v in (("female", head.group(1)), ("male", head.group(2)), ("total", tot.group(1))):
        rows.append({"age_label_printed": "column header", "age_lower": 1, "age_upper": 100, "sex": sex,
                     "measure": "n", "value": v.replace(",", ""),
                     **meta("Table 1, column header", f"Number of brains, {sex}, all ages (column header n)")})
    return pd.DataFrame(rows)


T2_ROWS = ["0", "a-c", "1a-1b", "I", "II", "III", "IV", "V", "VI"]
T2_COLS = ["Phase 0", "Phase 1", "Phase 2", "Phase 3", "Phase 4", "Total"]


def table2(txt: str) -> pd.DataFrame:
    rows = []
    start = txt.index("TABLE 2.")
    block = txt[start:start + 3000].splitlines()
    head = re.findall(r"\(n = ([\d,]+), ([\d.]+)%\)", " ".join(block[:4]))
    if len(head) != 6:
        raise ValueError(f"Table 2 header: expected 6 '(n = .., ..%)' groups, found {len(head)}")
    for col, (n, pct) in zip(T2_COLS, head):
        for meas, v in (("n", n.replace(",", "")), ("pct", pct)):
            rows.append({"at8_stage": "column total", "amyloid_phase": col, "measure": meas, "value": v,
                         **meta("Table 2, column header",
                                f"{'Number' if meas == 'n' else 'Percent of all 2,332'} of brains in amyloid {col} "
                                "(all AT8 stages)")})
    pat = re.compile(rf"^(0|a{DASH}c|1a{DASH}1b|VI|V|IV|III|II|I)\s+" + r"\s+".join([r"(\d+) \(([\d.]+)\)"] * 6))
    found = []
    for line in block:
        m = pat.match(line.strip())
        if not m:
            continue
        stage = re.sub(DASH, "-", m.group(1))
        found.append(stage)
        for k, col in enumerate(T2_COLS):
            n, pct = m.group(2 + 2 * k), m.group(3 + 2 * k)
            for meas, v in (("n", n), ("pct", pct)):
                what = ("Number of brains" if meas == "n" else "Row percent of brains")
                rows.append({"at8_stage": stage, "amyloid_phase": col, "measure": meas, "value": v,
                             **meta("Table 2", f"{what} with AT8 stage {stage} and amyloid {col}")})
    if found != T2_ROWS:
        raise ValueError(f"Table 2 rows: expected {T2_ROWS}, found {found}")
    return pd.DataFrame(rows)


# (quantity, regex with one group, locator, unit). Every regex must match exactly once.
TEXT = [
    ("Cases devoid of AD-related intraneuronal tau (no changes)", r"was found in (\d+) of 2,332 cases", "Results, Absence of Intraneuronal Lesions", "cases"),
    ("No-change cases younger than 10 years", r"categories; (\d+) persons were\s+younger than 10 years", "Results, Absence of Intraneuronal Lesions", "cases"),
    ("No-change cases younger than 20 years", r"younger than 10 years; (\d+) were younger than 20 years", "Results, Absence of Intraneuronal Lesions", "cases"),
    ("Age of the younger of the two eldest no-change cases", r"the eldest were (\d+) and \d+ years", "Results, Absence of Intraneuronal Lesions", "years"),
    ("Age of the eldest no-change case", r"the eldest were \d+ and (\d+) years", "Results, Absence of Intraneuronal Lesions", "years"),
    ("Pretangle stage a cases", r"(Forty) of\s+2,332 cases were classi", "Results, Pretangle Stage a", "cases"),
    ("Pretangle stage b cases", r"(Twelve) cases \(designated Stage b\)", "Results, Pretangle Stage b", "cases"),
    ("Pretangle stage c cases", r"(Six) of 2,332 cases displayed more widely", "Results, Pretangle Stage c", "cases"),
    ("Pretangle stages a-c cases (no cortical AT8-ir projection neurons)", r"In summary, (\d+) of 2,332 cases", "Results, Pretangle Stage c", "cases"),
    ("Pretangle stage 1a cases", r"ndings, (\d+) of\s+2,332 cases displayed mild cortical", "Results, Pretangle Stages 1a and 1b", "cases"),
    ("Pretangle stage 1b cases", r"previously, (\d+) of 2,332 cases showed", "Results, Pretangle Stages 1a and 1b", "cases"),
    ("Pretangle stages 1a and 1b cases", r"All 1a\s+and 1b cases \(n = (\d+)\)", "Results, Pretangle Stages 1a and 1b", "cases"),
    ("NFT stage I cases", r"The (\d+) NFT Stage I cases", "Results, Neurofibrillary Stages I to VI", "cases"),
    ("NFT stage II cases", r"The (\d+) NFT Stage\s+II cases", "Results, Neurofibrillary Stages I to VI", "cases"),
    ("NFT stage III cases", r"There were (\d+) NFT Stage III cases", "Results, Neurofibrillary Stages I to VI", "cases"),
    ("NFT stage IV cases", r"in (\d+) NFT Stage IV cases", "Results, Neurofibrillary Stages I to VI", "cases"),
    ("NFT stage V cases", r"The (\d+) NFT Stage V\s", "Results, Neurofibrillary Stages I to VI", "cases"),
    ("NFT stage VI cases", r"in the (\d+) NFT Stage VI cases", "Results, Neurofibrillary Stages I to VI", "cases"),
    ("NFT stages V and VI cases", r"Stages V and VI \((\d+) cases\)", "Results, Neurofibrillary Stages I to VI", "cases"),
    ("Cases with beta-amyloid plaques", r"All ([\d,]+) cases that displayed the presence of", "Results, Extracellular beta-Amyloid Protein Aggregation", "cases"),
    ("Percent males among ages 1 to 80 years", r"ing from 1 to 80 years \((\d+)% males\)", "Results, Statistical Analyses", "percent"),
    ("Percent females among ages older than 80 years", r"older than 80 years \(([\d.]+)% females\)", "Results, Statistical Analyses", "percent"),
    ("Younger-age cases with abnormal tau only in LC projection neurons (stage 'a'), as stated in the Discussion",
     r"We observed (\d+) cases in younger age groups", "Discussion, paragraph 3", "cases"),
    ("Cases with abnormal tau in pretangle or early NFT stages without beta-amyloid plaques (numerator)",
     r"plaques \(([\d,]+)/2,332", "Discussion, conclusion paragraph", "cases"),
]
WORDS = {"Forty": "40", "Twelve": "12", "Six": "6"}
# Printed in the Figure 1 panel titles and above the Figure 2 bars (the figures are images; read visually).
FIGURE = ([("Figure 1 title: all cases", "2332"), ("Figure 1 panel: cases devoid of changes", "10"),
           ("Figure 1 panel: stages a-c", "58"), ("Figure 1 panel: stages 1a-1b", "274"),
           ("Figure 1 panel: stages I-II", "1086"), ("Figure 1 panel: stages III-IV", "682"),
           ("Figure 1 panel: stages V-VI", "222")]
          + [(f"Figure 2A, n printed above the {lo}-{hi} y bar", n)
             for (lo, hi), n in zip(BANDS, ["8", "24", "61", "100", "188", "330", "492", "570", "484", "75"])]
          + [(f"Figure 2B, n printed above the {lo}-{hi} y bar", n)
             for (lo, hi), n in zip(BANDS, ["8", "24", "61", "100", "188", "330", "492", "570", "484", "75"])])


def text_counts(txt: str) -> pd.DataFrame:
    flat = re.sub(r"\s+", " ", txt)
    # the two-column layout interleaves lines; search each column-joined paragraph via the flattened
    # text first, then fall back to the raw text (which keeps each column's lines contiguous)
    raw = subprocess.run(["pdftotext", "-raw", str(RAW / "braak2011.pdf"), "-"], check=True,
                         capture_output=True, text=True).stdout
    raw_flat = re.sub(r"\s+", " ", raw.replace("-\n", "-"))
    rows = []
    for q, rx, loc, unit in TEXT:
        hits = re.findall(rx.replace(r"\s+", " ").replace(r"\s", " "), raw_flat)
        if len(hits) != 1:
            hits = re.findall(rx.replace(r"\s+", " ").replace(r"\s", " "), flat)
        if len(hits) != 1:
            raise ValueError(f"text count '{q}': pattern matched {len(hits)} times")
        v = WORDS.get(hits[0], hits[0]).replace(",", "")
        m = re.search(rx.replace(r"\s+", " ").replace(r"\s", " "), raw_flat) or re.search(
            rx.replace(r"\s+", " ").replace(r"\s", " "), flat)
        verb = (raw_flat if m.string is raw_flat else flat)[max(0, m.start() - 30):m.end() + 30].strip()
        rows.append({"quantity": q, "value": v, "unit": unit, "verbatim": verb, **meta(loc, q)})
    for q, v in FIGURE:
        loc = "Figure 1 (panel titles)" if q.startswith("Figure 1") else "Figure 2 (numbers above bars)"
        rows.append({"quantity": q, "value": v, "unit": "cases", "verbatim": "(figure image; printed number)",
                     **meta(loc, q)})
    df = pd.DataFrame(rows)
    return df[["quantity", "value", "unit", "source_id", "locator", "verbatim", "verified", "machine_rechecked"]]


# ------------------------------------------------------------- supplements ----
def sdc1() -> tuple[pd.DataFrame, pd.DataFrame]:
    t = doc_tables("braak2011_sdc1.doc")[0]
    hdr = next(r for r in t if len(r) > 1 and r[1] == "No changes")
    cols = {lab: i for i, lab in enumerate(hdr) if lab}
    stage_rows, group_rows = [], []
    for r in t:
        if not r or not re.fullmatch(r"\d\d-\d+|total", r[0].strip()):
            continue
        label = r[0].strip()
        lo, hi = (1, 100) if label == "total" else band_of(label)
        for lab, stage in SDC1_STAGES.items():
            if lab == "No changes":
                n, _ = n_pct(r[cols[lab]])
            else:
                n = int(r[cols[lab]])
            if label != "total":
                stage_rows.append({"age_lower": lo, "age_upper": hi, "stage": stage, "n": n,
                                   "age_label_printed": label, **meta(
                                       "Supplemental Digital Content 1 (Supplementary Table 1)",
                                       f"Number of brains at AT8 stage '{lab}', age category {label} y")})
            else:
                group_rows.append({"age_label_printed": label, "age_lower": lo, "age_upper": hi,
                                   "column": lab, "measure": "n", "value": n, **meta(
                                       "Supplemental Digital Content 1 (Supplementary Table 1), total row",
                                       f"Number of brains at AT8 stage '{lab}', all ages")})
        for g in SDC1_GROUPS:
            n, pct = n_pct(r[cols[g]])
            for meas, v in (("n", n), ("pct", pct)):
                if v is None:
                    continue
                if g == "No changes" and meas == "n":
                    continue                      # already stored (stage '0' in the stage file; total row above)
                group_rows.append({"age_label_printed": label, "age_lower": lo, "age_upper": hi, "column": g,
                                   "measure": meas, "value": v, **meta(
                                       "Supplemental Digital Content 1 (Supplementary Table 1)"
                                       + (", total row" if label == "total" else ""),
                                       f"{'Number' if meas == 'n' else 'Row percent'} of brains in AT8 group "
                                       f"'{g}', age category {label}")})
        tot = r[cols["Total"]]
        n, pct = n_pct(tot)
        group_rows.append({"age_label_printed": label, "age_lower": lo, "age_upper": hi, "column": "Total",
                           "measure": "n", "value": n, **meta(
                               "Supplemental Digital Content 1 (Supplementary Table 1)",
                               f"Number of brains, age category {label}")})
    st = pd.DataFrame(stage_rows)
    return (st[["age_lower", "age_upper", "stage", "n", "source_id", "locator", "verified", "quantity",
                "machine_rechecked", "age_label_printed"]], pd.DataFrame(group_rows))


def sdc2() -> pd.DataFrame:
    t = doc_tables("braak2011_sdc2.doc")[0]
    hdr = next(r for r in t if r and r[0].startswith("Age groups"))
    phases = ["Phase 0", "Phase 1", "Phase 2", "Phase 3", "Phase 4", "Total"]
    rows = []
    for col, cell in zip(phases, hdr[1:]):
        m = re.search(r"n = (\d+) \(([\d.]+)\)", cell)
        for meas, v in (("n", m.group(1)), ("pct", m.group(2))):
            rows.append({"age_label_printed": "column header", "age_lower": 1, "age_upper": 100,
                         "amyloid_phase": col, "measure": meas, "value": v, **meta(
                             "Supplemental Digital Content 2 (Supplementary Table 2), column header",
                             f"{'Number' if meas == 'n' else 'Percent of all 2,332'} of brains in amyloid {col}")})
    for r in t:
        if not r or not re.fullmatch(r"\d\d-\d+", r[0].strip()):
            continue
        label = r[0].strip()
        lo, hi = band_of(label)
        for col, cell in zip(phases, r[1:]):
            n, pct = n_pct(cell)
            for meas, v in (("n", n), ("pct", pct)):
                rows.append({"age_label_printed": label, "age_lower": lo, "age_upper": hi, "amyloid_phase": col,
                             "measure": meas, "value": v, **meta(
                                 "Supplemental Digital Content 2 (Supplementary Table 2)",
                                 f"{'Number' if meas == 'n' else 'Row percent'} of brains in amyloid {col}, "
                                 f"age category {label} y")})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ checks ----
def checks(t1, t2, tc, st, g1, s2) -> list[tuple[str, bool, str]]:
    out = []

    def add(name, ok, detail):
        # "count" checks (sums of counts) must pass; "printed" checks compare a printed percentage or a
        # number quoted in the prose with the tables, where the source can disagree with itself.
        kind = "printed" if ("percent" in name or name.startswith(("Discussion", "Statistics"))) else "count"
        out.append((name, bool(ok), detail, kind))

    def v(df, **kw):
        m = pd.Series(True, index=df.index)
        if "age_lower" in kw and "age_label_printed" not in kw and "age_label_printed" in df.columns:
            m &= ~df["age_label_printed"].isin(["column header", "total"])     # band rows only
        for k, x in kw.items():
            m &= df[k].astype(str) == str(x)
        r = df[m]
        if len(r) != 1:
            raise KeyError(f"{kw}: {len(r)} rows")
        return float(r["value"].iloc[0])

    # Table 1
    for sex in ("female", "male", "total"):
        s = sum(v(t1, age_lower=lo, sex=sex, measure="n") for lo, _ in BANDS)
        add(f"Table 1: {sex} column sums to its header n", s == v(t1, age_label_printed="column header", sex=sex),
            f"sum {s:.0f} vs header {v(t1, age_label_printed='column header', sex=sex):.0f}")
    bad = [lo for lo, _ in BANDS if v(t1, age_lower=lo, sex="female", measure="n") + v(t1, age_lower=lo, sex="male", measure="n")
           != v(t1, age_lower=lo, sex="total", measure="n")]
    add("Table 1: female + male = total in every age row", not bad, f"rows failing: {bad}")
    worst = max(abs(100 * v(t1, age_lower=lo, sex=s, measure="n") / v(t1, age_lower=lo, sex="total", measure="n")
                    - v(t1, age_lower=lo, sex=s, measure="pct")) for lo, _ in BANDS for s in ("female", "male"))
    add("Table 1: printed row percentages match n / total (within 0.05 rounding)", worst <= 0.05 + 1e-9,
        f"largest difference {worst:.3f} points")
    # Table 2
    for col in T2_COLS:
        s = sum(v(t2, at8_stage=r, amyloid_phase=col, measure="n") for r in T2_ROWS)
        add(f"Table 2: {col} column sums to its header n", s == v(t2, at8_stage="column total", amyloid_phase=col, measure="n"),
            f"sum {s:.0f} vs header {v(t2, at8_stage='column total', amyloid_phase=col, measure='n'):.0f}")
    bad = [r for r in T2_ROWS if sum(v(t2, at8_stage=r, amyloid_phase=c, measure="n") for c in T2_COLS[:5])
           != v(t2, at8_stage=r, amyloid_phase="Total", measure="n")]
    add("Table 2: Phase 0..4 sum to the row total in every AT8 row", not bad, f"rows failing: {bad}")
    worst = max(abs(100 * v(t2, at8_stage=r, amyloid_phase=c, measure="n") / v(t2, at8_stage=r, amyloid_phase="Total", measure="n")
                    - v(t2, at8_stage=r, amyloid_phase=c, measure="pct")) for r in T2_ROWS for c in T2_COLS[:5])
    add("Table 2: printed row percentages match n / row total (within 0.05)", worst <= 0.05 + 1e-9,
        f"largest difference {worst:.3f} points")
    worst = max(abs(100 * v(t2, at8_stage="column total", amyloid_phase=c, measure="n") / 2332
                    - v(t2, at8_stage="column total", amyloid_phase=c, measure="pct")) for c in T2_COLS)
    add("Table 2: header percentages match n / 2,332 (within 0.05)", worst <= 0.05 + 1e-9, f"largest difference {worst:.3f}")
    # SDC 1
    band_tot = {lo: v(g1, age_lower=lo, column="Total", measure="n") for lo, _ in BANDS}
    bad = [lo for lo, _ in BANDS if st[st["age_lower"] == lo]["n"].sum() != band_tot[lo]]
    add("SDC 1: the 12 stage counts sum to the printed band total in every age row", not bad, f"rows failing: {bad}")
    bad = [lo for lo, _ in BANDS if band_tot[lo] != v(t1, age_lower=lo, sex="total", measure="n")]
    add("SDC 1 band totals equal Table 1 totals", not bad, f"rows failing: {bad}")
    grp = {"a-c": ["a", "b", "c"], "1a-1b": ["1a", "1b"], "I-II": ["I", "II"], "III-IV": ["III", "IV"],
           "V-VI": ["V", "VI"]}
    bad = []
    for lo, _ in BANDS:
        for g, ss in grp.items():
            if st[(st["age_lower"] == lo) & st["stage"].isin(ss)]["n"].sum() != v(g1, age_lower=lo, column=g, measure="n"):
                bad.append((lo, g))
    add("SDC 1: printed group subtotals equal the sum of their stages", not bad, f"failing: {bad}")
    for lab, stage in SDC1_STAGES.items():
        s = st[st["stage"] == stage]["n"].sum()
        add(f"SDC 1: stage {stage} column sums to the printed total row", s == v(g1, age_label_printed="total", column=lab, measure="n"),
            f"sum {s} vs total row {v(g1, age_label_printed='total', column=lab, measure='n'):.0f}")
    diffs = []
    for lo, _ in BANDS:
        for g in SDC1_GROUPS:
            sel = g1[(g1["age_lower"] == lo) & (g1["column"] == g) & (g1["measure"] == "pct")]
            if len(sel):
                n = (st[(st["age_lower"] == lo) & (st["stage"] == "0")]["n"].sum() if g == "No changes"
                     else v(g1, age_lower=lo, column=g, measure="n"))
                diffs.append((abs(100 * n / band_tot[lo] - float(sel["value"].iloc[0])), lo, g, sel["value"].iloc[0]))
    over = [(lo, g, p, round(d, 2)) for d, lo, g, p in diffs if d > 0.05 + 1e-9]
    add("SDC 1: printed row percentages match n / band total (within 0.05)", not over,
        f"{len(over)} cell(s) off by more than rounding: {over}")
    # SDC 2
    for col in ["Phase 0", "Phase 1", "Phase 2", "Phase 3", "Phase 4", "Total"]:
        s = sum(v(s2, age_lower=lo, amyloid_phase=col, measure="n") for lo, _ in BANDS)
        add(f"SDC 2: {col} column sums to its header n", s == v(s2, age_label_printed="column header", amyloid_phase=col, measure="n"),
            f"sum {s:.0f} vs header {v(s2, age_label_printed='column header', amyloid_phase=col, measure='n'):.0f}")
    bad = [lo for lo, _ in BANDS if sum(v(s2, age_lower=lo, amyloid_phase=f"Phase {k}", measure="n") for k in range(5))
           != v(s2, age_lower=lo, amyloid_phase="Total", measure="n")]
    add("SDC 2: phases sum to the band total in every age row", not bad, f"rows failing: {bad}")
    bad = [lo for lo, _ in BANDS if v(s2, age_lower=lo, amyloid_phase="Total", measure="n") != band_tot[lo]]
    add("SDC 2 band totals equal Table 1 / SDC 1 totals", not bad, f"rows failing: {bad}")
    worst = max((abs(100 * v(s2, age_lower=lo, amyloid_phase=f"Phase {k}", measure="n")
                     / v(s2, age_lower=lo, amyloid_phase="Total", measure="n")
                     - v(s2, age_lower=lo, amyloid_phase=f"Phase {k}", measure="pct")), lo, k)
                for lo, _ in BANDS for k in range(5))
    add("SDC 2: printed row percentages match n / band total (within 0.05)", worst[0] <= 0.05 + 1e-9,
        f"largest difference {worst[0]:.3f} points (band {worst[1]}, phase {worst[2]})")
    # Across sources
    for c in ["Phase 0", "Phase 1", "Phase 2", "Phase 3", "Phase 4"]:
        a, b = v(t2, at8_stage="column total", amyloid_phase=c, measure="n"), v(s2, age_label_printed="column header", amyloid_phase=c, measure="n")
        add(f"Table 2 and SDC 2 agree on the {c} total", a == b, f"{a:.0f} vs {b:.0f}")
    t2map = {"0": ["0"], "a-c": ["a", "b", "c"], "1a-1b": ["1a", "1b"], "I": ["I"], "II": ["II"], "III": ["III"],
             "IV": ["IV"], "V": ["V"], "VI": ["VI"]}
    bad = [r for r, ss in t2map.items() if st[st["stage"].isin(ss)]["n"].sum() != v(t2, at8_stage=r, amyloid_phase="Total", measure="n")]
    add("Table 2 row totals equal SDC 1 stage totals", not bad, f"rows failing: {bad}")
    txt_map = {"Pretangle stage a cases": ["a"], "Pretangle stage b cases": ["b"], "Pretangle stage c cases": ["c"],
               "Pretangle stage 1a cases": ["1a"], "Pretangle stage 1b cases": ["1b"], "NFT stage I cases": ["I"],
               "NFT stage II cases": ["II"], "NFT stage III cases": ["III"], "NFT stage IV cases": ["IV"],
               "NFT stage V cases": ["V"], "NFT stage VI cases": ["VI"],
               "Cases devoid of AD-related intraneuronal tau (no changes)": ["0"],
               "Pretangle stages a-c cases (no cortical AT8-ir projection neurons)": ["a", "b", "c"],
               "Pretangle stages 1a and 1b cases": ["1a", "1b"], "NFT stages V and VI cases": ["V", "VI"]}
    bad = []
    for q, ss in txt_map.items():
        tv = float(tc.loc[tc["quantity"] == q, "value"].iloc[0])
        if tv != st[st["stage"].isin(ss)]["n"].sum():
            bad.append((q, tv))
    add("Counts in the Results text equal SDC 1 stage totals", not bad, f"failing: {bad}")
    plaques = float(tc.loc[tc["quantity"] == "Cases with beta-amyloid plaques", "value"].iloc[0])
    add("Results text: plaque cases = 2,332 - Phase 0 total", plaques == 2332 - v(t2, at8_stage="column total", amyloid_phase="Phase 0", measure="n"),
        f"{plaques:.0f} vs {2332 - v(t2, at8_stage='column total', amyloid_phase='Phase 0', measure='n'):.0f}")
    disc = float(tc.loc[tc["quantity"].str.startswith("Cases with abnormal tau in pretangle"), "value"].iloc[0])
    p0 = v(t2, at8_stage="column total", amyloid_phase="Phase 0", measure="n") - v(t2, at8_stage="0", amyloid_phase="Phase 0", measure="n")
    add("Discussion: 1,291 tau-positive cases without plaques = Phase 0 total minus the 10 no-change cases",
        disc == p0, f"{disc:.0f} vs {p0:.0f}")
    d52 = float(tc.loc[tc["quantity"].str.startswith("Younger-age cases with abnormal tau only"), "value"].iloc[0])
    sa = st[st["stage"] == "a"]["n"].sum()
    add("Discussion '52 cases ... (pretangle Stage a)' equals the stage a count (Results: 40)", d52 == sa,
        f"Discussion says {d52:.0f}; Results and SDC 1 give stage a = {sa}; stage a + b = "
        f"{sa + st[st['stage'] == 'b']['n'].sum()}. Recorded as printed, not used.")
    m80 = sum(v(t1, age_lower=lo, sex="male", measure="n") for lo, _ in BANDS[:8])
    a80 = sum(v(t1, age_lower=lo, sex="total", measure="n") for lo, _ in BANDS[:8])
    f80 = sum(v(t1, age_lower=lo, sex="female", measure="n") for lo, _ in BANDS[8:])
    o80 = sum(v(t1, age_lower=lo, sex="total", measure="n") for lo, _ in BANDS[8:])
    add("Statistics text: 58% males at ages 1-80 and 54.7% females above 80 match Table 1",
        round(100 * m80 / a80) == 58 and round(100 * f80 / o80, 1) == 54.7,
        f"males 1-80: {100 * m80 / a80:.1f}%; females 81-100: {100 * f80 / o80:.1f}%")
    fig = {f"Figure 2A, n printed above the {lo}-{hi} y bar": lo for lo, hi in BANDS}
    bad = [q for q, lo in fig.items() if float(tc.loc[tc["quantity"] == q, "value"].iloc[0]) != band_tot[lo]]
    add("Figure 2 bar labels equal the band totals", not bad, f"failing: {bad}")
    return out


def main() -> int:
    txt = pdf_text()
    t1, t2, tc = table1(txt), table2(txt), text_counts(txt)
    st, g1 = sdc1()
    s2 = sdc2()
    for df, name in ((t1, "braak2011_table1.csv"), (t2, "braak2011_table2.csv"), (tc, "braak2011_text_counts.csv"),
                     (st, "braak2011_stage_counts.csv"), (g1, "braak2011_sdc1_groups.csv"),
                     (s2, "braak2011_sdc2_amyloid_by_age.csv")):
        df.to_csv(OUT / name, index=False)
        print(f"{name}: {len(df)} rows")
    res = checks(t1, t2, tc, st, g1, s2)
    lines = ["# Braak et al. 2011: internal consistency checks of the extracted values", "",
             "Generated by `python -m src.extraction.braak2011`. Each check compares numbers printed in "
             "different places of the paper and its supplements. A failure means the source disagrees with "
             "itself or the extraction is wrong; it is reported, never corrected silently.", "",
             "Kinds: **count** = sums of counts (must all pass); **printed** = a printed percentage or a number "
             "quoted in the prose compared with the tables (a failure is the source disagreeing with itself; "
             "the value is stored as printed and the failure is listed in the report).", "",
             "| check | kind | result | detail |", "|---|---|---|---|"]
    lines += [f"| {n} | {k} | {'pass' if ok else ('**FAIL**' if k == 'count' else '**source disagrees**')} | {d} |"
              for n, ok, d, k in res]
    (OUT / "braak2011_sum_checks.md").write_text("\n".join(lines) + "\n")
    n_fail = sum(not ok for _, ok, _, _ in res)
    n_count_fail = sum(not ok for _, ok, _, k in res if k == "count")
    print(f"{len(res) - n_fail}/{len(res)} consistency checks pass; count checks failing: {n_count_fail}")
    for n, ok, d, k in res:
        if not ok:
            print(f"{'FAIL' if k == 'count' else 'SOURCE DISAGREES'} {n}: {d}")
    return 1 if n_count_fail else 0


if __name__ == "__main__":
    sys.exit(main())
