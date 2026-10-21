"""Human verification of extracted values: a printable checklist and a tick script.

Two files are covered: data/extracted/braak2011_stage_counts.csv (Module 2) and
data/extracted/meta_extraction.csv (Module 5). Every row carries
``verified: false`` until the researcher has checked it against the source.

    python -m src.extraction.verification checklist
        writes, under data/extracted/verification/:
          checklist.html            printable checklist (one box per row)
          checklist.pdf             the same, as a PDF (if Chromium is available)
          ticks_braak2011_stage_counts.csv, ticks_meta_extraction.csv
                                    tick sheets: put an x in `tick` for each
                                    row you have checked

    python -m src.extraction.verification apply TICKS.csv [--apply]
        sets verified: true for the ticked rows only. Without --apply it only
        shows what would change.

Safeguards in `apply`:
  * only `verified` changes, and only from false to true. Nothing is ever set
    to false, and no other cell is touched (the CSV is rewritten byte for
    byte except for those cells);
  * each tick-sheet row carries a fingerprint of the row as printed. If the
    row has changed since the checklist was printed, nothing is written;
  * any tick value other than blank or one of x / yes / y / true / 1 / ✓
    stops the run, as does an unknown row ID. All or nothing: either every
    ticked row is applied or none is;
  * every applied sheet is logged in data/extracted/verification/applied_log.md.

Claude never runs `apply --apply` on the researcher's behalf.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from src.common.paths import REPO_ROOT
from src.common.provenance import iso_utc, sha256_file, utc_now

EXTRACTED = REPO_ROOT / "data" / "extracted"
OUT_DIR = EXTRACTED / "verification"
TICK_YES = {"x", "yes", "y", "true", "1", "✓", "✔"}
SKIP_IN_FINGERPRINT = {"verified", "machine_rechecked"}


@dataclass(frozen=True)
class Target:
    key: str
    path: Path
    prefix: str


TARGETS = {
    "braak2011_stage_counts": Target("braak2011_stage_counts", EXTRACTED / "braak2011_stage_counts.csv", "B"),
    "meta_extraction": Target("meta_extraction", EXTRACTED / "meta_extraction.csv", "M"),
}


# --------------------------------------------------------------------------- CSV I/O
def read_rows(path: Path) -> tuple[list[str], list[list[str]], str]:
    """Header, data rows and line ending, read so that write_rows reproduces the file exactly."""
    raw = path.read_text(encoding="utf-8")
    le = "\r\n" if "\r\n" in raw else "\n"
    rows = list(csv.reader(io.StringIO(raw, newline="")))
    if _serialise(rows, le) != raw:
        raise ValueError(f"{path}: cannot be rewritten byte for byte; refusing to touch it")
    return rows[0], rows[1:], le


def _serialise(rows: list[list[str]], le: str) -> str:
    buf = io.StringIO()
    csv.writer(buf, lineterminator=le, quoting=csv.QUOTE_MINIMAL).writerows(rows)
    return buf.getvalue()


def row_id(target: Target, i: int, n: int) -> str:
    width = max(2, len(str(n)))
    return f"{target.prefix}-{i + 1:0{width}d}"


def fingerprint(header: list[str], row: list[str]) -> str:
    """Short hash of the row's content, excluding the verification columns."""
    parts = [f"{h}={v}" for h, v in zip(header, row) if h not in SKIP_IN_FINGERPRINT]
    return hashlib.sha256("\x1f".join(parts).encode("utf-8")).hexdigest()[:10]


def describe(target: Target, rec: dict) -> str:
    if target.key == "braak2011_stage_counts":
        return f"age {rec['age_lower']}-{rec['age_upper']}, stage {rec['stage']}: n = {rec['n']}"
    return (f"{rec['study_id']}: {rec['drug_class']} vs {rec['comparator_type']} comparator, "
            f"{rec['effect_measure']} {rec['estimate']} ({rec['ci_lower']}-{rec['ci_upper']})")


# --------------------------------------------------------------------------- checklist
def write_tick_sheet(target: Target, out_dir: Path) -> Path:
    header, rows, _ = read_rows(target.path)
    out = out_dir / f"ticks_{target.key}.csv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["id", "fingerprint", "currently_verified", "what", "tick", "note"])
        for i, r in enumerate(rows):
            rec = dict(zip(header, r))
            w.writerow([row_id(target, i, len(rows)), fingerprint(header, r), rec["verified"],
                        describe(target, rec), "", ""])
    return out


def _pct(level: str) -> str:
    try:
        x = float(level)
    except ValueError:
        return level
    return f"{x * 100:g}%" if x <= 1 else f"{x:g}%"


def _e(s: object) -> str:
    return html.escape(str(s))


def braak_html(target: Target) -> str:
    header, rows, _ = read_rows(target.path)
    recs = [dict(zip(header, r)) for r in rows]
    ids = {(r["age_lower"], r["age_upper"], r["stage"]): (row_id(target, i, len(rows)), fingerprint(header, rows[i]))
           for i, r in enumerate(recs)}
    stages = list(dict.fromkeys(r["stage"] for r in recs))
    bands = list(dict.fromkeys((r["age_lower"], r["age_upper"], r.get("age_label_printed", "")) for r in recs))
    val = {(r["age_lower"], r["age_upper"], r["stage"]): int(r["n"]) for r in recs}
    stage_label = {"0": "No changes"}
    head = "".join(f"<th>{_e(stage_label.get(s, s))}</th>" for s in stages)
    body = []
    for lo, hi, printed in bands:
        cells = []
        for s in stages:
            rid, _fp = ids[(lo, hi, s)]
            cells.append(f'<td><span class="v">{val[(lo, hi, s)]}</span><span class="box">☐</span>'
                         f'<span class="id">{rid}</span></td>')
        total = sum(val[(lo, hi, s)] for s in stages)
        body.append(f"<tr><th>{_e(printed or f'{lo}-{hi}')}</th>{''.join(cells)}<td class='tot'>{total}</td></tr>")
    col_tot = "".join(f"<td class='tot'>{sum(val[(lo, hi, s)] for lo, hi, _ in bands)}</td>" for s in stages)
    grand = sum(val.values())
    src = recs[0]
    return f"""
<section class="file">
<h2>1. Braak et al. 2011 — tau stage counts by age band</h2>
<p class="meta"><b>File:</b> <code>data/extracted/braak2011_stage_counts.csv</code> ({len(rows)} rows) ·
<b>Source:</b> {_e(src['locator'])}, <code>http://links.lww.com/NEN/A273</code> ·
<b>Used by:</b> Module 2 (all fits)</p>
<p class="how">Open Supplementary Table 1 next to this sheet. For each cell, check that the number printed
here equals the number in the supplement (same age row, same stage column), then tick its box. The grey
code under each box (B-…) is the row ID on the tick sheet <code>ticks_braak2011_stage_counts.csv</code>.
Row and column totals are computed from the file, so you can also compare them with the supplement's
totals. Abstract totals for reference: 2,332 brains; 10 without AT8 tau; 58 at stages a–c; 38 at 1a; 236 at 1b.</p>
<div class="scroll"><table class="grid">
<thead><tr><th>Age (printed)</th>{head}<th>Total</th></tr></thead>
<tbody>{''.join(body)}
<tr><th>Total</th>{col_tot}<td class="tot">{grand}</td></tr></tbody></table></div>
<p class="sign">Checked by: ____________________ &nbsp; Date: ____________ &nbsp; Notes: ______________________________</p>
</section>"""


def meta_html(target: Target) -> str:
    header, rows, _ = read_rows(target.path)
    cards = []
    for i, r in enumerate(rows):
        rec = dict(zip(header, r))
        rid = row_id(target, i, len(rows))
        inv = " · <b>invert_ratio: true</b> (the code flips the printed ratio; check the printed numbers as printed)" \
            if rec.get("invert_ratio", "").lower() == "true" else ""
        excl = (f'<div class="excl">Excluded from the analysis: {_e(rec["exclusion_reason"])}</div>'
                if rec.get("include", "").lower() != "true" else "")
        cards.append(f"""
<div class="card">
  <div class="cardhead"><span class="box big">☐</span><b>{rid}</b> · {_e(rec['study_id'])} ·
  source <code>{_e(rec['source_id'])}</code> · {_e(rec['locator'])}</div>
  <table class="kv">
   <tr><th>Exposure</th><td>{_e(rec['drug_class'])} ({_e(rec['drug'])})</td></tr>
   <tr><th>Comparator</th><td>{_e(rec['comparator_type'])}: {_e(rec['comparator'])}</td></tr>
   <tr><th>Outcome</th><td>{_e(rec['outcome'])}</td></tr>
   <tr><th>Number to check</th><td class="num">{_e(rec['effect_measure'])} {_e(rec['estimate'])}
       ({_e(_pct(rec['ci_level']))} CI {_e(rec['ci_lower'])}–{_e(rec['ci_upper'])}){inv}</td></tr>
   <tr><th>Verbatim in file</th><td class="verb">{_e(rec['verbatim'])}</td></tr>
  </table>{excl}
  <div class="checks">☐ estimate &nbsp; ☐ lower CI &nbsp; ☐ upper CI &nbsp; ☐ measure (HR/RR/OR/IRR) &nbsp;
  ☐ exposure and comparator the right way round &nbsp; Notes: ____________________</div>
</div>""")
    return f"""
<section class="file">
<h2>2. Meta-analysis extraction</h2>
<p class="meta"><b>File:</b> <code>data/extracted/meta_extraction.csv</code> ({len(rows)} rows) ·
<b>Used by:</b> Module 5. Citations and licences: <code>data/extracted/meta_extraction.sources.md</code></p>
<p class="how">For each row, open the paper at the place given (PMCID and locator), compare the numbers with
the paper, and tick the big box only if all of them match. The small boxes are for your own working. The
row ID (M-…) is the line on the tick sheet <code>ticks_meta_extraction.csv</code>. Excluded rows are listed
too: their numbers are kept in the file, so they can be verified like the others.</p>
{''.join(cards)}
<p class="sign">Checked by: ____________________ &nbsp; Date: ____________</p>
</section>"""


CSS = """
:root { --ink:#1d2321; --muted:#5c6763; --rule:#c9d1cd; --bg:#ffffff; --accent:#2f5d50; --warn:#8a3b12; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --ink:#e3e8e5; --muted:#a3ada8; --rule:#3c4542; --bg:#151917; --accent:#8cc3b0; --warn:#f0a67a; color-scheme: dark } }
:root[data-theme="dark"] { --ink:#e3e8e5; --muted:#a3ada8; --rule:#3c4542; --bg:#151917; --accent:#8cc3b0; --warn:#f0a67a; color-scheme: dark }
body { font: 11pt/1.45 Georgia, 'Times New Roman', serif; color: var(--ink); background: var(--bg);
       max-width: 60rem; margin: 0 auto; padding-block: 1.5rem; padding-inline: 16px; }
h1 { font-size: 1.5rem; margin: 0 0 .3rem; } h2 { font-size: 1.15rem; margin: 1.6rem 0 .4rem; color: var(--accent); }
.meta, .how { font-size: .9rem; } .how { color: var(--muted); }
code { font: .85em ui-monospace, Menlo, Consolas, monospace; }
.scroll { overflow-x: auto; }
table.grid { border-collapse: collapse; font: .82rem/1.2 ui-monospace, Menlo, Consolas, monospace; font-variant-numeric: tabular-nums; }
table.grid th, table.grid td { border: 1px solid var(--rule); padding: .25rem .35rem; text-align: center; vertical-align: top; }
table.grid td .v { display: block; font-size: 1rem; font-weight: 600; }
table.grid td .box { display: block; font-size: 1.05rem; }
table.grid td .id { display: block; font-size: .62rem; color: var(--muted); }
table.grid .tot { color: var(--muted); font-weight: 400; }
.card { border: 1px solid var(--rule); border-radius: 4px; padding: .55rem .7rem; margin: .6rem 0; break-inside: avoid; }
.cardhead { font-size: .92rem; margin-bottom: .3rem; } .box.big { font-size: 1.35rem; margin-right: .4rem; vertical-align: -2px; }
table.kv { border-collapse: collapse; width: 100%; font-size: .85rem; }
table.kv th { text-align: left; color: var(--muted); font-weight: 400; width: 9.5rem; vertical-align: top; padding: .1rem .5rem .1rem 0; }
table.kv td { padding: .1rem 0; } td.num { font-weight: 600; font-variant-numeric: tabular-nums; }
td.verb { font-style: italic; } .excl { color: var(--warn); font-size: .85rem; margin-top: .25rem; }
.checks { font-size: .8rem; color: var(--muted); margin-top: .35rem; } .sign { margin-top: .8rem; font-size: .9rem; }
@media print { body { max-width: none; font-size: 9.5pt; } h2 { break-before: page; } section.file:first-of-type h2 { break-before: auto; } }
"""


def checklist(out_dir: Path = OUT_DIR, pdf: bool = True) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written = [write_tick_sheet(t, out_dir) for t in TARGETS.values()] + [write_combined_sheet(out_dir)]
    stamp = iso_utc(utc_now())
    hashes = " · ".join(f"<code>{t.path.name}</code> sha256 {sha256_file(t.path)[:12]}…" for t in TARGETS.values())
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Extraction Verification Checklist</title><style>{CSS}</style></head><body>
<h1>Verification checklist — extracted values</h1>
<p class="meta">Printed from the files as they stood at {stamp} ({hashes}).
Every value below was extracted with Claude's help and is marked <code>verified: false</code>.</p>
<p class="how"><b>How this works.</b> Check each value against the source and tick its box here. Then put an
<b>x</b> in the <code>tick</code> column of the matching tick sheet (same ID) and run
<code>python -m src.extraction.verification apply data/extracted/verification/ticks_&lt;file&gt;.csv</code>
to see what would change, then the same command with <code>--apply</code>. Only rows you tick become
<code>verified: true</code>. If a row has changed since this sheet was printed, the script refuses and
asks for a fresh checklist.</p>
{braak_html(TARGETS['braak2011_stage_counts'])}
{meta_html(TARGETS['meta_extraction'])}
</body></html>"""
    html_path = out_dir / "checklist.html"
    html_path.write_text(doc, encoding="utf-8")
    written.append(html_path)
    if pdf:
        p = print_pdf(html_path, out_dir / "checklist.pdf")
        if p:
            written.append(p)
    return written


def _chromium() -> str | None:
    for name in ("chromium", "chromium-browser", "google-chrome"):
        if shutil.which(name):
            return shutil.which(name)
    for p in sorted(Path("/opt/pw-browsers").glob("chromium*/chrome-linux*/chrome")):
        return str(p)
    return None


def print_pdf(html_path: Path, pdf_path: Path) -> Path | None:
    exe = _chromium()
    if not exe:
        print("no Chromium found: checklist.pdf not written (print checklist.html from a browser instead)")
        return None
    subprocess.run([exe, "--headless", "--no-sandbox", "--disable-gpu", "--no-pdf-header-footer",
                    f"--print-to-pdf={pdf_path}", html_path.as_uri()], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
    return pdf_path


# --------------------------------------------------------------------------- apply
class TickError(ValueError):
    pass


def _sheet_rows(tick_path: Path) -> list[tuple[int, dict]]:
    with tick_path.open(encoding="utf-8-sig", newline="") as fh:
        return list(enumerate(csv.DictReader(fh), start=2))


def plan(tick_path: Path, target: Target, sheet: list[tuple[int, dict]] | None = None) -> list[tuple[int, str]]:
    """(row index, row id) for every ticked row; raises TickError on anything doubtful."""
    header, rows, _ = read_rows(target.path)
    by_id = {row_id(target, i, len(rows)): i for i in range(len(rows))}
    problems, todo = [], []
    if sheet is None:
        sheet = _sheet_rows(tick_path)
    if True:
        for line_no, t in sheet:
            tick = (t.get("tick") or "").strip().lower()
            if not tick:
                continue
            if tick not in TICK_YES:
                problems.append(f"line {line_no}: tick value {t.get('tick')!r} is not one of {sorted(TICK_YES)}")
                continue
            rid = (t.get("id") or "").strip()
            if rid not in by_id:
                problems.append(f"line {line_no}: unknown row id {rid!r} for {target.path.name}")
                continue
            i = by_id[rid]
            if fingerprint(header, rows[i]) != (t.get("fingerprint") or "").strip().removeprefix("fp-"):
                problems.append(f"line {line_no}: {rid} has changed since the checklist was printed "
                                f"(fingerprint {t.get('fingerprint')!r} != {fingerprint(header, rows[i])!r}); "
                                "print a fresh checklist and check it again")
                continue
            todo.append((i, rid))
    if problems:
        raise TickError("nothing written:\n  " + "\n  ".join(problems))
    return todo


def target_for(tick_path: Path) -> Target:
    for t in TARGETS.values():
        if tick_path.name == f"ticks_{t.key}.csv" or tick_path.name.startswith(f"ticks_{t.key}"):
            return t
    raise TickError(f"{tick_path.name}: cannot tell which file it belongs to (expected ticks_<file>.csv, "
                    f"one of {[f'ticks_{k}.csv' for k in TARGETS]})")


def apply(tick_path: Path, *, write: bool, target: Target | None = None,
          log_path: Path | None = None) -> list[str]:
    target = target or target_for(tick_path)
    return _apply_planned(tick_path, target, plan(tick_path, target), write=write, log_path=log_path)


def _apply_planned(tick_path: Path, target: Target, todo: list[tuple[int, str]], *, write: bool,
                   log_path: Path | None = None) -> list[str]:
    header, rows, le = read_rows(target.path)
    col = header.index("verified")
    flipped, already = [], []
    for i, rid in todo:
        if rows[i][col].strip().lower() == "true":
            already.append(rid)
        else:
            rows[i][col] = "true"
            flipped.append(rid)
    print(f"{target.path.name}: {len(todo)} ticked; {len(flipped)} would become verified: true; "
          f"{len(already)} already true; {len(rows) - len(todo)} untouched")
    for rid in flipped:
        print(f"  {rid} -> verified: true")
    if not write:
        print("dry run: nothing written (add --apply to write)")
        return flipped
    if flipped:
        target.path.write_text(_serialise([header] + rows, le), encoding="utf-8", newline="")
        log_path = log_path or (OUT_DIR / "applied_log.md")
        log_path.parent.mkdir(parents=True, exist_ok=True)
        new = not log_path.exists()
        with log_path.open("a", encoding="utf-8") as fh:
            if new:
                fh.write("# Verification ticks applied\n\n| when (UTC) | file | tick sheet sha256 | rows set to verified: true |\n|---|---|---|---|\n")
            fh.write(f"| {iso_utc(utc_now())} | {target.path.name} | {sha256_file(tick_path)[:16]} | "
                     f"{len(flipped)}: {', '.join(flipped)} |\n")
        print(f"wrote {target.path} and logged in {log_path}")
    return flipped


# --------------------------------------------------------------------------- one combined sheet
COMBINED = "verification_sheet.csv"
COMBINED_COLS = ["id", "file", "fingerprint", "currently_verified", "what_to_check", "value_in_file",
                 "source", "where_in_source", "link", "verbatim_or_description", "tick", "note"]


def _braak_rows(target: Target):
    header, rows, _ = read_rows(target.path)
    for i, r in enumerate(rows):
        rec = dict(zip(header, r))
        stage = "No changes (stage 0)" if rec["stage"] == "0" else f"stage {rec['stage']}"
        yield {"id": row_id(target, i, len(rows)), "file": target.path.name, "fingerprint": "fp-" + fingerprint(header, r),
               "currently_verified": rec["verified"],
               "what_to_check": f"number of brains, age {rec.get('age_label_printed') or rec['age_lower'] + '-' + rec['age_upper']}, {stage}",
               "value_in_file": rec["n"],
               "source": "Braak H, Thal DR, Ghebremedhin E, Del Tredici K. J Neuropathol Exp Neurol 2011;70:960-969",
               "where_in_source": f"{rec['locator']}: row 'age {rec.get('age_label_printed', '')}', column '{rec['stage'] if rec['stage'] != '0' else 'No changes'}'",
               "link": "http://links.lww.com/NEN/A273",
               "verbatim_or_description": rec.get("quantity", ""), "tick": "", "note": ""}


def _meta_rows(target: Target):
    header, rows, _ = read_rows(target.path)
    for i, r in enumerate(rows):
        rec = dict(zip(header, r))
        sid = rec["source_id"]
        link = f"https://pmc.ncbi.nlm.nih.gov/articles/{sid}/" if sid.upper().startswith("PMC") else ""
        inv = "; invert_ratio = true (the code flips it; check the numbers as printed)" \
            if rec.get("invert_ratio", "").lower() == "true" else ""
        excl = "" if rec.get("include", "").lower() == "true" else f"; EXCLUDED from analysis: {rec.get('exclusion_reason', '')}"
        yield {"id": row_id(target, i, len(rows)), "file": target.path.name, "fingerprint": "fp-" + fingerprint(header, r),
               "currently_verified": rec["verified"],
               "what_to_check": (f"{rec['study_id']}: {rec['drug_class']} ({rec['drug']}) vs {rec['comparator_type']} "
                                 f"comparator ({rec['comparator']}); outcome {rec['outcome']}; also check the measure "
                                 f"and that exposure/comparator are the right way round{inv}{excl}"),
               "value_in_file": (f"{rec['effect_measure']} {rec['estimate']} ({_pct(rec['ci_level'])} CI "
                                 f"{rec['ci_lower']} to {rec['ci_upper']})"),
               "source": sid, "where_in_source": rec["locator"], "link": link,
               "verbatim_or_description": rec["verbatim"], "tick": "", "note": ""}


def write_combined_sheet(out_dir: Path = OUT_DIR) -> Path:
    """One CSV for checking on a laptop: the 120 Braak counts first, then the
    37 meta-analysis rows. Put x in `tick` for each row checked."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / COMBINED
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COMBINED_COLS, lineterminator="\n")
        w.writeheader()
        w.writerows(_braak_rows(TARGETS["braak2011_stage_counts"]))
        w.writerows(_meta_rows(TARGETS["meta_extraction"]))
    return out


def apply_combined(tick_path: Path, *, write: bool, targets: dict[str, Target] | None = None,
                   log_path: Path | None = None) -> dict[str, list[str]]:
    """Apply a sheet holding rows for several files (routed by row-ID prefix).
    Every file is planned before any is written: all or nothing."""
    targets = targets or TARGETS
    sheet = _sheet_rows(tick_path)
    by_prefix = {t.prefix: t for t in targets.values()}
    groups: dict[str, list] = {k: [] for k in targets}
    bad = []
    for line_no, t in sheet:
        pre = (t.get("id") or "").split("-")[0]
        if pre not in by_prefix:
            if (t.get("tick") or "").strip():
                bad.append(f"line {line_no}: row id {t.get('id')!r} belongs to no known file")
            continue
        groups[by_prefix[pre].key].append((line_no, t))
    if bad:
        raise TickError("nothing written:\n  " + "\n  ".join(bad))
    plans = {k: plan(tick_path, targets[k], groups[k]) for k in targets}   # raises before any write
    out = {}
    for k, todo in plans.items():
        out[k] = _apply_planned(tick_path, targets[k], todo, write=write, log_path=log_path)
    return out


def _main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m src.extraction.verification", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("checklist", help="write the printable checklist and the tick sheets")
    c.add_argument("--no-pdf", action="store_true")
    a = sub.add_parser("apply", help="set verified: true for ticked rows")
    a.add_argument("ticks", type=Path)
    a.add_argument("--apply", action="store_true", help="write the change (default: dry run)")
    args = ap.parse_args(argv)
    if args.cmd == "checklist":
        for p in checklist(pdf=not args.no_pdf):
            print(f"wrote {p.relative_to(REPO_ROOT)}")
        return 0
    try:
        if "id" in (_sheet_rows(args.ticks)[0][1] if _sheet_rows(args.ticks) else {}) and \
                "file" in _sheet_rows(args.ticks)[0][1]:
            apply_combined(args.ticks, write=args.apply)
        else:
            apply(args.ticks, write=args.apply)
    except TickError as e:
        print(str(e), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(_main())
