"""The tick script sets verified: true for ticked rows only, and nothing else."""
import csv
import shutil

import pytest

from src.extraction import verification as v


@pytest.fixture(params=["braak2011_stage_counts", "meta_extraction"])
def setup(request, tmp_path):
    real = v.TARGETS[request.param]
    data = tmp_path / real.path.name
    shutil.copy(real.path, data)
    target = v.Target(real.key, data, real.prefix)
    sheet = v.write_tick_sheet(target, tmp_path)
    return target, sheet, tmp_path


def _tick(sheet, ids, value="x", fingerprint=None):
    rows = list(csv.DictReader(sheet.open(encoding="utf-8", newline="")))
    for r in rows:
        if r["id"] in ids:
            r["tick"] = value
            if fingerprint is not None:
                r["fingerprint"] = fingerprint
    with sheet.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return [r["id"] for r in rows]


def _verified(target):
    header, rows, _ = v.read_rows(target.path)
    i = header.index("verified")
    return [r[i] for r in rows]


def test_only_ticked_rows_flip_and_nothing_else_changes(setup):
    target, sheet, tmp = setup
    before_raw = target.path.read_text(encoding="utf-8")
    header, before_rows, _ = v.read_rows(target.path)
    ids = _tick(sheet, set())
    chosen = {ids[0], ids[2]}
    _tick(sheet, chosen)
    flipped = v.apply(sheet, write=True, target=target, log_path=tmp / "log.md")
    assert set(flipped) == chosen
    after = _verified(target)
    for k, rid in enumerate(ids):
        assert after[k] == ("true" if rid in chosen else "false")
    _, after_rows, _ = v.read_rows(target.path)
    col = header.index("verified")
    for b, a in zip(before_rows, after_rows):
        assert [x for j, x in enumerate(b) if j != col] == [x for j, x in enumerate(a) if j != col]
    # byte-for-byte: undoing the flipped cells restores the original file exactly
    for k, rid in enumerate(ids):
        if rid in chosen:
            after_rows[k][col] = "false"
    assert v._serialise([header] + after_rows, "\r\n" if "\r\n" in before_raw else "\n") == before_raw
    assert (tmp / "log.md").exists()


def test_dry_run_writes_nothing(setup):
    target, sheet, tmp = setup
    raw = target.path.read_bytes()
    ids = _tick(sheet, set())
    _tick(sheet, {ids[1]})
    assert v.apply(sheet, write=False, target=target, log_path=tmp / "log.md") == [ids[1]]
    assert target.path.read_bytes() == raw and not (tmp / "log.md").exists()


def test_changed_row_or_odd_tick_writes_nothing(setup):
    target, sheet, tmp = setup
    raw = target.path.read_bytes()
    ids = _tick(sheet, set())
    _tick(sheet, {ids[0]}, fingerprint="0000000000")
    _tick(sheet, {ids[1]})
    with pytest.raises(v.TickError, match="changed since the checklist was printed"):
        v.apply(sheet, write=True, target=target, log_path=tmp / "log.md")
    assert target.path.read_bytes() == raw
    _tick(sheet, {ids[0]}, value="maybe", fingerprint=v.fingerprint(*_hdr_row(target, 0)))
    with pytest.raises(v.TickError, match="not one of"):
        v.apply(sheet, write=True, target=target, log_path=tmp / "log.md")
    assert target.path.read_bytes() == raw


def _hdr_row(target, i):
    header, rows, _ = v.read_rows(target.path)
    return header, rows[i]


def test_blank_sheet_changes_nothing_and_never_sets_false(setup):
    target, sheet, tmp = setup
    raw = target.path.read_bytes()
    assert v.apply(sheet, write=True, target=target, log_path=tmp / "log.md") == []
    assert target.path.read_bytes() == raw



def test_combined_sheet_routes_ticks_to_each_file_and_is_all_or_nothing(tmp_path):
    targets = {}
    for key, real in v.TARGETS.items():
        dst = tmp_path / real.path.name
        shutil.copy(real.path, dst)
        targets[key] = v.Target(real.key, dst, real.prefix)
    sheet = tmp_path / v.COMBINED
    with sheet.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=v.COMBINED_COLS, lineterminator="\n")
        w.writeheader()
        w.writerows(v._braak_rows(targets["braak2011_stage_counts"]))
        w.writerows(v._meta_rows(targets["meta_extraction"]))
    rows = list(csv.DictReader(sheet.open(encoding="utf-8")))
    assert [r["id"][0] for r in rows].count("B") == 120 and [r["id"][0] for r in rows].count("M") == 37
    assert rows[0]["id"] == "B-001" and rows[120]["id"] == "M-01"
    for r in rows:
        if r["id"] in ("B-001", "B-050", "M-02"):
            r["tick"] = "x"
    raw = {k: t.path.read_bytes() for k, t in targets.items()}
    bad = [dict(r) for r in rows]
    bad[1]["tick"] = "maybe"                     # one unclear tick: nothing written anywhere
    with sheet.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=v.COMBINED_COLS, lineterminator="\n"); w.writeheader(); w.writerows(bad)
    with pytest.raises(v.TickError):
        v.apply_combined(sheet, write=True, targets=targets, log_path=tmp_path / "log.md")
    assert {k: t.path.read_bytes() for k, t in targets.items()} == raw
    with sheet.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=v.COMBINED_COLS, lineterminator="\n"); w.writeheader(); w.writerows(rows)
    out = v.apply_combined(sheet, write=True, targets=targets, log_path=tmp_path / "log.md")
    assert out == {"braak2011_stage_counts": ["B-001", "B-050"], "meta_extraction": ["M-02"]}
    assert _verified(targets["braak2011_stage_counts"]).count("true") == 2
    assert _verified(targets["meta_extraction"]).count("true") == 1
