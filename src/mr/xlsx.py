"""Minimal, standard-library reader for the sheets of an .xlsx workbook.

The supplementary tables of the pQTL / brain-eQTL papers are distributed as
.xlsx files. openpyxl is not in the locked Python environment, and adding a
dependency only to read a few published tables is not worth it, so this reads
the Office Open XML directly: shared strings, then each row of one sheet,
streamed with iterparse so a large sheet does not have to fit in memory as a
tree. Cell values are returned as strings exactly as stored (numbers in their
stored decimal form); empty cells are "".

    from src.mr.xlsx import sheet_names, iter_rows
    for row in iter_rows(path, "ST9"): ...
"""

from __future__ import annotations

import re
import zipfile
from collections.abc import Iterator
from pathlib import Path
from xml.etree.ElementTree import iterparse

_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _col_index(ref: str) -> int:
    letters = re.match(r"[A-Z]+", ref).group(0)
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def _shared_strings(z: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    out = []
    with z.open("xl/sharedStrings.xml") as fh:
        for _, el in iterparse(fh):
            if el.tag == _NS + "si":
                out.append("".join(t.text or "" for t in el.iter(_NS + "t")))
                el.clear()
    return out


def _sheet_paths(z: zipfile.ZipFile) -> dict[str, str]:
    rels = {}
    with z.open("xl/_rels/workbook.xml.rels") as fh:
        for _, el in iterparse(fh):
            if el.tag.endswith("Relationship"):
                rels[el.get("Id")] = el.get("Target")
    out = {}
    with z.open("xl/workbook.xml") as fh:
        for _, el in iterparse(fh):
            if el.tag == _NS + "sheet":
                target = rels[el.get(_REL + "id")].lstrip("/")
                out[el.get("name")] = target if target.startswith("xl/") else "xl/" + target
    return out


def sheet_names(path: Path) -> list[str]:
    with zipfile.ZipFile(path) as z:
        return list(_sheet_paths(z))


def iter_rows(path: Path, sheet: str) -> Iterator[list[str]]:
    """Yield each row of `sheet` as a list of strings (gaps filled with "")."""
    with zipfile.ZipFile(path) as z:
        paths = _sheet_paths(z)
        if sheet not in paths:
            raise KeyError(f"{path}: no sheet {sheet!r} (sheets: {list(paths)})")
        ss = _shared_strings(z)
        with z.open(paths[sheet]) as fh:
            for _, el in iterparse(fh):
                if el.tag != _NS + "row":
                    continue
                row: list[str] = []
                for c in el.iter(_NS + "c"):
                    i = _col_index(c.get("r")) if c.get("r") else len(row)
                    while len(row) < i:
                        row.append("")
                    t = c.get("t")
                    v = c.find(_NS + "v")
                    if t == "s" and v is not None:
                        val = ss[int(v.text)]
                    elif t == "inlineStr":
                        val = "".join(x.text or "" for x in c.iter(_NS + "t"))
                    else:
                        val = v.text if v is not None and v.text is not None else ""
                    row.append(val)
                el.clear()
                yield row
