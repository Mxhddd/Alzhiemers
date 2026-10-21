"""Turn a pre-registration (prereg/<name>.md) into a PDF with a real text layer.

    python -m src.registration.pdf <name>        -> prereg/<name>.pdf

Markdown -> standalone HTML (pandoc) -> PDF (headless Chromium "print to PDF").
Chromium writes real text, not an image, so the PDF is searchable and its text
can be copied. After printing, the text layer is checked with pdftotext: the
title and a sample of sentences from the document must be found in it,
otherwise the PDF is deleted and the command fails.

This package is deliberately outside src/common: files there are part of every
module's code fingerprint.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PREREG = REPO_ROOT / "prereg"

CSS = """
@page { size: A4; margin: 20mm 18mm; }
body { font-family: 'DejaVu Serif', Georgia, serif; font-size: 10.5pt; line-height: 1.45; color: #111; }
h1 { font-size: 17pt; line-height: 1.25; margin: 0 0 4pt; }
h2 { font-size: 13pt; margin: 16pt 0 4pt; }
h3 { font-size: 11pt; margin: 12pt 0 3pt; }
p, li { orphans: 3; widows: 3; }
code { font-family: 'DejaVu Sans Mono', monospace; font-size: 8.8pt; background: #f2f2f2; padding: 0 2pt; }
table { border-collapse: collapse; width: 100%; font-size: 9pt; margin: 6pt 0; }
th, td { border-bottom: 1px solid #bbb; padding: 3pt 4pt; text-align: left; vertical-align: top; }
th { border-bottom: 1.5px solid #555; }
tr { break-inside: avoid; }
"""


def _chromium() -> str:
    for name in ("chromium", "chromium-browser", "google-chrome"):
        if shutil.which(name):
            return shutil.which(name)
    found = sorted(Path("/opt/pw-browsers").glob("chromium-*/chrome-linux*/chrome"))
    if found:
        return str(found[0])
    raise SystemExit("no Chromium found")


def _plain(s: str) -> str:
    s = re.sub(r"`([^`]*)`", r"\1", s)
    s = re.sub(r"\*\*|__|\*|_(?=\w)|(?<=\w)_", "", s)
    s = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", s)
    return re.sub(r"\s+", " ", s).strip()


def _norm(s: str) -> str:
    return re.sub(r"[^0-9a-z]+", "", s.lower())


def text_layer_ok(md: Path, pdf: Path, n_samples: int = 12) -> tuple[bool, list[str]]:
    """Title plus a spread of sentences from the document must be in the PDF text."""
    out = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], capture_output=True, text=True, check=True).stdout
    hay = _norm(out)
    lines = [ln for ln in md.read_text(encoding="utf-8").splitlines()]
    title = next((ln[2:] for ln in lines if ln.startswith("# ")), "")
    prose = [_plain(ln.lstrip("*-0123456789. ")) for ln in lines
             if ln.strip() and not ln.startswith(("#", "|", "```")) and len(ln) > 60]
    step = max(1, len(prose) // n_samples)
    samples = [title] + [p[:60] for p in prose[::step]][:n_samples]
    missing = [s for s in samples if _norm(s) and _norm(s) not in hay]
    return not missing and len(hay) > 500, missing


def build(name: str) -> Path:
    md = PREREG / f"{name}.md"
    pdf = PREREG / f"{name}.pdf"
    with tempfile.TemporaryDirectory() as td:
        css = Path(td) / "style.css"
        css.write_text(CSS, encoding="utf-8")
        html = Path(td) / f"{name}.html"
        title = next((ln[2:] for ln in md.read_text(encoding="utf-8").splitlines() if ln.startswith("# ")), name)
        subprocess.run(["pandoc", str(md), "-f", "gfm", "-t", "html5", "-s", "--embed-resources",
                        "--css", str(css), "--metadata", f"pagetitle={title}", "-o", str(html)], check=True)
        subprocess.run([_chromium(), "--headless", "--no-sandbox", "--disable-gpu", "--no-pdf-header-footer",
                        f"--print-to-pdf={pdf}", html.as_uri()], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=180)
    ok, missing = text_layer_ok(md, pdf)
    if not ok:
        pdf.unlink(missing_ok=True)
        raise SystemExit(f"{pdf.name}: text layer check failed; not found in the PDF text: {missing}")
    return pdf


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python -m src.registration.pdf <name>")
    p = build(sys.argv[1])
    print(f"wrote {p.relative_to(REPO_ROOT)} (text layer checked)")
