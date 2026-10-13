"""What the researcher still has to supply — computed, never hand-maintained.

    python -m src.common.todo        (or `make todo`)

Lists every config / validation-target field that is null, every
hand-supplied data file referenced but absent, and every pre-registration
that is still a draft. Exit status is 0; this is a report, not a gate.
"""

from __future__ import annotations

import sys

import yaml

from src.common.paths import CONFIG_DIR, PREREG_DIR, REPO_ROOT
from src.common.prereg import BLOCK_MARKER

TARGETS_DIR = REPO_ROOT / "tests" / "validation" / "targets"
MODULES = ("mr", "baserate", "metaanalysis", "dynamics", "mediation", "network")


def null_fields(node, prefix="") -> list[str]:
    out = []
    if isinstance(node, dict):
        for k, v in node.items():
            out += null_fields(v, f"{prefix}.{k}" if prefix else str(k))
    elif node is None:
        out.append(prefix)
    return out


def data_paths(node) -> list[str]:
    out = []
    if isinstance(node, dict):
        for v in node.values():
            out += data_paths(v)
    elif isinstance(node, str) and node.startswith("data/raw/"):
        out.append(node)
    return out


def main() -> int:
    lines = []
    files = sorted(CONFIG_DIR.glob("*.yaml")) + sorted(TARGETS_DIR.glob("*.yaml"))
    missing_files = set()
    for f in files:
        doc = yaml.safe_load(f.read_text(encoding="utf-8")) or {}
        nulls = null_fields(doc)
        if nulls:
            lines.append(f"\n{f.relative_to(REPO_ROOT)}  ({len(nulls)} to fill)")
            lines += [f"  - {n}" for n in nulls]
        missing_files |= {p for p in data_paths(doc) if not (REPO_ROOT / p).exists()}
    if missing_files:
        lines.append(f"\nData files referenced but absent ({len(missing_files)}):")
        lines += [f"  - {p}" for p in sorted(missing_files)]
    extracted = sorted((REPO_ROOT / "data" / "extracted").glob("*.csv"))
    if extracted:
        lines.append("\nExtracted values (data/extracted/) — verify against the PDFs, flip verified to true:")
        for f in extracted:
            import csv
            with f.open(encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh))
            if not rows:
                lines.append(f"  - {f.name}: EMPTY (not yet extractable — see SOURCES.md)")
            else:
                n_unv = sum(1 for r in rows if str(r.get("verified", "")).strip().lower() != "true")
                lines.append(f"  - {f.name}: {n_unv} of {len(rows)} rows unverified")
    drafts = sorted((REPO_ROOT / "manuscript" / "derivations").glob("*.md"))
    if drafts:
        lines.append("\nDerivations to check by hand (manuscript/derivations/):")
        lines += [f"  - {d.name}" for d in drafts]
    params = yaml.safe_load((CONFIG_DIR / "dynamics_params.yaml").read_text(encoding="utf-8"))
    if not params.get("parameters"):
        lines.append("\nconfig/dynamics_params.yaml: no parameters yet; src/dynamics/lc_model.py: equations not written")
    lines.append("\nPre-registrations:")
    for m in MODULES:
        frozen = (PREREG_DIR / f"{m}.freeze.yaml").exists()
        complete = (PREREG_DIR / f"{m}.md").exists()
        draft = (PREREG_DIR / f"{m}.DRAFT.md").exists()
        blocked = complete and BLOCK_MARKER in (PREREG_DIR / f"{m}.md").read_text(encoding="utf-8")
        state = ("FROZEN" if frozen
                 else f"BLOCKED — open issue in prereg/{m}.md ({BLOCK_MARKER}); needs the researcher's decision"
                 if blocked
                 else "complete — review, upload to OSF, `make freeze-prereg`" if complete
                 else "partial draft (Claude to complete)" if draft else "not written")
        lines.append(f"  - {m}: {state}")
    print("Researcher inputs still needed (computed from the repository):" + "\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
