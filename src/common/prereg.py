"""The pre-registration gate (ground rule 7), enforced in code.

A pre-registration is two files in ``prereg/``:

* the document itself, e.g. ``prereg/mediation.md`` (the local copy), and
* a freeze record ``prereg/mediation.freeze.yaml`` written by ``freeze``::

      name: mediation
      document: mediation.md
      sha256: <hash of the document at freeze time>
      frozen_at_utc: 2026-11-02T10:15:00.000000Z
      osf_url: https://osf.io/abcde/

``check`` (and ``require_preregistration``, which modules call before reading
real target data) passes only if ALL of these hold:

1. the freeze record exists and has every field;
2. the document exists and its SHA-256 still equals the frozen hash — i.e. it
   has not been edited since freezing;
3. ``frozen_at_utc`` is a valid UTC timestamp that is not in the future;
4. ``osf_url`` points to osf.io — the external, third-party timestamp. A
   local git timestamp alone is under our own control and is not enough;
5. both files are committed to git with no uncommitted changes.

Command line:

    python -m src.common.prereg freeze NAME --document prereg/NAME.md --osf-url URL
    python -m src.common.prereg check NAME
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from src.common.paths import PREREG_DIR, REPO_ROOT
from src.common.provenance import iso_utc, sha256_file, utc_now

REQUIRED_FIELDS = ("name", "document", "sha256", "frozen_at_utc", "osf_url")
OSF_PREFIX = "https://osf.io/"
# A document containing this marker has an open issue that must be resolved
# (and the marker removed) before it can be frozen.
BLOCK_MARKER = "DO NOT FREEZE"


class PreregistrationError(PermissionError):
    pass


def freeze_path(name: str, prereg_dir: Path = PREREG_DIR) -> Path:
    return Path(prereg_dir) / f"{name}.freeze.yaml"


def freeze(name: str, document: Path, osf_url: str, *,
           prereg_dir: Path = PREREG_DIR) -> Path:
    """Write the freeze record. Refuses to overwrite an existing one."""
    document = Path(document).resolve()
    prereg_dir = Path(prereg_dir).resolve()
    if document.parent != prereg_dir:
        raise PreregistrationError(f"{document} must live directly in {prereg_dir}")
    if not osf_url.startswith(OSF_PREFIX):
        raise PreregistrationError(f"osf_url must start with {OSF_PREFIX}: {osf_url!r}")
    if BLOCK_MARKER in document.read_text(encoding="utf-8"):
        raise PreregistrationError(f"{document.name} is marked {BLOCK_MARKER!r}: resolve the open issue it "
                                   "describes and remove the marker first")
    out = freeze_path(name, prereg_dir)
    record = {
        "name": name,
        "document": document.name,
        "sha256": sha256_file(document),
        "frozen_at_utc": iso_utc(utc_now()),
        "osf_url": osf_url,
    }
    with out.open("x", encoding="utf-8") as fh:  # "x": never overwrite a freeze
        fh.write("# Written by `python -m src.common.prereg freeze`. Do not edit by hand.\n")
        fh.write(yaml.safe_dump(record, sort_keys=False))
    return out


def _git_clean_and_tracked(path: Path, repo_root: Path) -> str | None:
    rel = str(path.resolve().relative_to(Path(repo_root).resolve()))
    tracked = subprocess.run(
        ["git", "-C", str(repo_root), "ls-files", "--error-unmatch", "--", rel],
        capture_output=True, text=True, check=False,
    )
    if tracked.returncode != 0:
        return f"{rel} is not committed to git"
    status = subprocess.run(
        ["git", "-C", str(repo_root), "status", "--porcelain=v1", "--", rel],
        capture_output=True, text=True, check=False,
    )
    if status.returncode != 0:
        return f"git status failed for {rel}: {status.stderr.strip()}"
    if status.stdout.strip():
        return f"{rel} has uncommitted changes"
    return None


def check(name: str, *, prereg_dir: Path = PREREG_DIR, repo_root: Path = REPO_ROOT,
          now: datetime | None = None) -> dict:
    """Return the freeze record if the gate passes; raise listing every failure."""
    problems: list[str] = []
    fpath = freeze_path(name, prereg_dir)
    if not fpath.is_file():
        raise PreregistrationError(
            f"no frozen pre-registration for '{name}': {fpath} does not exist. "
            "Hypothesis-testing modules refuse to run on real target data until one does."
        )
    record = yaml.safe_load(fpath.read_text(encoding="utf-8")) or {}
    missing = [f for f in REQUIRED_FIELDS if not record.get(f)]
    if missing:
        raise PreregistrationError(f"{fpath}: missing fields {missing}")
    if record["name"] != name:
        problems.append(f"freeze record is for '{record['name']}', not '{name}'")

    doc = Path(prereg_dir) / record["document"]
    if not doc.is_file():
        problems.append(f"pre-registered document {doc} is missing")
    elif sha256_file(doc) != record["sha256"]:
        problems.append(f"{doc.name} was edited after freezing (SHA-256 mismatch)")

    try:
        frozen = datetime.strptime(record["frozen_at_utc"], "%Y-%m-%dT%H:%M:%S.%fZ")
        frozen = frozen.replace(tzinfo=timezone.utc)
        if frozen > (now or utc_now()):
            problems.append(f"frozen_at_utc {record['frozen_at_utc']} is in the future")
    except (TypeError, ValueError):
        problems.append(f"frozen_at_utc is not a valid UTC timestamp: {record['frozen_at_utc']!r}")

    if not str(record["osf_url"]).startswith(OSF_PREFIX):
        problems.append(f"osf_url must start with {OSF_PREFIX}")

    for p in (fpath, doc):
        if p.is_file():
            issue = _git_clean_and_tracked(p, repo_root)
            if issue:
                problems.append(issue)

    if problems:
        raise PreregistrationError(
            f"pre-registration gate for '{name}' FAILED:\n  - " + "\n  - ".join(problems)
        )
    return record


def require_preregistration(name: str) -> dict:
    """Call before any hypothesis test touches real target data.

    Returns the freeze record, which the module passes to ``start_run`` so the
    run header shows which pre-registration it ran under.
    """
    return check(name)


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m src.common.prereg")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p_f = sub.add_parser("freeze")
    p_f.add_argument("name")
    p_f.add_argument("--document", type=Path, required=True)
    p_f.add_argument("--osf-url", required=True)
    p_c = sub.add_parser("check")
    p_c.add_argument("name")
    args = parser.parse_args(argv)
    if args.cmd == "freeze":
        out = freeze(args.name, args.document, args.osf_url)
        print(f"wrote {out}. Commit it together with {args.document} now.")
    else:
        record = check(args.name)
        print(f"pre-registration '{args.name}' OK (frozen {record['frozen_at_utc']}, "
              f"{record['osf_url']})")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
