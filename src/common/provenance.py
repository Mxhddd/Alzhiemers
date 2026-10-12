"""Provenance: who/what/when produced a result.

Every run header records the git commit, whether the working tree had
uncommitted changes (and which files), the SHA-256 of the config file, the
SHA-256 of the lockfiles, the UTC timestamp, and the versions of the packages
the module names. Anything that cannot be determined raises; nothing is
recorded as "unknown".

Command-line use (also called from R so both languages hash identically):

    python -m src.common.provenance sha256 FILE
    python -m src.common.provenance git
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

from src.common.paths import LOCKFILES, REPO_ROOT


class ProvenanceError(RuntimeError):
    """Raised when provenance cannot be established. Never caught internally."""


@dataclass(frozen=True)
class GitState:
    commit: str
    dirty: bool
    uncommitted_files: tuple[str, ...]


def _git(args: list[str], repo_root: Path) -> str:
    proc = subprocess.run(
        ["git", "-C", str(repo_root), *args],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise ProvenanceError(
            f"git {' '.join(args)} failed in {repo_root}: {proc.stderr.strip()}"
        )
    return proc.stdout


def git_state(repo_root: Path = REPO_ROOT) -> GitState:
    """Current commit and working-tree status.

    Raises if the directory is not a git repository or has no commits: a
    result that cannot be tied to a commit has no provenance.
    """
    commit = _git(["rev-parse", "HEAD"], repo_root).strip()
    # --porcelain=v1 is stable across git versions; untracked files count as
    # dirty because untracked code can change a result. Ignored files
    # (data/raw, results) do not appear here by design.
    status = _git(["status", "--porcelain=v1", "--untracked-files=all"], repo_root)
    files = tuple(line[3:] for line in status.splitlines() if line.strip())
    return GitState(commit=commit, dirty=bool(files), uncommitted_files=files)


def sha256_file(path: Path) -> str:
    """SHA-256 of a file's raw bytes. Raises FileNotFoundError if absent."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"cannot hash {path}: not a file")
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def package_versions(names: list[str]) -> dict[str, str]:
    """Installed versions of the named distributions, plus Python itself.

    A package that is not installed raises PackageNotFoundError rather than
    being reported as missing — a run should not start in an environment it
    cannot describe.
    """
    versions = {"python": platform.python_version()}
    for name in names:
        versions[name] = metadata.version(name)
    return versions


def lockfile_hashes(lockfiles: tuple[Path, ...] = LOCKFILES) -> dict[str, str]:
    """SHA-256 of each lockfile. Both lockfiles must exist (ground rule 5)."""
    return {p.name: sha256_file(p) for p in lockfiles}


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def iso_utc(ts: datetime) -> str:
    if ts.tzinfo is None:
        raise ProvenanceError("refusing to format a timezone-naive timestamp")
    return ts.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def build_header(
    *,
    module: str,
    config_path: Path,
    seed: int | None,
    packages: list[str],
    started: datetime,
    repo_root: Path = REPO_ROOT,
    lockfiles: tuple[Path, ...] = LOCKFILES,
    preregistration: dict | None = None,
) -> dict:
    """Assemble the run header as an ordered dict (serialised as YAML)."""
    git = git_state(repo_root)
    config_path = Path(config_path).resolve()
    try:
        config_display = str(config_path.relative_to(repo_root))
    except ValueError:
        config_display = str(config_path)
    return {
        "lc_first_run_header": 1,
        "module": module,
        "language": "python",
        "started_utc": iso_utc(started),
        "git": {
            "commit": git.commit,
            "dirty": git.dirty,
            "uncommitted_files": list(git.uncommitted_files),
        },
        "config": {"path": config_display, "sha256": sha256_file(config_path)},
        "lockfiles": lockfile_hashes(lockfiles),
        "seed": seed,
        "packages": package_versions(packages),
        "platform": platform.platform(),
        "command": sys.argv,
        "preregistration": preregistration,
    }


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m src.common.provenance")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p_hash = sub.add_parser("sha256", help="print the SHA-256 of a file")
    p_hash.add_argument("file", type=Path)
    sub.add_parser("git", help="print git state as JSON")
    args = parser.parse_args(argv)
    if args.cmd == "sha256":
        print(sha256_file(args.file))
    elif args.cmd == "git":
        s = git_state()
        print(json.dumps({"commit": s.commit, "dirty": s.dirty,
                          "uncommitted_files": list(s.uncommitted_files)}))
    return 0


if __name__ == "__main__":
    sys.exit(_main())
