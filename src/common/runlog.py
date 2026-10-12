"""Run logging with a provenance header and an explicit drop ledger.

Usage (every Python module entry point looks like this):

    with start_run(module="baserate", config_path=cfg, seed=cfg_seed,
                   packages=["numpy", "scipy"]) as run:
        run.info("fitting model A")
        run.drop(stage="load", item="row 17", reason="age_upper < age_lower")
        run.record_output(path_to_table)

The log file is ``results/logs/<module>/<UTC timestamp>_<commit7>.log``. It
opens with a YAML header between ``---`` lines (machine-readable; see
``read_header``) and closes with a footer giving the status, the drop counts
and the hashes of every output. If the body raises, the footer records
``status: FAILED`` with the exception and the exception propagates — the
logger never swallows an error.
"""

from __future__ import annotations

import logging
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from types import TracebackType

import yaml

from src.common.paths import LOCKFILES, LOGS_DIR, REPO_ROOT
from src.common.provenance import build_header, iso_utc, sha256_file, utc_now

HEADER_DELIM = "---"


def _display(path: Path, repo_root: Path) -> str:
    """Repo-relative path when inside the repo, absolute otherwise."""
    path = Path(path).resolve()
    try:
        return str(path.relative_to(Path(repo_root).resolve()))
    except ValueError:
        return str(path)


@dataclass(frozen=True)
class Drop:
    stage: str
    item: str
    reason: str


@dataclass
class Run:
    module: str
    log_path: Path
    header: dict
    logger: logging.Logger
    repo_root: Path
    drops: list[Drop] = field(default_factory=list)
    outputs: dict[str, str] = field(default_factory=dict)
    _closed: bool = False

    # -- logging -------------------------------------------------------------
    def info(self, msg: str) -> None:
        self.logger.info(msg)

    def warning(self, msg: str) -> None:
        self.logger.warning(msg)

    # -- the drop ledger -----------------------------------------------------
    def drop(self, *, stage: str, item: str, reason: str) -> None:
        """Record one excluded record/instrument/study/target, with a reason.

        Every exclusion in every module must pass through here (ground rule:
        dropped records are counted and logged with reasons).
        """
        for name, value in (("stage", stage), ("item", item), ("reason", reason)):
            if not str(value).strip():
                raise ValueError(f"drop() needs a non-empty {name}")
        d = Drop(stage=stage, item=item, reason=reason)
        self.drops.append(d)
        self.logger.warning("DROPPED [%s] %s — %s", stage, item, reason)

    def drop_counts(self) -> dict[str, dict[str, int]]:
        counts: dict[str, Counter] = {}
        for d in self.drops:
            counts.setdefault(d.stage, Counter())[d.reason] += 1
        return {stage: dict(c) for stage, c in counts.items()}

    def write_drop_report(self, path: Path) -> Path:
        """Write every drop to a Markdown table (e.g. results/logs/mr_dropped.md)."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        lines = [
            f"# Dropped items — {self.module}",
            "",
            f"Run log: `{_display(self.log_path, self.repo_root)}`  ",
            f"Commit: `{self.header['git']['commit']}`  ",
            f"Total dropped: **{len(self.drops)}**",
            "",
        ]
        if self.drops:
            lines += ["| stage | item | reason |", "|---|---|---|"]
            lines += [
                f"| {d.stage} | {d.item} | {d.reason.replace('|', '/')} |"
                for d in self.drops
            ]
        else:
            lines.append("Nothing was dropped.")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        self.record_output(path)
        return path

    # -- outputs -------------------------------------------------------------
    def record_output(self, path: Path) -> None:
        """Log the SHA-256 of a file this run produced."""
        digest = sha256_file(path)
        key = _display(path, self.repo_root)
        self.outputs[key] = digest
        self.logger.info("OUTPUT %s sha256=%s", key, digest)

    # -- footer --------------------------------------------------------------
    def _close(self, status: str, error: str | None) -> None:
        if self._closed:
            raise RuntimeError("run already closed")
        self._closed = True
        footer = {
            "finished_utc": iso_utc(utc_now()),
            "status": status,
            "error": error,
            "n_dropped": len(self.drops),
            "dropped_by_stage_and_reason": self.drop_counts(),
            "outputs_sha256": self.outputs,
        }
        for h in list(self.logger.handlers):
            h.flush()
        with self.log_path.open("a", encoding="utf-8") as fh:
            fh.write(HEADER_DELIM + " footer\n")
            fh.write(yaml.safe_dump(footer, sort_keys=False, allow_unicode=True))
            fh.write(HEADER_DELIM + "\n")
        for h in list(self.logger.handlers):
            h.close()
            self.logger.removeHandler(h)

    def __enter__(self) -> "Run":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> bool:
        if exc is None:
            self._close("completed", None)
        else:
            self.logger.error("RUN FAILED: %s: %s", exc_type.__name__, exc)
            self._close("FAILED", f"{exc_type.__name__}: {exc}")
        return False  # never suppress the exception


def start_run(
    *,
    module: str,
    config_path: Path,
    seed: int | None,
    packages: list[str],
    preregistration: dict | None = None,
    logs_dir: Path = LOGS_DIR,
    repo_root: Path = REPO_ROOT,
    lockfiles: tuple[Path, ...] = LOCKFILES,
) -> Run:
    """Open a run log and write its provenance header.

    ``seed`` has no default on purpose: a module must say either which seed it
    used (taken from its config) or ``None`` to declare itself deterministic.
    """
    if not module or "/" in module:
        raise ValueError(f"bad module name: {module!r}")
    started = utc_now()
    header = build_header(
        module=module,
        config_path=config_path,
        seed=seed,
        packages=packages,
        started=started,
        repo_root=repo_root,
        lockfiles=lockfiles,
        preregistration=preregistration,
    )
    module_dir = Path(logs_dir) / module
    module_dir.mkdir(parents=True, exist_ok=True)
    stamp = started.strftime("%Y%m%dT%H%M%S%fZ")
    log_path = module_dir / f"{stamp}_{header['git']['commit'][:7]}.log"
    # "x" mode: two runs can never silently share (and overwrite) one log.
    with log_path.open("x", encoding="utf-8") as fh:
        fh.write(HEADER_DELIM + "\n")
        fh.write(yaml.safe_dump(header, sort_keys=False, allow_unicode=True))
        fh.write(HEADER_DELIM + "\n")

    logger = logging.getLogger(f"lcfirst.{module}.{stamp}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    fmt.converter = time.gmtime  # log lines in UTC, like the header
    file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    file_handler.setFormatter(fmt)
    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(fmt)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    run = Run(module=module, log_path=log_path, header=header, logger=logger,
              repo_root=Path(repo_root))
    run.info(f"log: {log_path}")
    if header["git"]["dirty"]:
        run.warning(
            "WORKING TREE IS DIRTY — this result is not reproducible from commit "
            f"{header['git']['commit'][:7]} alone. Uncommitted: "
            + ", ".join(header["git"]["uncommitted_files"])
        )
    return run


def read_header(log_path: Path) -> dict:
    """Parse the YAML header of a run log. Raises if the header is malformed."""
    lines = Path(log_path).read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != HEADER_DELIM:
        raise ValueError(f"{log_path}: no provenance header")
    try:
        end = lines.index(HEADER_DELIM, 1)
    except ValueError as e:
        raise ValueError(f"{log_path}: unterminated provenance header") from e
    return yaml.safe_load("\n".join(lines[1:end]))


def read_footer(log_path: Path) -> dict:
    """Parse the footer of a closed run log. Raises if the run never closed."""
    lines = Path(log_path).read_text(encoding="utf-8").splitlines()
    marker = HEADER_DELIM + " footer"
    if marker not in lines:
        raise ValueError(f"{log_path}: no footer (run did not close)")
    start = len(lines) - 1 - lines[::-1].index(marker)
    end = lines.index(HEADER_DELIM, start + 1)
    return yaml.safe_load("\n".join(lines[start + 1:end]))
