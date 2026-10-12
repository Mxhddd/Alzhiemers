"""Input guards.

Two ground rules are enforced here rather than trusted to habit:

* Simulated data (anything under ``sim/`` or ``tests/``) must never be read by
  an analysis pipeline. ``analysis_input`` refuses such paths.
* Data the researcher supplies by hand (e.g. ``data/raw/braak_counts.csv``)
  must never be created by code. When it is missing, ``analysis_input`` raises
  a ``MissingInputError`` that says who has to provide it — it never falls
  back to anything.
"""

from __future__ import annotations

from pathlib import Path

from src.common.paths import DATA_DIR, REPO_ROOT, SIM_DIR, TESTS_DIR


class SimulatedDataError(PermissionError):
    pass


class MissingInputError(FileNotFoundError):
    pass


def _is_under(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def analysis_input(path: Path, *, provided_by: str, data_dir: Path = DATA_DIR,
                   forbidden_dirs: tuple[Path, ...] = (SIM_DIR, TESTS_DIR)) -> Path:
    """Validate a path an analysis pipeline is about to read, and return it.

    ``provided_by`` says who produces the file (e.g. "researcher, hand-
    transcribed" or "download script src/mr/fetch.R"); it goes into the error
    if the file is missing so the message tells you what to do.
    """
    resolved = Path(path).resolve()
    for forbidden in forbidden_dirs:
        if _is_under(resolved, forbidden.resolve()):
            raise SimulatedDataError(
                f"{path} is under {forbidden.name}/ — simulated/test data must never "
                "be read by an analysis pipeline"
            )
    if not _is_under(resolved, Path(data_dir).resolve()):
        raise SimulatedDataError(
            f"{path} is outside data/ — analysis pipelines read only from data/"
        )
    if not resolved.is_file():
        try:
            shown = resolved.relative_to(REPO_ROOT)
        except ValueError:
            shown = resolved
        raise MissingInputError(
            f"required input {shown} does not exist. It is provided by: {provided_by}. "
            "Code must not create, scrape or fill in this file — ask for it."
        )
    return resolved
