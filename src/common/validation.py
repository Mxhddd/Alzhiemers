"""The validation gate (ground rule 1), enforced in code.

"A module that cannot recover a known answer must not be used on a new
question." ``make validate`` runs each module's reproductions of published
results and records the outcome in ``results/validation/<module>.status.yaml``
together with a **code fingerprint**: the SHA-256 over every file the result
depends on (shared code, the module's code, the validation harness and the
module's method config). A module calls ``require_validated(module)`` before
any analysis; it passes only if the last validation passed AND the code
fingerprint is unchanged since. Editing an estimator, a harmonisation rule or a
method setting therefore forces re-validation.

Command line (the R validation runner uses the first two):

    python -m src.common.validation fingerprint mr
    python -m src.common.validation record mr --result-file F --log L
    python -m src.common.validation check mr
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import yaml

from src.common.paths import REPO_ROOT, RESULTS_DIR
from src.common.provenance import git_state, iso_utc, sha256_file, utc_now

STATUS_DIR = RESULTS_DIR / "validation"

# What each module's validated result depends on. A new module is added here
# in the session that builds it.
FINGERPRINT_PATHS: dict[str, tuple[str, ...]] = {
    "mr": ("src/common", "src/mr", "tests/validation", "config/mr_methods.yaml"),
    "metaanalysis": ("src/common", "src/metaanalysis", "tests/validation/harness.R",
                     "tests/validation/run_validation.R", "tests/validation/targets",
                     "tests/validation/meta_published_pooled.R", "tests/validation/meta_hand_three_study.R",
                     "tests/validation/meta_sim_homogeneous.R", "config/metaanalysis.yaml"),
    "dynamics": ("src/common", "src/dynamics", "tests/validation/test_dynamics.py",
                 "tests/validation/test_dynamics_lc_model.py", "tests/validation/targets/dynamics_analytic_fixed_point.yaml",
                 "tests/validation/run_validation_py.py", "config/dynamics.yaml", "config/dynamics_params.yaml"),
    "mediation": ("src/common", "src/mediation", "sim/mediation_sim.py", "tests/validation/test_mediation.py",
                  "tests/validation/targets/mediation_published.yaml", "tests/validation/run_validation_py.py",
                  "config/mediation.yaml"),
    "baserate": ("src/common", "src/baserate", "sim/baserate_sim.py", "tests/validation/test_baserate.py",
                 "tests/validation/run_validation_py.py", "config/baserate.yaml"),
}

_EXCLUDE_PARTS = {"__pycache__", ".pytest_cache"}


class ValidationGateError(PermissionError):
    pass


def _files_for(module: str, repo_root: Path) -> list[Path]:
    if module not in FINGERPRINT_PATHS:
        raise ValidationGateError(f"no validation is defined for module '{module}'")
    files: list[Path] = []
    for rel in FINGERPRINT_PATHS[module]:
        p = repo_root / rel
        if p.is_file():
            files.append(p)
        elif p.is_dir():
            files.extend(f for f in p.rglob("*")
                         if f.is_file() and not (_EXCLUDE_PARTS & set(f.parts))
                         and f.suffix != ".pyc")
        else:
            raise FileNotFoundError(f"fingerprint path {rel} does not exist")
    return sorted(files)


def code_fingerprint(module: str, repo_root: Path = REPO_ROOT) -> str:
    h = hashlib.sha256()
    for f in _files_for(module, repo_root):
        h.update(f"{f.relative_to(repo_root)}\0{sha256_file(f)}\n".encode())
    return h.hexdigest()


def status_path(module: str, status_dir: Path = STATUS_DIR) -> Path:
    return Path(status_dir) / f"{module}.status.yaml"


def record(module: str, *, passed: bool, targets: list[dict], log: str,
           status_dir: Path = STATUS_DIR, repo_root: Path = REPO_ROOT) -> Path:
    if not targets:
        raise ValidationGateError(f"module '{module}' has no validation targets; "
                                  "an empty validation cannot pass")
    if passed and not all(t.get("passed") is True for t in targets):
        raise ValidationGateError("overall pass recorded while a target failed")
    git = git_state(repo_root)
    status = {
        "module": module,
        "passed": bool(passed),
        "targets": targets,
        "code_fingerprint": code_fingerprint(module, repo_root),
        "commit": git.commit,
        "dirty": git.dirty,
        "validated_at_utc": iso_utc(utc_now()),
        "log": log,
    }
    out = status_path(module, status_dir)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(status, sort_keys=False), encoding="utf-8")
    return out


def check(module: str, *, status_dir: Path = STATUS_DIR, repo_root: Path = REPO_ROOT) -> dict:
    p = status_path(module, status_dir)
    if not p.is_file():
        raise ValidationGateError(
            f"module '{module}' has never been validated ({p} missing). Run `make validate`.")
    status = yaml.safe_load(p.read_text(encoding="utf-8"))
    if status.get("passed") is not True:
        failed = [t.get("id") for t in status.get("targets", []) if t.get("passed") is not True]
        raise ValidationGateError(
            f"module '{module}' FAILED validation (targets: {failed}); it must not be used.")
    current = code_fingerprint(module, repo_root)
    if status.get("code_fingerprint") != current:
        raise ValidationGateError(
            f"module '{module}' code changed since it was validated "
            f"({status.get('validated_at_utc')}). Re-run `make validate`.")
    return status


def require_validated(module: str) -> dict:
    """Call at the top of every analysis entry point."""
    return check(module)


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m src.common.validation")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p_f = sub.add_parser("fingerprint"); p_f.add_argument("module")
    p_r = sub.add_parser("record"); p_r.add_argument("module")
    p_r.add_argument("--result-file", type=Path, required=True)
    p_r.add_argument("--log", required=True)
    p_c = sub.add_parser("check"); p_c.add_argument("module")
    args = parser.parse_args(argv)
    if args.cmd == "fingerprint":
        print(code_fingerprint(args.module))
    elif args.cmd == "record":
        res = yaml.safe_load(args.result_file.read_text(encoding="utf-8"))
        out = record(args.module, passed=res["passed"], targets=res["targets"], log=args.log)
        print(f"validation status for '{args.module}': "
              f"{'PASSED' if res['passed'] else 'FAILED'} -> {out}")
    else:
        s = check(args.module)
        print(f"module '{args.module}' validated {s['validated_at_utc']} (commit {s['commit'][:7]})")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
