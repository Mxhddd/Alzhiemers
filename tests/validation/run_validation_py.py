"""Run a Python module's validation tests and record the verdict.

    python -m tests.validation.run_validation_py <module>

Runs tests/validation/test_<module>*.py under pytest inside a logged run,
reads the per-test outcomes from pytest's JUnit XML, and records them (with
the code fingerprint) via src.common.validation.record. Any failed, errored
or skipped test makes the module FAIL; zero tests is a failure too.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

from src.common.paths import CONFIG_DIR, REPO_ROOT
from src.common.runlog import start_run
from src.common.validation import record

CONFIGS = {"baserate": CONFIG_DIR / "baserate.yaml", "dynamics": CONFIG_DIR / "dynamics.yaml",
           "mediation": CONFIG_DIR / "mediation.yaml"}


def main(argv: list[str]) -> int:
    if len(argv) != 1 or argv[0] not in CONFIGS:
        raise SystemExit(f"usage: python -m tests.validation.run_validation_py {{{','.join(CONFIGS)}}}")
    module = argv[0]
    test_files = sorted((REPO_ROOT / "tests" / "validation").glob(f"test_{module}*.py"))
    if not test_files:
        raise RuntimeError(f"no tests/validation/test_{module}*.py files")
    with start_run(module=f"validation_{module}", config_path=CONFIGS[module], seed=None,
                   packages=["pytest", "numpy", "scipy"]) as run:
        with tempfile.TemporaryDirectory() as tmp:
            xml = Path(tmp) / "junit.xml"
            proc = subprocess.run([sys.executable, "-m", "pytest", "-m", "validation", *map(str, test_files),
                                   f"--junitxml={xml}", "-q", "-p", "no:cacheprovider"],
                                  cwd=REPO_ROOT, capture_output=True, text=True)
            run.info(proc.stdout[-4000:])
            if not xml.is_file():
                raise RuntimeError(f"pytest produced no report:\n{proc.stderr[-2000:]}")
            targets = []
            for case in ET.parse(xml).getroot().iter("testcase"):
                bad = [c.tag for c in case if c.tag in ("failure", "error", "skipped")]
                detail = "passed" if not bad else f"{bad[0]}: " + (case.find(bad[0]).get("message") or "")[:500]
                targets.append({"id": case.get("name"), "passed": not bad, "detail": detail})
        passed = bool(targets) and all(t["passed"] for t in targets) and proc.returncode == 0
        out = record(module, passed=passed, targets=targets,
                     log=str(run.log_path.relative_to(REPO_ROOT)))
        for t in targets:
            (run.info if t["passed"] else run.warning)(f"{'PASS' if t['passed'] else 'FAIL'} {t['id']}: {t['detail']}")
        run.info(f"validation status -> {out}")
        if not passed:
            raise RuntimeError(f"validation of module '{module}' FAILED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
