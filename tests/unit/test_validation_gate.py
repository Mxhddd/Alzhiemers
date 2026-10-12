import subprocess

import pytest

from src.common.validation import ValidationGateError, check, code_fingerprint, record


def _git(repo, *args):
    subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t", *args],
                   check=True, capture_output=True)


@pytest.fixture
def repo(tmp_path, monkeypatch):
    for d in ("src/common", "src/mr", "tests/validation", "config"):
        (tmp_path / d).mkdir(parents=True)
    (tmp_path / "src/common/a.py").write_text("x = 1\n")
    (tmp_path / "src/mr/est.R").write_text("f <- 1\n")
    (tmp_path / "tests/validation/t.yaml").write_text("id: t\n")
    (tmp_path / "config/mr_methods.yaml").write_text("seed: 1  # why\n")
    _git(tmp_path, "init", "-q")
    _git(tmp_path, "add", ".")
    _git(tmp_path, "commit", "-qm", "init")
    return tmp_path


def _kw(repo):
    return {"status_dir": repo / "status", "repo_root": repo}


PASS = [{"id": "t", "passed": True, "detail": "ok"}]


def test_never_validated_blocks(repo):
    with pytest.raises(ValidationGateError, match="never been validated"):
        check("mr", **_kw(repo))


def test_passed_and_unchanged_passes(repo):
    record("mr", passed=True, targets=PASS, log="l", **_kw(repo))
    assert check("mr", **_kw(repo))["passed"] is True


def test_failed_validation_blocks(repo):
    record("mr", passed=False, targets=[{"id": "t", "passed": False, "detail": "outside CI"}],
           log="l", **_kw(repo))
    with pytest.raises(ValidationGateError, match="FAILED validation"):
        check("mr", **_kw(repo))


@pytest.mark.parametrize("edit", ["src/mr/est.R", "src/common/a.py", "config/mr_methods.yaml",
                                  "tests/validation/t.yaml"])
def test_any_code_or_setting_change_forces_revalidation(repo, edit):
    record("mr", passed=True, targets=PASS, log="l", **_kw(repo))
    (repo / edit).write_text((repo / edit).read_text() + "# changed\n")
    with pytest.raises(ValidationGateError, match="changed since"):
        check("mr", **_kw(repo))


def test_new_file_changes_fingerprint(repo):
    before = code_fingerprint("mr", repo)
    (repo / "src/mr/extra.R").write_text("")
    assert code_fingerprint("mr", repo) != before


def test_empty_or_inconsistent_records_refused(repo):
    with pytest.raises(ValidationGateError, match="no validation targets"):
        record("mr", passed=True, targets=[], log="l", **_kw(repo))
    with pytest.raises(ValidationGateError, match="target failed"):
        record("mr", passed=True, targets=[{"id": "t", "passed": False}], log="l", **_kw(repo))


def test_unknown_module_refused(repo):
    with pytest.raises(ValidationGateError, match="no validation is defined"):
        code_fingerprint("astrology", repo)
