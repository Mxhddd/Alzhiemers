import hashlib
import subprocess
from importlib.metadata import PackageNotFoundError

import pytest

from src.common.paths import REPO_ROOT
from src.common.provenance import (
    ProvenanceError,
    git_state,
    package_versions,
    sha256_file,
)


def test_sha256_matches_hashlib(tmp_path):
    f = tmp_path / "x.yaml"
    f.write_bytes(b"alpha: 1  # why\n")
    assert sha256_file(f) == hashlib.sha256(b"alpha: 1  # why\n").hexdigest()


def test_sha256_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        sha256_file(tmp_path / "absent.yaml")


def test_git_state_matches_git():
    head = subprocess.run(["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()
    assert git_state().commit == head


def test_git_state_detects_dirty_tree(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    (tmp_path / "a.txt").write_text("1")
    subprocess.run(["git", "-C", str(tmp_path), "add", "a.txt"], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "-c", "user.name=t", "-c", "user.email=t@t",
                    "commit", "-qm", "init"], check=True)
    assert git_state(tmp_path).dirty is False
    (tmp_path / "a.txt").write_text("2")
    (tmp_path / "new.py").write_text("")
    s = git_state(tmp_path)
    assert s.dirty is True
    assert set(s.uncommitted_files) == {"a.txt", "new.py"}


def test_git_state_outside_repo_raises(tmp_path):
    with pytest.raises(ProvenanceError):
        git_state(tmp_path)


def test_package_versions_raise_for_missing_package():
    assert "python" in package_versions([])
    with pytest.raises(PackageNotFoundError):
        package_versions(["definitely-not-installed-pkg"])
