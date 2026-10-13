import subprocess
from datetime import datetime, timedelta, timezone

import pytest
import yaml

from src.common.prereg import PreregistrationError, check, freeze


def _git(repo, *args):
    subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t", *args],
                   check=True, capture_output=True)


@pytest.fixture
def repo(tmp_path):
    _git(tmp_path, "init", "-q")
    (tmp_path / "prereg").mkdir()
    doc = tmp_path / "prereg" / "mediation.md"
    doc.write_text("# DAG\nX -> M -> Y\n")
    return tmp_path


def _freeze_and_commit(repo):
    freeze("mediation", repo / "prereg" / "mediation.md", "https://osf.io/abcde/",
           prereg_dir=repo / "prereg")
    _git(repo, "add", "prereg")
    _git(repo, "commit", "-qm", "freeze")


def _check(repo, **kw):
    return check("mediation", prereg_dir=repo / "prereg", repo_root=repo, **kw)


def test_no_freeze_record_blocks(repo):
    with pytest.raises(PreregistrationError, match="no frozen pre-registration"):
        _check(repo)


def test_frozen_and_committed_passes(repo):
    _freeze_and_commit(repo)
    assert _check(repo)["osf_url"] == "https://osf.io/abcde/"


def test_uncommitted_freeze_blocks(repo):
    freeze("mediation", repo / "prereg" / "mediation.md", "https://osf.io/abcde/",
           prereg_dir=repo / "prereg")
    with pytest.raises(PreregistrationError, match="not committed"):
        _check(repo)


def test_editing_document_after_freeze_blocks(repo):
    _freeze_and_commit(repo)
    (repo / "prereg" / "mediation.md").write_text("# DAG\nX -> Y\n")
    with pytest.raises(PreregistrationError, match="edited after freezing"):
        _check(repo)
    _git(repo, "commit", "-qam", "sneaky edit")  # committing the edit does not help
    with pytest.raises(PreregistrationError, match="edited after freezing"):
        _check(repo)


def test_future_timestamp_blocks(repo):
    _freeze_and_commit(repo)
    past = datetime.now(timezone.utc) - timedelta(days=3650)
    with pytest.raises(PreregistrationError, match="in the future"):
        _check(repo, now=past)


def test_non_osf_url_refused(repo):
    with pytest.raises(PreregistrationError, match="osf"):
        freeze("mediation", repo / "prereg" / "mediation.md", "https://example.com/x",
               prereg_dir=repo / "prereg")
    fpath = repo / "prereg" / "mediation.freeze.yaml"
    _freeze_and_commit(repo)
    rec = yaml.safe_load(fpath.read_text())
    rec["osf_url"] = "https://example.com/x"
    fpath.write_text(yaml.safe_dump(rec))
    _git(repo, "commit", "-qam", "swap url")
    with pytest.raises(PreregistrationError, match="osf_url"):
        _check(repo)


def test_freeze_never_overwrites(repo):
    _freeze_and_commit(repo)
    with pytest.raises(FileExistsError):
        freeze("mediation", repo / "prereg" / "mediation.md", "https://osf.io/zzzzz/",
               prereg_dir=repo / "prereg")


def test_document_marked_do_not_freeze_is_refused(repo):
    doc = repo / "prereg" / "mediation.md"
    doc.write_text("# DAG\n> **DO NOT FREEZE** — open issue\n")
    with pytest.raises(PreregistrationError, match="DO NOT FREEZE"):
        freeze("mediation", doc, "https://osf.io/abcde/", prereg_dir=repo / "prereg")
    assert not (repo / "prereg" / "mediation.freeze.yaml").exists()
