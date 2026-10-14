import pytest

from src.common.extracted import load_extracted
from src.common.paths import REPO_ROOT


class FakeRun:
    def __init__(self):
        self.marks = []

    def mark_unverified(self, d):
        self.marks.append(d)

    def warning(self, msg):
        self.marks.append(msg)


def test_shipped_extracted_files_carry_sources_and_unverified_flags():
    for name in ("dementia_incidence_cfas.csv", "amyloid_prevalence_jansen2015.csv",
                 "amyloid_hr_ossenkoppele2022.csv", "braak2011_abstract_totals.csv"):
        run = FakeRun()
        df = load_extracted(REPO_ROOT / "data/extracted" / name, provided_by="test", run=run)
        assert len(df) > 0 and run.marks, name     # all currently unverified -> run marked


def test_empty_counts_file_is_refused_not_imputed(tmp_path, monkeypatch):
    # (until 2026-10-14 the real Braak counts file was itself empty and served as this test's input)
    import src.common.extracted as ex
    monkeypatch.setattr(ex, "analysis_input", lambda p, provided_by: p)
    p = tmp_path / "braak2011_stage_counts.csv"
    p.write_text("age_lower,age_upper,stage,n,source_id,locator,verified\n")
    with pytest.raises(ValueError, match="no rows"):
        load_extracted(p, provided_by="test")


def test_missing_provenance_columns_are_refused(tmp_path, monkeypatch):
    import src.common.extracted as ex
    monkeypatch.setattr(ex, "analysis_input", lambda p, provided_by: p)
    p = tmp_path / "x.csv"
    p.write_text("a,verified\n1,false\n")
    with pytest.raises(ValueError, match="source"):
        load_extracted(p, provided_by="t")
    p.write_text("a,source_id,locator,verified\n1,s,T1,maybe\n")
    with pytest.raises(ValueError, match="true or false"):
        load_extracted(p, provided_by="t")


def test_machine_rechecked_must_be_valid_and_discrepancies_are_logged(tmp_path, monkeypatch):
    import src.common.extracted as ex

    monkeypatch.setattr(ex, "analysis_input", lambda p, provided_by: p)
    f = tmp_path / "x.csv"
    f.write_text("value,source_id,locator,verified,machine_rechecked\n1,p,Table 1,false,discrepancy\n")
    run = FakeRun()
    ex.load_extracted(f, provided_by="t", run=run)
    assert any("discrepancy" in m for m in run.marks)
    f.write_text("value,source_id,locator,verified,machine_rechecked\n1,p,Table 1,false,maybe\n")
    with pytest.raises(ValueError, match="machine_rechecked"):
        ex.load_extracted(f, provided_by="t", run=run)
