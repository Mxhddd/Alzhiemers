import pytest

from src.common.extracted import load_extracted
from src.common.paths import REPO_ROOT


class FakeRun:
    def __init__(self):
        self.marks = []

    def mark_unverified(self, d):
        self.marks.append(d)


def test_shipped_extracted_files_carry_sources_and_unverified_flags():
    for name in ("dementia_incidence_cfas.csv", "amyloid_prevalence_jansen2015.csv",
                 "amyloid_hr_ossenkoppele2022.csv", "braak2011_abstract_totals.csv"):
        run = FakeRun()
        df = load_extracted(REPO_ROOT / "data/extracted" / name, provided_by="test", run=run)
        assert len(df) > 0 and run.marks, name     # all currently unverified -> run marked


def test_empty_braak_counts_file_is_refused_not_imputed():
    with pytest.raises(ValueError, match="no rows"):
        load_extracted(REPO_ROOT / "data/extracted/braak2011_stage_counts.csv", provided_by="test")


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
