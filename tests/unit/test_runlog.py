import re

import pytest

from src.common.paths import CONFIG_DIR
from src.common.provenance import git_state, sha256_file
from src.common.runlog import read_footer, read_header, start_run

CFG = CONFIG_DIR / "provenance_check.yaml"


def test_header_records_provenance(tmp_path):
    with start_run(module="unit", config_path=CFG, seed=7,
                   packages=["pyyaml"], logs_dir=tmp_path) as run:
        run.info("hello")
    h = read_header(run.log_path)
    assert h["git"]["commit"] == git_state().commit
    assert re.fullmatch(r"[0-9a-f]{40}", h["git"]["commit"])
    assert isinstance(h["git"]["dirty"], bool)
    assert h["config"] == {"path": "config/provenance_check.yaml", "sha256": sha256_file(CFG)}
    assert set(h["lockfiles"]) == {"uv.lock", "renv.lock"}
    assert h["seed"] == 7
    assert set(h["packages"]) == {"python", "pyyaml"}
    assert re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{6}Z", h["started_utc"])
    assert run.log_path.parent == tmp_path / "unit"
    assert read_footer(run.log_path)["status"] == "completed"


def test_seed_is_mandatory(tmp_path):
    with pytest.raises(TypeError):
        start_run(module="unit", config_path=CFG, packages=[], logs_dir=tmp_path)  # type: ignore[call-arg]


def test_failure_is_logged_and_propagates(tmp_path):
    with pytest.raises(ZeroDivisionError):
        with start_run(module="unit", config_path=CFG, seed=None,
                       packages=[], logs_dir=tmp_path) as run:
            1 / 0
    footer = read_footer(run.log_path)
    assert footer["status"] == "FAILED"
    assert "ZeroDivisionError" in footer["error"]


def test_drops_are_counted_and_reported(tmp_path):
    report = tmp_path / "dropped.md"
    with start_run(module="unit", config_path=CFG, seed=None,
                   packages=[], logs_dir=tmp_path) as run:
        run.drop(stage="instruments", item="rs1", reason="F < 10")
        run.drop(stage="instruments", item="rs2", reason="F < 10")
        run.drop(stage="harmonise", item="rs3", reason="palindromic, MAF near 0.5")
        run.write_drop_report(report)
        with pytest.raises(ValueError):
            run.drop(stage="x", item="", reason="no item")
    footer = read_footer(run.log_path)
    assert footer["n_dropped"] == 3
    assert footer["dropped_by_stage_and_reason"] == {
        "instruments": {"F < 10": 2},
        "harmonise": {"palindromic, MAF near 0.5": 1},
    }
    assert "Total dropped: **3**" in report.read_text()
    assert str(report.resolve()) in footer["outputs_sha256"]


def test_missing_config_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        start_run(module="unit", config_path=tmp_path / "nope.yaml", seed=None,
                  packages=[], logs_dir=tmp_path)


def test_unverified_inputs_stamp_tables_and_footer(tmp_path):
    import pandas as pd
    with start_run(module="unit", config_path=CFG, seed=None, packages=[], logs_dir=tmp_path) as run:
        run.mark_unverified("x.csv: 3 of 3 values unverified")
        out = run.write_table(pd.DataFrame({"a": [1, 2]}), tmp_path / "t.csv")
    df = pd.read_csv(out)
    assert list(df.columns)[0] == "UNVERIFIED_INPUTS"
    assert read_footer(run.log_path)["unverified_inputs"] == ["x.csv: 3 of 3 values unverified"]


def test_verified_inputs_leave_tables_unstamped(tmp_path):
    import pandas as pd
    with start_run(module="unit", config_path=CFG, seed=None, packages=[], logs_dir=tmp_path) as run:
        out = run.write_table(pd.DataFrame({"a": [1]}), tmp_path / "t.csv")
    assert list(pd.read_csv(out).columns) == ["a"]
