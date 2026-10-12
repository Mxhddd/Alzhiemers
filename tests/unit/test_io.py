import pytest

from src.common.io import MissingInputError, SimulatedDataError, analysis_input


@pytest.fixture
def layout(tmp_path):
    for d in ("data/raw", "sim", "tests"):
        (tmp_path / d).mkdir(parents=True)
    return tmp_path


def _call(layout, rel):
    return analysis_input(layout / rel, provided_by="the researcher",
                          data_dir=layout / "data",
                          forbidden_dirs=(layout / "sim", layout / "tests"))


def test_simulated_paths_are_refused(layout):
    (layout / "sim" / "fake.csv").write_text("x\n1\n")
    (layout / "tests" / "fake.csv").write_text("x\n1\n")
    with pytest.raises(SimulatedDataError):
        _call(layout, "sim/fake.csv")
    with pytest.raises(SimulatedDataError):
        _call(layout, "tests/fake.csv")


def test_paths_outside_data_are_refused(layout):
    (layout / "elsewhere.csv").write_text("x\n")
    with pytest.raises(SimulatedDataError):
        _call(layout, "elsewhere.csv")


def test_missing_input_names_who_provides_it(layout):
    with pytest.raises(MissingInputError, match="the researcher"):
        _call(layout, "data/raw/braak_counts.csv")
    assert not (layout / "data/raw/braak_counts.csv").exists()


def test_real_input_passes(layout):
    (layout / "data/raw/x.csv").write_text("a\n")
    assert _call(layout, "data/raw/x.csv").name == "x.csv"
