import pytest

from src.dynamics.params import Parameter, ParameterSet


def P(**kw):
    base = dict(name="k", value=1.0, units="1/day", source="assumption", confidence="moderate", range=None)
    base.update(kw)
    return Parameter(**base)


def test_provenance_is_mandatory():
    with pytest.raises(ValueError, match="source"):
        P(source="")
    with pytest.raises(ValueError, match="confidence"):
        P(confidence="high")


def test_weak_parameters_need_a_range_and_are_never_reported():
    with pytest.raises(ValueError, match="weak parameter needs a plausible range"):
        P(confidence="weak")
    ps = ParameterSet.from_list([P(name="a"), P(name="b", confidence="weak", range=(0.5, 2.0))])
    assert ps.reportable_values() == {"a": 1.0}
    assert ps.weak() == ["b"]
    assert ps.ranges(["b"]) == {"b": (0.5, 2.0)}
    with pytest.raises(ValueError, match="no plausible range"):
        ps.ranges(["a"])


def test_value_must_lie_in_its_range_and_overrides_must_exist():
    with pytest.raises(ValueError, match="outside"):
        P(range=(2.0, 3.0))
    ps = ParameterSet.from_list([P(name="a")])
    assert ps.values({"a": 5.0}) == {"a": 5.0}
    with pytest.raises(KeyError):
        ps.values({"zz": 1.0})


def test_log_scale_needs_a_positive_range():
    with pytest.raises(ValueError, match="scale"):
        P(confidence="weak", range=(0.5, 2.0), scale="cubic")
    with pytest.raises(ValueError, match="positive range"):
        P(value=0.0, confidence="weak", range=(0.0, 2.0), scale="log")


def test_shipped_parameters_have_provenance_and_published_values_have_an_extracted_source():
    """Every shipped parameter states units, source and confidence; anything not
    an assumption must point at an extraction file under data/extracted/ (rule:
    published values only with a source file, verified: false until the
    researcher checks it). At present every parameter is an assumption."""
    import re

    from src.common.paths import CONFIG_DIR, REPO_ROOT
    from src.dynamics.lc_model import PARAMETERS

    ps = ParameterSet.from_yaml(CONFIG_DIR / "dynamics_params.yaml")
    assert sorted(ps.params) == sorted(PARAMETERS)
    for name, p in ps.params.items():
        if p.source.startswith("assumption:"):
            assert p.confidence == "weak" and p.range is not None, f"{name}: an assumption must be weak and swept"
        else:
            m = re.search(r"data/extracted/[\w.\-/]+", p.source)
            assert m and (REPO_ROOT / m.group(0)).exists(), f"{name}: a published value needs its data/extracted/ file"
    assert ps.reportable_values() == {}, "no parameter may be reported as a point value"
