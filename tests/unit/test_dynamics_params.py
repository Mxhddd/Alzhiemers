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


def test_shipped_parameter_file_is_empty_until_the_researcher_writes_it():
    from src.common.paths import CONFIG_DIR
    assert ParameterSet.from_yaml(CONFIG_DIR / "dynamics_params.yaml").params == {}
