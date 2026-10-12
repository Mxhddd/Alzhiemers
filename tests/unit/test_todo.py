from src.common.todo import main, null_fields


def test_null_fields_walks_nested_mappings():
    assert null_fields({"a": None, "b": {"c": None, "d": 1}}) == ["a", "b.c"]


def test_report_runs(capsys):
    assert main() == 0
    assert "Pre-registrations:" in capsys.readouterr().out
