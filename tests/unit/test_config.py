import pytest

from src.common.config import ConfigError, load_config, uncommented_keys
from src.common.paths import CONFIG_DIR


def _write(tmp_path, text):
    p = tmp_path / "c.yaml"
    p.write_text(text)
    return p


def test_get_nested_and_missing_key_raises(tmp_path):
    cfg = load_config(_write(tmp_path, "a:\n  b: 3  # why\n"))
    assert cfg.get("a.b") == 3
    with pytest.raises(ConfigError, match="a.c"):
        cfg.get("a.c")


def test_empty_value_raises(tmp_path):
    cfg = load_config(_write(tmp_path, "a:  # why\n"))
    with pytest.raises(ConfigError, match="empty"):
        cfg.get("a")


def test_duplicate_keys_raise(tmp_path):
    with pytest.raises(ConfigError, match="duplicate key 'f_min'"):
        load_config(_write(tmp_path, "f_min: 10  # why\nf_min: 5  # oops\n"))


def test_uncommented_detection(tmp_path):
    p = _write(tmp_path, "# why a\na: 1\nb: 2\nc: 3  # why c\nblock:\n  # why d\n  d: 4\n  e: 5\n")
    assert [line for line, _ in uncommented_keys(p)] == [3, 8]


@pytest.mark.parametrize("path", sorted(CONFIG_DIR.glob("*.yaml")), ids=lambda p: p.name)
def test_every_config_value_has_a_reason(path):
    """Ground rule 3: every setting in config/ carries a comment saying why."""
    assert uncommented_keys(path) == []
    load_config(path)  # and parses strictly
