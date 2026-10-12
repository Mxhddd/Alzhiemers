"""Strict YAML config loading.

Two failure modes of ordinary YAML handling are closed here:

1. **Silent defaults.** ``Config.get`` has no default argument. A missing key
   raises ``ConfigError`` naming the key and the file.
2. **Silent overwrites.** Plain PyYAML keeps the *last* of two duplicate keys
   without warning, so a pasted block could quietly change a threshold. This
   loader raises on any duplicate key.

``uncommented_keys`` backs a unit test that every ``key: value`` line in
``config/*.yaml`` carries a comment saying why (ground rule 3).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import yaml


class ConfigError(KeyError):
    pass


class _StrictLoader(yaml.SafeLoader):
    pass


def _construct_mapping(loader: _StrictLoader, node: yaml.MappingNode, deep: bool = False):
    seen: dict[Any, yaml.Node] = {}
    for key_node, _ in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in seen:
            raise ConfigError(
                f"duplicate key {key!r} at line {key_node.start_mark.line + 1} "
                f"(first defined at line {seen[key].start_mark.line + 1})"
            )
        seen[key] = key_node
    return loader.construct_mapping(node, deep=deep)


_StrictLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping
)


class Config:
    """Read-only view of a YAML config with dotted-path access and no defaults."""

    def __init__(self, data: dict, path: Path):
        if not isinstance(data, dict):
            raise ConfigError(f"{path}: top level must be a mapping")
        self._data = data
        self.path = Path(path)

    def get(self, dotted: str) -> Any:
        node: Any = self._data
        walked: list[str] = []
        for part in dotted.split("."):
            walked.append(part)
            if not isinstance(node, dict) or part not in node:
                raise ConfigError(
                    f"{self.path}: required key '{'.'.join(walked)}' is missing "
                    "(there are no defaults — add it with a comment saying why)"
                )
            node = node[part]
        if node is None:
            raise ConfigError(
                f"{self.path}: key '{dotted}' is empty; give it an explicit value"
            )
        return node

    def as_dict(self) -> dict:
        return yaml.safe_load(yaml.safe_dump(self._data))  # deep copy


def load_config(path: Path) -> Config:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"config file not found: {path}")
    with path.open(encoding="utf-8") as fh:
        data = yaml.load(fh, Loader=_StrictLoader)  # noqa: S506 (SafeLoader subclass)
    return Config(data, path)


_KEY_VALUE = re.compile(r"^\s*(?:-\s+)?[A-Za-z0-9_.\-]+\s*:\s*(?P<value>.*)$")


def uncommented_keys(path: Path) -> list[tuple[int, str]]:
    """Lines of the form ``key: scalar`` with no comment on the line or directly above.

    Keys that open a nested block (``key:`` with nothing after it) are exempt;
    their children are checked instead. List items (``- value``) are covered
    by the comment on their parent key.
    """
    problems: list[tuple[int, str]] = []
    prev = ""
    for lineno, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        line = raw.rstrip()
        stripped = line.strip()
        if not stripped or stripped == "---":
            continue
        if stripped.startswith("#"):
            prev = stripped
            continue
        m = _KEY_VALUE.match(line)
        if m:
            value = m.group("value").strip()
            has_inline = "#" in value and not value.startswith(("'", '"'))
            opens_block = value == "" or value in ("|", ">")
            if not opens_block and not has_inline and not prev.startswith("#"):
                problems.append((lineno, stripped))
        prev = stripped
    return problems
