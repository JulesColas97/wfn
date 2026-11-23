"""Shared helpers to read the single YAML input file.

The whole simulation is described in one YAML file, but every module only
reads its own top-level section (``receivers``, ``flow``, ``source``,
``propagation``...). This module therefore knows nothing about any section:
it only loads the file and converts plain dictionaries into typed dataclasses.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, TypeVar

import cattrs
import numpy as np
import yaml

T = TypeVar("T")

converter = cattrs.Converter(forbid_extra_keys=True)
converter.register_structure_hook(
    np.ndarray, lambda v, _: None if v is None else np.asarray(v, dtype=float)
)
converter.register_unstructure_hook(np.ndarray, lambda a: a.tolist())


def load_yaml(path: str | Path) -> dict[str, Any]:
    """Load a YAML file and return its content as a dictionary."""
    with open(path, "r") as f:
        data = yaml.safe_load(f)
    return data or {}


def get_section(path_or_dict: str | Path | dict, section: str, required: bool = True) -> dict | None:
    """Return one top-level section of the input file (or of an already loaded dict)."""
    data = path_or_dict if isinstance(path_or_dict, dict) else load_yaml(path_or_dict)
    if section not in data:
        if required:
            raise KeyError(f"Section '{section}' not found in input file")
        return None
    return data[section] or {}


def structure(data: dict, cls: type[T]) -> T:
    """Convert a plain dictionary into the dataclass ``cls`` (unknown keys raise an error)."""
    return converter.structure(data, cls)


def unstructure(obj: Any) -> Any:
    """Convert a dataclass into plain python types (YAML serialisable)."""
    return converter.unstructure(obj)


def resolve_path(path: str | Path | None, base_dir: str | Path | None) -> Path | None:
    """Resolve ``path`` relative to ``base_dir`` (typically the input file directory)."""
    if path is None:
        return None
    path = Path(path).expanduser()
    if path.is_absolute() or base_dir is None:
        return path
    return Path(base_dir) / path
