"""Read and write the Fortran namelist ``input.nml`` of the WAPE code."""
from __future__ import annotations

import dataclasses
import re
from pathlib import Path

import numpy as np

from .config import WapeInputs, backend_fields


def _format(key: str, value) -> list[str]:
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, str):
        return [f"{key}='{value}'"]
    if isinstance(value, (bool, np.bool_)):
        return [f"{key}={'.true.' if value else '.false.'}"]
    if isinstance(value, (list, tuple)):
        return [f"{key}(1)={str(list(value))[1:-1]}"]
    return [f"{key}={value}"]


def namelist_text(inputs: WapeInputs, backend: str = "piecewape") -> str:
    """Namelist content (``$inputs ... $end inputs``). ``None`` values (and empty lists) are skipped.

    Raises:
        ValueError: if a field set in ``inputs`` is not a key of ``backend``.
    """
    allowed = backend_fields(backend)
    lines = ["$inputs"]
    for f in dataclasses.fields(inputs):
        value = getattr(inputs, f.name)
        if value is None or (isinstance(value, (list, tuple)) and len(value) == 0):
            continue
        if f.name not in allowed:
            raise ValueError(f"{f.name!r} is not an input of the {backend} backend")
        lines += _format(f.name, value)
    lines.append("$end inputs")
    return "\n".join(lines) + "\n"


def write_namelist(inputs: WapeInputs, path: str | Path, backend: str = "piecewape") -> Path:
    """Write ``inputs`` to ``path`` (a file, or a directory where ``input.nml`` is created)."""
    path = Path(path)
    if path.is_dir():
        path = path / "input.nml"
    path.write_text(namelist_text(inputs, backend))
    return path


def _parse_value(val: str):
    val = val.strip().strip("'").strip('"')
    if val.lower() == ".true.":
        return True
    if val.lower() == ".false.":
        return False
    try:
        if "." in val or "e" in val.lower():
            return float(val)
        return int(val)
    except ValueError:
        return val


def read_namelist(path: str | Path) -> WapeInputs:
    """Build :class:`WapeInputs` from an ``input.nml`` file. Unknown keys are ignored."""
    path = Path(path)
    if path.is_dir():
        path = path / "input.nml"
    content = path.read_text()
    fields = {f.name: f for f in dataclasses.fields(WapeInputs)}
    data = {}
    for key, value in re.findall(r"(\w+)(?:\(1\))?\s*=\s*([^!\n$]+)", content):
        if key not in fields:
            continue
        parts = [_parse_value(p) for p in value.split(",") if p.strip()]
        if "list" in str(fields[key].type):
            data[key] = parts
        else:
            data[key] = parts[0] if len(parts) == 1 else parts
    return WapeInputs(**data)
