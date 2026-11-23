"""Read and write the Fortran namelist ``input.nml`` of the WAPE code."""
from __future__ import annotations

import dataclasses
import re
from pathlib import Path

import numpy as np

from .config import WapeInputs, backend_fields


_DEFAULTS = dataclasses.asdict(WapeInputs())


# keys of the old WAPE namelists, renamed in the current codes
LEGACY_KEYS = {"dout": "dx_out", "c0": "c_air", "rho0": "rho_air", "u0": "u_ref"}


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
    """Namelist content (``$inputs ... $end inputs``) for ``backend``. Empty lists are skipped.

    Every key of the backend is written (with its default value if not changed). A key
    of the other backend is ignored if it has its default value.

    Raises:
        ValueError: if a key that does not exist in ``backend`` was changed from its default.
    """
    allowed = backend_fields(backend)
    lines = ["$inputs"]
    for f in dataclasses.fields(inputs):
        value = getattr(inputs, f.name)
        if f.name not in allowed:
            if value != _DEFAULTS[f.name]:
                raise ValueError(f"{f.name!r} is not an input of the {backend} backend")
            continue
        if value is None or (isinstance(value, (list, tuple)) and len(value) == 0):
            continue
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
        key = LEGACY_KEYS.get(key, key)  # namelists written by the old code
        if key not in fields:
            continue
        parts = [_parse_value(p) for p in value.split(",") if p.strip()]
        if "list" in str(fields[key].type):
            data[key] = parts
        else:
            data[key] = parts[0] if len(parts) == 1 else parts
    if "cfl" in re.findall(r"^\s*(\w+)\s*=", content, flags=re.M):  # old key: dx = cfl * lambda
        data["dx_frac"] = 1 / float(re.search(r"^\s*cfl\s*=\s*([^!\n,]+)", content, flags=re.M).group(1))
    return WapeInputs(**data)
