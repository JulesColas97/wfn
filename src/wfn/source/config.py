"""Content of the ``source`` section of the input file."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np

from ..config import get_section, load_yaml, structure

R_AIR = 287.0


@dataclass
class AtmosphereConfig:
    """Constant atmospheric properties used by the source model.

    ``c0`` is computed from ``T_ref`` if not given.
    """

    T_ref: float = 10 + 273.15
    gamma: float = 1.4
    rho0: float = 1.225
    nu0: float = 1.45e-5
    c0: Optional[float] = None

    @property
    def sound_speed(self) -> float:
        if self.c0 is not None:
            return float(self.c0)
        return float(np.sqrt(self.gamma * R_AIR * self.T_ref))


@dataclass
class SourceConfig:
    """``source`` section.

    Example::

        source:
          turbine: NREL_5MW          # default name, path to a YAML file, {file: ...} or inline dict
          airfoil_database: NREL5MW_14seg.nc
          frequencies: [50, 100, 200]
          n_beta: 12                 # optional, overrides the turbine value
          positions: [[0, 0, 90]]    # optional hub positions, else from the flow, else (0, 0, hub_height)
          atmosphere: {T_ref: 283.15, rho0: 1.225}
          n_jobs: 1
          output: results/source     # optional prefix: <output>_<i>.nc

    Relative paths are relative to the directory of the input file.
    """

    turbine: Any
    frequencies: list[float]
    airfoil_database: Optional[str] = None
    n_beta: Optional[int] = None
    positions: Optional[list[list[float]]] = None
    atmosphere: AtmosphereConfig = field(default_factory=AtmosphereConfig)
    n_jobs: int = 1
    output: Optional[str] = None
    base_dir: Optional[str] = None

    @classmethod
    def from_dict(cls, data: dict, base_dir: str | Path | None = None) -> "SourceConfig":
        cfg = structure(dict(data), cls)
        if base_dir is not None:
            cfg.base_dir = str(base_dir)
        return cfg

    @classmethod
    def from_yaml(cls, path: str | Path) -> "SourceConfig":
        return cls.from_dict(get_section(load_yaml(path), "source"), base_dir=Path(path).parent)

    def path(self, p: str | None) -> Path | None:
        if p is None:
            return None
        p = Path(p).expanduser()
        if p.is_absolute() or self.base_dir is None:
            return p
        return Path(self.base_dir) / p
