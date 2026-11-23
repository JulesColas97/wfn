"""Content of the ``propagation`` section of the input file."""
from __future__ import annotations

import dataclasses
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from ..config import get_section, load_yaml, structure


BACKENDS = ("piecewape", "ewape")


@dataclass
class WapeInputs:
    """Parameters of the Fortran WAPE codes (``piecewape`` and ``EWAPE``).

    The field names are the keys of the Fortran namelist ``input.nml``: do not
    rename them. Fields left to ``None`` are not written to the namelist. Keys
    shared by both codes have the same name; the others are listed in
    :data:`EWAPE_ONLY` and :data:`PIECEWAPE_ONLY` and are refused by the other
    backend (see :func:`backend_fields`).
    """

    # ---- shared by both codes
    case_name: Optional[str] = "c0"
    c_air: Optional[float] = 343
    rho_air: Optional[float] = 1.2
    T0: Optional[float] = None  # temperature [K]
    p0: Optional[float] = None  # pressure [Pa]
    frequencies: Optional[list[float]] = field(default_factory=list)
    tau: Optional[list[float]] = None
    taumin: Optional[float] = None
    taumax: Optional[float] = None
    nb_tau: Optional[int] = None
    Lx1: Optional[float] = None
    Lx2: Optional[float] = None
    Ly1: Optional[float] = None
    Ly2: Optional[float] = None
    Lz: Optional[float] = None
    source_height: Optional[float] = None
    dx_frac: float = 10.0  # points per wavelength: dx = lambda / dx_frac
    pml_size: int = 30
    pml_param: float = 5e4
    pml_n: float = 2.5
    # flow
    external_flow: Optional[bool] = None
    uniform: Optional[bool] = None
    logarithmic: Optional[bool] = None
    u_ref: Optional[float] = None
    fname_flow: Optional[str] = None
    Xs_flow: Optional[float] = None
    Ys_flow: Optional[float] = None
    # output
    receiver_heights: list[float] = field(default_factory=list)
    side: bool = True

    # ---- piecewape only
    gamma0: Optional[float] = None
    pe_model: Optional[int] = None
    multi_domain: Optional[bool] = None
    domain_step: Optional[float] = None
    impedance_sigmae: Optional[float] = None
    impedance_alphae: Optional[float] = None
    impedance_mode: Optional[int] = None
    impedance_Rs: Optional[float] = None
    impedance_Om: Optional[float] = None
    impedance_q2: Optional[float] = None
    impedance_e: Optional[float] = None
    external_terrain: Optional[bool] = None
    fname_terrain: Optional[str] = None
    Xs_terrain: Optional[float] = None
    Ys_terrain: Optional[float] = None
    gaussian_hill: Optional[bool] = None
    hill_length: Optional[float] = None
    hill_pos: Optional[float] = None
    hill_height: Optional[float] = None
    dout: Optional[int] = None
    continuation: Optional[bool] = None
    write_flow: Optional[bool] = None
    xplane: Optional[list[float]] = None
    yplane: Optional[list[float]] = None

    # ---- EWAPE only
    H: Optional[float] = None
    c_water: Optional[float] = None
    rho_water: Optional[float] = None
    eta_air: Optional[float] = None
    eta_water: Optional[float] = None
    seabed_type: Optional[int] = None
    rho_sb: Optional[float] = None
    cp_sb: Optional[float] = None
    cs_sb: Optional[float] = None
    alpha_p: Optional[float] = None
    alpha_s: Optional[float] = None
    theta_sb: Optional[float] = None
    d_pml_bot: Optional[float] = None
    c_ref: Optional[float] = None
    dz_frac: Optional[float] = None
    n_pade: Optional[int] = None
    n_pade_start: Optional[int] = None
    theta_pade: Optional[float] = None
    power: Optional[bool] = None
    z_ref: Optional[float] = None
    alpha: Optional[float] = None
    z_roughness: Optional[float] = None
    external_wave: Optional[bool] = None
    sinus: Optional[bool] = None
    A_wave: Optional[float] = None
    lambda_wave: Optional[float] = None
    phi_wave: Optional[float] = None
    fname_wave: Optional[str] = None
    Xs_wave: Optional[float] = None
    Ys_wave: Optional[float] = None
    abs_air: Optional[bool] = None
    abs_water: Optional[bool] = None
    rh: Optional[float] = None
    dx_out: Optional[float] = None
    dz_out: Optional[float] = None


def _fields_from(first: str, last: str) -> frozenset[str]:
    names = [f.name for f in dataclasses.fields(WapeInputs)]
    return frozenset(names[names.index(first) : names.index(last) + 1])


PIECEWAPE_ONLY = _fields_from("gamma0", "yplane")
EWAPE_ONLY = _fields_from("H", "dz_out")


def backend_fields(backend: str) -> frozenset[str]:
    """Names of the :class:`WapeInputs` fields understood by ``backend``."""
    if backend not in BACKENDS:
        raise ValueError(f"unknown backend {backend!r}, expected one of {BACKENDS}")
    other = EWAPE_ONLY if backend == "piecewape" else PIECEWAPE_ONLY
    return frozenset(f.name for f in dataclasses.fields(WapeInputs)) - other


@dataclass
class DomainConfig:
    """Global acoustic domain (replaces the old ``Mesh`` xmin/xmax/... values)."""

    xmin: float
    xmax: float
    ymin: float
    ymax: float
    zmax: float


@dataclass
class RotorConfig:
    """Discretise a rotor in ``n`` source heights (see :func:`discretize_rotor`)."""

    hub_height: float
    diameter: float
    n: int


def discretize_rotor(hub_height: float, rotor_diameter: float, nb_source_heights: int) -> np.ndarray:
    """Source heights at the centre of ``nb_source_heights`` slices of the rotor, rounded to integers."""
    dh = rotor_diameter / nb_source_heights
    hmin = hub_height - rotor_diameter / 2 + dh / 2
    hmax = hub_height + rotor_diameter / 2
    heights = np.arange(hmin, hmax, dh)
    return np.rint(heights).astype(int)


@dataclass
class PropagationConfig:
    """``propagation`` section.

    Example::

        propagation:
          backend: piecewape              # or ewape
          executable: ~/wape/2D_WAPE      # else environment variable WFN_WAPE
          work_dir: wape_runs             # one folder t<i>/<int(height)>/ per case
          positions: [[0, 0]]             # source positions (x, y)
          source_heights: [50, 90, 120]   # or rotor: {hub_height: 90, diameter: 120, n: 7}
          tau: [0, 90, 180, 270]          # optional, else wape.tau
          frequencies: [50, 100]          # optional, else wape.frequencies
          domain: {xmin: -1000, xmax: 1000, ymin: -1000, ymax: 1000, zmax: 300}
          output: results/deltaL          # <output>_<i>.nc
          wape:                           # Fortran namelist parameters (WapeInputs)
            case_name: c0
            dout: 5                       # piecewape only
            receiver_heights: [2]
            side: true

    Without ``domain`` the ``Lx1/Lx2/Ly1/Ly2/Lz`` values of ``wape`` are used as
    is for every source. Relative paths are relative to the input file directory.
    """

    backend: str = "piecewape"  # "piecewape" (onshore, topography) or "ewape" (offshore)
    wape: WapeInputs = field(default_factory=WapeInputs)
    executable: Optional[str] = None
    work_dir: str = "."
    positions: Optional[list[list[float]]] = None
    source_heights: Optional[list[float]] = None
    rotor: Optional[RotorConfig] = None
    tau: Optional[list[float]] = None
    frequencies: Optional[list[float]] = None
    domain: Optional[DomainConfig] = None
    output: Optional[str] = None
    base_dir: Optional[str] = None

    @classmethod
    def from_dict(cls, data: dict, base_dir: str | Path | None = None) -> "PropagationConfig":
        cfg = structure(dict(data), cls)
        if base_dir is not None:
            cfg.base_dir = str(base_dir)
        backend_fields(cfg.backend)  # validate
        return cfg

    @classmethod
    def from_yaml(cls, path: str | Path) -> "PropagationConfig":
        return cls.from_dict(get_section(load_yaml(path), "propagation"), base_dir=Path(path).parent)

    def path(self, p: str | None) -> Path | None:
        if p is None:
            return None
        p = Path(p).expanduser()
        if p.is_absolute() or self.base_dir is None:
            return p
        return Path(self.base_dir) / p

    # ------------------------------------------------------------ resolved values
    @property
    def work_path(self) -> Path:
        return self.path(self.work_dir)

    def executable_path(self) -> Path:
        exe = self.executable or os.environ.get("WFN_WAPE")
        if exe is None:
            raise FileNotFoundError("no WAPE executable: set propagation.executable or WFN_WAPE")
        exe = self.path(exe) if self.executable else Path(exe).expanduser()
        if not exe.is_file():
            raise FileNotFoundError(f"WAPE executable not found: {exe} (maybe it needs to be compiled)")
        return exe

    def source_positions(self) -> list[list[float]]:
        return self.positions if self.positions else [[0.0, 0.0]]

    def heights(self) -> list[float]:
        if self.source_heights:
            return list(self.source_heights)
        if self.rotor is not None:
            return discretize_rotor(self.rotor.hub_height, self.rotor.diameter, self.rotor.n).tolist()
        if self.wape.source_height is not None:
            return [self.wape.source_height]
        raise ValueError("no source height: set source_heights, rotor or wape.source_height")

    def angles(self) -> list[float]:
        if self.tau:
            return list(self.tau)
        if self.wape.tau:
            return list(self.wape.tau)
        w = self.wape
        if None not in (w.taumin, w.taumax, w.nb_tau):
            return np.linspace(w.taumin, w.taumax, w.nb_tau).tolist()  # same as the Fortran codes
        raise ValueError("no propagation angle: set tau, wape.tau or wape.taumin/taumax/nb_tau")

    def freqs(self) -> list[float]:
        return list(self.frequencies) if self.frequencies else list(self.wape.frequencies)

    def output_file(self, isource: int) -> Path:
        out = self.output or "deltaL"
        return self.path(f"{out}_{isource}.nc")
