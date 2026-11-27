"""Content of the ``propagation`` section of the input file."""
from __future__ import annotations

import dataclasses
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar, Optional

import numpy as np

from ..config import get_section, load_yaml, structure


@dataclass
class WapeInputs:
    """Namelist parameters shared by the Fortran WAPE codes (``piecewape`` and ``EWAPE``).

    Use the subclass of the backend: :class:`PiecewapeInputs` or :class:`EwapeInputs`
    (``WapeInputs.for_backend(name)``). The field names are the keys of the Fortran
    namelist ``input.nml``: do not rename them. The default values are the ones of the
    Fortran codes (``EWAPE/input_params.f90`` and ``piecewape/read_write.f90``): keep them in sync.
    """

    backend: ClassVar[str] = ""
    _registry: ClassVar[dict[str, type["WapeInputs"]]] = {}

    def __init_subclass__(cls, backend: str | None = None, **kwargs):
        super().__init_subclass__(**kwargs)
        if backend is not None:
            cls.backend = backend
            WapeInputs._registry[backend] = cls

    @classmethod
    def for_backend(cls, backend: str) -> type["WapeInputs"]:
        """Inputs class of ``backend`` (``"piecewape"`` or ``"ewape"``)."""
        if backend not in cls._registry:
            raise ValueError(f"unknown backend {backend!r}, expected one of {tuple(cls._registry)}")
        return cls._registry[backend]

    def local_domain_extra(self, xs: float, ys: float) -> dict:
        """Backend specific keys that depend on the source position."""
        return {}

    def z_levels(self, Lz: float) -> np.ndarray:
        """Heights of the vertical plane written by the code, for a domain of height ``Lz``."""
        return np.arange(0, Lz, self.dx_out)

    # ---- shared by both codes
    case_name: str = "c0"
    c_air: float = 340.0
    rho_air: float = 1.2
    T0: float = 290.0  # temperature [K]
    p0: float = 101325.0  # pressure [Pa]
    frequencies: list[float] = field(default_factory=lambda: [100.0])
    tau: list[float] = field(default_factory=lambda: [0.0])
    taumin: Optional[float] = None  # unset: used with taumax and nb_tau only if all three are set
    taumax: Optional[float] = None
    nb_tau: Optional[int] = None
    Lx1: float = 300.0
    Lx2: float = 300.0
    Ly1: float = 300.0
    Ly2: float = 300.0
    Lz: float = 100.0
    source_height: float = 119.0
    dx_frac: float = 2.0  # points per wavelength: dx = lambda / dx_frac
    pml_size: int = 30
    pml_param: float = 5e4
    pml_n: float = 2.5
    # flow
    external_flow: bool = False
    uniform: bool = False
    logarithmic: bool = False
    u_ref: float = 0.0
    fname_flow: str = "./u.h5"
    Xs_flow: float = 0.0
    Ys_flow: float = 0.0
    # output
    dx_out: float = 1.0  # output step [m] (in x and z for piecewape, in x only for EWAPE)
    receiver_heights: list[float] = field(default_factory=lambda: [2.0])
    side: bool = True


@dataclass
class PiecewapeInputs(WapeInputs, backend="piecewape"):
    """Namelist parameters of ``piecewape`` (onshore, topography), on top of :class:`WapeInputs`."""

    gamma0: float = 1.4
    pe_model: int = 0
    multi_domain: bool = False
    domain_step: float = 0.0
    impedance_mode: int = 0  # 0: rigid ground
    impedance_sigmae: float = 50000.0
    impedance_alphae: float = 100.0
    impedance_Rs: float = 7000.0
    impedance_Om: float = 0.6
    impedance_q2: float = 1.3
    impedance_e: float = 0.13
    external_terrain: bool = False
    fname_terrain: str = "./terrain.h5"
    Xs_terrain: float = 0.0
    Ys_terrain: float = 0.0
    gaussian_hill: bool = False
    hill_length: float = 0.0
    hill_pos: float = 0.0
    hill_height: float = 0.0
    continuation: bool = False
    write_flow: bool = True  # write the mean flow (/flow group) in the output file
    xplane: list[float] = field(default_factory=list)
    yplane: list[float] = field(default_factory=list)

    def local_domain_extra(self, xs: float, ys: float) -> dict:
        return dict(Xs_terrain=xs, Ys_terrain=ys)


@dataclass
class EwapeInputs(WapeInputs, backend="ewape"):
    """Namelist parameters of ``EWAPE`` (offshore), on top of :class:`WapeInputs`."""

    H: float = 100.0
    c_water: float = 1500.0
    rho_water: float = 1000.0
    eta_air: float = 0.0
    eta_water: float = 0.0
    seabed_type: int = 0
    rho_sb: float = 1900.0
    cp_sb: float = 1700.0
    cs_sb: float = 0.0
    alpha_p: float = 0.0
    alpha_s: float = 0.0
    theta_sb: float = 0.0
    d_pml_bot: float = 2.0
    c_ref: float = 850.0
    dz_frac: float = 20.0
    n_pade: int = 10
    n_pade_start: int = 10
    theta_pade: float = float(np.pi / 8)
    power: bool = False
    z_ref: float = 90.0
    alpha: float = 0.2
    z_roughness: float = 0.01
    external_wave: bool = False
    sinus: bool = False
    A_wave: float = 2.0
    lambda_wave: float = 80.0
    phi_wave: float = 0.0
    fname_wave: str = "./wave.h5"
    Xs_wave: float = 0.0
    Ys_wave: float = 0.0
    abs_air: bool = False
    abs_water: bool = False
    rh: float = 80.0
    dz_out: float = 1.0

    def local_domain_extra(self, xs: float, ys: float) -> dict:
        return dict(Xs_wave=xs, Ys_wave=ys)

    def z_levels(self, Lz: float) -> np.ndarray:
        # side plane: z = -H, -H + dz_out, ... (water below z = 0), see EWAPE/mesh_module.f90
        return -self.H + self.dz_out * np.arange(int(np.floor((self.H + Lz) / self.dz_out)))


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
            dx_out: 5
            receiver_heights: [2]
            side: true

    Without ``domain`` the ``Lx1/Lx2/Ly1/Ly2/Lz`` values of ``wape`` are used as
    is for every source. Relative paths are relative to the input file directory.
    """

    backend: str = "piecewape"  # "piecewape" (onshore, topography) or "ewape" (offshore)
    wape: Optional[WapeInputs] = None  # None: defaults of the backend (see ``__post_init__``)
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

    def __post_init__(self):
        inputs_cls = WapeInputs.for_backend(self.backend)  # validates the backend
        if self.wape is None:
            self.wape = inputs_cls()
        elif not isinstance(self.wape, inputs_cls):
            raise TypeError(f"backend {self.backend!r} needs {inputs_cls.__name__}, got {type(self.wape).__name__}")

    @classmethod
    def from_dict(cls, data: dict, base_dir: str | Path | None = None) -> "PropagationConfig":
        data = dict(data)
        wape = data.pop("wape", None) or {}
        inputs_cls = WapeInputs.for_backend(data.get("backend", cls.backend))
        cfg = structure(data, cls)  # ``wape`` is built apart: its class depends on the backend
        cfg.wape = structure(dict(wape), inputs_cls)  # unknown keys (e.g. of the other backend) raise
        if base_dir is not None:
            cfg.base_dir = str(base_dir)
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
        return [self.wape.source_height]

    def angles(self) -> list[float]:
        if self.tau:
            return list(self.tau)
        w = self.wape
        if None not in (w.taumin, w.taumax, w.nb_tau):  # same priority as in the Fortran codes
            # taumax is excluded (full circle: 0, 360, 4 -> 0, 90, 180, 270), same operations as the Fortran
            return (w.taumin + np.arange(w.nb_tau) * ((w.taumax - w.taumin) / w.nb_tau)).tolist()
        return list(w.tau)

    def freqs(self) -> list[float]:
        return list(self.frequencies) if self.frequencies else list(self.wape.frequencies)

    def output_file(self, isource: int) -> Path:
        out = self.output or "deltaL"
        return self.path(f"{out}_{isource}.nc")
