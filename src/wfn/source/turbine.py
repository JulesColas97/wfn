"""Wind turbine geometry and operating condition.

``WindTurbine`` only describes the turbine (blade discretisation, rotor
parameters and controller). It does not read any flow and does not hold any
boundary layer data (see :mod:`wfn.source.airfoil`).
"""
from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field, replace
from importlib import resources
from pathlib import Path
from typing import Optional

import numpy as np
import scipy.interpolate as interp
import yaml

from ..config import load_yaml, structure, unstructure

logger = logging.getLogger(__name__)

TURBINE_DATA = resources.files("wfn.source") / "data" / "turbines"
AIRFOIL_DATA = resources.files("wfn.source") / "data" / "airfoil_geometry"

# keys of the old ``wind_turbine`` YAML section and their new name (None: ignored)
_LEGACY_KEYS = {
    "Nblade": "n_blade",
    "Nbeta": "n_beta",
    "Nseg": None,
    "controler": "controller",
    "airfoil_data_file": "airfoil_database",
    "boundaryLayerData": None,
    "yaml_fname": None,
}


@dataclass
class Controller:
    """Rotational speed and pitch as a function of the hub wind speed.

    Either tabulated (``uinf``, ``omega``, ``pitch``; nearest value is used) or
    defined by the three points ``cutin``, ``rated``, ``cutout`` given as
    ``[u, omega_rpm, pitch_deg]``.
    """

    uinf: Optional[np.ndarray] = None
    omega: Optional[np.ndarray] = None
    pitch: Optional[np.ndarray] = None
    cutin: Optional[np.ndarray] = None
    rated: Optional[np.ndarray] = None
    cutout: Optional[np.ndarray] = None

    def get_omega(self, u: float) -> float:
        """Rotational speed (rpm) for the hub wind speed ``u``."""
        if self.uinf is not None and self.omega is not None:
            return float(self.omega[np.argmin(np.abs(self.uinf - u))])
        if self.cutin is None or self.rated is None or self.cutout is None:
            return 0.0
        u_in, omega_in = self.cutin[0:2]
        u_rat, omega_rat = self.rated[0:2]
        u_out = self.cutout[0]
        if u < u_in or u > u_out:
            return 0.0
        if u >= u_rat:
            return float(omega_rat)
        return float(omega_in + (u - u_in) * (omega_rat - omega_in) / (u_rat - u_in))

    def get_pitch(self, u: float) -> float:
        """Pitch angle (deg) for the hub wind speed ``u``."""
        if self.uinf is not None and self.pitch is not None:
            return float(self.pitch[np.argmin(np.abs(self.uinf - u))])
        if self.rated is None or self.cutout is None:
            return 0.0
        u_rat, pitch_rat = self.rated[[0, 2]]
        u_out, pitch_out = self.cutout[[0, 2]]
        if u < u_rat or u > u_out:
            return 0.0
        return float(pitch_rat + (u - u_rat) * (pitch_out - pitch_rat) / (u_out - u_rat))


@dataclass
class WindTurbine:
    """Wind turbine blade discretisation and rotor parameters.

    Attributes:
        seg: radial position of the centre of each blade segment (m).
        span: span of each segment (m).
        chord: chord of each segment (m).
        twist: twist of each segment (deg).
        airfoil_index: index in ``airfoil_name`` of the airfoil of each segment.
        airfoil_name: airfoil names; the geometry is looked for in
            ``airfoil_path`` then in the package data (``<name>.dat``).
        airfoil_path: optional directory containing the airfoil ``.dat`` files.
        hub_height: hub height (m), used when no turbine position is given.
        rotor_diameter: if given and different from ``2*(seg[-1] + span[-1]/2)``,
            the blade is rescaled (chords unchanged).
        n_blade: number of blades.
        n_beta: number of azimuthal positions per ``1/n_blade`` of a turn.
        omega_rpm: rotational speed (rpm), used when there is no controller.
        pitch: pitch angle (deg), used when there is no controller.
        yaw_deg: yaw angle (deg) of the rotor axis from the x axis.
        controller: optional controller giving omega and pitch from the hub wind speed.
        airfoil_database: optional path of the boundary layer database of this turbine.
    """

    seg: np.ndarray
    chord: np.ndarray
    airfoil_index: np.ndarray
    span: Optional[np.ndarray] = None
    twist: Optional[np.ndarray] = None
    airfoil_name: list[str] = field(default_factory=list)
    airfoil_path: Optional[str] = None
    hub_height: Optional[float] = None
    rotor_diameter: Optional[float] = None
    n_blade: int = 3
    n_beta: int = 12
    omega_rpm: float = 0.0
    pitch: float = 0.0
    yaw_deg: float = 0.0
    controller: Optional[Controller] = None
    airfoil_database: Optional[str] = None

    def __post_init__(self):
        self.seg = np.asarray(self.seg, dtype=float)
        self.chord = np.asarray(self.chord, dtype=float)
        self.airfoil_index = np.asarray(self.airfoil_index, dtype=int)
        if self.span is not None:
            self.span = np.asarray(self.span, dtype=float)
        if self.twist is not None:
            self.twist = np.asarray(self.twist, dtype=float)
            _check_twist_units(self)
        self.airfoil_name = [str(n) for n in self.airfoil_name]
        n = self.seg.size
        for name in ("span", "chord", "twist", "airfoil_index"):
            value = getattr(self, name)
            if value is not None and value.size != n:
                raise ValueError(f"'{name}' must have the same size as 'seg' ({n})")
        half_last = 0.0 if self.span is None else 0.5 * self.span[-1]
        geometric_diameter = float(2 * (self.seg[-1] + half_last))
        if self.rotor_diameter is None:
            self.rotor_diameter = geometric_diameter
        elif not np.isclose(self.rotor_diameter, geometric_diameter):
            target = self.rotor_diameter
            self.rotor_diameter = geometric_diameter
            self.rescale(target)

    # --------------------------------------------------------------- props
    @property
    def n_seg(self) -> int:
        return self.seg.size

    @property
    def n_azimuth(self) -> int:
        """Total number of azimuthal positions over a full turn."""
        return self.n_beta * self.n_blade

    @property
    def omega(self) -> float:
        """Rotational speed (rad/s)."""
        return self.omega_rpm * 2 * np.pi / 60

    @property
    def yaw_rad(self) -> float:
        return np.deg2rad(self.yaw_deg)

    @property
    def twist_cotte_conv(self) -> np.ndarray:
        """Twist in the convention of Cotte/Tian (``90 - twist``), in degrees."""
        return 90.0 - self.twist

    @property
    def azimuths(self) -> np.ndarray:
        """Azimuth (rad) of the ``n_beta * n_blade`` blade positions over a full turn.

        Index ``k = iblade * n_beta + ibeta``. Azimuth 0 is a blade pointing
        upwards, increasing clockwise when looking downstream.
        """
        return 2 * np.pi * np.arange(self.n_azimuth) / self.n_azimuth

    def airfoil_files(self) -> list[Path]:
        """Paths of the airfoil geometry files (``airfoil_name`` order)."""
        out = []
        for name in self.airfoil_name:
            candidates = []
            if self.airfoil_path is not None:
                candidates.append(Path(self.airfoil_path) / f"{name}.dat")
            candidates.append(Path(str(AIRFOIL_DATA / f"{name}.dat")))
            for c in candidates:
                if c.is_file():
                    out.append(c)
                    break
            else:
                raise FileNotFoundError(f"airfoil geometry '{name}.dat' not found in {candidates}")
        return out

    # ------------------------------------------------------- modifications
    def rescale(self, new_diameter: float) -> None:
        """Rescale ``seg`` and ``span`` to a new rotor diameter (chords unchanged)."""
        ratio = new_diameter / self.rotor_diameter
        self.seg = ratio * self.seg
        if self.span is not None:
            self.span = ratio * self.span
        self.rotor_diameter = float(new_diameter)
        logger.info("rotor diameter set to %s m", new_diameter)

    def check_acoustic_geometry(self) -> None:
        """Raise if the geometry needed by the noise model is incomplete."""
        missing = [k for k in ("span", "twist") if getattr(self, k) is None]
        if missing:
            raise ValueError(f"turbine geometry needs {missing} for the noise computation")

    def with_operating_point(self, u_hub: float) -> "WindTurbine":
        """Copy of the turbine with omega and pitch given by the controller for ``u_hub``."""
        if self.controller is None:
            return self
        return replace(
            self,
            omega_rpm=self.controller.get_omega(u_hub),
            pitch=self.controller.get_pitch(u_hub),
        )

    def acoustic_strips(self, seg0: float | None = None, target_l_c: float = 3.0,
                        delta_span: float = 0.1) -> "WindTurbine":
        """New blade discretisation with segments of span/chord ratio ``target_l_c``.

        Port of the old ``createStripForAcoustic``: each segment span grows by
        ``2*delta_span`` until the aspect ratio target is reached, then a new
        segment starts. Chord and twist are interpolated with quadratic splines.
        """
        if self.twist is None:
            raise ValueError("twist is needed to create acoustic strips")
        if seg0 is None:
            seg0 = 0.0 if self.span is None else self.seg[0] - 0.5 * self.span[0]
        tck_c = interp.splrep(self.seg, self.chord, k=2)
        tck_twist = interp.splrep(self.seg, self.twist, k=2)
        rotor_radius = 0.5 * self.rotor_diameter

        seg, span, twist, chord, airfoil_index = [], [], [], [], []
        current_blade_length = seg0
        current_span = 0.0
        current_seg = seg0
        current_c = float(interp.splev(current_seg, tck_c))
        while current_blade_length < rotor_radius:
            if current_span / current_c < target_l_c:
                current_span += 2 * delta_span
                current_blade_length += 2 * delta_span
                current_seg += delta_span
                current_c = float(interp.splev(current_seg, tck_c))
            else:
                seg.append(current_seg)
                span.append(current_span)
                chord.append(current_c)
                twist.append(float(interp.splev(current_seg, tck_twist)))
                airfoil_index.append(self.airfoil_index[np.abs(current_seg - self.seg).argmin()])
                current_seg += 0.5 * current_span
                current_span = 0.0
                current_c = float(interp.splev(current_seg, tck_c))

        if not seg:
            raise ValueError("blade too short for the requested aspect ratio")
        seg[-1] += 0.5 * current_span
        span[-1] += current_span
        twist[-1] = float(interp.splev(seg[-1], tck_twist))
        chord[-1] = float(interp.splev(seg[-1], tck_c))
        airfoil_index[-1] = self.airfoil_index[np.abs(self.seg - current_seg).argmin()]

        return replace(
            self,
            seg=np.round(seg, 2),
            span=np.round(span, 2),
            twist=np.round(twist, 2),
            chord=np.round(chord, 2),
            airfoil_index=np.array(airfoil_index),
            rotor_diameter=None,
        )

    # ------------------------------------------------------------------ I/O
    @classmethod
    def from_dict(cls, data: dict) -> "WindTurbine":
        """Build from a dict; keys of the old ``wind_turbine`` section are accepted."""
        data = dict(data)
        for old, new in _LEGACY_KEYS.items():
            if old in data:
                value = data.pop(old)
                if new is not None and new not in data:
                    data[new] = value
        # old files could give the airfoils as paths instead of names
        if "airfoil_geometry" in data:
            geom = data.pop("airfoil_geometry")
            if "airfoil_name" not in data:
                data["airfoil_name"] = [Path(g).stem for g in geom]
                parents = {str(Path(g).parent) for g in geom}
                if len(parents) == 1 and "airfoil_path" not in data and parents != {"."}:
                    data["airfoil_path"] = parents.pop()
        return structure(data, cls)

    @classmethod
    def from_yaml(cls, path: str | Path, section: str | None = None) -> "WindTurbine":
        """Read a turbine file.

        The turbine may be at the root of the file, under ``wind_turbine`` (old
        format) or under ``section``. Relative ``airfoil_path`` and
        ``airfoil_database`` are resolved relative to the file.
        """
        path = Path(path)
        data = load_yaml(path)
        if section is not None:
            data = data[section]
        elif "wind_turbine" in data:
            data = data["wind_turbine"]
        wt = cls.from_dict(data)
        for attr in ("airfoil_path", "airfoil_database"):
            value = getattr(wt, attr)
            if value is not None and not Path(value).is_absolute():
                setattr(wt, attr, str(path.parent / value))
        return wt

    @classmethod
    def default(cls, name: str) -> "WindTurbine":
        """Turbine shipped with the package (``NREL_5MW``, ``IEA15MW``, ``SWT9.3MW``)."""
        f = TURBINE_DATA / f"{name}.yaml"
        if not f.is_file():
            available = sorted(p.name[:-5] for p in TURBINE_DATA.iterdir() if p.name.endswith(".yaml"))
            raise FileNotFoundError(f"no default turbine '{name}', available: {available}")
        return cls.from_yaml(Path(str(f)))

    @classmethod
    def load(cls, spec: str | dict, base_dir: str | Path | None = None) -> "WindTurbine":
        """Build from an inline dict, a ``{file: ...}`` dict, a file path or a default name."""
        if isinstance(spec, dict):
            if "file" in spec:
                return cls.load(spec["file"], base_dir)
            return cls.from_dict(spec)
        p = Path(spec)
        if not p.is_absolute() and base_dir is not None and (Path(base_dir) / p).is_file():
            p = Path(base_dir) / p
        if p.is_file():
            return cls.from_yaml(p)
        return cls.default(p.stem if p.suffix in (".yaml", ".yml") else str(spec))

    def to_dict(self) -> dict:
        return _plain({k: v for k, v in unstructure(self).items() if v is not None})

    def to_yaml(self, path: str | Path) -> None:
        with open(path, "w") as f:
            yaml.safe_dump({"wind_turbine": self.to_dict()}, f, sort_keys=False)


def _plain(obj):
    """Recursively convert numpy types into plain python (YAML safe)."""
    if isinstance(obj, dict):
        return {k: _plain(v) for k, v in obj.items() if v is not None}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def _check_twist_units(wt: WindTurbine) -> None:
    if np.all(np.abs(wt.twist) < 0.5) and np.any(wt.twist != 0):
        warnings.warn("twist values are all below 0.5: they must be given in degrees, not radians")
