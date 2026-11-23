"""Boundary layer database of a blade, stored as an ``xarray.Dataset``.

Variables (dims ``(uinf, aoa, seg)``): ``delta_star_{top,bot}``,
``theta_momen_{top,bot}``, ``cf_{top,bot}``, ``dpdx_{top,bot}``,
``UeUinf_{top,bot}``, ``cl``, ``cd``. ``top`` is the suction side, ``bot`` the
pressure side; quantities are taken at ``target_frac`` of the chord.
Coordinates: ``uinf`` (m/s), ``aoa`` (deg), ``seg`` (index) with ``chord(seg)``.

Interpolation in ``(uinf, aoa)`` is linear inside the convex hull of the valid
(non-NaN, non-zero) points and nearest-neighbour outside, as in the old code.
"""
from __future__ import annotations

import logging
import pickle
from pathlib import Path
from typing import Sequence

import numpy as np
import xarray as xr
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

from . import xfoil as _xfoil

logger = logging.getLogger(__name__)

BL_VARS = (
    "delta_star_top", "delta_star_bot",
    "theta_momen_top", "theta_momen_bot",
    "cf_top", "cf_bot",
    "dpdx_top", "dpdx_bot",
    "UeUinf_top", "UeUinf_bot",
)
AERO_VARS = ("cl", "cd")
ALL_VARS = BL_VARS + AERO_VARS

# xfoil output key -> database variable
_XFOIL_KEYS = {
    "Dstar_suction": "delta_star_top", "theta_suction": "theta_momen_top",
    "cf_suction": "cf_top", "dpdx_suction": "dpdx_top", "Ue_suction": "UeUinf_top",
    "Dstar_pressure": "delta_star_bot", "theta_pressure": "theta_momen_bot",
    "cf_pressure": "cf_bot", "dpdx_pressure": "dpdx_bot", "Ue_pressure": "UeUinf_bot",
}


class LinearNDInterpolatorExt:
    """Linear interpolation on scattered valid points, nearest neighbour outside.

    Port of ``prepost/interpolation2D.py``. NaN and exact zeros are treated as
    missing values.
    """

    def __init__(self, x: np.ndarray, y: np.ndarray, z: np.ndarray):
        xx, yy = np.meshgrid(x, y, indexing="ij")
        mask = ~np.isnan(z) & (z != 0)
        if not mask.any():
            raise ValueError("no valid data to interpolate")
        points = np.column_stack((xx[mask], yy[mask]))
        values = z[mask]
        self.linear = LinearNDInterpolator(points, values, rescale=True)
        self.nearest = NearestNDInterpolator(points, values, rescale=True)

    def __call__(self, *args):
        z = self.linear(*args)
        nan = np.isnan(z)
        if nan.any():
            return np.where(nan, self.nearest(*args), z)
        return z


class AirfoilDatabase:
    """Boundary layer and aerodynamic data of each blade segment."""

    def __init__(self, ds: xr.Dataset):
        missing = [v for v in BL_VARS if v not in ds]
        if missing:
            raise ValueError(f"airfoil database is missing variables {missing}")
        self.ds = ds
        self._interp: dict[tuple[str, int], LinearNDInterpolatorExt] = {}

    # ------------------------------------------------------------- creation
    @classmethod
    def from_arrays(cls, uinf, aoa, chord, data: dict[str, np.ndarray], **attrs) -> "AirfoilDatabase":
        """Build from arrays of shape ``(n_uinf, n_aoa, n_seg)``."""
        uinf = np.asarray(uinf, dtype=float)
        aoa = np.asarray(aoa, dtype=float)
        chord = np.asarray(chord, dtype=float)
        coords = {
            "uinf": ("uinf", uinf, {"units": "m/s"}),
            "aoa": ("aoa", aoa, {"units": "deg"}),
            "seg": ("seg", np.arange(chord.size)),
            "chord": ("seg", chord, {"units": "m"}),
        }
        dvars = {k: (("uinf", "aoa", "seg"), np.asarray(v, dtype=float))
                 for k, v in data.items() if k in ALL_VARS and v is not None}
        attrs = {k: v for k, v in attrs.items() if v is not None}
        return cls(xr.Dataset(dvars, coords=coords, attrs=attrs))

    @classmethod
    def from_legacy_pickle(cls, path: str | Path) -> "AirfoilDatabase":
        """Read an ``AirfoilData.save`` pickle file of the old code."""
        with open(path, "rb") as f:
            d = pickle.load(f)
        if d.get("INTERPOLAND"):
            raise ValueError("this legacy file was saved after interpolation, raw data is lost")
        if d.get("theta_momen_bot") is None and d.get("theta_mom_bot") is not None:
            d["theta_momen_bot"] = d["theta_mom_bot"]
        attrs = {
            "target_frac": d.get("target_frac"),
            "t_inf": d.get("t_inf"),
            "rho0": d.get("rho0"),
            "nu": d.get("nu"),
            "source": f"converted from legacy file {Path(path).name}",
        }
        db = cls.from_arrays(d["uinf"], d["aoa"], d["chord"], d, **attrs)
        if d.get("airfoil_index") is not None:
            db.ds["airfoil_index"] = ("seg", np.asarray(d["airfoil_index"], dtype=int))
        if d.get("airfoil_name") is not None:
            db.ds.attrs["airfoil_name"] = [str(n) for n in d["airfoil_name"]]
        return db

    @classmethod
    def open(cls, path: str | Path) -> "AirfoilDatabase":
        """Open a netCDF database, or an old pickle ``.dat`` file."""
        path = Path(path)
        if path.suffix in (".nc", ".h5", ".nc4"):
            return cls(xr.load_dataset(path))
        return cls.from_legacy_pickle(path)

    def save(self, path: str | Path) -> None:
        self.ds.to_netcdf(path)

    @classmethod
    def build_with_xfoil(
        cls,
        chord: Sequence[float],
        airfoil_index: Sequence[int],
        airfoil_files: Sequence[str | Path],
        uinf: Sequence[float],
        aoa: Sequence[float],
        nu: float = 1.48e-5,
        rho0: float = 1.2,
        t_inf: float = 273.0,
        target_frac: float = 0.95,
        trip_top: float | None = None,
        trip_bot: float | None = None,
        virtual_buffer: bool = False,
        display_id: int = 99,
    ) -> "AirfoilDatabase":
        """Compute the database with XFoil (port of ``AirfoilData.createBLData``).

        For each segment and each velocity, XFoil runs over all ``aoa`` (in the
        given order; starting from 0 and moving outwards helps convergence).
        The Reynolds number per unit chord is ``uinf / nu``. Duplicated AoA are
        removed and AoA are sorted at the end.

        Use :meth:`for_turbine` to build it directly from a ``WindTurbine``.
        """
        chord = np.asarray(chord, dtype=float)
        airfoil_index = np.asarray(airfoil_index, dtype=int)
        uinf = np.asarray(uinf, dtype=float)
        aoa = np.asarray(aoa, dtype=float)
        n_u, n_a, n_s = uinf.size, aoa.size, chord.size
        data = {v: np.zeros((n_u, n_a, n_s)) for v in BL_VARS}
        data["cl"] = np.full((n_u, n_a, n_s), np.nan)
        data["cd"] = np.full((n_u, n_a, n_s), np.nan)
        c0 = (_xfoil.GAMMA * _xfoil.R_AIR * t_inf) ** 0.5

        def run_all(env):
            for iseg in range(n_s):
                foil = airfoil_files[airfoil_index[iseg]]
                for iu in range(n_u):
                    logger.info("xfoil: segment %d/%d, velocity %d/%d", iseg + 1, n_s, iu + 1, n_u)
                    inputs = dict(chord=chord[iseg], re=uinf[iu] / nu, mach=uinf[iu] / c0, aoa=aoa,
                                  trip_top=trip_top, trip_bot=trip_bot, target_frac=target_frac,
                                  rho=rho0, T=t_inf)
                    try:
                        out = _xfoil.run_xfoil(inputs, foil, env=env)
                    except _xfoil.XfoilError as e:
                        logger.warning("xfoil failed for segment %d, uinf %s: %s", iseg, uinf[iu], e)
                        continue
                    for key, var in _XFOIL_KEYS.items():
                        data[var][iu, :, iseg] = out[key]
                    # kept from the old code: XFoil polars are divided by the chord
                    data["cl"][iu, :, iseg] = out["cl"] / chord[iseg]
                    data["cd"][iu, :, iseg] = out["cd"] / chord[iseg]

        if virtual_buffer:
            with _xfoil.VirtualDisplay(display_id) as disp:
                run_all(disp.env)
        else:
            run_all(None)

        db = cls.from_arrays(uinf, aoa, chord, data, target_frac=target_frac, t_inf=t_inf,
                             rho0=rho0, nu=nu, trip_top=trip_top, trip_bot=trip_bot,
                             source="xfoil")
        db.ds["airfoil_index"] = ("seg", airfoil_index)
        db.ds.attrs["airfoil_name"] = [Path(f).stem for f in airfoil_files]
        return db.sorted_aoa()

    @classmethod
    def for_turbine(cls, turbine, uinf, aoa, **kwargs) -> "AirfoilDatabase":
        """Run :meth:`build_with_xfoil` on the segments of a :class:`WindTurbine`."""
        return cls.build_with_xfoil(turbine.chord, turbine.airfoil_index, turbine.airfoil_files(),
                                    uinf, aoa, **kwargs)

    def sorted_aoa(self) -> "AirfoilDatabase":
        """Remove duplicated angles of attack (first kept) and sort them."""
        _, first = np.unique(self.ds.aoa.values, return_index=True)
        return AirfoilDatabase(self.ds.isel(aoa=np.sort(first)).sortby("aoa"))

    # -------------------------------------------------------- interpolation
    @property
    def n_seg(self) -> int:
        return self.ds.sizes["seg"]

    @property
    def chord(self) -> np.ndarray:
        return self.ds.chord.values

    def interpolator(self, var: str, iseg: int) -> LinearNDInterpolatorExt:
        key = (var, iseg)
        if key not in self._interp:
            self._interp[key] = LinearNDInterpolatorExt(
                self.ds.uinf.values, self.ds.aoa.values,
                self.ds[var].isel(seg=iseg).transpose("uinf", "aoa").values)
        return self._interp[key]

    def evaluate(self, var: str, iseg: int, uinf, aoa) -> np.ndarray:
        """Interpolate ``var`` of segment ``iseg`` at ``(uinf, aoa)`` (arrays allowed)."""
        return np.asarray(self.interpolator(var, iseg)(uinf, aoa), dtype=float)

    def boundary_layer(self, uinf: np.ndarray, aoa: np.ndarray) -> dict[str, np.ndarray]:
        """All BL variables at the operating points of every segment.

        Args:
            uinf, aoa: arrays of shape ``(..., n_seg)``.

        Returns:
            dict ``var -> array`` with the same shape as ``uinf``.
        """
        uinf = np.asarray(uinf, dtype=float)
        aoa = np.asarray(aoa, dtype=float)
        if uinf.shape[-1] != self.n_seg:
            raise ValueError(f"last axis must be the {self.n_seg} segments, got {uinf.shape}")
        out = {}
        for var in BL_VARS:
            vals = np.empty_like(uinf)
            for iseg in range(self.n_seg):
                vals[..., iseg] = self.evaluate(var, iseg, uinf[..., iseg], aoa[..., iseg])
            out[var] = vals
        return out

    def check_turbine(self, turbine, rtol: float = 1e-2) -> None:
        """Raise if the database segments do not match the turbine segments."""
        if self.n_seg != turbine.n_seg:
            raise ValueError(f"database has {self.n_seg} segments, turbine has {turbine.n_seg}")
        if not np.allclose(self.chord, turbine.chord, rtol=rtol):
            raise ValueError("database chords differ from turbine chords")
