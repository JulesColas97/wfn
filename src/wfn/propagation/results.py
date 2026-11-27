"""Read WAPE results and assemble them into ``xr.Dataset`` objects.

* :func:`read_wape`: one WAPE output file (one source height, one angle) with
  ``deltaL(r, z, freq)``.
* :func:`assemble_polar`: all the angles and source heights of one source, with
  ``deltaL(r, tau, z, freq, source_height)`` and ``terrain(r, tau)``.
* :func:`to_grid`: interpolate a polar result onto any receiver grid.
* :func:`third_octave`: average the frequencies into third-octave bands.

The classes :class:`PeResults` (one WAPE file) and :class:`DeltaL` (all the runs of one source)
wrap these functions and the resulting ``xr.Dataset`` (attribute ``ds``).
"""
from __future__ import annotations

import logging
from pathlib import Path
import json

import h5py
import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator

from .. import receivers
from ..units import third_octave_frequencies
from .cases import case_dir, local_domain
from .config import PropagationConfig
from wfn.receivers import ReceiversConfig
from wfn.config import unstructure

logger = logging.getLogger(__name__)

DL_MAX = 120.0  # |deltaL| above this value is set to NaN (as in the old code)
R_MARGIN = 50.0


def wape_filename(directory: str | Path, case_name: str, tau: float) -> Path:
    """Name of the WAPE output file: ``<case_name>_<tau>.h5`` (``tau`` as ``%04d`` if integer)."""
    tau_str = format(int(tau), "04d") if int(tau) == tau else str(tau)
    return Path(directory) / f"{case_name}_{tau_str}.h5"


def read_wape(fname: str | Path, side: bool = True) -> xr.Dataset:
    """Read one WAPE output file.

    Args:
        fname: ``.h5`` file written by WAPE.
        side: read the full vertical plane (``solution`` group) if True, else
            only the receiver heights (``receiver`` group).

    Returns:
        Dataset with ``deltaL(r, z, freq)``, the ground height ``h(r)``, the
        absolute height ``z_abs(r, z)`` of every point and, for the side plane, the
        mean flow ``u(r, z)`` and sound speed ``c(r, z)`` when the file has them.
    """
    flow = {}
    with h5py.File(fname, "r") as f:
        mesh = f["mesh"]
        r = np.array(mesh["x"])
        h = np.array(mesh["h"]) if "h" in mesh else np.zeros_like(r)
        s = np.array(mesh["s"]) if "s" in mesh else np.zeros_like(r)

        # same layout for EWAPE and piecewape: /side/<freq>/TL, /receiver/<freq>/TL, /flow/{u,c}
        # (old files: /solution/<freq>/deltaL and /receiver/<freq>/deltaL)
        if side:
            group = f["side"] if "side" in f else f["solution"]
        else:
            group = f["receiver"]
        keys = [k for k in group.keys() if k != "heights"]
        name = "TL" if "TL" in group[keys[0]] else "deltaL"
        if side:
            z = np.array(mesh["z"])
            z_abs = np.array(mesh["z_grid"]).T if "z_grid" in mesh else z[None, :] + h[:, None]
            if "flow" in f:
                flow = {k: (("r", "z"), np.array(f["flow"][k]).T, {"units": "m/s"}) for k in f["flow"]}
        else:
            z = np.array(group["heights"])
            z_abs = z[None, :] + h[:, None]
        dl = np.stack([np.array(group[k][name]).T for k in keys], axis=-1)

    order = np.argsort([float(k) for k in keys])
    freq = np.array([float(k) for k in keys])[order]
    dl = dl[..., order]
    dl[np.abs(dl) > DL_MAX] = np.nan
    return xr.Dataset(
        {
            "deltaL": (("r", "z", "freq"), dl, {"units": "dB"}),
            "h": ("r", h, {"units": "m"}),
            "s": ("r", s, {"units": "m"}),
            "z_abs": (("r", "z"), z_abs, {"units": "m"}),
            **flow,
        },
        coords={"r": r, "z": z, "freq": freq},
        attrs={"side": int(side), "field": name},
    )


def assemble_polar(cfg: PropagationConfig, isource: int = 0) -> xr.Dataset:
    """Gather all WAPE runs of one source on a common polar grid centred on the source.

    The grid has ``r = 0, dx_out, ... rmax`` with ``rmax`` the largest distance to
    a corner of the domain plus 50 m, and ``z = 0, dx_out, ..., < Lz`` (EWAPE: ``z = -H, -H + dz_out, ...``) (``side``) or
    the receiver heights. Beyond the end of a WAPE run, the last computed values
    are copied up to ``rmax``.

    Returns:
        Dataset with ``deltaL(r, tau, z, freq, source_height)`` and
        ``terrain(r, tau)``, plus the polar grid coordinates (``x``, ``y``).
    """
    w = cfg.wape
    side = bool(w.side)
    dom = {"Lx1": w.Lx1, "Lx2": w.Lx2, "Ly1": w.Ly1, "Ly2": w.Ly2, "Lz": w.Lz}
    dom.update(local_domain(cfg, isource))
    dr = w.dx_out
    rmax = np.sqrt(max(dom["Lx1"], dom["Lx2"]) ** 2 + max(dom["Ly1"], dom["Ly2"]) ** 2) + R_MARGIN
    r = np.arange(0, rmax + dr, dr)
    if not side:
        z = np.asarray(w.receiver_heights, dtype=float)
    else:
        z = w.z_levels(dom["Lz"])
    tau = np.asarray(cfg.angles(), dtype=float)
    heights = np.asarray(cfg.heights(), dtype=float)
    freq = np.asarray(cfg.freqs(), dtype=float)

    field = None
    dl = np.zeros((r.size, tau.size, z.size, freq.size, heights.size))
    terrain = np.full((r.size, tau.size), np.nan)
    surface = np.full((r.size, tau.size), np.nan)
    for ih, height in enumerate(heights):
        for it, t in enumerate(tau):
            res = read_wape(wape_filename(case_dir(cfg, isource, height), w.case_name, t), side)
            if field is None:
                field = res.attrs["field"]
            elif res.attrs["field"] != field:
                raise ValueError(f"WAPE run (height {height}, tau {t}) holds {res.attrs['field']}, others hold {field}")
            if side:
                # the WAPE side plane also contains the PML rows above Lz: crop them
                res = res.isel(z=slice(0, z.size))
            if res.sizes["z"] != z.size or not np.allclose(res.z.values, z):
                raise ValueError(
                    f"WAPE run (height {height}, tau {t}) has heights {res.z.values[[0, -1]]} "
                    f"({res.sizes['z']} values), expected {z.size} values from {z[0]} to {z[-1]}"
                )
            res_dl = res.deltaL.sel(freq=freq, method="nearest").values
            i0 = int(np.argmin(np.abs(r - res.r.values[0])))
            i1 = int(np.argmin(np.abs(r - res.r.values[-1])))
            dl[i0 : i1 + 1, it, :, :, ih] = res_dl
            terrain[i0 : i1 + 1, it] = res.h.values
            surface[i0 : i1 + 1, it] = res.s.values
            # copy the last values up to rmax
            last = i0 + np.where(res_dl > -180)[0][-1]
            dl[last:, it, :, :, ih] = dl[last - 1, it, :, :, ih]
            terrain[last:, it] = terrain[last - 1, it]
            surface[last:, it] = surface[last - 1, it]

    xs, ys = cfg.source_positions()[isource][:2]
    ds = receivers.polar(r, tau, z, center=(xs, ys))
    ds = ds.assign_coords(freq=("freq", freq, {"units": "Hz"}), source_height=("source_height", heights, {"units": "m"}))
    ds["deltaL"] = (("r", "tau", "z", "freq", "source_height"), dl, {"units": "dB"})
    ds["terrain"] = (("r", "tau"), terrain, {"units": "m"})
    ds["surface"] = (("r", "tau"), surface, {"units": "m"})
    topo = bool(np.nanmax(terrain) - np.nanmin(terrain) > 1e-6)
    ds.attrs.update(side=int(side), topo=int(topo), topo_following=int(topo and not side), isource=isource,
                    field=field)
    return ds


def _set_polar_xy(ds: xr.Dataset) -> xr.Dataset:
    grid = receivers.polar(ds.r.values, ds.tau.values, ds.z.values, center=(ds.attrs["center_x"], ds.attrs["center_y"]))
    return ds.assign_coords(x=grid.x, y=grid.y)


def loop_angle(ds: xr.Dataset) -> xr.Dataset:
    """Append ``tau = 360`` (copy of ``tau = 0``) so that interpolation is continuous around the source."""
    last = ds.isel(tau=[0]).assign_coords(tau=[360.0])
    return xr.concat([ds, last], dim="tau", data_vars="minimal", coords="minimal", compat="override")


def refine_angle(ds: xr.Dataset) -> xr.Dataset:
    """Insert the mid angle between each pair of angles (linear interpolation)."""
    tau = ds.tau.values
    new_tau = np.sort(np.concatenate([tau, 0.5 * (tau[1:] + tau[:-1])]))
    out = ds.drop_vars(["x", "y"]).interp(tau=new_tau)
    return _set_polar_xy(out.assign_attrs(ds.attrs))


def to_grid(ds: xr.Dataset, grid: xr.Dataset) -> xr.Dataset:
    """Interpolate a polar result (:func:`assemble_polar`) onto a receiver grid.

    Use :func:`loop_angle` first if the angles cover the full circle, otherwise
    points between the last angle and 360 deg get ``deltaL = 0``.

    ``deltaL`` uses the nearest value in ``(r, tau, z)``. If the result follows
    the terrain (``topo_following``) and the grid has a single ``z`` level, the
    receiver level closest to the mean height above the terrain is used instead.

    Returns:
        Dataset with ``deltaL(*grid_dims, freq, source_height)`` and ``terrain(*grid_dims)``.
    """
    ds = ensure_source_height_dim(ds)
    X, Y, Z, dims, coords = receivers.receiver_positions(grid)
    dx = X - ds.attrs["center_x"]
    dy = Y - ds.attrs["center_y"]
    new_r = np.sqrt(dx**2 + dy**2)
    new_tau = np.degrees(np.arctan2(dy, dx)) % 360
    rt = np.stack([new_r, new_tau], axis=-1)
    polar_rt = (ds.r.values, ds.tau.values)

    terrain = RegularGridInterpolator(polar_rt, ds.terrain.values, bounds_error=False, fill_value=0.0)(rt)
    surface = RegularGridInterpolator(polar_rt, ds.surface.values, bounds_error=False, fill_value=0.0)(rt)

    single_level = "z" in grid.dims and grid.sizes["z"] == 1
    if ds.attrs.get("topo_following") and single_level:
        z_rel = np.mean(Z - terrain)
        iz = int(np.argmin(np.abs(ds.z.values - z_rel)))
        dl = RegularGridInterpolator(
            polar_rt, ds.deltaL.values[:, :, iz], bounds_error=False, fill_value=0.0
        )(rt)
    else:
        dl = RegularGridInterpolator(
            (ds.r.values, ds.tau.values, ds.z.values),
            ds.deltaL.values,
            bounds_error=False,
            fill_value=0.0,
            method="nearest",
        )(np.stack([new_r, new_tau, Z], axis=-1))

    out = xr.Dataset(coords=coords).assign_coords(freq=ds.freq, source_height=ds.source_height)
    out["deltaL"] = (dims + ("freq", "source_height"), dl, {"units": "dB"})
    out["terrain"] = (dims, terrain, {"units": "m"})
    out["surface"] = (dims, surface, {"units": "m"})
    out.attrs.update(ds.attrs)
    out.attrs["grid_type"] = grid.attrs.get("grid_type", "")
    return out


def tl_to_deltaL(ds: xr.Dataset) -> xr.Dataset:
    """Convert a transmission loss field (EWAPE, ``attrs["field"] == "TL"``) to an effective deltaL.

    ``deltaL = -TL + 20 log10(Rff / 1 m) + 26 delta_water`` with ``Rff`` the distance between the source
    (at ``source_height``) and the receiver, and ``delta_water = 1`` for receivers below the surface
    ``ds.surface``. Works on the polar result and on the output of :func:`to_grid`.
    A dataset that is already a deltaL is returned unchanged.

    Convert before :func:`third_octave`: the energetic average does not commute with the sign of TL.
    """
    if ds.attrs.get("field", "deltaL") != "TL":
        return ds
    ds = ensure_source_height_dim(ds)
    X, Y, Z, grid_dims, _ = receivers.receiver_positions(ds)
    r = np.hypot(X - ds.attrs["center_x"], Y - ds.attrs["center_y"])
    hs = ds.source_height.values
    rff = np.maximum(np.sqrt(r[..., None] ** 2 + (Z[..., None] - hs) ** 2), 1.0)  # (*grid, source_height)
    water = Z < ds.surface.transpose(*grid_dims).values  # (*grid)
    corr = 20 * np.log10(rff)[..., None, :] + 26.0 * water[..., None, None]  # (*grid, 1, source_height)
    dl = ds.deltaL.transpose(*grid_dims, "freq", "source_height") + corr
    out = ds.copy()
    out["deltaL"] = dl.transpose(*ds.deltaL.dims).assign_attrs(units="dB")
    out.attrs.update(field="deltaL", converted_from_TL=1)
    return out


DEFAULT_FC = [50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000]
DEFAULT_NFC = [1, 1, 1, 1, 1, 1, 2, 2, 3, 4, 4, 4, 5, 5]


def third_octave(ds: xr.Dataset, fc=None, nfc=None) -> xr.Dataset:
    """Energetic average of ``deltaL`` over the frequencies of each third-octave band.

    ``ds.freq`` must be exactly :func:`wfn.units.third_octave_frequencies` ``(fc, nfc)``.
    The result has ``freq = fc``.
    """
    fc = DEFAULT_FC if fc is None else list(fc)
    nfc = DEFAULT_NFC if nfc is None else list(nfc)
    expected = third_octave_frequencies(fc, nfc)
    if ds.sizes["freq"] != expected.size or not np.all(ds.freq.values == expected):
        raise ValueError(f"frequencies {ds.freq.values} do not match the bands {expected}")
    band = np.repeat(np.asarray(fc, dtype=float), nfc)
    lin = 10 ** (ds.deltaL / 10)
    dl = 10 * np.log10(lin.assign_coords(freq=band).groupby("freq").mean())
    out = ds.drop_dims("freq")
    out["deltaL"] = dl.transpose(*ds.deltaL.dims)
    return out


def ensure_source_height_dim(ds: xr.Dataset) -> xr.Dataset:
    """Make ``source_height`` a dimension of ``deltaL`` (also when it is a scalar coordinate)."""
    if "deltaL" not in ds or "source_height" in ds.deltaL.dims:
        return ds
    if "source_height" not in ds.coords:
        ds = ds.assign_coords(source_height=[float(ds.attrs.get("source_height", np.nan))])
        return ds.assign(deltaL=ds.deltaL.expand_dims("source_height", axis=-1))
    h = np.atleast_1d(ds.source_height.values)
    ds = ds.drop_vars("source_height")
    ds = ds.assign(deltaL=ds.deltaL.expand_dims("source_height", axis=-1))
    return ds.assign_coords(source_height=("source_height", h, {"units": "m"}))


class _DatasetWrapper:
    """Holds an ``xr.Dataset`` in ``ds``; unknown attributes and ``[]`` are forwarded to it."""

    def __init__(self, ds: xr.Dataset):
        self.ds = ds

    def __getattr__(self, name: str):
        if name == "ds":  # not set yet (copy, pickle)
            raise AttributeError(name)
        return getattr(self.ds, name)

    def __getitem__(self, key):
        return self.ds[key]

    def __repr__(self) -> str:
        return f"{type(self).__name__}\n{self.ds!r}"

    def _new(self, ds: xr.Dataset):
        return type(self)(ds)

    def third_octave(self, fc=None, nfc=None):
        """See :func:`third_octave`."""
        return self._new(third_octave(self.ds, fc, nfc))


class PeResults(_DatasetWrapper):
    """Result of one WAPE run (one source height, one angle): ``deltaL(r, z, freq)``."""

    @classmethod
    def from_file(cls, fname: str | Path, side: bool = True) -> "PeResults":
        """See :func:`read_wape`."""
        return cls(read_wape(fname, side))

    def plot_side(self, freq: float, cmap: str = "RdBu_r", **kwargs):
        """Vertical plane of ``deltaL`` at the frequency closest to ``freq``."""
        da = self.ds.deltaL.sel(freq=freq, method="nearest")
        return da.assign_coords(z_abs=self.ds.z_abs).plot.pcolormesh(x="r", y="z", cmap=cmap, **kwargs)

    def plot_line(self, freq: float, z: float, ax=None, **kwargs):
        """``deltaL(r)`` at the height ``z`` and the frequency closest to ``freq``."""
        da = self.ds.deltaL.sel(freq=freq, z=z, method="nearest")
        return da.plot.line(x="r", ax=ax, **kwargs)


class DeltaL(_DatasetWrapper):
    """``deltaL(r, tau, z, freq, source_height)`` of one source on a polar grid (old ``DeltaLField``)."""
    def __init__(self, ds: xr.Dataset, cfg: PropagationConfig | None = None):
        super().__init__(ensure_source_height_dim(ds))
        self.cfg = cfg

    @classmethod
    def from_cases(cls, cfg: PropagationConfig, isource: int = 0) -> "DeltaL":
        return cls(assemble_polar(cfg, isource), cfg)

    def _new(self, ds):                      # keep cfg through loop_angle, to_grid, ...
        return type(self)(ds, self.cfg)    


    @classmethod
    def from_yaml(cls, fname: str, isource: int = 0) -> "DeltaL":
        """See :func:`assemble_polar`."""
        cfg = PropagationConfig.from_yaml(fname)
        ds = assemble_polar(cfg, isource)
        ds.attrs["yaml_fname"] = fname
        return cls(ds,cfg)

    def save(self, path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        ds = self.ds.copy()                  # shallow copy: do not touch self.ds.attrs
        if self.cfg is not None:
            ds.attrs["cfg"] = json.dumps(unstructure(self.cfg))
        ds.to_netcdf(path)
        return path

    @classmethod
    def load(cls, path) -> "DeltaL":
        with xr.open_dataset(path) as ds:
            ds = ds.load()
        cfg = None
        if "cfg" in ds.attrs:
            data = json.loads(ds.attrs.pop("cfg"))
            base_dir = data.pop("base_dir", None)
            cfg = PropagationConfig.from_dict(data, base_dir=base_dir)
        return cls(ds, cfg)

    @property
    def is_tl(self) -> bool:
        """True if the field is still a transmission loss (EWAPE), see :meth:`to_deltaL`."""
        return self.ds.attrs.get("field") == "TL"

    def to_deltaL(self) -> "DeltaL":
        """See :func:`tl_to_deltaL`."""
        return self._new(tl_to_deltaL(self.ds))

    def loop_angle(self) -> "DeltaL":
        """See :func:`loop_angle`."""
        return self._new(loop_angle(self.ds))

    def refine_angle(self) -> "DeltaL":
        """See :func:`refine_angle`."""
        return self._new(refine_angle(self.ds))

    def to_grid(self, grid: xr.Dataset=None) -> "DeltaL":
        """See :func:`to_grid`. The result is on the receiver grid, not polar."""
        if grid is None and self.ds.yaml_fname is not None:
            print(f"build receiver from {self.ds.yaml_fname}")
            grid = ReceiversConfig.from_yaml(self.ds.yaml_fname).build()
        return self._new(to_grid(self.ds, grid))
