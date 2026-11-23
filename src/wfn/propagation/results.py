"""Read WAPE results and assemble them into ``xr.Dataset`` objects.

* :func:`read_wape`: one WAPE output file (one source height, one angle) with
  ``deltaL(r, z, freq)``.
* :func:`assemble_polar`: all the angles and source heights of one source, with
  ``deltaL(r, tau, z, freq, source_height)`` and ``terrain(r, tau)``.
* :func:`to_grid`: interpolate a polar result onto any receiver grid.
* :func:`third_octave`: average the frequencies into third-octave bands.
"""
from __future__ import annotations

import logging
from pathlib import Path

import h5py
import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator

from .. import receivers
from ..units import third_octave_frequencies
from .cases import case_dir, local_domain
from .config import PropagationConfig

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
        Dataset with ``deltaL(r, z, freq)``, the ground height ``h(r)`` and the
        absolute height ``z_abs(r, z)`` of every point.
    """
    with h5py.File(fname, "r") as f:
        mesh = f["mesh"]
        r = np.array(mesh["x"])
        h = np.array(mesh["h"])
        group = f["solution"] if side else f["receiver"]
        keys = [k for k in group.keys() if k != "heights"]
        if side:
            z = np.array(mesh["z"])
            z_abs = np.array(mesh["z_grid"]).T
        else:
            z = np.array(group["heights"])
            z_abs = z[None, :] + h[:, None]
        dl = np.stack([np.array(group[k]["deltaL"]).T for k in keys], axis=-1)

    order = np.argsort([float(k) for k in keys])
    freq = np.array([float(k) for k in keys])[order]
    dl = dl[..., order]
    dl[np.abs(dl) > DL_MAX] = np.nan
    return xr.Dataset(
        {
            "deltaL": (("r", "z", "freq"), dl, {"units": "dB"}),
            "h": ("r", h, {"units": "m"}),
            "z_abs": (("r", "z"), z_abs, {"units": "m"}),
        },
        coords={"r": r, "z": z, "freq": freq},
        attrs={"side": int(side)},
    )


def assemble_polar(cfg: PropagationConfig, isource: int = 0) -> xr.Dataset:
    """Gather all WAPE runs of one source on a common polar grid centred on the source.

    The grid has ``r = 0, dout, ... rmax`` with ``rmax`` the largest distance to
    a corner of the domain plus 50 m, and ``z = 0, dout, ..., < Lz`` (``side``) or
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
    dr = w.dout
    rmax = np.sqrt(max(dom["Lx1"], dom["Lx2"]) ** 2 + max(dom["Ly1"], dom["Ly2"]) ** 2) + R_MARGIN
    r = np.arange(0, rmax + dr, dr)
    z = np.arange(0, dom["Lz"], dr) if side else np.asarray(w.receiver_heights, dtype=float)
    tau = np.asarray(cfg.angles(), dtype=float)
    heights = np.asarray(cfg.heights(), dtype=float)
    freq = np.asarray(cfg.freqs(), dtype=float)

    dl = np.zeros((r.size, tau.size, z.size, freq.size, heights.size))
    terrain = np.full((r.size, tau.size), np.nan)
    for ih, height in enumerate(heights):
        for it, t in enumerate(tau):
            res = read_wape(wape_filename(case_dir(cfg, isource, height), w.case_name, t), side)
            if res.sizes["z"] != z.size:
                raise ValueError(
                    f"WAPE run (height {height}, tau {t}) has {res.sizes['z']} heights, expected {z.size}"
                )
            res_dl = res.deltaL.sel(freq=freq, method="nearest").values
            i0 = int(np.argmin(np.abs(r - res.r.values[0])))
            i1 = int(np.argmin(np.abs(r - res.r.values[-1])))
            dl[i0 : i1 + 1, it, :, :, ih] = res_dl
            terrain[i0 : i1 + 1, it] = res.h.values
            # copy the last values up to rmax
            last = i0 + np.where(res_dl > -180)[0][-1]
            dl[last:, it, :, :, ih] = dl[last - 1, it, :, :, ih]
            terrain[last:, it] = terrain[last - 1, it]

    xs, ys = cfg.source_positions()[isource][:2]
    ds = receivers.polar(r, tau, z, center=(xs, ys))
    ds = ds.assign_coords(freq=("freq", freq, {"units": "Hz"}), source_height=("source_height", heights, {"units": "m"}))
    ds["deltaL"] = (("r", "tau", "z", "freq", "source_height"), dl, {"units": "dB"})
    ds["terrain"] = (("r", "tau"), terrain, {"units": "m"})
    topo = bool(np.nanmax(terrain) - np.nanmin(terrain) > 1e-6)
    ds.attrs.update(side=int(side), topo=int(topo), topo_following=int(topo and not side), isource=isource)
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
    X, Y, Z, dims, coords = receivers.receiver_positions(grid)
    dx = X - ds.attrs["center_x"]
    dy = Y - ds.attrs["center_y"]
    new_r = np.sqrt(dx**2 + dy**2)
    new_tau = np.degrees(np.arctan2(dy, dx)) % 360
    rt = np.stack([new_r, new_tau], axis=-1)
    polar_rt = (ds.r.values, ds.tau.values)

    terrain = RegularGridInterpolator(polar_rt, ds.terrain.values, bounds_error=False, fill_value=0.0)(rt)

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
    out.attrs.update(ds.attrs)
    out.attrs["grid_type"] = grid.attrs.get("grid_type", "")
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
