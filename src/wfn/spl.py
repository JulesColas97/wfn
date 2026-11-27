"""Coupling of the source model and the propagation model (old ``SplField.combine_linear``).

Typical use::

    from wfn import spl, units
    from wfn.propagation import DeltaL

    src = xr.open_dataset("results/nrel_0.nc")        # WindTurbineNoise.compute
    dl = DeltaL.load("results/deltaL_0.nc")            # deltaL(r, tau, z, freq, source_height)
    out = spl.combine(src, dl)                         # on the source grid
    level = units.db_sum(units.spl(out), ("blade", "seg"))

:func:`combine` multiplies the segment PSDs by ``10**(deltaL / 10)``, where ``deltaL``
is interpolated (in energy) between the propagation source heights at the
height ``z_seg`` of each blade segment. The output keeps the source layout
(linear PSDs on ``(*grid, beta, blade, seg, freq)``), so all :mod:`wfn.units`
helpers work on it as they do on the free-field result.

Difference with the old code: ``combine_linear`` swapped the two interpolation
weights (``ratio * dL_left + (1 - ratio) * dL_right``); here the interpolation
is the usual linear one.
"""
from __future__ import annotations

import numpy as np
import xarray as xr

from . import receivers as rcv
from .propagation.results import DeltaL, ensure_source_height_dim, loop_angle, to_grid


def _grid_dims(ds: xr.Dataset) -> list[str]:
    return list(dict.fromkeys((*ds.x.dims, *ds.y.dims, *ds.z.dims)))


def grid_of(ds: xr.Dataset) -> xr.Dataset:
    """Receiver grid (coordinates only) of a result dataset."""
    dims = _grid_dims(ds)
    grid = ds.drop_vars(list(ds.data_vars)).drop_dims([d for d in ds.dims if d not in dims])
    return grid.assign_attrs(ds.attrs)


def same_grid(a: xr.Dataset, b: xr.Dataset) -> bool:
    """True if ``a`` and ``b`` have the same grid dims and the same receiver positions."""
    Xa, Ya, Za, da, _ = rcv.receiver_positions(a)
    Xb, Yb, Zb, db, _ = rcv.receiver_positions(b)
    if da != db or Xa.shape != Xb.shape:
        return False
    return all(np.allclose(p, q) for p, q in ((Xa, Xb), (Ya, Yb), (Za, Zb)))


def _full_circle(tau: np.ndarray) -> bool:
    """True if the angles are evenly spaced over 360 deg without the last one (0, 90, 180, 270)."""
    if tau.size < 2:
        return False
    return bool(np.isclose(tau[0] + (tau[1] - tau[0]) * tau.size, tau[0] + 360.0))


def _is_cartesian(ds: xr.Dataset) -> bool:
    return all(d in ds.dims for d in ("x", "y", "z")) and ds.x.dims == ("x",)


def _interp_cartesian(ds: xr.Dataset, grid: xr.Dataset) -> xr.Dataset:
    """Linear interpolation cartesian -> cartesian (NaN outside the original grid)."""
    if not (_is_cartesian(ds) and _is_cartesian(grid)):
        raise ValueError(
            f"cannot interpolate a '{ds.attrs.get('grid_type')}' result onto a "
            f"'{grid.attrs.get('grid_type')}' grid (only cartesian -> cartesian)"
        )
    new = {}
    for d in ("x", "y", "z"):
        old_v, new_v = ds[d].values, grid[d].values
        if old_v.shape == new_v.shape and np.allclose(old_v, new_v):
            continue
        if old_v.size == 1:  # cannot interpolate along a single value
            raise ValueError(f"cannot interpolate along '{d}': the original grid has a single value {old_v}")
        new[d] = new_v
    return ds.interp(**new) if new else ds


def _deltaL_on_grid(dl: xr.Dataset, grid: xr.Dataset) -> xr.Dataset:
    if same_grid(dl, grid):
        return dl
    if dl.attrs.get("grid_type") == "polar" and "center_x" in dl.attrs:
        if _full_circle(dl.tau.values):
            dl = loop_angle(dl)
        return to_grid(dl, grid)
    return _interp_cartesian(dl, grid)


def _source_on_grid(src: xr.Dataset, grid: xr.Dataset) -> xr.Dataset:
    if same_grid(src, grid):
        return src
    out = _interp_cartesian(src, grid)
    if "distance" in out and "hub_x" in src.attrs:
        X, Y, Z, dims, _ = rcv.receiver_positions(out)
        d = np.sqrt((X - src.hub_x) ** 2 + (Y - src.hub_y) ** 2 + (Z - src.hub_z) ** 2)
        out["distance"] = (dims, d, src.distance.attrs)
    return out.assign_attrs(grid_type=grid.attrs.get("grid_type", "cartesian"))


def combine(src: xr.Dataset, dl: xr.Dataset | DeltaL, grid: xr.Dataset | None = None,
            keep_deltaL: bool = False) -> xr.Dataset:
    """Add the propagation effects to a source result.

    Args:
        src: output of :meth:`wfn.source.WindTurbineNoise.compute` (one turbine).
        dl: propagation result of the same turbine, either polar
            (:func:`wfn.propagation.assemble_polar`, :class:`DeltaL`) or already
            on a receiver grid (:func:`wfn.propagation.to_grid`).
        grid: common receiver grid. If None (default, as in the old code),
            ``deltaL`` is interpolated onto the source grid. Otherwise both are
            interpolated onto ``grid`` (the source only cartesian -> cartesian).
        keep_deltaL: also store ``deltaL`` of every segment (same size as ``Spp_TEN``).

    Returns:
        The source dataset on the common grid with ``Spp_TEN`` and ``Spp_TIN``
        multiplied by ``10**(deltaL(z_seg) / 10)``, plus ``deltaL``
        (dB, ``(*grid, beta, blade, seg, freq)``) and ``terrain`` if available.
    """
    if isinstance(dl, DeltaL):
        dl = dl.ds
    dl = ensure_source_height_dim(dl)
    if dl.attrs.get("field") == "TL":
        raise ValueError("deltaL is still a transmission loss (TL): call DeltaL.to_deltaL() first")
    if "isource" in src.attrs and "isource" in dl.attrs and int(src.isource) != int(dl.isource):
        raise ValueError(f"source isource={src.isource} but deltaL isource={dl.isource}")

    # 1. common grid
    grid = grid_of(src) if grid is None else grid
    src = _source_on_grid(src, grid)
    dl = _deltaL_on_grid(dl, grid)

    # 2. same frequencies (deltaL may contain more)
    try:
        dl = dl.sel(freq=src.freq.values, method="nearest", tolerance=1e-6)
    except KeyError:
        missing = np.setdiff1d(src.freq.values, dl.freq.values)
        raise ValueError(f"frequencies {missing} of the source are not in deltaL") from None
    dl = dl.assign_coords(freq=src.freq)

    # 3. deltaL at the height of each segment, interpolated in energy between the source heights
    lin = 10 ** (dl.deltaL / 10)
    h = dl.source_height
    if h.size == 1:
        lin = lin.isel(source_height=0, drop=True)
    else:
        z = src.z_seg.clip(float(h.min()), float(h.max()))
        lin = lin.interp(source_height=z).drop_vars("source_height")
    lin = xr.broadcast(lin, src.Spp_TEN)[0].transpose(*src.Spp_TEN.dims)  # single height: no (beta, blade, seg) dims yet

    # 4. output on the source layout
    out = src.copy()
    out["Spp_TEN"] = (src.Spp_TEN * lin).assign_attrs(src.Spp_TEN.attrs)
    out["Spp_TIN"] = (src.Spp_TIN * lin).assign_attrs(src.Spp_TIN.attrs)
    if keep_deltaL:
        out["deltaL"] = (10 * np.log10(lin)).assign_attrs(units="dB", long_name="deltaL at each segment")
    for name in ("terrain", "surface"):
        if name in dl:
            out[name] = dl[name]
    out.attrs.update(coupled=1, deltaL_source_heights=h.values,
                     deltaL_topo=int(dl.attrs.get("topo", 0)), deltaL_side=int(dl.attrs.get("side", 0)))
    return out
