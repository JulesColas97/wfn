"""Receiver grids, stored as ``xarray.Dataset`` objects (replaces the old ``Mesh`` class).

A receiver grid is an ``xr.Dataset`` with no data variables, only coordinates.
Whatever the layout, the coordinates ``x``, ``y`` and ``z`` are always defined,
so that ``xr.broadcast(ds.x, ds.y, ds.z)`` gives the full 3D receiver positions:

* cartesian: dims ``(x, y, z)``, 1D coordinates ``x``, ``y``, ``z``
* polar:     dims ``(r, tau, z)``, ``x(r, tau)``, ``y(r, tau)``, ``z(z)``; ``tau`` in degrees
* points:    dims ``(point,)``, ``x(point)``, ``y(point)``, ``z(point)``

Results computed on a grid keep the grid dims, so they can be selected,
interpolated (``ds.interp``) and saved (``ds.to_netcdf``) directly with xarray.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
import xarray as xr

from .config import get_section, structure


def cartesian(x: Sequence[float], y: Sequence[float], z: Sequence[float]) -> xr.Dataset:
    """Cartesian receiver grid with dims ``(x, y, z)``."""
    ds = xr.Dataset(
        coords={
            "x": ("x", np.atleast_1d(np.asarray(x, dtype=float)), {"units": "m"}),
            "y": ("y", np.atleast_1d(np.asarray(y, dtype=float)), {"units": "m"}),
            "z": ("z", np.atleast_1d(np.asarray(z, dtype=float)), {"units": "m"}),
        }
    )
    ds.attrs["grid_type"] = "cartesian"
    return ds


def polar(
    r: Sequence[float],
    tau: Sequence[float],
    z: Sequence[float],
    center: tuple[float, float] = (0.0, 0.0),
) -> xr.Dataset:
    """Polar receiver grid with dims ``(r, tau, z)`` centred on ``center``.

    ``tau`` is the angle in degrees measured from the x axis (anticlockwise).
    """
    r = np.atleast_1d(np.asarray(r, dtype=float))
    tau = np.atleast_1d(np.asarray(tau, dtype=float))
    tau_rad = np.deg2rad(tau)
    x2d = r[:, None] * np.cos(tau_rad)[None, :] + center[0]
    y2d = r[:, None] * np.sin(tau_rad)[None, :] + center[1]
    ds = xr.Dataset(
        coords={
            "r": ("r", r, {"units": "m"}),
            "tau": ("tau", tau, {"units": "deg"}),
            "z": ("z", np.atleast_1d(np.asarray(z, dtype=float)), {"units": "m"}),
            "x": (("r", "tau"), x2d, {"units": "m"}),
            "y": (("r", "tau"), y2d, {"units": "m"}),
        }
    )
    ds.attrs.update(grid_type="polar", center_x=float(center[0]), center_y=float(center[1]))
    return ds


def points(x: Sequence[float], y: Sequence[float], z: Sequence[float]) -> xr.Dataset:
    """Unstructured list of receivers with dim ``point``."""
    x, y, z = (np.atleast_1d(np.asarray(v, dtype=float)) for v in (x, y, z))
    if not (x.shape == y.shape == z.shape):
        raise ValueError("x, y and z must have the same length")
    ds = xr.Dataset(coords={"x": ("point", x), "y": ("point", y), "z": ("point", z)})
    ds.attrs["grid_type"] = "points"
    return ds


def receiver_positions(grid: xr.Dataset) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple, dict]:
    """Broadcast the grid coordinates.

    Returns:
        X, Y, Z arrays with the full grid shape, the grid dims and the grid
        coordinates (to rebuild a ``DataArray`` with the same layout).
    """
    X, Y, Z = xr.broadcast(grid.x, grid.y, grid.z)
    dims = X.dims
    Y = Y.transpose(*dims)
    Z = Z.transpose(*dims)
    return X.values, Y.values, Z.values, dims, dict(grid.coords)


def shift(grid: xr.Dataset, dx: float, dy: float) -> xr.Dataset:
    """Return a copy of the grid translated by ``(dx, dy)``."""
    out = grid.assign_coords(x=grid.x + dx, y=grid.y + dy)
    if "center_x" in out.attrs:
        out.attrs["center_x"] += dx
        out.attrs["center_y"] += dy
    return out


# -----------------------------------------------------------------------------
# YAML section ``receivers``
# -----------------------------------------------------------------------------
@dataclass
class AxisConfig:
    """1D axis, either given explicitly (``values``) or as ``min/max/n``."""

    values: Optional[list[float]] = None
    min: Optional[float] = None
    max: Optional[float] = None
    n: Optional[int] = None

    def array(self) -> np.ndarray:
        if self.values is not None:
            return np.asarray(self.values, dtype=float)
        if None in (self.min, self.max, self.n):
            raise ValueError("axis needs either 'values' or 'min', 'max' and 'n'")
        return np.linspace(self.min, self.max, self.n)


def _axis(value) -> AxisConfig:
    if isinstance(value, AxisConfig):
        return value
    if isinstance(value, dict):
        return structure(value, AxisConfig)
    return AxisConfig(values=list(np.atleast_1d(value).astype(float)))


@dataclass
class ReceiversConfig:
    """Content of the ``receivers`` section of the input file.

    Example::

        receivers:
          type: cartesian            # cartesian | polar | points
          x: {min: -200, max: 200, n: 70}
          y: {min: -200, max: 200, n: 70}
          z: [2]
          relative_to_source: false  # if true the grid is translated to each turbine hub
    """

    type: str = "cartesian"
    x: object = None
    y: object = None
    z: object = field(default_factory=lambda: [2.0])
    r: object = None
    tau: object = None
    relative_to_source: bool = False

    @classmethod
    def from_dict(cls, data: dict) -> "ReceiversConfig":
        return cls(**data)

    @classmethod
    def from_yaml(cls, path) -> "ReceiversConfig":
        return cls.from_dict(get_section(path, "receivers"))

    def build(self, center: tuple[float, float] = (0.0, 0.0)) -> xr.Dataset:
        """Build the receiver grid. ``center`` is only used if ``relative_to_source``
        is set or for polar grids."""
        c = center if (self.relative_to_source or self.type == "polar") else (0.0, 0.0)
        z = _axis(self.z).array()
        if self.type == "cartesian":
            return cartesian(_axis(self.x).array() + c[0], _axis(self.y).array() + c[1], z)
        if self.type == "polar":
            return polar(_axis(self.r).array(), _axis(self.tau).array(), z, center=c)
        if self.type == "points":
            return points(_axis(self.x).array() + c[0], _axis(self.y).array() + c[1], z)
        raise ValueError(f"unknown receiver grid type '{self.type}'")
