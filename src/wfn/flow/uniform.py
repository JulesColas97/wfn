"""Simplest flow: uniform (or power-law sheared) wind and constant dissipation."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .base import Flow


@dataclass
class UniformFlow(Flow, type_name="uniform"):
    """Analytical flow, mainly for tests and single turbine studies.

    ``u(z) = u * (z / z_ref) ** shear_alpha`` if ``shear_alpha`` is given, else ``u``.

    YAML::

        flow:
          type: uniform
          u: 8.0
          epsilon: 0.01
          shear_alpha: 0.2   # optional
          z_ref: 90          # needed with shear_alpha
          wind_direction: 0  # deg, direction the wind blows towards, from the x axis

    Attributes:
        u: wind speed (m/s) at ``z_ref``.
        epsilon: turbulent dissipation rate (m^2/s^3).
        shear_alpha: power-law exponent (None for a uniform profile).
        z_ref: reference height (m) of the power law.
        wind_direction: direction of the wind (deg) measured from the x axis.
        positions: optional list of turbine hub positions ``[[x, y, z], ...]``.
    """

    u: float = 8.0
    epsilon: float = 0.01
    shear_alpha: Optional[float] = None
    z_ref: Optional[float] = None
    wind_direction: float = 0.0
    positions: Optional[list[list[float]]] = None

    def __post_init__(self):
        if self.shear_alpha is not None and self.z_ref is None:
            raise ValueError("z_ref is required when shear_alpha is set")

    @property
    def turbine_positions(self) -> np.ndarray | None:
        return None if self.positions is None else np.asarray(self.positions, dtype=float)

    def speed(self, z: np.ndarray) -> np.ndarray:
        z = np.asarray(z, dtype=float)
        if self.shear_alpha is None:
            return np.full_like(z, self.u)
        return self.u * (np.maximum(z, 0.0) / self.z_ref) ** self.shear_alpha

    def at_points(self, x, y, z) -> dict[str, np.ndarray]:
        speed = self.speed(z)
        wd = np.deg2rad(self.wind_direction)
        return {
            "u": speed * np.cos(wd),
            "v": speed * np.sin(wd),
            "epsilon": np.full_like(speed, self.epsilon),
        }
