"""Abstract flow model.

A flow model provides the atmospheric fields to the other modules:

* the **source model** only needs point values of the wind (``u``, ``v``) and of
  the turbulent dissipation rate (``epsilon``) at the blade segment positions,
  plus optionally the turbine positions;
* the **propagation model** (later) will need the full fields on a regular
  grid (``u``, ``v``, ``w``, sound speed ``c``...), see :meth:`Flow.to_dataset`.

New flow types (Walberla, Meso-NH, FarmShadow, ...) subclass :class:`Flow` and
register themselves with a ``type`` name used in the YAML ``flow`` section::

    class FlowWalberla(Flow, type_name="walberla"):
        ...
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import ClassVar

import numpy as np
import xarray as xr

from ..config import get_section, structure


class Flow(ABC):
    """Interface of all flow models.

    Attributes expected on every implementation:
        turbine_positions: ``(N, 3)`` array of hub positions ``(x, y, z_hub)``
            of the turbines contained in the flow, or ``None`` if the flow does
            not define turbines.
    """

    _registry: ClassVar[dict[str, type["Flow"]]] = {}
    type_name: ClassVar[str | None] = None

    def __init_subclass__(cls, type_name: str | None = None, **kwargs):
        super().__init_subclass__(**kwargs)
        if type_name is not None:
            cls.type_name = type_name
            Flow._registry[type_name] = cls

    # ------------------------------------------------------------------ data
    @property
    def turbine_positions(self) -> np.ndarray | None:
        """Hub positions ``(N, 3)`` of the turbines defined by the flow, if any."""
        return None

    @abstractmethod
    def at_points(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict[str, np.ndarray]:
        """Flow quantities at arbitrary points.

        Args:
            x, y, z: arrays of identical shape (m).

        Returns:
            dict with at least the keys ``"u"``, ``"v"`` (m/s) and ``"epsilon"``
            (m^2/s^3), each with the same shape as ``x``.
        """

    def to_dataset(self) -> xr.Dataset:
        """Full fields on a regular grid, dims ``(x, y, z)``.

        Expected variables: ``u``, ``v``, ``w``, ``c`` (or ``temperature``),
        ``epsilon``. Needed by the propagation model only.
        """
        raise NotImplementedError(f"{type(self).__name__} does not provide gridded fields")

    # ------------------------------------------------------------- factories
    @classmethod
    def from_dict(cls, data: dict) -> "Flow":
        """Build the flow model registered under ``data["type"]``."""
        data = dict(data)
        type_name = data.pop("type", "uniform")
        if type_name not in Flow._registry:
            raise ValueError(
                f"unknown flow type '{type_name}', available: {sorted(Flow._registry)}"
            )
        return structure(data, Flow._registry[type_name])

    @classmethod
    def from_yaml(cls, path) -> "Flow":
        """Build the flow from the ``flow`` section of the input file."""
        return cls.from_dict(get_section(path, "flow"))
