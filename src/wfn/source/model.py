"""Wind turbine noise source model (TEN + TIN), results as ``xarray.Dataset``.

Typical use::

    from wfn.source import WindTurbineNoise
    from wfn import receivers

    src = WindTurbineNoise.from_yaml("input.yaml")
    ds = src.compute(receivers.cartesian(x, y, [2.0]))
    ds.to_netcdf("source.nc")

Output dataset (see :meth:`WindTurbineNoise.compute`):

* ``Spp_TEN``, ``Spp_TIN``: one-sided pressure PSD per rad/s of each segment,
  dims ``(*grid_dims, beta, blade, seg, freq)``. PSDs are linear, so the
  rotor total at one rotor position is ``(Spp_TEN + Spp_TIN).sum(("blade", "seg"))``
  and :func:`wfn.units.spl` converts to dB;
* ``distance``: hub-receiver distance (grid dims), used by :func:`wfn.units.swl`;
* rotor variables on ``(beta, blade, seg)``: ``x_seg, y_seg, z_seg, U_inf,
  U_rel, AoA, epsilon``. ``z_seg`` is what the coupling with the propagation
  model (ΔL depending on source height) will use.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

import numpy as np
import xarray as xr
from joblib import Parallel, delayed

from .. import receivers as rcv
from ..config import get_section, load_yaml
from ..flow import Flow, UniformFlow
from .airfoil import BL_VARS, AirfoilDatabase
from .amiet import SegmentState, kolmogorov_term, segment_spectra
from .config import AtmosphereConfig, SourceConfig
from .geometry import blade_frame_transforms, receiver_angles, segment_positions, to_segment_frame
from .turbine import WindTurbine

logger = logging.getLogger(__name__)


def _one_azimuth(k: int, azimuth: float, rec: np.ndarray, hub: np.ndarray, yaw: float,
                 R: np.ndarray, T: np.ndarray, rotor: dict, freq: np.ndarray,
                 c0: float, rho0: float, nu0: float) -> tuple[np.ndarray, np.ndarray]:
    """TEN and TIN of all segments at azimuth index ``k``: arrays ``(N, S, F)``."""
    n_seg = R.shape[0]
    ten = np.empty((rec.shape[0], n_seg, freq.size))
    tin = np.empty_like(ten)
    _, teta, psi = receiver_angles(rec, hub, azimuth, yaw)
    for s in range(n_seg):
        xyz = to_segment_frame(rec, R[s], T[s])
        st = SegmentState(
            chord=rotor["chord"][s],
            span=rotor["span"][s],
            u_rot=rotor["U_rot"][s],
            u_inf=rotor["U_inf"][k, s],
            kol_term=rotor["kol_term"][k, s],
            bl={v: rotor[v][k, s] for v in BL_VARS},
        )
        ten[:, s, :], tin[:, s, :] = segment_spectra(xyz, teta, psi, freq, st, c0, rho0, nu0, iseg=s)
    return ten, tin


class WindTurbineNoise:
    """Compute TEN and TIN of a wind turbine at a set of receivers.

    Args:
        turbine: blade geometry and operating parameters.
        airfoil: boundary layer database matching the turbine segments.
        frequencies: frequencies (Hz).
        flow: flow model giving ``u, v, epsilon`` at the blade segments
            (default: uniform 8 m/s, epsilon 0.01).
        atmosphere: constant sound speed / density / viscosity.
        positions: hub positions ``(N, 3)``; default: from the flow, else ``(0, 0, hub_height)``.
        n_jobs: number of processes (parallel over azimuths).
    """

    def __init__(self, turbine: WindTurbine, airfoil: AirfoilDatabase, frequencies: Sequence[float],
                 flow: Flow | None = None, atmosphere: AtmosphereConfig | None = None,
                 positions: np.ndarray | None = None, n_jobs: int = 1):
        turbine.check_acoustic_geometry()
        airfoil.check_turbine(turbine)
        self.turbine = turbine
        self.airfoil = airfoil
        self.frequencies = np.atleast_1d(np.asarray(frequencies, dtype=float))
        self.flow = flow if flow is not None else UniformFlow()
        self.atmosphere = atmosphere if atmosphere is not None else AtmosphereConfig()
        self._positions = None if positions is None else np.atleast_2d(np.asarray(positions, float))
        self.n_jobs = n_jobs
        self.output: Path | None = None

    # ---------------------------------------------------------------- build
    @classmethod
    def from_config(cls, cfg: SourceConfig, flow: Flow | None = None) -> "WindTurbineNoise":
        turbine = WindTurbine.load(cfg.turbine, cfg.base_dir)
        if cfg.n_beta is not None:
            turbine.n_beta = cfg.n_beta
        db_path = cfg.path(cfg.airfoil_database) or (
            Path(turbine.airfoil_database) if turbine.airfoil_database else None)
        if db_path is None:
            raise ValueError("no airfoil database given (source.airfoil_database)")
        model = cls(turbine, AirfoilDatabase.open(db_path), cfg.frequencies, flow=flow,
                    atmosphere=cfg.atmosphere, positions=cfg.positions, n_jobs=cfg.n_jobs)
        model.output = cfg.path(cfg.output)
        return model

    @classmethod
    def from_yaml(cls, path: str | Path) -> "WindTurbineNoise":
        """Build from the ``source`` (and, if present, ``flow``) sections of the input file."""
        data = load_yaml(path)
        cfg = SourceConfig.from_dict(get_section(data, "source"), base_dir=Path(path).parent)
        flow_section = get_section(data, "flow", required=False)
        flow = Flow.from_dict(flow_section) if flow_section is not None else None
        return cls.from_config(cfg, flow)

    # ------------------------------------------------------------ positions
    @property
    def hub_positions(self) -> np.ndarray:
        """``(N, 3)`` hub positions of all turbines to compute."""
        if self._positions is not None:
            return self._positions
        if self.flow.turbine_positions is not None:
            return np.atleast_2d(self.flow.turbine_positions)
        if self.turbine.hub_height is None:
            raise ValueError("no turbine position: give source.positions or turbine.hub_height")
        return np.array([[0.0, 0.0, self.turbine.hub_height]])

    # ---------------------------------------------------------------- rotor
    def rotor_state(self, hub: Sequence[float]) -> xr.Dataset:
        """Inflow, operating point and BL quantities of every segment at every azimuth.

        Returns:
            Dataset with dims ``(azimuth, seg)`` (``azimuth`` covers a full turn,
            index ``iblade * n_beta + ibeta``). The operating point (omega,
            pitch) is in the attributes.
        """
        hub = np.asarray(hub, dtype=float)
        wt = self.turbine
        az = wt.azimuths
        xs, ys, zs = segment_positions(az, wt.seg, hub, wt.yaw_rad)
        f = self.flow.at_points(xs, ys, zs)
        u_n = f["u"] * np.cos(wt.yaw_rad) + f["v"] * np.sin(wt.yaw_rad)

        wt = wt.with_operating_point(float(np.mean(u_n)))
        u_rot = wt.seg * wt.omega
        u_rel = np.sqrt(u_n**2 + u_rot[None, :] ** 2)
        aoa = (wt.twist_cotte_conv - wt.pitch)[None, :] - np.rad2deg(np.arctan(u_rot[None, :] / u_n))
        bl = self.airfoil.boundary_layer(u_rel, aoa)

        dims = ("azimuth", "seg")
        ds = xr.Dataset(
            {
                "x_seg": (dims, xs), "y_seg": (dims, ys), "z_seg": (dims, zs),
                "U_inf": (dims, u_n), "U_rel": (dims, u_rel), "AoA": (dims, aoa),
                "epsilon": (dims, f["epsilon"] * np.ones_like(u_n)),
                "kol_term": (dims, kolmogorov_term(f["epsilon"]) * np.ones_like(u_n)),
                "U_rot": ("seg", u_rot),
                **{v: (dims, bl[v]) for v in BL_VARS},
            },
            coords={
                "azimuth": ("azimuth", az, {"units": "rad"}),
                "seg": ("seg", np.arange(wt.n_seg)),
                "r_seg": ("seg", wt.seg, {"units": "m"}),
                "chord": ("seg", wt.chord, {"units": "m"}),
                "span": ("seg", wt.span, {"units": "m"}),
            },
            attrs={
                "hub_x": hub[0], "hub_y": hub[1], "hub_z": hub[2],
                "omega_rpm": wt.omega_rpm, "pitch": wt.pitch, "yaw_deg": wt.yaw_deg,
                "n_blade": wt.n_blade, "n_beta": wt.n_beta,
            },
        )
        ds["U_inf"].attrs["description"] = "inflow velocity normal to the rotor"
        ds["AoA"].attrs["units"] = "deg"
        return ds

    # -------------------------------------------------------------- compute
    def compute(self, receivers: xr.Dataset, hub: Sequence[float] | None = None) -> xr.Dataset:
        """TEN and TIN of one turbine at the receivers of ``receivers``.

        Args:
            receivers: receiver grid (see :mod:`wfn.receivers`).
            hub: hub position; default: first of :attr:`hub_positions`.
        """
        hub = self.hub_positions[0] if hub is None else np.asarray(hub, dtype=float)
        rotor = self.rotor_state(hub)
        wt = self.turbine
        X, Y, Z, grid_dims, grid_coords = rcv.receiver_positions(receivers)
        grid_shape = X.shape
        rec = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))

        chord_angle = -wt.twist_cotte_conv + rotor.attrs["pitch"]
        R, T = blade_frame_transforms(rotor.azimuth.values, wt.seg, chord_angle, hub, wt.yaw_rad)
        rotor_np = {k: rotor[k].values for k in ["U_inf", "kol_term", "U_rot", *BL_VARS]}
        rotor_np.update(chord=wt.chord, span=wt.span)
        c0 = self.atmosphere.sound_speed
        rho0, nu0 = self.atmosphere.rho0, self.atmosphere.nu0

        logger.info("computing %d receivers x %d azimuths x %d segments x %d frequencies",
                    rec.shape[0], wt.n_azimuth, wt.n_seg, self.frequencies.size)
        args = [(k, az, rec, hub, wt.yaw_rad, R[k], T[k], rotor_np, self.frequencies, c0, rho0, nu0)
                for k, az in enumerate(rotor.azimuth.values)]
        if self.n_jobs == 1:
            results = [_one_azimuth(*a) for a in args]
        else:
            results = Parallel(n_jobs=self.n_jobs, prefer="processes")(
                delayed(_one_azimuth)(*a) for a in args)

        # (K, N, S, F) -> grid + (beta, blade, seg, freq) ; k = iblade * n_beta + ibeta
        def assemble(i):
            arr = np.stack([r[i] for r in results], axis=0)
            arr = arr.reshape(wt.n_blade, wt.n_beta, *grid_shape, wt.n_seg, self.frequencies.size)
            return np.moveaxis(arr, (0, 1), (-3, -4))

        def rotor_var(name):
            return rotor[name].values.reshape(wt.n_blade, wt.n_beta, wt.n_seg).transpose(1, 0, 2)

        out_dims = (*grid_dims, "beta", "blade", "seg", "freq")
        rdims = ("beta", "blade", "seg")
        distance = np.sqrt((X - hub[0]) ** 2 + (Y - hub[1]) ** 2 + (Z - hub[2]) ** 2)
        beta = 2 * np.pi * np.arange(wt.n_beta) / wt.n_azimuth
        azimuth = beta[:, None] + 2 * np.pi * np.arange(wt.n_blade)[None, :] / wt.n_blade

        ds = xr.Dataset(
            {
                "Spp_TEN": (out_dims, assemble(0), {"units": "Pa^2 s/rad", "long_name": "trailing edge noise PSD"}),
                "Spp_TIN": (out_dims, assemble(1), {"units": "Pa^2 s/rad", "long_name": "turbulent inflow noise PSD"}),
                "distance": (grid_dims, distance, {"units": "m", "long_name": "hub-receiver distance"}),
                **{n: (rdims, rotor_var(n), rotor[n].attrs)
                   for n in ["x_seg", "y_seg", "z_seg", "U_inf", "U_rel", "AoA", "epsilon"]},
            },
            coords={
                **{k: v for k, v in grid_coords.items()},
                "freq": ("freq", self.frequencies, {"units": "Hz"}),
                "beta": ("beta", beta, {"units": "rad", "long_name": "azimuth of blade 0"}),
                "blade": ("blade", np.arange(wt.n_blade)),
                "seg": ("seg", np.arange(wt.n_seg)),
                "r_seg": ("seg", wt.seg, {"units": "m"}),
                "azimuth": (("beta", "blade"), azimuth, {"units": "rad"}),
            },
            attrs={**rotor.attrs, **receivers.attrs, "c0": c0, "rho0": rho0, "nu0": nu0,
                   "model": "wfn.source TEN (Lee/Goody + Amiet) and TIN (Amiet, Kolmogorov)"},
        )
        return ds

    def compute_all(self, receivers: xr.Dataset | rcv.ReceiversConfig,
                    save: bool | None = None) -> list[xr.Dataset]:
        """Compute every turbine of :attr:`hub_positions`.

        Args:
            receivers: a grid, or a :class:`ReceiversConfig` built for each hub
                (useful with ``relative_to_source`` or polar grids).
            save: write ``<output>_<i>.nc`` (default: if ``output`` is set).
        """
        save = self.output is not None if save is None else save
        out = []
        for i, hub in enumerate(self.hub_positions):
            logger.info("turbine %d at %s", i, hub)
            grid = receivers.build((hub[0], hub[1])) if isinstance(receivers, rcv.ReceiversConfig) else receivers
            ds = self.compute(grid, hub).assign_attrs(isource=i)
            if save:
                if self.output is None:
                    raise ValueError("no output prefix set")
                self.output.parent.mkdir(parents=True, exist_ok=True)
                ds.to_netcdf(f"{self.output}_{i}.nc")
            out.append(ds)
        return out
