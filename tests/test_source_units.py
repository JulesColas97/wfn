import shutil
from pathlib import Path

import numpy as np
import pytest
import scipy.special
import xarray as xr

from wfn import receivers, units
from wfn.flow import Flow, UniformFlow
from wfn.source import AirfoilDatabase, SourceConfig, WindTurbine, WindTurbineNoise
from wfn.source.geometry import blade_frame_transforms, segment_positions
from wfn.source.kernels import fresnel_cs

from conftest import NREL_DB

FREQ = [100.0, 500.0, 1000.0]


@pytest.fixture(scope="module")
def nrel():
    return WindTurbine.default("NREL_5MW"), AirfoilDatabase.open(NREL_DB)


# ------------------------------------------------------------------ kernels
def test_fresnel_matches_scipy():
    x = np.linspace(0, 12, 600).reshape(3, -1)
    c, s = fresnel_cs(x)
    s_ref, c_ref = scipy.special.fresnel(x)
    np.testing.assert_allclose(c, c_ref, atol=1e-12)
    np.testing.assert_allclose(s, s_ref, atol=1e-12)


# ----------------------------------------------------------------- geometry
def test_transforms_are_rotations_and_hub_is_origin_of_rotor():
    az = np.linspace(0, 2 * np.pi, 7, endpoint=False)
    seg = np.array([10.0, 30.0])
    hub = np.array([5.0, -3.0, 90.0])
    R, T = blade_frame_transforms(az, seg, np.array([-70.0, -80.0]), hub, yaw_rad=0.3)
    eye = np.broadcast_to(np.eye(3), R.shape)
    np.testing.assert_allclose(R @ np.swapaxes(R, -1, -2), eye, atol=1e-12)
    # each segment centre maps to the origin of its own frame
    xs, ys, zs = segment_positions(az, seg, hub, yaw_rad=0.3)
    centres = np.stack([xs, ys, zs], axis=-1)
    local = np.einsum("ksij,ksj->ksi", R, centres) + T
    np.testing.assert_allclose(local, 0.0, atol=1e-9)


def test_segment_positions_convention():
    xs, ys, zs = segment_positions(np.array([0.0, np.pi / 2]), np.array([10.0]), np.zeros(3))
    np.testing.assert_allclose(zs[:, 0], [10.0, 0.0], atol=1e-12)   # beta=0: blade up
    # clockwise looking downstream (+x, z up => y to the left): blade goes to -y
    np.testing.assert_allclose(ys[:, 0], [0.0, -10.0], atol=1e-12)


def test_yaw_is_a_rigid_rotation(nrel):
    """Rotating the whole scene (rotor yaw, wind direction, receivers) around the
    hub must not change the spectra."""
    wt, db = nrel
    yaw = 30.0
    hub = (0.0, 0.0, 90.0)
    x, y, z = np.array([400.0, -150.0, 50.0]), np.array([120.0, 300.0, -500.0]), np.array([2.0, 2.0, 40.0])
    ref = WindTurbineNoise(wt, db, FREQ, flow=UniformFlow(u=9.0), positions=[hub]).compute(
        receivers.points(x, y, z))
    a = np.deg2rad(yaw)
    xr_, yr_ = np.cos(a) * x - np.sin(a) * y, np.sin(a) * x + np.cos(a) * y
    wt_yaw = WindTurbine.default("NREL_5MW")
    wt_yaw.yaw_deg = yaw
    rot = WindTurbineNoise(wt_yaw, db, FREQ, flow=UniformFlow(u=9.0, wind_direction=yaw),
                           positions=[hub]).compute(receivers.points(xr_, yr_, z))
    np.testing.assert_allclose(rot.Spp_TEN, ref.Spp_TEN, rtol=1e-9)
    np.testing.assert_allclose(rot.Spp_TIN, ref.Spp_TIN, rtol=1e-9)


# ------------------------------------------------------------------ airfoil
def test_airfoil_netcdf_roundtrip(tmp_path, nrel):
    _, db = nrel
    db.save(tmp_path / "db.nc")
    db2 = AirfoilDatabase.open(tmp_path / "db.nc")
    xr.testing.assert_identical(db.ds, db2.ds)
    u, a = np.array([[20.0] * 14]), np.array([[5.0] * 14])
    for var, val in db.boundary_layer(u, a).items():
        np.testing.assert_array_equal(val, db2.boundary_layer(u, a)[var])


def test_airfoil_interpolation_hits_grid_points(nrel):
    _, db = nrel
    iu, ia, iseg = 1, 25, 6
    u, a = db.ds.uinf.values[iu], db.ds.aoa.values[ia]
    expected = db.ds.delta_star_top.values[iu, ia, iseg]
    assert db.evaluate("delta_star_top", iseg, u, a) == pytest.approx(expected, rel=1e-10)


def test_database_must_match_turbine(nrel):
    wt, db = nrel
    short = WindTurbine(seg=wt.seg[:5], span=wt.span[:5], chord=wt.chord[:5], twist=wt.twist[:5],
                        airfoil_index=wt.airfoil_index[:5], airfoil_name=wt.airfoil_name)
    with pytest.raises(ValueError):
        WindTurbineNoise(short, db, FREQ)


# ------------------------------------------------------------------ turbine
@pytest.mark.parametrize("name", ["NREL_5MW", "IEA15MW", "SWT9.3MW"])
@pytest.mark.filterwarnings("ignore:twist values")
def test_default_turbines_load(name):
    wt = WindTurbine.default(name)
    assert wt.n_seg == wt.chord.size
    assert len(wt.airfoil_files()) == len(wt.airfoil_name)


def test_turbine_yaml_roundtrip(tmp_path, nrel):
    wt, _ = nrel
    wt.to_yaml(tmp_path / "wt.yaml")
    wt2 = WindTurbine.from_yaml(tmp_path / "wt.yaml")
    for k in ("seg", "span", "chord", "twist", "airfoil_index"):
        np.testing.assert_array_equal(getattr(wt, k), getattr(wt2, k))
    assert wt2.omega_rpm == wt.omega_rpm and wt2.airfoil_name == wt.airfoil_name


def test_turbine_old_format_keys():
    wt = WindTurbine.from_dict({
        "Nblade": 2, "Nbeta": 6, "Nseg": 2, "seg": [10, 20], "span": [10, 10], "chord": [2, 1],
        "twist": [10, 5], "airfoil_index": [0, 0],
        "airfoil_geometry": ["./airfoil_geometry/NACA0012.dat"],
        "controler": {"cutin": [3, 6, 0], "rated": [11, 12, 0], "cutout": [25, 12, 20]},
    })
    assert (wt.n_blade, wt.n_beta, wt.airfoil_name) == (2, 6, ["NACA0012"])
    assert wt.controller.get_omega(7.0) == pytest.approx(9.0)
    assert wt.controller.get_pitch(18.0) == pytest.approx(10.0)
    assert wt.with_operating_point(30.0).omega_rpm == 0.0


def test_acoustic_strips_respect_aspect_ratio(nrel):
    wt, _ = nrel
    strips = wt.acoustic_strips(target_l_c=3.0)
    ratio = strips.span[:-1] / strips.chord[:-1]
    assert np.all(ratio >= 3.0 - 0.1)
    assert strips.rotor_diameter == pytest.approx(wt.rotor_diameter, rel=1e-2)


# ------------------------------------------------------------ config / flow
INPUT = """
receivers:
  type: polar
  r: {min: 100, max: 500, n: 3}
  tau: [0, 90, 180, 270]
  z: [2]
flow:
  type: uniform
  u: 10
  positions: [[0, 0, 90], [500, 0, 90]]
source:
  turbine: NREL_5MW
  airfoil_database: DB
  frequencies: [100, 1000]
  n_beta: 2
propagation:
  anything: ignored by the source
"""


def test_single_input_file(tmp_path):
    f = tmp_path / "input.yaml"
    f.write_text(INPUT.replace("DB", str(NREL_DB)))
    src = WindTurbineNoise.from_yaml(f)
    assert isinstance(src.flow, UniformFlow) and src.turbine.n_beta == 2
    np.testing.assert_array_equal(src.hub_positions, [[0, 0, 90], [500, 0, 90]])
    out = src.compute_all(receivers.ReceiversConfig.from_yaml(f))
    assert len(out) == 2
    assert dict(out[1].Spp_TIN.sizes) == {"r": 3, "tau": 4, "z": 1, "beta": 2, "blade": 3, "seg": 14, "freq": 2}
    # polar grid centred on each hub
    np.testing.assert_allclose(out[1].x.sel(r=100, tau=0), 600.0)
    np.testing.assert_allclose(out[1].distance.sel(r=100, tau=0, z=2), np.hypot(100, 88))


def test_unknown_keys_rejected():
    with pytest.raises(Exception):
        SourceConfig.from_dict({"turbine": "NREL_5MW", "frequencies": [100], "typo": 1})
    with pytest.raises(ValueError):
        Flow.from_dict({"type": "does_not_exist"})


def test_parallel_matches_serial(nrel):
    wt, db = nrel
    grid = receivers.cartesian([-100, 100], [0, 200], [2])
    a = WindTurbineNoise(wt, db, FREQ, n_jobs=1, positions=[(0, 0, 90)]).compute(grid)
    b = WindTurbineNoise(wt, db, FREQ, n_jobs=2, positions=[(0, 0, 90)]).compute(grid)
    xr.testing.assert_identical(a, b)


def test_netcdf_output_and_levels(tmp_path, nrel):
    wt, db = nrel
    ds = WindTurbineNoise(wt, db, FREQ, positions=[(0, 0, 90)]).compute(
        receivers.cartesian([200.0], [0.0], [2.0]))
    ds.to_netcdf(tmp_path / "out.nc")
    back = xr.load_dataset(tmp_path / "out.nc")
    xr.testing.assert_allclose(units.spl(back), units.spl(ds))
    rotor_total = units.db_sum(units.spl(ds), ["blade", "seg"])
    assert rotor_total.dims == ("x", "y", "z", "beta", "freq")


# -------------------------------------------------------------------- xfoil
@pytest.mark.xfoil
@pytest.mark.skipif(shutil.which("xfoil") is None, reason="xfoil not installed")
def test_xfoil_database_smoke():
    wt = WindTurbine.default("NREL_5MW")
    db = AirfoilDatabase.build_with_xfoil(
        chord=wt.chord[-1:], airfoil_index=[0], airfoil_files=[wt.airfoil_files()[-1]],
        uinf=[40.0], aoa=[0, 2, 4, 2], virtual_buffer=shutil.which("Xvfb") is not None)
    assert db.ds.aoa.values.tolist() == [0, 2, 4]
    assert np.all(db.ds.delta_star_top.values > 0)
    assert np.all(np.isfinite(db.ds.cl.values))
