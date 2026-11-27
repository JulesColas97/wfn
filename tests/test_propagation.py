"""Unit tests of ``wfn.propagation`` (no WAPE executable needed)."""
import h5py
import numpy as np
import pytest

from wfn import receivers
from wfn.propagation.cases import case_inputs
from wfn.propagation import (
    PropagationConfig,
    EwapeInputs,
    PiecewapeInputs,
    WapeInputs,
    assemble_polar,
    case_dir,
    discretize_rotor,
    loop_angle,
    namelist_text,
    prepare_cases,
    read_namelist,
    read_wape,
    refine_angle,
    third_octave,
    tl_to_deltaL,
    to_grid,
    wape_filename,
    write_namelist,
)
from wfn.propagation import DeltaL
from wfn.units import third_octave_frequencies


def write_fake_wape(fname, r, z, freqs, func, h=None, rec_heights=(2.0,), legacy=False):
    """Write a file with the layout of the WAPE outputs. ``func(r, z, f)`` gives deltaL."""
    h = np.zeros_like(r) if h is None else h
    side, name = ("solution", "deltaL") if legacy else ("side", "TL")
    R, Z = np.meshgrid(r, z, indexing="ij")
    with h5py.File(fname, "w") as f:
        f["mesh/x"] = r
        f["mesh/z"] = z
        f["mesh/h"] = h
        f["mesh/x_grid"] = R.T
        f["mesh/z_grid"] = (Z + h[:, None]).T
        f["receiver/heights"] = np.asarray(rec_heights)
        if not legacy:
            f["flow/u"] = (10 + 0 * R).T
            f["flow/c"] = (340 + 0 * R).T
        for fr in freqs:
            key = format(int(fr), "04d")
            f[f"{side}/{key}/{name}"] = func(R, Z, fr).T
            Rr, Zr = np.meshgrid(r, np.asarray(rec_heights), indexing="ij")
            f[f"receiver/{key}/{name}"] = func(Rr, Zr, fr).T


def base_cfg(tmp_path, **kw):
    data = dict(
        work_dir=str(tmp_path),
        positions=[[100.0, 50.0]],
        source_heights=[40, 80],
        tau=[0, 90, 180, 270],
        frequencies=[50, 100],
        domain=dict(xmin=0, xmax=300, ymin=0, ymax=200, zmax=30),
        wape=dict(case_name="c0", dx_out=10, receiver_heights=[2.0], side=True),
    )
    data.update(kw)
    return PropagationConfig.from_dict(data)


# ----------------------------------------------------------------- config
def test_config_from_yaml(tmp_path):
    (tmp_path / "in.yaml").write_text(
        "propagation:\n  work_dir: runs\n  rotor: {hub_height: 90, diameter: 120, n: 3}\n"
        "  wape: {frequencies: [100], tau: [0, 180], pml_param: 5e4}\n"
    )
    cfg = PropagationConfig.from_yaml(tmp_path / "in.yaml")
    assert cfg.work_path == tmp_path / "runs"
    assert cfg.heights() == [50, 90, 130]
    assert cfg.angles() == [0, 180]
    assert cfg.freqs() == [100]
    assert cfg.wape.pml_param == 5e4


def test_config_unknown_key():
    with pytest.raises(Exception):
        PropagationConfig.from_dict({"wape": {"not_a_key": 1}})


def test_discretize_rotor():
    np.testing.assert_array_equal(discretize_rotor(90, 126, 7), [36, 54, 72, 90, 108, 126, 144])


# --------------------------------------------------------------- namelist
def test_namelist_format():
    txt = namelist_text(PiecewapeInputs(case_name="a", frequencies=[50.0, 100.0], impedance_mode=0))
    lines = txt.splitlines()
    assert lines[0] == "$inputs" and lines[-1] == "$end inputs"
    assert "case_name='a'" in lines
    assert "frequencies(1)=50.0, 100.0" in lines
    assert not any(l.startswith("nb_") for l in lines)
    assert "impedance_mode=0" in lines
    assert "Lx1=300.0" in lines  # Fortran default, written explicitly
    assert not any(l.startswith(("H=", "dz_out")) for l in lines)  # EWAPE-only keys


def test_namelist_roundtrip(tmp_path):
    w = PiecewapeInputs(case_name="c0", frequencies=[100.0], tau=[0.0, 45.0], Lx1=1000.0, Lz=200.0,
                   impedance_mode=0, receiver_heights=[2.0, 10.0], fname_flow="/a/flow.h5")
    path = write_namelist(w, tmp_path)
    assert read_namelist(path) == w


def test_angles_from_taumin_taumax(tmp_path):
    cfg = base_cfg(tmp_path, tau=None, wape=dict(taumin=0, taumax=360, nb_tau=4))
    assert cfg.angles() == [0, 90, 180, 270]  # taumax excluded, as in the Fortran codes
    cfg = base_cfg(tmp_path, tau=None, wape=dict(taumin=30, taumax=360, nb_tau=1))
    assert cfg.angles() == [30]


def test_case_inputs_explicit_tau(tmp_path):
    cfg = base_cfg(tmp_path, tau=[0, 45], wape=dict(taumin=0, taumax=360, nb_tau=4))
    w = case_inputs(cfg, 0, 40)
    assert w.tau == [0, 45] and (w.taumin, w.taumax, w.nb_tau) == (None, None, None)


def test_namelist_backend(tmp_path):
    w = EwapeInputs(case_name="a", source_height=40.0, H=50.0, pml_size=30)
    txt = namelist_text(w)
    assert "source_height=40.0" in txt and "H=50.0" in txt and "pml_size=30" in txt
    assert "H=" not in namelist_text(PiecewapeInputs())
    assert "pe_model" not in txt  # piecewape only
    assert read_namelist(write_namelist(w, tmp_path), "ewape") == w
    assert WapeInputs.for_backend("ewape") is EwapeInputs
    with pytest.raises(ValueError):
        WapeInputs.for_backend("unknown")


def test_wrong_backend_key():
    with pytest.raises(Exception):
        PropagationConfig.from_dict({"backend": "ewape", "wape": {"pe_model": 2}})
    cfg = PropagationConfig.from_dict({"backend": "ewape", "wape": {"H": 50}})
    assert isinstance(cfg.wape, EwapeInputs) and cfg.wape.H == 50


# ------------------------------------------------------------------ cases
def test_prepare_cases(tmp_path):
    exe = tmp_path / "2D_WAPE_src"
    exe.write_text("#!/bin/bash\n")
    cfg = base_cfg(tmp_path / "runs", executable=str(exe))
    dirs = prepare_cases(cfg)
    assert dirs == [tmp_path / "runs/t0/40", tmp_path / "runs/t0/80"]
    w = read_namelist(dirs[1] / "input.nml")
    assert (w.Lx1, w.Lx2, w.Ly1, w.Ly2, w.Lz) == (200, 100, 150, 50, 30)
    assert (w.Xs_flow, w.Ys_flow) == (100, 50)
    assert w.source_height == 80 and w.tau == [0, 90, 180, 270] and w.frequencies == [50, 100]
    assert (dirs[0] / "2D_WAPE").is_file()
    run = (tmp_path / "runs/run.sh").read_text().splitlines()
    assert run[1:4] == ["cd t0/40/", "time ./2D_WAPE >out.out", "cd - >/dev/null"]


def test_missing_executable(tmp_path, monkeypatch):
    monkeypatch.delenv("WFN_WAPE", raising=False)
    with pytest.raises(FileNotFoundError):
        prepare_cases(base_cfg(tmp_path))


# ---------------------------------------------------------------- results
def test_wape_filename():
    assert wape_filename("d", "c0", 0).name == "c0_0000.h5"
    assert wape_filename("d", "c0", 22.5).name == "c0_22.5.h5"


def test_read_wape(tmp_path):
    r = np.arange(0, 100, 10.0)
    z = np.arange(0, 30, 10.0)
    func = lambda R, Z, f: R + Z + f
    write_fake_wape(tmp_path / "a.h5", r, z, [100, 50], func, rec_heights=(2.0, 5.0))
    ds = read_wape(tmp_path / "a.h5")
    assert ds.deltaL.dims == ("r", "z", "freq")
    np.testing.assert_array_equal(ds.freq, [50, 100])
    assert float(ds.deltaL.sel(r=20, z=10, freq=50)) == 80
    assert np.isnan(ds.deltaL.sel(r=90, z=20, freq=100))  # > 120 dB
    rec = read_wape(tmp_path / "a.h5", side=False)
    np.testing.assert_array_equal(rec.z, [2, 5])
    assert float(rec.deltaL.sel(r=30, z=5, freq=50)) == 85
    assert float(ds.u.sel(r=20, z=10)) == 10 and float(ds.c.sel(r=20, z=10)) == 340
    assert "u" not in rec


def test_read_wape_legacy_layout(tmp_path):
    """Files of the old code: /solution/<freq>/deltaL, no flow."""
    r = np.arange(0, 100, 10.0)
    z = np.arange(0, 30, 10.0)
    write_fake_wape(tmp_path / "a.h5", r, z, [50], lambda R, Z, f: R + Z + f, legacy=True)
    ds = read_wape(tmp_path / "a.h5")
    assert float(ds.deltaL.sel(r=20, z=10, freq=50)) == 80 and "u" not in ds


def make_runs(cfg, func, z=None):
    """Fake WAPE outputs for every (height, tau); the run along tau stops at r = 300."""
    r = np.arange(0, 310, 10.0)
    z = np.arange(0, 30, 10.0) if z is None else z
    for h in cfg.heights():
        d = case_dir(cfg, 0, h)
        d.mkdir(parents=True, exist_ok=True)
        for t in cfg.angles():
            write_fake_wape(wape_filename(d, "c0", t), r, z, cfg.freqs(),
                            lambda R, Z, f: func(R, Z, f, h, t))


def test_assemble_polar(tmp_path):
    cfg = base_cfg(tmp_path)
    make_runs(cfg, lambda R, Z, f, h, t: 0.01 * R + Z / 10 + f / 100 + h / 40 + t / 90)
    ds = assemble_polar(cfg, 0)
    assert ds.deltaL.dims == ("r", "tau", "z", "freq", "source_height")
    # rmax = hypot(200, 150) + 50 = 300
    np.testing.assert_allclose(ds.r, np.arange(0, 310, 10))
    np.testing.assert_allclose(ds.z, [0, 10, 20])
    assert ds.attrs["center_x"] == 100 and ds.attrs["topo"] == 0
    v = ds.deltaL.sel(r=100, tau=90, z=20, freq=100, source_height=80)
    assert float(v) == pytest.approx(1 + 2 + 1 + 2 + 1)
    # last value is overwritten by the previous one (as in the old code)
    np.testing.assert_array_equal(ds.deltaL.isel(r=-1), ds.deltaL.isel(r=-2))


def test_assemble_polar_crops_pml_rows(tmp_path):
    """piecewape also writes the PML rows above Lz: they must be dropped."""
    cfg = base_cfg(tmp_path)
    make_runs(cfg, lambda R, Z, f, h, t: 0.01 * R + Z / 10, z=np.arange(0, 50, 10.0))
    ds = assemble_polar(cfg, 0)
    np.testing.assert_allclose(ds.z, [0, 10, 20])


def test_assemble_polar_ewape(tmp_path):
    """EWAPE side plane starts at z = -H (water)."""
    cfg = base_cfg(tmp_path, backend="ewape", wape=dict(H=20.0, dz_out=10.0, dx_out=10.0))
    z = -20 + 10.0 * np.arange(5)  # floor((20 + 30) / 10) points
    make_runs(cfg, lambda R, Z, f, h, t: 0.01 * R + Z / 10, z=z)
    ds = assemble_polar(cfg, 0)
    np.testing.assert_allclose(ds.z, z)


def test_loop_refine_angle(tmp_path):
    cfg = base_cfg(tmp_path)
    make_runs(cfg, lambda R, Z, f, h, t: t / 90 + 0 * R)
    ds = loop_angle(assemble_polar(cfg, 0))
    np.testing.assert_array_equal(ds.tau, [0, 90, 180, 270, 360])
    assert float(ds.deltaL.sel(tau=360).max()) == 0
    ref = refine_angle(ds)
    assert ref.sizes["tau"] == 9
    assert float(ref.deltaL.sel(tau=135, r=50, z=0, freq=50, source_height=40)) == pytest.approx(1.5)
    np.testing.assert_allclose(ref.x.sel(tau=90), 100, atol=1e-9)


def test_to_grid(tmp_path):
    cfg = base_cfg(tmp_path)
    make_runs(cfg, lambda R, Z, f, h, t: t / 90 + 0 * R)
    ds = loop_angle(assemble_polar(cfg, 0))
    grid = receivers.cartesian([100, 150, 200], [50, 100], [0, 10])  # centre at (100, 50)
    out = to_grid(ds, grid)
    assert out.deltaL.dims == ("x", "y", "z", "freq", "source_height")
    assert float(out.deltaL.sel(x=200, y=50, z=10, freq=50, source_height=40)) == 0  # tau = 0
    assert float(out.deltaL.sel(x=100, y=100, z=0, freq=50, source_height=40)) == 1  # tau = 90
    pts = to_grid(ds, receivers.points([100], [0], [10]))  # tau = 270
    assert float(pts.deltaL.isel(point=0, freq=0, source_height=0)) == 3


def test_third_octave(tmp_path):
    fc, nfc = [100, 200], [1, 2]
    freqs = third_octave_frequencies(fc, nfc)
    cfg = base_cfg(tmp_path, frequencies=freqs.tolist(), source_heights=[40], tau=[0])
    vals = {freqs[0]: 3.0, freqs[1]: 0.0, freqs[2]: 10 * np.log10(3)}
    make_runs(cfg, lambda R, Z, f, h, t: vals[f] + 0 * R)
    out = third_octave(assemble_polar(cfg, 0), fc, nfc)
    np.testing.assert_array_equal(out.freq, fc)
    np.testing.assert_allclose(out.deltaL.sel(freq=100), 3.0)
    np.testing.assert_allclose(out.deltaL.sel(freq=200), 10 * np.log10(2))
    with pytest.raises(ValueError):
        third_octave(out, fc, nfc)


def test_field_attribute_is_inherited(tmp_path):
    r = np.arange(0, 100, 10.0)
    z = np.arange(0, 30, 10.0)
    write_fake_wape(tmp_path / "a.h5", r, z, [50], lambda R, Z, f: R + Z + f)
    write_fake_wape(tmp_path / "b.h5", r, z, [50], lambda R, Z, f: R + Z + f, legacy=True)
    assert read_wape(tmp_path / "a.h5").attrs["field"] == "TL"
    assert read_wape(tmp_path / "b.h5").attrs["field"] == "deltaL"

    cfg = base_cfg(tmp_path)
    make_runs(cfg, lambda R, Z, f, h, t: 0.01 * R)
    dl = DeltaL.from_cases(cfg, 0)
    assert dl.is_tl
    assert dl.loop_angle().is_tl and dl.refine_angle().is_tl
    grid = receivers.cartesian([100.0, 150.0], [50.0, 80.0], [10.0])
    assert dl.loop_angle().to_grid(grid).is_tl
    dl.save(tmp_path / "dl.nc")
    assert DeltaL.load(tmp_path / "dl.nc").is_tl


def test_tl_to_deltaL():
    r = np.array([0.0, 30.0, 60.0])
    z = np.array([-5.0, 5.0])
    ds = receivers.polar(r, [0.0, 90.0], z, center=(10.0, 20.0))
    ds = ds.assign_coords(freq=[100.0], source_height=[40.0])
    ds["deltaL"] = (("r", "tau", "z", "freq", "source_height"), np.full((3, 2, 2, 1, 1), -50.0))
    ds["surface"] = (("r", "tau"), np.zeros((3, 2)))
    ds.attrs.update(field="TL")
    out = tl_to_deltaL(ds)
    assert out.attrs["field"] == "deltaL" and out.deltaL.dims == ds.deltaL.dims
    # above the surface: 50 + 20 log10(Rff); below: + 26 dB
    rff = np.sqrt(30.0**2 + (5.0 - 40.0) ** 2)
    assert float(out.deltaL.sel(r=30, tau=0, z=5)) == pytest.approx(50 + 20 * np.log10(rff))
    rff = np.sqrt(30.0**2 + (-5.0 - 40.0) ** 2)
    assert float(out.deltaL.sel(r=30, tau=90, z=-5)) == pytest.approx(50 + 20 * np.log10(rff) + 26)
    # r = 0 is clipped to Rff >= 1 m (no -inf) and a second call does nothing
    assert np.isfinite(out.deltaL).all()
    assert tl_to_deltaL(out) is out
    assert DeltaL(ds).to_deltaL().ds.attrs["field"] == "deltaL"
