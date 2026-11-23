"""The new propagation module must reproduce the old ``prepost`` WAPE pre/post-processing."""
import numpy as np
import pytest

from wfn import receivers
from wfn.config import load_yaml
from wfn.propagation import PropagationConfig, assemble_polar, namelist_text, read_namelist, read_wape, to_grid

from conftest import OLD_CODE, ROOT

pytestmark = pytest.mark.legacy

OLD_EXAMPLES = ROOT / "old_project" / "examples"
NO_FLOW_H5 = OLD_EXAMPLES / "1_simple_propa" / "no_flow" / "noFlow_0000.h5"
WALBERLA = OLD_EXAMPLES / "2_multi_wape" / "walberla"


@pytest.fixture(scope="module")
def old(old_prepost):
    import inspect
    import re

    from prepost import deltaLfield, simu, wape
    from prepost.source import mesh

    # the old ``simu.read_wape_inputs`` uses ``re`` and ``inspect`` without importing them
    simu.re, simu.inspect = re, inspect
    return dict(simu=simu, wape=wape, deltaLfield=deltaLfield, mesh=mesh)


@pytest.mark.skip(reason="intentional deviation: keys renamed (c0->c_air, u0->u_ref, cfl->dx_frac), group `inputs`, no nb_* counters")
def test_namelist_same_as_old(old, tmp_path):
    yaml_file = OLD_CODE / "input_default.yaml"
    s = old["simu"].Simu(yaml_fname=str(yaml_file))
    s.writeWapeInputs(directory=str(tmp_path))
    cfg = PropagationConfig.from_dict({"wape": load_yaml(yaml_file)["wape_inputs"]})
    assert namelist_text(cfg.wape) == (tmp_path / "input.nml").read_text()


@pytest.mark.parametrize("side", [True, False])
def test_read_wape_same_as_old(old, side):
    if not NO_FLOW_H5.is_file():
        pytest.skip("old WAPE output not available")
    res = old["wape"].PeResults()
    res.read_carto(str(NO_FLOW_H5)) if side else res.read_receiver(str(NO_FLOW_H5))
    ds = read_wape(NO_FLOW_H5, side=side)
    np.testing.assert_array_equal(ds.freq, res.frequencies)
    np.testing.assert_array_equal(ds.r, res.x)
    np.testing.assert_array_equal(ds.z, res.z)
    np.testing.assert_array_equal(ds.h, res.h)
    np.testing.assert_array_equal(ds.z_abs, res.z_grid)
    if side:
        np.testing.assert_array_equal(ds.deltaL, res.deltaL)
    else:  # the old receiver reader does not clip |deltaL| > 120 dB
        valid = np.abs(res.deltaL) <= 120
        np.testing.assert_array_equal(ds.deltaL.values[valid], res.deltaL[valid])


@pytest.fixture(scope="module")
def walberla(old):
    """Old and new polar deltaL for the source 3, height 40 m of the Walberla example."""
    nml = WALBERLA / "t3" / "40" / "input.nml"
    if not (nml.parent / "c0_0000.h5").is_file():
        pytest.skip("Walberla WAPE outputs not available")
    dl = old["deltaLfield"].DeltaLField(nml_path=str(nml), side=True)
    dl.concatenateAllAngles()

    w = read_namelist(nml)
    cfg = PropagationConfig(wape=w, work_dir=str(WALBERLA), positions=[[0, 0]] * 3 + [[w.Xs_flow, w.Ys_flow]])
    ds = assemble_polar(cfg, isource=3)
    return dl, ds


def test_assemble_polar_same_as_old(walberla):
    dl, ds = walberla
    np.testing.assert_array_equal(ds.r, dl.mesh.r_array)
    np.testing.assert_array_equal(ds.tau, dl.mesh.tau_array)
    np.testing.assert_array_equal(ds.z, dl.mesh.z_array)
    np.testing.assert_array_equal(ds.deltaL, dl.deltaL)
    np.testing.assert_array_equal(ds.terrain, dl.terrain)
    np.testing.assert_allclose(ds.x.transpose("r", "tau"), dl.mesh.x_grid[:, :, 0])
    np.testing.assert_allclose(ds.y.transpose("r", "tau"), dl.mesh.y_grid[:, :, 0])


def test_to_grid_same_as_old(old, walberla):
    dl, ds = walberla
    x, y, z = np.arange(0, 9000, 50.0), np.array([600.0]), np.arange(0, 400, 25.0)
    new_mesh = old["mesh"].Mesh()
    new_mesh.x_array, new_mesh.y_array, new_mesh.z_array = x, y, z
    new_mesh.createCartesianGrid()
    dl.interpolate(new_mesh)

    out = to_grid(ds, receivers.cartesian(x, y, z))
    np.testing.assert_array_equal(out.deltaL, dl.deltaL)
    np.testing.assert_array_equal(out.terrain.isel(z=0), dl.terrain)
