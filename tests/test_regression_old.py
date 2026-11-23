"""The new source model must reproduce the old ``prepost`` WindTurbineNoise."""
import numpy as np
import pytest

from wfn import receivers
from wfn.flow import UniformFlow
from wfn.source import AirfoilDatabase, AtmosphereConfig, WindTurbine, WindTurbineNoise

from conftest import NREL_DB, OLD_NREL

pytestmark = pytest.mark.legacy

FREQ = [50, 80, 100, 200, 400, 500, 600, 800, 1000, 1200, 1500, 1800, 2000]
HUB = (0.0, 0.0, 90.0)


def fix_old_yaw_bug(module):
    """In the old ``WindTurbine``, ``yaw_rad`` is declared with ``@omega.setter``, so
    it returns omega (rad/s) instead of the yaw: the old code applied a yaw of
    ``omega`` radians in the blade frame transforms. Patch it to compare physics."""
    module.WindTurbine.yaw_rad = property(lambda self: np.deg2rad(self.yaw_deg))


def run_old(module, in_dir, x, y, z, u_ref, alpha):
    fix_old_yaw_bug(module)
    in_dir(OLD_NREL)
    src = module.WindTurbineNoise(yaml_fname="./input_NREL.yaml")
    src.inputs.frequencies = np.array(FREQ, dtype=float)
    src.wt.airfoil_data.createAllBlDataInterpoland()
    src.mesh.x_array = np.asarray(x, dtype=float)
    src.mesh.y_array = np.asarray(y, dtype=float)
    src.mesh.z_array = np.asarray(z, dtype=float)
    src.mesh.createCartesianGrid()
    src._compute_blade_segment_position()
    src.inputs.u_inf = u_ref if alpha is None else u_ref * (src.z_rotor / HUB[2]) ** alpha
    src.inputs.epsilon = 0.01
    src.setFlowAtRotorPosition()
    src.computeNoiseFromRotor(n_jobs=1)
    return src


def run_new(x, y, z, u_ref, alpha):
    flow = UniformFlow(u=u_ref, epsilon=0.01, shear_alpha=alpha, z_ref=None if alpha is None else HUB[2])
    model = WindTurbineNoise(WindTurbine.default("NREL_5MW"), AirfoilDatabase.open(NREL_DB), FREQ,
                             flow=flow, atmosphere=AtmosphereConfig(), positions=[HUB])
    return model.compute(receivers.cartesian(x, y, z))


def as_old_layout(da):
    """(x, y, z, beta, blade, seg, freq) -> (x, y, z, K, 1, seg, freq), K = iblade*n_beta + ibeta."""
    arr = da.transpose("x", "y", "z", "blade", "beta", "seg", "freq").values
    s = arr.shape
    return arr.reshape(s[0], s[1], s[2], s[3] * s[4], 1, s[5], s[6])


@pytest.mark.parametrize("u_ref, alpha", [(8.0, None), (15.0, 0.2)])
def test_single_receiver(old_prepost, in_dir, u_ref, alpha):
    pos = ([1000.0], [600.0], [90.0])
    old = run_old(old_prepost, in_dir, *pos, u_ref, alpha)
    new = run_new(*pos, u_ref, alpha)
    np.testing.assert_allclose(as_old_layout(new.Spp_TEN), old.Spp_TEN, rtol=1e-10)
    np.testing.assert_allclose(as_old_layout(new.Spp_TIN), old.Spp_TIN, rtol=1e-10)
    np.testing.assert_allclose(new.AoA.transpose("blade", "beta", "seg").values.reshape(-1, 14),
                               old.AoA_rotor, rtol=1e-12)


def test_grid(old_prepost, in_dir):
    pos = (np.linspace(-300, 300, 5), np.linspace(-200, 250, 4), [2.0, 50.0])
    old = run_old(old_prepost, in_dir, *pos, 15.0, 0.2)
    new = run_new(*pos, 15.0, 0.2)
    np.testing.assert_allclose(as_old_layout(new.Spp_TEN), old.Spp_TEN, rtol=1e-10)
    np.testing.assert_allclose(as_old_layout(new.Spp_TIN), old.Spp_TIN, rtol=1e-10)
    # SWL of the old code with the new helpers
    from wfn import units
    np.testing.assert_allclose(as_old_layout(units.swl(new)), old.SWL_tot, rtol=0, atol=1e-9)
