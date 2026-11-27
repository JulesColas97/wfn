"""Tests of ``wfn.spl.combine`` (synthetic datasets, no WAPE / XFoil needed)."""
import numpy as np
import pytest
import xarray as xr

from wfn import receivers, spl, units
from wfn.propagation import DeltaL

FREQ = np.array([100.0, 200.0])
HEIGHTS = np.array([50.0, 100.0])


def fake_source(grid, z_seg, freq=FREQ, seed=0):
    """Source result with the layout of ``WindTurbineNoise.compute``; ``z_seg`` is ``(beta, blade, seg)``."""
    z_seg = np.asarray(z_seg, dtype=float)
    X, Y, Z, dims, coords = receivers.receiver_positions(grid)
    shape = (*X.shape, *z_seg.shape, freq.size)
    rng = np.random.default_rng(seed)
    out_dims = (*dims, "beta", "blade", "seg", "freq")
    rdims = ("beta", "blade", "seg")
    ds = xr.Dataset(
        {
            "Spp_TEN": (out_dims, rng.uniform(1, 2, shape)),
            "Spp_TIN": (out_dims, rng.uniform(1, 2, shape)),
            "distance": (dims, np.sqrt(X**2 + Y**2 + (Z - 75) ** 2)),
            "z_seg": (rdims, z_seg),
        },
        coords={**coords, "freq": freq, "beta": np.arange(z_seg.shape[0]),
                "blade": np.arange(z_seg.shape[1]), "seg": np.arange(z_seg.shape[2])},
        attrs={**grid.attrs, "hub_x": 0.0, "hub_y": 0.0, "hub_z": 75.0},
    )
    return ds


def fake_deltaL_cart(grid, values, freq=FREQ, heights=HEIGHTS):
    """deltaL on a receiver grid, ``values[ih]`` is the (constant) deltaL of height ``ih``."""
    X, _, _, dims, coords = receivers.receiver_positions(grid)
    dl = np.broadcast_to(np.asarray(values, float), (*X.shape, freq.size, heights.size)).copy()
    ds = xr.Dataset(coords=coords).assign_coords(freq=freq, source_height=heights)
    ds["deltaL"] = ((*dims, "freq", "source_height"), dl)
    ds["terrain"] = (dims, np.zeros(X.shape))
    return ds.assign_attrs(grid.attrs)


def fake_deltaL_polar(func, tau=(0, 90, 180, 270), freq=FREQ, heights=HEIGHTS):
    """Polar deltaL like ``assemble_polar``; ``func(r, tau)`` gives the value (same for all z, f, h)."""
    r = np.arange(0, 400, 5.0)
    ds = receivers.polar(r, tau, [2.0])
    ds = ds.assign_coords(freq=freq, source_height=heights)
    R, T = np.meshgrid(r, np.asarray(tau, float), indexing="ij")
    val = func(R, T)[:, :, None, None, None] * np.ones((1, 1, 1, freq.size, heights.size))
    ds["deltaL"] = (("r", "tau", "z", "freq", "source_height"), val)
    ds["terrain"] = (("r", "tau"), np.zeros(R.shape))
    ds["surface"] = (("r", "tau"), np.zeros(R.shape))
    ds.attrs.update(side=1, topo=0, topo_following=0, isource=0)
    return ds


@pytest.fixture
def cart():
    return receivers.cartesian(np.linspace(-100, 100, 5), np.linspace(-100, 100, 4), [2.0])


def test_same_grid(cart):
    assert spl.same_grid(cart, cart)
    assert not spl.same_grid(cart, receivers.cartesian([0, 1], [0, 1], [2.0]))
    assert not spl.same_grid(cart, receivers.polar([0, 10], [0, 90], [2.0]))
    pts = receivers.points([0, 1], [0, 1], [2, 2])
    assert spl.same_grid(pts, receivers.points([0, 1], [0, 1], [2, 2]))
    assert not spl.same_grid(pts, receivers.points([0, 1], [0, 2], [2, 2]))


def test_constant_deltaL(cart):
    src = fake_source(cart, np.full((2, 3, 4), 70.0))
    out = spl.combine(src, fake_deltaL_cart(cart, [6.0, 6.0]))
    xr.testing.assert_allclose(out.Spp_TEN, src.Spp_TEN * 10**0.6)
    xr.testing.assert_allclose(out.Spp_TIN, src.Spp_TIN * 10**0.6)
    assert out.Spp_TEN.dims == src.Spp_TEN.dims
    np.testing.assert_allclose(out.deltaL, 6.0)
    # the units helpers work on the result
    lvl = units.db_sum(units.spl(out), ("blade", "seg"))
    ref = units.db_sum(units.spl(src), ("blade", "seg")) + 6.0
    xr.testing.assert_allclose(lvl, ref)


def test_height_interpolation(cart):
    # segment heights: exactly on the first height, in the middle, on the second, outside both ends
    z_seg = np.array([50.0, 75.0, 100.0, 20.0, 150.0]).reshape(1, 1, 5)
    src = fake_source(cart, z_seg)
    out = spl.combine(src, fake_deltaL_cart(cart, [0.0, 10.0]))
    factor = (out.Spp_TEN / src.Spp_TEN).isel(x=0, y=0, z=0, freq=0, beta=0, blade=0).values
    np.testing.assert_allclose(factor, [1.0, 5.5, 10.0, 1.0, 10.0])


@pytest.mark.parametrize("scalar", [False, True])
def test_single_source_height(cart, scalar):
    z_seg = np.array([50.0, 75.0, 100.0]).reshape(1, 1, 3)
    src = fake_source(cart, z_seg)
    dl = fake_deltaL_cart(cart, [4.0], heights=np.array([80.0]))
    if scalar:
        dl = dl.isel(source_height=0)  # source_height becomes a scalar coordinate
    out = spl.combine(src, dl)
    factor = (out.Spp_TEN / src.Spp_TEN).isel(x=0, y=0, z=0, freq=0, beta=0, blade=0).values
    np.testing.assert_allclose(factor, 10 ** 0.4)


def test_polar_deltaL_to_source_grid(cart):
    src = fake_source(cart, np.full((1, 1, 2), 75.0)).assign_attrs(isource=0)
    dl = fake_deltaL_polar(lambda r, t: np.full_like(r, 3.0))
    for d in (dl, DeltaL(dl)):
        out = spl.combine(src, d)
        assert out.Spp_TEN.dims == src.Spp_TEN.dims
        np.testing.assert_allclose(out.deltaL, 3.0)
        assert "terrain" in out and out.terrain.dims == ("x", "y")


def test_polar_full_circle_is_looped():
    # receivers between tau = 270 and 360: need the tau = 360 copy to be inside the polar grid
    grid = receivers.points([100.0], [-10.0], [2.0])
    src = fake_source(grid, np.full((1, 1, 1), 75.0))
    out = spl.combine(src, fake_deltaL_polar(lambda r, t: np.full_like(r, 4.0)))
    np.testing.assert_allclose(out.deltaL, 4.0)


def test_isource_mismatch(cart):
    src = fake_source(cart, np.full((1, 1, 1), 75.0)).assign_attrs(isource=1)
    with pytest.raises(ValueError, match="isource"):
        spl.combine(src, fake_deltaL_polar(lambda r, t: 0 * r))


def test_missing_frequency(cart):
    src = fake_source(cart, np.full((1, 1, 1), 75.0), freq=np.array([100.0, 300.0]))
    with pytest.raises(ValueError, match="300"):
        spl.combine(src, fake_deltaL_cart(cart, [0.0, 0.0]))


def test_extra_deltaL_frequencies(cart):
    src = fake_source(cart, np.full((1, 1, 1), 75.0), freq=np.array([200.0]))
    dl = fake_deltaL_cart(cart, [0.0, 0.0], freq=np.array([100.0, 200.0, 400.0]))
    dl["deltaL"] = dl.deltaL + xr.DataArray([1.0, 2.0, 4.0], dims="freq", coords={"freq": dl.freq})
    out = spl.combine(src, dl)
    np.testing.assert_allclose(out.deltaL, 2.0)


def test_common_grid(cart):
    src = fake_source(cart, np.full((1, 1, 1), 75.0))
    new = receivers.cartesian(np.linspace(-50, 50, 3), [0.0], [2.0])
    out = spl.combine(src, fake_deltaL_polar(lambda r, t: np.full_like(r, 2.0)), grid=new)
    assert out.sizes["x"] == 3 and out.sizes["y"] == 1
    assert spl.same_grid(out, new)
    np.testing.assert_allclose(out.distance.values[:, 0, 0], np.sqrt(np.array([50.0, 0, 50]) ** 2 + 73**2))
    ref = src.Spp_TEN.interp(x=new.x.values, y=new.y.values) * 10**0.2
    xr.testing.assert_allclose(out.Spp_TEN, ref)


def test_non_cartesian_source_to_other_grid():
    src = fake_source(receivers.points([0, 1], [0, 1], [2, 2]), np.full((1, 1, 1), 75.0))
    with pytest.raises(ValueError, match="cartesian"):
        spl.combine(src, fake_deltaL_polar(lambda r, t: 0 * r), grid=receivers.cartesian([0], [0], [2]))


def test_combine_refuses_tl(cart):
    src = fake_source(cart, np.full((1, 1, 1), 75.0))
    dl = fake_deltaL_cart(cart, [0.0, 0.0]).assign_attrs(field="TL")
    with pytest.raises(ValueError, match="to_deltaL"):
        spl.combine(src, dl)
