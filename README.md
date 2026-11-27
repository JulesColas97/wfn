# wfn — wind farm noise

Rewrite of `old_project/` (package `prepost`) as independent modules that share
one input file format and `xarray` data structures.

```
src/wfn/
  config.py        load the YAML file, dict -> dataclass (cattrs)
  receivers.py     receiver grids as xr.Dataset (cartesian / polar / points), replaces Mesh
  units.py         SPL/SWL from PSD, dB sums, A-weighting, third-octave frequencies
  flow/            Flow (abstract interface) + UniformFlow (uniform / power-law wind)
  source/          wind turbine noise source model (TEN + TIN)
    turbine.py     WindTurbine geometry, Controller, acoustic strips
    airfoil.py     AirfoilDatabase (BL data on (uinf, aoa, seg)), XFoil generation, legacy import
    xfoil.py       XFoil runner and post-processing
    kernels.py     numba kernels: Fresnel integrals, Lee and Goody wall pressure spectra
    geometry.py    rotor kinematics and blade segment frames
    amiet.py       TEN / TIN formulas for one segment
    model.py       WindTurbineNoise: orchestrates and returns an xr.Dataset
  propagation/     WAPE (next step)
```

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest                      # legacy comparison tests need old_project/
```

## One input file, one section per module

```yaml
receivers: {type: cartesian, x: {min: -200, max: 200, n: 41}, y: {min: -200, max: 200, n: 41}, z: [2]}
flow:      {type: uniform, u: 15, epsilon: 0.01, shear_alpha: 0.2, z_ref: 90}
source:
  turbine: NREL_5MW                   # package default, a YAML path, or an inline dict
  airfoil_database: NREL5MW_14seg.nc
  frequencies: [50, 100, 200, 500, 1000]
  positions: [[0, 0, 90]]
  output: results/nrel                # -> results/nrel_<i>.nc
propagation: {...}                    # read by wfn.propagation only
```

```python
from wfn.receivers import ReceiversConfig
from wfn.source import WindTurbineNoise
from wfn import units

src = WindTurbineNoise.from_yaml("input.yaml")
results = src.compute_all(ReceiversConfig.from_yaml("input.yaml"))   # one xr.Dataset per turbine
spl = units.spl(results[0])                                            # dims (x, y, z, beta, blade, seg, freq)
```

See `examples/source_nrel/`.

## Source output

`WindTurbineNoise.compute` returns an `xr.Dataset` with:

- `Spp_TEN`, `Spp_TIN`: PSD per segment, dims `(*grid, beta, blade, seg, freq)`
  (linear, so summing over `blade`/`seg` gives the rotor total);
- `distance`: hub-receiver distance, used by `units.swl`;
- `x_seg, y_seg, z_seg, U_inf, U_rel, AoA, epsilon` on `(beta, blade, seg)`.
  The coupling with propagation interpolates ΔL at `z_seg`.

## Old airfoil databases

```bash
python scripts/convert_legacy_airfoil.py old.dat new.nc     # or AirfoilDatabase.open("old.dat")
```

To build a database with XFoil and inspect it, see `examples/airfoil_nrel/`:
`0_create_bl_data.py` builds it, `1_plot_bl.py` makes static plots and the comparison with the
NREL polars, and `2_explore_bl.py` opens the interactive Dash explorer
(`pip install -e ".[explore]"`). The plotting helpers are in `wfn.source.explore`.

## Differences with the old code

- **Yaw bug fixed.** In the old `WindTurbine`, `yaw_rad` was declared with `@omega.setter`, so it
  returned omega (rad/s). The blade frame transforms therefore used a yaw of omega radians
  (~73° at 12.1 rpm). Rotor-averaged levels change by up to ~1.5 dB depending on the receiver.
- **Blade positions for the flow sampling fixed.** The old `_compute_blade_segment_position` put
  the blade at `y = +sin(beta) r`, the mirror image of the rotor used in the acoustic transforms
  (`y = -sin(beta) r`, clockwise seen from upwind). This matters only for flows that are not
  symmetric in y (wakes, veer).
- **Yaw is now handled consistently.** The segment positions, the Doppler angles and the
  transforms all use it.
- **Source/propagation coupling weights fixed.** The old `SplField.combine_linear` interpolated
  ΔL between the two nearest source heights as `ratio*dL_left + (1-ratio)*dL_right`, which swaps
  the weights. Its `searchsorted(side="left")` also picked the wrong neighbour when a segment sat
  exactly on a source height. `wfn.spl.combine` uses plain linear interpolation (in energy) and
  clips `z_seg` to the range of source heights, as before.
- **Third-octave bands from the linear PSD.** `units.spl` is a level per Hz (from the PSD), so
  the OASPL needs third-octave band levels first (`units.spl_third_octave`, band energy =
  bandwidth × mean PSD in the band). The old `SplField.computeThirdOctave` passed the SPL in dB to
  `integrateThirdOctave`, which multiplied dB values by the bandwidth.
- **Smaller changes.** `rho0`/`nu0` are now passed to the wall pressure models instead of being
  hard-coded. `c0` works when it is given. Every Clauser-parameter case is handled. XFoil runs in
  its own temp dir, so the 64-character path limit no longer applies.
- **Removed.** The pickle results files, `Mesh`, `convertToOneBlade/MultiBlade` (the output is
  always `(beta, blade)`) and `computeBEM`.

With these two bugs patched in the old code, `tests/test_regression_old.py` checks that the new
model reproduces it to `rtol=1e-10`.
