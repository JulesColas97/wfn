# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.


# Quickfix errors

When I ask about an error, a build failure, or "the quickfix list",
read `.qf.txt` in the project root first. It contains the latest
Neovim quickfix list (from :make, LSP, grep, etc.) in the format
`file:line:col: message`.

You can also read the output and error files as err.err and out.out in the folder
corresponding to he quickfix list folder.

Do not edit files or run commands. I only want explanations of where
the error comes from and how to fix it.

Dont propose a plan instead asked explicetely. 
Directly answer the questions with code snippet inside the chat
Dont try to run tests instead explicetely asked to.

When you propose code try to :
    - give a simple code snipet and dont try to avoid complex multi functions code 


## What this is
`wfn` (wind farm noise) is a rewrite of the legacy package `prepost` (in `old_project/code/prepost`, which is gitignored and present only locally). The goals are listed in `changes.md`: decoupled modules, one shared YAML input with one section per module, and xarray outputs instead of the old `Mesh` and pickle files. `wfn.source` is fully ported. `wfn.flow` has only `Flow` and `UniformFlow`. `wfn.propagation` (WAPE) is ported with the local launcher and fixed-grid assembly only (no Slurm launcher, no variable-grid/pe_model==5 mode yet). Next come the source/propagation coupling and the other flows (Walberla, Meso-NH, FarmShadow).

## Commands
```bash
source venv/bin/activate            # local venv is ./venv (README says .venv)
pip install -e ".[dev]"             # extras: plot, explore (Dash explorer)
pytest                              # all tests
pytest tests/test_source_units.py::test_name -q   # single test
pytest -m "not legacy"              # skip comparisons against old_project/
pytest -m "not xfoil"               # skip tests needing the xfoil executable
python scripts/convert_legacy_airfoil.py old.dat new.nc
```
There is no linter or formatter configured.

## Architecture
- **One YAML file, one section per module.** `wfn/config.py` knows nothing about the sections. It only provides `load_yaml`, `get_section`, and `structure`/`unstructure`, which use cattrs with `forbid_extra_keys=True`, so unknown keys raise an error. Each module owns its own dataclass config with `from_dict`/`from_yaml`, for example `receivers.ReceiversConfig` and `source.config.SourceConfig`. Paths are resolved relative to the input file with `resolve_path`. Do not add a global config class, and keep `source` and `propagation` from importing each other's configs.
- **Receivers** (`receivers.py`) are plain `xr.Dataset` grids (`cartesian`, `polar`, `points`) that replace the old `Mesh`.
- **Flow** (`flow/base.py`) is an abstract interface. Subclasses register themselves with `class X(Flow, type_name="...")`, and the YAML `flow.type` picks the class from that registry. The source model needs only `at_points(x, y, z) -> {u, v, epsilon}` and optionally `turbine_positions`. Propagation will need full fields through `to_dataset`.
- **Source** (`source/`): `model.WindTurbineNoise` orchestrates the computation. It uses `turbine.WindTurbine` (geometry/controller, with defaults in `source/data/turbines/*.yaml`), `airfoil.AirfoilDatabase` (BL data on `(uinf, aoa, seg)` stored as `.nc`; built with XFoil from `xfoil.py` or imported from the legacy `.dat`), `geometry.py` (rotor kinematics), `amiet.py` (TEN/TIN per segment) and the numba kernels in `kernels.py`. `compute()` returns an `xr.Dataset` with `Spp_TEN`/`Spp_TIN` on dims `(*grid, beta, blade, seg, freq)`. These are linear PSDs, so you sum over `blade`/`seg` before converting to dB with `units.py`. `units.spl` is a level per Hz and must not be summed over `freq`; for an OASPL use `units.spl_third_octave` (band levels) then `db_sum(..., "freq")`. It also returns the segment variables (`x_seg, y_seg, z_seg, U_inf, AoA, ...`) on `(beta, blade, seg)`. Coupling with propagation is in `spl.py` (ΔL interpolated at `z_seg`).
- **Propagation** (`propagation/`): drives the external Fortran WAPE code (path in `propagation.executable` or env `WFN_WAPE`; the Fortran source lives in `old_project/code/wape`). `config.WapeInputs` (shared keys) with subclasses `PiecewapeInputs`/`EwapeInputs` (registered by `backend=`; `PropagationConfig.from_dict` picks the class from `propagation.backend`, so a key of the other backend raises). Field names are the Fortran namelist keys, so do not rename them. `cases.prepare_cases` writes `<work_dir>/t<isource>/<int(height)>/input.nml` + `run.sh`, and `launcher.run_local` runs them. `results.read_wape` reads one `<case>_<tau:04d>.h5`. `results.PeResults` (one file) and `results.DeltaL` (polar result, `.ds` is the Dataset; `from_cases`, `loop_angle`, `refine_angle`, `to_grid`, `save/load`) wrap the free functions, which stay available. `assemble_polar` builds `deltaL(r, tau, z, freq, source_height)` + `terrain(r, tau)` on a `receivers.polar` grid, and `to_grid` interpolates it onto any receiver grid (call `loop_angle` first). Legacy checks are in `tests/test_regression_propagation.py`.
- **Coupling** (`spl.py`): `combine(src_ds, deltaL, grid=None)` brings ΔL onto the source grid (or both onto `grid`; source interpolation is cartesian only), checks that the frequencies match, interpolates ΔL in energy over `source_height` at `z_seg`, and returns the source Dataset with `Spp_* *= 10**(ΔL/10)` plus `deltaL` per segment. The old `combine_linear` weight swap is intentionally fixed (see README).
- **Blade convention:** rotation is clockwise seen from upwind (`y = -sin(beta) r`). Yaw is applied consistently to the positions, the Doppler angles and the frame transforms. The README lists the old-code bugs that were fixed on purpose (yaw_rad returned omega, mirrored blade positions). Do not "fix" the new code back to match the old one.

## Testing against the legacy code
`tests/test_regression_old.py` (marker `legacy`) imports the old `prepost` through a stubbed package in `tests/conftest.py::old_prepost`. That fixture avoids GUI imports, patches the old yaw bug and checks the new model to `rtol=1e-10`. Any change to the source physics must keep this test passing, unless the change is an intentional, documented deviation. Work in small steps that can be checked this way.

Note: the README still points to `examples/source_nrel/`, but the example folders are now numbered (`examples/0_airfoil_nrel`, `1_source_nrel`, `2_airfoil_from_data`, `3_propagation_simple`).
