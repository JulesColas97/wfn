"""WAPE propagation model.

Workflow::

    cfg = PropagationConfig.from_yaml("input.yaml")   # ``propagation`` section
    prepare_cases(cfg)                                # input.nml + 2D_WAPE in t<i>/<h>/, run.sh
    run_local(cfg.work_path)                          # run all cases
    ds = assemble_polar(cfg, isource=0)               # deltaL(r, tau, z, freq, source_height)
    ds.to_netcdf(cfg.output_file(0))

Coupling with the source is planned as interpolating ``deltaL`` at the blade
segment heights ``source_ds.z_seg``.
"""
from .cases import case_dir, case_inputs, local_domain, prepare_cases
from .config import DomainConfig, PropagationConfig, RotorConfig, WapeInputs, discretize_rotor
from .launcher import run_local, write_run_script
from .namelist import namelist_text, read_namelist, write_namelist
from .results import assemble_polar, loop_angle, read_wape, refine_angle, third_octave, to_grid, wape_filename

__all__ = [
    "PropagationConfig",
    "WapeInputs",
    "DomainConfig",
    "RotorConfig",
    "discretize_rotor",
    "namelist_text",
    "write_namelist",
    "read_namelist",
    "case_dir",
    "case_inputs",
    "local_domain",
    "prepare_cases",
    "write_run_script",
    "run_local",
    "read_wape",
    "wape_filename",
    "assemble_polar",
    "loop_angle",
    "refine_angle",
    "to_grid",
    "third_octave",
]
