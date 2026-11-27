"""WAPE propagation model.

Workflow::

    cfg = PropagationConfig.from_yaml("input.yaml")   # ``propagation`` section
    prepare_cases(cfg)                                # input.nml + 2D_WAPE in t<i>/<h>/, run.sh
    run_local(cfg.work_path)                          # run all cases
    dl = DeltaL.from_cases(cfg, isource=0)            # deltaL(r, tau, z, freq, source_height)
    dl.save(cfg.output_file(0))                       # dl.ds is the xr.Dataset

Coupling with the source: :func:`wfn.spl.combine` interpolates ``deltaL`` at the
blade segment heights ``source_ds.z_seg`` and applies it to the segment PSDs.
"""
from .cases import case_dir, case_inputs, local_domain, prepare_cases
from .config import (
    DomainConfig,
    EwapeInputs,
    PiecewapeInputs,
    PropagationConfig,
    RotorConfig,
    WapeInputs,
    discretize_rotor,
)
from .launcher import run_local, write_run_script
from .namelist import namelist_text, read_namelist, write_namelist
from .results import DeltaL, PeResults, assemble_polar, loop_angle, read_wape, refine_angle, third_octave, tl_to_deltaL, to_grid, wape_filename

__all__ = [
    "PropagationConfig",
    "WapeInputs",
    "PiecewapeInputs",
    "EwapeInputs",
    "PeResults",
    "DeltaL",
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
    "tl_to_deltaL",
]
