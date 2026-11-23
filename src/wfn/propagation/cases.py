"""Prepare one WAPE run folder per (source position, source height).

Folder layout (same as the old code)::

    <work_dir>/t<isource>/<int(height)>/input.nml
    <work_dir>/t<isource>/<int(height)>/2D_WAPE
    <work_dir>/t<isource>/<int(height)>/<case_name>_<tau>.h5   (written by WAPE)

In every run the source is at the origin of the WAPE domain.
"""
from __future__ import annotations

import dataclasses
import logging
import shutil
from pathlib import Path

from .config import PropagationConfig, WapeInputs
from .launcher import write_run_script
from .namelist import write_namelist

logger = logging.getLogger(__name__)


def case_dir(cfg: PropagationConfig, isource: int, height: float) -> Path:
    """Run folder of one (source, height) case."""
    return cfg.work_path / f"t{isource}" / f"{int(height)}"


def local_domain(cfg: PropagationConfig, isource: int) -> dict:
    """WAPE domain extents around source ``isource`` (empty dict if no ``domain`` is set)."""
    if cfg.domain is None:
        return {}
    d = cfg.domain
    xs, ys = cfg.source_positions()[isource][:2]
    return dict(
        Lx1=d.xmax - xs,
        Lx2=xs - d.xmin,
        Ly1=d.ymax - ys,
        Ly2=ys - d.ymin,
        Lz=d.zmax,
        Xs_flow=xs,
        Ys_flow=ys,
        **(
            dict(Xs_wave=xs, Ys_wave=ys)
            if cfg.backend == "ewape"
            else dict(Xs_terrain=xs, Ys_terrain=ys)
        ),
    )


def case_inputs(cfg: PropagationConfig, isource: int, height: float) -> WapeInputs:
    """Namelist parameters of one case.

    The resolved angles are written as an explicit ``tau`` list and ``taumin/taumax/nb_tau`` are
    cleared: the Fortran codes give priority to them, and would otherwise ignore the list.
    """
    return dataclasses.replace(
        cfg.wape,
        **local_domain(cfg, isource),
        source_height=height,
        tau=cfg.angles(),
        taumin=None,
        taumax=None,
        nb_tau=None,
        frequencies=cfg.freqs(),
    )


def prepare_cases(cfg: PropagationConfig) -> list[Path]:
    """Write ``input.nml`` and copy the executable in every run folder, then write ``run.sh``.

    Returns:
        the list of run folders.
    """
    exe = cfg.executable_path()
    dirs = []
    for isource in range(len(cfg.source_positions())):
        for height in cfg.heights():
            path = case_dir(cfg, isource, height)
            path.mkdir(parents=True, exist_ok=True)
            write_namelist(case_inputs(cfg, isource, height), path / "input.nml", cfg.backend)
            shutil.copy(exe, path / "2D_WAPE")
            dirs.append(path)
    write_run_script(dirs, cfg.work_path)
    logger.info("prepared %d WAPE runs in %s", len(dirs), cfg.work_path)
    return dirs
