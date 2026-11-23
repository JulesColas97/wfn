"""Run XFoil to compute boundary layer quantities near the trailing edge.

Port of ``old_project/code/prepost/source/xfoil.py``. The XFoil script and the
post-processing (spline interpolation of the BL quantities at ``target_frac``
of the chord, pressure gradient from Cp) are unchanged. Differences:

* every run happens in its own working directory, the airfoil file is copied
  there under a short name, so the 64 characters path limit of XFoil no longer
  matters;
* XFoil is called without a shell and the virtual frame buffer (Xvfb) is a
  managed subprocess that is always terminated (:class:`VirtualDisplay`).
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.interpolate as interp

logger = logging.getLogger(__name__)

XFOIL_EXE = "xfoil"
GAMMA = 1.4
R_AIR = 287.0


class XfoilError(RuntimeError):
    pass


class VirtualDisplay:
    """Context manager starting ``Xvfb`` so that XFoil does not open windows."""

    def __init__(self, display_id: int = 99):
        self.display_id = display_id
        self.proc: subprocess.Popen | None = None

    def __enter__(self) -> "VirtualDisplay":
        if shutil.which("Xvfb") is None:
            raise XfoilError("Xvfb not found, install it or use virtual_buffer=False")
        lock = Path(f"/tmp/.X{self.display_id}-lock")
        if lock.exists():
            lock.unlink()
        self.proc = subprocess.Popen(
            ["Xvfb", f":{self.display_id}"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        time.sleep(0.5)
        return self

    @property
    def env(self) -> dict:
        return {**os.environ, "DISPLAY": f":{self.display_id}"}

    def __exit__(self, *exc):
        if self.proc is not None:
            self.proc.terminate()
            try:
                self.proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.proc.kill()
        return False


# -----------------------------------------------------------------------------
# input file
# -----------------------------------------------------------------------------
def write_input_file(workdir: Path, inputs: dict, foil_file: str = "foil.dat") -> Path:
    """Write the XFoil command file. All paths are relative to ``workdir``.

    ``inputs`` keys: ``chord``, ``re``, ``mach``, ``aoa``, ``trip_top``,
    ``trip_bot`` (trips may be None).
    """
    lines = ["PLOP", "W 0.05", "", f"LOAD {foil_file}"]
    lines += ["GDES", "GSET", "SCAL", f"{inputs['chord']}", "EXEC", "GSET", ""]
    lines += ["PPAR", "N 200", "", "", "OPER"]
    if inputs.get("trip_top") is not None and inputs.get("trip_bot") is not None:
        lines += ["VPAR", "XTR", f"{inputs['trip_top']}", f"{inputs['trip_bot']}", ""]
    else:
        lines += ["VPAR", ""]
    lines += ["ITER", "800"]
    lines += [f"VISC", f"{inputs['re']:.3f}", "M", f"{inputs['mach']:.5f}"]
    lines += ["pacc", "polar.dat", ""]
    for ii, aoa in enumerate(inputs["aoa"]):
        lines += [f"ALFA {aoa:.3f}", f"DUMP bl_{ii}.dat", f"CPWRT Cp_{ii}.dat"]
    lines += ["", "pacc", "QUIT", ""]
    fname = workdir / "xfoil_input.txt"
    fname.write_text("\n".join(lines))
    return fname


def run_xfoil_process(workdir: Path, env: dict | None = None, timeout: float = 300) -> None:
    """Run XFoil on ``workdir/xfoil_input.txt`` (stdout in ``out.out``)."""
    if shutil.which(XFOIL_EXE) is None:
        raise XfoilError(f"'{XFOIL_EXE}' executable not found")
    with open(workdir / "xfoil_input.txt") as fin, open(workdir / "out.out", "w") as fout:
        try:
            res = subprocess.run(
                [XFOIL_EXE], stdin=fin, stdout=fout, stderr=subprocess.PIPE, text=True,
                cwd=workdir, env=env, timeout=timeout,
            )
        except subprocess.TimeoutExpired as e:
            raise XfoilError(f"XFoil timed out after {timeout} s") from e
    if res.returncode != 0:
        if "SIGFPE" in (res.stderr or ""):
            raise XfoilError("floating-point exception in XFoil")
        raise XfoilError(f"XFoil failed with exit code {res.returncode}: {res.stderr}")


# -----------------------------------------------------------------------------
# post-processing (unchanged from the old code)
# -----------------------------------------------------------------------------
def c_round(value: float) -> float:
    """Round to 6 decimals."""
    return np.round(value * 10**6) / 10**6


def _suction_pressure_splines(x: np.ndarray, data: np.ndarray, le_index: int) -> tuple:
    suction = interp.splrep(x[1:le_index][::-1], data[1:le_index][::-1], k=2)
    pressure = interp.splrep(x[le_index + 1:], data[le_index + 1:], k=2)
    return suction, pressure


def read_bl_data(workdir: Path, aoa_index: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read the BL dump and Cp files of one angle of attack (``*****`` are NaN)."""
    colspecs = [(1, 10), (10, 19), (19, 28), (28, 37), (37, 47), (47, 57), (57, 67)]
    bl = pd.read_fwf(workdir / f"bl_{aoa_index}.dat", colspecs=colspecs, skiprows=1,
                     header=None, na_values=["**********"])
    bl.columns = ["s", "x", "y", "Ue", "Dstar", "theta", "cf"]
    cp = pd.read_csv(workdir / f"Cp_{aoa_index}.dat", sep=r"\s+", comment="#", header=None,
                     skiprows=1)
    cp.columns = ["x", "cp"]
    return bl, cp


def _bl_splines(bl: pd.DataFrame, cp: pd.DataFrame) -> dict:
    x = np.array(bl["x"])
    le_index = int(np.argmin(x))
    out = {}
    for key in ["Ue", "Dstar", "theta", "cf"]:
        out[f"{key}_suction"], out[f"{key}_pressure"] = _suction_pressure_splines(
            x, np.array(bl[key]), le_index)
    out["cp_suction"], out["cp_pressure"] = _suction_pressure_splines(
        np.array(cp["x"]), np.array(cp["cp"]), le_index)
    return out


def _at(pos: float, tck) -> float:
    return c_round(interp.splev(pos, tck, der=0))


def _empty_output(aoa: np.ndarray) -> dict:
    n = len(aoa)
    out = {"aoa": np.array(aoa)}
    for side in ("suction", "pressure"):
        for key in ("Dstar", "theta", "Ue", "dpdx", "cf", "dpdx1"):
            out[f"{key}_{side}"] = np.zeros(n)
    out["cl"] = np.full(n, np.nan)
    out["cd"] = np.full(n, np.nan)
    return out


def _fill_aoa(iaoa: int, target_frac: float, chord: float, splines: dict, out: dict,
              dyn_press: float) -> None:
    pos = target_frac * chord
    for key in ["Ue", "Dstar", "theta", "cf"]:
        out[f"{key}_suction"][iaoa] = _at(pos, splines[f"{key}_suction"])
        out[f"{key}_pressure"][iaoa] = _at(pos, splines[f"{key}_pressure"])
    pos_m = (target_frac - 0.05) * chord
    for side in ("suction", "pressure"):
        cp_p1 = _at(pos, splines[f"cp_{side}"])
        cp_m1 = _at(pos_m, splines[f"cp_{side}"])
        out[f"dpdx1_{side}"][iaoa] = cp_p1 * dyn_press
        out[f"dpdx_{side}"][iaoa] = (cp_p1 - cp_m1) / (0.05 * chord) * dyn_press


def _read_polar(workdir: Path, out: dict) -> None:
    polar = pd.read_csv(workdir / "polar.dat", sep=r"\s+", skiprows=11, usecols=[0, 1, 2])
    polar.columns = ["aoa", "cl", "cd"]
    polar = polar.drop_duplicates(subset=["aoa"])
    merged = pd.merge(pd.DataFrame({"aoa": out["aoa"]}), polar, on="aoa", how="left")
    out["cl"] = merged["cl"].to_numpy()
    out["cd"] = merged["cd"].to_numpy()


def post_process(workdir: Path, inputs: dict, n_aoa: int | None = None) -> dict:
    """Interpolate the BL quantities at ``target_frac * chord`` for each AoA."""
    T = inputs["T"]
    dyn_press = 0.5 * inputs["rho"] * (inputs["mach"] * (GAMMA * R_AIR * T) ** 0.5) ** 2
    out = _empty_output(inputs["aoa"])
    n_aoa = len(inputs["aoa"]) if n_aoa is None else n_aoa
    for ii in range(n_aoa):
        bl, cp = read_bl_data(workdir, ii)
        _fill_aoa(ii, inputs["target_frac"], inputs["chord"], _bl_splines(bl, cp), out, dyn_press)
    _read_polar(workdir, out)
    return out


def run_xfoil(inputs: dict, foil_path: str | Path, env: dict | None = None,
              workdir: str | Path | None = None, keep_files: bool = False) -> dict:
    """Run XFoil for one airfoil, one Reynolds/Mach number and a list of AoA.

    Args:
        inputs: ``chord, re, mach, aoa, trip_top, trip_bot, target_frac, rho, T``.
        foil_path: airfoil geometry file (XFoil format).
        env: environment for XFoil (e.g. ``VirtualDisplay.env``).
        workdir: directory for XFoil files (a temporary one if None).
        keep_files: keep the XFoil files (only for a temporary workdir).

    Returns:
        dict with ``Dstar_*, theta_*, Ue_*, cf_*, dpdx_*`` for ``suction`` and
        ``pressure`` sides and ``cl``, ``cd``, one value per AoA. If XFoil
        crashes, only the AoA computed before the crash are post-processed.
    """
    tmp = None
    if workdir is None:
        tmp = tempfile.mkdtemp(prefix="wfn_xfoil_")
        workdir = Path(tmp)
    else:
        workdir = Path(workdir)
        if workdir.exists():
            shutil.rmtree(workdir)
        workdir.mkdir(parents=True)
    try:
        shutil.copy(foil_path, workdir / "foil.dat")
        write_input_file(workdir, inputs)
        n_aoa = None
        try:
            run_xfoil_process(workdir, env=env)
        except XfoilError as e:
            done = [int(f.stem.split("_")[-1]) for f in workdir.glob("Cp_*.dat")]
            if not done:
                raise
            n_aoa = max(done)  # the last dump may be incomplete
            logger.warning("%s; keeping AoA up to %s", e, inputs["aoa"][n_aoa - 1] if n_aoa else None)
        return post_process(workdir, inputs, n_aoa=n_aoa)
    finally:
        if tmp is not None and not keep_files:
            shutil.rmtree(tmp, ignore_errors=True)
