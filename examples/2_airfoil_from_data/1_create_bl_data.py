"""Create the boundary layer database of the NREL 5MW blade with XFoil.

Port of old_project/examples/3_wind_turbine_aero/NREL/0_create_bl_data.py.

    python 0_create_bl_data.py           # 14 segments x 4 velocities x 41 AoA (a few minutes)
    python 0_create_bl_data.py --quick   # 2 segments x 2 velocities x 9 AoA (quick check)

Needs ``xfoil`` (and ``Xvfb`` to hide the XFoil windows).
"""
import logging
import shutil
import sys
import time

import numpy as np

from wfn.source import AirfoilDatabase, WindTurbine

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

wt = WindTurbine.from_yaml("input.yaml")

# Velocity range seen by the blade segments over the operating range.
# NB: the old script used rpm * pi / 60 for omega, i.e. half the real rotation speed.
u_cutin, u_cutout = 2.0, 14.0
omega_min, omega_max = 4 * 2 * np.pi / 60, 12.1 * 2 * np.pi / 60   # rad/s
u_min = np.hypot(u_cutin, wt.seg[0] * omega_min)
u_max = np.hypot(u_cutout, wt.seg[-1] * omega_max)
print(f"relative velocity range: {u_min:.1f} - {u_max:.1f} m/s")

# AoA from 0 to -20 then from 0 to 20: XFoil starts each point from the previous
# converged solution, so moving away from 0 helps convergence. Duplicates
# (0 twice) are removed and AoA sorted by the database.
aoa = np.concatenate((np.arange(0, -10, -1), np.arange(0, 10, 1)))
uinf = np.linspace(u_min, u_max, 4)
output = "MySE3_13seg_xfoil.nc"

t0 = time.perf_counter()
db = AirfoilDatabase.for_turbine(
    wt, uinf, aoa,
    target_frac=0.95,                             # BL quantities at 95 % of the chord
    virtual_buffer=False,#shutil.which("Xvfb") is not None,
)
db.save(output)
print(f"done in {time.perf_counter() - t0:.0f} s -> {output}")
print(db.ds)
