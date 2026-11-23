"""NREL 5MW noise spectrum at a single receiver (new version of
old_project/examples/3_wind_turbine_aero/NREL/3_compute_noise_at_one_receiver.py)."""
import time

import matplotlib.pyplot as plt
import numpy as np

from wfn import receivers, units
from wfn.source import WindTurbineNoise

# turbine, airfoil database and flow from the shared input file
src = WindTurbineNoise.from_yaml("input.yaml")

xR, yR, zR = 1000.0, 600.0, 90.0
t0 = time.perf_counter()
ds = src.compute(receivers.points([xR], [yR], [zR]))
print(f"computed in {time.perf_counter() - t0:.1f} s")
ds.to_netcdf("single_receiver.nc")

# rotor maps (blade 0 over beta, all blades -> full turn)
rotor = ds.isel(point=0)
az = rotor.azimuth.values.ravel()
order = np.argsort(az)


def full_turn(da):
    return da.stack(k=("beta", "blade")).transpose("k", "seg").values[order]


Y = full_turn(rotor.y_seg)
Z = full_turn(rotor.z_seg)
swl_seg = units.swl(rotor).isel(freq=0)
fig, axes = plt.subplots(1, 3, figsize=(10, 3.5))
for ax, (name, data) in zip(axes, [("U_rel (m/s)", full_turn(rotor.U_rel)),
                                   ("AoA (deg)", full_turn(rotor.AoA)),
                                   (f"SWL {ds.freq.values[0]:.0f} Hz (dB)", full_turn(swl_seg))]):
    m = ax.scatter(Y, Z, c=data, s=12)
    ax.set_aspect("equal")
    ax.set_title(name)
    ax.axis("off")
    fig.colorbar(m, ax=ax, shrink=0.8)
fig.tight_layout()
fig.savefig("polars_NREL.png", dpi=150)

# spectra averaged over one rotation (energetic mean over beta, sum over blades/segments)
def rotor_level(level):
    return units.db_mean(units.db_sum(level, ["blade", "seg"]), "beta").isel(point=0)


f = ds.freq.values
aw = units.a_weighting(f)
plt.figure()
plt.plot(f, rotor_level(units.swl(ds)) + aw, "gray", label="SWL_tot")
for c in ("TIN", "TEN"):
    plt.plot(f, rotor_level(units.swl(ds, c)) + aw, label=c)
plt.xscale("log")
plt.xlabel("f (Hz)")
plt.ylabel("SWL (dBA)")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig("spectrum_NREL.png", dpi=150)
print("saved polars_NREL.png, spectrum_NREL.png, single_receiver.nc")
