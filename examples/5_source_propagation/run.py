"""Add the propagation effects (deltaL) to the free-field source result.

Uses the results of the previous examples:
* ``examples/1_source_nrel/results/nrel_0.nc``         (run ``1_source_nrel/run_plane.py``)
* ``examples/3_propagation_simple/results/deltaL_0.nc`` (run ``3_propagation_simple/run.py``)

The propagation example only computes 200 Hz, so the source is reduced to that frequency
(``combine`` raises an error if a source frequency is missing from deltaL).
"""
import matplotlib.pyplot as plt
import xarray as xr

from wfn import spl, units
from wfn.propagation import DeltaL

src = xr.open_dataset("../1_source_nrel/results/nrel_0.nc").sel(freq=[200.0])
dl = DeltaL.load("../3_propagation_simple/results/deltaL_0.nc")

# deltaL (polar) is interpolated onto the source grid, then onto the height of each segment
out = spl.combine(src, dl)
out.to_netcdf("spl_0.nc")

# rotor total, energetic average over the rotation
free = units.db_mean(units.db_sum(units.spl(src), ("blade", "seg")), "beta")
prop = units.db_mean(units.db_sum(units.spl(out), ("blade", "seg")), "beta")

fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharey=True)
kw = dict(x="x", y="y")
free.isel(z=0, freq=0).plot(ax=axes[0], cmap="viridis", **kw)
prop.isel(z=0, freq=0).plot(ax=axes[1], cmap="viridis", **kw)
(prop - free).isel(z=0, freq=0).plot(ax=axes[2], cmap="RdBu_r", vmin=-15, vmax=15, **kw)
for ax, title in zip(axes, ["free field SPL", "SPL with propagation", "difference"]):
    ax.set_aspect("equal")
    ax.set_title(f"{title}, f = 200 Hz")
plt.tight_layout()
plt.savefig("spl_plane.png", dpi=150)
print("saved spl_0.nc and spl_plane.png")
