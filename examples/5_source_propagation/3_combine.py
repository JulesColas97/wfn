"""Add the propagation effects (deltaL) to the free-field source result.

Uses the results of the previous examples:
* ``examples/1_source_nrel/results/nrel_0.nc``         (run ``1_source_nrel/run_plane.py``)
* ``examples/3_propagation_simple/results/deltaL_0.nc`` (run ``3_propagation_simple/run.py``)

The propagation example only computes 200 Hz, so the source is reduced to that frequency
(``combine`` raises an error if a source frequency is missing from deltaL).
"""
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.animation import FuncAnimation, PillowWriter

from wfn import spl, units
from wfn.propagation import DeltaL
IFREQ = 1
src = xr.open_dataset("./results/MySe3_0.nc")
dl = DeltaL.load("./results/deltaL_0.nc")

# deltaL (polar) is interpolated onto the source grid, then onto the height of each segment
# -----------------------------------------------------------------------------
out = spl.combine(src, dl)
out.to_netcdf("results/spl_0.nc")

# rotor total, energetic average over the rotation
# -----------------------------------------------------------------------------
free = units.db_mean(units.db_sum(units.spl(src), ("blade", "seg")), "beta")
prop = units.db_mean(units.db_sum(units.spl(out), ("blade", "seg")), "beta")
SPL_beta = units.db_sum(units.spl(src), ("blade", "seg"))
AM = SPL_beta.max("beta") - SPL_beta.min("beta")

# Plot 
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharey=True)
kw = dict(x="x", y="y")
free.isel(z=0, freq=IFREQ).plot(ax=axes[0], cmap="viridis", **kw)
prop.isel(z=0, freq=IFREQ).plot(ax=axes[1], cmap="viridis", **kw)
AM.isel(z=0, freq=IFREQ).plot(ax=axes[2], cmap="inferno", vmin=0, vmax=10, **kw)
for ax, title in zip(axes, ["free field SPL", "SPL with propagation", "AM"]):
    ax.set_aspect("equal")
    ax.set_title(f"{title}, f = {int(out.freq[IFREQ])} Hz")
plt.tight_layout()
# plt.show()
# plt.savefig("spl_plane.png", dpi=150)
# print("saved spl_0.nc and spl_plane.png")


# animate
lvl = SPL_beta.isel(z=0).sel(freq=100, method="nearest")    # (x, y, beta)

vmin, vmax = float(lvl.min()), float(lvl.max())        # fixed color scale for all frames
fig, ax = plt.subplots(figsize=(6, 5))
im = ax.pcolormesh(lvl.x, lvl.y, lvl.isel(beta=0).T, vmin=vmin, vmax=vmax, cmap="viridis")
ax.set_aspect("equal")
fig.colorbar(im, label="SPL (dB)")
title = ax.set_title("")

def update(i):
    im.set_array(lvl.isel(beta=i).T.values.ravel())
    title.set_text(f"beta = {float(lvl.beta[i]):.2f} rad")
    return im, title

anim = FuncAnimation(fig, update, frames=lvl.sizes["beta"], interval=100, blit=False)
anim.save("spl_beta.gif", writer=PillowWriter(fps=10))   # loops forever by default
plt.show()
