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

import numpy as np
from wfn import spl, units
from wfn.propagation import DeltaL
IFREQ = 1
src = xr.open_dataset("./results/MySe3_0.nc")
dl = DeltaL.load("./results/deltaL_side_0.nc")
dl = dl.to_deltaL()
# dl.deltaL.sel(y=0,freq=100, source_height=90).plot(x="x", y="z", 
#                                                    cmap="RdBu_r",vmin=-10,vmax=10)
# plt.gca().set_aspect("equal")
# plt.show()

# deltaL (polar) is interpolated onto the source grid, then onto the height of each segment
# -----------------------------------------------------------------------------
out = spl.combine(src, dl)
out.to_netcdf("results/spl_0.nc")

def one_blade(ds):
    """(beta, blade) -> a single blade doing a full turn (dim beta of size n_beta * n_blade)."""
    nb = ds.sizes["blade"]
    st = ds.drop_vars("azimuth").stack(rot=("blade", "beta"))
    az = st.beta.values + 2 * np.pi * st.blade.values / nb
    st = (st.drop_vars(["rot", "blade", "beta"])
            .rename(rot="beta")
            .assign_coords(beta=("beta", az, {"units": "rad"})))
    return st.transpose(..., "beta", "seg", "freq")

one = one_blade(src)                                   # Spp_TEN: (x, y, z, beta=36, seg, freq)



# rotor total, energetic average over the rotation
# -----------------------------------------------------------------------------
free = units.db_mean(units.db_sum(units.spl(src), ("blade", "seg")), "beta")
prop = units.db_mean(units.db_sum(units.spl(out), ("blade", "seg")), "beta")
# SPL_beta = units.db_sum(units.spl(src), ("blade", "seg"))
SPL_beta = units.db_sum(units.spl(one), ("seg"))
AM = SPL_beta.max("beta") - SPL_beta.min("beta")

# lvl = units.db_sum(units.spl(one), "seg")              # one blade over 0..2π
# AM = lvl.max("beta") - lvl.min("beta")
# Plot 
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(3,1, figsize=(5,5), sharey=True)
kw = dict(x="x", y="z")
free.isel(y=0, freq=IFREQ).plot(ax=axes[0], cmap="viridis",vmin=30,vmax=80, **kw)
prop.isel(y=0, freq=IFREQ).plot(ax=axes[1], cmap="viridis",vmin=30,vmax=80, **kw)
AM.isel(y=0, freq=IFREQ).plot(ax=axes[2], cmap="inferno", vmin=0, vmax=10, **kw)
for ax, title in zip(axes, ["free field SPL", "SPL with propagation", "AM"]):
    ax.set_aspect("equal")
    ax.set_title(f"{title}, f = {int(out.freq[IFREQ])} Hz")
plt.tight_layout()
plt.savefig("spl_side.png", dpi=300)
# print("saved spl_0.nc and spl_plane.png")


# animate
lvl = SPL_beta.isel(y=0).sel(freq=100, method="nearest")    # (x, y, beta)

vmin, vmax = float(lvl.min()), float(lvl.max())        # fixed color scale for all frames
fig, ax = plt.subplots(figsize=(8,4))
im = ax.pcolormesh(lvl.x, lvl.z, lvl.isel(beta=0).T, vmin=30, vmax=80, cmap="viridis")
ax.set_aspect("equal")
fig.colorbar(im, label="SPL (dB)")
title = ax.set_title("")
ax.set_xlabel("x (m)")
ax.set_ylabel("z (m)")

plt.tight_layout()

def update(i):
    im.set_array(lvl.isel(beta=i).T.values.ravel())
    title.set_text(f"beta = {float(lvl.beta[i]):.2f} rad")
    return im, title

anim = FuncAnimation(fig, update, frames=lvl.sizes["beta"], interval=100, blit=False)
anim.save("spl_ff_1_blade.gif", writer=PillowWriter(fps=10))   # loops forever by default
plt.show()
