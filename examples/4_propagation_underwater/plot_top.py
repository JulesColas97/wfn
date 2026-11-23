
"""WAPE runs around one source (no flow, flat impedance ground), then deltaL on the receiver plane.

Needs the compiled WAPE code: set ``propagation.executable`` in input.yaml or the
environment variable ``WFN_WAPE``.
"""
import xarray as xr
import matplotlib.pyplot as plt

from wfn.propagation import PropagationConfig, assemble_polar, refine_angle, loop_angle, prepare_cases, run_local, to_grid
from wfn.receivers import ReceiversConfig

ds = xr.open_dataset("./results/deltaL_0.nc")   # deltaL(r, tau, z, freq, source_height)
ds = refine_angle(ds)
ds = refine_angle(ds)

plane = to_grid(loop_angle(ds), ReceiversConfig.from_yaml("input.yaml").build())



fig, axes = plt.subplots(1, plane.sizes["source_height"], figsize=(10, 4), sharey=True)

for ax, h in zip(axes, plane.source_height.values):
    plane.deltaL.sel(source_height=h).sel(z=-15).isel(freq=0).plot(
        ax=ax, x="x", y="y", cmap="RdBu_r", vmin=-80, vmax=-20, cbar_kwargs={"label": "deltaL (dB)"}
    )
    ax.set_aspect("equal")
    ax.set_title(f"source at {h:.0f} m, f = {plane.freq.values[0]:.0f} Hz")
plt.tight_layout()
# plt.savefig("deltaL_plane.png", dpi=150)
plt.show()
