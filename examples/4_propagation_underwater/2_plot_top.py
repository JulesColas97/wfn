
"""WAPE runs around one source (no flow, flat impedance ground), then deltaL on the receiver plane.

Needs the compiled WAPE code: set ``propagation.executable`` in input.yaml or the
environment variable ``WFN_WAPE``.
"""
import matplotlib.pyplot as plt
from wfn.propagation import DeltaL

dl = DeltaL.load("./results/deltaL_0.nc")

fig, axes = plt.subplots(1, dl.sizes["source_height"], figsize=(10, 4), sharey=True)

for ax, h in zip(axes, dl.source_height.values):
    dl.deltaL.sel(source_height=h).sel(z=-10).isel(freq=0).plot(
        ax=ax, x="x", y="y", cmap="RdBu_r", vmin=-100, vmax=-20, cbar_kwargs={"label": "deltaL (dB)"}
    )
    ax.set_aspect("equal")
    ax.set_title(f"source at {h:.0f} m, f = {dl.freq.values[0]:.0f} Hz")
plt.tight_layout()
# plt.savefig("deltaL_plane.png", dpi=150)
plt.show()
