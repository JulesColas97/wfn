"""WAPE runs around one source (no flow, flat impedance ground), then deltaL on the receiver plane.

Needs the compiled WAPE code: set ``propagation.executable`` in input.yaml or the
environment variable ``WFN_WAPE``.
"""
import matplotlib.pyplot as plt
import numpy as np

from wfn.propagation import DeltaL, PropagationConfig
from wfn.receivers import ReceiversConfig

# cfg = PropagationConfig.from_yaml("input_top.yaml")
dl = DeltaL.from_yaml("input_top.yaml")
dl = dl.refine_angle()
dl = dl.refine_angle()
dl = dl.loop_angle()
dl = dl.to_grid() #ReceiversConfig.from_yaml("input_top.yaml").build())
dl.save(dl.cfg.output_file(0))

fig, axes = plt.subplots(1, dl.sizes["source_height"], figsize=(10, 4), sharey=True)
for ax, h in zip(axes, dl.source_height.values):
    dl.sel(source_height=h).isel(z=0, freq=0).deltaL.plot(
        ax=ax, x="x", y="y", cmap="RdBu_r", vmin=-10, vmax=10, cbar_kwargs={"label": "deltaL (dB)"}
    )
    ax.set_aspect("equal")
    ax.set_title(f"source at {h:.0f} m, f = {dl.freq.values[0]:.0f} Hz")
plt.tight_layout()
plt.show()
plt.savefig("deltaL_plane.png", dpi=150)
