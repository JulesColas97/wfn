"""WAPE runs around one source (no flow, flat impedance ground), then deltaL on the receiver plane.

Needs the compiled WAPE code: set ``propagation.executable`` in input.yaml or the
environment variable ``WFN_WAPE``.
"""
import matplotlib.pyplot as plt

from wfn.propagation import PropagationConfig, assemble_polar, loop_angle, prepare_cases, run_local, to_grid
from wfn.receivers import ReceiversConfig

cfg = PropagationConfig.from_yaml("input_side.yaml")
prepare_cases(cfg)            # input.nml + 2D_WAPE in results/wape/t0/<h>/ and run.sh
run_local(cfg.work_path)      # runs every case one after the other

ds = assemble_polar(cfg, isource=0)   # deltaL(r, tau, z, freq, source_height)
ds.to_netcdf(cfg.output_file(0))

plane = to_grid(loop_angle(ds), ReceiversConfig.from_yaml("input_side.yaml").build())
fig, axes = plt.subplots( plane.sizes["source_height"],1, figsize=(10, 4), sharey=True)
for ax, h in zip(axes, plane.source_height.values):
    plane.deltaL.sel(source_height=h).isel(y=0, freq=0).plot(
        ax=ax, x="x", y="z", cmap="RdBu_r", vmin=-10, vmax=10, cbar_kwargs={"label": "deltaL (dB)"}
    )
    ax.set_aspect("equal")
    ax.set_title(f"source at {h:.0f} m, f = {plane.freq.values[0]:.0f} Hz")
plt.tight_layout()
plt.savefig("deltaL_plane.png", dpi=150)
print("saved deltaL_plane.png and", cfg.output_file(0))
plt.show()
