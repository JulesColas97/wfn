"""Plot the NREL 5MW boundary layer database.

    python 1_plot_bl.py [database.nc] [segment]

Default database: NREL5MW_14seg_xfoil.nc (from 0_create_bl_data.py), else the
converted legacy database of examples/source_nrel.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from wfn.source import AirfoilDatabase
from wfn.source.explore import plot_database, segment_label

here = Path(__file__).parent
candidates = [Path(sys.argv[1])] if len(sys.argv) > 1 else [
    here / "NREL5MW_14seg_xfoil.nc", here / "../source_nrel/NREL5MW_14seg.nc"]
path = next(p for p in candidates if p.is_file())
seg = int(sys.argv[2]) if len(sys.argv) > 2 else -1

db = AirfoilDatabase.open(path)
seg = seg % db.n_seg
print(f"{path.name}: {segment_label(db, seg)}")

# 1. all BL variables vs AoA (markers: database, lines: interpolation used by the source model)
fig = plot_database(db, seg)
fig.savefig(f"bl_seg{seg}.png", dpi=150)

# 2. cl / cd against the NREL 5MW reference polars (fastest velocity)
reference = {"DU40": "DU40_A17", "DU35": "DU35_A17", "DU30": "DU30_A17", "DU25": "DU25_A17",
             "DU21": "DU21_A17", "NACA64": "NA64_A17"}
names = db.ds.attrs.get("airfoil_name")
name = list(names)[int(db.ds.airfoil_index.values[seg])] if names is not None else None
iu = db.ds.sizes["uinf"] - 1
fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
for ax, var in zip(axes, ("cl", "cd")):
    # NB: as in the old code, the stored XFoil coefficients are divided by the chord
    ax.plot(db.ds.aoa, db.ds[var].isel(uinf=iu, seg=seg), "o-", ms=3, label="xfoil database")
    if name in reference:
        ref = np.loadtxt(here / "comparison_data" / reference[name], skiprows=31)[1:]
        ref = ref[np.abs(ref[:, 0]) <= 30]
        ax.plot(ref[:, 0], ref[:, 1 if var == "cl" else 2], label=f"NREL {reference[name]}")
    ax.set_xlabel("AoA (deg)")
    ax.set_ylabel(var)
    ax.set_xlim(-30, 30)
    ax.grid(True, alpha=0.3)
axes[1].legend()
fig.suptitle(f"{segment_label(db, seg)}, uinf={db.ds.uinf.values[iu]:.1f} m/s")
fig.tight_layout()
fig.savefig(f"cl_cd_seg{seg}.png", dpi=150)
print(f"saved bl_seg{seg}.png and cl_cd_seg{seg}.png")
plt.show()
