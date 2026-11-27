"""SPL map of the NREL 5MW on the receiver grid of input.yaml, for every turbine position.

Each turbine result is saved as ``<source.output>_<i>.nc`` (an xarray Dataset).
"""
import matplotlib.pyplot as plt

from wfn import units
from wfn.receivers import ReceiversConfig
from wfn.source import WindTurbineNoise

src = WindTurbineNoise.from_yaml("input.yaml")
results = src.compute_all(ReceiversConfig.from_yaml("input.yaml"))

ds = results[0]
# rotor-averaged A-weighted overall level
# third-octave bands first: units.spl is a level per Hz and cannot be summed over freq
spl = units.db_mean(units.db_sum(units.spl_third_octave(ds), ["blade", "seg"]), "beta")
oaspl = units.db_sum(spl + units.a_weighting(spl.freq), "freq")

oaspl.isel(z=0).plot(x="x", y="y", cmap="viridis", cbar_kwargs={"label": "OASPL (dBA)"})
plt.gca().set_aspect("equal")
plt.title("NREL 5MW, free field, z = %.0f m" % ds.z.values[0])
plt.savefig("oaspl_plane.png", dpi=150)
print("saved oaspl_plane.png and", src.output)
