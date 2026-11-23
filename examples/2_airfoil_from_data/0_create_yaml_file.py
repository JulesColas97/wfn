"""Build input.yaml for the MySE3.0-135 proxy turbine from the data in blade/."""
from pathlib import Path
import yaml

import numpy as np
import pandas as pd

from wfn.source import WindTurbine

HERE = Path(__file__).parent
BLADE = HERE / "blade"
TARGET_L_C = 2          # target span/chord ratio of the acoustic segments
OUT = HERE / "input.yaml"

# ---------------------------------------------------------------- turbine parameters
par = pd.read_csv(HERE / "turbine_parameters.csv").set_index("parameter")["value"]
hub_radius = float(par["hub radius"])
hub_height = float(par["hub height (rotor apex) above MSL"])
omega_rpm = float(par["rated rotor speed"])
pitch = float(par["fine pitch below rated"])

# ---------------------------------------------------------------- blade geometry
geo = pd.read_csv(BLADE / "blade_geometry.csv")
r = geo["r_from_apex_m"].to_numpy()
chord = geo["chord_m"].to_numpy()
pitch = float(par["fine pitch below rated"])

# ---------------------------------------------------------------- blade geometry
geo = pd.read_csv(BLADE / "blade_geometry.csv")
r = geo["r_from_apex_m"].to_numpy()
chord = geo["chord_m"].to_numpy()
twist = geo["twist_deg"].to_numpy()
span_norm = (geo["span_from_root_m"] / geo["span_from_root_m"].iloc[-1]).to_numpy()

# ---------------------------------------------------------------- airfoils
pos = pd.read_csv(BLADE / "airfoil_positions.csv")
pos.columns = ["airfoil", "s"]
pos["airfoil"] = pos["airfoil"].replace({"cylinder": "circular"})   # packaged circular.dat

# XFoil .dat files (name line + x y) next to the script
dat_dir = HERE / "airfoils_dat"
dat_dir.mkdir(exist_ok=True)
for f in (BLADE / "airfoils").glob("*.csv"):
    xy = pd.read_csv(f)[["x_over_c", "y_over_c"]].to_numpy()
    with open(dat_dir / f"{f.stem}.dat", "w") as fout:
        fout.write(f"{f.stem}\n")
        np.savetxt(fout, xy, fmt="%.6f")

names = list(dict.fromkeys(pos["airfoil"]))          # unique, root -> tip order
# each station takes the airfoil whose listed position is the closest
nearest = np.abs(span_norm[:, None] - pos["s"].to_numpy()[None, :]).argmin(axis=1)
airfoil_index = np.array([names.index(pos["airfoil"].iloc[i]) for i in nearest])

# ---------------------------------------------------------------- turbine + strips
wt = WindTurbine(
    seg=r, chord=chord, twist=twist,
    airfoil_index=airfoil_index, airfoil_name=names, airfoil_path=str(dat_dir.resolve()),
    hub_height=hub_height, n_blade=3, n_beta=12, omega_rpm=omega_rpm, pitch=pitch,
)
wt = wt.acoustic_strips(seg0=hub_radius, target_l_c=TARGET_L_C)

with open(OUT, "w") as f:
    yaml.safe_dump({"wind_turbine": wt.to_dict()}, f,
                   sort_keys=False, default_flow_style=None, width=10_000)
print(f"{wt.n_seg} segments, D = {wt.rotor_diameter:.2f} m -> {OUT}")
