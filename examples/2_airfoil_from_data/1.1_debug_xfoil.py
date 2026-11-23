from wfn.source.xfoil import run_xfoil

inputs = dict(chord=2.7, re=10 / 1.48e-5, mach=10 / 340, aoa=[0.0],
              trip_top=None, trip_bot=None, target_frac=0.95, rho=1.2, T=273.0)
try:
    run_xfoil(inputs, "airfoils_dat/FX77-W-500.dat", workdir="xfoil_debug")
except ValueError:
    pass
