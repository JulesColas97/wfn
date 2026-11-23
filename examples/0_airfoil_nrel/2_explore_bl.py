"""Interactive explorer of a boundary layer database (Dash app).

    pip install -e ".[explore]"
    python 2_explore_bl.py [database.nc]     # then open http://127.0.0.1:8050
"""
import sys
from pathlib import Path

from wfn.source import AirfoilDatabase
from wfn.source.explore import explore

here = Path(__file__).parent
candidates = [Path(sys.argv[1])] if len(sys.argv) > 1 else [
    here / "NREL5MW_14seg_xfoil.nc", here / "../source_nrel/NREL5MW_14seg.nc"]
path = next(p for p in candidates if p.is_file())
print(f"exploring {path}")
explore(AirfoilDatabase.open(path))
