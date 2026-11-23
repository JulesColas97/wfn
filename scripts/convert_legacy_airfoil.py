"""Convert an old ``AirfoilData.save`` pickle (.dat) into a netCDF airfoil database.

Usage:
    python scripts/convert_legacy_airfoil.py OLD.dat [NEW.nc]
"""
import sys
from pathlib import Path

from wfn.source import AirfoilDatabase


def main(argv):
    if len(argv) not in (2, 3):
        print(__doc__)
        return 1
    src = Path(argv[1])
    dst = Path(argv[2]) if len(argv) == 3 else src.with_suffix(".nc")
    db = AirfoilDatabase.from_legacy_pickle(src)
    db.save(dst)
    print(f"{src} -> {dst}")
    print(db.ds)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
