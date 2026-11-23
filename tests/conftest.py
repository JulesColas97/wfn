import os
import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
OLD_CODE = ROOT / "old_project" / "code" / "prepost"
OLD_NREL = ROOT / "old_project" / "examples" / "3_wind_turbine_aero" / "NREL"
NREL_DB = ROOT / "examples" / "1_source_nrel" / "NREL5MW_14seg.nc"


@pytest.fixture(scope="session")
def old_prepost():
    """Import the old source model without running ``prepost/__init__.py``
    (which pulls GUI/plotting dependencies that are not needed here)."""
    if not OLD_CODE.is_dir():
        pytest.skip("old_project not available")
    pkg = types.ModuleType("prepost")
    pkg.__path__ = [str(OLD_CODE)]
    sys.modules["prepost"] = pkg
    for name in ["pyvista", "typer", "dash", "dash_bootstrap_components", "plotly",
                 "plotly.graph_objects", "requests"]:
        sys.modules.setdefault(name, types.ModuleType(name))
    ft = sys.modules.setdefault("func_timeout", types.ModuleType("func_timeout"))
    ft.func_timeout = None
    ft.FunctionTimedOut = Exception
    from prepost.source import wind_turbine_noise
    return wind_turbine_noise


@pytest.fixture
def in_dir():
    """Temporarily change directory (old code resolves paths from the cwd)."""
    cwd = os.getcwd()

    def _cd(path):
        os.chdir(path)

    yield _cd
    os.chdir(cwd)
