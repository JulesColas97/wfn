import matplotlib

matplotlib.use("Agg")

import pytest

from wfn.source import AirfoilDatabase
from wfn.source.explore import plot_database

from conftest import NREL_DB


@pytest.fixture(scope="module")
def db():
    return AirfoilDatabase.open(NREL_DB)


def test_plot_database(db):
    fig = plot_database(db, seg=13)
    assert len([ax for ax in fig.axes if ax.get_visible()]) == 12
    assert "NACA64" in fig._suptitle.get_text()


def test_dash_figure_and_app(db):
    pytest.importorskip("dash")
    from wfn.source.explore import bl_figure, build_app

    fig = bl_figure(db, "cf_top", iu=2, seg=5)
    assert len(fig.data) == 2
    assert fig.data[1].x.size == db.ds.sizes["aoa"]
    app = build_app(db)
    assert app.layout is not None
