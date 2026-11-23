"""Visual inspection of an :class:`~wfn.source.airfoil.AirfoilDatabase`.

* :func:`plot_database`: static matplotlib figure (needs ``matplotlib``);
* :func:`explore`: interactive Dash app, port of the old
  ``prepost.source.explore_airfoil_db.explore`` (needs ``pip install -e .[explore]``).

Both show the raw database points (markers) and the interpolation actually
used by the source model (lines, :meth:`AirfoilDatabase.evaluate`).
"""
from __future__ import annotations

from typing import Sequence

import numpy as np

from .airfoil import ALL_VARS, AirfoilDatabase

DEFAULT_VARS = ("delta_star_top", "delta_star_bot", "theta_momen_top", "theta_momen_bot",
                "cf_top", "cf_bot", "dpdx_top", "dpdx_bot", "UeUinf_top", "UeUinf_bot",
                "cl", "cd")


def _fine_aoa(db: AirfoilDatabase, n: int = 200) -> np.ndarray:
    return np.linspace(db.ds.aoa.min().item(), db.ds.aoa.max().item(), n)


def segment_label(db: AirfoilDatabase, seg: int) -> str:
    """``"seg 13: NACA64, c=1.42 m"``."""
    label = f"seg {seg}"
    names = db.ds.attrs.get("airfoil_name")
    if names is not None and "airfoil_index" in db.ds:
        names = [names] if isinstance(names, str) else list(names)
        label += f": {names[int(db.ds.airfoil_index.values[seg])]}"
    return label + f", c={db.chord[seg]:.2f} m"


def chord_reynolds(db: AirfoilDatabase, iu: int, seg: int) -> float | None:
    nu = db.ds.attrs.get("nu")
    return None if nu is None else float(db.ds.uinf.values[iu] * db.chord[seg] / nu)


# -----------------------------------------------------------------------------
# matplotlib
# -----------------------------------------------------------------------------
def plot_database(db: AirfoilDatabase, seg: int, variables: Sequence[str] | None = None,
                  ncols: int = 4):
    """Grid of ``variables`` vs AoA for segment ``seg``, one colour per velocity.

    Returns:
        the matplotlib ``Figure``.
    """
    import matplotlib.pyplot as plt

    variables = [v for v in (variables or DEFAULT_VARS) if v in db.ds]
    nrows = int(np.ceil(len(variables) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.6 * nrows), squeeze=False,
                             sharex=True)
    aoa = db.ds.aoa.values
    aoa_fine = _fine_aoa(db)
    colors = plt.cm.viridis(np.linspace(0, 0.9, db.ds.sizes["uinf"]))
    for ax, var in zip(axes.flat, variables):
        for iu, u in enumerate(db.ds.uinf.values):
            ax.plot(aoa, db.ds[var].isel(uinf=iu, seg=seg).values, "o", ms=2.5, color=colors[iu])
            ax.plot(aoa_fine, db.evaluate(var, seg, np.full_like(aoa_fine, u), aoa_fine), "-",
                    lw=1, color=colors[iu], label=f"{u:.1f} m/s")
        ax.set_title(var, fontsize=10)
        ax.grid(True, alpha=0.3)
    for ax in axes.flat[len(variables):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("AoA (deg)")
    axes.flat[0].legend(fontsize=7, title="uinf", title_fontsize=7)
    fig.suptitle(segment_label(db, seg))
    fig.tight_layout()
    return fig


# -----------------------------------------------------------------------------
# Dash
# -----------------------------------------------------------------------------
def bl_figure(db: AirfoilDatabase, var: str, iu: int, seg: int):
    """Plotly figure of one variable for one velocity and one segment."""
    import plotly.graph_objects as go

    u = float(db.ds.uinf.values[iu])
    aoa_fine = _fine_aoa(db)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=aoa_fine,
                             y=db.evaluate(var, seg, np.full_like(aoa_fine, u), aoa_fine),
                             mode="lines", name="interpolation"))
    fig.add_trace(go.Scatter(x=db.ds.aoa.values, y=db.ds[var].isel(uinf=iu, seg=seg).values,
                             mode="markers", name="database"))
    re = chord_reynolds(db, iu, seg)
    re_txt = f", Re_c={re:.3g}" if re is not None else ""
    fig.update_layout(title=f"{segment_label(db, seg)} | uinf={u:.2f} m/s{re_txt}",
                      xaxis_title="AoA (deg)", yaxis_title=var)
    return fig


def build_app(db: AirfoilDatabase):
    """Create the Dash app (without running it)."""
    try:
        import dash
        import dash_bootstrap_components as dbc
        from dash import Input, Output, dcc
    except ImportError as e:
        raise ImportError("the explorer needs Dash: pip install -e .[explore]") from e

    variables = [v for v in ALL_VARS if v in db.ds]
    n_u, n_seg = db.ds.sizes["uinf"], db.n_seg
    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
    app.layout = dbc.Container(dbc.Row([
        dbc.Col([
            dbc.Label("Variable"),
            dcc.Dropdown(id="variable", options=variables, value=variables[0],
                         style={"margin": "20px 0"}),
            dbc.Label("Velocity uinf (m/s)"),
            dcc.Slider(id="uinf", min=0, max=n_u - 1, step=1, value=0,
                       marks={i: f"{u:.1f}" for i, u in enumerate(db.ds.uinf.values)}),
            dbc.Label("Segment"),
            dcc.Slider(id="seg", min=0, max=n_seg - 1, step=1, value=0,
                       marks={i: str(i) for i in range(n_seg)}),
        ], width=4),
        dbc.Col([dcc.Graph(id="plot", style={"height": "600px"})], width=8),
    ]), fluid=True)

    @app.callback(Output("plot", "figure"),
                  [Input("variable", "value"), Input("uinf", "value"), Input("seg", "value")])
    def update(var, iu, seg):
        return bl_figure(db, var, iu, seg)

    return app


def explore(db: AirfoilDatabase, port: int = 8050, debug: bool = False) -> None:
    """Run the interactive explorer at ``http://127.0.0.1:<port>``."""
    build_app(db).run(port=port, debug=debug)
