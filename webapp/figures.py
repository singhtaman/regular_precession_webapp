"""Plotly figures for the three result panels.

Presentation only: values from physics.Result are plotted unchanged. Default
axis windows match the original app; a window is widened only when data would
otherwise fall outside it.
"""
import math
from fractions import Fraction

import numpy as np
import plotly.graph_objects as go

# RP blue (as in the original app), J red, non-precessing reference (zero line) black dashed.
RP_COLOR = "#2457d6"
J_COLOR = "#d93a2b"
REF_COLOR = "#1d2433"
TEXT = "#1d2433"
MUTED = "#5d6677"
GRID = "#eceff4"
AXIS = "#c9cfd9"
FONT = "Inter, system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif"

FREQ_TICKS = [20, 30, 40, 60, 80, 100, 150, 200, 300, 400]

# B_RP/B_NP − 1 diverges where the non-precessing amplitude B_NP vanishes (e.g. an edge-on
# source whose + polarization the detector cannot see). Beyond this value the axis stops
# growing and the panel is labelled instead, so a 10^15 spike doesn't flatten the curve.
AMP_DISPLAY_LIMIT = 20.0
POINT_TOL = 1e-9


def is_point(x, y):
    """True when a curve has (numerically) zero extent, i.e. it is a single point."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    return x.size > 0 and np.ptp(x) < POINT_TOL and np.ptp(y) < POINT_TOL


def is_flat_zero(*arrays):
    """True when every array is (numerically) zero, i.e. the curves sit on the zero line."""
    return all(np.size(a) and np.max(np.abs(a)) < POINT_TOL for a in arrays)


def amp_peak(r):
    """Largest |B_RP/B_NP − 1| if it exceeds the display limit, else None."""
    peak = float(np.nanmax(np.abs(r.amp))) if r.f.size else 0.0
    return peak if peak > AMP_DISPLAY_LIMIT else None


def sci(v):
    """Format a number as e.g. 3.9×10<sup>15</sup> for Plotly/HTML text."""
    if v == 0 or 1e-3 <= abs(v) < 1e4:
        return f"{v:.3g}"
    exp = int(math.floor(math.log10(abs(v))))
    return f"{v / 10 ** exp:.2g}×10<sup>{exp}</sup>"


def _base_layout(xaxis_title, yaxis_title):
    axis = dict(showline=True, linecolor=AXIS, linewidth=1, mirror=True,
                ticks="outside", tickcolor=AXIS, ticklen=4,
                gridcolor=GRID, zeroline=False,
                title=dict(font=dict(size=13, color=TEXT), standoff=8),
                tickfont=dict(size=11.5, color=MUTED))
    return dict(
        font=dict(family=FONT, size=12, color=TEXT),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#ffffff",
        margin=dict(l=62, r=18, t=14, b=50),
        showlegend=False,
        autosize=True,
        hoverlabel=dict(bgcolor="#ffffff", bordercolor=AXIS,
                        font=dict(family=FONT, size=12, color=TEXT)),
        xaxis=dict(axis, title=dict(axis["title"], text=xaxis_title)),
        yaxis=dict(axis, title=dict(axis["title"], text=yaxis_title)),
        uirevision="keep",
    )


def _fit(default_lo, default_hi, *arrays, pad_frac=0.04):
    """Default window, widened to include all finite data plus a small margin."""
    finite = [a[np.isfinite(a)] for a in arrays if a is not None and np.size(a)]
    finite = [a for a in finite if a.size]
    if not finite:
        return [default_lo, default_hi]
    lo = min(float(a.min()) for a in finite)
    hi = max(float(a.max()) for a in finite)
    pad = pad_frac * max(hi - lo, 1e-9)
    return [min(default_lo, lo - pad), max(default_hi, hi + pad)]


def _pi_label(k: Fraction) -> str:
    if k == 0:
        return "0"
    sign = "−" if k < 0 else ""
    k = abs(k)
    num = "" if k.numerator == 1 else str(k.numerator)
    return f"{sign}{num}π" if k.denominator == 1 else f"{sign}{num}π/{k.denominator}"


def _pi_ticks(lo, hi, max_ticks=7):
    for denom in (8, 4, 2, 1):
        start = math.ceil(lo / (math.pi / denom))
        stop = math.floor(hi / (math.pi / denom))
        if stop - start + 1 <= max_ticks:
            break
    ks = [Fraction(i, denom) for i in range(start, stop + 1)]
    return [float(k) * math.pi for k in ks], [_pi_label(k) for k in ks]


def _break_wraps(phi, y):
    """Insert gaps where φ_L jumps across ±π so lines don't streak across the panel."""
    if phi.size < 2:
        return phi, y
    jumps = np.where(np.abs(np.diff(phi)) > np.pi)[0] + 1
    if not jumps.size:
        return phi, y
    return np.insert(phi.astype(float), jumps, np.nan), np.insert(y.astype(float), jumps, np.nan)


def modulation_figure(r):
    phi_J, cos_J = r.J
    # φ_L comes from arctan2 and lies in (−π, π]; show J at the same angle in that convention.
    phi_J_plot = float(np.arctan2(np.sin(phi_J), np.cos(phi_J)))

    fig = go.Figure(layout=_base_layout("φ<sub>L</sub>", "cos θ<sub>L</sub>"))
    fig.add_trace(go.Scatter(x=[phi_J_plot], y=[cos_J], mode="markers", name="J",
                             marker=dict(color=J_COLOR, size=14, line=dict(color="#ffffff", width=1.5)),
                             hovertemplate="J<br>φ<sub>J</sub> = %{x:.4f}<br>cos θ<sub>J</sub> = %{y:.4f}<extra></extra>"))
    # A stationary L is a single point; a line through identical points would draw nothing.
    point = is_point(r.phi_L, r.cos_L)
    x, y = (r.phi_L[:1], r.cos_L[:1]) if point else _break_wraps(r.phi_L, r.cos_L)
    fig.add_trace(go.Scatter(x=x, y=y, mode="markers" if point else "lines", name="RP",
                             line=dict(color=RP_COLOR, width=2),
                             marker=dict(color=RP_COLOR, size=8),
                             hovertemplate="RP<br>φ<sub>L</sub> = %{x:.4f}<br>"
                                           "cos θ<sub>L</sub> = %{y:.4f}<extra></extra>"))

    xr = _fit(phi_J_plot - np.pi / 4, phi_J_plot + np.pi / 4, r.phi_L)
    yr = _fit(cos_J - 0.8, cos_J + 0.8, r.cos_L)
    yr = [max(yr[0], -1.05), min(yr[1], 1.05)]
    tickvals, ticktext = _pi_ticks(*xr)
    fig.update_xaxes(range=xr, tickvals=tickvals, ticktext=ticktext)
    fig.update_yaxes(range=yr)
    fig.update_layout(hovermode="closest")
    return fig


def _frequency_figure(r, y, y_title, default_y, digits, limit=None):
    fig = go.Figure(layout=_base_layout("f (Hz)", y_title))
    # layer="below" keeps an RP curve that coincides with the reference visible on top of it.
    fig.add_hline(y=0, line=dict(color=REF_COLOR, width=1, dash="dash"), opacity=0.7, layer="below")
    fig.add_trace(go.Scatter(x=r.f, y=y, mode="lines", name="RP",
                             line=dict(color=RP_COLOR, width=2),
                             hovertemplate=f"RP: %{{y:.{digits}f}}<extra></extra>"))

    f_hi = max(150.0, float(r.f[-1]) * 1.03) if r.f.size else 150.0
    ticks = [t for t in FREQ_TICKS if t <= f_hi]
    fig.update_xaxes(type="log", range=[math.log10(19), math.log10(f_hi)],
                     tickvals=ticks, ticktext=[str(t) for t in ticks],
                     hoverformat=".2f")
    yr = _fit(*default_y, y)
    if limit is not None and (yr[0] < -limit or yr[1] > limit):
        yr = [max(yr[0], -limit), min(yr[1], limit)]
        peak = float(np.nanmax(np.abs(y)))
        fig.add_annotation(
            text=f"Off scale: values reach {sci(peak)}<br>(B<sub>NP</sub> ≈ 0 here, so the ratio diverges)",
            xref="paper", yref="paper", x=0.5, y=0.97, xanchor="center", yanchor="top", showarrow=False,
            align="center", font=dict(size=11.5, color="#6d4500"),
            bgcolor="#fff8eb", bordercolor="#f1d49a", borderwidth=1, borderpad=5)
    fig.update_yaxes(range=yr)
    fig.update_layout(hovermode="x unified")
    return fig


def amplitude_figure(r):
    return _frequency_figure(r, r.amp, "B<sub>RP</sub> / B<sub>NP</sub> − 1", (-1.05, 0.95), 4,
                             limit=AMP_DISPLAY_LIMIT)


def phase_figure(r):
    return _frequency_figure(r, r.phase, "Φ<sub>RP</sub> − Φ<sub>NP</sub> (rad)", (-3.25, 0.75), 4)


def placeholder_figure(message=None):
    """Axis-less figure used before the first computation and in empty/error states.

    The graph keeps its full size in these states so it can never collapse or be hidden.
    """
    fig = go.Figure(layout=dict(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                                margin=dict(l=0, r=0, t=0, b=0), autosize=True,
                                xaxis=dict(visible=False), yaxis=dict(visible=False)))
    if message:
        fig.add_annotation(text=message, xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False,
                           font=dict(family=FONT, size=13, color=MUTED))
    return fig


def graph_config(filename):
    return {
        "displaylogo": False,
        "responsive": True,
        "modeBarButtonsToRemove": ["select2d", "lasso2d"],
        "toImageButtonOptions": {"format": "png", "filename": filename, "scale": 3},
    }
