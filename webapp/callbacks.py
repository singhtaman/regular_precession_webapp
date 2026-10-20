"""Dash callbacks."""
import io
import logging
import time

import numpy as np
from dash import ALL, ClientsideFunction, Input, Output, State, ctx, html, no_update

from . import figures
from .layout import unit_label
from .params import BY_ID, CALL_ORDER, GROUPS, PARAMS, PRESETS
from .physics import DELTA_F, F_MIN, compute, derived_quantities

log = logging.getLogger(__name__)

SLIDER_INPUTS = [Input(pid, "value") for pid in CALL_ORDER]
SLIDER_STATES = [State(pid, "value") for pid in CALL_ORDER]


def _validate(values):
    """Return an error message for missing or out-of-range values, else None."""
    for pid, v in zip(CALL_ORDER, values):
        p = BY_ID[pid]
        if v is None:
            return f"Enter a value for {p.name.lower()}."
        if not (p.min - 1e-9 <= v <= p.max + 1e-9):
            return f"{p.name} must be between {p.min:g} and {p.max:g}{' ' + p.unit if p.unit else ''}."
    return None


def _notice(kind, title, body):
    icon = {"warn": "!", "error": "×", "info": "i"}[kind]
    if isinstance(body, list) and body and isinstance(body[0], list):  # several notes
        body = html.Ul([html.Li(item) for item in body])
    else:
        body = html.P(body)
    return [html.Span(icon, className="notice-icon", **{"aria-hidden": "true"}),
            html.Div([html.Strong(title), body])], f"notice notice-{kind}"


def _sci(v):
    exp = int(np.floor(np.log10(abs(v))))
    return [f"{v / 10 ** exp:.2g}×10", html.Sup(str(exp))] if exp >= 4 else [f"{v:.3g}"]


def _notes(r):
    """Plain-language notes for results that are correct but easy to misread."""
    notes = []
    bad = sum(int(np.count_nonzero(~np.isfinite(a))) for a in (r.phi_L, r.cos_L, r.amp, r.phase))
    if bad:
        notes.append([f"{bad:,} computed values are not finite (NaN or ∞) and are left out of the plots."])
    if figures.is_point(r.phi_L, r.cos_L):
        notes.append(["L does not move for these parameters, so the modulation of L is a single point."])
    if figures.is_flat_zero(r.amp, r.phase):
        notes.append(["RP coincides with the non-precessing signal here, so its curves lie on the zero line."])
    peak = figures.amp_peak(r)
    if peak:
        notes.append(["The non-precessing amplitude B", html.Sub("NP"), " is ≈ 0 for this orientation, so B",
                      html.Sub("RP"), "/B", html.Sub("NP"), " − 1 diverges (peak ", *_sci(peak),
                      f"). The amplitude axis stops at ±{figures.AMP_DISPLAY_LIMIT:g}; use autoscale in the "
                      "plot toolbar to see the full range."])
    return notes


def _stat(label, value, unit, note):
    return html.Div(className="stat", children=[
        html.Span(label, className="stat-label"),
        html.Span([html.Span(value, className="stat-value"),
                   html.Span(unit, className="stat-unit") if unit else None]),
        html.Span(note, className="stat-note"),
    ])


def _stats(d):
    sub = html.Sub
    return [
        _stat("Total mass", f"{d.total_mass:.2f}", unit_label("M☉"), "source frame"),
        _stat("Redshifted mass", f"{d.total_mass_redshifted:.2f}", unit_label("M☉"), ["M", sub("z"), " = M(1+z)"]),
        _stat("Cut-off frequency", f"{d.f_cut:.2f}", "Hz", ["f", sub("cut"), ", upper end of the band"]),
        _stat(["cos ι", sub("JN")], f"{d.cos_i_JN:.2f}", "", "angle between J and N"),
        _stat("Luminosity distance", f"{d.luminosity_distance:.3f}", "Gpc", "flat ΛCDM, Planck 2018"),
    ]


def register(app):
    @app.callback(
        Output("mod-graph", "figure"),
        Output("amp-graph", "figure"),
        Output("phase-graph", "figure"),
        Output("stats", "children"),
        Output("status", "children"),
        Output("notice", "children"),
        Output("notice", "className"),
        Output("export-csv", "disabled"),
        *SLIDER_INPUTS,
    )
    def update(*values):
        """Returns three figures, stats, status, notice (children, class) and export-disabled."""
        msg = _validate(values)
        if msg:
            return (no_update,) * 5 + _notice("warn", "Check the parameters", msg) + (True,)

        d = derived_quantities(*values[:7])
        t0 = time.perf_counter()
        try:
            r = compute(*values)
        except Exception as exc:  # report model failures instead of a silent callback error
            log.exception("Model evaluation failed for %s", dict(zip(CALL_ORDER, values)))
            body, title = f"{type(exc).__name__}: {exc}", "The model could not be evaluated for these parameters"
            blank = figures.placeholder_figure("No plot for these parameters: see the message above")
            return (blank, blank, blank, _stats(d), "", *_notice("error", title, body), True)
        ms = (time.perf_counter() - t0) * 1000

        if r.f.size == 0:
            blank = figures.placeholder_figure(
                f"No frequencies in band: f<sub>cut</sub> = {r.f_cut:.2f} Hz < f<sub>min</sub> = {F_MIN} Hz")
            return (blank, blank, blank, _stats(d), "",
                    *_notice("warn", "No signal in the detector band",
                             ["The cut-off frequency f", html.Sub("cut"), f" = {r.f_cut:.2f} Hz is below f",
                              html.Sub("min"), f" = {F_MIN} Hz, so there are no frequencies to plot. "
                              "Lower the chirp mass or redshift, or raise the mass ratio."]), True)

        status = f"{r.f.size:,} samples · {F_MIN}–{r.f[-1]:.1f} Hz, Δf = {DELTA_F} Hz · {ms:.0f} ms"
        notes = _notes(r)
        notice, notice_cls = (_notice("info", "About these results", notes) if notes
                              else (None, "notice is-empty"))
        return (figures.modulation_figure(r), figures.amplitude_figure(r), figures.phase_figure(r),
                _stats(d), status, notice, notice_cls, False)

    # Angle readouts as fractions of π, formatted in the browser (no server round trip).
    for p in PARAMS:
        if p.angle:
            app.clientside_callback(
                ClientsideFunction("rp", "angle"),
                Output(f"{p.id}-value", "children"),
                Input(p.id, "drag_value"),
                Input(p.id, "value"),
            )

    group_ids = {key: [p.id for p in params] for key, _, _, params in GROUPS}

    @app.callback(
        [Output(p.id, "value") for p in PARAMS],
        Input("reset-all", "n_clicks"),
        Input({"type": "reset-group", "group": ALL}, "n_clicks"),
        Input({"type": "preset", "name": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def reset(*_):
        trig = ctx.triggered_id
        if not trig or not ctx.triggered[0]["value"]:
            return [no_update] * len(PARAMS)
        if trig == "reset-all":
            new = {p.id: p.value for p in PARAMS}
        elif trig["type"] == "reset-group":
            new = {pid: BY_ID[pid].value for pid in group_ids[trig["group"]]}
        else:
            new = PRESETS[trig["name"]]["values"]
        return [new.get(p.id, no_update) for p in PARAMS]

    @app.callback(
        Output("download", "data"),
        Input("export-csv", "n_clicks"),
        *SLIDER_STATES,
        prevent_initial_call=True,
    )
    def export(n_clicks, *values):
        if not n_clicks or _validate(values):
            return no_update
        try:
            r = compute(*values)
        except Exception:  # the button is disabled in this state; never fail the request
            log.exception("Export failed for %s", dict(zip(CALL_ORDER, values)))
            return no_update
        if r.f.size == 0:
            return no_update
        buf = io.StringIO()
        buf.write("# Regular precession explorer export\n# Parameters:\n")
        for pid, v in zip(CALL_ORDER, values):
            p = BY_ID[pid]
            buf.write(f"#   {p.name} ({p.symbol}{('_' + p.sub) if p.sub else ''}) = {v!r}"
                      f"{' ' + p.unit if p.unit else ''}\n")
        buf.write(f"# f_cut = {r.f_cut!r} Hz; frequency grid {F_MIN} Hz to f_cut, step {DELTA_F} Hz\n")
        cols = ["f_Hz", "phi_L_RP", "cos_theta_L_RP", "amp_RP_over_NP_minus_1", "phase_RP_minus_NP_rad"]
        data = np.column_stack([r.f, r.phi_L, r.cos_L, r.amp, r.phase])
        np.savetxt(buf, data, delimiter=",", header=",".join(cols), comments="", fmt="%.17g")
        name = "regular_precession_Mc{:g}_q{:g}_z{:g}.csv".format(*values[:3])
        return {"content": buf.getvalue(), "filename": name, "type": "text/csv"}
