"""Page layout. Slider ids are the ones used by RP_map_all_app.py."""
from dash import dcc, get_asset_url, html

from .figures import graph_config, placeholder_figure
from .params import GROUPS, PRESETS

REPO_URL = "https://github.com/singhtaman/regular_precession_webapp"
GITHUB_URL = "https://github.com/singhtaman"
WEBSITE_URL = "https://singhtaman.github.io/"


def symbol(p):
    return html.Span([p.symbol, html.Sub(p.sub) if p.sub else None], className="sym")


def unit_label(unit):
    return ["M", html.Sub("⊙")] if unit == "M☉" else unit


def control(p):
    unit = html.Span(unit_label(p.unit), className="control-unit") if p.unit and not p.angle else None
    value = html.Span(id=f"{p.id}-value", className="control-value") if p.angle else None
    return html.Div(className="control", children=[
        html.Div(className="control-head", children=[
            html.Div(className="control-label", children=[symbol(p), html.Span(p.name, className="control-name"), unit]),
            value,
        ]),
        dcc.Slider(
            id=p.id, min=p.min, max=p.max, step=p.step, value=p.value, marks=p.marks,
            tooltip={"placement": "bottom", "always_visible": False},
            updatemode="mouseup",
            # Angles snap to multiples of π/16 and are shown as fractions of π instead.
            allow_direct_input=not p.angle,
            className="slider",
        ),
        html.P(p.hint, className="control-hint"),
    ])


def stats_card():
    return html.Article(className="plot-card stats-card", children=[
        html.Header(className="plot-head", children=[
            html.H3("System properties", className="plot-title"),
            html.P("Derived from the source and orientation parameters", className="plot-subtitle"),
        ]),
        html.Div(id="stats", className="stats"),
    ])


def group(key, title, desc, params):
    extra = None
    if key == "orientation":
        extra = html.Div(className="presets", children=[
            html.Span("Presets", className="presets-label"),
            *[html.Button([html.Span(v["label"]),
                           html.Small([v["detail"][0], html.Sub(v["detail"][1]), v["detail"][2]])],
                          id={"type": "preset", "name": name}, className="chip", n_clicks=0,
                          title=f"Set θ_S, φ_S, θ_J, φ_J to the {v['label'].lower()} system in systems_lib.py")
              for name, v in PRESETS.items()],
        ])
    return html.Details(open=True, className="group", children=[
        html.Summary(className="group-summary", children=[
            html.Span(title, className="group-title"),
            html.Span(desc, className="group-desc"),
        ]),
        html.Div(className="group-body", children=[
            extra,
            *[control(p) for p in params],
            html.Button("Reset " + title.lower(), id={"type": "reset-group", "group": key},
                        className="link-button", n_clicks=0),
        ]),
    ])


def plot_card(graph_id, title, subtitle, filename):
    return html.Article(className="plot-card", children=[
        html.Header(className="plot-head", children=[
            html.H3(title, className="plot-title"),
            html.P(subtitle, className="plot-subtitle"),
        ]),
        # Fixed-height frame: the plot's size never depends on Plotly measuring its container.
        html.Div(className="plot-frame", children=dcc.Graph(
            id=graph_id, figure=placeholder_figure(), config=graph_config(filename),
            style={"height": "100%", "width": "100%"})),
    ])


def legend():
    item = lambda cls, text: html.Span([html.I(className=f"key {cls}", **{"aria-hidden": "true"}), text], className="legend-item")
    return html.Div(className="legend", children=[
        item("key-rp", "Regular precession (RP)"),
        item("key-np", "Non-precessing (NP) reference"),
        item("key-j", "J direction"),
    ])


def about():
    return html.Details(className="about", children=[
        html.Summary("About this figure and model"),
        html.Div(className="about-body", children=[
            html.P([
                "Stellar binary black holes coalescing in the LIGO band can have misaligned spins. The orbital "
                "angular momentum ", html.B("L"), " then precesses and nutates about the total angular momentum ",
                html.B("J"), ", modulating the gravitational-wave signal. In ", html.B("regular precession"),
                " L moves on a cone around J with no nutation, and the direction of J is nearly constant over "
                "the radiation-reaction time.",
            ]),
            html.P([
                "Two dimensionless parameters describe the cone in the post-Newtonian waveform: the precession "
                "amplitude θ̃ sets its opening angle, θ", html.Sub("LJ"), " = 0.1 θ̃ (f/f", html.Sub("cut"),
                ")", html.Sup("1/3"), ", and the precession frequency Ω̃ sets how fast L moves around it. "
                "The precessional phase γ", html.Sub("P"), " is a nuisance parameter. The panels compare a "
                "regularly precessing signal (RP, blue) with the non-precessing signal (NP) of the same source.",
            ]),
            html.Ul([
                html.Li([html.B("Modulation of L: "), "the direction of L in the detector frame, cos θ", html.Sub("L"),
                         " against φ", html.Sub("L"), ". The red point marks the direction of J."]),
                html.Li([html.B("Amplitude: "), "fractional change of the GW amplitude, B", html.Sub("RP"),
                         "/B", html.Sub("NP"), " − 1."]),
                html.Li([html.B("Phase: "), "phase difference Φ", html.Sub("RP"), " − Φ", html.Sub("NP"),
                         " in radians: the precessional polarization phase plus twice the precession "
                         "correction δφ, with the sign used by the original app."]),
            ]),
            html.P([
                "Frequencies run from 20 Hz up to the cut-off f", html.Sub("cut"),
                " = 1/(6", html.Sup("3/2"), " π M", html.Sub("z"), ") in steps of 0.05 Hz. "
                "Distances use flat ΛCDM with Planck 2018 parameters (H", html.Sub("0"),
                " = 67.4, Ω", html.Sub("m"), " = 0.315). Plots open on the original app's axis windows "
                "and extend them when the curves would otherwise run outside.",
            ]),
            html.P(["Drag a slider or type a value; plots update when you release the handle. "
                    "Use the plot toolbar to zoom, pan or save a PNG, and double-click a plot to reset its view."],
                   className="muted"),
        ]),
    ])


def build():
    return html.Div(className="app-shell", children=[
        html.A("Skip to results", href="#results", className="skip-link"),
        html.Header(className="topbar", children=[
            html.Div(className="brand", children=[
                html.Img(src=get_asset_url("bbh_spin.jpg"), alt="", className="brand-logo"),
                html.Div([
                    html.H1("Regular Precession", className="brand-title"),
                    html.P("How a precessing orbit imprints on gravitational waves from binary black holes",
                           className="brand-tagline"),
                ]),
            ]),
            html.Nav(className="top-links", children=[
                html.A([html.I(className="icon icon-code"), html.Span("Source")], href=REPO_URL,
                       target="_blank", rel="noopener", className="top-link"),
                html.A([html.I(className="icon icon-user"), html.Span("GitHub")], href=GITHUB_URL,
                       target="_blank", rel="noopener", className="top-link", title="GitHub profile"),
                html.A([html.I(className="icon icon-globe"), html.Span("Website")], href=WEBSITE_URL,
                       target="_blank", rel="noopener", className="top-link", title="Personal website"),
            ]),
        ]),

        html.Main(className="workspace", children=[
            html.Aside(className="controls", **{"aria-label": "Parameters"}, children=[
                html.Div(className="controls-head", children=[
                    html.H2("Parameters", className="section-title"),
                    html.Button("Reset all", id="reset-all", className="button button-ghost", n_clicks=0),
                ]),
                *[group(*g) for g in GROUPS],
            ]),

            html.Section(id="results", className="results", **{"aria-label": "Results"}, children=[
                html.Div(className="results-bar", children=[
                    legend(),
                    html.Div(className="results-actions", children=[
                        html.Span(id="status", className="status", role="status"),
                        html.Button([html.I(className="icon icon-download"), "Export CSV"],
                                    id="export-csv", className="button", n_clicks=0),
                        dcc.Download(id="download"),
                    ]),
                ]),
                html.Div(id="notice", role="alert"),
                dcc.Loading(
                    type="dot", color="#2457d6", delay_show=250,
                    overlay_style={"visibility": "visible", "opacity": 0.55},
                    parent_className="plot-loading",
                    children=html.Div(id="plot-grid", className="plot-grid", children=[
                        plot_card("mod-graph", "Modulation of L",
                                  "Direction of the orbital angular momentum about J", "rp_modulation"),
                        stats_card(),
                        plot_card("amp-graph", "Amplitude",
                                  "Fractional change relative to the non-precessing signal", "rp_amplitude"),
                        plot_card("phase-graph", "Phase",
                                  "Phase difference relative to the non-precessing signal", "rp_phase"),
                    ]),
                ),
                about(),
            ]),
        ]),

        html.Footer(className="footer", children=[
            html.Span("© Tamanjyot Singh"),
            html.Span("Supported by NSF grant PHY-2309320", className="muted"),
        ]),
    ])
