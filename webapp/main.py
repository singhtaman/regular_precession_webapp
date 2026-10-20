"""Dash application for the redesigned Regular Precession explorer."""
import warnings

warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")

from dash import Dash  # noqa: E402
from flask import request  # noqa: E402

from . import callbacks, layout  # noqa: E402

FONTS = ("https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700"
         "&family=JetBrains+Mono:wght@400;500&display=swap")

INDEX = """<!DOCTYPE html>
<html lang="en">
  <head>
    {%metas%}
    <title>{%title%}</title>
    <link rel="icon" type="image/jpeg" href="/assets/bbh_spin.jpg">
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    {%css%}
  </head>
  <body>
    {%app_entry%}
    <footer>{%config%}{%scripts%}{%renderer%}</footer>
  </body>
</html>"""

app = Dash(
    __name__,
    title="Regular Precession",
    update_title=None,
    external_stylesheets=[FONTS],
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1"},
        {"name": "description", "content": "Interactive explorer for regularly precessing "
                                           "binary black hole gravitational waveforms."},
        {"name": "theme-color", "content": "#ffffff"},
    ],
)
app.index_string = INDEX
app.layout = layout.build
callbacks.register(app)

server = app.server


@server.after_request
def _cdn_cache(response):
    # Dash component bundles are fingerprinted and assets are versioned per deployment, so
    # they can be cached at the CDN edge instead of invoking the function on every page load.
    if response.status_code == 200 and request.path.startswith(("/_dash-component-suites/", "/assets/")):
        response.headers["CDN-Cache-Control"] = "public, max-age=86400"
    return response
