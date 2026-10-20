"""Browser tests (optional): `pip install playwright && playwright install chromium`, then `pytest -m ui`."""
import threading

import pytest

sync_api = pytest.importorskip("playwright.sync_api", reason="Playwright not installed")
pytestmark = pytest.mark.ui

GRAPHS = ["mod-graph", "amp-graph", "phase-graph"]

# Each plot's SVG must fill its graph frame and stay inside its card (never collapsed or clipped).
VISIBILITY = """() => %s.map(id => {
  const g = document.getElementById(id), card = g.closest('.plot-card'), svg = g.querySelector('.main-svg');
  const gr = g.getBoundingClientRect(), cr = card.getBoundingClientRect();
  const sr = svg ? svg.getBoundingClientRect() : {height: 0, bottom: 0};
  const ok = gr.height > 250 && Math.abs(sr.height - gr.height) < 2 && sr.bottom <= cr.bottom + 1;
  return ok ? null : `${id}: graph ${Math.round(gr.height)}px, svg ${Math.round(sr.height)}px`;
}).filter(Boolean)""" % GRAPHS


@pytest.fixture(scope="module")
def server_url(wsgi):
    from werkzeug.serving import make_server
    srv = make_server("127.0.0.1", 0, wsgi, threaded=True)  # threaded, like production servers
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_port}/"
    srv.shutdown()


@pytest.fixture(scope="module")
def browser():
    with sync_api.sync_playwright() as p:
        try:
            b = p.chromium.launch()
        except Exception as exc:  # browser binary not installed
            pytest.skip(f"Chromium not available: {exc}")
        yield b
        b.close()


@pytest.fixture()
def page(browser, server_url):
    pg = browser.new_page(viewport={"width": 1600, "height": 1000})
    pg.errors = []
    pg.on("pageerror", lambda e: pg.errors.append(str(e)))
    pg.on("console", lambda m: m.type == "error" and pg.errors.append(m.text))
    pg.goto(server_url)
    pg.wait_for_selector("#mod-graph .scatterlayer .trace")
    yield pg
    pg.close()


def _drag(pg, slider_id, frac):
    track = pg.locator(f"#{slider_id} .dash-slider-track").first
    track.scroll_into_view_if_needed()
    bb, th = track.bounding_box(), pg.locator(f"#{slider_id} .dash-slider-thumb").first.bounding_box()
    pg.mouse.move(th["x"] + th["width"] / 2, th["y"] + th["height"] / 2)
    pg.mouse.down()
    pg.mouse.move(bb["x"] + bb["width"] * frac, bb["y"] + bb["height"] / 2, steps=6)
    pg.mouse.up()
    pg.wait_for_timeout(900)


def test_renders_without_errors(page):
    for gid in GRAPHS:
        assert page.locator(f"#{gid} .scatterlayer .trace").count() > 0
    assert page.evaluate(VISIBILITY) == []
    assert page.errors == []


def test_two_by_two_arrangement(page):
    box = {k: page.locator(sel).bounding_box() for k, sel in
           {"mod": "#mod-graph", "stats": "#stats", "amp": "#amp-graph", "phase": "#phase-graph"}.items()}
    assert box["mod"]["x"] < box["stats"]["x"] and abs(box["mod"]["y"] - box["stats"]["y"]) < 120
    assert box["amp"]["x"] < box["phase"]["x"] and abs(box["amp"]["y"] - box["phase"]["y"]) < 1
    assert box["amp"]["y"] > box["mod"]["y"] + box["mod"]["height"] - 1


def test_plots_stay_visible_through_empty_states_and_back(page):
    """Regression: leaving an empty/error state used to leave plots collapsed behind their cards."""
    for slider, frac in [("theta_S-slider", 0), ("theta_J-slider", 0), ("theta_J-slider", 0.5),
                         ("chirp_mass-slider", 1), ("redshift-slider", 1), ("mass_ratio-slider", 1),
                         ("chirp_mass-slider", 0.1)]:
        _drag(page, slider, frac)
        assert page.evaluate(VISIBILITY) == [], f"after {slider} -> {frac}"
    assert page.errors == []


def test_enter_commits_typed_value(page):
    box = page.locator(".control:has(#chirp_mass-slider) input[type=number]").first
    box.click()
    box.fill("40")
    box.press("Enter")
    page.wait_for_timeout(1200)
    assert "95.52" in page.inner_text("#stats")  # total mass for M_c = 40, q = 0.6


def test_face_on_preset(page):
    page.click("text=Face-on")
    page.wait_for_timeout(1200)
    assert "1.00" in page.inner_text("#stats")  # cos ι_JN
    assert page.inner_text("#theta_J-slider-value").startswith("π/4")
