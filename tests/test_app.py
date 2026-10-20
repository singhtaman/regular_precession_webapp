"""The Dash app over HTTP: page, assets, every callback and the exported data."""
import io
import json
import threading

import numpy as np
import pytest

from conftest import callback_outputs, notice_text
from webapp import physics
from webapp.params import BY_ID, CALL_ORDER, GROUPS, PARAMS, PRESETS

PI = np.pi
DEFAULTS = [BY_ID[k].value for k in CALL_ORDER]
VERCEL_MAX_PAYLOAD = int(4.5 * 1024 * 1024)


@pytest.mark.parametrize("path", ["/", "/_dash-layout", "/_dash-dependencies",
                                  "/assets/style.css", "/assets/app.js", "/assets/bbh_spin.jpg"])
def test_endpoints_respond(client, path):
    assert client.get(path).status_code == 200


def test_static_assets_cacheable_at_cdn(client):
    assert client.get("/assets/style.css").headers.get("CDN-Cache-Control") == "public, max-age=86400"
    assert "CDN-Cache-Control" not in client.get("/_dash-layout").headers


def test_index_has_title_and_viewport(client):
    html = client.get("/").get_data(as_text=True)
    assert "<title>Regular Precession</title>" in html
    assert 'name="viewport"' in html


def test_layout_contains_every_control_and_panel(client):
    layout = client.get("/_dash-layout").get_data(as_text=True)
    for pid in [p.id for p in PARAMS] + ["mod-graph", "amp-graph", "phase-graph", "stats", "status",
                                          "notice", "export-csv", "download", "reset-all"]:
        assert f'"{pid}"' in layout, pid


def test_default_callback_returns_three_figures(call_main):
    resp, size = call_main(DEFAULTS)
    for gid in ("mod-graph", "amp-graph", "phase-graph"):
        assert resp[gid]["figure"]["data"], gid
    assert "2,433 samples" in resp["status"]["children"]
    assert resp["notice"]["className"] == "notice is-empty"
    assert size < VERCEL_MAX_PAYLOAD


def test_largest_payload_fits_vercel_limit(call_main):
    resp, size = call_main([5, 1, 0.05] + DEFAULTS[3:])  # widest band → most samples
    assert "6,891 samples" in resp["status"]["children"]
    assert size < VERCEL_MAX_PAYLOAD


@pytest.mark.parametrize("values,kind,phrase", [
    ([50, 1, 1] + DEFAULTS[3:], "notice-warn", "No signal in the detector band"),
    ([None] + DEFAULTS[1:], "notice-warn", "Enter a value for chirp mass"),
    ([500] + DEFAULTS[1:], "notice-warn", "must be between 5 and 50"),
    ([10, 0.6, 0.3, 0, 0, PI / 2, PI / 4] + DEFAULTS[7:], "notice-info", "diverges"),
    (DEFAULTS[:7] + [0, 2, 0], "notice-info", "single point"),
])
def test_states_are_explained(call_main, values, kind, phrase):
    resp, _ = call_main(values)
    assert kind in resp["notice"]["className"]
    assert phrase in notice_text(resp)


@pytest.fixture()
def model_fails(monkeypatch):
    """Simulate a model failure (no real input currently triggers one)."""
    def boom(*args):
        raise FloatingPointError("invalid value encountered in scalar divide")
    monkeypatch.setattr("webapp.callbacks.compute", boom)


def test_model_failure_is_explained(call_main, model_fails):
    resp, _ = call_main(DEFAULTS)
    assert "notice-error" in resp["notice"]["className"]
    assert "could not be evaluated" in notice_text(resp)
    assert "FloatingPointError" in notice_text(resp)
    assert resp["export-csv"]["disabled"] is True


@pytest.mark.parametrize("state", ["out_of_band", "model_failure"])
def test_plots_are_never_hidden(call_main, monkeypatch, state):
    """Empty and error states keep every graph (placeholder figure + message), never remove it."""
    if state == "model_failure":
        monkeypatch.setattr("webapp.callbacks.compute", lambda *a: (_ for _ in ()).throw(FloatingPointError("x")))
    resp, _ = call_main([50, 1, 1] + DEFAULTS[3:] if state == "out_of_band" else DEFAULTS)
    assert "plot-grid" not in resp
    for gid in ("mod-graph", "amp-graph", "phase-graph"):
        ann = resp[gid]["figure"]["layout"].get("annotations", [])
        assert ann and "No " in ann[0]["text"], gid


@pytest.mark.parametrize("values,disabled", [
    (DEFAULTS, False),
    ([50, 1, 1] + DEFAULTS[3:], True),
    ([10, 0.6, 0.3, 0, 0, 0, PI / 2] + DEFAULTS[7:], False),  # J parallel to N: plots normally
])
def test_export_button_disabled_without_data(call_main, values, disabled):
    resp, _ = call_main(values)
    assert resp["export-csv"]["disabled"] is disabled


def test_export_does_not_fail_when_model_fails(client, dash_deps, model_fails):
    dep = next(d for d in dash_deps if "download.data" in d["output"])
    body = {"output": dep["output"], "outputs": {"id": "download", "property": "data"},
            "inputs": [{"id": "export-csv", "property": "n_clicks", "value": 1}],
            "state": [{"id": k, "property": "value", "value": v} for k, v in zip(CALL_ORDER, DEFAULTS)],
            "changedPropIds": ["export-csv.n_clicks"]}
    assert client.post("/_dash-update-component", json=body).status_code in (200, 204)


def test_callback_in_worker_thread_plots_pole_orientation(wsgi):
    """Production servers handle requests in threads; J parallel to N must plot with finite data there."""
    dep = next(d for d in wsgi.test_client().get("/_dash-dependencies").get_json() if "mod-graph.figure" in d["output"])
    values = [10, 0.6, 0.3, 0, 0, 0, PI / 2] + DEFAULTS[7:]
    body = {"output": dep["output"], "outputs": callback_outputs(dep),
            "inputs": [{"id": k, "property": "value", "value": v} for k, v in zip(CALL_ORDER, values)],
            "changedPropIds": ["theta_J-slider.value"]}
    out = {}
    t = threading.Thread(target=lambda: out.update(r=wsgi.test_client().post("/_dash-update-component", json=body)))
    t.start()
    t.join()
    resp = out["r"].get_json()["response"]
    assert resp["notice"]["className"] in ("notice is-empty", "notice notice-info")
    assert "not finite" not in notice_text(resp)
    assert resp["amp-graph"]["figure"]["data"]


def _post_reset(client, dash_deps, trigger_id, groups_clicked=(), presets_clicked=()):
    dep = next(d for d in dash_deps if '"reset-all"' in json.dumps(d["inputs"]))
    groups = [key for key, *_ in GROUPS]
    body = {"output": dep["output"], "outputs": callback_outputs(dep),
            "inputs": [{"id": "reset-all", "property": "n_clicks", "value": int(trigger_id == "reset-all")},
                       [{"id": {"group": g, "type": "reset-group"}, "property": "n_clicks",
                         "value": int(g in groups_clicked)} for g in groups],
                       [{"id": {"name": n, "type": "preset"}, "property": "n_clicks",
                         "value": int(n in presets_clicked)} for n in sorted(PRESETS)]],
            "changedPropIds": [trigger_id + ".n_clicks"]}
    r = client.post("/_dash-update-component", json=body)
    assert r.status_code == 200
    return {k: v["value"] for k, v in r.get_json()["response"].items()}


def test_reset_all_restores_every_default(client, dash_deps):
    out = _post_reset(client, dash_deps, "reset-all")
    assert out == {p.id: p.value for p in PARAMS}


def test_reset_group_touches_only_that_group(client, dash_deps):
    trig = json.dumps({"group": "precession", "type": "reset-group"}, separators=(",", ":"))
    out = _post_reset(client, dash_deps, trig, groups_clicked=("precession",))
    precession = next(params for key, _, _, params in GROUPS if key == "precession")
    assert out == {p.id: p.value for p in precession}


@pytest.mark.parametrize("name", sorted(PRESETS))
def test_presets_set_orientation(client, dash_deps, name):
    trig = json.dumps({"name": name, "type": "preset"}, separators=(",", ":"))
    out = _post_reset(client, dash_deps, trig, presets_clicked=(name,))
    assert out == PRESETS[name]["values"]


def test_csv_export_matches_computation(client, dash_deps):
    dep = next(d for d in dash_deps if "download.data" in d["output"])
    body = {"output": dep["output"], "outputs": {"id": "download", "property": "data"},
            "inputs": [{"id": "export-csv", "property": "n_clicks", "value": 1}],
            "state": [{"id": k, "property": "value", "value": v} for k, v in zip(CALL_ORDER, DEFAULTS)],
            "changedPropIds": ["export-csv.n_clicks"]}
    data = client.post("/_dash-update-component", json=body).get_json()["response"]["download"]["data"]
    assert data["filename"] == "regular_precession_Mc10_q0.6_z0.3.csv"
    text = data["content"]
    for p in PARAMS:
        assert p.name in text
    lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    cols = lines[0].split(",")
    table = np.loadtxt(io.StringIO("\n".join(lines[1:])), delimiter=",")
    rows = dict(zip(cols, table.T))
    r = physics.compute(*DEFAULTS)
    np.testing.assert_array_equal(rows["f_Hz"], r.f)
    np.testing.assert_array_equal(rows["phi_L_RP"], r.phi_L)
    np.testing.assert_array_equal(rows["cos_theta_L_RP"], r.cos_L)
    np.testing.assert_array_equal(rows["amp_RP_over_NP_minus_1"], r.amp)
    np.testing.assert_array_equal(rows["phase_RP_minus_NP_rad"], r.phase)
