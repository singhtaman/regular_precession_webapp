"""Shared fixtures. Tests run from the repository root: `pytest`."""
import json
import os
import sys
import warnings

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")


@pytest.fixture(scope="session")
def original():
    """The deployed Render app, imported unchanged, as the scientific reference."""
    return pytest.importorskip("RP_map_all_app",
                               reason="needs matplotlib (requirements-dev.txt) to import the original app")


@pytest.fixture(scope="session")
def wsgi():
    from app import app
    return app


@pytest.fixture()
def client(wsgi):
    return wsgi.test_client()


@pytest.fixture(scope="session")
def dash_deps(wsgi):
    return wsgi.test_client().get("/_dash-dependencies").get_json()


def callback_outputs(dep):
    return [{"id": o.split(".")[0], "property": o.split(".")[1]} for o in dep["output"].strip(".").split("...")]


@pytest.fixture()
def call_main(client, dash_deps):
    """POST the main plotting callback with slider values in CALL_ORDER; returns the response dict."""
    from webapp.params import CALL_ORDER

    dep = next(d for d in dash_deps if "mod-graph.figure" in d["output"])

    def call(values):
        body = {"output": dep["output"], "outputs": callback_outputs(dep),
                "inputs": [{"id": k, "property": "value", "value": v} for k, v in zip(CALL_ORDER, values)],
                "changedPropIds": ["chirp_mass-slider.value"]}
        r = client.post("/_dash-update-component", json=body)
        assert r.status_code == 200, r.data[:300]
        return r.get_json()["response"], len(r.data)

    return call


def notice_text(resp):
    return json.dumps(resp["notice"]["children"], ensure_ascii=False)
