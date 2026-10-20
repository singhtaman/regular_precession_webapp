"""Every scientific control keeps the id, range, step, marks and default of the original app."""
import numpy as np
import pytest
from dash import dcc

from webapp.params import BY_ID, CALL_ORDER, PARAMS, PRESETS


def _walk(component):
    yield component
    children = getattr(component, "children", None)
    if children is None:
        return
    for child in children if isinstance(children, (list, tuple)) else [children]:
        if hasattr(child, "to_plotly_json"):
            yield from _walk(child)


@pytest.fixture(scope="module")
def original_sliders(original):
    return {c.id: c for c in _walk(original.app.layout) if isinstance(c, dcc.Slider)}


def test_same_set_of_sliders(original_sliders):
    assert set(original_sliders) == {p.id for p in PARAMS} == set(CALL_ORDER)


@pytest.mark.parametrize("pid", sorted(BY_ID))
def test_range_step_and_marks_preserved(original_sliders, pid):
    o, p = original_sliders[pid], BY_ID[pid]
    assert (o.min, o.max) == (p.min, p.max)
    assert o.step == p.step
    assert {float(k): v for k, v in o.marks.items()} == {float(k): v for k, v in p.marks.items()}


@pytest.mark.parametrize("pid", sorted(BY_ID))
def test_defaults_unchanged(original_sliders, pid):
    assert BY_ID[pid].value == original_sliders[pid].value


def test_call_order_matches_original_callback(original):
    """CALL_ORDER is the argument order of the original update_plot callback."""
    deps = original.app.callback_map["combined-plot.figure"]["inputs"]
    assert [d["id"] for d in deps] == CALL_ORDER


@pytest.mark.parametrize("name", sorted(PRESETS))
def test_presets_lie_on_slider_grid(name):
    for pid, v in PRESETS[name]["values"].items():
        p = BY_ID[pid]
        k = (v - p.min) / p.step
        assert p.min <= v <= p.max and np.isclose(k, round(k)), (pid, v)


def test_presets_match_systems_lib():
    import systems_lib as lib
    for name, system in [("edge-on", lib.default_precession_params_sys2_RP), ("face-on", lib.default_precession_params_sys1_RP)]:
        vals = PRESETS[name]["values"]
        for key in ("theta_S", "phi_S", "theta_J", "phi_J"):
            assert vals[f"{key}-slider"] == system[key], (name, key)
