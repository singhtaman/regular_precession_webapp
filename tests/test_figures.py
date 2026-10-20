"""Presentation rules: data is plotted unchanged and never pushed out of view."""
import math

import numpy as np
import pytest

from webapp import figures, physics
from webapp.params import BY_ID, CALL_ORDER

PI = np.pi
DEFAULTS = [BY_ID[k].value for k in CALL_ORDER]


def _finite(a):
    a = np.asarray(a, float)
    return a[np.isfinite(a)]


def _within(values, rng):
    v = _finite(values)
    return v.size == 0 or (v.min() >= rng[0] - 1e-12 and v.max() <= rng[1] + 1e-12)


@pytest.fixture(scope="module")
def default_result():
    return physics.compute(*DEFAULTS)


def test_default_windows_equal_original_app(default_result):
    """With the original app's defaults the axes are exactly the original fixed windows."""
    mod = figures.modulation_figure(default_result)
    assert list(mod.layout.xaxis.range) == pytest.approx([PI / 4, 3 * PI / 4])
    assert list(mod.layout.yaxis.range) == pytest.approx([-0.8, 0.8])
    amp = figures.amplitude_figure(default_result)
    assert list(amp.layout.yaxis.range) == pytest.approx([-1.05, 0.95])
    assert list(amp.layout.xaxis.range) == pytest.approx([math.log10(19), math.log10(150)])
    assert list(figures.phase_figure(default_result).layout.yaxis.range) == pytest.approx([-3.25, 0.75])


@pytest.mark.parametrize("args", [
    DEFAULTS,
    [10, 0.6, 0.3, PI / 4, 0, PI / 4, 0] + DEFAULTS[7:],          # face-on: L far from the original window
    [5, 1, 0.05] + DEFAULTS[3:],                                  # f_cut ≈ 364 Hz, beyond the 150 Hz window
    [10, 0.6, 0.3, PI / 2, 0, PI / 16, PI] + DEFAULTS[7:],         # φ_L wraps across ±π
    [10, 0.6, 0.3, PI / 4, 0, PI / 2, PI / 2, 10, 4, PI],          # widest cone, fastest precession
])
def test_all_data_inside_axes(args):
    r = physics.compute(*args)
    mod = figures.modulation_figure(r)
    for t in mod.data:
        assert _within(t.x, mod.layout.xaxis.range) and _within(t.y, mod.layout.yaxis.range), t.name
    for fig in (figures.amplitude_figure(r), figures.phase_figure(r)):
        xr = [10 ** v for v in fig.layout.xaxis.range]
        for t in fig.data:
            assert _within(t.x, xr), t.name
            if figures.amp_peak(r) is None:
                assert _within(t.y, fig.layout.yaxis.range), t.name


def test_plotted_values_are_the_computed_arrays(default_result):
    r = default_result
    amp = {t.name: t for t in figures.amplitude_figure(r).data}
    phase = {t.name: t for t in figures.phase_figure(r).data}
    np.testing.assert_array_equal(amp["RP"].y, r.amp)
    np.testing.assert_array_equal(phase["RP"].y, r.phase)
    np.testing.assert_array_equal(phase["RP"].x, r.f)
    mod = {t.name: t for t in figures.modulation_figure(r).data}
    # Only NaN gaps may be added (at φ_L wraps); the points themselves are unchanged.
    np.testing.assert_array_equal(_finite(mod["RP"].x), r.phi_L)
    np.testing.assert_array_equal(_finite(mod["RP"].y), r.cos_L)


def test_phi_wrap_is_broken_not_streaked():
    phi = np.array([3.0, 3.1, -3.1, -3.0])
    x, y = figures._break_wraps(phi, np.arange(4.0))
    assert np.isnan(x[2]) and np.isnan(y[2]) and _finite(x).tolist() == phi.tolist()


@pytest.mark.parametrize("phi_J,expected", [(PI / 2, PI / 2), (3 * PI / 2, -PI / 2), (2 * PI, 0.0)])
def test_J_marker_uses_phi_L_convention(phi_J, expected):
    r = physics.compute(*(DEFAULTS[:6] + [phi_J] + DEFAULTS[7:]))
    j = next(t for t in figures.modulation_figure(r).data if t.name == "J")
    assert j.x[0] == pytest.approx(expected, abs=1e-12)


def test_stationary_L_is_drawn_as_a_point():
    r = physics.compute(*(DEFAULTS[:7] + [0, 2, 0]))  # θ̃ = 0: L stays along J
    assert figures.is_point(r.phi_L, r.cos_L)
    rp = next(t for t in figures.modulation_figure(r).data if t.name == "RP")
    assert rp.mode == "markers"


def test_J_marker_is_drawn_under_the_curve(default_result):
    assert [t.name for t in figures.modulation_figure(default_result).data] == ["J", "RP"]


def test_zero_line_is_drawn_below_the_curves(default_result):
    for fig in (figures.amplitude_figure(default_result), figures.phase_figure(default_result)):
        assert [s.layer for s in fig.layout.shapes] == ["below"]


def test_divergent_amplitude_is_capped_and_labelled():
    r = physics.compute(*([10, 0.6, 0.3, 0, 0, PI / 2, PI / 4] + DEFAULTS[7:]))  # B_NP ≈ 0
    assert figures.amp_peak(r) > 1e10
    fig = figures.amplitude_figure(r)
    lo, hi = fig.layout.yaxis.range
    assert -figures.AMP_DISPLAY_LIMIT <= lo and hi <= figures.AMP_DISPLAY_LIMIT
    assert any("Off scale" in a.text for a in fig.layout.annotations)


def test_rp_is_distinct_from_j():
    assert figures.RP_COLOR.lower() != figures.J_COLOR.lower()
