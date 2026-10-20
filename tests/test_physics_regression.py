"""The redesigned app must reproduce RP_map_all_app.py exactly.

Every array plotted by the original plot_combined() is compared bit-for-bit with
webapp.physics.compute(), over named cases, slider extremes, a seeded random sample
and the grid of special orientation angles where the model has singular points.
"""
import itertools
import random
import threading

import numpy as np
import pytest

from webapp import physics
from webapp.params import BY_ID, CALL_ORDER

PI = np.pi
DEFAULTS = [BY_ID[k].value for k in CALL_ORDER]

NAMED = {
    "defaults": DEFAULTS,
    "face_on": [10, 0.6, 0.3, PI / 4, 0, PI / 4, 0, 5, 2, 0],
    "lowest_mass_widest_band": [5, 1, 0.05] + DEFAULTS[3:],
    "out_of_band": [50, 1, 1] + DEFAULTS[3:],
    "no_precession": DEFAULTS[:7] + [0, 0, 0],
    "theta_tilde_zero": DEFAULTS[:7] + [0, 2, 0],
    "omega_tilde_zero": DEFAULTS[:7] + [5, 0, 0],
    "all_maxima": [50, 1, 1, PI, 2 * PI, PI, 2 * PI, 10, 4, 2 * PI],
    "all_minima": [5, 0.05, 0.05, 0, 0, 0, 0, 0, 0, 0],
    "pole_J_parallel_N": [10, 0.6, 0.3, 0, 0, 0, PI / 2] + DEFAULTS[7:],
    "bnp_vanishes": [10, 0.6, 0.3, 0, 0, PI / 2, PI / 4] + DEFAULTS[7:],
}


def _grid(lo, hi, step, rng):
    return lo + step * rng.randint(0, round((hi - lo) / step))


def _random_cases(n, seed):
    rng = random.Random(seed)
    return [[_grid(p.min, p.max, p.step, rng) for p in (BY_ID[k] for k in CALL_ORDER)] for _ in range(n)]


def _special_angle_cases():
    th = [0, PI / 4, PI / 2, 3 * PI / 4, PI]
    ph = [k * PI / 4 for k in range(9)]
    cases = [DEFAULTS[:3] + [a, b, c, d] + DEFAULTS[7:] for a, b, c, d in itertools.product(th, ph, th, ph)]
    return cases[::3]  # every third of the 2025 combinations keeps the run short


CASES = list(NAMED.items()) \
    + [(f"random_{i}", c) for i, c in enumerate(_random_cases(150, seed=2024))] \
    + [(f"angles_{i}", c) for i, c in enumerate(_special_angle_cases())]


def _evaluate(fn, args):
    try:
        return fn(*args), None
    except Exception as exc:  # compared below: both versions must fail the same way
        return None, type(exc)


@pytest.mark.parametrize("name,args", CASES, ids=[n for n, _ in CASES])
def test_compute_matches_original_bit_for_bit(original, name, args):
    fig, orig_err = _evaluate(original.plot_combined, args)
    r, new_err = _evaluate(physics.compute, args)

    assert orig_err == new_err, f"original raised {orig_err}, new raised {new_err}"
    if orig_err:
        return

    mine = [(r.phi_L, r.cos_L), ([r.J[0]], [r.J[1]]), (r.f, r.amp), (r.f, r.phase)]
    assert len(fig.data) == len(mine)
    for (x, y), trace in zip(mine, fig.data):
        np.testing.assert_array_equal(np.asarray(x, float), np.asarray(trace.x, float))
        np.testing.assert_array_equal(np.asarray(y, float), np.asarray(trace.y, float))


def test_pole_orientation_evaluates():
    """theta_S = theta_J = 0 (J parallel to N on the detector z-axis) is handled by the face-on branch."""
    r = physics.compute(*NAMED["pole_J_parallel_N"])
    for a in (r.cos_L, r.phi_L, r.amp, r.phase):
        assert np.all(np.isfinite(a))


@pytest.mark.parametrize("args", [NAMED["defaults"], NAMED["face_on"], NAMED["out_of_band"], NAMED["all_minima"]])
def test_derived_quantities_match_original_readouts(original, args):
    mc, q, z, tS, pS, tJ, pJ = args[:7]
    d = physics.derived_quantities(mc, q, z, tS, pS, tJ, pJ)
    assert original.update_total_mass(mc, q, z) == (
        f"Total mass of the binary is, M = {d.total_mass:.2f} solar masses (redshifted, "
        f"M_z = {d.total_mass_redshifted:.2f} solar masses), and the cut off frequency for this "
        f"signal is, f_cut = {d.f_cut:.2f} Hz")
    assert original.update_cos_i_JN(tJ, tS, pJ, pS) == \
        f"Cosine of angle between J and N, cos(ι_JN) : {d.cos_i_JN:.2f}"
    assert d.luminosity_distance == original.redshift_to_luminosity_distance(z)


def test_build_params_does_not_mutate_library_defaults():
    from systems_lib import default_precession_params_sys2_NP as lib
    before = dict(lib)
    physics.compute(*NAMED["defaults"])
    assert lib == before


def test_model_error_state_is_enforced_in_worker_threads(monkeypatch):
    """np.seterr is per-thread; web servers evaluate callbacks in worker threads."""
    assert physics.MODEL_ERRSTATE["invalid"] == "raise"
    monkeypatch.setattr(physics, "_compute", lambda *args: np.geterr())
    seen = {}

    def worker():
        seen["thread_default"] = np.geterr()["invalid"]
        seen["during_compute"] = physics.compute(*NAMED["defaults"])["invalid"]

    t = threading.Thread(target=worker)
    t.start()
    t.join()
    assert seen["thread_default"] != "raise"  # otherwise this test proves nothing
    assert seen["during_compute"] == "raise"
