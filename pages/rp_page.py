"""Bridge between the browser build and the evaluation layer (runs in Pyodide).

The physics is the same code the Dash app uses: webapp.physics and webapp.params,
which call the unmodified Regular_precession class. This module only converts
inputs and outputs for JavaScript. _validate and export_csv mirror
webapp/callbacks.py, which can't be imported here because it depends on Dash.
"""
import io
import json

import numpy as np

from webapp.params import BY_ID, CALL_ORDER, GROUPS, PRESETS

physics = None  # webapp.physics, imported by load_model() once SciPy and Astropy are loaded


def load_model():
    """Import the model; returns the frequency-grid constants."""
    global physics
    from webapp import physics as module
    physics = module
    return json.dumps({"fMin": physics.F_MIN, "deltaF": physics.DELTA_F})


def params_json():
    """Slider definitions and presets, for building the controls (needs only NumPy)."""
    groups = [{
        "key": key, "title": title, "desc": desc,
        "params": [{
            "id": p.id, "name": p.name, "symbol": p.symbol, "sub": p.sub, "hint": p.hint,
            "min": p.min, "max": p.max, "step": p.step, "value": p.value, "unit": p.unit,
            "angle": p.angle, "marks": [[float(v), label] for v, label in p.marks.items()],
        } for p in params],
    } for key, title, desc, params in GROUPS]
    return json.dumps({"groups": groups, "callOrder": CALL_ORDER, "presets": PRESETS})


def _validate(values):
    """Return an error message for missing or out-of-range values, else None."""
    for pid, v in zip(CALL_ORDER, values):
        p = BY_ID[pid]
        if v is None:
            return f"Enter a value for {p.name.lower()}."
        if not (p.min - 1e-9 <= v <= p.max + 1e-9):
            return f"{p.name} must be between {p.min:g} and {p.max:g}{' ' + p.unit if p.unit else ''}."
    return None


def _f64(a):
    # Raw float64 bytes keep NaN, ±inf and full precision on the way to JavaScript.
    return np.ascontiguousarray(a, dtype=np.float64).tobytes()


def evaluate(values):
    """Validate, then compute derived quantities and the RP/NP comparison for values (CALL_ORDER)."""
    msg = _validate(values)
    if msg:
        return {"invalid": msg}

    d = physics.derived_quantities(*values[:7])
    out = {"derived": {"total_mass": float(d.total_mass),
                       "total_mass_redshifted": float(d.total_mass_redshifted),
                       "f_cut": float(d.f_cut), "cos_i_JN": float(d.cos_i_JN),
                       "luminosity_distance": float(d.luminosity_distance)}}
    try:
        r = physics.compute(*values)
    except Exception as exc:  # reported in the page, as the Dash app does
        out["error"] = f"{type(exc).__name__}: {exc}"
        return out

    out.update(f_cut=float(r.f_cut), J=[float(r.J[0]), float(r.J[1])], n=int(r.f.size),
               f=_f64(r.f), phi_L=_f64(r.phi_L), cos_L=_f64(r.cos_L),
               amp=_f64(r.amp), phase=_f64(r.phase))
    return out


def export_csv(values):
    """CSV of every plotted array, in the same format as the Dash app's export."""
    if _validate(values):
        return None
    r = physics.compute(*values)
    if r.f.size == 0:
        return None
    buf = io.StringIO()
    buf.write("# Regular precession explorer export\n# Parameters:\n")
    for pid, v in zip(CALL_ORDER, values):
        p = BY_ID[pid]
        buf.write(f"#   {p.name} ({p.symbol}{('_' + p.sub) if p.sub else ''}) = {v!r}"
                  f"{' ' + p.unit if p.unit else ''}\n")
    buf.write(f"# f_cut = {r.f_cut!r} Hz; frequency grid {physics.F_MIN} Hz to f_cut, step {physics.DELTA_F} Hz\n")
    cols = ["f_Hz", "phi_L_RP", "cos_theta_L_RP", "amp_RP_over_NP_minus_1", "phase_RP_minus_NP_rad"]
    data = np.column_stack([r.f, r.phi_L, r.cos_L, r.amp, r.phase])
    np.savetxt(buf, data, delimiter=",", header=",".join(cols), comments="", fmt="%.17g")
    name = "regular_precession_Mc{:g}_q{:g}_z{:g}.csv".format(*values[:3])
    return [name, buf.getvalue()]
