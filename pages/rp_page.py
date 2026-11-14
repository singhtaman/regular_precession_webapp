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

physics = None      # webapp.physics, webapp.time_domain and webapp.sphere,
time_domain = None  # imported by load_model() once SciPy and Astropy are loaded
sphere = None


def load_model():
    """Import the model; returns the frequency-grid constants."""
    global physics, time_domain, sphere
    from webapp import physics as p, time_domain as td, sphere as sp
    physics, time_domain, sphere = p, td, sp
    return json.dumps({"fMin": physics.F_MIN, "deltaF": physics.DELTA_F,
                       "taperAlpha": time_domain.TAPER_ALPHA})


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


def evaluate(values, domain="frequency"):
    """Validate, then compute derived quantities and the RP/NP comparison for values (CALL_ORDER).

    In the time domain the result also holds h(t) for RP and NP (webapp.time_domain) and the
    J-frame directions of L and N with the time map t(f) (webapp.sphere).
    """
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
    if domain == "time" and r.f.size:
        try:
            td = time_domain.strain(*values)
            sp = sphere.directions(*values)
        except Exception as exc:
            out["error"] = f"{type(exc).__name__}: {exc}"
            return out
        out.update(t=_f64(td.t), h_RP=_f64(td.h_RP), h_NP=_f64(td.h_NP), fs=td.fs, nt=int(td.t.size))
        out.update(sp_t=_f64(sp.t), sp_theta=_f64(sp.theta_LJ), sp_phi=_f64(sp.phi_LJ),
                   sp_N=[float(v) for v in sp.N], sp_ldotn_error=sp.ldotn_error,
                   sp_t_monotonic=sp.t_monotonic)
    return out


def _header(buf, values):
    buf.write("# Regular precession explorer export\n# Parameters:\n")
    for pid, v in zip(CALL_ORDER, values):
        p = BY_ID[pid]
        buf.write(f"#   {p.name} ({p.symbol}{('_' + p.sub) if p.sub else ''}) = {v!r}"
                  f"{' ' + p.unit if p.unit else ''}\n")


def export_csv(values, domain="frequency"):
    """CSV of every plotted array; the frequency-domain file matches the Dash app's export."""
    if _validate(values):
        return None
    r = physics.compute(*values)
    if r.f.size == 0:
        return None
    if domain == "time":
        return _export_time(values)
    buf = io.StringIO()
    _header(buf, values)
    buf.write(f"# f_cut = {r.f_cut!r} Hz; frequency grid {physics.F_MIN} Hz to f_cut, step {physics.DELTA_F} Hz\n")
    cols = ["f_Hz", "phi_L_RP", "cos_theta_L_RP", "amp_RP_over_NP_minus_1", "phase_RP_minus_NP_rad"]
    data = np.column_stack([r.f, r.phi_L, r.cos_L, r.amp, r.phase])
    np.savetxt(buf, data, delimiter=",", header=",".join(cols), comments="", fmt="%.17g")
    name = "regular_precession_Mc{:g}_q{:g}_z{:g}.csv".format(*values[:3])
    return [name, buf.getvalue()]


def _export_time(values):
    td = time_domain.strain(*values)
    buf = io.StringIO()
    _header(buf, values)
    buf.write(f"# Time domain: inverse FFT of precessing_strain(f) on {physics.F_MIN} Hz to f_cut = {td.f_cut!r} Hz, "
              f"step {physics.DELTA_F} Hz, Tukey taper alpha = {time_domain.TAPER_ALPHA}; "
              f"sampling rate {td.fs:g} Hz; t_c = 0\n")
    data = np.column_stack([td.t, td.h_RP, td.h_NP])
    np.savetxt(buf, data, delimiter=",", header="t_s,h_RP,h_NP", comments="", fmt="%.17g")
    name = "regular_precession_time_Mc{:g}_q{:g}_z{:g}.csv".format(*values[:3])
    return [name, buf.getvalue()]
