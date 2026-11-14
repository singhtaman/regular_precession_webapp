"""Time-domain strain: inverse Fourier transform of the frequency-domain model.

h̃(f) is Regular_precession.precessing_strain(f), unmodified, on the app's frequency
grid (F_MIN to f_cut in steps of DELTA_F), for the RP and NP models built exactly as
in physics.compute. Nothing here changes the model; the only additions are
presentation choices:

* The band edges are tapered with a Tukey window (TAPER_ALPHA) so the sharp start
  at F_MIN and cut-off at f_cut don't ring in the time domain.
* The transform uses the model's Fourier convention (Cutler & Flanagan 1994,
  h̃(f) = ∫ h(t) e^{2πift} dt), so h(t) = 2 Re ∫₀^∞ h̃(f) e^{−2πift} df.
* With t_c = 0 the signal ends at t = 0 (the cut-off frequency); the window covers
  1/DELTA_F = 20 s, from −(20 s − POST_TC) to POST_TC.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.signal.windows import tukey

from regular_precession_no_match import Regular_precession

from .physics import DELTA_F, F_MIN, MODEL_ERRSTATE, build_params

TAPER_ALPHA = 0.1   # fraction of the band tapered (half at each edge)
POST_TC = 0.25      # s shown after t_c
LEAD_IN = 0.25      # s shown before the signal enters the band
VISIBLE = 1e-3      # the plot starts where |h| first exceeds this fraction of its peak
SAMPLES_PER_CYCLE = 6


@dataclass
class TimeDomain:
    t: np.ndarray       # [s], t_c = 0
    h_RP: np.ndarray    # detector strain, dimensionless
    h_NP: np.ndarray
    fs: float           # sampling rate [Hz]
    f_cut: float        # [Hz]


def strain(chirp_mass, mass_ratio, redshift, theta_S, phi_S, theta_J, phi_J,
           theta_tilde, omega_tilde, gamma_P) -> TimeDomain:
    with np.errstate(**MODEL_ERRSTATE):
        params_NP, params_RP = build_params(
            chirp_mass, mass_ratio, redshift, theta_S, phi_S, theta_J, phi_J,
            theta_tilde, omega_tilde, gamma_P)
        RP = Regular_precession(params_RP)
        NP = Regular_precession(params_NP)

        f_cut = RP.get_f_cut()
        f = np.arange(F_MIN, f_cut, DELTA_F)
        if f.size == 0:
            empty = np.empty(0)
            return TimeDomain(empty, empty, empty, 0.0, f_cut)

        # Power-of-two sampling rate with at least SAMPLES_PER_CYCLE samples per cycle at f_cut, so
        # the curves are smooth when zoomed in. Above the band h̃ is zero, so a higher rate only
        # interpolates the same waveform.
        fs = float(2 ** math.ceil(math.log2(SAMPLES_PER_CYCLE * f_cut)))
        n = int(round(fs / DELTA_F))
        bins = np.rint(f / DELTA_F).astype(int)
        window = tukey(f.size, TAPER_ALPHA)

        def to_time(model):
            hf = np.zeros(n // 2 + 1, dtype=complex)
            hf[bins] = model.precessing_strain(f) * window
            # numpy's irfft uses e^{+2πi...}/n; conjugating gives the model's e^{−2πift} convention.
            return np.fft.irfft(np.conj(hf), n) * fs

        h_RP, h_NP = to_time(RP), to_time(NP)

    # Reorder the periodic window so time runs from −(T − POST_TC) to POST_TC:
    # sample i is at t_i = (i − (n − post))/fs, which is index (i + post) mod n of the transform.
    post = int(round(POST_TC * fs))
    h_RP, h_NP = np.roll(h_RP, -post), np.roll(h_NP, -post)
    t = (np.arange(n) - (n - post)) / fs

    # Drop the empty stretch before the signal enters the band.
    env = np.maximum(np.abs(h_RP), np.abs(h_NP))
    peak = np.nanmax(env)
    start = int(np.argmax(env > VISIBLE * peak)) if peak > 0 else 0
    start = max(0, start - int(LEAD_IN * fs))
    return TimeDomain(t[start:], h_RP[start:], h_NP[start:], fs, f_cut)
