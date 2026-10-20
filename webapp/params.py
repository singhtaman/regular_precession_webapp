"""Slider definitions for the explorer.

IDs, ranges, steps, marks and defaults are copied from RP_map_all_app.py so the
redesigned UI drives the calculation with exactly the same inputs.
"""
from dataclasses import dataclass, field

import numpy as np

PI = np.pi

ANGLE_MARKS_PI = {0: "0", PI / 4: "π/4", PI / 2: "π/2", 3 * PI / 4: "3π/4", PI: "π"}
ANGLE_MARKS_2PI = {0: "0", PI / 2: "π/2", PI: "π", 3 * PI / 2: "3π/2", 2 * PI: "2π"}


@dataclass(frozen=True)
class Param:
    id: str             # Dash component id of the slider (unchanged from the original app)
    name: str           # human-readable name
    symbol: str         # base symbol
    sub: str            # subscript ("" if none)
    hint: str           # one-line explanation shown under the label
    min: float
    max: float
    step: float
    value: float        # default
    marks: dict = field(default_factory=dict)
    unit: str = ""
    angle: bool = False


def _marks(values):
    return {v: str(v) for v in values}


GROUPS = [
    ("precession", "Precession", "Motion of L in a cone around J", [
        Param("theta_tilde-slider", "Precession amplitude", "θ̃", "",
              "Opening angle of the cone; it reaches θ̃/10 rad at the cut-off frequency",
              0, 10, 0.1, 5, _marks(range(11))),
        Param("omega_tilde-slider", "Precession frequency", "Ω̃", "",
              "Sets how fast L circles J; at 0, L does not circle J",
              0, 4.0, 0.1, 2.0, _marks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4])),
        Param("gamma_P-slider", "Precessional phase", "γ", "P",
              "Phase of the precession when the signal enters the band at 20 Hz",
              0, 2 * PI, PI / 16, 0.0, ANGLE_MARKS_2PI, "rad", True),
    ]),
    ("source", "Source", "Binary masses and distance", [
        Param("chirp_mass-slider", "Chirp mass", "M", "c",
              "Source frame; the waveform uses (1+z) times this value",
              5, 50, 1, 10, _marks(range(5, 51, 5)), "M☉"),
        Param("mass_ratio-slider", "Mass ratio", "q", "",
              "Sets the symmetric mass ratio η = q/(1+q)²",
              0.05, 1, 0.05, 0.6, _marks([0.2, 0.4, 0.6, 0.8, 1])),
        Param("redshift-slider", "Redshift", "z", "",
              "Sets the distance via flat ΛCDM (Planck 2018)",
              0.05, 1, 0.05, 0.3, _marks([0.2, 0.4, 0.6, 0.8, 1])),
    ]),
    ("orientation", "Orientation", "Detector-frame directions of N and J", [
        Param("theta_S-slider", "N polar angle", "θ", "S",
              "Polar angle of the line of sight N",
              0, PI, PI / 16, PI / 4, ANGLE_MARKS_PI, "rad", True),
        Param("phi_S-slider", "N azimuthal angle", "φ", "S",
              "Azimuthal angle of the line of sight N",
              0, 2 * PI, PI / 16, 0, ANGLE_MARKS_2PI, "rad", True),
        Param("theta_J-slider", "J polar angle", "θ", "J",
              "Polar angle of the total angular momentum J",
              0, PI, PI / 16, PI / 2, ANGLE_MARKS_PI, "rad", True),
        Param("phi_J-slider", "J azimuthal angle", "φ", "J",
              "Azimuthal angle of the total angular momentum J",
              0, 2 * PI, PI / 16, PI / 2, ANGLE_MARKS_2PI, "rad", True),
    ]),
]

PARAMS = [p for _, _, _, group in GROUPS for p in group]
BY_ID = {p.id: p for p in PARAMS}

# Argument order of the original update_plot / plot_combined callbacks.
CALL_ORDER = [
    "chirp_mass-slider", "mass_ratio-slider", "redshift-slider",
    "theta_S-slider", "phi_S-slider", "theta_J-slider", "phi_J-slider",
    "theta_tilde-slider", "omega_tilde-slider", "gamma_P-slider",
]

# Orientation presets taken from systems_lib.py (both lie on the slider grid; system 3,
# with θ_J = 8π/9, does not, so it is not offered).
PRESETS = {
    "edge-on": {"label": "Edge-on", "detail": ("System 2 · cos ι", "JN", " = 0"),
                "values": {"theta_S-slider": PI / 4, "phi_S-slider": 0.0,
                           "theta_J-slider": PI / 2, "phi_J-slider": PI / 2}},
    "face-on": {"label": "Face-on", "detail": ("System 1 · cos ι", "JN", " = 1"),
                "values": {"theta_S-slider": PI / 4, "phi_S-slider": 0.0,
                           "theta_J-slider": PI / 4, "phi_J-slider": 0.0}},
}
