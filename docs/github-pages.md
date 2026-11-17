# Regular Precession explorer: browser version on GitHub Pages

This branch (`feature/github-pages`) adds a version of the explorer that runs entirely in
the browser, so it can be hosted on GitHub Pages next to the Vercel app.

| | Vercel | GitHub Pages (this branch) |
|---|---|---|
| Code | `app.py` → `webapp/` (Dash) | `pages/` (static page) + the same Python model |
| Where the model runs | Vercel Function (Python 3.12) | the visitor's browser, in [Pyodide](https://pyodide.org) 0.29.0 (Python 3.13) |
| URL | https://regularprecession.vercel.app | https://singhtaman.github.io/regular_precession_webapp/ |

## How it works

- `pages/worker.js` starts Pyodide in a Web Worker, copies in `regular_precession_no_match.py`,
  `systems_lib.py`, `webapp/params.py` and `webapp/physics.py` unchanged (plus the browser-only
  `webapp/time_domain.py` and `webapp/sphere.py`), and calls them
  through `pages/rp_page.py`. Files are fetched one at a time with retries: Python's
  `http.server` refuses bursts of parallel connections, which showed up as "Failed to fetch".
- `pages/app.js` builds the controls from `webapp/params.py` and mirrors `webapp/callbacks.py`
  (validation, notices, system properties, presets, resets, CSV export in the same format).
- `pages/figures.js` is a JavaScript port of `webapp/figures.py`; Plotly is not available in
  Pyodide. Keep the two in step when changing a plot.
- The page reuses `webapp/assets/style.css`; `pages/pages.css` styles the native sliders and
  the first-load panel.

## Time-domain view (browser version only, for now)

A Frequency / Time switch in the results bar. In the time view the amplitude and phase panels
are replaced by one strain panel, h(t), with the RP and NP signals overlaid; the modulation of L
and the system properties are unchanged. `webapp/time_domain.py` computes it:

- h̃(f) = `Regular_precession.precessing_strain(f)` (unmodified) for RP and NP, on the app's grid
  (20 Hz to f_cut, Δf = 0.05 Hz, so a 20 s window).
- Tukey taper, α = 0.1 (5% of the band at each edge), so the sharp start and cut-off don't ring.
- Inverse FFT in the model's Fourier convention (Cutler & Flanagan 1994, h̃(f) = ∫ h(t) e^{2πift} dt):
  h(t) = 2 Re ∫₀^∞ h̃(f) e^{−2πift} df. t_c = 0, so t = 0 is the coalescence time; the signal reaches f_cut a few ms earlier
  (4.4 ms at the defaults, from t(f) = (1/2π) dΨ/df).
- Sampling rate: the smallest power of two with ≥ 6 samples per cycle at f_cut (1024 Hz at the
  defaults). The plot opens on the last 0.5 s; a range bar shows the whole signal.
- The export in the time view has columns t_s, h_RP, h_NP.

Checks done in Pyodide: the instantaneous frequency rises from ≈20 Hz to near f_cut at t ≈ 0;
signal lengths are 3.5 s at the defaults and 12.7 s for the longest slider setting (M_c = 5,
q = 1, z = 0.05), inside the 20 s window; a forward FFT of h(t) recovers the tapered h̃(f) in band
to 0.1–0.4% (the remainder is the near-zero stretch before the signal, which the plot omits).
There is no Dash/Vercel version of this view yet.

### Sphere: L, the orbital plane and N about J

A rotatable unit sphere, shown in both views. Frequency view: beside the modulation of L, above
the amplitude and phase panels, with the system properties below. Time view: a 2 × 2 grid of
modulation of L | system properties over sphere | strain, the strain filling the row's height. `webapp/sphere.py` computes it from the unmodified
class on the same frequency grid:

- Frame: z = J (north pole), x = the axis from which the model measures φ_LJ.
  L(f) = (sin θ_LJ cos φ_LJ, sin θ_LJ sin φ_LJ, cos θ_LJ) with `get_theta_LJ(f)`, `get_phi_LJ(f)`.
- Blue band: polar angles from min to max θ_LJ (20 Hz to f_cut); blue curve: the path of L.
- Orange band: latitudes within ±max θ_LJ of the equator, swept by the orbital plane (the great
  circle perpendicular to L, tilted by θ_LJ). Where the bands overlap (max θ_LJ > π/4, i.e.
  θ̃ ≳ 7.9) the overlap is purple and a "Both bands" key appears.
- N̂ = (0, sin ι_JN, cos ι_JN) with cos/sin ι_JN from `precession_angles()`: a dot on the sphere
  and an arrow. This is the direction implied by the model's `LdotN(f)` =
  sin θ_LJ sin ι_JN sin φ_LJ + cos θ_LJ cos ι_JN; `directions()` returns max |L·N − LdotN| on every
  call (`sp_ldotn_error`).
- Linked to time: hovering over (or tapping) the strain plot, or the amplitude or phase plot at a
  frequency f (converted to t(f)), moves the L marker, its arrow and the
  orbital-plane great circle to that time. Time → frequency uses the stationary-phase relation of
  the model's orbital phase, t(f) = (1/2π) dΨ/df with Ψ = `Regular_precession.Psi` (t_c = 0),
  differentiated numerically on the grid (`np.gradient`); the page inverts t(f) by binary search
  and interpolates f, θ_LJ, φ_LJ linearly between grid points. Times outside
  [t(20 Hz), t(f_cut)] show no marker. Before any hover the marker is at 20 Hz.
- Time evolution: **Play** runs the time from t(20 Hz) to t(f_cut) at 1×, ½×, ¼× (default) or
  ⅒× real time. While playing (and when paused or scrubbed afterwards) L leaves a trail on the
  sphere (the full path is faded), and the strain after the current time is covered by a
  translucent white band behind a cursor; if the strain view is zoomed in (narrower than 90% of
  the signal) it scrolls with the cursor, keeping it at 80% of the width. A time slider scrubs
  (and pauses); a click on the strain plot pauses at that time; hovering scrubs when paused.
  Reaching the end clears the trail and fade. In the frequency view the same cursor and fade
  mark the amplitude and phase plots at the current frequency (the frequency axis does not scroll).
  Switching views pauses.
- The strain panel now has a fixed bottom margin (118 px, `margin.autoexpand: false`): with
  automatic expansion Plotly 3.1 grew the range bar's margin by ~1 px on every relayout, which
  shrank the plot visibly during playback (and slightly on every zoom before).
- Plotly's 3D traces are not in plotly-basic, so `app.js` loads the cdnjs plotly.js 3.1.0 gl3d
  bundle (1.6 MB) once the model is ready, keeps it as a separate instance for the
  sphere only, and restores `window.Plotly` to the basic bundle.

Checks (14 Nov): L·N − LdotN ≤ 2.2×10⁻¹⁶ over 45 parameter sets (defaults, extremes, 40 random;
one random set had no band); t(f) strictly increasing in all of them. The strain's own zero-crossing
frequency at t(f) is 24.95, 40.18 and 81.76 Hz for f = 25, 40 and 80 Hz (defaults). In Chromium and
WebKit, the marker after hovering at four in-band times matches the model's L at the marker's
frequency to ≤ 4×10⁻⁶ and sits within half a strain sample (< 0.5 ms) of the hovered point; out-of-band
times hide the marker.

Playback checks (Chromium and WebKit): play at 1× reaches the end in ≈3 s at the defaults and
clears the trace; pause holds the marker and cursor; scrubbing pans the strain view to the cursor;
a click on the strain while playing pauses there; a slider change during playback continues with
the new signal; plot height stays constant during playback (320 px desktop, 268 px phone);
keyboard Space toggles Play. Headless frame rate during playback: 121 fps Chromium, 33 fps WebKit
(software GL).

The first visit downloads about 32 MB (Python, NumPy, SciPy, Astropy) and is ready in about
8 seconds on a fast connection; the browser caches it for later visits. Each update then takes
tens of milliseconds.

## Agreement with the Vercel app

The CSV export of 26 inputs (defaults, presets, slider extremes and 20 random slider
positions) was compared with the Vercel app's export:

- f, φ_L, cos θ_L and B_RP/B_NP − 1 agree to within 7×10⁻¹³ (relative): floating-point rounding.
- Φ_RP − Φ_NP agrees to within 5×10⁻⁶ rad (relative 2×10⁻⁷) in the worst case (θ̃ = 10,
  Ω̃ = 4) and to 10⁻¹³ rad in typical cases. The phase correction δφ is integrated with
  `scipy.integrate.odeint` at its default tolerances, and Pyodide ships SciPy 1.14.1 while the
  Vercel app pins 1.18.1; the difference is at the level of the integrator's tolerance.
- Inputs with no frequencies in band give no export in both.

## Build and preview locally

```bash
pages/build.sh                       # assembles _site/
python3 -m http.server -d _site 8000 # http://localhost:8000
```

## Deploy

`.github/workflows/pages.yml` runs `pages/build.sh` and publishes `_site/` on every push to
`feature/github-pages`, or by hand from the Actions tab. One-time setup in the repository
settings: Pages → Source: **GitHub Actions**, and allow `feature/github-pages` under
Environments → github-pages → Deployment branches.

Sphere in the frequency view (added 14 Nov, after the time view): `rp_page.evaluate` now returns
the sphere arrays in both views. Checked in headless Chromium: hovering the amplitude plot at 60 Hz
and the phase plot at 100 Hz places L at 59.6 and 99.6 Hz (pixel rounding); playback advances the
cursor on both plots on the log frequency axis; click to pause, scrub, play to the end, a slider
change during playback, switching views both ways, 390 px width without horizontal scroll, no
console errors.

## Vercel (browser version, separate project)

On this branch `vercel.json` describes the static browser build (no framework, no install,
`sh pages/build.sh`, output `_site/`), so a separate Vercel project with production branch
`feature/github-pages` serves the same site as GitHub Pages. The Dash app's Vercel project is
unaffected: it deploys `feature/vercel-ui`, which keeps its own Flask `vercel.json`. Note that the
Dash project also builds preview deployments of this branch, which now show the browser build.
