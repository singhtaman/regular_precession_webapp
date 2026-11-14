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
  `systems_lib.py`, `webapp/params.py` and `webapp/physics.py` unchanged, and calls them
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
  h(t) = 2 Re ∫₀^∞ h̃(f) e^{−2πift} df. t_c = 0, so the signal reaches f_cut at t = 0.
- Sampling rate: the smallest power of two with ≥ 6 samples per cycle at f_cut (1024 Hz at the
  defaults). The plot opens on the last 0.5 s; a range bar shows the whole signal.
- The export in the time view has columns t_s, h_RP, h_NP.

Checks done in Pyodide: the instantaneous frequency rises from ≈20 Hz to near f_cut at t ≈ 0;
signal lengths are 3.5 s at the defaults and 12.7 s for the longest slider setting (M_c = 5,
q = 1, z = 0.05), inside the 20 s window; a forward FFT of h(t) recovers the tapered h̃(f) in band
to 0.1–0.4% (the remainder is the near-zero stretch before the signal, which the plot omits).
There is no Dash/Vercel version of this view yet.

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
