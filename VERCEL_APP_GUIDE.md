# Implementation guide: an interactive scientific Dash app on Vercel

A reusable blueprint for building, from scratch in a new repository, the kind of app
this repository contains: an interactive scientific explorer, where sliders drive a
Python model and Plotly figures update, designed to read as a polished tool, served
from a single Vercel Python function on the free Hobby plan, and able to run alongside
an existing deployment such as Render.

It distils what was built and learned on branch `feature/vercel-ui` of
`singhtaman/generic_precession_app`. The live reference is
https://genericprecession.vercel.app.

**How to read the labels.** Every recommendation is marked:

- **[Essential]**: needed for any app of this kind. Skipping it caused real bugs,
  failed builds or slow pages in this project.
- **[Project-specific]**: a choice or value particular to the generic-precession app.
  Replace it with your own.

**How to use this guide in a new repository.** Copy this file to `docs/` of the new
repo. Write that repo's `CLAUDE.md` from the template in section 14, then follow the
phases in section 13.

---

## Contents

1. [What you are building](#1-what-you-are-building)
2. [Architecture](#2-architecture)
3. [Repository layout](#3-repository-layout)
4. [The parameter registry](#4-the-parameter-registry)
5. [The model adapter](#5-the-model-adapter)
6. [Figures](#6-figures)
7. [Page structure](#7-page-structure)
8. [Callbacks and interaction](#8-callbacks-and-interaction)
9. [Styling system](#9-styling-system)
10. [API integrations, authentication and environment variables](#10-api-integrations-authentication-and-environment-variables)
11. [Vercel deployment](#11-vercel-deployment)
12. [Testing and verification](#12-testing-and-verification)
13. [Build plan](#13-build-plan)
14. [Adapting to another app, and a CLAUDE.md template](#14-adapting-to-another-app)
15. [Pitfalls](#15-pitfalls)

---

## 1. What you are building

A single-page app with:

- a **parameter panel** of grouped sliders with symbols, units, explanations, resets
  and presets;
- a **results area** with figures that update when a slider is released, a card of
  derived quantities, a legend, a status line and data export;
- **explicit states**: invalid input, "no data in range", model error, and results that
  are correct but easy to misread;
- **scientific fidelity**: the model code is unchanged, and the app's numbers are
  tested bit-for-bit against the original app;
- **deployment** as one Vercel Python function, with dependencies bundled for fast cold
  starts, beside the existing host.

| Aspect | Essential pattern | This project's instance [Project-specific] |
|---|---|---|
| Model | Pure Python functions or classes, imported unchanged | `Generic_precession` in `generic_precession_no_match.py` |
| Inputs | Registry of parameters (id, range, step, default, marks, unit) | 14 sliders in 4 groups |
| Outputs | Figures plus derived scalars | Modulation of L, amplitude, phase; 5 system properties |
| Per-request cost | Must stay far below the function limit | 10–50 ms per update |
| Existing deployment | Must keep working untouched | Render serves `scripts/GP_map_all_app.py` from `main` |

---

## 2. Architecture

### 2.1 Decision: keep Dash in a single WSGI function [Essential]

Choose in this order:

| Approach | Use when | Notes |
|---|---|---|
| **A. Dash as one Vercel Python function** | Each callback finishes in well under 300 s, and responses are under 4.5 MB | Used here. Smallest change, and the model stays in Python. |
| B. Static frontend plus a Python API | The UI needs features Dash can't provide, or the API serves other clients | Twice the moving parts. Only worth it with a clear reason. |
| C. Separate compute host | Single calculations take minutes, or need GPUs or long-running jobs | Vercel serves the UI only, and calls the compute host. |

A Dash app is a Flask (WSGI) app. Vercel's Python runtime (generally available, with
Python 3.12 the default) loads a top-level `app` WSGI callable from `app.py`, and
routes every path to it.

### 2.2 Request flow

```
Browser                                   Vercel                      Function (Python 3.12)
───────                                   ──────                      ──────────────────────
GET /                ───────────────────► CDN ──(miss)──────────────► Flask: Dash index HTML
GET /_dash-layout    ─────────────────────────────────────────────► layout.build()
GET /_dash-dependencies ──────────────────────────────────────────► callback graph
GET /_dash-component-suites/*, /assets/* ► CDN (HIT after the first request, CDN-Cache-Control)
slider released ─► POST /_dash-update-component ──────────────────► update(): validate → model → figures
angle slider dragged ─► clientside callback (runs in the browser, no request)
Export CSV ─► POST /_dash-update-component ───────────────────────► export(): model → CSV text → dcc.Download
```

### 2.3 Layering [Essential]

```
app.py                 entry point: exposes WSGI `app` for Vercel; `python app.py` runs locally
└─ webapp/main.py      creates Dash(), index template, meta tags, fonts, CDN cache header
   ├─ layout.py        page structure (pure functions that return components)
   ├─ callbacks.py     all server and client callbacks, validation, notices, export
   ├─ figures.py       Plotly figure builders and display rules (presentation only)
   ├─ physics.py       adapter: inputs → model calls → arrays (no presentation)
   ├─ params.py        single source of truth for every input control
   └─ assets/          style.css (tokens and layout), app.js (client helpers), images
model modules          the original scientific code, imported, never modified for the UI
```

Rules:

- Each layer depends only on layers below it. `figures.py` never calls the model;
  `physics.py` never builds UI.
- The model is not edited to suit the UI. If the model has to change, change it in its
  own commit, test it, and deploy it to every site that uses it (section 11.8).

---

## 3. Repository layout

```
.
├── app.py                    [Essential] Vercel entry point (top-level `app`)
├── webapp/                   [Essential] the new app, isolated from any old app
│   ├── __init__.py
│   ├── main.py  layout.py  callbacks.py  figures.py  physics.py  params.py
│   └── assets/ style.css  app.js  <logo image>
├── tests/                    [Essential] regression, params, figures, HTTP, optional UI
├── requirements.txt          [Essential] pinned runtime dependencies (Vercel installs these)
├── requirements-dev.txt      [Essential] -r requirements.txt + pytest, gunicorn, test-only packages
├── vercel.json               [Essential] framework preset + bundle exclusions
├── pytest.ini                [Essential] test paths, markers, warning filters
├── docs/                     deployment notes and this guide
├── <model modules>.py        [Project-specific] the science code
└── scripts/, Procfile        [Project-specific] the original app and its host's config, untouched
```

**[Essential]** Add the new app beside the old one; don't edit the old one. The branch
then only *adds* files, so the original deployment can't break. The only shared edit
here was one line appended to `.gitignore`.

**[Essential]** Don't add files that the old host reads automatically. For example,
Render reads `.python-version`, so pin Python elsewhere or rely on Vercel's default.

---

## 4. The parameter registry

One table describes every control. The layout, validation, resets, presets, export
headers and tests all read from it.

```python
# webapp/params.py  [Essential pattern]
from dataclasses import dataclass, field

@dataclass(frozen=True)
class Param:
    id: str            # Dash component id; keep the ORIGINAL app's id
    name: str          # human-readable name
    symbol: str        # base symbol, e.g. "θ"
    sub: str           # subscript, e.g. "S" (rendered with html.Sub)
    hint: str          # one-line explanation under the label
    min: float
    max: float
    step: float
    value: float       # default
    marks: dict = field(default_factory=dict)
    unit: str = ""
    angle: bool = False   # display as a fraction of π, snap to the step

GROUPS = [  # (key, title, description, [Param, ...])
    ("precession", "Precession", "Motion of L in a cone around J", [...]),
    ...
]
PARAMS = [p for *_, group in GROUPS for p in group]
BY_ID = {p.id: p for p in PARAMS}
CALL_ORDER = [...]   # argument order of the model adapter (== the original callback's Input order)
PRESETS = {"edge-on": {"label": ..., "detail": ..., "values": {...}}}
```

- **[Essential]** Copy ids, ranges, steps and marks from the original app, and test that
  they match (section 12). Record any changed default explicitly; here the nutation
  defaults Δθ̃ = 1.2, ΔΩ̃ = 0.2, ω̃ = 1.3 differ from the original app.
- **[Essential]** Presets must lie on the slider grid (`(v − min)/step` is an integer),
  or the slider and the value disagree.
- **[Project-specific]** The 14 parameters, the four groups, and the edge-on and
  face-on presets taken from `systems_lib.py`.

---

## 5. The model adapter

`webapp/physics.py` turns slider values into arrays. It is the only file that calls the
model.

### 5.1 Reproduce the original assembly exactly [Essential]

Copy the original app's parameter assembly line for line: unit conversions, derived
quantities, special-case substitutions such as `ω̃ = 0 → 0.001`, and order of
operations. Then **test it bit-for-bit** against the original (section 12.1). Make
only non-numeric changes, for example building each model object once instead of once
per panel.

```python
@dataclass
class Result:          # everything the figures and export need
    f: np.ndarray
    ...                # one field per plotted array, plus scalars such as f_cut

def compute(*values) -> Result: ...
def derived_quantities(*values) -> Derived: ...   # the original app's readouts
```

### 5.2 Enforce the model's NumPy error state in every request [Essential]

If the model (or anything it imports) calls `np.seterr(...)` at import time, note that
**NumPy error state is per-thread**. Web servers (Flask's development server, gunicorn
with threads, Vercel) handle requests in worker threads, where that call never ran. In
those threads, invalid operations silently produce NaN instead of raising, and NaN
plots as nothing. In this project, that caused blank panels with no error message.

```python
MODEL_ERRSTATE = np.geterr()          # captured right after importing the model modules

def compute(*values):
    with np.errstate(**MODEL_ERRSTATE):
        return _compute(*values)
```

Test it by calling `compute` from a `threading.Thread` (section 12.1).

### 5.3 Keep presentation out [Essential]

The adapter returns raw arrays. Anything visual (angle wrapping, NaN gaps, axis ranges)
belongs in `figures.py`, so exported data always equals the computation.

---

## 6. Figures

`webapp/figures.py`: one builder per panel, `f(result) -> go.Figure`.

### 6.1 Display rules [Essential]

1. **Open on the original app's axis windows, and widen only to fit the data.** Fixed
   windows hid data for most orientations here: the modulation curve was partly
   off-screen in 87 % of random inputs.

   ```python
   def _fit(default_lo, default_hi, *arrays, pad_frac=0.04):
       finite = [a[np.isfinite(a)] for a in arrays if np.size(a)]
       finite = [a for a in finite if a.size]
       if not finite:
           return [default_lo, default_hi]
       lo = min(a.min() for a in finite); hi = max(a.max() for a in finite)
       pad = pad_frac * max(hi - lo, 1e-9)
       return [min(default_lo, lo - pad), max(default_hi, hi + pad)]
   ```

2. **Cap axes where a quantity diverges, and say so.** If the ratio B_GP/B_NP − 1
   reaches 10¹⁵, stretching the axis flattens everything. Clamp the axis to a sensible
   limit (here ±20), and add an "Off scale: values reach …" annotation plus a note.
3. **Never let a curve disappear.**
   - Draw reference lines underneath the curves (`add_hline(..., layer="below")`).
   - Draw a stationary curve (zero extent) as a marker.
   - Draw the background marker (here the J direction) before the points that may sit
     on top of it.
4. **Periodic coordinates:** insert NaN where an angle wraps (|Δφ| > π), and express
   reference markers in the same convention (here `arctan2` gives (−π, π]).
5. **Add only gaps, never change values.** Test that each plotted array equals the
   computed one.
6. **Placeholder figures, not hidden figures.** In empty or error states, return a blank
   figure with a centred message. Never hide the graph container (see 9.4).

### 6.2 Presentation defaults [Essential]

- A shared `_base_layout()`: transparent paper background, white plot area, light grid,
  mirrored axis lines, outside ticks, 12–13 px Inter font, compact margins,
  `uirevision="keep"` (so a user's zoom survives updates), `showlegend=False`.
- Hover: `"x unified"` for curves against frequency; `"closest"` with formatted values
  for parametric plots.
- Axis labels with HTML subscripts (`B<sub>GP</sub>`) instead of MathJax, which avoids
  loading a large extra library.
- Graph config:
  `{"displaylogo": False, "responsive": True, "modeBarButtonsToRemove": ["select2d", "lasso2d"], "toImageButtonOptions": {"format": "png", "filename": ..., "scale": 3}}`.
- Colours carry meaning and stay distinct without colour vision. Here GP is blue and
  solid, RP is amber (`#d97a1c`) and dashed, and the J marker is red. Line colours
  should have at least 3:1 contrast against white.

**[Project-specific]** The three panels, their default windows, π-fraction ticks, the
log-frequency axis and its tick values, and the ±20 amplitude cap.

---

## 7. Page structure

### 7.1 Wireframe (desktop) [Essential pattern]

```
┌──────────────────────────────────────────────────────────────────────────────┐
│ [logo] App title · one-line description                Source  GitHub  Web  │  sticky header
├───────────────────────┬──────────────────────────────────────────────────────┤
│ PARAMETERS  Reset all │ ── GP   - - RP   ● J         status · Export CSV     │  results bar
│ ▸ Group 1             │ ┌ notice (hidden when empty) ───────────────────────┐ │
│   sym  Name   [value] │ └───────────────────────────────────────────────────┘ │
│   ──●────────── [box] │ ┌ Figure 1 ──────────────┐ ┌ System properties ────┐ │
│   hint text           │ │                         │ │ Quantity      value u │ │
│   Reset group         │ │                         │ │ ...                   │ │
│ ▸ Group 2 ...         │ └─────────────────────────┘ └───────────────────────┘ │
│ ▸ Group 4  [presets]  │ ┌ Figure 2 ──────────────┐ ┌ Figure 3 ─────────────┐ │
│ (sticky, scrolls)     │ └─────────────────────────┘ └───────────────────────┘ │
│                       │ ▸ About this figure and model                         │
├───────────────────────┴──────────────────────────────────────────────────────┤
│ © Author                                                 Funding acknowledgement│
└──────────────────────────────────────────────────────────────────────────────┘
```

### 7.2 Responsive behaviour [Essential]

| Viewport | Layout |
|---|---|
| ≥ 1180 px | Sidebar (360 px; 320 px below 1180 px) plus a 2×2 results grid |
| 880–1179 px | Sidebar plus results in one column |
| < 880 px | Single column: **results first**, then parameters. Sidebar no longer sticky; header links become icons. |
| < 560 px | Smaller plot frames (310 px), no frame padding |

**[Essential]** Layout choices learned from user feedback:

- A full-width plot of a parametric curve stretches it into a misleading shape. Use a
  2×2 grid and put the derived-quantity card beside the first figure.
- Put the legend once, in the results bar, not in each figure.

### 7.3 Accessibility [Essential]

- A skip link ("Skip to results").
- `aria-label` on regions, and `role="status"` / `role="alert"` on the status line and
  the notice.
- Visible `:focus-visible` rings.
- Native `<details>/<summary>` for collapsible groups and the About section; these are
  keyboard-accessible without JavaScript.
- Keyboard slider control (Dash 4 sliders support the arrow keys).
- `prefers-reduced-motion` support.

### 7.4 Component catalogue (`layout.py`) [Essential pattern]

| Function | Produces |
|---|---|
| `symbol(p)` | Symbol with HTML subscript |
| `control(p)` | Label row (symbol, name, unit pill, or a π-fraction badge for angles), `dcc.Slider`, hint |
| `group(key, title, desc, params)` | `html.Details` with summary, optional presets row, controls, "Reset {group}" link button |
| `plot_card(id, title, subtitle, filename)` | Card with header and a **fixed-height frame** containing `dcc.Graph(style={"height": "100%"})` |
| `stats_card()` | Card holding `#stats` (rows of label, note and value with unit) |
| `legend()` | Line and marker keys drawn in CSS |
| `about()` | `html.Details` explaining each panel, the frequency grid, cosmology and how to interact |
| `build()` | The whole page; assigned as `app.layout = build` (a function, so it is rebuilt per page load) |

Slider settings: `updatemode="mouseup"` (one request per release, not per pixel);
`tooltip={"placement": "bottom", "always_visible": False}`;
`allow_direct_input=not p.angle` (a number box for scalar sliders; angles use a
π-fraction badge instead).

---

## 8. Callbacks and interaction

### 8.1 The main callback [Essential]

```python
@app.callback(
    Output("fig-1", "figure"), Output("fig-2", "figure"), Output("fig-3", "figure"),
    Output("stats", "children"), Output("status", "children"),
    Output("notice", "children"), Output("notice", "className"),
    Output("export-csv", "disabled"),
    *[Input(pid, "value") for pid in CALL_ORDER],
)
def update(*values):
    if (msg := _validate(values)):                         # missing or out-of-range input
        return (no_update,) * 5 + _notice("warn", "Check the parameters", msg) + (True,)
    d = derived_quantities(*values)
    try:
        r = compute(*values)
    except Exception as exc:                               # never fail silently
        log.exception(...)
        blank = placeholder_figure("No plot for these parameters: see the message above")
        return (blank,) * 3 + (stats(d), "", *_notice("error", "...", f"{type(exc).__name__}: {exc}"), True)
    if r.f.size == 0:                                      # valid input, nothing to plot
        blank = placeholder_figure("No frequencies in band ...")
        return (blank,) * 3 + (stats(d), "", *_notice("warn", "No signal in the detector band", ...), True)
    notes = _notes(r)                                      # correct-but-surprising results
    return (fig1(r), fig2(r), fig3(r), stats(d), status_line(r),
            *(_notice("info", "About these results", notes) if notes else (None, "notice is-empty")),
            False)
```

- **[Essential]** One server round trip per slider release. The figures, stats, status,
  notice and export state all come from that one request.
- **[Essential]** Notes for results that look broken but aren't: non-finite values left
  out, a stationary curve, curves equal to the reference, a divergent ratio. Generate
  them from the data, not from assumptions about which parameters cause them.
- **[Project-specific]** The note texts, the 20 Hz lower frequency, and the status line
  format ("2,433 samples · 20–141.6 Hz, Δf = 0.05 Hz · 10 ms").

### 8.2 Resets and presets [Essential]

A single callback writes every slider value, triggered by `reset-all`, by
pattern-matching ids `{"type": "reset-group", "group": ALL}` and
`{"type": "preset", "name": ALL}`, with `prevent_initial_call=True`. Use
`ctx.triggered_id` to choose which values to set, and return `no_update` for the rest.

### 8.3 Export [Essential]

A `dcc.Download` fed by a callback that recomputes from `State` values and writes:

- a commented header with every parameter, its value and unit;
- one column per plotted array, formatted with `%.17g` so values round-trip exactly.

It never raises: on an exception it logs and returns `no_update`. The button is disabled
whenever there is nothing to export. Size check: the largest export here is 1.24 MB (measured),
under Vercel's 4.5 MB response limit.

### 8.4 Client-side helpers (`assets/app.js`) [Essential]

- **Display formatting without a server call:** a `dash_clientside` namespace function
  formats angle sliders as π-fractions, registered once per angle slider with
  `ClientsideFunction("gp", "angle")`. Read `callback_context.triggered` to prefer
  `value` over a stale `drag_value` after a reset.
- **Enter commits typed values:** Dash 4 slider number boxes commit on blur. Add a
  document-level `keydown` handler that blurs the box on Enter.

---

## 9. Styling system

One stylesheet, `webapp/assets/style.css`; Dash loads every `.css` file in `assets/`.

### 9.1 Design tokens [Essential]

```css
:root {
  --bg:#f4f6f9; --surface:#fff; --surface-2:#f8f9fb; --border:#e3e7ed; --border-strong:#cfd5de;
  --text:#1d2433; --text-2:#3b4456; --muted:#5d6677;
  --accent:#2457d6; --accent-strong:#1b44ad; --accent-weak:#eaf0fd;   /* also the GP colour */
  --rp:#d97a1c; --j:#d93a2b;
  --warn-bg:#fff8eb; --warn-border:#f1d49a; --warn-text:#6d4500;
  --error-bg:#fdf0ef; --error-border:#f2c3bd; --error-text:#8a2117;
  --font: Inter, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
  --mono: "JetBrains Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
  --radius:12px; --radius-sm:8px;
  --shadow: 0 1px 2px rgba(16,24,40,.04), 0 1px 3px rgba(16,24,40,.06);
  --space-1:4px; --space-2:8px; --space-3:12px; --space-4:16px; --space-5:24px; --space-6:32px;
  --header-h:64px; --sidebar-w:360px;
}
```

The Python figure colours (`GP_COLOR`, `RP_COLOR`, …) must match these tokens.
**[Project-specific]**: the palette values, the fonts (Inter and JetBrains Mono from
Google Fonts), and the logo.

### 9.2 Theming Dash 4 components [Essential]

Dash 4 components read `--Dash-*` CSS variables, which Dash injects on `:root` at
runtime, after your stylesheet loads. Override them with `html:root`: it is more
specific, so it wins without `!important`.

```css
html:root {
  --Dash-Fill-Interactive-Strong: var(--accent);   /* slider track and thumb, focus */
  --Dash-Fill-Interactive-Weak: var(--accent-weak);
  --Dash-Fill-Disabled: #e3e7ee;
  --Dash-Stroke-Strong: var(--border-strong);
  --Dash-Text-Primary: var(--text);
  --Dash-Text-Strong: var(--text-2);
  --Dash-Text-Weak: var(--muted);
}
```

Fine-tune with stable class names: `.dash-slider-mark`, `.dash-input-container`, and
`.dash-slider-root.has-marks`.

### 9.3 Component styles [Essential pattern]

Cards (`.plot-card`, `.about`, `.controls`), buttons (`.button`, `.button-ghost`,
`.link-button`, `.chip`, `:disabled`), notices (`.notice-warn`, `.notice-error`,
`.notice-info`, `.is-empty { display: none }`), stats rows (CSS grid: label/note on the
left, value on the right), and legend keys (CSS borders for solid and dashed lines, and
a dot). Icons are inline SVG data URIs used as CSS masks with `currentColor`
(`.icon-code`, `.icon-user`, `.icon-globe`, `.icon-download`), so there are no external
image requests.

### 9.4 Fixed-height plot frames [Essential]

```css
.plot-frame { height: 380px; padding: 0 var(--space-2) var(--space-2); }
/* dcc.Graph inside gets style={"height": "100%", "width": "100%"} */
```

**Why:** if a Plotly graph's container is ever `display:none`, Plotly measures it as
0 px. When it is shown again, the container stays collapsed (8 px here) while the SVG
overflows it, and the card clips it, so the plot sits hidden behind its own card. A
fixed-height frame, plus never hiding graphs (6.1 rule 6), prevents this. Test it with
real browser interaction (12.5).

---

## 10. API integrations, authentication and environment variables

### 10.1 API integrations

**This app has no external API integrations.** All computation happens in the function.
Its only HTTP surface is Dash's own internal endpoints:

| Endpoint | Method | Purpose |
|---|---|---|
| `/` | GET | Index HTML (custom `index_string`: favicon, font preconnects) |
| `/_dash-layout` | GET | Serialised layout |
| `/_dash-dependencies` | GET | Callback graph |
| `/_dash-update-component` | POST | Every server-side callback (JSON in, JSON out) |
| `/_dash-component-suites/...` | GET | Dash and React bundles (fingerprinted; CDN-cached) |
| `/assets/...` | GET | CSS, JS, images (CDN-cached) |

The only external resource is **Google Fonts** (stylesheet plus preconnect). It's
optional; the font stacks fall back to system fonts.

**If your app needs an external API [Essential guidance]:**

- Call it from a server-side callback, never from browser JavaScript with a secret.
- Keep keys in Vercel environment variables (Sensitive type), read with
  `os.environ[...]`, never committed.
- Set explicit timeouts well under the 300 s function limit, handle failures with a
  notice (8.1), and remember that cold starts add about 3 s.
- Don't rely on in-memory caches persisting between requests; instances are ephemeral.
  Use a managed store if results must persist.

### 10.2 Authentication

**This app has no authentication.** The production site is public by design; the model
and the page hold nothing private.

| Need | Option | Plan |
|---|---|---|
| Keep previews private | Vercel Authentication (default Deployment Protection for preview URLs) | Hobby |
| Check a protected preview from a terminal | `vercel curl <url>` | Hobby |
| Password-protect production | Vercel Password Protection | Paid add-on on Pro |
| Per-user login in the app | Flask/Dash auth middleware, or an identity provider; secrets in Vercel env vars | Any; adds complexity, so justify it first |

### 10.3 Environment variables

| Name | Where | Value | Required | Purpose |
|---|---|---|---|---|
| `PORT` | local | e.g. `8050` | no | Port for `python app.py` |
| `HOST` | local | e.g. `127.0.0.1` | no | Bind address for `python app.py` |
| `DASH_DEBUG` | local | `1` | no | Dash debug mode locally (never in production) |
| `VERCEL_SUPPORT_LARGE_FUNCTIONS` | Vercel project, all environments | `1` | **yes, if bundle > 250 MB** | Allows Large Functions (beta) |
| `VERCEL_PYTHON_DISABLE_RUNTIME_INSTALL` | Vercel project, Production and Preview | `1` | recommended | Bundles every dependency instead of installing at cold start (11.4) |

None of these are secrets. `vercel link` writes a `.env.local` holding a short-lived
`VERCEL_OIDC_TOKEN`. The app doesn't use it; delete it before deploying from an export,
and keep `.env*` out of git.

---

## 11. Vercel deployment

### 11.1 Entry point [Essential]

```python
# app.py
import os
from webapp.main import app as dash_app, server

app = server   # WSGI callable Vercel loads

if __name__ == "__main__":
    dash_app.run(host=os.environ.get("HOST", "127.0.0.1"),
                 port=int(os.environ.get("PORT", 8050)),
                 debug=os.environ.get("DASH_DEBUG") == "1")
```

Use `app.run` (`run_server` was removed in Dash 3). Don't manipulate `sys.path`: Vercel
runs from the project root, and relative paths resolve from there.

### 11.2 `vercel.json` [Essential]

```json
{
  "$schema": "https://openapi.vercel.sh/vercel.json",
  "framework": "flask",
  "functions": {
    "app.py": {
      "excludeFiles": "{scripts/**,tests/**,docs/**,**/*.md,**/__pycache__/**,<unused-heavy-modules>.py,Procfile,LICENSE}"
    }
  }
}
```

- `"framework": "flask"` is required. A project created as "Other" looks only for
  `api/` functions, and the build fails with "The pattern "app.py" defined in
  `functions` doesn't match any Serverless Functions".
- Exclude the old app, tests, docs and any module that imports packages you don't
  install.

### 11.3 Dependencies [Essential]

```
# requirements.txt: pinned to the versions you tested
dash==4.4.1
plotly==7.1.0
flask==3.1.3        # list Flask explicitly so Vercel detects the Flask preset
numpy==2.5.4
scipy==1.18.1
astropy==8.0.2      # [Project-specific]: needed only for the cosmology
```

- Install only what runtime code imports. Here matplotlib was imported but never used,
  so it was dropped.
- Check that Linux wheels exist and measure the real bundle:
  `uv pip install --python-version 3.12 --python-platform x86_64-manylinux_2_28 --target /tmp/bundle -r requirements.txt && du -sh /tmp/bundle`.
  Here that was 316 MB.
- Python version: Vercel's default (3.12) unless you need another; set it in a way the
  old host doesn't read.

### 11.4 Bundle size and cold starts [Essential]

Vercel's Python builder (`@vercel/python`) behaves like this. Its source was read
because this behaviour isn't in the public docs:

| Bundle (uncompressed) | Default behaviour |
|---|---|
| ≤ 225 MB | Everything bundled |
| > 225 MB | **Packages installed into `/tmp` at cold start** (the log shows "Runtime dependencies installed in ~3 s") |
| > 225 MB with `VERCEL_PYTHON_DISABLE_RUNTIME_INSTALL=1` and `VERCEL_SUPPORT_LARGE_FUNCTIONS=1` | Everything bundled as a Large Function (up to 5 GB, beta) |

Measured on fresh deployments of this app:

| Configuration | First request after idle | Warm request |
|---|---|---|
| Runtime install (default) | 6.9 s | 0.17 s |
| Everything bundled | 3.2 s | 0.14 s |

Removing test suites, stubs and source maps from packages only reached 263 MB here, so
trimming couldn't avoid the threshold. If Vercel drops the undocumented variable, builds
fall back to runtime install, which is slower but still works.

### 11.5 CDN caching [Essential]

```python
@server.after_request
def _cdn_cache(response):
    if response.status_code == 200 and request.path.startswith(("/_dash-component-suites/", "/assets/")):
        response.headers["CDN-Cache-Control"] = "public, max-age=86400"
    return response
```

Static bundles are then served from Vercel's CDN (`x-vercel-cache: HIT`) instead of
invoking the function, which saves invocations and CPU on Hobby.

### 11.6 Deploying [Essential]

**From the CLI, without pushing.** Deploy exactly what is committed:

```bash
vercel login                                             # in your own terminal (browser OAuth)
rm -rf /tmp/app-deploy && mkdir /tmp/app-deploy
git archive <branch> | tar -x -C /tmp/app-deploy         # committed files only
cd /tmp/app-deploy
vercel link --yes --project <project>                    # first time; keep .vercel/project.json for reuse
rm -f .env.local                                         # don't upload the OIDC token
vercel deploy            # preview (protected)   |   vercel deploy --prod   # production
```

The first deployment of a new project is assigned to production automatically. If Node
isn't installed, download the official Node LTS tarball (verify `SHASUMS256.txt`) into a
scratch folder and `npm install --prefix <dir> vercel`; nothing is installed
system-wide.

**Validate before deploying:** `vercel build` on a copy of the repo, with
`.vercel/project.json` set to `{"settings": {"framework": "flask"}}`. Note that a build
on macOS bundles macOS wheels; only Vercel's builders produce the Linux bundle.

### 11.7 Project, domains and settings

```bash
vercel project add <name>                                 # or create on first link
printf "1" | vercel env add VERCEL_SUPPORT_LARGE_FUNCTIONS production    # repeat: preview, development
printf "1" | vercel env add VERCEL_PYTHON_DISABLE_RUNTIME_INSTALL production   # and preview
vercel project rename <old> <new>
vercel domains add <name>.vercel.app <project>            # attach: future prod deploys use it
vercel alias set <deployment-url> <name>.vercel.app       # point the current deployment at it now
```

**[Project-specific]** Project `generic-precession`; main address
`genericprecession.vercel.app` (the older `generic-precession.vercel.app` and
`generic-precession-ui.vercel.app` still point to the same app); region `iad1`.

### 11.8 Running beside an existing host [Essential]

| | Existing host (e.g. Render) | Vercel |
|---|---|---|
| Branch | `main` | `feature/vercel-ui` |
| Code | old app in `scripts/` | `app.py` → `webapp/` |
| Shared | model modules (identical on both branches) | same |

- **Model changes go to both sites:** commit (or cherry-pick) the model change to
  `main`, which redeploys the old host. Then merge `main` into the feature branch (don't
  rebase) and redeploy Vercel. Before pushing `main`, run the old app under gunicorn
  from the directory its host uses, and make one real callback request.
- **Git integration (optional):** connect the repository to the Vercel project, then:
  - set **Production Branch** to the feature branch;
  - add an **Ignored Build Step**, so pushes to `main` (which has no `app.py`) don't
    create failing builds:
    `if [ "$VERCEL_GIT_COMMIT_REF" = "feature/vercel-ui" ]; then exit 1; else exit 0; fi`.
- **Check for other Vercel projects already connected to the repository.** Every push
  triggers their builds. Here an older project failed on every push.

### 11.9 Hobby limits that matter

| Limit | Hobby | Design for |
|---|---|---|
| Function duration | 300 s | Callbacks in milliseconds; time the worst case |
| Memory | 2 GB / 1 vCPU | Measure peak RSS (about 200 MB here) |
| Request/response body | 4.5 MB | Measure the largest callback response (1.06 MB here) |
| Bundle | 500 MB Python; Large Functions beyond | See 11.4 |
| Active CPU | 4 CPU-hours / month | About 0.05 s per update plus about 3 s per cold start |
| Invocations | 1,000,000 / month | One per slider release; static assets CDN-cached |
| Runtime logs | 1 hour retention | `vercel logs --deployment <url> --json` without `--since` |
| Use | Non-commercial, personal | |

---

## 12. Testing and verification

`pytest` runs everything; `pytest -m ui` runs the browser tests.

### 12.1 Scientific regression (`test_physics_regression.py`) [Essential]

- Import the **original app unchanged** as the reference. This needs its dependencies
  in `requirements-dev.txt`, here matplotlib.
- Parametrise over:
  - named cases (defaults, original defaults, presets, slider minima and maxima,
    out-of-range, known special cases);
  - a **seeded** random sample on the slider grid;
  - the grid of special angles (multiples of π/4). Degenerate geometry lives exactly
    on these grid points, and random sampling rarely hits them.
- Compare **every plotted array bit-for-bit** (`np.testing.assert_array_equal`). Where
  the original raises, require the same exception type from the new code.
- Compare derived readouts against the original callbacks' formatted strings.
- Test the error state in a worker thread (5.2) and that library defaults aren't mutated.

### 12.2 Controls (`test_params.py`) [Essential]

Walk the original `app.layout`, collect its `dcc.Slider`s, and assert the same ids,
ranges, steps and marks. List the intentional default changes explicitly. Check that
`CALL_ORDER` equals the original callback's Input order, and that presets lie on the
grid and match their source.

### 12.3 Figures (`test_figures.py`) [Essential]

- With the original defaults, the axes equal the original windows.
- For extreme cases, no finite data lies outside the axes.
- Plotted values equal the computed arrays (only NaN gaps added).
- Wrap gaps, the J marker convention, stationary-as-marker, reference line `layer="below"`.
- The capped axis and its off-scale annotation.

### 12.4 HTTP (`test_app.py`) [Essential]

Use the Flask test client:

- every endpoint returns 200; CDN headers are present on assets and absent on layout;
- the index has a title and viewport meta; the layout contains every id;
- every notice state;
- plots are never removed (placeholders with messages);
- the export button is disabled without data, and export equals `compute` exactly;
- the largest payload is under 4.5 MB;
- simulated model failures, using `monkeypatch` on the callback module's `compute`.

### 12.5 Browser (`test_ui.py`, optional) [Essential for UI changes]

Use Playwright with a threaded `werkzeug.serving.make_server` on port 0:

- the 2×2 arrangement;
- **after every real mouse drag, each graph's SVG fills its frame and stays inside its
  card**, including drags into and out of empty and error states;
- Enter commits typed values; presets update the readouts; no page or console errors.

Skip automatically when Playwright isn't installed.

### 12.6 Live checks after deploying [Essential]

- Request the page, layout and assets from the real URL.
- Send a **parameter sweep over HTTP** (about 1,300 requests at concurrency 4). Classify
  each response as plots, plots with a note, or the specific empty or error message. The
  target is zero unexplained responses.
- Compare **live figure arrays against local computation**. Decode Plotly's `bdata`
  typed arrays with base64 and the given `dtype`. Across platforms (Linux x86-64 vs
  macOS arm64), expect agreement near machine precision for closed-form quantities, and
  within the integrator tolerance for ODE-integrated ones (here phase, about 10⁻⁷).
  Bit-for-bit equality only holds on the same platform.
- Run the browser visibility test against the live URL at three widths.
- Search the function logs for "Runtime dependencies installed" (it should be absent;
  see 11.4).

---

## 13. Build plan

Follow these phases, which match the CLAUDE.md workflow:

1. **Inspect.** Read every app variant, callback and model module; find which file is
   actually deployed (the `Procfile` may not say: check the host's settings); list
   dependencies, filesystem use and threads. Record working-tree state.
2. **Baseline.** Create the feature branch. Set up a Python 3.12 environment outside the
   repo. Run the original under its real start command. Record reference outputs and
   timings. Run a seeded random sweep **and the special-value grid** through the
   original, recording crashes, empty results and data the fixed axes hide.
3. **Assess Vercel.** Fetch the current Python runtime and limits documentation. Measure
   the Linux bundle size, worst-case callback time and the largest response. Choose an
   approach (2.1).
4. **Build:** `params.py` → `physics.py` (bit-for-bit test immediately) → `figures.py`
   → `layout.py` + `style.css` → `callbacks.py` → `app.py`. Review screenshots at
   1600, 1280, 820 and 390 px after each visual change.
5. **Deploy config:** `requirements.txt`, `requirements-dev.txt`, `vercel.json`,
   `.gitignore` (`.vercel/`).
6. **Tests:** section 12, then run in a fresh environment built only from
   `requirements-dev.txt`.
7. **Review:** diff against `main` should show only additions. Old app and model files
   identical; secret scan; lint (`pyflakes`); unused CSS tokens; `!important` only for
   reduced motion.
8. **Deploy and verify:** 11.6–11.7, then 12.6. Document in `docs/` what was verified
   and what wasn't.

---

## 14. Adapting to another app

### 14.1 What carries over unchanged [Essential]

The layering (2.3); the parameter registry (4); the model adapter with
`np.errstate` (5.2); the display rules (6.1); the page skeleton and responsive rules
(7); the main callback shape, resets, export and Enter-to-commit (8); the tokens, Dash
theming and fixed plot frames (9); everything in 11; the test suite structure (12).

### 14.2 What you replace [Project-specific]

| Item | Here |
|---|---|
| Model modules and the adapter's assembly code | `Generic_precession`, redshifted chirp mass, Planck 2018 ΛCDM |
| Parameters, groups, defaults, presets | 14 sliders; edge-on and face-on presets |
| Panels, default axis windows, ticks, caps | Modulation, amplitude (±20 cap), phase |
| Derived quantities | Total mass, redshifted mass, f_cut, cos ι_JN, d_L |
| Notes and messages | B_NP ≈ 0, stationary L, curves equal to the reference, out of band |
| Copy, branding, palette, fonts, logo | "Generic Precession", Inter, `bbh_spin.jpg`, "© Tamanjyot Singh" |
| Hosting names | Project `generic-precession`, `genericprecession.vercel.app`, Render on `main` |

### 14.3 Feasibility checks for a new app [Essential]

- **Dependencies:** `lal`, `PyCBC` or similar large binary packages may not fit or have
  wheels. Check the Linux bundle first (11.3).
- **Compute time:** matches, parameter sweeps or fits can exceed 300 s. Time the worst
  case; consider approach C.
- **State:** files written per request, uploads or caches don't persist. Redesign or add
  storage, and get approval for any paid service.

### 14.4 CLAUDE.md template for the new repository

```markdown
# CLAUDE.md: <App name>: UI redesign and Vercel deployment

Follow docs/implementation-guide.md (copied from generic_precession_app).

## Goal
Redesign <live app file> into a polished interactive app on a new branch
`feature/vercel-ui`, deployed to Vercel Hobby alongside <existing host> on `main`.

## Hard rules
- Scientific correctness first: do not modify <model modules> for the UI; every plotted
  array must equal the original app bit-for-bit (same platform) in tests.
- Only add files on the branch; never touch <old app path>, its requirements or its
  start configuration. Model changes go to `main` first, then merge into the branch.
- Never hide or collapse a figure; explain empty, error and surprising states in the UI.
- Don't push, merge, deploy to production or add paid services without approval.

## Project facts
- Live app file: <path>; host: <host>; start command: <cmd>
- Model entry points: <functions/classes>; special substitutions: <list>
- Parameters: <id, range, step, default, unit for each>
- Panels and their original axis windows: <list>
- Branding: <title, copyright, logo, palette>

## Definition of done
Same as guide sections 12–13: regression, params, figures, HTTP and UI tests pass;
deployed; live sweep, numeric comparison and browser checks pass; docs/ updated with
what was and wasn't verified.
```

---

## 15. Pitfalls

Each of these happened in this project.

| Pitfall | Symptom | Fix |
|---|---|---|
| `np.seterr` is per-thread | Blank panels with no error, only on the server | `np.errstate(**MODEL_ERRSTATE)` around every model call (5.2) |
| Hiding a Plotly container (`display:none`) | Plots shrink to 8 px and sit hidden behind their cards after leaving an error state | Fixed-height frames; placeholders instead of hiding (9.4) |
| Fixed axis windows | Data off-screen for most inputs | Default window widened to fit (6.1) |
| Reference line drawn over curves | Curves equal to the reference look missing | `layer="below"`; notes (6.1, 8.1) |
| Divergent ratio | Axis stretched to 10¹⁵, panel looks flat | Cap the axis and label it (6.1) |
| Special values on the slider grid | Exact singular geometry (0/0) reached by ordinary drags | Sweep the special-value grid; handle in the model or explain in the UI |
| Dash 4 number box | Typed value ignored on Enter | Blur on Enter (8.4) |
| `run_server` | `python app.py` crashes on Dash ≥ 3 | `app.run(...)` |
| Framework preset "Other" | Build error about `api/` functions | `"framework": "flask"` in `vercel.json` |
| Bundle > 225 MB | Packages installed at every cold start (6.9 s) | 11.4 variables (3.2 s) |
| `vercel link` writes `.env.local` | OIDC token could be uploaded with the source | Delete before deploy; deploy from `git archive` |
| Untracked local files | Personal files uploaded | Deploy from `git archive` |
| Old Vercel project connected to the repo | Failed builds on every push | Delete or disconnect it |
| Unpinned requirements on the old host | Rebuild silently upgrades to new major versions | Pin on the new branch; test the old app with current versions before pushing `main` |
| Cross-platform numerics | Live ODE results differ from local at about 10⁻⁷ | Bit-for-bit only on the same platform; use tolerances live |
| Hobby log retention | `vercel logs --since 1h` errors | Omit `--since` |
| Element screenshots with a sticky header | The header overlaps the captured region | Use full-page screenshots for review |
