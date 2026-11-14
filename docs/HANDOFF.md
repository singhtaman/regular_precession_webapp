# Handoff — 14 Nov 2026

Record of the work done with Claude Code from 26 Oct to 14 Nov 2026 on the personal
website and the Regular Precession app, what state each is in, and what is left.

## 1. Regular Precession app (this repository)

### Branches

| Branch | Purpose | Deploys | State |
|---|---|---|---|
| `main` | Original Render app (`RP_map_all_app.py`) | Render | untouched |
| `feature/vercel-ui` | Redesigned Dash app (`app.py` → `webapp/`) | **Vercel production on every push** | untouched |
| `feature/github-pages` | Browser version (`pages/`) + time-domain view | GitHub Pages, once set up | **3 commits, local only, not pushed** |

`feature/github-pages` is based on `feature/vercel-ui`:

- `660cad1` Add a browser version for GitHub Pages
- `ca03426` Rename the Pyodide bridge module (`bridge.py` → `rp_page.py`)
- `64fbcbe` Add a frequency/time switch with a time-domain strain panel

None of the existing files (the physics, `webapp/*.py` except the new `time_domain.py`, the
Dash app, `vercel.json`, the tests) were changed; `.gitignore` gained `_site/`.

### What the browser version is

Static page in `pages/`, assembled into `_site/` by `pages/build.sh` and published by
`.github/workflows/pages.yml`. The unmodified model runs in Pyodide 0.29.0 in a Web Worker;
`pages/figures.js` is a JavaScript port of `webapp/figures.py` (Plotly is not available in
Pyodide), so plot changes must be made in both. Full description: `docs/github-pages.md`.

First load ≈32 MB and ≈8 s; each update 20–45 ms.

### Verified

- **Agreement with the live Vercel app**, CSV export of 26 inputs (defaults, presets, slider
  extremes, 20 seeded random): f, φ_L, cos θ_L, amplitude within 7×10⁻¹³ relative; phase within
  5×10⁻⁶ rad worst case (θ̃ = 10, Ω̃ = 4), ~10⁻¹³ typically. The phase difference comes from
  `odeint` in SciPy 1.14.1 (Pyodide) vs 1.18.1 (Vercel pins). No-signal inputs agree.
- **Interface** in headless Chromium and WebKit: presets, resets, typed and out-of-range values,
  keyboard slider control, no-signal state, CSV download in both domains, frequency↔time switch,
  phone width (390 px) without horizontal scroll, no console errors.
- **Time domain** in Pyodide: chirp runs forward (≈20 Hz → near f_cut at t ≈ 0); durations
  3.5 s (defaults) to 12.7 s (longest slider setting) fit the 20 s window; forward FFT recovers
  the tapered h̃(f) to 0.1–0.4% (loss = trimmed pre-signal stretch).

### Not verified

- The repository's pytest suite (needs Python 3.12; not installed on this Mac). No file it tests
  was changed.
- Real phones/tablets and browsers with extensions.
- The GitHub Actions workflow has never run.

### Correction to record

The "Failed to fetch" seen on 14 Nov was first blamed on a content blocker and `bridge.py` was
renamed. The real cause was Python's `http.server` refusing parallel connections (backlog of 5)
during local preview; the worker now fetches files sequentially with retries. The rename is
harmless and was kept.

### Time-domain choices made without explicit sign-off (review these)

`webapp/time_domain.py`: Tukey taper α = 0.1; t_c = 0 at f_cut; 0.25 s shown after t_c and
0.25 s before the signal; sampling rate = smallest power of two with ≥ 6 samples/cycle at f_cut;
the panel opens on the last 0.5 s. These are presentation choices; the model is unchanged.

### Next steps

1. Review locally: `pages/build.sh && python3 -m http.server -d _site 8000`.
2. Push `feature/github-pages` (Vercel ignores non-production branches).
3. GitHub → Settings → Pages → Source: **GitHub Actions**; Settings → Environments →
   github-pages → allow `feature/github-pages`. Site: https://singhtaman.github.io/regular_precession_webapp/
4. Optional: link it from `Precession_GWs.html` on the website next to the Vercel link.
5. Optional: bring the time view to the Vercel app — reuse `webapp/time_domain.py`, add a
   strain figure to `webapp/figures.py`, a switch in `layout.py`/`callbacks.py`, and tests.
   Pushing `feature/vercel-ui` deploys to production.
6. Optional: same browser build for Generic Precession (local `../generic_precession_app`; no
   public repo was found for it).

## 2. Personal website (`singhtaman.github.io`)

All on `main`, pushed and live; `main` is the only branch.

- Pages was serving `site-redesign`, so commits on `site-classic` weren't live. Merged both into
  `main` (`6c36528`, `5f485e8`, `24132be`), removed `/classic/`, Pages switched to `main`,
  deleted `site-classic` and `site-redesign`. PR #1 (from `site-classic`) is contained in `main`.
- 5 Nov: the whole site returned 404 because the repository had become private; it was made
  public again.
- `63513a7` light/dark toggle (`assets/js/theme.js`, header button, saved per browser).
- `a90d98c` poetry page moved to `assets/hidden/kabit.html` (noindex; nothing links to it).
- `7f65c1e` poetry page restyled in the redesign; reading-shelf covers enlarged and linked;
  Omar El Akkad cover replaced (Barnes & Noble blocks hotlinking).
- `CLAUDE.md` in the website repo is still untracked.
