# Regular Precession explorer: redesigned UI on Vercel

This branch (`feature/vercel-ui`) adds a redesigned version of the interactive
regular-precession app and the configuration to host it on Vercel's Hobby plan.
The Render app (`RP_map_all_app.py`), the `Procfile` and the scientific modules
(`regular_precession_no_match.py`, `systems_lib.py`) are byte-identical to `main`.

| | Render (original) | Vercel (this branch) |
|---|---|---|
| Code | `RP_map_all_app.py` | `app.py` → `webapp/` |
| Branch | `main` | `feature/vercel-ui` |
| Dependencies | `requirements.txt` on `main` (unpinned) | `requirements.txt` on this branch (pinned) |
| Server | `gunicorn RP_map_all_app:server` | Vercel Function (Python 3.12) |
| URL | your Render service | https://regularprecession.vercel.app |

**Do not merge this branch into `main` as is.** Render builds from `main` and installs
its `requirements.txt`; on this branch that file is the pinned Vercel runtime set, which
lacks `gunicorn` and `matplotlib` that the Render app needs.

## What changed

- Parameters grouped into Precession, Source and Orientation, with symbols, units and
  hints, per-group and global reset, and edge-on / face-on presets (systems 2 and 1 in
  `systems_lib.py`). All slider ids, ranges, steps, marks and defaults are unchanged.
- Results in a 2×2 grid (modulation of L | system properties, amplitude | phase), with a
  status line, CSV export of every plotted array at full precision, and PNG export.
- Explicit messages for invalid input, no frequencies in band (f_cut < 20 Hz), model
  errors, a divergent amplitude ratio (B_NP ≈ 0) and a stationary L.
- Axes open on the original windows and widen only when data would be off-screen (the
  fixed windows hid part of the data in 2,343 of 2,409 sampled inputs). The amplitude
  axis is capped at ±20 with a label when the ratio diverges.

`webapp/physics.py` reproduces the original parameter assembly line for line and calls
the unmodified `Regular_precession` class. Each model (RP, NP) is built once per request
instead of once per panel, and evaluation runs under the model's NumPy error state
(`np.seterr` is per-thread). Neither changes any value.

The phase panel keeps the original definition: it plots −(φ_P + 2 δφ) of the RP model.
The NP polarization phase is constant in f and is not subtracted, as in the original.

## Run and test locally

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements-dev.txt
python app.py                     # http://127.0.0.1:8050
pytest                            # 913 tests; browser tests: pip install playwright, then pytest -m ui
```

`test_physics_regression.py` compares every plotted array with `RP_map_all_app.py`
bit-for-bit over named cases, slider extremes, 150 seeded random inputs and the
special-angle grid.

## Deploy

The GitHub repository is connected to the Vercel project `regularprecession`:

- Pushing to `feature/vercel-ui` (the production branch) deploys to production.
- Pushes to any other branch, including `main`, are skipped by the Ignored Build Step
  `if [ "$VERCEL_GIT_COMMIT_REF" = "feature/vercel-ui" ]; then exit 1; else exit 0; fi`.
- `VERCEL_SUPPORT_LARGE_FUNCTIONS=1` and `VERCEL_PYTHON_DISABLE_RUNTIME_INSTALL=1` are set
  on the project because the dependency bundle is above 225 MB.

To deploy by hand instead, follow `VERCEL_APP_GUIDE.md` section 11.6
(`git archive feature/vercel-ui`, then `vercel deploy --prod`).
