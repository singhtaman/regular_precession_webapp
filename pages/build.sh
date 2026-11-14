#!/bin/sh
# Assemble the static browser build in _site/ (used by .github/workflows/pages.yml).
# The Python model files are copied unchanged from the repository.
set -eu
cd "$(dirname "$0")/.."

rm -rf _site
mkdir -p _site/py/webapp
cp pages/index.html pages/app.js pages/figures.js pages/worker.js pages/pages.css _site/
cp webapp/assets/style.css webapp/assets/bbh_spin.jpg _site/
cp regular_precession_no_match.py systems_lib.py pages/rp_page.py _site/py/
cp webapp/__init__.py webapp/params.py webapp/physics.py webapp/time_domain.py webapp/sphere.py _site/py/webapp/
echo "Built _site/ ($(find _site -type f | wc -l | tr -d ' ') files)"
