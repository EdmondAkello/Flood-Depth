#!/bin/bash
# FloodDepth full local check: flake8, Bandit, the pure suites (with and without SciPy) and the
# QGIS smoke test on QGIS 3.44 and 4.2 (Docker images qgis/qgis, as in CI).
#   bash tools/dev/run_tests.sh                 (PY=python3 by default)
#   QGIS_IMAGES="4.2.3-trixie" bash tools/dev/run_tests.sh
# Without Docker, the smoke test runs with $PY when it can import qgis (e.g. the OSGeo4W shell).
set -u
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
PY=${PY:-python3}
status=0
cd "$REPO"
echo "== flake8";  $PY -m flake8 . && echo "flake8: ok" || status=1
echo "== bandit";  $PY -m bandit -q -r . -x ./tests,./tools,./docs,./benchmarks && echo "bandit: ok" || status=1
echo "== pure suites (SciPy if installed)"; PY=$PY bash tools/dev/ci_pure.sh || status=1
echo "== pure suites (fallbacks)"; FLOODDEPTH_NO_SCIPY=1 PY=$PY bash tools/dev/ci_pure.sh | tail -1 || status=1

WORK="$(mktemp -d)"
mkdir -p "$WORK/flooddepth"
tar --exclude=.git --exclude=dist -cf - . | tar -x -C "$WORK/flooddepth"
if command -v docker >/dev/null 2>&1; then
  for img in ${QGIS_IMAGES:-3.44.15-trixie 4.2.3-trixie}; do
    docker run --rm -v "$WORK:/w" -w /w -e QT_QPA_PLATFORM=offscreen \
      -e PYTHONPATH=/usr/share/qgis/python:/usr/share/qgis/python/plugins "qgis/qgis:$img" \
      python3 -m flooddepth.tests.qgis_smoke > "$WORK/smoke_$img.log" 2>&1 || status=1
    echo "QGIS $img smoke: pass=$(grep -c '\[PASS\]' "$WORK/smoke_$img.log") fail=$(grep -c '\[FAIL\]' "$WORK/smoke_$img.log")"
    grep '\[FAIL\]' "$WORK/smoke_$img.log"
  done
elif $PY -c "import qgis.core" 2>/dev/null; then
  (cd "$WORK" && QT_QPA_PLATFORM=offscreen $PY -m flooddepth.tests.qgis_smoke > smoke.log 2>&1) || status=1
  echo "QGIS smoke: pass=$(grep -c '\[PASS\]' "$WORK/smoke.log") fail=$(grep -c '\[FAIL\]' "$WORK/smoke.log")"
else
  echo "QGIS smoke: skipped (no Docker and no QGIS Python)"
fi
rm -rf "$WORK"
exit $status
