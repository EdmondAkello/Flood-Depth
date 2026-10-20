#!/bin/bash
# FloodDepth full test run: pure suites + the QGIS smoke test.
#   bash tools/dev/run_tests.sh            (PY=python3.12 by default; needs QGIS's Python for the smoke test)
set -u
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
PY=${PY:-python3.12}
PY=$PY bash "$REPO/tools/dev/ci_pure.sh"
WORK="$(mktemp -d)"; ln -s "$REPO" "$WORK/flooddepth"; cd "$WORK"
QT_QPA_PLATFORM=offscreen $PY -m flooddepth.tests.qgis_smoke > "$REPO/smoke.log" 2>&1
echo "QGIS smoke: pass=$(grep -c '\[PASS\]' "$REPO/smoke.log") fail=$(grep -c '\[FAIL\]' "$REPO/smoke.log") done=$(grep -c 'ALL QGIS SMOKE CHECKS PASSED' "$REPO/smoke.log")"
rm -rf "$WORK"
if command -v bandit >/dev/null; then
  b=$(bandit -q -r "$REPO" -x "$REPO/tools,$REPO/tests,$REPO/docs" -f csv 2>/dev/null | tail -n +2 | grep -c .)
  echo "bandit (shipped files): issues=$b"
fi
