#!/bin/bash
# Run the pure (QGIS-free) test suites; fail if any suite errors or reports a [FAIL].
# Used by .github/workflows/tests.yml; also runs locally: bash tools/dev/ci_pure.sh
# Tests import `flooddepth.*`, so the repository is linked as a package named
# `flooddepth` in a temporary folder (the clone may be called Flood-Depth).
set -u
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
WORK="$(mktemp -d)"
ln -s "$REPO" "$WORK/flooddepth"
cd "$WORK"
PY=${PY:-python3}
status=0
total=0
for f in flooddepth/tests/test_*.py; do
  m=$(basename "$f" .py)
  out=$($PY -m "flooddepth.tests.$m" 2>&1); rc=$?
  p=$(printf '%s\n' "$out" | grep -c '\[PASS\]')
  x=$(printf '%s\n' "$out" | grep -c '\[FAIL\]')
  total=$((total + p))
  echo "$m: pass=$p fail=$x rc=$rc"
  if [ "$rc" -ne 0 ] || [ "$x" -ne 0 ]; then
    status=1
    printf '%s\n' "$out" | grep -E '\[FAIL\]|Error|Traceback' | head -20
  fi
done
echo "total checks passed: $total"
rm -rf "$WORK"
exit $status
