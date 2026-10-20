#!/bin/bash
# Build the installable plugin zip (folder flooddepth/ inside) from the committed tree.
# Files marked export-ignore in .gitattributes (tests, tools, docs, .github) are left out.
#   bash tools/dev/make_zip.sh            -> dist/FloodDepth_<version>.zip
set -eu
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
VER=$(grep '^version=' metadata.txt | cut -d= -f2)
mkdir -p dist
git archive --format=zip --prefix=flooddepth/ -o "dist/FloodDepth_${VER}.zip" HEAD
echo "dist/FloodDepth_${VER}.zip"
