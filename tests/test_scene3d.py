# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""3D scene: meshes, packing, mud and debris, the HTML page (no browser needed).

    python -m flooddepth.tests.test_scene3d
"""

import json
import os
import re
import tempfile

import numpy as np

from .common import check, finish
from ..core import raster as R
from ..core import scene3d as S3


def main():
    d = tempfile.mkdtemp()
    print("\n1. Packing is lossless")
    a = np.random.default_rng(1).normal(size=(500, 3)).astype("float32")
    check("float positions", np.array_equal(S3.unpack(S3._pack(a, "<f4")), a))
    idx = np.random.default_rng(2).integers(0, 10_000, 3000).astype("uint32")
    check("delta-coded indices", np.array_equal(S3.unpack(S3._pack(idx, "<u4", delta=True)), idx))

    print("\n2. Grid mesh")
    z = np.arange(12, dtype=float).reshape(3, 4)
    keep = np.ones_like(z, bool)
    keep[0, 0] = False
    xyz, faces, _ = S3.grid_mesh(z, keep, 0.0, 30.0, 10.0)
    check("11 vertices, 2 triangles per complete square", len(xyz) == 11 and len(faces) == 2 * 5, f"{len(faces)}")
    check("cell centres in map units", np.allclose(xyz[0, :2], [15.0, 25.0]))

    print("\n3. Scene from arrays: water, mud and debris, lines")
    rows, cols = 120, 160
    yy, xx = np.mgrid[0:rows, 0:cols].astype("float32")
    dem = 1500.0 + 0.02 * yy + 0.002 * (xx - 80) ** 2
    flood = np.abs(xx - 80) < 12
    wse = np.where(flood, dem.min(axis=1, keepdims=True) + 1.5, np.nan)
    depth = np.where(flood, np.maximum(wse - dem, 0), np.nan)
    evid = np.where(flood, 1, 255).astype("uint8")
    evid[:40][flood[:40]] = 3                       # the upper third is mud and debris
    g = R.Grid(cols, rows, (500000.0, 30.0, 0.0, 9900000.0, 0.0, -30.0), R.srs_from_epsg(32737).ExportToWkt())
    road = ["LINESTRING(500300 9898000, 504500 9898200)"]
    info = S3.write_scene(os.path.join(d, "s.html"), g, dem, flood, wse, depth, evidence=evid,
                          lines=[("Road", "Road", road, g.wkt, "#555555")], max_vertices=5000,
                          style={"mud": "#7a5230"})
    page = open(info["path"], encoding="utf-8").read()
    data = json.loads(re.search(r'<script id="scene-data" type="application/json">(.*?)</script>', page, re.S).group(1))
    ids = [s["id"] for s in data["surfaces"]]
    check("terrain, water and mud surfaces", ids == ["terrain", "water", "mud"], str(ids))
    check("level of detail respected", len(S3.unpack(data["meshes"]["terrain"]["positions"])) <= 5000,
          f"step {info['step']}")
    mud_rgb = S3.unpack(data["meshes"]["mud"]["colors"])
    k = mud_rgb[:, 0].astype(float) / 0x7a
    check("mud in the user's colour (shaded by elevation)",
          np.allclose(mud_rgb[:, 2], k * 0x30, atol=2) and k.max() <= 1.01, str(mud_rgb[0]))
    wd = S3.unpack(data["meshes"]["water"]["depth"])
    check("water depths carried for the click read-out", np.nanmax(wd) > 0.5 and np.nanmin(wd) >= 0)
    check("road drawn on the terrain", any(ln["name"] == "Road" for ln in data["lines"]))
    check("depth legend on the water", len(data["surfaces"][1]["style"]["classes"]) == 5)
    check("no network in the page", "<script src" not in page and "__SCENE__" not in page)
    check("vertical exaggeration chosen", 1.0 <= info["vertical_exaggeration"] <= 8.0)
    try:
        S3.write_scene(os.path.join(d, "e.html"), g, dem, np.zeros_like(flood), wse, depth)
        check("no flood is refused", False)
    except ValueError:
        check("no flood is refused", True)


if __name__ == "__main__":
    main()
    finish("test_scene3d")
