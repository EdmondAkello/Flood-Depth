# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""The whole quick workflow on the synthetic site (no QGIS).

    python -m flooddepth.tests.test_pipeline
"""

import json
import os
import tempfile

import numpy as np

from .common import check, finish
from . import synth as S
from ..core.pipeline import QuickOptions, run_quick, permanent_from_raster
from ..core import raster as R


def depth_error(path):
    a, g, _ = R.read_band(path)
    x = g.gt[0] + (np.arange(g.cols) + 0.5) * g.gt[1]
    y = g.gt[3] + (np.arange(g.rows) + 0.5) * g.gt[5]
    X, Y = np.meshgrid(x, y)
    z, wse, lake = S.truth(X, Y)
    t = np.where(~lake & (wse > z), wse - z, np.nan)
    both = np.isfinite(a) & np.isfinite(t)
    return (a - t)[both], both.sum(), a, lake


def main():
    d = tempfile.mkdtemp(prefix="fd_pipe_")
    dem = S.write_dem(os.path.join(d, "dem_fabdem.tif"))
    post = S.write_s2_raw(os.path.join(d, "post"), "2024-05-05", True)
    pre = S.write_s2_raw(os.path.join(d, "pre"), "2024-04-07", False, seed=2)
    q = lambda **k: run_quick(QuickOptions(dem_path=dem, **k), log=None)

    print("\n1. Image -> depth, lake removed with the pre-event image")
    info = q(image_paths=post, pre_paths=pre, out_dir=os.path.join(d, "a"), uncertainty_members=8)
    err, n, a, lake = depth_error(info["outputs"]["depth"])
    check("depth bias |b| < 0.1 m", abs(err.mean()) < 0.1, f"bias {err.mean():+.3f} m over {n} cells")
    check("depth MAE < 0.15 m", np.abs(err).mean() < 0.15, f"MAE {np.abs(err).mean():.3f} m")
    check("no depth on the lake", not np.isfinite(a[lake]).any())
    check("FABDEM preset picked from the file name", info["uncertainty"]["sigma_dem_m"] == 1.5)
    with open(os.path.join(d, "a", "run.json")) as f:
        js = json.load(f)
    check("run.json is valid JSON with the threshold", "threshold" in js["classification"])
    check("evidence: observed water", any(r["code"] == 1 for r in info["evidence"]))

    print("\n1b. Pre-event image without SWIR1 (Browser 'False color' case): NDWI fallback")
    pre_noswir = [p for p in pre if "B11" not in os.path.basename(p) and "B12" not in os.path.basename(p)]
    info = q(image_paths=post, pre_paths=pre_noswir, out_dir=os.path.join(d, "a2"), uncertainty_members=0)
    ci = info["classification"]
    check("pre-event used for permanent water", ci.get("permanent_share", 0) > 0, str(ci.get("pre_bands")))
    err, n, a, lake = depth_error(info["outputs"]["depth"])
    check("no depth on the lake (NDWI pre-event)", not np.isfinite(a[lake]).any())

    print("\n2. Without the pre-event image the lake is not flood")
    info = q(image_paths=post, out_dir=os.path.join(d, "b"), uncertainty_members=0)
    err, n, a, lake = depth_error(info["outputs"]["depth"])
    check("lake edges are waterline only where dry ground is next to them; flood MAE < 0.2 m",
          np.abs(err).mean() < 0.2, f"MAE {np.abs(err).mean():.3f}")

    print("\n3. Polygon extent with a permanent-water polygon")
    info = q(extent_wkts=[S.flood_polygon_wkt()], extent_crs=S.utm_wkt(), out_dir=os.path.join(d, "c"),
             uncertainty_members=0)
    err, n, _, _ = depth_error(info["outputs"]["depth"])
    check("polygon: MAE < 0.1 m", np.abs(err).mean() < 0.1, f"MAE {np.abs(err).mean():.3f}")

    print("\n4. Flood trace (vegetation stripped, water gone)")
    tpre = S.write_s2_trace(os.path.join(d, "tpre"), "2024-04-07", False)
    tpost = S.write_s2_trace(os.path.join(d, "tpost"), "2024-05-12", True, seed=4)
    try:
        q(image_paths=tpost, out_dir=os.path.join(d, "d0"), uncertainty_members=0)
        check("water mode on a drained scene is refused", False)
    except ValueError as exc:
        check("water mode on a drained scene is refused", "flood-trace" in str(exc), str(exc)[:70])
    info = q(image_paths=tpost, pre_paths=tpre, mode="trace", out_dir=os.path.join(d, "d"), uncertainty_members=6,
             seed_wkts=[S.centreline_wkt()], seed_crs=S.utm_wkt())
    err, n, _, _ = depth_error(info["outputs"]["depth"])
    check("trace recovers the flood strip", n > 2000, f"{n} cells")
    check("trace depth MAE < 0.2 m", np.abs(err).mean() < 0.2, f"MAE {np.abs(err).mean():.3f}")
    codes = {r["code"]: r for r in info["evidence"]}
    check("bare sediment classed as erosion / deposition", 3 in codes and codes[3]["cells"] > 0.8 * n)
    check("erosion cells carry extra terrain error", info["uncertainty"]["sigma_erosion_m"] > 0)

    print("\n5. Permanent water from a raster")
    wc = np.array([[10, 80, 30], [80, 0, 40]], float)
    m, rule = permanent_from_raster(wc)
    check("WorldCover -> class 80", m.tolist() == [[False, True, False], [True, False, False]] and "WorldCover" in rule)
    occ = np.linspace(0, 100, 30).reshape(5, 6)
    m, rule = permanent_from_raster(occ)
    check("JRC occurrence -> >= 50 %", m.sum() == (occ >= 50).sum() and "occurrence" in rule)
    m, rule = permanent_from_raster(np.array([[1.0, 2.0, 5.0]]), "2,5")
    check("explicit class list", m.tolist() == [[False, True, True]])


if __name__ == "__main__":
    main()
    finish("test_pipeline")
