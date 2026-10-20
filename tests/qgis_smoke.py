# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""QGIS-level smoke test: every FloodDepth Processing tool, run inside QGIS.

Needs a QGIS installation. Builds the synthetic site (tests/synth.py) and
runs the real Processing chain, checking that every tool completes, its
outputs open, and the depths match the known answer.

Run from the folder that CONTAINS the `flooddepth` plugin folder:

    QT_QPA_PLATFORM=offscreen python -m flooddepth.tests.qgis_smoke

or paste into the QGIS Python console:

    from flooddepth.tests import qgis_smoke; qgis_smoke.main(in_qgis=True)
"""

import os
import sys
import tempfile

import numpy as np

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def _start_qgis():
    from qgis.core import QgsApplication
    app = QgsApplication([], False)
    app.initQgis()
    for p in ("/usr/share/qgis/python/plugins",
              os.path.join(QgsApplication.prefixPath(), "python", "plugins")):
        if os.path.isdir(p) and p not in sys.path:
            sys.path.append(p)
    from processing.core.Processing import Processing
    Processing.initialize()
    return app


def _read(path):
    from osgeo import gdal
    ds = gdal.Open(path)
    b = ds.GetRasterBand(1)
    a = b.ReadAsArray().astype("float64")
    nd = b.GetNoDataValue()
    if nd is not None:
        a[a == nd] = np.nan
    return a, ds.GetGeoTransform(), ds.GetProjection()


def _true_depth_on(path):
    from . import synth as S
    a, gt, _ = _read(path)
    rows, cols = a.shape
    x = gt[0] + (np.arange(cols) + 0.5) * gt[1]
    y = gt[3] + (np.arange(rows) + 0.5) * gt[5]
    X, Y = np.meshgrid(x, y)
    z, wse, lake = S.truth(X, Y)
    return a, np.where(~lake & (wse > z), wse - z, np.nan)


def main(in_qgis=False):
    app = None if in_qgis else _start_qgis()
    import processing
    from qgis.core import QgsApplication
    from ..processing_provider.provider import FloodDepthProvider
    from . import synth as S
    prov = FloodDepthProvider()
    QgsApplication.processingRegistry().addProvider(prov)
    ids = sorted(a.id() for a in prov.algorithms())
    print("\n1. Provider")
    check("7 algorithms registered", len(ids) == 7, ", ".join(ids))
    for old in ("flooddepth:classify_inundation", "flooddepth:compute_flood_depth",
                "flooddepth:compute_constant_stage_depth", "flooddepth:estimate_stage_dem_boundary"):
        check(f"0.3 id kept: {old}", old in ids)

    tmp = tempfile.mkdtemp(prefix="flooddepth_smoke_")
    dem = S.write_dem(os.path.join(tmp, "dem_tandemx.tif"))
    post = S.write_s2_raw(os.path.join(tmp, "post"), "2024-05-05", flooded=True)
    pre = S.write_s2_raw(os.path.join(tmp, "pre"), "2024-04-07", flooded=False, seed=2)
    ext = S.write_vector(os.path.join(tmp, "extent.gpkg"), [S.flood_polygon_wkt()])
    line = S.write_vector(os.path.join(tmp, "river.gpkg"), [S.centreline_wkt()], "line")

    print("\n2. Quick flood depth from the image files (pre-event removes the lake)")
    out = os.path.join(tmp, "quick")
    r = processing.run("flooddepth:quickflooddepth", {
        "IMAGES": post, "PRE_IMAGES": pre, "DEM": dem, "MODE": 0, "METHOD": 0, "MEMBERS": 6,
        "OUT_DIR": out})
    for f in ("water_mask.tif", "depth.tif", "wse.tif", "evidence.tif", "flood_extent.gpkg", "report.html",
              "run.json", "depth_p10.tif", "depth_p90.tif", "flood_probability.tif"):
        check(f"writes {f}", os.path.exists(os.path.join(out, f)))
    d, t = _true_depth_on(r["DEPTH"])
    both = np.isfinite(d) & np.isfinite(t)
    err = (d - t)[both]
    check("flooded cells found", both.sum() > 1000, f"{both.sum()} cells")
    check("depth bias within 0.15 m", abs(err.mean()) < 0.15, f"bias {err.mean():+.3f} m")
    check("depth MAE within 0.25 m", np.abs(err).mean() < 0.25, f"MAE {np.abs(err).mean():.3f} m")
    lake_x, lake_y = S.LAKE[0], S.LAKE[1]
    a, gt, _ = _read(r["DEPTH"])
    rr, cc = int((lake_y - gt[3]) / gt[5]), int((lake_x - gt[0]) / gt[1])
    check("lake (permanent water) has no flood depth", not np.isfinite(a[rr, cc]))
    p10, _, _ = _read(os.path.join(out, "depth_p10.tif"))
    p90, _, _ = _read(os.path.join(out, "depth_p90.tif"))
    ok = np.isfinite(p10) & np.isfinite(p90)
    check("P10 <= P90 everywhere", bool(np.all(p10[ok] <= p90[ok] + 1e-5)))
    with open(r["REPORT"], encoding="utf-8") as fh:
        check("report mentions the uncertainty", "Uncertainty" in fh.read())

    print("\n3. Quick flood depth from a polygon extent")
    r2 = processing.run("flooddepth:quickflooddepth", {"EXTENT": ext, "DEM": dem, "MEMBERS": 0,
                                                       "OUT_DIR": os.path.join(tmp, "quick_poly")})
    d, t = _true_depth_on(r2["DEPTH"])
    both = np.isfinite(d) & np.isfinite(t)
    err = (d - t)[both]
    check("polygon extent: MAE within 0.15 m", np.abs(err).mean() < 0.15, f"MAE {np.abs(err).mean():.3f} m")

    print("\n4. Step 1: classify")
    r3 = processing.run("flooddepth:classify_inundation", {
        "IMAGES": post, "MASK": os.path.join(tmp, "mask.tif"), "POLYGONS": os.path.join(tmp, "water.gpkg")})
    m, _, _ = _read(r3["MASK"])
    check("mask has water and dry", (m == 1).sum() > 0 and (m == 0).sum() > 0)
    check("alpha strip is NoData", np.isnan(m[:, -1]).all())

    print("\n5. Step 2: surface and depth from the extent (with waterline points)")
    r4 = processing.run("flooddepth:estimate_stage_dem_boundary", {
        "EXTENT": ext, "DEM": dem, "WSE": os.path.join(tmp, "wse2.tif"), "DEPTH": os.path.join(tmp, "depth2.tif"),
        "WATERLINE": os.path.join(tmp, "waterline.gpkg")})
    d, t = _true_depth_on(r4["DEPTH"])
    both = np.isfinite(d) & np.isfinite(t)
    check("step 2 MAE within 0.15 m", np.abs(d - t)[both].mean() < 0.15)
    check("waterline points written", os.path.exists(r4["WATERLINE"]))
    r4b = processing.run("flooddepth:estimate_stage_dem_boundary", {
        "MASK": r3["MASK"], "DEM": dem, "WSE": os.path.join(tmp, "wse3.tif"), "DEPTH": os.path.join(tmp, "depth3.tif")})
    check("step 2 from a mask raster", os.path.exists(r4b["DEPTH"]))

    print("\n6. Step 3: depth from a water-surface raster, and below a level")
    r5 = processing.run("flooddepth:compute_flood_depth", {"DEM": dem, "STAGE": r4["WSE"],
                                                           "OUTPUT": os.path.join(tmp, "depth4.tif")})
    a5, _, _ = _read(r5["OUTPUT"])
    check("depth from WSE raster", np.isfinite(a5).sum() > 1000)
    r6 = processing.run("flooddepth:compute_constant_stage_depth", {"DEM": dem, "LEVEL": 1199.0,
                                                                    "OUTPUT": os.path.join(tmp, "depth5.tif")})
    a6, _, _ = _read(r6["OUTPUT"])
    check("constant level: max depth = level - min DEM", abs(np.nanmax(a6) - (1199.0 - (S.Z0 - S.S_LONG * S.ROWS * S.CS
                                                                                       + S.S_CROSS * S.CS / 2))) < 0.2,
          f"{np.nanmax(a6):.2f} m")

    print("\n7. Compare extents")
    r7 = processing.run("flooddepth:compareextents", {"PRED_RASTER": r3["MASK"], "REF_POLYS": ext, "CELL": 10,
                                                      "AGREEMENT": os.path.join(tmp, "agree.tif")})
    check("CSI of the classification vs the true extent > 0.8", r7["CSI"] > 0.8, f"CSI {r7['CSI']:.3f}")

    print("\n8. Supervised (training polygons) on the green, NIR, SWIR files")
    green = [p for p in post if "_B03_" in p][0]
    w = S.write_vector(os.path.join(tmp, "wt.gpkg"), [
        f"POLYGON(({S.XC - 100} {S.Y1 - 3000}, {S.XC + 100} {S.Y1 - 3000}, {S.XC + 100} {S.Y1 - 3600}, "
        f"{S.XC - 100} {S.Y1 - 3600}, {S.XC - 100} {S.Y1 - 3000}))"])
    l_ = S.write_vector(os.path.join(tmp, "lt.gpkg"), [
        f"POLYGON(({S.XC + 1500} {S.Y1 - 3000}, {S.XC + 2500} {S.Y1 - 3000}, {S.XC + 2500} {S.Y1 - 3600}, "
        f"{S.XC + 1500} {S.Y1 - 3600}, {S.XC + 1500} {S.Y1 - 3000}))"])
    nir = [p for p in post if "_B08_" in p][0]
    from osgeo import gdal
    stack = os.path.join(tmp, "stack.vrt")
    gdal.BuildVRT(stack, [green, nir], separate=True)
    r8 = processing.run("flooddepth:supervisedwater", {"IMAGE": stack, "WATER_TRAIN": w, "LAND_TRAIN": l_,
                                                       "MASK": os.path.join(tmp, "sup.tif"),
                                                       "POLYGONS": os.path.join(tmp, "sup.gpkg")})
    r9 = processing.run("flooddepth:compareextents", {"PRED_RASTER": r8["MASK"], "REF_POLYS": ext, "CELL": 10,
                                                      "AGREEMENT": os.path.join(tmp, "agree2.tif")})
    check("supervised CSI > 0.8", r9["CSI"] > 0.8, f"CSI {r9['CSI']:.3f}")

    print("\n9. Refusals with a reason")
    try:
        processing.run("flooddepth:quickflooddepth", {"IMAGES": [green], "DEM": dem, "MEMBERS": 0,
                                                      "OUT_DIR": os.path.join(tmp, "bad")})
        check("one band only is refused", False)
    except Exception as exc:
        check("one band only is refused", "bands" in str(exc).lower(), str(exc)[:90])

    print()
    if FAILURES:
        print(f"{len(FAILURES)} QGIS SMOKE CHECK(S) FAILED: {FAILURES}")
    else:
        print("ALL QGIS SMOKE CHECKS PASSED")
    if app is not None:
        app.exitQgis()
    return not FAILURES


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
