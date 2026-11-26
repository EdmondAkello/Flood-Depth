# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""DEM tiles, DEM audit and Prepare DEM (GDAL).

    python -m flooddepth.tests.test_dem
"""

import json
import os
import tempfile

import numpy as np

from .common import check, finish
from ..core import dem as DM
from ..core import raster as R


def _tile(path, arr, gt, epsg):
    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.GetDriverByName("GTiff").Create(path, arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32,
                                              ["TILED=YES", "COMPRESS=LZW"])
    ds.SetGeoTransform(gt)
    ds.SetProjection(R.srs_from_epsg(epsg).ExportToWkt())
    ds.GetRasterBand(1).WriteArray(arr.astype("float32"))
    ds = None
    return path


def main():
    d = tempfile.mkdtemp()
    yy, xx = np.mgrid[0:200, 0:200].astype("float32")
    z = 1200.0 + 0.3 * xx + 0.1 * yy + 2.0 * np.sin(xx / 9.0)

    print("\n1. Tiles: mosaic and choice of coordinate system")
    a = _tile(os.path.join(d, "gem15_a.tif"), z[:, :100], (500000, 15, 0, 9900000, 0, -15), 32737)
    b = _tile(os.path.join(d, "gem15_b.tif"), z[:, 100:], (501500, 15, 0, 9900000, 0, -15), 32737)
    other = _tile(os.path.join(d, "gem15_z36.tif"), z, (166000, 15, 0, 9900000, 0, -15), 32736)
    src, notes = DM.dem_source([a])
    check("one DEM is used as it is", src == a and not notes)
    ext = (501300, 9898000, 501800, 9899000)
    src, notes = DM.dem_source([a, b, other], ext, R.srs_from_epsg(32737).ExportToWkt())
    check("tiles of the CRS covering the area are mosaicked", src.startswith("/vsimem/") and len(notes) == 2,
          "; ".join(notes))
    g = R.Grid(40, 20, (501200, 15, 0, 9898900, 0, -15), R.srs_from_epsg(32737).ExportToWkt())
    m, _ = R.read_on_grid(src, g, 1, "near")
    check("the mosaic reads across the tile seam", np.allclose(m, z[73:93, 80:120]))

    print("\n2. DEM audit")
    w, fig = DM.audit(z)
    check("clean continuous DEM: no warning", not w and fig["vertical_step_m"] is None, str(fig))
    rep = np.repeat(z[:, :100], 2, axis=1)
    w, fig = DM.audit(rep)
    check("nearest-neighbour repeats are flagged", len(w) == 1 and fig["repeated_col_share"] > 0.4)
    flat = np.round(1100.0 + 0.002 * xx)
    w, fig = DM.audit(flat)
    check("whole-metre steps on flat ground are flagged", any("whole steps" in x for x in w), str(fig))

    print("\n3. Prepare DEM")
    ll = R.srs_from_epsg(32737).ExportToWkt()
    out = os.path.join(d, "prepared.tif")
    info = DM.prepare_dem([a, b], out, clip_bounds=ext, clip_wkt=ll, buffer_m=300.0)
    pg = R.grid_of(out)
    check("prepared DEM is in UTM 37S at the native 15 m", "37S" in info["crs"] and pg.cell_x == 15.0,
          f"{info['crs']}, {pg.cell_x}")
    check("clipped to the area plus the buffer", abs(pg.bounds[0] - (ext[0] - 300)) <= 15 and
          abs(pg.bounds[2] - (ext[2] + 300)) <= 15, str(pg.bounds))
    check("log written next to it", json.load(open(os.path.splitext(out)[0] + ".json"))["cell_m"] == 15.0)
    try:
        DM.prepare_dem([a], out, resampling="near")
        check("nearest neighbour is refused", False)
    except ValueError:
        check("nearest neighbour is refused", True)

    print("\n4. Memory check")
    need, w = DM.memory_check(10_000_000, available_gb=8.0)
    check("10 M cells fit in 8 GB without a warning", w is None and 0.3 < need < 1.0, f"{need:.2f} GB")
    need, w = DM.memory_check(80_000_000, available_gb=8.0)
    check("80 M cells in 8 GB: warning", w is not None)
    try:
        DM.memory_check(1_073_741_824, available_gb=8.0)     # a whole 32 768 x 32 768 GEM tile
        check("a whole GEM tile is refused with advice", False)
    except ValueError as exc:
        check("a whole GEM tile is refused with advice", "area of interest" in str(exc))


if __name__ == "__main__":
    main()
    finish("test_dem")
