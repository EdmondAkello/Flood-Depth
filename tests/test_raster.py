# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Raster reading, NoData conventions and the analysis grid (GDAL).

    python -m flooddepth.tests.test_raster
"""

import os
import tempfile

import numpy as np

from .common import check, finish
from ..core import raster as R


def _write(path, arrays, gt, epsg, dtype, nodata=None):
    from osgeo import gdal, osr
    gdal.UseExceptions()
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    t = {"Byte": gdal.GDT_Byte, "Float32": gdal.GDT_Float32, "UInt16": gdal.GDT_UInt16}[dtype]
    ds = gdal.GetDriverByName("GTiff").Create(path, arrays[0].shape[1], arrays[0].shape[0], len(arrays), t)
    ds.SetGeoTransform(gt)
    ds.SetProjection(srs.ExportToWkt())
    for i, a in enumerate(arrays, 1):
        ds.GetRasterBand(i).WriteArray(a)
        if nodata is not None:
            ds.GetRasterBand(i).SetNoDataValue(nodata)
    ds = None
    return path


def main():
    d = tempfile.mkdtemp()
    print("\n1. NoData conventions")
    v = np.full((50, 60), 120, "uint8")
    al = np.full((50, 60), 255, "uint8")
    al[:, :10] = 0
    p = _write(os.path.join(d, "raw.tif"), [v, al], (36.0, 1e-4, 0, -1.0, 0, -1e-4), 4326, "Byte")
    a, g, notes = R.read_band(p, 1)
    check("alpha band detected (value + dataMask)", R.alpha_band_index(R.open_dataset(p)) == 2)
    check("alpha 0 -> NaN", np.isnan(a[:, :10]).all() and np.isfinite(a[:, 10:]).all())
    z = np.full((40, 40), 1150.0, "float32")
    z[:5] = -3.4028235e38
    p = _write(os.path.join(d, "hand.tif"), [z], (34.9, 2.8e-4, 0, -0.1, 0, -2.8e-4), 4326, "Float32")
    a, _, notes = R.read_band(p)
    check("undeclared -3.4e38 sentinel -> NaN", np.isnan(a[:5]).all() and np.isfinite(a[5:]).all(), notes[0][:50])
    z = np.full((40, 40), 1150.0, "float32")
    z[:, :8] = 0.0
    p = _write(os.path.join(d, "dem0.tif"), [z], (500000, 30, 0, 9990000, 0, -30), 32736, "Float32")
    a, _, notes = R.read_band(p, zero_nodata_auto=True)
    check("DEM 0 fill with no NoData -> NaN (auto)", np.isnan(a[:, :8]).all(), notes[0][:60] if notes else "")
    a, _, _ = R.read_band(p)
    check("... but only when asked", np.isfinite(a).all())
    lo = np.full((40, 40), 3.0, "float32")
    lo[:, :8] = 0.0
    p = _write(os.path.join(d, "hand0.tif"), [lo], (500000, 30, 0, 9990000, 0, -30), 32736, "Float32")
    a, _, _ = R.read_band(p, zero_nodata_auto=True)
    check("zeros kept when other values are low (HAND, depth)", np.isfinite(a).all())

    print("\n2. Analysis grid")
    z = np.random.default_rng(0).normal(1200, 5, (300, 400)).astype("float32")
    p = _write(os.path.join(d, "glo30.tif"), [z], (34.8, 1 / 3600, 0, -0.1, 0, -1 / 3600), 4326, "Float32")
    g, notes = R.analysis_grid(p)
    check("geographic DEM -> UTM 36S", "32736" in notes[0], notes[0][:70])
    check("cell ~ DEM resolution in metres (31 m)", abs(g.cell_x - 31.0) < 0.6, f"{g.cell_x}")
    p2 = _write(os.path.join(d, "utm.tif"), [z], (600000, 30, 0, 9990000, 0, -30), 32736, "Float32")
    g2, _ = R.analysis_grid(p2, (603015, 9985000, 606000, 9988000), R.srs_from_epsg(32736).ExportToWkt())
    check("projected DEM: cut to the extent on the DEM lattice", (g2.gt[0] - 600000) % 30 == 0 and (9990000 - g2.gt[3]) % 30 == 0)
    check("cut size (snapped outward to whole DEM cells)", g2.cols == 100 and g2.rows == 101, f"{g2.cols} x {g2.rows}")
    try:
        R.analysis_grid(p2, (0, 0, 10, 10), R.srs_from_epsg(32736).ExportToWkt())
        check("no overlap refused", False)
    except R.RasterError:
        check("no overlap refused", True)

    print("\n3. Masks between grids")
    src = R.Grid(30, 30, (600000, 10, 0, 9990300, 0, -10), R.srs_from_epsg(32736).ExportToWkt())
    dst = R.Grid(10, 10, (600000, 30, 0, 9990300, 0, -30), src.wkt)
    m = np.zeros((30, 30), bool)
    m[:, :15] = True
    valid = np.ones_like(m)
    mm, vv = R.mask_to_grid(m, valid, src, dst, 0.5)
    check("area fraction >= 0.5 -> water", mm[:, :5].all() and not mm[:, 5:].any())
    valid[:, 27:] = False
    mm, vv = R.mask_to_grid(m, valid, src, dst, 0.5)
    check("all-NoData destination cell is invalid", not vv[:, 9].any())


if __name__ == "__main__":
    main()
    finish("test_raster")
