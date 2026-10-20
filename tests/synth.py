# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Synthetic test site with a known answer (GDAL only, no QGIS).

A valley 9 km long draining south on a 30 m DEM (WGS 84 / UTM 36S):
    z(x, y) = Z0 - S_LONG * d_south + S_CROSS * |x - XC|
A flood whose water surface is the sloping plane
    wse(y)  = Z0 - S_LONG * d_south + H
so it is exactly the cells within H / S_CROSS of the centreline and the
true depth is wse - z. A small lake (flat bed, water in both dates) sits
outside the valley as permanent water.

Sentinel-2 'Raw' bands are written the way Copernicus Browser exports
them: one Byte file per band, reflectance x 255, a second alpha band, in
EPSG:4326 at ~10 m, file names ..._Sentinel-2_L2A_B03_(Raw).tiff.
"""

import os

import numpy as np

Z0, S_LONG, S_CROSS, H = 1200.0, 0.0005, 0.004, 2.0
ROWS, COLS, CS = 300, 300, 30.0
X0, Y1 = 600000.0, 9950000.0          # UTM 36S
XC = X0 + COLS * CS / 2.0
LAKE = (X0 + 900.0, Y1 - 900.0, 400.0)   # centre x, y, radius (m)
EPSG = 32736


def _gdal():
    from osgeo import gdal, osr
    gdal.UseExceptions()
    return gdal, osr


def truth(x, y):
    d = Y1 - y
    z = Z0 - S_LONG * d + S_CROSS * np.abs(x - XC)
    wse = Z0 - S_LONG * d + H
    lake = (x - LAKE[0]) ** 2 + (y - LAKE[1]) ** 2 < LAKE[2] ** 2
    z = np.where(lake, Z0 + 5.0, z)       # lake bed (flat, DEM shows the water surface)
    return z, wse, lake


def dem_grid():
    cols = np.arange(COLS) * CS + X0 + CS / 2
    rows = Y1 - (np.arange(ROWS) * CS + CS / 2)
    return np.meshgrid(cols, rows)


def write_dem(path, nodata_strip=True):
    gdal, osr = _gdal()
    x, y = dem_grid()
    z, _, _ = truth(x, y)
    z = z.astype("float32")
    if nodata_strip:
        z[:, :3] = -9999.0                 # a NoData strip along the west edge
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    ds = gdal.GetDriverByName("GTiff").Create(path, COLS, ROWS, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((X0, CS, 0, Y1, 0, -CS))
    ds.SetProjection(srs.ExportToWkt())
    b = ds.GetRasterBand(1)
    b.WriteArray(z)
    b.SetNoDataValue(-9999.0)
    ds = None
    return path


def _lonlat_grid(res_deg=9e-5):
    gdal, osr = _gdal()
    utm = osr.SpatialReference()
    utm.ImportFromEPSG(EPSG)
    ll = osr.SpatialReference()
    ll.ImportFromEPSG(4326)
    for s in (utm, ll):
        s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    to_ll = osr.CoordinateTransformation(utm, ll)
    to_utm = osr.CoordinateTransformation(ll, utm)
    pts = [to_ll.TransformPoint(X0 + 300, Y1 - 300), to_ll.TransformPoint(X0 + COLS * CS - 300, Y1 - ROWS * CS + 300)]
    lon0, lat1 = pts[0][0], pts[0][1]
    lon1, lat0 = pts[1][0], pts[1][1]
    nx = int((lon1 - lon0) / res_deg)
    ny = int((lat1 - lat0) / res_deg)
    lon = lon0 + (np.arange(nx) + 0.5) * res_deg
    lat = lat1 - (np.arange(ny) + 0.5) * res_deg
    LON, LAT = np.meshgrid(lon, lat)
    xy = np.array(to_utm.TransformPoints(np.c_[LON.ravel(), LAT.ravel()].tolist()))
    return (lon0, res_deg, 0, lat1, 0, -res_deg), nx, ny, xy[:, 0].reshape(ny, nx), xy[:, 1].reshape(ny, nx)


REFL = {  # (water, land) reflectance
    "B02": (0.07, 0.06), "B03": (0.07, 0.08), "B04": (0.05, 0.07),
    "B08": (0.03, 0.30), "B11": (0.02, 0.24), "B12": (0.01, 0.15),
}


def write_s2_raw(folder, date="2024-05-05", flooded=True, seed=1, dtype="Byte"):
    """Browser-style one-file-per-band export. Returns the file list."""
    gdal, osr = _gdal()
    os.makedirs(folder, exist_ok=True)
    gt, nx, ny, x, y = _lonlat_grid()
    z, wse, lake = truth(x, y)
    water = lake | ((z < wse) if flooded else np.zeros_like(lake))
    rng = np.random.default_rng(seed)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    paths = []
    for band, (rw, rl) in REFL.items():
        r = np.where(water, rw, rl) + rng.normal(0, 0.006, water.shape)
        r = np.clip(r, 0.001, 1.0)
        name = f"{date}-00_00_{date}-23_59_Sentinel-2_L2A_{band}_(Raw).tiff"
        p = os.path.join(folder, name)
        if dtype == "Byte":
            ds = gdal.GetDriverByName("GTiff").Create(p, nx, ny, 2, gdal.GDT_Byte)
            ds.GetRasterBand(1).WriteArray(np.round(r * 255).astype("uint8"))
            alpha = np.full(r.shape, 255, "uint8")
            alpha[:, -5:] = 0                  # a no-data strip on the east edge
            ds.GetRasterBand(2).WriteArray(alpha)
        else:
            ds = gdal.GetDriverByName("GTiff").Create(p, nx, ny, 1, gdal.GDT_Float32)
            ds.GetRasterBand(1).WriteArray(r.astype("float32"))
        ds.SetGeoTransform(gt)
        ds.SetProjection(srs.ExportToWkt())
        ds = None
        paths.append(p)
    return paths


VEG = {"B02": 0.03, "B03": 0.06, "B04": 0.04, "B08": 0.36, "B11": 0.20, "B12": 0.10}
BARE = {"B02": 0.10, "B03": 0.12, "B04": 0.16, "B08": 0.19, "B11": 0.32, "B12": 0.26}


def write_s2_trace(folder, date, stripped, seed=3):
    """Vegetated valley; with `stripped`, the flood strip is bare sediment
    (a receded flash flood). Float32 reflectance, one file per band."""
    gdal, osr = _gdal()
    os.makedirs(folder, exist_ok=True)
    gt, nx, ny, x, y = _lonlat_grid()
    z, wse, lake = truth(x, y)
    strip = (z < wse) & ~lake if stripped else np.zeros_like(lake)
    rng = np.random.default_rng(seed)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    paths = []
    for band in VEG:
        r = np.where(strip, BARE[band], VEG[band]) + rng.normal(0, 0.008, strip.shape)
        p = os.path.join(folder, f"{date}-00_00_{date}-23_59_Sentinel-2_L2A_{band}_(Raw).tiff")
        ds = gdal.GetDriverByName("GTiff").Create(p, nx, ny, 1, gdal.GDT_Float32)
        ds.GetRasterBand(1).WriteArray(np.clip(r, 0.001, 1).astype("float32"))
        ds.SetGeoTransform(gt)
        ds.SetProjection(srs.ExportToWkt())
        ds = None
        paths.append(p)
    return paths


def flood_polygon_wkt():
    """The true flood strip as a UTM polygon (clipped to the DEM)."""
    half = H / S_CROSS
    return (f"POLYGON(({XC - half} {Y1 - 5}, {XC + half} {Y1 - 5}, {XC + half} {Y1 - ROWS * CS + 5}, "
            f"{XC - half} {Y1 - ROWS * CS + 5}, {XC - half} {Y1 - 5}))")


def utm_wkt():
    _, osr = _gdal()
    s = osr.SpatialReference()
    s.ImportFromEPSG(EPSG)
    return s.ExportToWkt()


def write_vector(path, wkts, geom="polygon"):
    gdal, osr = _gdal()
    from osgeo import ogr
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    drv = ogr.GetDriverByName("GPKG")
    if os.path.exists(path):
        drv.DeleteDataSource(path)
    ds = drv.CreateDataSource(path)
    gt = {"polygon": ogr.wkbPolygon, "line": ogr.wkbLineString, "point": ogr.wkbPoint}[geom]
    lyr = ds.CreateLayer("layer", srs=srs, geom_type=gt)
    for w in wkts:
        f = ogr.Feature(lyr.GetLayerDefn())
        f.SetGeometry(ogr.CreateGeometryFromWkt(w))
        lyr.CreateFeature(f)
    ds = None
    return path


def centreline_wkt():
    return f"LINESTRING({XC} {Y1 - 10}, {XC} {Y1 - ROWS * CS + 10})"
