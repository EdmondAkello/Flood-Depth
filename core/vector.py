# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Polygons <-> grids with OGR (C++ inside the QGIS process). No QGIS imports."""

import numpy as np

from . import raster as R


def _ogr():
    from osgeo import ogr
    ogr.UseExceptions()
    return ogr


def rasterize_wkts(wkts, src_wkt, grid, all_touched=False):
    """Burn WKT geometries (in CRS `src_wkt`) onto `grid`. Returns a bool mask."""
    gdal = R._gdal()
    ogr = _ogr()
    from osgeo import osr
    mem = gdal.GetDriverByName("MEM").Create("", grid.cols, grid.rows, 1, gdal.GDT_Byte)
    mem.SetGeoTransform(grid.gt)
    mem.SetProjection(grid.wkt)
    if not wkts:
        return np.zeros(grid.shape, bool)
    dst = R.srs_from_wkt(grid.wkt)
    src = R.srs_from_wkt(src_wkt) if src_wkt else dst
    ct = None if src.IsSame(dst) else osr.CoordinateTransformation(src, dst)
    vds = ogr.GetDriverByName("Memory").CreateDataSource("r")
    lyr = vds.CreateLayer("r", srs=dst, geom_type=ogr.wkbUnknown)
    for w in wkts:
        g = ogr.CreateGeometryFromWkt(w)
        if g is None or g.IsEmpty():
            continue
        if ct is not None:
            g.Transform(ct)
        f = ogr.Feature(lyr.GetLayerDefn())
        f.SetGeometry(g)
        lyr.CreateFeature(f)
    opts = ["ALL_TOUCHED=TRUE"] if all_touched else []
    gdal.RasterizeLayer(mem, [1], lyr, burn_values=[1], options=opts)
    return mem.GetRasterBand(1).ReadAsArray() > 0


def read_wkts(path):
    """(list of WKT, layer CRS WKT) from any OGR-readable file (first layer)."""
    ogr = _ogr()
    ds = ogr.Open(str(path))
    if ds is None:
        raise R.RasterError(f"OGR could not open {path}.")
    lyr = ds.GetLayer(0)
    srs = lyr.GetSpatialRef()
    out = []
    for f in lyr:
        g = f.GetGeometryRef()
        if g is not None and not g.IsEmpty():
            out.append(g.ExportToWkt())
    return out, (srs.ExportToWkt() if srs else "")


def _open_for_write(path):
    """(driver, datasource) by extension: .shp, .geojson/.json, else GeoPackage.
    A GeoPackage is opened for update (layers are replaced); other formats
    are recreated."""
    import os
    ogr = _ogr()
    ext = os.path.splitext(str(path))[1].lower()
    name = {".shp": "ESRI Shapefile", ".geojson": "GeoJSON", ".json": "GeoJSON"}.get(ext, "GPKG")
    drv = ogr.GetDriverByName(name)
    if os.path.exists(str(path)):
        if name == "GPKG":
            return drv, drv.Open(str(path), 1)
        drv.DeleteDataSource(str(path))
    return drv, drv.CreateDataSource(str(path))


def polygonize(mask, grid, path, layer="flood_extent", attrs=None, min_cells=0):
    """Write the True regions of `mask` as polygons to a GeoPackage layer.

    Each polygon gets area_m2 (in grid units) and the extra constant
    attributes in `attrs`. Returns the number of polygons written.
    """
    gdal = R._gdal()
    ogr = _ogr()
    mem = gdal.GetDriverByName("MEM").Create("", grid.cols, grid.rows, 1, gdal.GDT_Byte)
    mem.SetGeoTransform(grid.gt)
    mem.SetProjection(grid.wkt)
    b = mem.GetRasterBand(1)
    b.WriteArray(np.asarray(mask, bool).astype("uint8"))
    srs = R.srs_from_wkt(grid.wkt)
    tmp = ogr.GetDriverByName("Memory").CreateDataSource("p")
    tl = tmp.CreateLayer("p", srs=srs, geom_type=ogr.wkbPolygon)
    tl.CreateField(ogr.FieldDefn("v", ogr.OFTInteger))
    gdal.Polygonize(b, b, tl, 0, ["8CONNECTED=8"], callback=None)
    drv, ds = _open_for_write(path)
    if ds.GetLayerByName(layer) is not None:
        ds.DeleteLayer(layer)
    out = ds.CreateLayer(layer, srs=srs, geom_type=ogr.wkbMultiPolygon)
    out.CreateField(ogr.FieldDefn("area_m2", ogr.OFTReal))
    attrs = attrs or {}
    for k, v in attrs.items():
        t = ogr.OFTReal if isinstance(v, float) else (ogr.OFTInteger64 if isinstance(v, int) else ogr.OFTString)
        out.CreateField(ogr.FieldDefn(k, t))
    n = 0
    min_area = min_cells * grid.cell_area()
    out.StartTransaction()
    tl.ResetReading()
    for f in tl:
        if f.GetField("v") != 1:
            continue
        g = f.GetGeometryRef()
        a = g.GetArea()
        if a < min_area:
            continue
        o = ogr.Feature(out.GetLayerDefn())
        o.SetGeometry(ogr.ForceToMultiPolygon(g.Clone()))
        o.SetField("area_m2", float(a))
        for k, v in attrs.items():
            o.SetField(k, v)
        out.CreateFeature(o)
        n += 1
    out.CommitTransaction()
    ds = None
    return n


def write_points(path, layer, grid, rows, cols, fields):
    """Cell-centre points with attribute columns {name: array}."""
    ogr = _ogr()
    import os
    srs = R.srs_from_wkt(grid.wkt)
    drv, ds = _open_for_write(path)
    if ds.GetLayerByName(layer) is not None:
        ds.DeleteLayer(layer)
    lyr = ds.CreateLayer(layer, srs=srs, geom_type=ogr.wkbPoint)
    for name, arr in fields.items():
        t = ogr.OFTInteger if np.asarray(arr).dtype.kind in "iub" else ogr.OFTReal
        lyr.CreateField(ogr.FieldDefn(name, t))
    lyr.StartTransaction()
    for i, (r, c) in enumerate(zip(rows, cols)):
        x, y = grid.rowcol_to_xy(int(r), int(c))
        f = ogr.Feature(lyr.GetLayerDefn())
        g = ogr.Geometry(ogr.wkbPoint)
        g.AddPoint_2D(x, y)
        f.SetGeometry(g)
        for name, arr in fields.items():
            v = arr[i]
            f.SetField(name, int(v) if np.asarray(arr).dtype.kind in "iub" else float(v))
        lyr.CreateFeature(f)
    lyr.CommitTransaction()
    ds = None
    return len(rows)
