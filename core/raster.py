# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""GDAL raster I/O on a common analysis grid. No QGIS imports.

Everything the algorithms see is a float32 array with NaN where the data
are not valid, on a `Grid`. "Not valid" collects every convention found in
real downloads, because missing one silently turns NoData into terrain or
into water:

* the band's declared NoData value;
* non-finite values;
* undeclared float sentinels such as -3.4028e38 (|v| > 1e30), common in
  Copernicus DEM and GLO-30 HAND clips;
* an alpha / data-mask band (GDAL colour interpretation "Alpha", or the
  last band of a Sentinel Hub / Copernicus Browser export whose values are
  only 0 and 255 or 0 and 1).

Reprojection uses gdal.Warp on in-memory datasets: a C++ call inside the
QGIS process, never a subprocess.
"""

import math
from dataclasses import dataclass

import numpy as np

SENTINEL_ABS = 1e30


def _gdal():
    from osgeo import gdal
    gdal.UseExceptions()
    return gdal


def _osr():
    from osgeo import osr
    osr.UseExceptions()
    return osr


def srs_from_wkt(wkt):
    osr = _osr()
    s = osr.SpatialReference()
    s.ImportFromWkt(wkt)
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return s


def srs_from_epsg(code):
    osr = _osr()
    s = osr.SpatialReference()
    s.ImportFromEPSG(int(code))
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return s


def is_geographic(wkt):
    return bool(wkt) and bool(srs_from_wkt(wkt).IsGeographic())


@dataclass
class Grid:
    """A north-up raster grid: size, geotransform and CRS (WKT)."""
    cols: int
    rows: int
    gt: tuple
    wkt: str

    @property
    def cell_x(self):
        return abs(self.gt[1])

    @property
    def cell_y(self):
        return abs(self.gt[5])

    @property
    def bounds(self):
        x0, y0 = self.gt[0], self.gt[3]
        x1 = x0 + self.cols * self.gt[1]
        y1 = y0 + self.rows * self.gt[5]
        return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))

    @property
    def shape(self):
        return (self.rows, self.cols)

    def cell_area(self):
        return self.cell_x * self.cell_y

    def rowcol_to_xy(self, row, col):
        return (self.gt[0] + (col + 0.5) * self.gt[1],
                self.gt[3] + (row + 0.5) * self.gt[5])

    def xy_to_rowcol(self, x, y):
        return (int(math.floor((y - self.gt[3]) / self.gt[5])),
                int(math.floor((x - self.gt[0]) / self.gt[1])))

    def same_as(self, other, tol=1e-6):
        if (self.cols, self.rows) != (other.cols, other.rows):
            return False
        if any(abs(a - b) > tol * max(1.0, abs(a)) for a, b in zip(self.gt, other.gt)):
            return False
        return srs_from_wkt(self.wkt).IsSame(srs_from_wkt(other.wkt)) == 1


class RasterError(Exception):
    """An input raster cannot be used; the message says why and what to do."""


class quiet(object):
    """Silence GDAL warnings (e.g. TIFF ExtraSamples notices in Browser exports)."""

    def __enter__(self):
        _gdal().PushErrorHandler("CPLQuietErrorHandler")
        return self

    def __exit__(self, *a):
        _gdal().PopErrorHandler()
        return False


def open_dataset(path):
    gdal = _gdal()
    try:
        with quiet():
            ds = gdal.Open(str(path), gdal.GA_ReadOnly)
    except RuntimeError as exc:
        raise RasterError(f"GDAL could not open {path}: {exc}")
    if ds is None:
        raise RasterError(f"GDAL could not open {path}.")
    return ds


def grid_of(ds_or_path):
    ds = open_dataset(ds_or_path) if isinstance(ds_or_path, str) else ds_or_path
    gt = ds.GetGeoTransform()
    if gt[2] != 0 or gt[4] != 0:
        raise RasterError("Rotated rasters are not supported. Warp the raster to a north-up grid first.")
    wkt = ds.GetProjection()
    if not wkt:
        raise RasterError(f"{ds.GetDescription()} has no coordinate reference system.")
    return Grid(ds.RasterXSize, ds.RasterYSize, tuple(gt), wkt)


def alpha_band_index(ds):
    """1-based index of an alpha / data-mask band, or None.

    Declared alpha (colour interpretation) wins. Otherwise the LAST band of
    a 2- or 4-band integer raster counts as a data mask when it holds only
    {0, 255} or {0, 1} - the layout of Sentinel Hub / Copernicus Browser
    "Raw" exports (value + dataMask) and of RGBA screenshots.
    """
    gdal = _gdal()
    n = ds.RasterCount
    for i in range(1, n + 1):
        if ds.GetRasterBand(i).GetColorInterpretation() == gdal.GCI_AlphaBand:
            return i
    if n in (2, 4):
        b = ds.GetRasterBand(n)
        if b.DataType in (gdal.GDT_Byte, gdal.GDT_UInt16, gdal.GDT_Int16):
            step = max(1, int(math.sqrt(ds.RasterXSize * ds.RasterYSize / 250000.0)))
            with quiet():
                a = b.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize,
                                  buf_xsize=max(1, ds.RasterXSize // step),
                                  buf_ysize=max(1, ds.RasterYSize // step))
            vals = np.unique(a)
            if vals.size <= 2 and set(vals.tolist()) <= {0, 1, 255} and (vals.max() > 0):
                return n
    return None


def data_band_count(ds):
    a = alpha_band_index(ds)
    return ds.RasterCount - (1 if a else 0)


def read_band(path, band=1, ds=None, zero_nodata_auto=False):
    """(float32 array with NaN for invalid cells, Grid, notes list).

    zero_nodata_auto: for elevation rasters with no declared NoData, treat
    exact zeros as NoData when they cover > 0.5 % of the grid and every
    other value is above 50 (a DEM clipped with 0 fill, a common export)."""
    ds = ds or open_dataset(path)
    if band < 1 or band > ds.RasterCount:
        raise RasterError(f"{path} has {ds.RasterCount} band(s); band {band} does not exist.")
    b = ds.GetRasterBand(band)
    try:
        with quiet():
            arr = b.ReadAsArray()
    except RuntimeError as exc:
        raise RasterError(f"Could not read band {band} of {path} (file damaged or truncated?): {exc}")
    arr = arr.astype("float32")
    notes = []
    nd = b.GetNoDataValue()
    invalid = ~np.isfinite(arr)
    if nd is not None and math.isfinite(nd):
        invalid |= arr == np.float32(nd)
    sentinel = np.abs(arr) > SENTINEL_ABS
    if sentinel.any() and (nd is None or abs(nd) < SENTINEL_ABS):
        notes.append(f"{int(sentinel.sum()):,} cells hold an undeclared float sentinel (|v| > 1e30) and were treated as NoData.")
    invalid |= sentinel
    ai = alpha_band_index(ds)
    if ai and ai != band:
        with quiet():
            alpha = ds.GetRasterBand(ai).ReadAsArray()
        invalid |= alpha == 0
        notes.append(f"band {ai} used as the data mask (alpha).")
    if zero_nodata_auto and nd is None:
        z = (arr == 0) & ~invalid
        rest = arr[~invalid & ~z]
        if z.mean() > 0.005 and rest.size and float(rest.min()) > 50.0:
            invalid |= z
            notes.append(f"no NoData value declared; {int(z.sum()):,} cells of exactly 0 (all other "
                         f"values >= {float(rest.min()):.0f}) treated as NoData.")
    arr[invalid] = np.nan
    return arr, grid_of(ds), notes


def band_info(path):
    """Light description used for band identification and error messages."""
    gdal = _gdal()
    ds = open_dataset(path)
    alpha = alpha_band_index(ds)
    bands = []
    for i in range(1, ds.RasterCount + 1):
        b = ds.GetRasterBand(i)
        bands.append({"index": i, "description": b.GetDescription() or "",
                      "dtype": gdal.GetDataTypeName(b.DataType),
                      "colorinterp": gdal.GetColorInterpretationName(b.GetColorInterpretation()),
                      "metadata": b.GetMetadata() or {},
                      "alpha": i == alpha})
    return {"path": str(path), "count": ds.RasterCount, "alpha": alpha,
            "bands": bands, "metadata": ds.GetMetadata() or {}, "grid": grid_of(ds)}


def transform_bounds(bounds, src_wkt, dst_wkt, densify=21):
    """Bounding box of `bounds` (xmin, ymin, xmax, ymax) reprojected."""
    osr = _osr()
    src, dst = srs_from_wkt(src_wkt), srs_from_wkt(dst_wkt)
    if src.IsSame(dst):
        return tuple(bounds)
    ct = osr.CoordinateTransformation(src, dst)
    xmin, ymin, xmax, ymax = bounds
    xs, ys = [], []
    for t in np.linspace(0.0, 1.0, densify):
        for (x, y) in ((xmin + t * (xmax - xmin), ymin), (xmin + t * (xmax - xmin), ymax),
                       (xmin, ymin + t * (ymax - ymin)), (xmax, ymin + t * (ymax - ymin))):
            px, py, _ = ct.TransformPoint(float(x), float(y))
            if math.isfinite(px) and math.isfinite(py):
                xs.append(px)
                ys.append(py)
    if not xs:
        raise RasterError("Could not transform the raster extent between coordinate systems.")
    return (min(xs), min(ys), max(xs), max(ys))


def intersect(a, b):
    out = (max(a[0], b[0]), max(a[1], b[1]), min(a[2], b[2]), min(a[3], b[3]))
    if out[0] >= out[2] or out[1] >= out[3]:
        return None
    return out


def utm_epsg(lon, lat):
    zone = int(math.floor((lon + 180.0) / 6.0)) % 60 + 1
    return (32600 if lat >= 0 else 32700) + zone


def metres_per_degree(lat):
    return 111320.0 * math.cos(math.radians(lat))


def analysis_grid(dem_path, extent=None, extent_wkt=None, cell=0.0, max_cells=None):
    """The projected grid every analysis step shares.

    * DEM in a projected CRS: the DEM's own grid (cells are never shifted),
      cut to `extent` (any CRS) when given, optionally re-celled to `cell`.
    * DEM in geographic coordinates (Copernicus GLO-30, HAND): the WGS 84 /
      UTM zone of the extent centre, at the DEM resolution in metres rounded
      to 0.5 m unless `cell` is given.
    Returns (Grid, notes).
    """
    ds = open_dataset(dem_path)
    g = grid_of(ds)
    notes = []
    dem_bounds = g.bounds
    if is_geographic(g.wkt):
        cx = (dem_bounds[0] + dem_bounds[2]) / 2.0
        cy = (dem_bounds[1] + dem_bounds[3]) / 2.0
        if extent is not None:
            eb = transform_bounds(extent, extent_wkt or g.wkt, g.wkt)
            ib = intersect(eb, dem_bounds)
            if ib is None:
                raise RasterError("The analysis extent does not overlap the DEM.")
            cx, cy = (ib[0] + ib[2]) / 2.0, (ib[1] + ib[3]) / 2.0
        epsg = utm_epsg(cx, cy)
        wkt = srs_from_epsg(epsg).ExportToWkt()
        native = min(g.cell_x * metres_per_degree(cy), g.cell_y * 111320.0)
        cs = float(cell) if cell and cell > 0 else max(0.5, round(native * 2.0) / 2.0)
        b = transform_bounds(dem_bounds, g.wkt, wkt)
        if extent is not None:
            b2 = transform_bounds(extent, extent_wkt or g.wkt, wkt)
            b = intersect(b, b2)
            if b is None:
                raise RasterError("The analysis extent does not overlap the DEM.")
        x0 = math.floor(b[0] / cs) * cs
        y1 = math.ceil(b[3] / cs) * cs
        cols = int(math.ceil((b[2] - x0) / cs))
        rows = int(math.ceil((y1 - b[1]) / cs))
        out = Grid(cols, rows, (x0, cs, 0.0, y1, 0.0, -cs), wkt)
        notes.append(f"DEM is geographic; analysis grid is EPSG:{epsg} (WGS 84 / UTM) at {cs:g} m "
                     f"(DEM native about {native:.1f} m).")
    else:
        cs_x, cs_y = g.cell_x, g.cell_y
        b = dem_bounds
        if extent is not None:
            b2 = transform_bounds(extent, extent_wkt or g.wkt, g.wkt)
            b = intersect(b, b2)
            if b is None:
                raise RasterError("The analysis extent does not overlap the DEM.")
        if cell and cell > 0:
            cs_x = cs_y = float(cell)
        # snap to the DEM lattice so DEM cells are not shifted
        gx0, gy0 = g.gt[0], g.gt[3]
        c0 = math.floor((b[0] - gx0) / g.cell_x)
        r0 = math.floor((gy0 - b[3]) / g.cell_y)
        x0 = gx0 + c0 * g.cell_x
        y1 = gy0 - r0 * g.cell_y
        cols = int(math.ceil((b[2] - x0) / cs_x))
        rows = int(math.ceil((y1 - b[1]) / cs_y))
        out = Grid(cols, rows, (x0, cs_x, 0.0, y1, 0.0, -cs_y), g.wkt)
        notes.append(f"Analysis grid: DEM grid, {cs_x:g} x {cs_y:g} map units.")
    if out.cols < 3 or out.rows < 3:
        raise RasterError("The analysis area is smaller than 3 x 3 cells.")
    if max_cells and out.cols * out.rows > max_cells:
        raise RasterError(f"The analysis grid would have {out.cols * out.rows:,} cells (limit {max_cells:,}). "
                          "Give an area of interest or a larger cell size.")
    notes.append(f"Analysis grid {out.cols:,} x {out.rows:,} cells.")
    return out, notes


def _mem_from_array(arr, grid):
    gdal = _gdal()
    mem = gdal.GetDriverByName("MEM").Create("", grid.cols, grid.rows, 1, gdal.GDT_Float32)
    mem.SetGeoTransform(grid.gt)
    mem.SetProjection(grid.wkt)
    b = mem.GetRasterBand(1)
    b.SetNoDataValue(float("nan"))
    b.WriteArray(np.asarray(arr, dtype="float32"))
    return mem


def warp_array(arr, src_grid, dst_grid, resample="bilinear"):
    """Warp a NaN-masked float array between grids (NaN stays NaN)."""
    gdal = _gdal()
    if src_grid.same_as(dst_grid):
        return np.array(arr, dtype="float32", copy=True)
    src = _mem_from_array(arr, src_grid)
    dst = gdal.GetDriverByName("MEM").Create("", dst_grid.cols, dst_grid.rows, 1, gdal.GDT_Float32)
    dst.SetGeoTransform(dst_grid.gt)
    dst.SetProjection(dst_grid.wkt)
    db = dst.GetRasterBand(1)
    db.SetNoDataValue(float("nan"))
    db.Fill(float("nan"))
    gdal.Warp(dst, src, resampleAlg=resample, srcNodata=float("nan"), dstNodata=float("nan"),
              multithread=False)
    return db.ReadAsArray().astype("float32")


def read_on_grid(path, grid, band=1, resample="bilinear", zero_nodata_auto=False):
    """Read one band, mask it, and warp it onto `grid`. Returns (array, notes)."""
    arr, g, notes = read_band(path, band, zero_nodata_auto=zero_nodata_auto)
    return warp_array(arr, g, grid, resample), notes


def mask_to_grid(mask, valid, src_grid, dst_grid, min_fraction=0.5):
    """Binary mask between grids by area fraction: a destination cell is
    True when at least `min_fraction` of its valid area is True.
    Returns (mask, valid) on dst_grid."""
    m = np.where(valid, mask.astype("float32"), np.nan)
    if src_grid.same_as(dst_grid):
        return mask & valid, valid.copy()
    frac = warp_array(m, src_grid, dst_grid, "average")
    v = np.isfinite(frac)
    return v & (np.nan_to_num(frac) >= min_fraction), v


def write_raster(path, arr, grid, nodata=-9999.0, dtype="Float32", tags=None, color_table=None):
    gdal = _gdal()
    types = {"Float32": gdal.GDT_Float32, "Byte": gdal.GDT_Byte, "Int16": gdal.GDT_Int16,
             "UInt16": gdal.GDT_UInt16, "Int32": gdal.GDT_Int32}
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(str(path), grid.cols, grid.rows, 1, types[dtype],
                    options=["TILED=YES", "COMPRESS=DEFLATE", "BIGTIFF=IF_SAFER"])
    if ds is None:
        raise RasterError(f"Could not create {path}.")
    ds.SetGeoTransform(grid.gt)
    ds.SetProjection(grid.wkt)
    b = ds.GetRasterBand(1)
    if color_table:  # must precede the first write
        ct = gdal.ColorTable()
        for k, rgba in color_table.items():
            ct.SetColorEntry(int(k), tuple(rgba))
        b.SetRasterColorTable(ct)
        b.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
    out = np.asarray(arr)
    if dtype == "Float32":
        out = np.where(np.isfinite(out), out, nodata).astype("float32")
    b.WriteArray(out)
    if nodata is not None:
        b.SetNoDataValue(float(nodata))
    if tags:
        ds.SetMetadata({str(k): str(v) for k, v in tags.items()})
    b.FlushCache()
    ds = None
    return str(path)
