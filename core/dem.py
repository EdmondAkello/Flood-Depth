# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""DEM preparation: tiles, checks and a prepared GeoTIFF. No QGIS imports.

National DEMs come as tiles, often several CRSs at once (the Google GEM
15 m DTM for Kenya ships the whole country in four UTM zones, 24 files).
`dem_source` turns a list of tiles into one virtual raster (gdal.BuildVRT
in /vsimem: an XML file, nothing is copied) of the tiles that share the
CRS best suited to the area, so a user can select the whole folder.

`audit` flags the DEM defects that move a waterline elevation and so the
depth: repeated rows / columns from nearest-neighbour resampling, and
elevations stored in coarse steps (whole metres) on flat ground. The
checks follow QEHT's DEM audit (same author), tuned to flood depth.

`prepare_dem` merges, clips (with a buffer), reprojects to the local UTM
zone with bilinear resampling and audits, as QEHT's "Prepare DEM for
hydrology" does. There is deliberately no depression filling: filling
raises the ground in exactly the sinks where flood water stands, and the
depth is water surface minus the real ground.
"""

import json
import math
import os

import numpy as np

from . import raster as R

AUDIT_MAX_CELLS = 4_000_000
OUT_NODATA = -9999.0
_vrt_count = [0]


def _gdal():
    return R._gdal()


def _group_by_crs(paths):
    groups = []   # [(srs, [paths], bounds)]
    for p in paths:
        ds = R.open_dataset(p)
        g = R.grid_of(ds)
        srs = R.srs_from_wkt(g.wkt)
        for grp in groups:
            if grp[0].IsSame(srs):
                grp[1].append(p)
                b = grp[2]
                grp[2] = (min(b[0], g.bounds[0]), min(b[1], g.bounds[1]),
                          max(b[2], g.bounds[2]), max(b[3], g.bounds[3]))
                break
        else:
            groups.append([srs, [p], g.bounds])
    return groups


def _epsg(srs):
    try:
        srs.AutoIdentifyEPSG()
    except RuntimeError:
        pass
    code = srs.GetAuthorityCode(None)
    return int(code) if code and code.isdigit() else None


def _overlap_share(bounds, extent, extent_wkt, wkt):
    if extent is None:
        return 0.0
    try:
        e = R.transform_bounds(extent, extent_wkt or wkt, wkt)
    except R.RasterError:
        return 0.0
    i = R.intersect(bounds, e)
    if i is None:
        return 0.0
    return (i[2] - i[0]) * (i[3] - i[1]) / max(1e-9, (e[2] - e[0]) * (e[3] - e[1]))


def dem_source(paths, extent=None, extent_wkt=None):
    """(path GDAL can open, notes): the single DEM, or a VRT of the tiles.

    With tiles in several CRSs, the group used is the one covering the
    extent best; ties go to the UTM zone of the extent centre."""
    paths = [str(p) for p in paths if p]
    if not paths:
        raise ValueError("Give the DEM.")
    if len(paths) == 1:
        return paths[0], []
    notes = []
    groups = _group_by_crs(paths)
    if len(groups) > 1:
        centre_utm = None
        if extent is not None:
            ll = R.srs_from_epsg(4326).ExportToWkt()
            try:
                b = R.transform_bounds(extent, extent_wkt, ll)
                centre_utm = R.utm_epsg((b[0] + b[2]) / 2.0, (b[1] + b[3]) / 2.0)
            except R.RasterError:
                pass

        def score(grp):
            share = _overlap_share(grp[2], extent, extent_wkt, grp[0].ExportToWkt())
            return (round(share, 3), _epsg(grp[0]) == centre_utm, len(grp[1]))
        groups.sort(key=score, reverse=True)
        best = groups[0]
        notes.append(f"DEM tiles in {len(groups)} coordinate systems: the {len(best[1])} tile(s) in "
                     f"{best[0].GetName()} were used (best cover of the area); the other "
                     f"{sum(len(g[1]) for g in groups[1:])} were ignored.")
    use = groups[0][1]
    if len(use) == 1:
        return use[0], notes
    gdal = _gdal()
    _vrt_count[0] += 1
    vrt = f"/vsimem/flooddepth_dem_{os.getpid()}_{_vrt_count[0]}.vrt"
    with R.quiet():
        ds = gdal.BuildVRT(vrt, use)
    if ds is None:
        raise R.RasterError("The DEM tiles could not be combined (gdal.BuildVRT).")
    ds = None
    notes.append(f"{len(use)} DEM tiles combined into one virtual mosaic (only the cells under the "
                 "analysis area are read).")
    return vrt, notes


# -- memory --------------------------------------------------------------

BYTES_PER_CELL = 64   # measured peak growth per analysis cell (Nyando, 15 m vs 30 m runs: ~52 B)


def available_memory_gb():
    """Free physical memory in GB (best effort; None when unknown). From QEHT."""
    try:
        if os.name == "nt":
            import ctypes

            class _MS(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                            ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                            ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                            ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                            ("sullAvailExtendedVirtual", ctypes.c_ulonglong)]
            ms = _MS()
            ms.dwLength = ctypes.sizeof(_MS)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms)):
                return ms.ullAvailPhys / 1e9
            return None
        return os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 1e9
    except (AttributeError, ValueError, OSError):
        return None


def memory_check(cells, available_gb=None):
    """(estimate GB, warning or None); raises ValueError when the run would
    clearly not fit in the free memory."""
    need = cells * BYTES_PER_CELL / 1e9
    free = available_memory_gb() if available_gb is None else available_gb
    if free is None:
        return need, None
    advice = ("Give an area of interest around the flood, or a larger analysis cell size "
              "(advanced parameter), or close other programs.")
    if need > 0.9 * free:
        raise ValueError(f"The analysis grid has {cells / 1e6:,.0f} M cells and needs about {need:.1f} GB of "
                         f"memory, but only {free:.1f} GB is free. " + advice)
    if need > 0.5 * free:
        return need, (f"The run needs about {need:.1f} GB of the {free:.1f} GB of free memory; QGIS may slow "
                      "down. " + advice)
    return need, None


# -- audit ---------------------------------------------------------------

def _central(a, max_cells=AUDIT_MAX_CELLS):
    rows, cols = a.shape
    if rows * cols <= max_cells:
        return a
    side = int(math.sqrt(max_cells))
    r0, c0 = (rows - side) // 2, (cols - side) // 2
    return a[r0:r0 + side, c0:c0 + side]


def _repeats(z, min_overlap=50):
    """(rows repeating the previous row, rows considered)."""
    v = np.isfinite(z)
    ov = v[1:] & v[:-1]
    n = ov.sum(axis=1)
    eq = ((z[1:] == z[:-1]) & ov).sum(axis=1)
    with np.errstate(invalid="ignore"):
        lo = np.where(ov, z[1:], np.inf).min(axis=1)
        hi = np.where(ov, z[1:], -np.inf).max(axis=1)
    dup = (n >= min_overlap) & (eq >= 0.99 * n) & (hi > lo)
    return int(dup.sum()), int((n >= min_overlap).sum())


def vertical_step(z, candidates=(1.0, 0.5, 0.25, 0.1, 0.05, 0.01), need=0.999, sample=200_000):
    """Storage step of the elevations (m), or None when continuous."""
    x = z[np.isfinite(z)].astype("float64")
    if x.size == 0:
        return None
    if x.size > sample:
        x = x[np.random.default_rng(0).choice(x.size, sample, replace=False)]
    eps = 4.0 * float(np.spacing(np.float32(np.max(np.abs(x)))))
    for q in candidates:
        if (np.abs(x - np.round(x / q) * q) <= max(0.02 * q, eps)).mean() >= need:
            return float(q)
    return None


def audit(dem, flood=None):
    """Warnings about the DEM on the analysis grid (list of str) and a dict
    of the figures. `flood` (optional) restricts the flat-ground check to
    the flood and its surroundings."""
    z = _central(np.asarray(dem, dtype="float32"))
    out = {"repeated_row_share": 0.0, "repeated_col_share": 0.0, "vertical_step_m": None}
    warns = []
    if np.isfinite(z).sum() < 100:
        return warns, out
    dr, nr = _repeats(z)
    dc, nc = _repeats(z.T)
    out["repeated_row_share"] = dr / max(1, nr)
    out["repeated_col_share"] = dc / max(1, nc)
    if max(out["repeated_row_share"], out["repeated_col_share"]) >= 0.002:
        warns.append(f"{100 * out['repeated_row_share']:.1f} % of DEM rows and "
                     f"{100 * out['repeated_col_share']:.1f} % of columns exactly repeat their neighbour: the DEM "
                     "was probably resampled with nearest neighbour (or exported finer than its native grid). "
                     "The repeats make steps in the terrain that the waterline elevations pick up. Use the "
                     "native product, or the Prepare DEM tool (bilinear).")
    q = vertical_step(z)
    out["vertical_step_m"] = q
    if q is not None and q >= 0.5:
        area = z if flood is None else np.where(_central(np.asarray(flood, bool)), z, np.nan)
        gy, gx = np.gradient(area.astype("float64"))
        drop = float(np.nanmedian(np.hypot(gx, gy))) if np.isfinite(gy).any() else 0.0
        out["median_drop_per_cell_m"] = drop
        if drop < q:
            warns.append(f"DEM elevations are stored in whole steps of {q:g} m, larger than the typical "
                         f"change from one cell to the next on the flood ({drop:.2f} m). Waterline elevations "
                         "on flat ground are then rounded to the step and the depth with them; the "
                         "uncertainty bands include this only through the DEM error. Prefer a "
                         "floating-point DEM (FABDEM, GEM 15 m, LiDAR).")
    return warns, out


# -- prepare ---------------------------------------------------------------

def prepare_dem(paths, out_path, clip_bounds=None, clip_wkt=None, buffer_m=1000.0, cell_m=0.0,
                target_wkt=None, resampling="bilinear", progress=None):
    """Merge DEM tiles, clip to `clip_bounds` (+ buffer), reproject to
    `target_wkt` (default: the WGS 84 / UTM zone of the area) and audit.
    Returns a dict (also written next to the output as <name>.json)."""
    gdal = _gdal()
    if resampling not in ("bilinear", "cubic"):
        raise ValueError("Resampling must be bilinear or cubic: nearest neighbour repeats rows and columns "
                         "and puts steps in the terrain.")
    src, notes = dem_source(paths, clip_bounds, clip_wkt)
    ds = R.open_dataset(src)
    g = R.grid_of(ds)
    ll = R.srs_from_epsg(4326).ExportToWkt()
    area = clip_bounds if clip_bounds is not None else g.bounds
    area_wkt = (clip_wkt or g.wkt) if clip_bounds is not None else g.wkt
    b = R.transform_bounds(area, area_wkt, ll)
    lon, lat = (b[0] + b[2]) / 2.0, (b[1] + b[3]) / 2.0
    if not target_wkt:
        epsg = R.utm_epsg(lon, lat)
        target_wkt = R.srs_from_epsg(epsg).ExportToWkt()
    dst = R.srs_from_wkt(target_wkt)
    if dst.IsGeographic():
        raise ValueError("The target CRS is geographic; FloodDepth works in metres (choose a projected CRS).")
    if R.is_geographic(g.wkt):
        native = min(g.cell_x * R.metres_per_degree(lat), g.cell_y * 111320.0)
    else:
        native = min(g.cell_x, g.cell_y) * (R.srs_from_wkt(g.wkt).GetLinearUnits() or 1.0)
    cell = float(cell_m) if cell_m and cell_m > 0 else max(0.5, round(native * 2.0) / 2.0)
    opts = dict(format="GTiff", dstSRS=target_wkt, xRes=cell, yRes=cell, targetAlignedPixels=True,
                resampleAlg=resampling, outputType=gdal.GDT_Float32, dstNodata=OUT_NODATA,
                creationOptions=["COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=3", "BIGTIFF=IF_SAFER"],
                multithread=True, warpMemoryLimit=512)
    nd = ds.GetRasterBand(1).GetNoDataValue()
    opts["srcNodata"] = nd if nd is not None else float("nan")
    if clip_bounds is not None:
        cb = R.transform_bounds(clip_bounds, clip_wkt or g.wkt, target_wkt)
        bf = float(buffer_m or 0.0)
        opts["outputBounds"] = (cb[0] - bf, cb[1] - bf, cb[2] + bf, cb[3] + bf)
    same_grid = (not R.is_geographic(g.wkt) and dst.IsSame(R.srs_from_wkt(g.wkt))
                 and abs(cell - g.cell_x) < 1e-6 and abs(cell - g.cell_y) < 1e-6)
    if same_grid:
        # already in the target CRS at the target cell: keep the source lattice, so no cell is resampled
        b = opts.get("outputBounds", g.bounds)
        x0 = g.gt[0] + math.floor((b[0] - g.gt[0]) / cell) * cell
        x1 = g.gt[0] + math.ceil((b[2] - g.gt[0]) / cell) * cell
        y1 = g.gt[3] - math.floor((g.gt[3] - b[3]) / cell) * cell
        y0 = g.gt[3] - math.ceil((g.gt[3] - b[1]) / cell) * cell
        opts["outputBounds"] = (x0, y0, x1, y1)
        opts["targetAlignedPixels"] = False
        notes.append("The DEM is already in the target CRS at the target cell size: its cells were copied, "
                     "not resampled.")
    if progress is not None:
        opts["callback"] = lambda frac, msg, data: (progress(frac, "Warping the DEM") or 1)
    if os.path.exists(out_path):
        gdal.GetDriverByName("GTiff").Delete(out_path)
    with R.quiet():
        out = gdal.Warp(str(out_path), ds, options=gdal.WarpOptions(**opts))
    if out is None:
        raise R.RasterError("gdal.Warp failed to write the prepared DEM.")
    gt = out.GetGeoTransform()
    band = out.GetRasterBand(1)
    try:
        band.ComputeStatistics(False)
    except RuntimeError:
        pass
    a = band.ReadAsArray()
    a = np.where(a == OUT_NODATA, np.nan, a).astype("float32")
    warns, figures = audit(a)
    out.SetMetadata({"FLOODDEPTH_DEM_PREP": f"{len(paths)} input file(s), {dst.GetName()}, {cell:g} m, "
                                            f"{resampling}" + (f", buffer {buffer_m:g} m" if clip_bounds else ""),
                     "FLOODDEPTH_DEM_SOURCE": os.path.basename(str(paths[0]))})
    info = {"output": str(out_path), "inputs": [str(p) for p in paths], "crs": dst.GetName(),
            "cell_m": cell, "native_cell_m": round(native, 2), "cols": out.RasterXSize, "rows": out.RasterYSize,
            "bounds": [gt[0], gt[3] + out.RasterYSize * gt[5], gt[0] + out.RasterXSize * gt[1], gt[3]],
            "valid_share": float(np.isfinite(a).mean()), "resampling": resampling, "notes": notes,
            "warnings": warns, "audit": figures}
    out.FlushCache()
    out = None
    with open(os.path.splitext(str(out_path))[0] + ".json", "w", encoding="utf-8") as f:
        json.dump(info, f, indent=1, default=str)
    return info
