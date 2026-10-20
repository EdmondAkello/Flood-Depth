# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Water-surface elevation from a flood extent, and flood depth.

Method (boundary-elevation method; Cohen et al. 2018, 2019 "FwDET")
------------------------------------------------------------------
1. Waterline cells: flood cells with at least one 4-neighbour that is
   valid DRY ground. A flood cell whose only non-flood neighbours are
   NoData, the grid edge, cloud or other water (a lake, the permanent
   river) is NOT a waterline: there the extent stops for reasons that say
   nothing about the water level. Sampling those edges is the main bias
   of simple boundary methods (FwDET v2.0 handles the domain edge the same
   way).
2. The water-surface elevation at a waterline cell is the terrain
   elevation there (depth ~ 0 at the edge of the water).
3. Outliers: each waterline value is compared with the mean of the
   waterline values around it (a window of `outlier_radius` cells).
   Values more than `outlier_k` robust standard deviations (1.4826 x MAD)
   and more than `outlier_min` (m) away are dropped. This removes edge
   cells that sit on embankments, buildings or tree canopy in a DSM.
   Waterline cells on ground steeper than `steep_slope_pct` are never
   dropped: there the edge elevation changes by metres per cell for real.
4. The kept values are spread over the flood cells. Default "nearest":
   each flood cell takes the elevation of its nearest waterline cell
   (Euclidean allocation, as in FwDET), computed by jump flooding in
   NumPy; then a focal mean over flood cells (radius `smooth_radius`,
   FwDET v2) removes the allocation steps. Nearest allocation keeps a
   steep, narrow valley from borrowing edge elevations from far up- or
   downstream. Smoothing is applied only where it changes the surface by
   less than `smooth_tolerance` (0.5 m): on steep ground (an escarpment,
   a gully) the waterline elevations change by metres per cell and an
   average across them would invent water depth. Option "idw": GDAL's
   inverse-distance fill
   (gdal.FillNodata), smoother on wide flat floodplains.
5. depth = WSE - terrain inside the extent. Cells where the terrain is
   above the estimated surface get depth 0 and are counted.

The same procedure applied to a HAND raster instead of a DEM gives the
water height above the drainage at the edge, and depth = that height -
HAND. That is the correct use of HAND for depth; HAND values inside the
extent are heights above the channel, not water depths.

All steps are whole-array NumPy or GDAL calls: there is no per-cell
Python loop.
"""

import numpy as np

from . import raster as R


def neighbours_any(mask, connectivity=4):
    """True where at least one neighbour (inside the grid) is True."""
    m = np.asarray(mask, bool)
    out = np.zeros_like(m)
    out[1:, :] |= m[:-1, :]
    out[:-1, :] |= m[1:, :]
    out[:, 1:] |= m[:, :-1]
    out[:, :-1] |= m[:, 1:]
    if connectivity == 8:
        out[1:, 1:] |= m[:-1, :-1]
        out[1:, :-1] |= m[:-1, 1:]
        out[:-1, 1:] |= m[1:, :-1]
        out[:-1, :-1] |= m[1:, 1:]
    return out


def waterline(flood, valid, other_water=None, blocked=None):
    """Flood cells that touch valid dry ground (see module notes).

    flood        : cells in the flood extent
    valid        : cells where terrain (and imagery) are valid
    other_water  : water that is not part of the flood (permanent water);
                   counts as water, not as dry ground
    blocked      : cells to treat as unknown (cloud, outside the AOI)
    """
    flood = np.asarray(flood, bool) & valid
    wet = flood.copy()
    if other_water is not None:
        wet |= np.asarray(other_water, bool)
    dry = valid & ~wet
    if blocked is not None:
        dry &= ~np.asarray(blocked, bool)
    return flood & neighbours_any(dry, 4)


def box_sum(a, r):
    """Sum over a (2r+1) x (2r+1) window, edges padded with zeros."""
    a = np.asarray(a, dtype="float64")
    if r <= 0:
        return a.copy()
    k = 2 * r + 1
    c = np.pad(a, ((r + 1, r), (r + 1, r)))
    c = np.cumsum(np.cumsum(c, axis=0), axis=1)
    return c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]


def focal_mean(values, mask, r):
    """Mean of `values` over `mask` cells in a (2r+1)^2 window; NaN if none."""
    v = np.where(mask, np.nan_to_num(values), 0.0)
    n = box_sum(mask.astype("float64"), r)
    s = box_sum(v, r)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = s / n
    out[n <= 0] = np.nan
    return out


def remove_outliers(z, line, radius=15, k=3.0, min_abs=0.5, protect=None):
    """Drop waterline values far from their neighbourhood. Returns (keep, n_dropped, info).

    `protect`: waterline cells that are never dropped (steep ground, where
    edge elevations legitimately change by metres from cell to cell)."""
    keep = line.copy()
    if radius <= 0 or k <= 0 or line.sum() < 10:
        return keep, 0, {}
    local = focal_mean(z, line, radius)
    # exclude the cell itself from its reference when it has neighbours
    n = box_sum(line.astype("float64"), radius)
    s = local * n
    with np.errstate(invalid="ignore", divide="ignore"):
        loo = np.where(n > 1, (s - np.where(line, z, 0.0)) / (n - 1), local)
    res = (z - loo)[line]
    res = res[np.isfinite(res)]
    if res.size < 10:
        return keep, 0, {}
    med = float(np.median(res))
    mad = float(np.median(np.abs(res - med))) * 1.4826
    lim = max(k * mad, float(min_abs))
    bad = line & np.isfinite(loo) & (np.abs(z - loo - med) > lim)
    if protect is not None:
        bad &= ~protect
    keep &= ~bad
    return keep, int(bad.sum()), {"residual_median": med, "robust_sd": mad, "limit": lim}


def proximity_cells(seeds):
    """Distance (cells) to the nearest seed, by GDAL ComputeProximity (C++)."""
    gdal = R._gdal()
    rows, cols = seeds.shape
    src = gdal.GetDriverByName("MEM").Create("", cols, rows, 1, gdal.GDT_Byte)
    src.GetRasterBand(1).WriteArray(seeds.astype("uint8"))
    dst = gdal.GetDriverByName("MEM").Create("", cols, rows, 1, gdal.GDT_Float32)
    gdal.ComputeProximity(src.GetRasterBand(1), dst.GetRasterBand(1), ["VALUES=1", "DISTUNITS=PIXEL"])
    return dst.GetRasterBand(1).ReadAsArray()


def nearest_seed(seeds, cell_x=1.0, cell_y=1.0, max_dist_cells=None):
    """Row/col of the nearest seed cell for every cell (Euclidean).

    Jump flooding (Rong & Tan 2006) with a final 1-step pass: log2(n)
    rounds of 8 shifted comparisons, all whole-array NumPy on the
    overlapping slabs (no full-size temporaries). Verified exact against
    brute force in tests/test_surface.py.
    Returns (seed_row, seed_col, distance) with -1 / inf where no seed.
    """
    seeds = np.asarray(seeds, bool)
    rows, cols = seeds.shape
    sr = np.full(seeds.shape, -1, dtype="int32")
    sc = np.full(seeds.shape, -1, dtype="int32")
    rr, cc = np.nonzero(seeds)
    sr[rr, cc] = rr
    sc[rr, cc] = cc
    best = np.full(seeds.shape, np.inf, dtype="float32")
    best[rr, cc] = 0.0
    if rr.size == 0:
        return sr, sc, best
    gr = np.arange(rows, dtype="float32")[:, None]
    gc = np.arange(cols, dtype="float32")[None, :]
    wy, wx = np.float32(float(cell_y) ** 2), np.float32(float(cell_x) ** 2)
    lim = max(rows, cols) if max_dist_cells is None else min(max(rows, cols), int(max_dist_cells) + 2)
    k = 1
    while k * 2 < lim:
        k *= 2
    steps = []
    while k >= 1:
        steps.append(k)
        k //= 2
    steps.append(1)
    for k in steps:
        for dy in (-k, 0, k):
            for dx in (-k, 0, k):
                if (dy == 0 and dx == 0) or abs(dy) >= rows or abs(dx) >= cols:
                    continue
                # target cells t = (r, c); candidate = seed stored at (r + dy, c + dx)
                ty = slice(max(0, -dy), rows - max(0, dy))
                tx = slice(max(0, -dx), cols - max(0, dx))
                sy = slice(max(0, dy), rows - max(0, -dy))
                sx = slice(max(0, dx), cols - max(0, -dx))
                cr = sr[sy, sx]
                cc_ = sc[sy, sx]
                d = wy * (gr[ty] - cr) ** 2 + wx * (gc[:, tx] - cc_) ** 2
                d[cr < 0] = np.inf
                bt = best[ty, tx]
                better = d < bt
                if better.any():
                    bt[better] = d[better]
                    sr[ty, tx][better] = cr[better]
                    sc[ty, tx][better] = cc_[better]
    return sr, sc, np.sqrt(best)


def fill_nearest(values, known, cell_x=1.0, cell_y=1.0, targets=None):
    """Every cell takes the value of its nearest known cell (Euclidean allocation).

    targets: cells that must be exact (the flood); the jump schedule starts
    just above their largest distance to a known cell (ComputeProximity),
    so cells farther away than that may be approximate - they are not used.
    """
    md = None
    if targets is not None and targets.any():
        aniso = max(cell_x, cell_y) / min(cell_x, cell_y)
        md = float(proximity_cells(known)[targets].max()) * aniso
    sr, sc, dist = nearest_seed(known, cell_x, cell_y, md)
    out = np.full(values.shape, np.nan, dtype="float32")
    ok = sr >= 0
    out[ok] = np.asarray(values, dtype="float32")[sr[ok], sc[ok]]
    return out, dist


def fill_from_points(values, known, grid, max_search=None, smoothing_iterations=0):
    """Inverse-distance fill of a grid from known cells (gdal.FillNodata)."""
    gdal = R._gdal()
    rows, cols = values.shape
    mem = gdal.GetDriverByName("MEM").Create("", cols, rows, 1, gdal.GDT_Float32)
    mem.SetGeoTransform(grid.gt)
    mem.SetProjection(grid.wkt)
    b = mem.GetRasterBand(1)
    nd = -3.0e38
    b.SetNoDataValue(nd)
    b.WriteArray(np.where(known, values, nd).astype("float32"))
    mk = gdal.GetDriverByName("MEM").Create("", cols, rows, 1, gdal.GDT_Byte)
    mk.GetRasterBand(1).WriteArray(known.astype("uint8"))
    if max_search is None:
        max_search = float(rows + cols)
    gdal.FillNodata(b, mk.GetRasterBand(1), maxSearchDist=float(max_search),
                    smoothingIterations=int(smoothing_iterations))
    out = b.ReadAsArray().astype("float32")
    out[out <= -1e38] = np.nan
    return out


def estimate_surface(elev, flood, valid, grid, other_water=None, blocked=None,
                     outlier_radius=15, outlier_k=3.0, outlier_min=0.5,
                     smooth_radius=2, smooth_tolerance=0.5, interpolation="nearest",
                     steep_slope_pct=10.0, max_search_cells=None, progress=None):
    """Water-surface (or HAND water-height) raster over the flood cells.

    Returns (surface float32 with NaN outside the flood, info dict, waterline
    mask of the cells used, waterline mask of the cells dropped).
    """
    def step(f, msg):
        if progress:
            progress(f, msg)

    elev = np.asarray(elev, dtype="float32")
    valid = np.asarray(valid, bool) & np.isfinite(elev)
    flood = np.asarray(flood, bool) & valid
    info = {"flood_cells": int(flood.sum())}
    if info["flood_cells"] == 0:
        raise ValueError("The flood extent has no cells with valid terrain on the analysis grid.")
    step(0.05, "Finding the waterline")
    line = waterline(flood, valid, other_water, blocked)
    info["waterline_cells"] = int(line.sum())
    border = np.zeros_like(valid)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    edge_like = flood & (border | neighbours_any(~valid | (blocked if blocked is not None else np.zeros_like(valid)), 4))
    info["edge_cells_excluded"] = int((edge_like & ~line).sum())
    if line.sum() < 3:
        raise ValueError(
            f"Only {int(line.sum())} waterline cell(s) were found: the extent does not border dry ground "
            "with valid terrain. Check that the DEM covers the area around the flood and that the extent "
            "is not clipped to the DEM or image edge.")
    step(0.15, "Screening waterline outliers")
    slope = slope_percent(elev, grid.cell_x, grid.cell_y)
    steep = line & (slope > float(steep_slope_pct))
    info["waterline_steep"] = int(steep.sum())
    info["flood_steep_share"] = float((flood & (slope > float(steep_slope_pct))).sum() / max(1, flood.sum()))
    keep, dropped, oinfo = remove_outliers(elev, line, outlier_radius, outlier_k, outlier_min, protect=steep)
    info.update({"waterline_used": int(keep.sum()), "waterline_dropped": dropped})
    info.update({"outlier_" + k: v for k, v in oinfo.items()})
    if keep.sum() < 3:
        keep = line
        info["waterline_used"] = int(line.sum())
        info["waterline_dropped"] = 0
    step(0.3, "Interpolating the water surface")
    # work inside the flood's bounding box (+ margin): the waterline lies in it
    rr, cc = np.nonzero(flood)
    m = int(max(2, smooth_radius or 0))
    r0, r1 = max(0, rr.min() - m), min(flood.shape[0], rr.max() + m + 1)
    c0, c1 = max(0, cc.min() - m), min(flood.shape[1], cc.max() + m + 1)
    win = (slice(r0, r1), slice(c0, c1))
    e_w, k_w, f_w = elev[win], keep[win], flood[win]
    if interpolation == "idw":
        sub = R.Grid(c1 - c0, r1 - r0, (grid.gt[0] + c0 * grid.gt[1], grid.gt[1], 0.0,
                                         grid.gt[3] + r0 * grid.gt[5], 0.0, grid.gt[5]), grid.wkt)
        s_w = fill_from_points(e_w, k_w, sub, max_search=max_search_cells)
    else:
        s_w, dist = fill_nearest(e_w, k_w, grid.cell_x, grid.cell_y, targets=f_w)
        info["max_distance_to_waterline_m"] = float(np.max(dist[f_w])) if f_w.any() else 0.0
    info["interpolation"] = interpolation
    if smooth_radius and smooth_radius > 0:
        step(0.8, "Smoothing")
        sm = focal_mean(s_w, f_w & np.isfinite(s_w), int(smooth_radius))
        ok = np.isfinite(sm) & (np.abs(sm - s_w) <= float(smooth_tolerance))
        info["smoothing_skipped_cells"] = int((f_w & np.isfinite(sm) & ~ok).sum())
        s_w = np.where(ok, sm, s_w).astype("float32")
    surf = np.full(flood.shape, np.nan, dtype="float32")
    surf[win] = s_w
    surf = np.where(flood, surf, np.nan).astype("float32")
    info["unfilled_cells"] = int((flood & ~np.isfinite(surf)).sum())
    step(1.0, "")
    return surf, info, keep, line & ~keep


def slope_percent(elev, cell_x, cell_y):
    z = np.asarray(elev, dtype="float64")
    gy, gx = np.gradient(np.where(np.isfinite(z), z, np.nanmean(z)), cell_y, cell_x)
    return (100.0 * np.hypot(gx, gy)).astype("float32")


def depth_from_surface(surface, elev, flood):
    """(depth with NaN outside the flood, number of cells clipped to 0)."""
    d = np.asarray(surface, dtype="float32") - np.asarray(elev, dtype="float32")
    inside = np.asarray(flood, bool) & np.isfinite(d)
    neg = inside & (d < 0)
    d = np.where(inside, np.maximum(d, 0.0), np.nan).astype("float32")
    return d, int(neg.sum())


def depth_constant_stage(level, elev, flood=None, valid=None):
    """Bathtub depth below one water level; restricted to `flood` if given."""
    e = np.asarray(elev, dtype="float32")
    v = np.isfinite(e) if valid is None else (np.asarray(valid, bool) & np.isfinite(e))
    d = float(level) - e
    inside = v & (d > 0)
    if flood is not None:
        inside &= np.asarray(flood, bool)
    return np.where(inside, d, np.nan).astype("float32")
