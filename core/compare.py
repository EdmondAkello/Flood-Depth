# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Compare a flood extent with a reference extent (validation)."""

import math

import numpy as np

from . import raster as R
from . import vector as V
from . import stats as ST

AGREEMENT = {0: "dry in both", 1: "flooded in both (hit)", 2: "false alarm", 3: "miss"}
AGREEMENT_COLORS = {0: (240, 240, 240, 255), 1: (33, 113, 181, 255), 2: (230, 85, 13, 255),
                    3: (117, 107, 177, 255), 255: (0, 0, 0, 0)}


def _bounds(src):
    if "raster" in src:
        g = R.grid_of(R.open_dataset(src["raster"]))
        return g.bounds, g.wkt
    from .pipeline import _wkts_bounds
    return _wkts_bounds(src["wkts"], src["crs"])


def _mask(src, grid):
    """(water, valid) of a raster (value == water value) or polygons on grid."""
    if "raster" in src:
        a, _ = R.read_on_grid(src["raster"], grid, src.get("band", 1), "near")
        valid = np.isfinite(a) & (a != 255) & (a >= 0)
        return valid & (a == float(src.get("value", 1))), valid
    m = V.rasterize_wkts(src["wkts"], src["crs"], grid)
    return m, np.ones(grid.shape, bool)


def compare_extents(pred, ref, aoi_wkts=None, aoi_crs="", cell=0.0):
    """pred / ref: {'raster': path, 'value': 1} or {'wkts': [...], 'crs': wkt}.

    Grid: the reference raster's grid, else the predicted raster's, else a
    WGS 84 / UTM grid of `cell` m (default 10) over both extents.
    Returns (scores dict, agreement uint8 array, grid)."""
    if "raster" in ref and not cell:
        grid = R.grid_of(R.open_dataset(ref["raster"]))
    elif "raster" in pred and not cell:
        grid = R.grid_of(R.open_dataset(pred["raster"]))
    else:
        b1, w1 = _bounds(pred)
        b2, w2 = _bounds(ref)
        b2 = R.transform_bounds(b2, w2, w1)
        b = (min(b1[0], b2[0]), min(b1[1], b2[1]), max(b1[2], b2[2]), max(b1[3], b2[3]))
        ll = R.transform_bounds(b, w1, R.srs_from_epsg(4326).ExportToWkt())
        epsg = R.utm_epsg((ll[0] + ll[2]) / 2, (ll[1] + ll[3]) / 2)
        wkt = R.srs_from_epsg(epsg).ExportToWkt()
        bb = R.transform_bounds(b, w1, wkt)
        cs = float(cell) if cell else 10.0
        x0, y1 = math.floor(bb[0] / cs) * cs, math.ceil(bb[3] / cs) * cs
        grid = R.Grid(int(math.ceil((bb[2] - x0) / cs)), int(math.ceil((y1 - bb[1]) / cs)),
                      (x0, cs, 0.0, y1, 0.0, -cs), wkt)
    p, pv = _mask(pred, grid)
    r, rv = _mask(ref, grid)
    valid = pv & rv
    if aoi_wkts:
        valid &= V.rasterize_wkts(aoi_wkts, aoi_crs, grid)
    s = ST.extent_agreement(p, r, valid)
    ca = grid.cell_area() if not R.is_geographic(grid.wkt) else None
    if ca:
        s.update({"pred_area_km2": (s["tp"] + s["fp"]) * ca / 1e6, "ref_area_km2": (s["tp"] + s["fn"]) * ca / 1e6})
    ag = np.full(grid.shape, 255, dtype="uint8")
    ag[valid] = 0
    ag[valid & p & r] = 1
    ag[valid & p & ~r] = 2
    ag[valid & ~p & r] = 3
    return s, ag, grid
