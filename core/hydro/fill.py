# Vendored unchanged (imports adapted) from QEHT 1.0.1 core - QGIS Engineering Hydrology
# Toolkit, Copyright (C) 2026 Edmond Akello, GPL-2.0-or-later, used here under GPL-3.0-or-later.
# Source: https://github.com/EdmondAkello/qeht (commit 5ab0322). Validated there against
# RichDEM / TauDEM (tests/test_flats.py); keep in step with QEHT rather than editing here.
# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Depression filling by priority-flood.

Deliberately separated from D8 routing. DDM HydroLogic fuses the two -
its priority-flood pass emits the receiver graph directly, so the flow
direction it produces is the traversal parent rather than the steepest
descent neighbour. That is fine for connectivity tracing but it is not
compatible with standard D8 routing and it ignores diagonal distance entirely.

Here, conditioning produces ONLY a conditioned elevation surface. Routing
is a separate decision made in `qeht.core.flow.direction`. That separation
is what lets a user document which conditioning method was applied - which
is the difference between a defensible drainage design and a black box.

Reference: Barnes, Lehman & Mulla (2014), "Priority-flood: an optimal
depression-filling and watershed-labeling algorithm for digital elevation
models", Computers & Geosciences 62:117-127. O(n log n).
"""

import numpy as np

from .grid import DROW, DCOL


def valid_boundary_mask(valid):
    """Valid cells that touch the grid edge or touch a NoData cell.

    A clipped DEM usually arrives with a NoData collar, so the outer row
    and column are not the hydrological boundary. Seeding only the raster
    edge would force the whole grid to drain to one global low point.
    """
    rows, cols = valid.shape
    boundary = np.zeros_like(valid, dtype=bool)

    boundary[0, :] |= valid[0, :]
    boundary[-1, :] |= valid[-1, :]
    boundary[:, 0] |= valid[:, 0]
    boundary[:, -1] |= valid[:, -1]

    # Valid cells orthogonally or diagonally adjacent to NoData.
    invalid = ~valid
    padded = np.zeros((rows + 2, cols + 2), dtype=bool)
    padded[1:-1, 1:-1] = invalid
    touches_nodata = np.zeros_like(valid, dtype=bool)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            touches_nodata |= padded[1 + dr:rows + 1 + dr, 1 + dc:cols + 1 + dc]
    boundary |= valid & touches_nodata
    return boundary


def fill_depressions(dem, valid, min_slope=0.0, cell_width=1.0, cell_height=1.0,
                     progress=None):
    """Fill depressions, returning (filled_dem, n_filled_cells, seed_count).

    Parameters
    ----------
    dem : 2-D float array, elevations. Values under ~valid are ignored.
    valid : 2-D bool array, True where the cell carries real elevation.
    min_slope : optional gradient (m per m) imposed across filled flats.
        0.0 gives a classic flat fill (matches common `Fill` tool behaviour). A small
        positive value (1e-4 is typical) tilts filled surfaces so that
        downstream routing has a defined direction instead of relying on
        flat-resolution. Applied as an epsilon increment per traversal
        step scaled by cell size.
    progress : optional callable(fraction_0_to_1, message).

    Returns
    -------
    filled : 2-D float array. Cells outside `valid` are set to NaN.
    n_filled : int, count of cells whose elevation was raised.
    seed_count : int, number of boundary outlet cells seeded.

    Method (v0.12)
    --------------
    Priority-flood (Barnes et al. 2014a) - with or without epsilon -
    produces the unique solution of

        F = z                                   on the boundary seeds
        F = max(z, min over valid neighbours of F + epsilon)   elsewhere

    (the lowest "spill" surface; with epsilon > 0 every filled step rises by
    epsilon). v0.12 computes that solution by frontier relaxation, fully
    vectorised: start from +inf inside, and repeatedly lower the neighbours
    of the cells that changed, until nothing changes. The result is
    identical, bit for bit, to the v0.8.3 heap implementation (kept in
    core/flow/_reference.py as the test oracle) and runs in a small
    fraction of the time.
    """
    dem = np.asarray(dem, dtype=np.float64)
    rows, cols = dem.shape
    mean_cell = 0.5 * (abs(cell_width) + abs(cell_height))
    epsilon = float(min_slope) * mean_cell

    boundary = valid_boundary_mask(valid)
    if not boundary.any():
        # Degenerate case: fully interior valid region with no edge and no
        # NoData contact. Seed the single lowest valid cell.
        ids = np.flatnonzero(valid.reshape(-1))
        if ids.size == 0:
            return np.full((rows, cols), np.nan), 0, 0
        boundary = np.zeros_like(valid)
        boundary.reshape(-1)[ids[np.argmin(dem.reshape(-1)[ids])]] = True
    seed_count = int(boundary.sum())

    # Pad by one cell so neighbour shifts need no bounds checks.
    P = np.full((rows + 2, cols + 2), np.inf)
    Z = np.full((rows + 2, cols + 2), np.inf)
    Z[1:-1, 1:-1] = np.where(valid, dem, np.inf)
    V = np.zeros((rows + 2, cols + 2), dtype=bool)
    V[1:-1, 1:-1] = valid
    S = np.zeros((rows + 2, cols + 2), dtype=bool)
    S[1:-1, 1:-1] = boundary
    P[S] = Z[S]
    W = cols + 2
    offs = np.array([int(DROW[k]) * W + int(DCOL[k]) for k in range(8)], dtype=np.int64)
    Pf, Zf, Vf, Sf = P.reshape(-1), Z.reshape(-1), V.reshape(-1), S.reshape(-1)
    updatable = Vf & ~Sf

    front = np.flatnonzero(Sf)
    slot = np.zeros(Pf.size, dtype=np.int64)       # sort-free de-duplication
    it = 0
    while front.size:
        it += 1
        cand = (front[:, None] + offs[None, :]).reshape(-1)
        cand = cand[updatable[cand]]
        pos = np.arange(cand.size)
        slot[cand] = pos
        cand = cand[slot[cand] == pos]             # keep one copy of each cell
        if cand.size == 0:
            break
        nb = cand[:, None] + offs[None, :]
        m = Pf[nb].min(axis=1)
        new = np.maximum(Zf[cand], m + epsilon)
        ch = new < Pf[cand]
        Pf[cand[ch]] = new[ch]
        front = cand[ch]
        if progress is not None and it % 200 == 0:
            progress(min(0.99, it / float(rows + cols)), "Filling depressions")

    filled = P[1:-1, 1:-1].copy()
    filled[~valid] = np.nan
    n_filled = int((valid & (filled > dem)).sum())
    if progress is not None:
        progress(1.0, "Filling depressions")
    return filled, n_filled, seed_count


def depression_depth(dem, filled, valid):
    """Per-cell fill depth - a QA/QC product, not an intermediate.

    Large contiguous depths in a satellite DEM usually mark a road
    embankment, a dam wall, or a culvert the DEM does not know about.
    Reviewing this raster before accepting the conditioning is the step
    that catches a catchment draining the wrong side of a road.
    """
    depth = np.full(dem.shape, np.nan, dtype=np.float64)
    depth[valid] = filled[valid] - dem[valid]
    depth[valid & (depth < 0)] = 0.0
    return depth
