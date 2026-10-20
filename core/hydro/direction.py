# Vendored unchanged (imports adapted) from QEHT 1.0.1 core - QGIS Engineering Hydrology
# Toolkit, Copyright (C) 2026 Edmond Akello, GPL-2.0-or-later, used here under GPL-3.0-or-later.
# Source: https://github.com/EdmondAkello/qeht (commit 5ab0322). Validated there against
# RichDEM / TauDEM (tests/test_flats.py); keep in step with QEHT rather than editing here.
# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""D8 flow direction - steepest descent weighted by distance.

This is the module that most needed writing from scratch. The reference
plugin we studied assigns each cell's receiver during a priority-flood
traversal, which means:

  * diagonal and orthogonal neighbours are weighted identically, and
  * the receiver is the traversal parent, not the steepest descent.

Standard D8 (as used by commercial GIS/hydrology toolsets) selects the neighbour maximising

    drop / distance     where distance = cellsize      (orthogonal)
                                       = cellsize*sqrt(2)  (diagonal)

Getting this wrong biases flow paths toward diagonals and changes both
accumulation values and catchment boundaries. We implement the weighted
form.

Flat handling
-------------
After filling, flats are common and have no steepest descent. We resolve
them with a gradient-toward-lower-terrain BFS: cells adjacent to a flat
outlet get distance 1, their flat neighbours 2, and so on; each flat cell
then routes to the flat neighbour with the smaller BFS distance.

This is the "gradient towards lower terrain" half of Garbrecht & Martz
(1997). We do NOT implement the "gradient away from higher terrain" half
in v0.1. Consequence: drainage across wide flats is plausible and fully
connected, but converges to fewer channels than commercial reference
toolsets typically do. This is a
KNOWN AND DOCUMENTED DIVERGENCE - it must be stated in any validation
report comparing our output against a reference hydrology toolset's result over flat terrain.
Avoid this by using min_slope > 0 during filling, which largely removes
flats before routing.
"""

import numpy as np
from collections import deque

from .grid import DROW, DCOL, IS_DIAGONAL, NO_RECEIVER, neighbour_distances

# Standard D8 tie-breaking priority. When two or more neighbours give an
# equal steepest drop/distance, a naive scan-order pick does not match what
# commercial reference toolsets produce - they apply a fixed preference.
# Recovered empirically from a reference hydrology toolset's flow
# direction grid over low-relief coastal terrain (21% of cells tied), the
# order among the cardinals is S > W > N > E, with diagonals ranked below
# all cardinals. Applying it lifted off-flat agreement with the reference
# toolset from ~84% to ~91% across four independent windows. Earlier index = wins.
#   index: 0=E 1=SE 2=S 3=SW 4=W 5=NW 6=N 7=NE
_TIE_PRIORITY = [2, 4, 6, 0, 1, 3, 5, 7]   # S, W, N, E, SE, SW, NW, NE
_TIE_RANK = [0] * 8
for _p, _k in enumerate(_TIE_PRIORITY):
    _TIE_RANK[_k] = _p
_TIE_RANK = np.array(_TIE_RANK, dtype=np.int64)

# Relative tolerance for treating two slopes as tied. Low-relief DEMs
# produce many near-equal drops that are "the same" within floating point
# and within DEM vertical precision; 1e-6 captures these without merging
# genuinely different gradients.
_TIE_TOL = 1e-6


def d8_direction(elevation, valid, cell_width=1.0, cell_height=1.0,
                 resolve_flats=True, flat_method="barnes", progress=None):
    """Compute D8 direction indices (0-7) from a conditioned DEM.

    flat_method : "barnes" (default since v0.13; Barnes, Lehman & Mulla
        2014, identical to Barnes' own RichDEM implementation) or "toward"
        (shortest path to the flat's outlet; the default up to v0.12, kept
        for matching outputs of reference platforms that route flats that
        way). The WP-G benchmark (v0.13, 80 Kenyan AOIs x ALOS/FABDEM)
        chose Barnes; the v0.12 hybrid family was removed from this API
        because it reproduced Barnes (core/flow/flats.resolve_flats keeps w).

    Returns
    -------
    direction : 2-D int64 array of internal indices, NO_RECEIVER where the
        cell is NoData, is an outlet draining off-grid, or is an
        unresolvable sink.
    stats : dict with diagnostic counts.

    Ties (neighbours within a relative 1e-6 of the steepest drop/distance)
    are broken by the fixed priority S, W, N, E, SE, SW, NW, NE - the rule
    recovered empirically from a reference GIS platform (off-flat agreement
    84 % -> 91 %). Deterministic and reproducible.
    """
    elev = np.asarray(elevation, dtype=np.float64)
    rows, cols = elev.shape
    dist = neighbour_distances(cell_width, cell_height)

    work = np.where(valid, elev, np.inf)

    # Streamed over the 8 directions (v0.12): no rows x cols x 8 arrays, so
    # peak memory is a few copies of the grid instead of ~17. Identical to
    # the v0.8.3 two-array version (kept in core/flow/_reference.py).
    def _slope_k(k):
        dr, dc = int(DROW[k]), int(DCOL[k])
        r0_s, r1_s = max(0, -dr), rows - max(0, dr)
        c0_s, c1_s = max(0, -dc), cols - max(0, dc)
        r0_n, r1_n = max(0, dr), rows - max(0, -dr)
        c0_n, c1_n = max(0, dc), cols - max(0, -dc)
        with np.errstate(invalid="ignore"):
            sl = (work[r0_s:r1_s, c0_s:c1_s] - work[r0_n:r1_n, c0_n:c1_n]) / dist[k]
        return (slice(r0_s, r1_s), slice(c0_s, c1_s)), np.where(np.isfinite(sl), sl, -np.inf)

    # Pass 1: the steepest drop/distance available at each cell.
    best_slope = np.full((rows, cols), -np.inf, dtype=np.float64)
    for k in range(8):
        win, sl = _slope_k(k)
        np.maximum(best_slope[win], sl, out=best_slope[win])
        if progress is not None:
            progress((k + 1) / 8.0 * 0.3, "Computing D8 flow direction")

    # Pass 2: among neighbours tied within tolerance at the steepest drop,
    # pick by the fixed tie priority S, W, N, E, SE, SW, NW, NE (recovered
    # empirically from a reference GIS platform), not by scan order.
    thresh = best_slope * (1.0 - _TIE_TOL)
    best_rank = np.full((rows, cols), 99, dtype=np.int8)
    direction = np.zeros((rows, cols), dtype=np.int64)
    for k in range(8):
        win, sl = _slope_k(k)
        with np.errstate(invalid="ignore"):
            take = (sl >= thresh[win]) & (_TIE_RANK[k] < best_rank[win])
        best_rank[win][take] = _TIE_RANK[k]
        direction[win][take] = k
        if progress is not None:
            progress(0.3 + (k + 1) / 8.0 * 0.3, "Computing D8 flow direction")

    # Cells with no downhill neighbour: no positive best slope.
    no_flow = ~(valid & np.isfinite(best_slope) & (best_slope > 0))
    direction[no_flow] = NO_RECEIVER
    direction[~valid] = NO_RECEIVER
    if progress is not None:
        progress(0.6, "Computing D8 flow direction")

    n_no_descent = int((valid & (direction < 0)).sum())
    flat_stats = {}

    # Boundary cells (grid edge or next to NoData) with no descent are
    # outlets that drain off the grid. v0.13: flats touching them drain to
    # them (TauDEM and RichDEM do the same); up to v0.12 such flats stayed
    # unrouted and became false sinks at clip edges and coastlines.
    from .fill import valid_boundary_mask
    drain = valid & (direction < 0) & valid_boundary_mask(valid)

    if resolve_flats and n_no_descent:
        # Garbrecht & Martz: impose a two-gradient surface on the flats,
        # then re-run the ordinary steepest-descent pass over it. One
        # routing rule, applied twice - not two competing rules.
        if flat_method not in ("barnes", "toward"):
            raise ValueError(f"flat_method must be 'barnes' or 'toward', not {flat_method!r}")
        if flat_method == "barnes":
            # Iterated Barnes 2014 convergent flat resolution. Repeats the
            # Barnes pass until no further flats resolve, so flats drain in
            # hierarchy order (a flat whose outlet is a lower flat resolves
            # once that lower flat has drained).
            from .flats import resolve_flats_barnes
            if progress is not None:
                progress(0.75, "Resolving flats (Barnes 2014, iterated)")
            direction, flat_stats = resolve_flats_barnes(elev, valid, direction, drain=drain)
        else:
            direction, n_res = _resolve_flats_toward(work, valid, direction, drain=drain)
            flat_stats = {"flat_cells": n_res, "method": "toward-lower only"}

    stats = {
        "cells_valid": int(valid.sum()),
        "cells_without_descent": n_no_descent,
        "flat_method": flat_stats.get("method", "none"),
        "flat_iterations": flat_stats.get("iterations", 1),
        "flat_cells": flat_stats.get("flat_cells", 0),
        "flat_outlet_seeds": flat_stats.get("flat_outlet_seeds", 0),
        "flat_high_edge_seeds": flat_stats.get("flat_high_edge_seeds", 0),
        "boundary_outlets": int((drain & (direction < 0)).sum()),
        "cells_still_unrouted": int((valid & (direction < 0) & ~drain).sum()),
    }
    return direction, stats


def _resolve_flats_toward(work, valid, direction, progress=None, drain=None):
    """Route flat cells by BFS distance toward the nearest flat outlet.

    Vectorised in v0.12 (whole-frontier BFS, neighbour shifts). Reproduces
    the v0.8.3 per-cell algorithm (core/flow/_reference.py) exactly:
      seeds  : flat cells with a routed neighbour at equal or lower elevation;
      BFS    : 1 at the seeds, +1 per step across equal-elevation flat cells;
      assign : each flat cell drains to the not-higher neighbour with the
               smallest distance (a routed non-flat neighbour counts as 0),
               first in E, SE, S ... order on ties.
    """
    from .flats import _neighbours, _csr, _bfs
    rows, cols = work.shape
    drain = np.zeros(work.shape, dtype=bool) if drain is None else drain
    routed = (direction >= 0) | drain
    flat = valid & (direction < 0) & ~drain
    if not flat.any():
        return direction, 0
    fr, fc = np.nonzero(flat)
    n = fr.size
    idx = np.full((rows, cols), -1, dtype=np.int64)
    idx[fr, fc] = np.arange(n)
    wv = work[fr, fc]
    seed = np.zeros(n, dtype=bool)
    ea, eb = [], []
    for k in range(8):
        nr, nc, inside = _neighbours(fr, fc, k, rows, cols)
        ok = inside & valid[nr, nc]
        seed |= ok & routed[nr, nc] & (work[nr, nc] <= wv)
        j = idx[nr, nc]
        link = ok & (j >= 0) & (work[nr, nc] == wv)
        ea.append(np.flatnonzero(link)); eb.append(j[link])
    indptr, indices = _csr(n, np.concatenate(ea), np.concatenate(eb))
    dist = _bfs(n, indptr, indices, np.flatnonzero(seed))
    INF = np.iinfo(np.int32).max
    bfs = np.where(dist > 0, dist, INF)
    bfs_grid = np.full((rows, cols), INF, dtype=np.int64)
    bfs_grid[fr, fc] = bfs

    best_d = np.full(n, INF, dtype=np.int64)
    best_k = np.full(n, -1, dtype=np.int64)
    for k in range(8):
        nr, nc, inside = _neighbours(fr, fc, k, rows, cols)
        ok = inside & valid[nr, nc] & (work[nr, nc] <= wv)
        nd = np.where(routed[nr, nc] & ~flat[nr, nc], 0, bfs_grid[nr, nc])
        better = ok & (nd < best_d)
        best_d = np.where(better, nd, best_d)
        best_k = np.where(better, k, best_k)
    good = (best_k >= 0) & (best_d < INF)
    direction[fr[good], fc[good]] = best_k[good]
    if progress is not None:
        progress(1.0, "Resolving flats")
    return direction, int(good.sum())
