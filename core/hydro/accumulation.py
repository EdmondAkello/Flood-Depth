# Vendored unchanged (imports adapted) from QEHT 1.0.1 core - QGIS Engineering Hydrology
# Toolkit, Copyright (C) 2026 Edmond Akello, GPL-2.0-or-later, used here under GPL-3.0-or-later.
# Source: https://github.com/EdmondAkello/qeht (commit 5ab0322). Validated there against
# RichDEM / TauDEM (tests/test_flats.py); keep in step with QEHT rather than editing here.
# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Flow accumulation and Strahler ordering by topological traversal.

O(N) after flow direction is established. No recursion - a recursive
downstream walk overflows the stack on any real DEM, and Python's
recursion limit would be hit long before the memory limit.

Algorithm: count how many upstream cells drain into each cell (indegree),
queue every cell with indegree 0 (a ridge/source cell), then repeatedly
pop a cell, push its load downstream, and decrement the receiver's
indegree - enqueueing it once it reaches 0. Each cell is popped exactly
once and its total is complete when popped, so there is no double
counting.
"""

import numpy as np

from .grid import NO_RECEIVER, receivers_from_direction


def flow_accumulation(direction, valid, weights=None, progress=None):
    """Accumulated upstream contribution per cell.

    Parameters
    ----------
    direction : 2-D internal direction index array (from d8_direction).
    valid : 2-D bool array.
    weights : optional 2-D float array of per-cell contribution. Default
        is 1.0 per cell, giving the standard D8 convention where accumulation is
        a count of upstream cells and EXCLUDES the cell itself. Pass a
        cell-area array for contributing area, or a rainfall/runoff grid
        for weighted accumulation.

    Returns
    -------
    accum : 2-D float array. NaN outside `valid`.
    stats : dict.
    """
    rows, cols = direction.shape
    n = rows * cols
    valid_flat = valid.reshape(-1)

    receiver = receivers_from_direction(direction, (rows, cols))
    # A receiver that is itself invalid is treated as off-grid.
    has_recv = receiver != NO_RECEIVER
    bad = has_recv.copy()
    bad[has_recv] = ~valid_flat[receiver[has_recv]]
    receiver[bad] = NO_RECEIVER

    if weights is None:
        weight_flat = np.ones(n, dtype=np.float64)
    else:
        weight_flat = np.asarray(weights, dtype=np.float64).reshape(-1).copy()
    weight_flat[~valid_flat] = 0.0

    # Standard D8 convention: a cell's accumulation excludes its own weight.
    total = np.zeros(n, dtype=np.float64)

    indegree = np.zeros(n, dtype=np.int32)
    live = valid_flat & (receiver != NO_RECEIVER)
    np.add.at(indegree, receiver[live], 1)

    # Topological order in waves (v0.12, vectorised): every cell whose
    # upstream is complete passes its total on in one step. Identical to the
    # v0.8.3 queue for cell counts; with fractional weights the sums may
    # differ in the last binary digit (different addition order).
    front = np.flatnonzero(valid_flat & (indegree == 0))
    processed = 0
    n_valid = max(1, int(valid_flat.sum()))
    while front.size:
        processed += front.size
        j = receiver[front]
        ok = j != NO_RECEIVER
        src, dst = front[ok], j[ok]
        np.add.at(total, dst, total[src] + weight_flat[src])
        np.subtract.at(indegree, dst, 1)
        dst = np.unique(dst)
        front = dst[indegree[dst] == 0]
        if progress is not None:
            progress(min(0.99, processed / n_valid), "Accumulating flow")

    unresolved = int((indegree > 0).sum())

    accum = np.full(n, np.nan, dtype=np.float64)
    accum[valid_flat] = total[valid_flat]
    accum = accum.reshape(rows, cols)

    stats = {
        "cells_processed": processed,
        "cells_in_cycles": unresolved,   # >0 means the graph has a loop
        "max_accumulation": float(np.nanmax(accum)) if processed else 0.0,
    }
    if progress is not None:
        progress(1.0, "Accumulating flow")
    return accum, stats


def strahler_order(direction, valid, stream_mask, progress=None):
    """Strahler order computed over an EXTRACTED STREAM NETWORK.

    This distinction is deliberate and important. Strahler order computed
    over every cell of a DEM is not the same quantity as Strahler order of
    the channel network, and the two are routinely confused. Orders here
    are defined only on cells where `stream_mask` is True; everything else
    is 0.

    Rule: a cell with no upstream stream cells is order 1. Otherwise, if
    the maximum upstream order occurs more than once, order = max + 1;
    if it occurs exactly once, order = max.
    """
    rows, cols = direction.shape
    n = rows * cols
    valid_flat = valid.reshape(-1)
    stream_flat = np.asarray(stream_mask).reshape(-1) & valid_flat

    receiver = receivers_from_direction(direction, (rows, cols))
    has_recv = receiver != NO_RECEIVER
    bad = has_recv.copy()
    bad[has_recv] = ~stream_flat[receiver[has_recv]]
    receiver = receiver.copy()
    receiver[bad] = NO_RECEIVER

    order = np.zeros(n, dtype=np.int32)
    max_up = np.zeros(n, dtype=np.int32)
    count_at_max = np.zeros(n, dtype=np.int32)

    indegree = np.zeros(n, dtype=np.int32)
    live = stream_flat & (receiver != NO_RECEIVER)
    np.add.at(indegree, receiver[live], 1)

    # Waves (vectorised in v0.12). Strahler order does not depend on the
    # processing order, so the result equals the v0.8.3 queue exactly.
    front = np.flatnonzero(stream_flat & (indegree == 0))
    wave_max = np.zeros(n, dtype=np.int32)     # scratch, reset per wave
    wave_cnt = np.zeros(n, dtype=np.int32)
    while front.size:
        o = np.where(max_up[front] == 0, 1,
                     np.where(count_at_max[front] >= 2, max_up[front] + 1, max_up[front]))
        order[front] = o
        j = receiver[front]
        ok = j != NO_RECEIVER
        dst, oo = j[ok], o[ok]
        if dst.size == 0:
            break
        np.maximum.at(wave_max, dst, oo)
        at = oo == wave_max[dst]
        np.add.at(wave_cnt, dst[at], 1)
        ud = np.unique(dst)
        wm, wc = wave_max[ud].copy(), wave_cnt[ud].copy()
        wave_max[ud] = 0
        wave_cnt[ud] = 0
        higher = wm > max_up[ud]
        equal = wm == max_up[ud]
        count_at_max[ud] = np.where(higher, wc, np.where(equal, count_at_max[ud] + wc,
                                                          count_at_max[ud]))
        max_up[ud] = np.maximum(max_up[ud], wm)
        np.subtract.at(indegree, dst, 1)
        front = ud[indegree[ud] == 0]

    if progress is not None:
        progress(1.0, "Strahler ordering")
    return order.reshape(rows, cols)
