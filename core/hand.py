# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Height Above Nearest Drainage (HAND) from a DEM.

HAND (Renno et al. 2008; Nobre et al. 2011) is the height of a cell above
the drainage cell its flow path reaches. Computed here with the routing
core of QEHT (core/hydro, vendored): priority-flood depression filling
(Barnes et al. 2014), D8 flow direction with iterated Barnes flat
resolution, flow accumulation. Streams are cells whose contributing area
reaches `stream_km2`. The drainage elevation is carried from every stream
cell up its contributing area in reverse topological order (whole-array
waves, no per-cell Python loop):

    drain_z(stream cell) = filled z
    drain_z(cell)        = drain_z(receiver)      otherwise
    HAND                 = z - drain_z            (raw DEM z; >= 0 except in pits)

Cells whose flow path leaves the grid or ends in an unresolved sink before
reaching a stream take the elevation of that outlet. Near the grid edge,
where part of the catchment is cut off, HAND is less reliable: give the
DEM a margin around the flood.
"""

import numpy as np

from .hydro.fill import fill_depressions
from .hydro.direction import d8_direction
from .hydro.grid import NO_RECEIVER, receivers_from_direction


def compute_hand(dem, valid, cell_x, cell_y, stream_km2=1.0, progress=None):
    """(hand, info) for a DEM array (NaN where invalid)."""
    def step(f, m):
        if progress:
            progress(f, m)

    z = np.asarray(dem, dtype="float64")
    valid = np.asarray(valid, bool) & np.isfinite(z)
    rows, cols = z.shape
    step(0.0, "Filling depressions (HAND)")
    filled, n_filled, _ = fill_depressions(np.where(valid, z, 0.0), valid, cell_width=cell_x, cell_height=cell_y)
    step(0.35, "Flow direction (HAND)")
    direction, dstats = d8_direction(np.where(valid, filled, np.nan), valid, cell_x, cell_y)
    step(0.6, "Flow accumulation (HAND)")
    receiver = receivers_from_direction(direction, (rows, cols))
    vf = valid.reshape(-1)
    has = receiver != NO_RECEIVER
    bad = has.copy()
    bad[has] = ~vf[receiver[has]]
    receiver[bad] = NO_RECEIVER

    # forward topological waves (upstream -> downstream), recording the order
    n = rows * cols
    indeg = np.zeros(n, dtype=np.int32)
    live = vf & (receiver != NO_RECEIVER)
    np.add.at(indeg, receiver[live], 1)
    acc = np.zeros(n, dtype=np.float64)
    front = np.flatnonzero(vf & (indeg == 0))
    waves = []
    while front.size:
        waves.append(front)
        j = receiver[front]
        ok = j != NO_RECEIVER
        np.add.at(acc, j[ok], acc[front[ok]] + 1.0)
        np.subtract.at(indeg, j[ok], 1)
        d = np.unique(j[ok])
        front = d[indeg[d] == 0]
    area_km2 = (acc + 1.0) * cell_x * cell_y / 1e6
    stream = vf & (area_km2 >= float(stream_km2))

    step(0.8, "Drainage elevation (HAND)")
    ff = filled.reshape(-1)
    drain = np.full(n, np.nan)
    outlet = vf & (receiver == NO_RECEIVER)
    for w in reversed(waves):            # downstream first
        src = np.where(stream[w] | outlet[w], ff[w], drain[np.where(receiver[w] >= 0, receiver[w], 0)])
        drain[w] = src
    hand = (z.reshape(-1) - drain).reshape(rows, cols)
    hand[~valid] = np.nan
    info = {"stream_km2": float(stream_km2), "stream_cells": int(stream.sum()), "filled_cells": int(n_filled),
            "unrouted_cells": int(dstats.get("cells_still_unrouted", 0)), "waves": len(waves)}
    step(1.0, "")
    return hand.astype("float32"), info
