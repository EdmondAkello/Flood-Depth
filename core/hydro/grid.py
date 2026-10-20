# Vendored unchanged (imports adapted) from QEHT 1.0.1 core - QGIS Engineering Hydrology
# Toolkit, Copyright (C) 2026 Edmond Akello, GPL-2.0-or-later, used here under GPL-3.0-or-later.
# Source: https://github.com/EdmondAkello/qeht (commit 5ab0322). Validated there against
# RichDEM / TauDEM (tests/test_flats.py); keep in step with QEHT rather than editing here.
# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Copyright (C) 2026
# Licensed under the GNU General Public License v2 or later.
"""Grid conventions shared by every core module.

NO QGIS IMPORTS IN THIS PACKAGE. Ever. `qeht.core` must remain importable
from a bare Python interpreter with only NumPy present, so the hydrology
can be tested without launching QGIS.

Neighbour convention
--------------------
Eight neighbours are indexed 0-7 internally. All computation uses these
indices; the standard 1/2/4/.../128 D8 codes appear only at the I/O boundary
(see `encode_d8` / `decode_d8`). This avoids log2 lookups in hot loops.

    index : 0=E 1=SE 2=S 3=SW 4=W 5=NW 6=N 7=NE
    D8    : 1    2    4    8    16   32   64   128

which matches the standard clockwise-from-east D8 ordering used across GIS
raster hydrology tools, so encoding is a simple table lookup rather than an
angular calculation.
"""

import numpy as np

# (drow, dcol) for internal indices 0-7. Row increases southward
# (north-up raster), so S = +1 row.
DROW = np.array([0, 1, 1, 1, 0, -1, -1, -1], dtype=np.int64)
DCOL = np.array([1, 1, 0, -1, -1, -1, 0, 1], dtype=np.int64)

# Standard D8 code for each internal index.
D8_CODE = np.array([1, 2, 4, 8, 16, 32, 64, 128], dtype=np.int32)

# True where the neighbour is diagonal (distance factor sqrt(2)).
IS_DIAGONAL = np.array([False, True, False, True, False, True, False, True])

# Sentinel used in receiver arrays for "no downstream cell" (an outlet,
# an unresolvable sink, or a NoData cell).
NO_RECEIVER = -1


def neighbour_distances(cell_width, cell_height):
    """Planimetric distance to each of the 8 neighbours, in map units.

    Handles non-square cells correctly, which matters for DEMs delivered
    in geographic-derived projections where x and y resolution differ
    slightly after reprojection.
    """
    dist = np.empty(8, dtype=np.float64)
    for k in range(8):
        dx = DCOL[k] * cell_width
        dy = DROW[k] * cell_height
        dist[k] = np.hypot(dx, dy)
    return dist


def encode_d8(direction_index, valid=None):
    """Internal 0-7 index array -> standard 1/2/4/.../128 D8 raster.

    Cells with no receiver (NO_RECEIVER) encode to 0, matching the standard
    convention for sinks and cells that flow off the grid edge.
    """
    out = np.zeros(direction_index.shape, dtype=np.int32)
    has = direction_index >= 0
    if valid is not None:
        has &= valid
    out[has] = D8_CODE[direction_index[has]]
    return out


def decode_d8(d8_array):
    """Standard D8-coded raster -> internal 0-7 index array (NO_RECEIVER for 0).

    Lets us ingest a flow direction grid produced by another GIS/hydrology
    toolset and run the rest of our pipeline on it, which is how we
    cross-validate against a trusted production result.
    """
    out = np.full(d8_array.shape, NO_RECEIVER, dtype=np.int64)
    for k in range(8):
        out[d8_array == D8_CODE[k]] = k
    return out


def receivers_from_direction(direction_index, shape):
    """Direction index grid -> flat receiver-id array.

    receiver[i] is the flat index of the cell that i drains to, or
    NO_RECEIVER. Flat indexing keeps the topological traversal in
    accumulation cheap.
    """
    rows, cols = shape
    idx = np.arange(rows * cols, dtype=np.int64).reshape(rows, cols)
    rr, cc = np.divmod(idx, cols)

    receiver = np.full(rows * cols, NO_RECEIVER, dtype=np.int64)
    has = direction_index >= 0
    k = direction_index[has]
    nr = rr[has] + DROW[k]
    nc = cc[has] + DCOL[k]

    inside = (nr >= 0) & (nr < rows) & (nc >= 0) & (nc < cols)
    flat_src = idx[has][inside]
    receiver[flat_src] = nr[inside] * cols + nc[inside]
    return receiver
