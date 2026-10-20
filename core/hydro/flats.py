# Vendored unchanged (imports adapted) from QEHT 1.0.1 core - QGIS Engineering Hydrology
# Toolkit, Copyright (C) 2026 Edmond Akello, GPL-2.0-or-later, used here under GPL-3.0-or-later.
# Source: https://github.com/EdmondAkello/qeht (commit 5ab0322). Validated there against
# RichDEM / TauDEM (tests/test_flats.py); keep in step with QEHT rather than editing here.
# -*- coding: utf-8 -*-
# QEHT - QGIS Engineering Hydrology Toolkit
# Licensed under the GNU General Public License v2 or later.
"""Flat resolution: Barnes (2014), the QEHT default since v0.13.

Barnes, Lehman & Mulla (2014) resolve a flat with two BFS gradients:
    toward : distance from the flat's low edges (its outlets)
    away   : distance from its high edges (adjacent higher terrain)
combined as   flat_mask = w * toward + (flat_height - away),  w = 2.
Each flat cell flows to the same-flat neighbour with the lowest flat_mask.

v0.13: this reproduces Barnes' own implementation (RichDEM) cell for cell -
low-edge cells leave the flat directly and neighbours are scanned W, NW,
N, NE, E, SE, S, SW (oracle: _reference.reference_barnes_richdem). Flats
touching boundary drains (grid edge / NoData) drain to them.

The weight w is kept for research; the Processing tools use w = 2 (the
v0.12 'hybrid' option was removed after the WP-G benchmark showed it
reproduces Barnes). Properties of w:
  * w = 2        exactly Barnes 2014;
  * w -> inf     the away term only breaks ties: flow follows a shortest
                 path to the outlet (tested), like toward-lower - but among
                 EQUALLY short paths the Barnes rules (same flat, cardinal
                 first) choose differently from toward-lower (first in E,
                 SE, S ... order). On broad real flats that choice dominates:
                 on a 30 m coastal clip every w from 1.5 to 1e6 gave nearly
                 the same network (stream IoU vs toward-lower 0.394-0.399,
                 ~44 % of flat cells routed differently). w is therefore a
                 weak control there; the tie rule is the strong one;
  * 1 < w < 2    stronger push off the high edges.
Why any w > 1 is safe: a cell only ever flows to a neighbour with a
STRICTLY lower flat_mask, so no cycle is possible for any w. Every cell
with toward = t > 1 has a neighbour with toward = t - 1 whose flat_mask is
lower by at least w - 1 > 0 (away changes by at most 1 per step), so only
low-edge cells (toward = 1) can be local minima, and those drain out of the
flat by the exit rule. For w <= 1 interior false sinks become possible.

Size switch (optional speed setting): flats with fewer than
`small_flat_cells` cells use the toward gradient only (the away BFS is
skipped); larger flats use the full gradient. Per flat, so there are no
seams.

Implementation: fully vectorised - edge finding by neighbour shifts,
labelling by union-find with pointer jumping over equal-elevation links,
BFS by whole frontiers over a compact adjacency list - except the exit
rule for low-edge cells, which exits only to cells routed before the pass
(order-independent and cycle-free; v0.8.3's rule could form 2-cell cycles).

Iteration: a pass only resolves flats whose outlet is already routed;
passes repeat until nothing more resolves (flats drain in hierarchy order).

Reference: Barnes, R., Lehman, C., Mulla, D. (2014). Computers &
Geosciences 62, 128-135. doi:10.1016/j.cageo.2013.01.009
"""

import numpy as np

from .grid import DROW, DCOL

_IS_CARDINAL = np.array([True, False, True, False, True, False, True, False])
_TOL = 1e-9
# Neighbour scan order for assignment and exits: W, NW, N, NE, E, SE, S, SW
# (QEHT indices), the order of Barnes' reference implementation (RichDEM).
# With cardinal-over-diagonal on ties this reproduces RichDEM's Barnes flat
# resolution cell for cell (WP-G benchmark, v0.13).
_ORDER = (4, 5, 6, 7, 0, 1, 2, 3)


def _neighbours(r, c, k, rows, cols):
    nr = r + int(DROW[k]); nc = c + int(DCOL[k])
    inside = (nr >= 0) & (nr < rows) & (nc >= 0) & (nc < cols)
    return np.where(inside, nr, 0), np.where(inside, nc, 0), inside


def _components(n, a, b):
    """Connected components of n nodes with undirected edges (a, b)."""
    parent = np.arange(n, dtype=np.int64)
    if a.size == 0:
        return parent
    while True:
        pa, pb = parent[a], parent[b]
        m = np.minimum(pa, pb)
        before = parent.copy()
        np.minimum.at(parent, pa, m)
        np.minimum.at(parent, pb, m)
        while True:                               # pointer jumping
            pp = parent[parent]
            if np.array_equal(pp, parent):
                break
            parent = pp
        if np.array_equal(parent, before):
            return parent


def _csr(n, src, dst):
    order = np.lexsort((np.arange(src.size), src))
    src, dst = src[order], dst[order]
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.add.at(indptr, src + 1, 1)
    np.cumsum(indptr, out=indptr)
    return indptr, dst


def _bfs(n, indptr, indices, seeds):
    """Layered multi-source BFS: 1 at the seeds, +1 per 8-neighbour step; 0 = unreached."""
    dist = np.zeros(n, dtype=np.int64)
    front = np.unique(seeds)
    if front.size == 0:
        return dist
    dist[front] = 1
    loops = 1
    while front.size:
        cnt = indptr[front + 1] - indptr[front]
        tot = int(cnt.sum())
        if tot == 0:
            break
        offs = np.repeat(indptr[front] - (np.cumsum(cnt) - cnt), cnt) + np.arange(tot)
        nb = indices[offs]
        nb = np.unique(nb[dist[nb] == 0])
        loops += 1
        dist[nb] = loops
        front = nb
    return dist


def _one_pass(elev, valid, direction, w, small_flat_cells, drain):
    rows, cols = elev.shape
    fr, fc = np.nonzero(valid & (direction < 0) & ~drain)
    n = fr.size
    if n == 0:
        return 0, 0
    idx = np.full((rows, cols), -1, dtype=np.int64)
    idx[fr, fc] = np.arange(n)
    e = elev[fr, fc]

    # -- edges and equal-elevation links -----------------------------------
    low = np.zeros(n, dtype=bool); high = np.zeros(n, dtype=bool)
    ea, eb = [], []
    for k in range(8):
        nr, nc, inside = _neighbours(fr, fc, k, rows, cols)
        ok = inside & valid[nr, nc]
        ne = elev[nr, nc]
        routed = (direction[nr, nc] >= 0) | drain[nr, nc]
        low |= ok & ((ne < e - _TOL) | (routed & (np.abs(ne - e) <= _TOL)))
        high |= ok & (ne > e + _TOL)
        j = idx[nr, nc]
        link = ok & (j >= 0) & (np.abs(ne - e) <= _TOL)
        ea.append(np.flatnonzero(link)); eb.append(j[link])
    if not low.any():
        return 0, 0
    a = np.concatenate(ea); b = np.concatenate(eb)

    # -- labels: equal-elevation components that contain a low edge --------
    root = _components(n, a, b)
    has_low = np.zeros(n, dtype=bool)
    has_low[root[low]] = True
    labelled = has_low[root]
    uroots, lab = np.unique(root, return_inverse=True)
    lab = np.where(labelled, lab + 1, 0)                   # 0 = not labelled
    n_labels = int(np.unique(lab[labelled]).size)

    keep = labelled[a]                                       # links inside labelled flats
    indptr, indices = _csr(n, a[keep], b[keep])

    # -- gradients ------------------------------------------------------------
    toward = _bfs(n, indptr, indices, np.flatnonzero(low & labelled))
    size = np.bincount(lab, minlength=int(lab.max()) + 1)
    big = size[lab] >= int(small_flat_cells) if small_flat_cells else np.ones(n, dtype=bool)
    away = _bfs(n, indptr, indices, np.flatnonzero(high & labelled & big))
    fheight = np.zeros(int(lab.max()) + 1, dtype=np.int64)
    np.maximum.at(fheight, lab, away)
    wt = float(w) * toward
    fmask = np.where(away > 0, (fheight[lab] - away) + wt, wt).astype(np.float64)
    fmask = np.where(big, fmask, toward.astype(np.float64))
    fmask[~labelled] = 0.0

    # -- assignment: lowest flat_mask among same-flat neighbours -------------
    L = np.flatnonzero(labelled)
    lr, lc = fr[L], fc[L]
    lab_grid = np.zeros((rows, cols), dtype=np.int64)
    lab_grid[fr, fc] = lab
    fm_grid = np.zeros((rows, cols), dtype=np.float64)
    fm_grid[fr, fc] = fmask
    best_v = fmask[L].copy()
    best_k = np.full(L.size, -1, dtype=np.int64)
    for k in _ORDER:
        nr, nc, inside = _neighbours(lr, lc, k, rows, cols)
        same = inside & (lab_grid[nr, nc] == lab[L])
        v = fm_grid[nr, nc]
        better = same & ((v < best_v) | ((v == best_v) & (best_k >= 0)
                                         & ~_IS_CARDINAL[np.maximum(best_k, 0)]
                                         & _IS_CARDINAL[k]))
        best_v = np.where(better, v, best_v)
        best_k = np.where(better, k, best_k)

    # Low-edge cells (toward = 1) leave the flat directly, to a neighbour
    # that was routed before this pass or is a boundary drain. In Barnes
    # 2014 the routed low-edge cells carry the smallest flat_mask of all, so
    # every flat cell next to one drains into it (v0.13; up to v0.12 such a
    # cell could first run along the rim of the flat, which is where QEHT
    # departed from the reference implementation). The remaining local
    # minima, if any, use the same exit. Exits only ever go to cells routed
    # BEFORE the pass, so they are order-independent and cannot form cycles.
    # Cardinal preferred over diagonal, then W, NW, N, NE ... order.
    best_k[low[L]] = -1
    d0 = direction.copy()
    direction[lr[best_k >= 0], lc[best_k >= 0]] = best_k[best_k >= 0]
    resolved = int((best_k >= 0).sum())
    lm = np.flatnonzero(best_k < 0)
    if lm.size:
        mr, mc = lr[lm], lc[lm]
        ev = elev[mr, mc]
        exit_k = np.full(lm.size, -1, dtype=np.int64)
        for k in _ORDER:
            nr, nc, inside = _neighbours(mr, mc, k, rows, cols)
            ok = inside & valid[nr, nc] & (
                (elev[nr, nc] < ev - _TOL)
                | (((d0[nr, nc] >= 0) | drain[nr, nc]) & (np.abs(elev[nr, nc] - ev) <= _TOL)))
            take = ok & ((exit_k < 0) | (~_IS_CARDINAL[np.maximum(exit_k, 0)] & _IS_CARDINAL[k]))
            exit_k = np.where(take, k, exit_k)
        got = exit_k >= 0
        direction[mr[got], mc[got]] = exit_k[got]
        resolved += int(got.sum())
    return resolved, n_labels


def resolve_flats(elev, valid, direction, w=2.0, small_flat_cells=0, max_iterations=50,
                  drain=None):
    """Iterated Barnes-family flat resolution. Modifies `direction` in place.

    w : toward-gradient weight (w = 2 is Barnes 2014; must be > 1).
    small_flat_cells : flats smaller than this use the toward gradient only
        (0 = off).
    drain : optional bool grid of boundary cells without descent that drain
        off the grid or into NoData (v0.13). They are outlets: a flat that
        touches one drains to it, as in TauDEM and RichDEM. None = no drains
        (the v0.12 behaviour, where such flats stayed unrouted).
    """
    if not w > 1.0:
        raise ValueError("w must be greater than 1 (w <= 1 can create false sinks)")
    elev = np.asarray(elev, dtype=np.float64)
    drain = np.zeros(elev.shape, dtype=bool) if drain is None else np.asarray(drain, dtype=bool)
    total = passes = first = 0
    for it in range(max_iterations):
        resolved, n = _one_pass(elev, valid, direction, w, small_flat_cells, drain)
        if it == 0:
            first = n
        passes += 1
        total += resolved
        if resolved == 0:
            break
    method = "Barnes 2014 (iterated)" if w == 2.0 and not small_flat_cells else \
        f"hybrid w={w:g}" + (f", toward-only below {small_flat_cells} cells" if small_flat_cells else "")
    return direction, {
        "flat_cells_resolved": total,
        "flat_regions_first_pass": first,
        "iterations": passes,
        "outletless_flat_cells": int((valid & (direction < 0) & ~drain).sum()),
        "method": method,
        "w": float(w),
        "small_flat_cells": int(small_flat_cells),
    }


def resolve_flats_barnes(elev, valid, direction, max_iterations=50, drain=None):
    """Iterated Barnes 2014 (w = 2). Kept for API compatibility."""
    return resolve_flats(elev, valid, direction, w=2.0, max_iterations=max_iterations, drain=drain)
