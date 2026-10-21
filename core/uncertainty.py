# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Uncertainty of the flood extent and depth.

Two parts.

Evidence class (per flood cell, what the extent is based on)
    1 observed water        water seen in the image at the overpass
                            (ponding, as on the Kano plains, Nyando 2024)
    2 trace - inundation    vegetation lost between two dates but still
                            vegetated / not bare afterwards: water passed
    3 trace - erosion or    bare soil or fresh sediment after the event:
      deposition            scour, bank collapse, debris (Mai Mahiu 2024).
                            The trace edge marks the erosion width, not the
                            water width, and the pre-event DEM no longer
                            describes the ground: depths are least certain
    4 steep flow path       ground steeper than the steep-slope limit:
                            water flowed, it did not pond
    0 polygon (given)       the extent came from a user layer

Ensemble (Monte Carlo)
    Every member re-runs the surface and depth with
    * a spatially correlated DEM error: Gaussian field, standard deviation
      sigma_dem, correlation length L (Gaussian kernel, FFT);
    * an extra independent terrain error sigma_erosion in class-3 cells;
    * a water-surface (method) error: the residual SD of a spatially
      blocked cross-validation of the waterline (method_error), applied as
      a correlated field - how far the real surface departs from what the
      flood edge can tell, which the DEM term does not cover;
    * a perturbed extent: for image extents the classification threshold
      is drawn across the near-optimal Otsu band (between-class variance
      >= 95 % of its maximum) and the area-fraction cut-off used to move
      from the image grid to the analysis grid is drawn in [0.3, 0.7];
      for polygon extents and traces the edge is shifted by whole cells
      (+-1; +-2 in trace mode).
    Outputs: flood probability (share of members flooded), depth P10 / P50
    / P90 (a dry member counts as 0 m), and the spread of area, mean depth
    and volume across members.
"""

import math
import os

import numpy as np

from . import surface as SF
from . import classify as C

EVIDENCE = {0: "extent from polygons", 1: "observed water", 2: "trace - inundation",
            3: "trace - erosion / deposition", 4: "steep flow path"}
EVIDENCE_COLORS = {0: (150, 150, 150, 255), 1: (33, 113, 181, 255), 2: (116, 196, 118, 255),
                   3: (217, 95, 14, 255), 4: (203, 24, 29, 255), 255: (0, 0, 0, 0)}

DEM_PRESETS = (  # (key, label, sigma m, file-name patterns)
    ("lidar", "LiDAR / drone survey (0.15 m)", 0.15, ("lidar", "drone", "uav", "dtm_1m", "survey")),
    ("fabdem", "FABDEM bare-earth 30 m (1.5 m)", 1.5, ("fabdem",)),
    ("copernicus", "Copernicus GLO-30 DSM (2.0 m)", 2.0, ("copernicus", "glo30", "glo-30", "cop30", "cop-dem")),
    ("tandemx", "TanDEM-X 30 m (2.0 m)", 2.0, ("tandem", "tdx")),
    ("alos", "ALOS AW3D30 DSM (2.5 m)", 2.5, ("alos", "aw3d")),
    ("srtm", "SRTM 30 m (4.0 m)", 4.0, ("srtm", "nasadem")),
)


def dem_sigma_from_name(path, default=2.0):
    """(sigma, label) guessed from the DEM file name; `default` if unknown."""
    n = os.path.basename(str(path)).lower()
    for key, label, sig, pats in DEM_PRESETS:
        if any(p in n for p in pats):
            return sig, label
    return default, f"unknown DEM - {default:g} m assumed"


def correlated_noise(shape, sigma, corr_cells, rng):
    """Zero-mean Gaussian random field with standard deviation `sigma` and a
    Gaussian correlation of `corr_cells` (std of the kernel, in cells)."""
    rows, cols = shape
    white = rng.standard_normal(shape)
    if corr_cells <= 0.5:
        return (sigma * white).astype("float32")
    ky = np.fft.fftfreq(rows)[:, None]
    kx = np.fft.rfftfreq(cols)[None, :]
    h = np.exp(-2.0 * (math.pi ** 2) * (corr_cells ** 2) * (kx ** 2 + ky ** 2))
    f = np.fft.irfft2(np.fft.rfft2(white) * h, s=shape)
    sd = f.std()
    if sd <= 0:
        return np.zeros(shape, "float32")
    return (sigma * f / sd).astype("float32")


def evidence_classes(flood, source, bare=None, slope=None, steep_pct=10.0):
    """uint8 evidence raster (see module notes); 255 outside the flood."""
    out = np.full(flood.shape, 255, dtype="uint8")
    base = {"image": 1, "trace": 2, "polygons": 0}[source]
    out[flood] = base
    if source == "trace" and bare is not None:
        out[flood & bare] = 3
    if slope is not None:
        out[flood & (slope > steep_pct)] = 4
    return out


def block_mean(a, f):
    """Mean over f x f blocks (NaN-aware), edges padded."""
    a = np.asarray(a, dtype="float32")
    r, c = a.shape
    pr, pc = (-r) % f, (-c) % f
    p = np.pad(a, ((0, pr), (0, pc)), constant_values=np.nan)
    b = p.reshape((r + pr) // f, f, (c + pc) // f, f)
    with np.errstate(invalid="ignore"):
        n = np.sum(np.isfinite(b), axis=(1, 3))
        s = np.nansum(b, axis=(1, 3))
        out = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    return out.astype("float32")


def expand(a, f, shape):
    return np.repeat(np.repeat(a, f, axis=0), f, axis=1)[:shape[0], :shape[1]]


def method_error(elev, flood, valid, grid, other_water=None, interpolation="idw", folds=3, rng=None):
    """Spatially blocked cross-validation of the water surface along the edge.

    Waterline cells are held out in square blocks whose side is twice the
    median distance from a flood cell to the waterline (so a held-out cell
    is predicted from about as far away as an interior cell), predicted
    from the remaining waterline, and compared with their elevation.
    Returns (sd of the residuals in m, block size in cells, n held out).
    The residual SD measures how far the real water surface departs from
    what the edge can tell - the method error the DEM term does not cover.
    """
    rng = rng or np.random.default_rng(0)
    line = SF.waterline(flood, valid, other_water)
    if line.sum() < 30:
        return 0.0, 0, 0
    d = SF.proximity_cells(line)
    b = int(max(4, min(200, 2 * np.median(d[flood]))))
    rows, cols = flood.shape
    br, bc = np.indices((rows, cols))
    tile = (br // b) * ((cols + b - 1) // b) + (bc // b)
    ids = np.unique(tile[line])
    res = []
    for _ in range(folds):
        held_ids = rng.choice(ids, max(1, int(0.2 * ids.size)), replace=False)
        held = line & np.isin(tile, held_ids)
        train = line & ~held
        if train.sum() < 10 or held.sum() < 5:
            continue
        if interpolation == "idw":
            pred = SF.fill_from_points(elev, train, grid, max_search=SF.search_distance(train, held))
        else:
            pred, _ = SF.fill_nearest(elev, train, grid.cell_x, grid.cell_y, targets=held)
        r = (elev - pred)[held]
        res.append(r[np.isfinite(r)])
    if not res:
        return 0.0, b, 0
    r = np.concatenate(res)
    med = np.median(r)
    sd = float(1.4826 * np.median(np.abs(r - med)))
    return sd, b, int(r.size)


def shift_edge(mask, valid, k):
    """Grow (k > 0) or shrink (k < 0) a mask by |k| cells (4-neighbour)."""
    if k == 0:
        return mask.copy()
    if k > 0:
        return C.dilate(mask, k) & valid
    return mask & ~C.dilate(~mask, -k)


def run_ensemble(dem, valid, grid, base_flood, members=20, sigma_dem=2.0, corr_m=150.0,
                 sigma_erosion=1.0, erosion=None, frac_options=None, edge_shift=1,
                 other_water=None, seed=42, surface_kwargs=None, depth_base=None,
                 seed_mask=None, max_cells=2_000_000, method_cv=True, progress=None):
    """Monte Carlo depth ensemble on the flood's bounding box.

    frac_options: list of area-fraction rasters (image extents at the low,
        optimal and high ends of the threshold band) on the analysis grid;
        a member draws one and a cut-off in [0.3, 0.7]. None -> edge shifts.
    depth_base: 'dem' (depth = surface - DEM) or a HAND array of the same
        shape (depth = surface - HAND; the DEM error is applied to HAND).
    Returns dict of full-grid arrays (prob, p10, p50, p90) and summary.
    """
    rng = np.random.default_rng(seed)
    surface_kwargs = dict(surface_kwargs or {})
    surface_kwargs.pop("progress", None)
    flood_any = base_flood.copy()
    if frac_options:
        for f in frac_options:
            flood_any |= np.nan_to_num(f) >= 0.3
    m = 3 + 2 * int(edge_shift)
    rr, cc = np.nonzero(flood_any)
    r0, r1 = max(0, rr.min() - m), min(dem.shape[0], rr.max() + m + 1)
    c0, c1 = max(0, cc.min() - m), min(dem.shape[1], cc.max() + m + 1)
    win = (slice(r0, r1), slice(c0, c1))
    elev = (depth_base if depth_base is not None else dem)[win].astype("float32")
    v = valid[win] & np.isfinite(elev)
    ow = other_water[win] if other_water is not None else None
    ero = erosion[win] if erosion is not None else None
    bf = base_flood[win]
    fo = [f[win] for f in frac_options] if frac_options else None
    sm = seed_mask[win] if seed_mask is not None else None
    # coarsen large windows: the bands do not need full resolution
    full_shape = elev.shape
    f = int(math.ceil(math.sqrt(full_shape[0] * full_shape[1] / float(max_cells)))) if max_cells else 1
    if f > 1:
        elev = block_mean(elev, f)
        v = block_mean(v.astype("float32"), f) >= 0.5
        v &= np.isfinite(elev)
        bf = block_mean(bf.astype("float32"), f) >= 0.5
        ow = (block_mean(ow.astype("float32"), f) >= 0.5) if ow is not None else None
        ero = (block_mean(ero.astype("float32"), f) >= 0.5) if ero is not None else None
        sm = (block_mean(sm.astype("float32"), f) > 0) if sm is not None else None
        fo = [block_mean(np.nan_to_num(x), f) for x in fo] if fo else None
    cx, cy = grid.cell_x * f, grid.cell_y * f
    sub_gt = (grid.gt[0] + c0 * grid.gt[1], cx, 0.0, grid.gt[3] + r0 * grid.gt[5], 0.0, -cy)
    from .raster import Grid
    sub = Grid(elev.shape[1], elev.shape[0], sub_gt, grid.wkt)
    corr_cells = corr_m / max(cx, 1e-9)
    shape = elev.shape
    m_sd, m_block, m_n = (0.0, 0, 0)
    if method_cv:
        m_sd, m_block, m_n = method_error(elev, bf, v, sub, ow, surface_kwargs.get("interpolation", "idw"),
                                          rng=np.random.default_rng(seed + 1))
    store = []
    wet = np.zeros(shape, dtype="uint16")
    areas, means, vols, failed = [], [], [], 0
    ca = cx * cy
    for i in range(members):
        if progress:
            progress(i / float(members), f"Uncertainty member {i + 1} of {members}")
        e = elev + correlated_noise(shape, sigma_dem, corr_cells, rng)
        if ero is not None and sigma_erosion > 0:
            e = np.where(ero, e + rng.normal(0.0, sigma_erosion, shape).astype("float32"), e)
        if fo:
            fr = fo[int(rng.integers(0, len(fo)))]
            fl = v & (np.nan_to_num(fr) >= rng.uniform(0.3, 0.7))
        else:
            k = int(rng.integers(-edge_shift, edge_shift + 1))
            fl = shift_edge(bf, v, k)
        if ow is not None:
            fl &= ~ow
        if sm is not None:
            fl, _, _ = C.keep_connected(fl, sm)
        if fl.sum() < 3:
            failed += 1
            continue
        try:
            s, _, _, _ = SF.estimate_surface(e, fl, v, sub, other_water=ow, **surface_kwargs)
        except ValueError:
            failed += 1
            continue
        if m_sd > 0:
            s = s + correlated_noise(shape, m_sd, max(1.0, m_block / 2.0), rng)
        d, _ = SF.depth_from_surface(s, e, fl)
        dd = np.nan_to_num(d, nan=0.0)
        store.append(dd)
        wet += fl.astype("uint16")
        n = int(fl.sum())
        areas.append(n * ca / 1e6)
        means.append(float(dd[fl].mean()))
        vols.append(float(dd[fl].sum() * ca))
    ok = members - failed
    if ok < 2:
        raise ValueError("Fewer than two uncertainty members succeeded.")
    depths = np.stack(store)
    prob = wet.astype("float32") / float(ok)
    q = np.percentile(depths, [10, 50, 90], axis=0).astype("float32")
    if f > 1:
        prob = expand(prob, f, full_shape)
        q = np.stack([expand(x, f, full_shape) for x in q])
        v = valid[win] & np.isfinite((depth_base if depth_base is not None else dem)[win])
    full = {}
    for name, arr in (("prob", prob), ("p10", q[0]), ("p50", q[1]), ("p90", q[2])):
        a = np.full(dem.shape, np.nan, dtype="float32")
        a[win] = arr
        if name != "prob":
            a[win] = np.where(prob > 0, arr, np.nan)
        else:
            a[win] = np.where(v, arr, np.nan)
        full[name] = a

    def pct(x):
        x = np.asarray(x, dtype="float64")
        return [float(np.percentile(x, 10)), float(np.percentile(x, 50)), float(np.percentile(x, 90))]
    summary = {"members": members, "succeeded": ok, "sigma_dem_m": sigma_dem, "corr_length_m": corr_m,
               "sigma_erosion_m": sigma_erosion if ero is not None else 0.0,
               "extent_perturbation": "threshold band + area fraction" if fo else f"edge shift +-{edge_shift} cell(s)",
               "area_km2_p10_p50_p90": pct(areas), "mean_depth_m_p10_p50_p90": pct(means),
               "volume_m3_p10_p50_p90": pct(vols), "seed": seed, "coarsening": f,
               "method_sd_m": m_sd, "method_block_m": m_block * cx, "method_cv_cells": m_n}
    return full, summary
