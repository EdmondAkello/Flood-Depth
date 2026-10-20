# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Water classification: indices, automatic thresholds, cleanup.

Mask convention everywhere: a pair of boolean arrays (water, valid).
Written rasters use 1 = water, 0 = dry, 255 = NoData.

Automatic thresholds
--------------------
Default: edge-based Otsu (Donchyts et al. 2016) seeded at the physical
threshold (0 for NDWI / MNDWI, -16 dB VV, -22 dB VH) - see edge_threshold.

Otsu (1979) picks the histogram split that maximises between-class
variance. Its separability eta = sigma_between^2 / sigma_total^2 (0-1)
says how bimodal the histogram really is. A single Gaussian gives
eta = 2/pi ~ 0.64, so FloodDepth calls eta >= 0.75 good, 0.65-0.75
moderate and warns below 0.65. Flood scenes are usually dominated by dry land, which drags a
global Otsu threshold into the land mode. The default "tile" method
therefore follows the split-based approach of Martinis et al. (2009) and
Chini et al. (2017): the image is cut into tiles, the tiles that contain
both classes (by the global threshold) and have the largest spread are
pooled, and Otsu is applied to that pooled, balanced sample. When fewer
than five such tiles exist the global threshold is kept and the report
says so.
"""

import math

import numpy as np

from . import raster as R


# -- indices ----------------------------------------------------------

def normalized_difference(a, b):
    """(a - b) / (a + b); NaN where either is NaN or a + b is ~0."""
    a = np.asarray(a, dtype="float32")
    b = np.asarray(b, dtype="float32")
    den = a + b
    out = np.full(a.shape, np.nan, dtype="float32")
    ok = np.isfinite(a) & np.isfinite(b) & (np.abs(den) > 1e-6)
    out[ok] = (a[ok] - b[ok]) / den[ok]
    return np.clip(out, -1.0, 1.0)


def to_db(x):
    x = np.asarray(x, dtype="float32")
    out = np.full(x.shape, np.nan, dtype="float32")
    ok = np.isfinite(x) & (x > 0)
    out[ok] = 10.0 * np.log10(x[ok])
    return out


# -- thresholds --------------------------------------------------------

class ThresholdError(ValueError):
    pass


def otsu(values, bins=512, value_range=None):
    """Otsu threshold of a 1-D sample. Returns (threshold, eta)."""
    v = np.asarray(values, dtype="float64")
    v = v[np.isfinite(v)]
    if v.size < 50:
        raise ThresholdError("fewer than 50 valid values - cannot choose a threshold automatically")
    lo, hi = value_range if value_range else (np.percentile(v, 0.1), np.percentile(v, 99.9))
    if not (hi - lo) > 1e-9:
        raise ThresholdError("the values are (almost) constant - there is nothing to separate")
    hist, edges = np.histogram(np.clip(v, lo, hi), bins=bins, range=(lo, hi))
    return otsu_from_histogram(hist, edges)


def otsu_from_histogram(hist, edges):
    counts = np.asarray(hist, dtype="float64")
    edges = np.asarray(edges, dtype="float64")
    if counts.ndim != 1 or edges.size != counts.size + 1:
        raise ThresholdError("histogram counts and edges do not match")
    total = counts.sum()
    if total <= 0:
        raise ThresholdError("empty histogram")
    centers = 0.5 * (edges[:-1] + edges[1:])
    w1 = np.cumsum(counts)
    w2 = total - w1
    s1 = np.cumsum(counts * centers)
    mu_t = s1[-1] / total
    ok = (w1 > 0) & (w2 > 0)
    if not ok.any():
        raise ThresholdError("the histogram has a single occupied bin")
    mu1 = np.where(w1 > 0, s1 / np.where(w1 > 0, w1, 1), 0.0)
    mu2 = np.where(w2 > 0, (s1[-1] - s1) / np.where(w2 > 0, w2, 1), 0.0)
    between = np.where(ok, w1 * w2 * (mu1 - mu2) ** 2 / total ** 2, -1.0)
    k = int(np.argmax(between))
    # A flat maximum (empty bins between the classes - common with 8-bit
    # data or well separated classes) makes every threshold in the gap
    # equally optimal; take the middle of the gap, not its first edge.
    top = between[k] * (1.0 - 1e-9)
    k1 = k
    while k1 + 1 < between.size and between[k1 + 1] >= top:
        k1 += 1
    mid = 0.5 * (edges[k + 1] + edges[k1 + 1])   # = edges[k + 1] when there is no plateau
    var_t = float(np.sum(counts * (centers - mu_t) ** 2) / total)
    eta = float(between[k] / var_t) if var_t > 0 else 0.0
    return float(mid), max(0.0, min(1.0, eta))


def otsu_band(values, frac=0.95, bins=512):
    """Thresholds whose between-class variance is >= `frac` of the Otsu
    maximum: the range over which the split is almost equally good. Used as
    the threshold uncertainty. Returns (lo, hi)."""
    v = np.asarray(values, dtype="float64")
    v = v[np.isfinite(v)]
    lo, hi = np.percentile(v, 0.1), np.percentile(v, 99.9)
    hist, edges = np.histogram(np.clip(v, lo, hi), bins=bins, range=(lo, hi))
    c = hist.astype("float64")
    tot = c.sum()
    ctr = 0.5 * (edges[:-1] + edges[1:])
    w1 = np.cumsum(c)
    w2 = tot - w1
    s1 = np.cumsum(c * ctr)
    ok = (w1 > 0) & (w2 > 0)
    mu1 = np.where(w1 > 0, s1 / np.where(w1 > 0, w1, 1), 0)
    mu2 = np.where(w2 > 0, (s1[-1] - s1) / np.where(w2 > 0, w2, 1), 0)
    b = np.where(ok, w1 * w2 * (mu1 - mu2) ** 2, -1)
    k = int(np.argmax(b))
    good = np.nonzero(b >= frac * b[k])[0]
    # contiguous run around the maximum
    i0 = k
    while i0 - 1 >= 0 and (i0 - 1) in set(good.tolist()):
        i0 -= 1
    i1 = k
    gs = set(good.tolist())
    while i1 + 1 < len(b) and (i1 + 1) in gs:
        i1 += 1
    return float(edges[i0 + 1]), float(edges[i1 + 1])


def tile_threshold(index, water_above=True, tile=None, min_tiles=5, max_tiles=64):
    """Split-based Otsu (see module notes). Returns a dict with threshold,
    eta, method, tiles_used and the global values for the report."""
    idx = np.asarray(index, dtype="float32")
    finite = idx[np.isfinite(idx)]
    t0, eta0 = otsu(finite)
    out = {"threshold": t0, "eta": eta0, "method": "global Otsu", "tiles_used": 0,
           "global_threshold": t0, "global_eta": eta0, "band": otsu_band(finite)}
    rows, cols = idx.shape
    n = rows * cols
    if tile is None:
        tile = int(max(16, min(256, math.sqrt(n / 400.0))))
    pr, pc = (-rows) % tile, (-cols) % tile
    padded = np.pad(idx, ((0, pr), (0, pc)), constant_values=np.nan)
    tiles = padded.reshape((rows + pr) // tile, tile, (cols + pc) // tile, tile).swapaxes(1, 2)
    tiles = tiles.reshape(-1, tile * tile)
    nvalid = np.sum(np.isfinite(tiles), axis=1)
    keep = nvalid >= 0.5 * tile * tile
    if not keep.any():
        return out
    t = tiles[keep]
    with np.errstate(invalid="ignore"):
        wet = (t > t0) if water_above else (t < t0)
        frac = np.sum(wet & np.isfinite(t), axis=1) / np.sum(np.isfinite(t), axis=1)
        sd = np.nanstd(t, axis=1)
    mixed = (frac >= 0.1) & (frac <= 0.9)
    if mixed.sum() < min_tiles:
        out["note"] = (f"only {int(mixed.sum())} tile(s) contain both water and land; "
                       "the global Otsu threshold was kept")
        return out
    order = np.argsort(-sd[mixed])[:max_tiles]
    pool = t[mixed][order].ravel()
    try:
        t1, eta1 = otsu(pool)
    except ThresholdError:
        return out
    out.update({"threshold": t1, "eta": eta1, "method": "tile-based Otsu", "tiles_used": int(order.size),
                "tile_size": tile})
    try:
        out["band"] = otsu_band(pool)
    except (ValueError, IndexError) as exc:
        out["band_note"] = f"threshold band not computed: {exc}"
    return out


UNIMODAL_ETA = 2.0 / 3.141592653589793   # eta of a single Gaussian = 2/pi ~ 0.637


PHYSICAL_THRESHOLD = {"ndwi": 0.0, "mndwi": 0.0, "sar_vv": -16.0, "sar_vh": -22.0}


def edge_threshold(index, water_above=True, t0=0.0, buf=2, min_pixels=200):
    """Edge-based Otsu (Donchyts et al. 2016).

    Classify with the physical threshold t0, take the pixels within `buf`
    cells of the resulting water / land edges - a sample in which both
    classes are present in similar numbers - and apply Otsu to it. On the
    446 Sen1Floods11 chips this scores CSI 0.61 for MNDWI against 0.39 for
    the tile-based method (benchmarks/sen1floods11.py). Falls back to t0
    (and says so) when one class is absent or the edge is too short.
    """
    x = np.asarray(index, dtype="float32")
    fin = np.isfinite(x)
    with np.errstate(invalid="ignore"):
        w = ((x > t0) if water_above else (x < t0)) & fin
    out = {"threshold": float(t0), "eta": None, "method": "physical threshold (fallback)", "tiles_used": 0,
           "prior": float(t0)}
    if w.sum() < 20 or (fin & ~w).sum() < 20:
        out["note"] = "one class (almost) absent at the physical threshold; it was kept"
        return out
    edge = dilate(w, buf) & dilate(fin & ~w, buf) & fin
    if edge.sum() < min_pixels:
        out["note"] = f"water / land edge too short ({int(edge.sum())} pixels); the physical threshold was kept"
        return out
    t, eta = otsu(x[edge])
    out.update({"threshold": t, "eta": eta, "method": "edge-based Otsu", "edge_pixels": int(edge.sum())})
    try:
        out["band"] = otsu_band(x[edge])
    except (ValueError, IndexError) as exc:
        out["band_note"] = f"threshold band not computed: {exc}"
    return out


def confidence_label(eta):
    """good >= 0.75; moderate 0.65-0.75; low below 0.65, i.e. no better
    than splitting one Gaussian population (eta = 2/pi)."""
    if eta >= 0.75:
        return "good"
    if eta >= 0.65:
        return "moderate"
    return "low"


def classify(index, threshold, water_above=True):
    """(water, valid) from an index array and a threshold."""
    valid = np.isfinite(index)
    with np.errstate(invalid="ignore"):
        water = (index > threshold) if water_above else (index < threshold)
    return water & valid, valid


# -- connected components (in C, via GDAL) ------------------------------

def label_components(mask, connectivity=8):
    """Label 4- or 8-connected True regions. Returns (labels int32, count).

    Uses gdal.Polygonize then rasterizes each polygon with its id: both
    are C++ calls, so this is fast on large grids with no SciPy.
    """
    gdal = R._gdal()
    from osgeo import ogr
    mask = np.asarray(mask, dtype=bool)
    rows, cols = mask.shape
    if not mask.any():
        return np.zeros(mask.shape, "int32"), 0
    mem = gdal.GetDriverByName("MEM").Create("", cols, rows, 1, gdal.GDT_Byte)
    mem.SetGeoTransform((0.0, 1.0, 0.0, float(rows), 0.0, -1.0))
    b = mem.GetRasterBand(1)
    b.WriteArray(mask.astype("uint8"))
    vds = ogr.GetDriverByName("Memory").CreateDataSource("cc")
    lyr = vds.CreateLayer("cc", geom_type=ogr.wkbPolygon)
    lyr.CreateField(ogr.FieldDefn("v", ogr.OFTInteger))
    opts = ["8CONNECTED=8"] if connectivity == 8 else []
    with R.quiet():
        gdal.Polygonize(b, b, lyr, 0, opts, callback=None)
    lab = gdal.GetDriverByName("MEM").Create("", cols, rows, 1, gdal.GDT_Int32)
    lab.SetGeoTransform((0.0, 1.0, 0.0, float(rows), 0.0, -1.0))
    lyr2 = vds.CreateLayer("ids", geom_type=ogr.wkbPolygon)
    lyr2.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    n = 0
    lyr.ResetReading()
    for f in lyr:
        if f.GetField("v") != 1:
            continue
        n += 1
        g = ogr.Feature(lyr2.GetLayerDefn())
        g.SetGeometry(f.GetGeometryRef().Clone())
        g.SetField("id", n)
        lyr2.CreateFeature(g)
    with R.quiet():
        gdal.RasterizeLayer(lab, [1], lyr2, options=["ATTRIBUTE=id"])
    labels = lab.GetRasterBand(1).ReadAsArray().astype("int32")
    labels[~mask] = 0
    return labels, n


def clean_mask(water, valid, min_water_cells=0, max_hole_cells=0):
    """Remove water patches smaller than `min_water_cells` and fill dry
    holes smaller than `max_hole_cells` (holes touching NoData or the grid
    edge are never filled). Each operation is separate, so islands, roads
    and buildings are not filled unless asked. Returns (mask, stats)."""
    out = water & valid
    stats = {"water_patches_removed": 0, "water_cells_removed": 0,
             "holes_filled": 0, "hole_cells_filled": 0}
    if min_water_cells and min_water_cells > 1:
        lab, n = label_components(out)
        if n:
            size = np.bincount(lab.ravel(), minlength=n + 1)
            small = np.zeros(n + 1, bool)
            small[1:] = size[1:] < min_water_cells
            rm = small[lab]
            stats["water_patches_removed"] = int(small.sum())
            stats["water_cells_removed"] = int(rm.sum())
            out = out & ~rm
    if max_hole_cells and max_hole_cells > 0:
        dry = valid & ~out
        lab, n = label_components(dry, connectivity=4)
        if n:
            size = np.bincount(lab.ravel(), minlength=n + 1)
            # holes touching NoData or the edge are open land, not holes
            edge = np.zeros_like(valid)
            edge[0, :] = edge[-1, :] = edge[:, 0] = edge[:, -1] = True
            near_invalid = _dilate(~valid, 1) | edge
            open_ids = np.unique(lab[near_invalid & dry])
            fill = np.zeros(n + 1, bool)
            fill[1:] = size[1:] <= max_hole_cells
            fill[open_ids] = False
            fill[0] = False
            add = fill[lab]
            stats["holes_filled"] = int(fill.sum())
            stats["hole_cells_filled"] = int(add.sum())
            out = out | add
    return out, stats


def keep_connected(water, seed):
    """Keep only the water regions that contain at least one seed cell."""
    lab, n = label_components(water)
    if n == 0:
        return water.copy(), 0, 0
    ids = np.unique(lab[seed & water])
    ids = ids[ids > 0]
    keep = np.isin(lab, ids)
    return keep, int(ids.size), int(n - ids.size)


def _dilate(mask, r):
    out = mask.copy()
    for _ in range(int(r)):
        m = out.copy()
        m[1:, :] |= out[:-1, :]
        m[:-1, :] |= out[1:, :]
        m[:, 1:] |= out[:, :-1]
        m[:, :-1] |= out[:, 1:]
        out = m
    return out


def dilate(mask, cells):
    return _dilate(np.asarray(mask, bool), cells)


# -- supervised (advanced) -----------------------------------------------

def fit_centroids(water_samples, land_samples, min_samples=10):
    """Standardised two-class nearest-centroid model (n_samples x n_bands)."""
    w = np.asarray(water_samples, dtype="float64")
    l = np.asarray(land_samples, dtype="float64")
    if w.ndim != 2 or l.ndim != 2 or w.shape[1] != l.shape[1]:
        raise ValueError("Training samples must be 2-D arrays with matching band counts.")
    w = w[np.all(np.isfinite(w), axis=1)]
    l = l[np.all(np.isfinite(l), axis=1)]
    if len(w) < min_samples or len(l) < min_samples:
        raise ValueError(f"Not enough valid training pixels (water {len(w)}, land {len(l)}; "
                         f"at least {min_samples} each).")
    allv = np.vstack((w, l))
    mean = allv.mean(axis=0)
    scale = allv.std(axis=0)
    scale[~np.isfinite(scale) | (scale < 1e-8)] = 1.0
    return mean, scale, ((w - mean) / scale).mean(axis=0), ((l - mean) / scale).mean(axis=0)


def classify_centroids(stack, model):
    """stack (bands, rows, cols) -> (water, valid)."""
    mean, scale, cw, cl = model
    a = np.asarray(stack, dtype="float32")
    valid = np.all(np.isfinite(a), axis=0)
    z = (a - mean[:, None, None]) / scale[:, None, None]
    dw = np.sum((z - cw[:, None, None]) ** 2, axis=0)
    dl = np.sum((z - cl[:, None, None]) ** 2, axis=0)
    return (dw < dl) & valid, valid
