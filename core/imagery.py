# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""From satellite files to a water mask, automatically.

classify_image() is the whole classification step: identify the bands,
read them onto one grid, build the water index, mask bright cloud, pick a
threshold, clean the mask, and optionally remove permanent water seen in a
pre-event image. Every decision is recorded in `info` for the report.

Modes
-----
water  : water at the time of the image (NDWI, MNDWI, or dark SAR).
trace  : flood TRACE after the water has gone - vegetation lost between a
         pre-event and a post-event image (NDVI drop). For flash floods
         and debris flows (e.g. Mai Mahiu, 29 April 2024) the water drains
         within hours and no optical image shows it; the stripped
         vegetation and fresh sediment do. A trace is not a water extent:
         depths derived from it describe the flood that left the trace and
         carry more uncertainty.
"""

import os

import numpy as np

from . import raster as R
from . import classify as C
from . import sensors as S

INDEX_NAMES = {"ndwi": "NDWI (green, NIR) - McFeeters 1996",
               "mndwi": "MNDWI (green, SWIR1) - Xu 2006",
               "sar": "SAR backscatter (dark = water)",
               "dndvi": "NDVI change post - pre (vegetation loss)"}


def scene_from_paths(paths, overrides=None):
    paths = [str(p) for p in paths if p]
    if not paths:
        raise ValueError("No image file was given.")
    for p in paths:
        if not os.path.exists(p):
            raise ValueError(f"{p} does not exist.")
    return S.identify(paths, R.band_info, overrides)


def _ref_grid(scene, roles):
    for r in roles:
        if r in scene.sources:
            p, b = scene.sources[r]
            return R.grid_of(p)
    p, b = next(iter(scene.sources.values()))
    return R.grid_of(p)


def read_role(scene, role, grid, notes):
    p, b = scene.sources[role]
    a, g, nt = R.read_band(p, b)
    for t in nt:
        if t not in notes:
            notes.append(f"{os.path.basename(p)}: {t}")
    off = scene.offset.get(role, 0.0)
    if off:
        a = a + np.float32(off)
    return R.warp_array(a, g, grid)


def decide_offset(scene, notes):
    """Apply the Sentinel-2 -1000 BOA offset only when the data show it.

    Since processing baseline 04.00 (25 Jan 2022) L2A digital numbers are
    reflectance x 10000 + 1000. Exports that harmonise the values (EO /
    Copernicus Browser, most cloud archives) or rescale them (Browser uint16
    'Raw' = reflectance x 65535) must NOT be shifted. The offset is applied
    when every water-index band has its 1st percentile >= 1000 and its
    99th percentile <= 12000 (the offset DN window) - otherwise not.
    """
    if not scene.offset_candidate:
        return
    stats = []
    for role in ("green", "nir", "swir1"):
        if role not in scene.sources:
            continue
        p, b = scene.sources[role]
        a, _, _ = R.read_band(p, b)
        v = a[np.isfinite(a) & (a > 0)]
        if v.size == 0:
            return
        stats.append((float(np.percentile(v, 1)), float(np.percentile(v, 99))))
    if stats and all(lo >= 1000.0 and hi <= 12000.0 for lo, hi in stats):
        for role in scene.sources:
            scene.offset[role] = -1000.0
        notes.append("Sentinel-2 L2A digital numbers carry the +1000 BOA offset (baseline 04.00+; every band's "
                     "values lie in the offset DN range); it is removed before computing indices.")
    else:
        notes.append("Sentinel-2 L2A after Jan 2022 in integer form, but the values are not in the offset DN "
                     "range (harmonised or rescaled export): no BOA offset applied.")


def choose_index(scene, index="auto"):
    index = (index or "auto").lower()
    if index == "auto":
        if scene.kind == "sar":
            return "sar"
        if scene.rendered_sar:
            raise ValueError(_cannot_classify(scene))
        if scene.has("green", "swir1"):
            return "mndwi"
        if scene.has("green", "nir"):
            return "ndwi"
        raise ValueError(_cannot_classify(scene))
    need = {"ndwi": ("green", "nir"), "mndwi": ("green", "swir1"), "sar": ()}[index]
    if index == "sar" and scene.kind != "sar":
        raise ValueError("SAR thresholding needs a Sentinel-1 product (VV band or a Sentinel-1 rendered image).")
    if need and not scene.has(*need):
        raise ValueError(f"{index.upper()} needs the {' and '.join(need)} bands; found: {scene.describe()}. "
                         + _cannot_classify(scene))
    return index


def _cannot_classify(scene):
    msg = ("Could not find the bands for an automatic water index. Found: " + scene.describe() + ". ")
    if scene.rendered_sar:
        msg += (f"{scene.rendered_sar} is a rendered Sentinel-1 visualisation: its colours are a display "
                "choice, not backscatter, so water cannot be thresholded reliably. Download the VV band as "
                "linear gamma0/sigma0 or dB (e.g. Copernicus Browser 'VV - decibel gamma0', or Sentinel-1 RTC). ")
    if scene.kind == "rgb":
        msg += ("This is a true-colour (RGB) image: water indices need NIR or SWIR. Use the Raw bands "
                "(B03 + B08 or B11), the 'False color' product, or supervised mode with training polygons. ")
    msg += ("If the bands are in one multiband file without names, set the band numbers "
            "(green, NIR, SWIR1) in the advanced parameters.")
    for n in scene.notes:
        msg += " " + n
    return msg


def compute_index(scene, index, grid=None, notes=None):
    """(index array, grid, water_above, notes)."""
    notes = [] if notes is None else notes
    if index == "sar":
        role = "vv" if "vv" in scene.sources else "vh"
        grid = grid or _ref_grid(scene, [role])
        a = read_role(scene, role, grid, notes)
        finite = a[np.isfinite(a)]
        if finite.size and np.nanpercentile(finite, 99) < 5.0 and np.nanmin(finite) >= 0:
            a = C.to_db(a)
            notes.append("SAR values look linear; converted to dB.")
        return a, grid, False, notes
    pair = {"ndwi": ("green", "nir"), "mndwi": ("green", "swir1")}[index]
    grid = grid or _ref_grid(scene, pair)
    a = read_role(scene, pair[0], grid, notes)
    b = read_role(scene, pair[1], grid, notes)
    return C.normalized_difference(a, b), grid, True, notes


def plausible_range(iname, scene, sar_band=None):
    """Range in which an automatic threshold is accepted (None = no check).

    NDWI / MNDWI from reflectance: open water separates from land between
    about -0.3 and 0.6 (McFeeters 1996; Xu 2006; Ji et al. 2009). An Otsu
    threshold outside that range means the histogram was split inside the
    land (or water) population. SAR (dB): water VV about -26 to -10 dB,
    VH -32 to -16 dB. Rendered composites are not reflectance: no check.
    """
    if iname in ("ndwi", "mndwi") and scene.kind == "bands":
        return (-0.3, 0.6)
    if iname == "sar":
        return (-26.0, -10.0) if sar_band == "vv" else (-32.0, -16.0)
    return None


def bright_cloud(scene, grid, threshold, notes):
    """Cells whose blue reflectance exceeds `threshold` (approximate)."""
    if not threshold or threshold <= 0 or "blue" not in scene.sources:
        return None
    blue = read_role(scene, "blue", grid, notes)
    fin = blue[np.isfinite(blue)]
    if fin.size == 0:
        return None
    scale = S.reflectance_scale(scene.dtype.get("blue"), float(np.percentile(fin, 99)))
    with np.errstate(invalid="ignore"):
        return np.isfinite(blue) & (blue * scale > threshold)


def ndvi(scene, grid, notes):
    if not scene.has("red", "nir"):
        raise ValueError("Flood-trace mode needs red and NIR bands (e.g. B04 and B08) in both images; "
                         f"found: {scene.describe()}.")
    return C.normalized_difference(read_role(scene, "nir", grid, notes), read_role(scene, "red", grid, notes))


def classify_image(paths, overrides=None, index="auto", threshold=None, method="tile",
                   min_water_cells=0, max_hole_cells=0, cloud_threshold=0.25,
                   pre_paths=None, pre_overrides=None, mode="water", cloud_buffer=2, progress=None):
    """Classify water (or flood trace) from image files.

    Returns dict: water, valid, grid, permanent (or None), cloud (or None),
    index (array), info (dict with every decision).
    """
    def step(f, m):
        if progress:
            progress(f, m)

    notes = []
    info = {"mode": mode}
    step(0.02, "Identifying bands")
    scene = scene_from_paths(paths, overrides)
    info.update({"scene_kind": scene.kind, "sensor": scene.family or "unknown", "date": scene.date,
                 "bands": scene.describe(), "files": [os.path.basename(p) for p in paths]})
    notes.extend(scene.notes)
    decide_offset(scene, notes)

    pre = None
    if pre_paths:
        pre = scene_from_paths(pre_paths, pre_overrides)
        info["pre_bands"] = pre.describe()
        info["pre_date"] = pre.date
        notes.extend("pre-event: " + n for n in pre.notes)
        decide_offset(pre, notes)

    if mode == "trace":
        if pre is None:
            raise ValueError("Flood-trace mode compares two dates: give the pre-event image files as well.")
        step(0.1, "NDVI before and after")
        grid = _ref_grid(scene, ["red", "nir"])
        idx = ndvi(scene, grid, notes) - ndvi(pre, grid, notes)
        water_above = False
        iname = "dndvi"
    else:
        iname = choose_index(scene, index)
        step(0.1, f"Computing {iname.upper()}")
        idx, grid, water_above, notes = compute_index(scene, iname, None, notes)
    info["index"] = INDEX_NAMES[iname]
    info["index_key"] = iname
    info["grid"] = f"{grid.cols} x {grid.rows} cells of {grid.cell_x:.6g} x {grid.cell_y:.6g} map units"

    cloud = None
    if scene.kind != "sar":
        cloud = bright_cloud(scene, grid, cloud_threshold, notes)
        if pre is not None and mode == "trace":
            pc = bright_cloud(pre, grid, cloud_threshold, notes)
            if pc is not None:
                cloud = pc if cloud is None else (cloud | pc)
        if cloud is not None:
            if cloud_buffer:
                cloud = C.dilate(cloud, int(cloud_buffer))
            info["cloud_share"] = float(cloud.mean())
            idx = np.where(cloud, np.nan, idx)
        elif scene.kind != "sar":
            notes.append("No blue band: bright cloud was not masked. Check the image for cloud.")

    valid = np.isfinite(idx)
    info["valid_share"] = float(valid.mean())
    if valid.sum() < 100:
        raise ValueError("Fewer than 100 valid pixels in the image (all NoData or cloud?).")

    step(0.35, "Choosing the threshold")
    if threshold is not None:
        t = {"threshold": float(threshold), "eta": None, "method": "user threshold", "tiles_used": 0}
    else:
        try:
            if method == "global":
                th, eta = C.otsu(idx[valid])
                t = {"threshold": th, "eta": eta, "method": "global Otsu", "tiles_used": 0}
            else:
                t = C.tile_threshold(idx, water_above)
        except C.ThresholdError as exc:
            raise ValueError(f"Automatic threshold failed: {exc}. Set a manual threshold "
                             "(e.g. 0 for NDWI/MNDWI) or use supervised mode.")
    info["threshold"] = t
    if t.get("eta") is not None:
        info["confidence"] = C.confidence_label(t["eta"])
        if t["eta"] < 0.65:
            notes.append(f"Low separability (eta = {t['eta']:.2f}; a single Gaussian population gives 0.64): the "
                         "histogram is not clearly two-class. Inspect the mask; consider a manual threshold or a "
                         "smaller area of interest.")
    info["sar_band"] = ("vv" if "vv" in scene.sources else "vh") if iname == "sar" else None
    lim = plausible_range(iname, scene, info["sar_band"])
    if threshold is None and lim and not (lim[0] <= t["threshold"] <= lim[1]):
        raise ValueError(
            f"The automatic {iname.upper()} threshold ({t['threshold']:.3f}) is outside the range where open "
            f"water is normally separated ({lim[0]:g} to {lim[1]:g}). The image probably contains little or "
            "no water at that date (for a flash flood the water may have drained before the overpass - use "
            "flood-trace mode with a pre-event image), or it is dominated by cloud or shadow. If you are "
            "sure, set the threshold manually.")

    water, valid = C.classify(idx, t["threshold"], water_above)
    raw_share = float(water.sum() / max(1, valid.sum()))
    info["water_share_raw"] = raw_share

    permanent = None
    if pre is not None and mode == "water":
        step(0.55, "Pre-event water (permanent water)")
        try:
            pname = choose_index(pre, iname if iname != "sar" else "auto")
            pidx, pgrid, pabove, _ = compute_index(pre, pname, None, notes)
            pv = np.isfinite(pidx)
            pt = C.tile_threshold(pidx, pabove) if threshold is None else {"threshold": float(threshold)}
            pw, pv = C.classify(pidx, pt["threshold"], pabove)
            permanent, _ = R.mask_to_grid(pw, pv, pgrid, grid, 0.5)
            info["permanent_threshold"] = pt["threshold"]
            info["permanent_share"] = float(permanent.sum() / max(1, valid.sum()))
        except ValueError as exc:
            notes.append(f"Pre-event image not used: {exc}")
            permanent = None

    step(0.7, "Cleaning the mask")
    cleaned, cstats = C.clean_mask(water, valid, min_water_cells, max_hole_cells)
    info["cleanup"] = cstats
    info["cleanup"].update({"min_water_cells": int(min_water_cells), "max_hole_cells": int(max_hole_cells)})
    info["water_share"] = float(cleaned.sum() / max(1, valid.sum()))
    if info["water_share"] < 0.001:
        notes.append("Less than 0.1 % of the valid image is classified as water. If a flood is expected, the "
                     "image may postdate it (water already drained - consider flood-trace mode with a "
                     "pre-event image) or the threshold is too strict.")
    if info["water_share"] > 0.6:
        notes.append("More than 60 % of the image is classified as water: check the threshold and the bands.")
    # threshold uncertainty: the extents at both ends of the near-optimal band
    incl = excl = None
    band = t.get("band")
    if threshold is None and band and band[1] > band[0]:
        step(0.85, "Threshold uncertainty band")
        lo_t, hi_t = (band[0], band[1]) if water_above else (band[1], band[0])
        w1, _ = C.classify(idx, lo_t, water_above)
        w2, _ = C.classify(idx, hi_t, water_above)
        incl, _ = C.clean_mask(w1, valid, min_water_cells, max_hole_cells)
        excl, _ = C.clean_mask(w2, valid, min_water_cells, max_hole_cells)
        info["threshold_band"] = [float(band[0]), float(band[1])]
        info["water_share_band"] = [float(excl.sum() / max(1, valid.sum())), float(incl.sum() / max(1, valid.sum()))]

    # flood trace: bare soil / fresh sediment after the event (erosion- or
    # deposition-dominated) versus vegetation damaged but still present
    bare = None
    if mode == "trace" and scene.has("swir1", "red", "nir", "blue"):
        sw = read_role(scene, "swir1", grid, notes)
        rd = read_role(scene, "red", grid, notes)
        nr = read_role(scene, "nir", grid, notes)
        bl = read_role(scene, "blue", grid, notes)
        bsi = C.normalized_difference(sw + rd, nr + bl)   # Rikimaru et al. 2002
        with np.errstate(invalid="ignore"):
            bare = cleaned & np.isfinite(bsi) & (bsi > 0.0)
        info["trace_bare_share"] = float(bare.sum() / max(1, cleaned.sum()))
    info["notes"] = notes
    step(1.0, "")
    return {"water": cleaned, "water_raw": water, "valid": valid, "grid": grid, "permanent": permanent,
            "cloud": cloud, "index": idx, "info": info, "scene": scene,
            "water_incl": incl, "water_excl": excl, "bare": bare}


MASK_COLORS = {0: (235, 228, 210, 255), 1: (30, 110, 200, 255), 2: (120, 80, 160, 255), 255: (0, 0, 0, 0)}


def mask_raster_values(water, valid, permanent=None):
    """uint8: 1 flood water, 2 permanent water, 0 dry, 255 NoData."""
    out = np.full(water.shape, 255, dtype="uint8")
    out[valid] = 0
    out[water & valid] = 1
    if permanent is not None:
        out[permanent & valid] = 2
    return out
