# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Quick flood depth: image (or extent) + DEM -> depth, in one call.

This is the whole workflow without QGIS, so it can be tested and scripted:

    from flooddepth.core.pipeline import QuickOptions, run_quick
    res = run_quick(QuickOptions(image_paths=[...], dem_path="dem.tif",
                                 out_dir="run1"))

Output folder (fixed layout):
    water_mask.tif     classification on the image grid
                       (1 flood water, 2 permanent water, 0 dry, 255 NoData)
    flood_extent.gpkg  flood_extent (polygons), waterline (points used and
                       dropped) - on the analysis grid
    wse.tif            water-surface elevation (DEM method) or water height
                       above drainage (HAND method)
    depth.tif          flood depth (m); NoData outside the flood
    report.html        inputs, every automatic decision, statistics, limits
    run.json           the same as machine-readable JSON
"""

import json
import os
import time
from dataclasses import dataclass, field, asdict

import numpy as np

from . import raster as R
from . import imagery as I
from . import classify as C
from . import surface as SF
from . import stats as ST
from . import vector as V
from . import uncertainty as U

VERSION = "0.4.1"


@dataclass
class QuickOptions:
    dem_path: str
    out_dir: str
    image_paths: list = field(default_factory=list)
    pre_paths: list = field(default_factory=list)
    mode: str = "water"                 # water | trace
    index: str = "auto"                 # auto | ndwi | mndwi | sar
    threshold: float = None             # None = automatic
    threshold_method: str = "edge"      # edge | tile | global (benchmarks/sen1floods11.py)
    band_overrides: dict = field(default_factory=dict)
    extent_wkts: list = field(default_factory=list)     # flood polygons instead of an image
    extent_crs: str = ""
    aoi_wkts: list = field(default_factory=list)
    aoi_crs: str = ""
    seed_wkts: list = field(default_factory=list)       # keep flood connected to these
    seed_crs: str = ""
    seed_buffer_m: float = 30.0
    permanent_wkts: list = field(default_factory=list)  # known permanent water
    permanent_crs: str = ""
    permanent_raster: str = ""          # WorldCover, JRC occurrence or any water mask
    permanent_values: str = "auto"      # auto | comma list of class values | >=N
    method: str = "dem"                 # dem | hand
    hand_path: str = ""
    cell: float = 0.0
    min_water_cells: int = 10
    max_hole_cells: int = 0
    cloud_threshold: float = 0.25
    outlier_radius: int = 10
    outlier_k: float = 3.0
    outlier_min: float = 0.5
    smooth_radius: int = 2
    smooth_tolerance: float = 0.5
    edge_step: float = 0.5              # water level at the edge: fraction of the step to the first dry cell
    edge_step_cap: float = 1.0
    steep_slope_pct: float = 10.0
    interpolation: str = "idw"          # idw | nearest (benchmarks/dem_recovery.py)
    min_fraction: float = 0.5
    max_cells: int = 200_000_000
    uncertainty_members: int = 20      # 0 = no ensemble
    dem_sigma: float = 0.0              # 0 = from the DEM file name (preset)
    dem_corr_m: float = 150.0
    erosion_sigma: float = 1.0
    seed: int = 42


def _progress(cb, lo, hi):
    def f(frac, msg=""):
        if cb:
            cb(lo + (hi - lo) * frac, msg)
    return f


def run_quick(opt, progress=None, log=print):
    t0 = time.time()
    os.makedirs(opt.out_dir, exist_ok=True)
    info = {"version": VERSION, "started": time.strftime("%Y-%m-%d %H:%M:%S"), "options": _opt_dict(opt),
            "warnings": [], "steps": []}
    W = info["warnings"]

    def say(msg):
        info["steps"].append(msg)
        if log:
            log(msg)

    if not opt.image_paths and not opt.extent_wkts:
        raise ValueError("Give either satellite image file(s) or a flood-extent layer.")
    if opt.method == "hand" and not opt.hand_path:
        raise ValueError("The HAND method needs a HAND raster.")

    # 1. classification on the image's own grid
    cls = None
    if opt.image_paths:
        cls = I.classify_image(opt.image_paths, opt.band_overrides, opt.index, opt.threshold,
                               opt.threshold_method, opt.min_water_cells, opt.max_hole_cells,
                               opt.cloud_threshold, opt.pre_paths or None, None, opt.mode,
                               progress=_progress(progress, 0.0, 0.35))
        info["classification"] = cls["info"]
        W.extend(cls["info"].get("notes", []))
        t = cls["info"]["threshold"]
        say(f"Classified {cls['info']['index']} with {t['method']} threshold {t['threshold']:.3f}"
            + (f" (eta {t['eta']:.2f}, {cls['info'].get('confidence')})" if t.get("eta") is not None else "")
            + f"; water {100 * cls['info']['water_share']:.2f} % of valid pixels.")
        mpath = os.path.join(opt.out_dir, "water_mask.tif")
        R.write_raster(mpath, I.mask_raster_values(cls["water"], cls["valid"], cls["permanent"]), cls["grid"],
                       nodata=255, dtype="Byte", color_table=I.MASK_COLORS,
                       tags={"FLOODDEPTH": VERSION, "VALUES": "1 flood water, 2 permanent water, 0 dry, 255 NoData",
                             "INDEX": cls["info"]["index"], "THRESHOLD": t["threshold"]})
        info["water_mask"] = mpath

    # 2. analysis grid
    if cls is not None:
        ext, ext_wkt = cls["grid"].bounds, cls["grid"].wkt
    else:
        ext, ext_wkt = _wkts_bounds(opt.extent_wkts, opt.extent_crs, pad_frac=0.1)
    if opt.aoi_wkts:
        ab, aw = _wkts_bounds(opt.aoi_wkts, opt.aoi_crs)
        ab = R.transform_bounds(ab, aw, ext_wkt)
        ext = R.intersect(ext, ab)
        if ext is None:
            raise ValueError("The area of interest does not overlap the image / flood extent.")
    grid, gnotes = R.analysis_grid(opt.dem_path, ext, ext_wkt, opt.cell, opt.max_cells)
    for n in gnotes:
        say(n)
    info["grid"] = {"cols": grid.cols, "rows": grid.rows, "cell_x": grid.cell_x, "cell_y": grid.cell_y,
                    "crs": R.srs_from_wkt(grid.wkt).GetName()}
    if R.is_geographic(grid.wkt):
        raise ValueError("Internal error: the analysis grid is geographic.")
    prog = _progress(progress, 0.35, 0.95)
    prog(0.0, "Reading the DEM")
    dem, dnotes = R.read_on_grid(opt.dem_path, grid, 1, "bilinear", zero_nodata_auto=True)
    W.extend("DEM: " + n for n in dnotes)
    valid = np.isfinite(dem)
    if valid.mean() < 0.5:
        W.append(f"The DEM covers only {100 * valid.mean():.0f} % of the analysis area.")

    # 3. flood (and permanent water) on the analysis grid
    other = None
    if cls is not None:
        flood_src = cls["water"] if cls["permanent"] is None else (cls["water"] & ~cls["permanent"])
        flood, img_valid = R.mask_to_grid(flood_src, cls["valid"], cls["grid"], grid, opt.min_fraction)
        if cls["permanent"] is not None:
            other, _ = R.mask_to_grid(cls["permanent"], cls["valid"], cls["grid"], grid, opt.min_fraction)
        valid &= img_valid
    else:
        flood = V.rasterize_wkts(opt.extent_wkts, opt.extent_crs, grid)
    if opt.permanent_wkts:
        pw = V.rasterize_wkts(opt.permanent_wkts, opt.permanent_crs, grid)
        other = pw if other is None else (other | pw)
    if opt.permanent_raster:
        pr, _ = R.read_on_grid(opt.permanent_raster, grid, 1, "near")
        pw, rule = permanent_from_raster(pr, opt.permanent_values)
        say(f"Permanent water from {os.path.basename(opt.permanent_raster)}: {rule}; "
            f"{int(pw.sum()):,} cells ({pw.sum() * grid.cell_area() / 1e6:.2f} km2).")
        other = pw if other is None else (other | pw)
    if other is not None:
        flood &= ~other
        info["permanent_cells"] = int(other.sum())
    if opt.aoi_wkts:
        valid &= V.rasterize_wkts(opt.aoi_wkts, opt.aoi_crs, grid)
    flood &= valid
    if opt.seed_wkts:
        buf = max(1, int(round(opt.seed_buffer_m / grid.cell_x)))
        seed = C.dilate(V.rasterize_wkts(opt.seed_wkts, opt.seed_crs, grid, all_touched=True), buf)
        flood, kept, dropped = C.keep_connected(flood, seed)
        say(f"Connectivity: kept {kept} flood region(s) touching the reference features, dropped {dropped}.")
        info["connectivity"] = {"kept": kept, "dropped": dropped}
    else:
        seed = None
    if flood.sum() == 0:
        raise ValueError("No flood cells remain on the analysis grid (after masking, permanent-water removal "
                         "and connectivity). Check the inputs and the area of interest.")
    say(f"Flood extent: {int(flood.sum()):,} cells = {flood.sum() * grid.cell_area() / 1e6:.3f} km2.")

    # 4. surface and depth
    hand = None
    if opt.method == "hand":
        prog(0.15, "Reading HAND")
        hand, hn = R.read_on_grid(opt.hand_path, grid, 1, "bilinear", zero_nodata_auto=False)
        W.extend("HAND: " + n for n in hn)
        hvalid = valid & np.isfinite(hand)
        surf, sinfo, used, dropped = SF.estimate_surface(
            hand, flood, hvalid, grid, other_water=other, outlier_radius=opt.outlier_radius,
            outlier_k=opt.outlier_k, outlier_min=opt.outlier_min, smooth_radius=opt.smooth_radius,
            smooth_tolerance=opt.smooth_tolerance, interpolation=opt.interpolation,
            steep_slope_pct=opt.steep_slope_pct, edge_step=opt.edge_step, edge_step_cap=opt.edge_step_cap,
            progress=_progress(prog, 0.2, 0.8))
        depth, neg = SF.depth_from_surface(surf, hand, flood & hvalid)
        base = hand
        surf_label = "Water height above nearest drainage (m)"
    else:
        surf, sinfo, used, dropped = SF.estimate_surface(
            dem, flood, valid, grid, other_water=other, outlier_radius=opt.outlier_radius,
            outlier_k=opt.outlier_k, outlier_min=opt.outlier_min, smooth_radius=opt.smooth_radius,
            smooth_tolerance=opt.smooth_tolerance, interpolation=opt.interpolation,
            steep_slope_pct=opt.steep_slope_pct, edge_step=opt.edge_step, edge_step_cap=opt.edge_step_cap,
            progress=_progress(prog, 0.2, 0.8))
        depth, neg = SF.depth_from_surface(surf, dem, flood)
        base = dem
        surf_label = "Water-surface elevation (DEM units)"
    info["surface"] = sinfo
    info["surface"]["zero_clipped_cells"] = neg
    say(f"Waterline: {sinfo['waterline_cells']:,} cells, {sinfo['waterline_used']:,} used, "
        f"{sinfo['waterline_dropped']:,} dropped as outliers; {sinfo['edge_cells_excluded']:,} flood-edge cells "
        "at NoData / grid edge not used.")
    if sinfo["waterline_cells"] and sinfo["waterline_dropped"] / sinfo["waterline_cells"] > 0.25:
        W.append("More than 25 % of the waterline was dropped as outliers: the DEM is noisy along the flood "
                 "edge (DSM canopy or buildings?) or the extent is poor.")
    if sinfo.get("flood_steep_share", 0) > 0.2:
        W.append(f"{100 * sinfo['flood_steep_share']:.0f} % of the flood cells lie on ground steeper than "
                 f"{opt.steep_slope_pct:g} %: water does not pond there. In a flash flood or debris flow these "
                 "are flow paths; depths there are the least certain (see the uncertainty bands).")
    if neg / max(1, flood.sum()) > 0.1:
        W.append(f"{100 * neg / flood.sum():.0f} % of flood cells lie above the estimated water surface "
                 "(depth set to 0): the extent includes high ground or the DEM is a surface model there.")

    # 5. evidence classes and uncertainty
    slope = SF.slope_percent(dem, grid.cell_x, grid.cell_y)
    bare_g = None
    if cls is not None and cls.get("bare") is not None:
        bare_g, _ = R.mask_to_grid(cls["bare"], cls["valid"], cls["grid"], grid, 0.5)
    source = "polygons" if cls is None else ("trace" if opt.mode == "trace" else "image")
    evid = U.evidence_classes(flood, source, bare_g, slope, opt.steep_slope_pct)
    ca = grid.cell_area()
    unc = None
    if opt.uncertainty_members and opt.uncertainty_members >= 2:
        sig, siglab = (opt.dem_sigma, "user value") if opt.dem_sigma > 0 else U.dem_sigma_from_name(opt.dem_path)
        if opt.dem_sigma <= 0 and siglab.startswith("unknown"):
            W.append(f"DEM vertical error not known from the file name; {sig:g} m (1 sigma) assumed for the "
                     "uncertainty bands. Set it if you know it (LiDAR ~0.15 m, FABDEM ~1.5 m).")
        fracs = None
        if cls is not None and cls.get("water_incl") is not None:
            perm = cls["permanent"] if cls["permanent"] is not None else np.zeros_like(cls["water"])
            fracs = []
            for mk in (cls["water_excl"], cls["water"], cls["water_incl"]):
                fr = R.warp_array(np.where(cls["valid"], (mk & ~perm).astype("float32"), np.nan),
                                  cls["grid"], grid, "average")
                fracs.append(np.where(valid, fr, np.nan))
        kw = dict(outlier_radius=opt.outlier_radius, outlier_k=opt.outlier_k, outlier_min=opt.outlier_min,
                  smooth_radius=opt.smooth_radius, smooth_tolerance=opt.smooth_tolerance,
                  interpolation=opt.interpolation, steep_slope_pct=opt.steep_slope_pct,
                  edge_step=opt.edge_step, edge_step_cap=opt.edge_step_cap)
        try:
            unc, usum = U.run_ensemble(
                dem, valid if hand is None else (valid & np.isfinite(hand)), grid, flood,
                members=int(opt.uncertainty_members), sigma_dem=sig, corr_m=opt.dem_corr_m,
                sigma_erosion=opt.erosion_sigma, erosion=(evid == 3), frac_options=fracs,
                edge_shift=2 if source == "trace" else 1, other_water=other, seed=opt.seed,
                surface_kwargs=kw, depth_base=hand, seed_mask=seed,
                progress=_progress(prog, 0.8, 0.95))
            usum["sigma_source"] = siglab
            a = usum["area_km2_p10_p50_p90"]
            if a[0] > 0 and a[2] / a[0] > 3.0:
                W.append(f"The flooded area changes {a[2] / a[0]:.0f}-fold across the ensemble (P10 {a[0]:.2f} to "
                         f"P90 {a[2]:.2f} km2): the extent, not the DEM, dominates the uncertainty. The "
                         "classification threshold is weakly constrained (see the threshold band); review the "
                         "mask, give an area of interest, or set the threshold manually.")
            info["uncertainty"] = usum
            say(f"Uncertainty: {usum['succeeded']} members; area {usum['area_km2_p10_p50_p90'][0]:.2f}-"
                f"{usum['area_km2_p10_p50_p90'][2]:.2f} km2, mean depth {usum['mean_depth_m_p10_p50_p90'][0]:.2f}-"
                f"{usum['mean_depth_m_p10_p50_p90'][2]:.2f} m (P10-P90).")
        except ValueError as exc:
            W.append(f"Uncertainty ensemble not computed: {exc}")
    rows = []
    for code in sorted(set(np.unique(evid).tolist()) - {255}):
        sel = evid == code
        n = int(sel.sum())
        row = {"code": code, "label": U.EVIDENCE[code], "cells": n, "area_km2": n * ca / 1e6,
               "mean_depth_m": float(np.nanmean(depth[sel])) if n else float("nan")}
        if unc is not None:
            row["mean_band_m"] = float(np.nanmean((unc["p90"] - unc["p10"])[sel]))
            row["mean_prob"] = float(np.nanmean(unc["prob"][sel]))
        rows.append(row)
    info["evidence"] = rows
    if source == "trace":
        e3 = sum(r["cells"] for r in rows if r["code"] == 3)
        W.append(f"Flood-trace extent: {100 * e3 / max(1, flood.sum()):.0f} % of it is bare soil or fresh sediment "
                 "after the event (erosion / deposition). There the trace edge is the erosion width and the "
                 "pre-event DEM no longer describes the ground; treat those depths as indicative.")

    prog(0.95, "Writing outputs")
    info["depth"] = ST.depth_stats(depth, ca)
    # interior: the flood without its waterline ring (depth there is ~0 by construction);
    # in fragmented extents the all-cell statistics are dominated by edges
    interior = flood & ~(used | dropped)
    info["depth_interior"] = ST.depth_stats(np.where(interior, depth, np.nan), ca)
    info["depth_interior"]["share_of_flood"] = float(interior.sum() / max(1, flood.sum()))
    wpath = os.path.join(opt.out_dir, "wse.tif")
    dpath = os.path.join(opt.out_dir, "depth.tif")
    tags = {"FLOODDEPTH": VERSION, "METHOD": opt.method}
    R.write_raster(wpath, surf, grid, -9999.0, "Float32", dict(tags, CONTENT=surf_label))
    R.write_raster(dpath, depth, grid, -9999.0, "Float32", dict(tags, CONTENT="flood depth (m)"))
    gpkg = os.path.join(opt.out_dir, "flood_extent.gpkg")
    if os.path.exists(gpkg):
        os.remove(gpkg)
    npoly = V.polygonize(flood, grid, gpkg, "flood_extent",
                         {"source": "image" if cls is not None else "polygons"})
    rr, cc = np.nonzero(used | dropped)
    if rr.size > 200000:
        sel = np.linspace(0, rr.size - 1, 200000).astype(int)
        rr, cc = rr[sel], cc[sel]
    V.write_points(gpkg, "waterline", grid, rr, cc,
                   {"used": used[rr, cc].astype("int32"), "z": base[rr, cc].astype("float64"),
                    "surface": np.nan_to_num(surf[rr, cc]).astype("float64")})
    epath = os.path.join(opt.out_dir, "evidence.tif")
    R.write_raster(epath, evid, grid, 255, "Byte", dict(tags, VALUES="; ".join(f"{k} {v}" for k, v in U.EVIDENCE.items())),
                   color_table=U.EVIDENCE_COLORS)
    outs = {"wse": wpath, "depth": dpath, "extent": gpkg, "water_mask": info.get("water_mask"), "evidence": epath}
    if unc is not None:
        for key, name, lab in (("prob", "flood_probability.tif", "share of ensemble members flooded"),
                               ("p10", "depth_p10.tif", "depth 10th percentile (m)"),
                               ("p50", "depth_p50.tif", "depth median (m)"),
                               ("p90", "depth_p90.tif", "depth 90th percentile (m)")):
            pth = os.path.join(opt.out_dir, name)
            R.write_raster(pth, unc[key], grid, -9999.0, "Float32", dict(tags, CONTENT=lab))
            outs[key] = pth
    info.update({"outputs": outs,
                 "polygons": npoly, "seconds": round(time.time() - t0, 1)})
    from .report import write_report
    rpath = os.path.join(opt.out_dir, "report.html")
    write_report(rpath, info, depth=depth, flood=flood, dem=dem, grid=grid)
    info["outputs"]["report"] = rpath
    with open(os.path.join(opt.out_dir, "run.json"), "w", encoding="utf-8") as f:
        json.dump(_jsonable(info), f, indent=2)
    prog(1.0, "Done")
    say(f"Done in {info['seconds']} s.")
    return info


WORLDCOVER_CODES = {10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100}


def permanent_from_raster(arr, rule="auto"):
    """Permanent-water mask from a raster. Returns (mask, rule used).

    auto: ESA WorldCover (all values are WorldCover codes) -> class 80;
          values within 0-100 with many levels -> JRC Global Surface Water
          occurrence (%) >= 50; otherwise any value > 0.
    '80' or '80,90': those class values; '>=75': a threshold.
    """
    a = np.asarray(arr, dtype="float64")
    fin = np.isfinite(a)
    rule = (rule or "auto").strip().lower()
    if rule == "auto":
        vals = set(np.unique(a[fin]).astype(int).tolist()) if fin.any() else set()
        if vals and vals <= (WORLDCOVER_CODES | {0}):
            return fin & (a == 80), "ESA WorldCover class 80 (permanent water bodies)"
        if fin.any() and a[fin].min() >= 0 and a[fin].max() <= 100 and len(vals) > 12:
            return fin & (a >= 50), "occurrence >= 50 % (JRC Global Surface Water convention)"
        return fin & (a > 0), "any value > 0"
    if rule.startswith(">="):
        t = float(rule[2:])
        return fin & (a >= t), f"value >= {t:g}"
    codes = [float(x) for x in rule.replace(";", ",").split(",") if x.strip()]
    return fin & np.isin(a, codes), "values " + ", ".join(f"{c:g}" for c in codes)


def _wkts_bounds(wkts, crs, pad_frac=0.0):
    from osgeo import ogr
    xs0, ys0, xs1, ys1 = [], [], [], []
    for w in wkts:
        g = ogr.CreateGeometryFromWkt(w)
        if g is None or g.IsEmpty():
            continue
        a, b, c, d = g.GetEnvelope()
        xs0.append(a)
        xs1.append(b)
        ys0.append(c)
        ys1.append(d)
    if not xs0:
        raise ValueError("The layer has no usable geometry.")
    b = [min(xs0), min(ys0), max(xs1), max(ys1)]
    px, py = (b[2] - b[0]) * pad_frac, (b[3] - b[1]) * pad_frac
    return (b[0] - px, b[1] - py, b[2] + px, b[3] + py), crs


def _opt_dict(opt):
    d = asdict(opt)
    for k in ("extent_wkts", "aoi_wkts", "seed_wkts", "permanent_wkts"):
        d[k] = f"{len(d[k])} geometries" if d[k] else None
    for k in ("extent_crs", "aoi_crs", "seed_crs", "permanent_crs"):
        d[k] = R.srs_from_wkt(d[k]).GetName() if d[k] else None
    return d


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o
