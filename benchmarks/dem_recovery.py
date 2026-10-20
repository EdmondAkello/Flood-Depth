# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Depth-recovery benchmark on real DEMs (not shipped in the plugin zip).

Truth: a water surface that follows the valley trend - the DEM smoothed
over ~1.5 km plus an offset h chosen so that a target share of the area is
flooded; extent = cells below it connected to the largest wet region;
true depth = surface - DEM. The engine sees only the extent and a DEM.

(a) same DEM   : method error alone (interpolation of the surface from
                 the edge), for several settings;
(b) other DEMs : truth on the bare-earth DEM (FABDEM), depth estimated
                 on ALOS / TanDEM-X / Copernicus - the DEM-choice error,
                 and whether the ensemble P10-P90 covers it.

    python -m flooddepth.benchmarks.dem_recovery <site> <dem_dir> <glo30_dir> [out.json]

dem_dir holds <site>_fabdem.tif, <site>_alos.tif, <site>_tandemx.tif clips;
glo30_dir holds <site>_copernicus.vrt and <site>_hand.vrt (GLO-30 DEM and
HAND tiles from the AWS open-data buckets).
"""

import json
import math
import os
import sys
import time

import numpy as np

from ..core import raster as R, surface as SF, classify as C, uncertainty as U, stats as ST

INTERP = os.environ.get("FD_INTERP", "idw")

SITES = {  # bbox lon0, lat0, lon1, lat1 ; DEM file stems in dem_dir
    "nyando": ((34.82, -0.33, 35.02, -0.13), ["fabdem", "alos", "tandemx"]),
    "maimahiu": ((36.55, -1.04, 36.68, -0.92), ["fabdem", "alos", "tandemx"]),
}


def smooth(z, radius_cells):
    m = np.isfinite(z)
    return SF.focal_mean(np.where(m, z, 0.0), m, radius_cells)


def truth(dem, share, radius_cells):
    trend = smooth(dem, radius_cells)
    rel = dem - trend
    v = rel[np.isfinite(rel)]
    h = float(np.percentile(v, 100 * share))
    wse = trend + h
    wet = np.isfinite(dem) & (dem < wse)
    lab, n = C.label_components(wet)
    if n:
        sizes = np.bincount(lab.ravel())
        sizes[0] = 0
        wet = lab == int(np.argmax(sizes))
    return wse.astype("float32"), wet


def truth_channel(dem, h, grid):
    """HAND-type inundation: water h m above the nearest valley-bottom cell.
    Valley bottom = lowest 3 % of DEM - 300 m local mean (connected to the
    largest bottom network). Extent = cells below that level connected to
    the bottom cells."""
    r = max(2, int(round(300.0 / grid.cell_x)))
    rel = dem - smooth(dem, r)
    v = rel[np.isfinite(rel)]
    bottom = np.isfinite(rel) & (rel < np.percentile(v, 3))
    wse, _ = SF.fill_nearest(np.where(bottom, dem, np.nan).astype("float32"), bottom, grid.cell_x, grid.cell_y)
    wse = wse + h
    wet = np.isfinite(dem) & (dem < wse)
    keep, _, _ = C.keep_connected(wet, bottom)
    return wse.astype("float32"), keep


def truth_hand(hand, dem, h):
    """HAND inundation (Nobre et al. 2011): water h m above the drainage.
    Wet = HAND < h, connected to the drainage (HAND < 0.5 m); the water
    surface is DEM + (h - HAND), true depth = h - HAND."""
    wet = np.isfinite(hand) & np.isfinite(dem) & (hand < h)
    chan = np.isfinite(hand) & (hand < 0.5)
    wet, _, _ = C.keep_connected(wet, chan)
    return (dem + (h - hand)).astype("float32"), wet


def run(site, dem_dir, hand_dir, members=12):
    bbox, dems = SITES[site]
    ll = R.srs_from_epsg(4326).ExportToWkt()
    paths = {k: os.path.join(dem_dir, f"{site}_{k}.tif") for k in dems}
    paths["copernicus"] = os.path.join(hand_dir, f"{site}_copernicus.vrt")
    grid, _ = R.analysis_grid(paths["fabdem"], bbox, ll)
    arr = {k: R.read_on_grid(p, grid, 1, "bilinear", zero_nodata_auto=True)[0] for k, p in paths.items()}
    hand = R.read_on_grid(os.path.join(hand_dir, f"{site}_hand.vrt"), grid, 1, "bilinear")[0]
    rad = int(round(1500.0 / grid.cell_x))
    out = {"site": site, "grid": f"{grid.cols}x{grid.rows} @ {grid.cell_x:.1f} m", "interp_default": INTERP,
           "cases": []}
    settings = {
        "nearest edge, smooth 2": {"interpolation": "nearest"},
        "nearest edge, no outlier screen": {"interpolation": "nearest", "outlier_k": 0},
        "inverse distance, smooth 2": {"interpolation": "idw"},
        "inverse distance, no smoothing": {"interpolation": "idw", "smooth_radius": 0},
        "inverse distance, no outlier screen": {"interpolation": "idw", "outlier_k": 0},
    }
    scenarios = [("valley trend, 10 % flooded", "fabdem", lambda: truth(arr["fabdem"], 0.10, rad)),
                 ("valley trend, 25 % flooded", "fabdem", lambda: truth(arr["fabdem"], 0.25, rad)),
                 ("HAND, water 1 m above drainage", "copernicus", lambda: truth_hand(hand, arr["copernicus"], 1.0)),
                 ("HAND, water 3 m above drainage", "copernicus", lambda: truth_hand(hand, arr["copernicus"], 3.0))]
    for name, tdem, make in scenarios:
        wse_t, flood = make()
        true_d = np.where(flood, wse_t - arr[tdem], np.nan)
        case = {"scenario": name, "truth_dem": tdem, "flood_km2": float(flood.sum() * grid.cell_area() / 1e6),
                "true_mean_m": float(np.nanmean(true_d)), "true_p90_m": float(np.nanpercentile(true_d, 90)),
                "same_dem": {}, "other_dem": {}}
        valid = np.isfinite(arr[tdem])
        for sname, kw in settings.items():
            t = time.time()
            s_, info, _, _ = SF.estimate_surface(arr[tdem], flood, valid, grid, **kw)
            d, _ = SF.depth_from_surface(s_, arr[tdem], flood)
            e = ST.error_stats(d, true_d)
            e["seconds"] = round(time.time() - t, 2)
            case["same_dem"][sname] = e
        for k, dem in arr.items():
            if k == tdem:
                continue
            v = np.isfinite(dem)
            s_, info, _, _ = SF.estimate_surface(dem, flood & v, v, grid, interpolation=INTERP)
            d, _ = SF.depth_from_surface(s_, dem, flood & v)
            e = ST.error_stats(d, true_d)
            sig, lab = U.dem_sigma_from_name(paths[k])
            full, us = U.run_ensemble(dem, v, grid, flood & v, members=members, sigma_dem=sig, corr_m=150.0,
                                      edge_shift=1, seed=1, surface_kwargs={"interpolation": INTERP})
            inside = flood & v & np.isfinite(full["p10"]) & np.isfinite(true_d)
            cov = float(np.mean((full["p10"][inside] <= true_d[inside]) & (true_d[inside] <= full["p90"][inside])))
            e.update({"sigma_m": sig, "p10_p90_coverage": cov,
                      "mean_band_m": float(np.nanmean((full["p90"] - full["p10"])[inside])),
                      "method_sd_m": us.get("method_sd_m")})
            case["other_dem"][k] = e
        out["cases"].append(case)
    return out


if __name__ == "__main__":
    site, dem_dir, hand_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    res = run(site, dem_dir, hand_dir)
    print(json.dumps(res, indent=1))
    if len(sys.argv) > 4:
        with open(sys.argv[4], "w") as f:
            json.dump(res, f, indent=1)
