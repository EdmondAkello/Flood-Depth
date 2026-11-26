# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Google GEM 15 m DTM against the 30 m DEMs (not shipped in the plugin zip).

Same truth as dem_recovery.py (a water surface following the valley trend on
FABDEM, 10 % and 25 % of the area flooded); depth estimated on ALOS, TanDEM-X
and GEM (warped to the same 30 m grid). Also the method error on GEM at its
own 15 m grid (truth and estimate on GEM), and each DEM's offset to FABDEM.

    python -m flooddepth.benchmarks.dem_gem15 <dem_dir> [out.json]

dem_dir holds <site>_fabdem.tif, <site>_alos.tif, <site>_tandemx.tif and
<site>_gem15.tif (a clip of the GEM tile covering the site) for nyando and maimahiu.
"""

import json
import os
import sys
import time

import numpy as np

from ..core import raster as R, surface as SF, stats as ST
from .dem_recovery import SITES, truth


def _r(e):
    return {k: (round(v, 3) if isinstance(v, float) else v) for k, v in e.items()}


def run(dem_dir):
    ll = R.srs_from_epsg(4326).ExportToWkt()
    out = {}
    for site in ("nyando", "maimahiu"):
        bbox = SITES[site][0]
        path = {k: os.path.join(dem_dir, f"{site}_{k}.tif") for k in ("fabdem", "alos", "tandemx", "gem15")}
        grid, _ = R.analysis_grid(path["fabdem"], bbox, ll)
        arr = {k: R.read_on_grid(p, grid, 1, "bilinear", zero_nodata_auto=True)[0] for k, p in path.items()}
        offset = {}
        for k in ("alos", "tandemx", "gem15"):
            d = (arr[k] - arr["fabdem"])
            d = d[np.isfinite(d)]
            offset[k] = {"mean_m": round(float(d.mean()), 2), "sd_m": round(float(d.std()), 2)}
        rad = int(round(1500.0 / grid.cell_x))
        cases = []
        for share in (0.10, 0.25):
            wse_t, flood = truth(arr["fabdem"], share, rad)
            true_d = np.where(flood, wse_t - arr["fabdem"], np.nan)
            row = {"scenario": f"valley trend, {int(share * 100)} % flooded",
                   "true_mean_m": round(float(np.nanmean(true_d)), 2)}
            for k in ("alos", "tandemx", "gem15"):
                v = np.isfinite(arr[k])
                s_, _, _, _ = SF.estimate_surface(arr[k], flood & v, v, grid)
                d, _ = SF.depth_from_surface(s_, arr[k], flood & v)
                row[k] = _r(ST.error_stats(d, true_d))
            cases.append(row)
        g15g = R.grid_of(path["gem15"])
        g15, _ = R.analysis_grid(path["gem15"], R.transform_bounds(bbox, ll, g15g.wkt), g15g.wkt)
        z = R.read_on_grid(path["gem15"], g15, 1, "bilinear")[0]
        own = {}
        for share in (0.10, 0.25):
            wse_t, flood = truth(z, share, int(round(1500.0 / g15.cell_x)))
            true_d = np.where(flood, wse_t - z, np.nan)
            t = time.time()
            s_, _, _, _ = SF.estimate_surface(z, flood, np.isfinite(z), g15)
            d, _ = SF.depth_from_surface(s_, z, flood)
            e = _r(ST.error_stats(d, true_d))
            e["seconds"] = round(time.time() - t, 2)
            own[f"{int(share * 100)} % flooded"] = e
        out[site] = {"grid_30m": f"{grid.cols} x {grid.rows}", "grid_15m": f"{g15.cols} x {g15.rows}",
                     "offset_to_fabdem": offset, "truth_on_fabdem": cases, "gem15_own_grid": own}
    return out


if __name__ == "__main__":
    res = run(sys.argv[1])
    print(json.dumps(res, indent=1))
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w") as f:
            json.dump(res, f, indent=1)
