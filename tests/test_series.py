# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Several image dates -> per-date runs, maximum depth, dates flooded.

    python -m flooddepth.tests.test_series
"""

import csv
import os
import tempfile

import numpy as np

from .common import check, finish
from . import synth as S
from ..core import raster as R
from ..core.pipeline import QuickOptions
from ..core.series import group_by_date, run_series


def main():
    d = tempfile.mkdtemp(prefix="fd_series_")
    dem = S.write_dem(os.path.join(d, "dem_fabdem.tif"))
    a = S.write_s2_raw(os.path.join(d, "img"), "2024-05-05", True, seed=1)
    b = S.write_s2_raw(os.path.join(d, "img"), "2024-05-12", True, seed=3)
    c = S.write_s2_raw(os.path.join(d, "img"), "2024-05-20", False, seed=4)   # drained: only the lake
    pre = S.write_s2_raw(os.path.join(d, "pre"), "2024-04-07", False, seed=2)

    print("\n1. Grouping by date")
    g = group_by_date(a + b + c)
    check("three dates, in order", list(g) == ["2024-05-05", "2024-05-12", "2024-05-20"])
    check("files stay with their date", all(len(v) == len(a) for v in g.values()))
    try:
        group_by_date(a + [os.path.join(d, "undated_B03.tif")])
        check("an undated file is refused with its name", False)
    except ValueError as exc:
        check("an undated file is refused with its name", "undated_B03.tif" in str(exc))

    print("\n2. The series")
    out = os.path.join(d, "series")
    res = run_series(QuickOptions(dem_path=dem, out_dir=out, image_paths=a + b + c, pre_paths=pre,
                                  uncertainty_members=4), log=None)
    st = {r["date"]: r["status"] for r in res["rows"]}
    check("flooded dates mapped", st["2024-05-05"] == "ok" and st["2024-05-12"] == "ok", str(st))
    check("drained date skipped with a reason, the series continues", st["2024-05-20"].startswith("skipped"),
          st["2024-05-20"][:80])
    for k in ("max_depth", "dates_flooded", "first_date_flooded", "csv", "html"):
        check(f"{k} written", os.path.exists(res["outputs"][k]))
    mx, gm, _ = R.read_band(res["outputs"]["max_depth"])
    cnt, _, _ = R.read_band(res["outputs"]["dates_flooded"])
    d1, g1, _ = R.read_band(os.path.join(out, "2024-05-05", "depth.tif"))
    d1u = R.warp_array(d1, g1, gm, "near")
    both = np.isfinite(d1u) & np.isfinite(mx)
    check("maximum depth >= each date's depth", bool(np.all(mx[both] >= d1u[both] - 1e-5)))
    check("dates flooded is 1 or 2 where flooded", set(np.unique(cnt[np.isfinite(mx)]).tolist()) <= {1.0, 2.0})
    with open(res["outputs"]["csv"], newline="") as fh:
        rows = list(csv.DictReader(fh))
    check("CSV has a row per date", [r["date"] for r in rows] == list(g))
    html = open(res["outputs"]["html"], encoding="utf-8").read()
    check("summary links each date's report", "2024-05-05/report.html" in html and "2024-05-12/report.html" in html)
    check("envelope area >= each date's area", res["envelope"]["area_km2"] + 1e-9 >=
          max(r["area_km2"] for r in res["rows"] if r["status"] == "ok"))
    finish("test_series")


if __name__ == "__main__":
    main()
