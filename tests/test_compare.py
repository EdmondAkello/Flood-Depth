# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Extent agreement scores and depth statistics.

    python -m flooddepth.tests.test_compare
"""

import numpy as np

from .common import check, finish
from ..core import stats as ST
from ..core.compare import compare_extents
from . import synth as S


def main():
    print("\n1. Contingency scores")
    pred = np.array([1, 1, 1, 0, 0, 0, 1, 0], bool)
    ref = np.array([1, 1, 0, 1, 0, 0, 1, 0], bool)
    s = ST.extent_agreement(pred, ref, np.ones_like(pred))
    check("hits 3, false alarms 1, misses 1", (s["tp"], s["fp"], s["fn"], s["tn"]) == (3, 1, 1, 3))
    check("CSI = 3/5", abs(s["csi"] - 0.6) < 1e-12)
    check("hit rate 0.75, FAR 0.25, bias 1, F1 0.75",
          (round(s["hit_rate"], 3), round(s["false_alarm_ratio"], 3), round(s["bias"], 3), round(s["f1"], 3))
          == (0.75, 0.25, 1.0, 0.75))
    e = ST.error_stats(np.array([1.0, 2.0, np.nan]), np.array([0.5, 2.5, 1.0]))
    check("bias 0, MAE 0.5, RMSE 0.5", e["cells"] == 2 and abs(e["bias_m"]) < 1e-12 and abs(e["mae_m"] - 0.5) < 1e-12)

    print("\n2. Depth statistics")
    d = np.array([[0.2, 0.7, 1.5], [2.5, 3.5, np.nan]])
    st = ST.depth_stats(d, 100.0)
    check("5 cells, 500 m2", st["cells"] == 5 and abs(st["area_km2"] - 0.0005) < 1e-12)
    check("one cell in each class", [c["cells"] for c in st["classes"]] == [1, 1, 1, 1, 1])
    check("volume = sum x cell area", abs(st["volume_m3"] - 840.0) < 1e-9)

    print("\n3. Polygons vs polygons on a UTM grid")
    w = S.flood_polygon_wkt()
    s, ag, g = compare_extents({"wkts": [w], "crs": S.utm_wkt()}, {"wkts": [w], "crs": S.utm_wkt()}, cell=20)
    check("identical extents: CSI 1", s["csi"] == 1.0 and s["fp"] == 0 and s["fn"] == 0)
    check("agreement raster codes 0/1 only", set(np.unique(ag).tolist()) <= {0, 1})


if __name__ == "__main__":
    main()
    finish("test_compare")
