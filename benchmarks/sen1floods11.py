# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Classification benchmark on Sen1Floods11 hand-labelled chips (not shipped).

Bonafilia et al. (2020), Sen1Floods11: 446 hand-labelled 512 x 512 chips of
11 flood events (Bolivia, Ghana, India, Mekong, Nigeria, Pakistan, Paraguay,
Somalia, Spain, Sri Lanka, USA) with co-registered Sentinel-1 (VV, VH dB)
and Sentinel-2 L1C (13 bands). Labels: 1 water, 0 not water, -1 no data.
Public: gs://sen1floods11/v1.1/data/flood_events/HandLabeled/.

Each chip goes through flooddepth.core.imagery.classify_image exactly as in
the plugin (automatic index, tile-based Otsu, plausibility guard, cleanup of
patches < 10 cells). Band numbers are given because the chips carry no band
names (S2: green 3, NIR 8, SWIR1 12, blue 2; S1: VV 1, VH 2).

    python -m flooddepth.benchmarks.sen1floods11 <chip_dir> [out.json]
"""

import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

from ..core import imagery as I, raster as R, stats as ST

CONFIGS = {
    "S2 fully automatic (edge Otsu)": ("S2Hand", {}, "auto", "edge"),
    "S1 fully automatic (edge Otsu, speckle 5x5)": ("S1Hand", {}, "auto", "edge"),
    "S2 NDWI edge Otsu": ("S2Hand", {"green": 3, "nir": 8, "blue": 2}, "ndwi", "edge"),
    "S2 MNDWI tile-Otsu": ("S2Hand", {"green": 3, "nir": 8, "swir1": 12, "blue": 2, "red": 4}, "mndwi", "tile"),
    "S2 NDWI tile-Otsu": ("S2Hand", {"green": 3, "nir": 8, "swir1": 12, "blue": 2, "red": 4}, "ndwi", "tile"),
    "S2 MNDWI global Otsu": ("S2Hand", {"green": 3, "nir": 8, "swir1": 12, "blue": 2, "red": 4}, "mndwi", "global"),
    "S1 VV tile-Otsu": ("S1Hand", {"vv": 1, "vh": 2}, "sar", "tile"),
    "S1 VV global Otsu": ("S1Hand", {"vv": 1, "vh": 2}, "sar", "global"),
}


def run(chip_dir):
    labels = sorted(glob.glob(os.path.join(chip_dir, "*_LabelHand.tif")))
    out = {}
    for name, (kind, ov, index, method) in CONFIGS.items():
        pooled = defaultdict(lambda: np.zeros(4, dtype="int64"))   # tp fp fn tn per country
        per_chip, refused, failed = [], [], []
        for lp in labels:
            stem = os.path.basename(lp).replace("_LabelHand.tif", "")
            country = stem.split("_")[0]
            ip = os.path.join(chip_dir, f"{stem}_{kind}.tif")
            lab, lg, _ = R.read_band(lp)
            valid_l = np.isfinite(lab) & (lab >= 0)
            truth_share = float((lab[valid_l] == 1).mean()) if valid_l.any() else 0.0
            try:
                res = I.classify_image([ip], ov, index, None, method, 10, 0, 0.25 if kind == "S2Hand" else 0)
            except ValueError as exc:
                msg = str(exc)
                (refused if ("outside the range" in msg or "No water was found" in msg) else failed).append(
                    {"chip": stem, "truth_water_share": truth_share, "why": msg[:80]})
                # a refusal is scored as "no water mapped"
                w = np.zeros(lab.shape, bool)
                v = valid_l
            else:
                w, v = res["water"], res["valid"] & valid_l
            s = ST.extent_agreement(w, lab == 1, v)
            pooled[country] += [s["tp"], s["fp"], s["fn"], s["tn"]]
            per_chip.append({"chip": stem, "csi": s["csi"], "truth_water_share": truth_share})
        tot = sum(pooled.values())
        res = {"pooled": _scores(tot), "by_country": {c: _scores(v) for c, v in sorted(pooled.items())},
               "chips": len(labels), "refused": len(refused), "failed": len(failed),
               "refused_truth_water_share_median": float(np.median([r["truth_water_share"] for r in refused]))
               if refused else None,
               "refused_with_over_5pct_water": int(sum(r["truth_water_share"] > 0.05 for r in refused)),
               "median_chip_csi": float(np.nanmedian([c["csi"] for c in per_chip]))}
        out[name] = res
        print(f"{name:44s} pooled CSI {res['pooled']['csi']:.3f}  F1 {res['pooled']['f1']:.3f}  "
              f"hit {res['pooled']['hit_rate']:.3f}  FAR {res['pooled']['false_alarm_ratio']:.3f}  "
              f"refused {len(refused)} (with >5% water: {res['refused_with_over_5pct_water']}), failed {len(failed)}")
    return out


def _scores(v):
    tp, fp, fn, tn = [int(x) for x in v]
    d = lambda a, b: float(a) / b if b else float("nan")
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "csi": d(tp, tp + fp + fn), "f1": d(2 * tp, 2 * tp + fp + fn),
            "hit_rate": d(tp, tp + fn), "false_alarm_ratio": d(fp, tp + fp), "bias": d(tp + fp, tp + fn)}


if __name__ == "__main__":
    res = run(sys.argv[1])
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w") as f:
            json.dump(res, f, indent=1)
