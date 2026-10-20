# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Depth statistics and flood-extent agreement scores."""

import numpy as np

DEPTH_CLASSES = ((0.0, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, None))


def depth_stats(depth, cell_area_m2, classes=DEPTH_CLASSES):
    d = np.asarray(depth, dtype="float64")
    v = d[np.isfinite(d)]
    out = {"cells": int(v.size), "area_km2": v.size * cell_area_m2 / 1e6}
    if v.size == 0:
        return out
    out.update({
        "mean_m": float(v.mean()), "median_m": float(np.median(v)),
        "p90_m": float(np.percentile(v, 90)), "p99_m": float(np.percentile(v, 99)),
        "max_m": float(v.max()), "volume_m3": float(v.sum() * cell_area_m2),
        "zero_depth_share": float(np.mean(v <= 0.0)),
    })
    rows = []
    for lo, hi in classes:
        sel = (v > lo if lo > 0 else v >= lo) & ((v <= hi) if hi is not None else True)
        n = int(np.sum(sel))
        rows.append({"label": f"{lo:g}-{hi:g} m" if hi is not None else f"> {lo:g} m",
                     "cells": n, "area_km2": n * cell_area_m2 / 1e6,
                     "share": n / v.size})
    out["classes"] = rows
    return out


def extent_agreement(pred, ref, valid):
    """Contingency scores of a predicted against a reference extent.

    CSI (critical success index) = TP / (TP + FP + FN); hit rate = TP /
    (TP + FN); false-alarm ratio = FP / (TP + FP); bias = (TP + FP) / (TP +
    FN); F1 = 2TP / (2TP + FP + FN). Only cells valid in both count.
    """
    p = np.asarray(pred, bool) & valid
    r = np.asarray(ref, bool) & valid
    tp = int(np.sum(p & r))
    fp = int(np.sum(p & ~r))
    fn = int(np.sum(~p & r))
    tn = int(np.sum(valid & ~p & ~r))

    def div(a, b):
        return float(a) / b if b else float("nan")
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "csi": div(tp, tp + fp + fn), "hit_rate": div(tp, tp + fn),
            "false_alarm_ratio": div(fp, tp + fp), "bias": div(tp + fp, tp + fn),
            "f1": div(2 * tp, 2 * tp + fp + fn)}


def error_stats(est, ref):
    """Bias, MAE and RMSE of est - ref over cells finite in both."""
    a = np.asarray(est, dtype="float64")
    b = np.asarray(ref, dtype="float64")
    ok = np.isfinite(a) & np.isfinite(b)
    if not ok.any():
        return {"cells": 0}
    e = a[ok] - b[ok]
    return {"cells": int(ok.sum()), "bias_m": float(e.mean()), "mae_m": float(np.abs(e).mean()),
            "rmse_m": float(np.sqrt(np.mean(e ** 2))),
            "corr": float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 2 else float("nan")}
