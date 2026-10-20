# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Water indices, automatic thresholds, cleanup and connectivity.

    python -m flooddepth.tests.test_classify
"""

import numpy as np

from .common import check, finish
from ..core import classify as C


def test_indices():
    print("\n1. Indices")
    nd = C.normalized_difference(np.array([3.0, 1.0, 0.0, np.nan]), np.array([1.0, 3.0, 0.0, 1.0]))
    check("ND = (a-b)/(a+b)", abs(nd[0] - 0.5) < 1e-6 and abs(nd[1] + 0.5) < 1e-6)
    check("zero denominator -> NaN", np.isnan(nd[2]))
    check("NaN input -> NaN", np.isnan(nd[3]))
    check("dB of 0.1 = -10", abs(C.to_db(np.array([0.1]))[0] + 10.0) < 1e-5)
    check("dB of 0 -> NaN", np.isnan(C.to_db(np.array([0.0]))[0]))


def test_otsu():
    print("\n2. Otsu")
    t, eta = C.otsu(np.r_[np.full(5000, -0.5), np.full(1000, 0.5)], value_range=(-1, 1))
    check("empty gap: threshold in the middle of the gap", abs(t) < 0.01, f"t={t:.4f}")
    check("separable classes: eta ~ 1", eta > 0.99, f"eta={eta:.3f}")
    rng = np.random.default_rng(0)
    t, eta = C.otsu(rng.normal(0.0, 1.0, 50000))
    check("unimodal Gaussian: eta = 2/pi (0.637)", abs(eta - 2 / np.pi) < 0.01, f"eta={eta:.3f}")
    check("unimodal Gaussian labelled low", C.confidence_label(eta) == "low")
    try:
        C.otsu(np.full(1000, 0.3))
        check("constant values refused", False)
    except C.ThresholdError:
        check("constant values refused", True)
    try:
        C.otsu(np.array([0.1, 0.2]))
        check("too few values refused", False)
    except C.ThresholdError:
        check("too few values refused", True)
    lo, hi = C.otsu_band(np.r_[rng.normal(-0.4, 0.08, 20000), rng.normal(0.3, 0.08, 4000)])
    check("near-optimal band brackets the valley", lo < 0.0 < hi or (-0.15 < lo < hi < 0.15), f"{lo:.3f}..{hi:.3f}")


def test_tile_threshold():
    print("\n3. Tile-based Otsu on an imbalanced scene (2 % water)")
    rng = np.random.default_rng(1)
    img = rng.normal(-0.45, 0.12, (600, 600)).astype("float32")      # land, broad
    water = np.zeros_like(img, dtype=bool)
    water[280:292, :] = True                                         # a 12-cell river, 2 %
    img[water] = rng.normal(0.35, 0.05, water.sum())
    t = C.tile_threshold(img, True)
    w, _ = C.classify(img, t["threshold"], True)
    tw, _ = C.classify(img, t["global_threshold"], True)
    err_tile = (w != water).mean()
    err_glob = (tw != water).mean()
    check("tile method uses mixed tiles", t["tiles_used"] >= 5, f"{t['tiles_used']} tiles")
    check("tile threshold errors <= global", err_tile <= err_glob + 1e-9, f"tile {err_tile:.4f} vs global {err_glob:.4f}")
    check("tile threshold classifies the river (error < 0.5 %)", err_tile < 0.005, f"{100 * err_tile:.3f} %")
    check("band reported", "band" in t and t["band"][0] <= t["threshold"] <= t["band"][1] + 1e-6)


def test_components_and_cleanup():
    print("\n4. Connected components and cleanup")
    m = np.zeros((40, 40), bool)
    m[2:5, 2:5] = True          # 9 cells
    m[10, 10] = True            # 1 cell
    m[11, 11] = True            # diagonal neighbour -> same 8-connected region
    m[20:35, 20:35] = True      # 225 cells with a hole
    m[26:28, 26:28] = False     # 4-cell hole (island)
    lab, n = C.label_components(m)
    check("8-connected regions counted", n == 3, f"n={n}")
    lab4, n4 = C.label_components(m, connectivity=4)
    check("4-connected splits the diagonal pair", n4 == 4, f"n={n4}")
    valid = np.ones_like(m)
    out, st = C.clean_mask(m, valid, min_water_cells=5)
    check("small patches removed", not out[10, 10] and out[3, 3], str(st))
    check("island kept by default", not out[26, 26])
    out2, st2 = C.clean_mask(m, valid, 0, max_hole_cells=4)
    check("hole filled only when asked", out2[26, 26] and st2["holes_filled"] == 1, str(st2))
    edge = np.zeros((20, 20), bool)
    edge[:, 5:15] = True
    edge[0:2, 8:10] = False      # notch open to the grid edge
    out3, st3 = C.clean_mask(edge, np.ones_like(edge), 0, max_hole_cells=10)
    check("a notch touching the grid edge is not a hole", not out3[0, 8], str(st3))
    seed = np.zeros_like(m)
    seed[3, 3] = True
    keep, k, d = C.keep_connected(m, seed)
    check("keep_connected keeps only the seeded region", keep[2, 2] and not keep[25, 25] and k == 1 and d == 2)


def test_supervised():
    print("\n5. Nearest-centroid model")
    w = np.tile([[0.1, 0.2, 0.8]], (20, 1))
    l_ = np.tile([[0.8, 0.7, 0.2]], (20, 1))
    model = C.fit_centroids(w, l_)
    wat, val = C.classify_centroids(np.array([[[0.1, 0.8]], [[0.2, 0.7]], [[0.8, 0.2]]]), model)
    check("separates the two classes", wat.tolist() == [[True, False]])
    try:
        C.fit_centroids(np.ones((2, 3)), np.ones((20, 3)))
        check("too few training pixels refused", False)
    except ValueError:
        check("too few training pixels refused", True)


if __name__ == "__main__":
    test_indices()
    test_otsu()
    test_tile_threshold()
    test_components_and_cleanup()
    test_supervised()
    finish("test_classify")
