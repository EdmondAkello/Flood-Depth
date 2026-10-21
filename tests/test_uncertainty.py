# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Correlated DEM error, evidence classes and the depth ensemble.

    python -m flooddepth.tests.test_uncertainty
"""

import numpy as np

from .common import check, finish
from ..core import uncertainty as U
from .test_surface import plane_case, grid


def main():
    print("\n1. Correlated noise")
    rng = np.random.default_rng(0)
    n = U.correlated_noise((400, 400), 2.0, 10.0, rng)
    check("standard deviation = sigma", abs(n.std() - 2.0) < 0.05, f"{n.std():.3f}")
    c1 = np.corrcoef(n[:, :-5].ravel(), n[:, 5:].ravel())[0, 1]
    c2 = np.corrcoef(n[:, :-60].ravel(), n[:, 60:].ravel())[0, 1]
    check("correlated at short lag, not at long lag", c1 > 0.8 and abs(c2) < 0.2, f"r(5)={c1:.2f}, r(60)={c2:.2f}")
    w = U.correlated_noise((100, 100), 1.0, 0.0, rng)
    check("zero correlation length = white noise", abs(np.corrcoef(w[:, :-1].ravel(), w[:, 1:].ravel())[0, 1]) < 0.05)
    check("DEM preset from name: FABDEM 1.5 m", U.dem_sigma_from_name("nyando_fabdem.tif")[0] == 1.5)
    check("DEM preset unknown -> default with a note", U.dem_sigma_from_name("dem.tif")[1].startswith("unknown"))
    b = U.block_mean(np.arange(16, dtype="float32").reshape(4, 4), 2)
    check("block mean", b.tolist() == [[2.5, 4.5], [10.5, 12.5]])
    check("expand back", U.expand(b, 2, (4, 3)).shape == (4, 3))

    print("\n2. Evidence classes")
    flood = np.ones((3, 3), bool)
    flood[0, 0] = False
    bare = np.zeros_like(flood)
    bare[1, 1] = True
    slope = np.zeros((3, 3))
    slope[2, 2] = 30.0
    ev = U.evidence_classes(flood, "trace", bare, slope, 10.0)
    check("outside 255, inundation 2, erosion 3, steep 4",
          ev[0, 0] == 255 and ev[0, 1] == 2 and ev[1, 1] == 3 and ev[2, 2] == 4)
    check("image extent -> observed water", U.evidence_classes(flood, "image")[1, 1] == 1)
    above = np.zeros_like(flood)
    above[0, 1] = above[1, 1] = above[2, 2] = True
    ev5 = U.evidence_classes(flood, "trace", bare, slope, 10.0, above=above)
    check("ground above the surface -> doubtful 5, but erosion 3 and steep 4 are kept",
          ev5[0, 1] == 5 and ev5[1, 1] == 3 and ev5[2, 2] == 4 and ev5[0, 0] == 255)

    print("\n3. Method error (blocked cross-validation of the waterline)")
    z, wse, fl = plane_case()
    g = grid(*z.shape)
    sd0, b0, n0 = U.method_error(z, fl, np.ones_like(fl), g)
    check("planar water surface: method error ~ 0", sd0 < 0.05 and n0 > 0, f"sd {sd0:.3f} m, block {b0} cells")
    rng2 = np.random.default_rng(3)
    zr = z + U.correlated_noise(z.shape, 0.5, 8.0, rng2)
    sd1, _, _ = U.method_error(zr, fl, np.ones_like(fl), g)
    check("rough edge elevations: method error grows", sd1 > 5 * max(sd0, 0.01), f"sd {sd1:.3f} m")

    print("\n4. Ensemble on a known sloping surface")
    z, wse, fl = plane_case()
    g = grid(*z.shape)
    full, s = U.run_ensemble(z, np.ones_like(fl), g, fl, members=20, sigma_dem=0.3, corr_m=100.0,
                             edge_shift=1, seed=7)
    true = np.where(fl, wse - z, np.nan)
    inside = fl & np.isfinite(full["p10"]) & np.isfinite(full["p90"])
    cover = np.mean((full["p10"][inside] - 0.15 <= true[inside]) & (true[inside] <= full["p90"][inside] + 0.15))
    check("P10-P90 (+-0.15 m) covers the true depth in >= 80 % of cells", cover >= 0.8, f"{100 * cover:.0f} %")
    check("P10 <= P50 <= P90", bool(np.all(full["p10"][inside] <= full["p50"][inside] + 1e-6)
                                    and np.all(full["p50"][inside] <= full["p90"][inside] + 1e-6)))
    check("probability 1 in the middle, < 1 at the edge", full["prob"][100, 80] == 1.0
          and full["prob"][fl & ~np.roll(fl, 1, 1)].mean() < 1.0)
    a = s["area_km2_p10_p50_p90"]
    check("area band contains the base area", a[0] <= fl.sum() * 100 / 1e6 <= a[2], str([round(x, 4) for x in a]))
    full2, s2 = U.run_ensemble(z, np.ones_like(fl), g, fl, members=6, sigma_dem=0.3, max_cells=8000, seed=7)
    check("coarsened ensemble runs and reports its factor", s2["coarsening"] >= 2 and full2["p90"].shape == z.shape,
          f"factor {s2['coarsening']}")
    full3, _ = U.run_ensemble(z, np.ones_like(fl), g, fl, members=6, sigma_dem=0.3, seed=7)
    check("same seed -> same result", np.allclose(np.nan_to_num(full3["p50"]),
                                                  np.nan_to_num(U.run_ensemble(z, np.ones_like(fl), g, fl, members=6,
                                                                               sigma_dem=0.3, seed=7)[0]["p50"])))


if __name__ == "__main__":
    main()
    finish("test_uncertainty")
