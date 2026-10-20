# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Water surface and depth on analytic terrain with known answers.

    python -m flooddepth.tests.test_surface
"""

import numpy as np

from .common import check, finish
from ..core import surface as SF
from ..core.raster import Grid, srs_from_epsg

WKT = srs_from_epsg(32736).ExportToWkt()


def grid(rows, cols, cs=10.0):
    return Grid(cols, rows, (500000.0, cs, 0.0, 9990000.0, 0.0, -cs), WKT)


def test_nearest():
    print("\n1. Nearest-seed transform vs brute force")
    rng = np.random.default_rng(4)
    ok = True
    for (r, c, n, cy) in ((70, 90, 25, 1.0), (120, 60, 7, 2.0), (50, 200, 1, 1.0), (33, 33, 200, 1.0)):
        s = np.zeros((r, c), bool)
        s.flat[rng.choice(r * c, n, replace=False)] = True
        sr, sc, d = SF.nearest_seed(s, 1.0, cy)
        rr, cc = np.nonzero(s)
        g = np.indices((r, c)).reshape(2, -1).T
        bd = np.sqrt(np.min((cy * (g[:, None, 0] - rr[None])) ** 2 + (g[:, None, 1] - cc[None]) ** 2, axis=1)).reshape(r, c)
        ok &= np.allclose(d, bd, atol=1e-3)
    check("exact distances (incl. anisotropic cells, one seed)", bool(ok))
    s = np.zeros((40, 40), bool)
    s[0, 0] = True
    flood = np.zeros_like(s)
    flood[:5, :5] = True
    v, d = SF.fill_nearest(np.arange(1600, dtype="float32").reshape(40, 40), s, 1, 1, targets=flood)
    check("bounded schedule exact on the target cells", np.all(v[flood] == 0))


def test_waterline():
    print("\n2. Waterline cells")
    valid = np.ones((10, 10), bool)
    flood = np.zeros_like(valid)
    flood[2:8, 0:6] = True               # touches the west grid edge
    valid[:, 9] = False
    line = SF.waterline(flood, valid)
    check("edge against dry ground is waterline", line[2, 3] and line[4, 5])
    check("grid edge is not waterline", not line[4, 0])
    other = np.zeros_like(valid)
    other[2:8, 6] = True                 # permanent water east of the flood
    line2 = SF.waterline(flood, valid, other_water=other)
    check("edge against permanent water is not waterline", not line2[4, 5])
    blocked = np.zeros_like(valid)
    blocked[1, :] = True                 # cloud north of the flood
    line3 = SF.waterline(flood, valid, blocked=blocked)
    check("edge against cloud is not waterline", not line3[2, 3])


def plane_case(slope_long=0.0005, cross=0.004, h=2.0, rows=200, cols=160, cs=10.0):
    y = np.arange(rows)[:, None] * cs
    x = np.arange(cols)[None, :] * cs
    xc = cols * cs / 2
    z = (1000.0 - slope_long * y + cross * np.abs(x - xc)).astype("float32")
    wse = (1000.0 - slope_long * y + h) * np.ones_like(z)
    flood = z < wse
    return z, wse.astype("float32"), flood


def test_plane():
    print("\n3. Sloping water surface in a V valley (known answer)")
    z, wse, flood = plane_case()
    g = grid(*z.shape)
    valid = np.ones_like(flood)
    for interp in ("nearest", "idw"):
        s, info, used, dropped = SF.estimate_surface(z, flood, valid, g, interpolation=interp)
        d, neg = SF.depth_from_surface(s, z, flood)
        true = (wse - z)[flood]
        err = d[flood] - true
        check(f"{interp}: depth MAE < 0.05 m", np.abs(err).mean() < 0.05, f"MAE {np.abs(err).mean():.4f} m, "
              f"max depth {np.nanmax(d):.2f} (true {true.max():.2f})")
    check("no waterline dropped on clean terrain", info["waterline_dropped"] == 0)
    z2 = z.copy()
    z2[:, :3] = np.nan
    valid2 = np.isfinite(z2)
    flood2 = flood.copy()
    flood2[:5, :] = False
    s, info, _, _ = SF.estimate_surface(z2, flood2, valid2, g)
    check("flood running off the grid edge: edge cells excluded", info["edge_cells_excluded"] > 0)


def test_outliers_and_steep():
    print("\n4. Outlier screening and steep ground")
    z, wse, flood = plane_case()
    g = grid(*z.shape)
    valid = np.ones_like(flood)
    line = SF.waterline(flood, valid)
    z_bad = z.copy()
    rr, cc = np.nonzero(line)
    pick = (rr[::25], cc[::25])
    z_bad[pick] += 6.0                     # embankment / canopy spikes on the edge
    s, info, used, dropped = SF.estimate_surface(z_bad, flood, valid, g)
    check("spikes on the edge dropped", dropped[pick].mean() > 0.9, f"{info['waterline_dropped']} dropped")
    d, _ = SF.depth_from_surface(s, z_bad, flood)
    err = (d - (wse - z_bad))[flood & ~line]
    check("depth unaffected by the spikes (MAE < 0.1 m)", np.nanmean(np.abs(err)) < 0.1, f"{np.nanmean(np.abs(err)):.3f}")
    # steep: an escarpment where the edge drops 3 m per cell
    zs, wses, fl = plane_case(slope_long=0.3, cross=0.05, h=1.0)
    s, info, used, dropped = SF.estimate_surface(zs, fl, np.ones_like(fl), grid(*zs.shape))
    check("steep edges never dropped as outliers", info["waterline_dropped"] == 0 and info["waterline_steep"] > 0)
    d, _ = SF.depth_from_surface(s, zs, fl)
    check("steep corridor: no invented depth (max <= true max + 0.5 m)", np.nanmax(d) <= 1.5, f"max {np.nanmax(d):.2f}")
    check("smoothing skipped where it would change the surface > tolerance", info.get("smoothing_skipped_cells", 0) > 0)


def test_edge_step():
    print("\n5. Edge water level between the last wet and first dry cell")
    z = np.full((60, 60), 2.0, dtype="float32")    # banks 2 m above a flat floor
    z[15:45, 15:45] = 0.0                            # floor; water 1 m deep fills it
    flood = z < 1.0
    g = grid(60, 60)
    v = np.ones_like(flood)
    s0, i0, _, _ = SF.estimate_surface(z, flood, v, g, edge_step=0.0)
    s1, i1, _, _ = SF.estimate_surface(z, flood, v, g, edge_step=0.5, edge_step_cap=1.0)
    s2, _, _, _ = SF.estimate_surface(z, flood, v, g, edge_step=0.5, edge_step_cap=0.4)
    d0, _ = SF.depth_from_surface(s0, z, flood)
    d1, _ = SF.depth_from_surface(s1, z, flood)
    d2, _ = SF.depth_from_surface(s2, z, flood)
    check("edge at the wet cell: depth 0 (underestimate)", abs(np.nanmean(d0)) < 1e-6)
    check("halfway up the 2 m step: the true 1 m", abs(np.nanmean(d1) - 1.0) < 1e-5, f"{np.nanmean(d1):.3f} m")
    check("the cap limits the lift", abs(np.nanmean(d2) - 0.4) < 1e-5, f"{np.nanmean(d2):.3f} m")
    check("lift reported", abs(i1["edge_lift_mean_m"] - 1.0) < 1e-6)


def test_depth_helpers():
    print("\n6. Depth helpers")
    d, neg = SF.depth_from_surface(np.array([[2.0, 1.0]]), np.array([[1.0, 1.5]]), np.array([[True, True]]))
    check("negative depth clipped to 0 and counted", d[0, 1] == 0 and neg == 1)
    c = SF.depth_constant_stage(10.0, np.array([[8.0, 12.0, np.nan]]))
    check("constant level: wet, dry, NoData", c[0, 0] == 2.0 and np.isnan(c[0, 1]) and np.isnan(c[0, 2]))
    b = SF.box_sum(np.ones((5, 5)), 1)
    check("box sum (corner 4, centre 9)", b[0, 0] == 4 and b[2, 2] == 9)


if __name__ == "__main__":
    test_nearest()
    test_waterline()
    test_plane()
    test_outliers_and_steep()
    test_edge_step()
    test_depth_helpers()
    finish("test_surface")
