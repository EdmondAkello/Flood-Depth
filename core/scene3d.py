# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Interactive 3D scene: one self-contained HTML file (WebGL 2, offline).

The terrain (DEM, draped with the satellite image or a hillshade, permanent
water and doubtful cells painted on it), the water surface coloured by depth,
the flood edge and any reference lines, opened in the 3D viewer of the
LandXML Road & Terrain GIS Tools (same author; `viewer3d.html`): rotate, pan,
zoom, layer switches and opacity, vertical exaggeration, click to read the
level and depth, saved views, PNG and A3 export.

Meshes are regular grids on the analysis grid, block-averaged so the scene
stays within `max_vertices` per surface; they are packed as the viewer
expects (byte-planar, delta-coded indices, gzip, base64). No QGIS imports.
"""

import base64
import gzip
import html
import json
import math
import os
import time

import numpy as np

from . import raster as R
from . import report as RP

TEMPLATE = os.path.join(os.path.dirname(__file__), "viewer3d.html")
MAX_VERTICES = 600_000
_TYPES = {"<f4": "f4", "<u4": "u4", "u1": "u1"}


def _pack(array, dtype, delta=False):
    """An array packed for the viewer (same encoding as the LandXML scene writer)."""
    data = np.ascontiguousarray(array, dtype=dtype)
    components = data.shape[1] if data.ndim == 2 else 1
    flat = data.T.ravel() if components > 1 else data.ravel()
    if delta:
        steps = np.diff(np.concatenate([[0], flat.astype(np.int64)]))
        flat = ((steps << 1) ^ (steps >> 63)).astype("<u4")
    size = flat.dtype.itemsize
    planes = flat.view(np.uint8).reshape(-1, size).T.tobytes() if size > 1 else flat.tobytes()
    return {"z": base64.b64encode(gzip.compress(planes, 6, mtime=0)).decode("ascii"), "t": _TYPES[dtype],
            "c": components, "d": int(bool(delta))}


def unpack(packed):
    """The array back from `_pack` (tests)."""
    dtype = {"f4": "<f4", "u4": "<u4", "u1": "u1"}[packed["t"]]
    raw = np.frombuffer(gzip.decompress(base64.b64decode(packed["z"])), np.uint8)
    size = np.dtype(dtype).itemsize
    flat = raw.reshape(size, -1).T.copy().view(dtype).ravel() if size > 1 else raw.copy()
    if packed["d"]:
        steps = flat.astype(np.int64)
        flat = np.cumsum((steps >> 1) ^ -(steps & 1)).astype(dtype)
    c = packed["c"]
    return flat.reshape(c, -1).T.copy() if c > 1 else flat


def _block_nanmean(a, step):
    rows, cols = a.shape
    br, bc = -(-rows // step), -(-cols // step)
    pad = np.full((br * step, bc * step), np.nan, dtype="float64")
    pad[:rows, :cols] = a
    blk = pad.reshape(br, step, bc, step)
    with np.errstate(invalid="ignore"):
        n = np.isfinite(blk).sum(axis=(1, 3))
        s = np.nansum(blk, axis=(1, 3))
        return np.where(n > 0, s / np.maximum(n, 1), np.nan)


def _block_frac(mask, step):
    return _block_nanmean(np.where(np.isfinite(mask.astype("float64")), mask.astype("float64"), np.nan), step)


def grid_mesh(z, keep, x0, y0, cell):
    """Vertices (x, y, z) at the cells of `z` where `keep`, and two triangles
    per square whose four corners are kept. Returns (xyz, faces, index map)."""
    rows, cols = z.shape
    ok = keep & np.isfinite(z)
    idx = np.full(z.shape, -1, dtype=np.int64)
    rr, cc = np.nonzero(ok)
    idx[rr, cc] = np.arange(rr.size)
    xyz = np.column_stack([x0 + (cc + 0.5) * cell, y0 - (rr + 0.5) * cell, z[rr, cc]])
    a, b = idx[:-1, :-1], idx[:-1, 1:]
    c, d = idx[1:, :-1], idx[1:, 1:]
    quad = (a >= 0) & (b >= 0) & (c >= 0) & (d >= 0)
    f1 = np.column_stack([a[quad], c[quad], b[quad]])
    f2 = np.column_stack([b[quad], c[quad], d[quad]])
    return xyz, np.vstack([f1, f2]), idx


def _classes(stops, unit="m"):
    out = []
    for i, (v, c) in enumerate(stops):
        hi = stops[i + 1][0] if i + 1 < len(stops) else None
        lab = f"{v:g}–{hi:g} {unit}" if hi is not None else f"> {v:g} {unit}"
        out.append({"low": float(v), "high": None if hi is None else float(hi),
                    "color": "#{:02x}{:02x}{:02x}".format(*[int(x) for x in c]), "label": lab})
    return out


def _ring_lines(wkts, crs_wkt, grid_wkt, simplify):
    """Polylines (list of Nx2 arrays) of the outlines of polygons / the lines given as WKT."""
    from osgeo import ogr, osr
    out = []
    ct = None
    if crs_wkt and grid_wkt:
        src, dst = R.srs_from_wkt(crs_wkt), R.srs_from_wkt(grid_wkt)
        if not src.IsSame(dst):
            ct = osr.CoordinateTransformation(src, dst)
    for w in wkts or []:
        g = ogr.CreateGeometryFromWkt(w)
        if g is None:
            continue
        if ct is not None:
            g.Transform(ct)
        if simplify:
            g = g.SimplifyPreserveTopology(simplify)
        stack = [g]
        while stack:
            h = stack.pop()
            n = h.GetGeometryCount()
            if n:
                stack.extend(h.GetGeometryRef(i) for i in range(n))
            elif h.GetPointCount() > 1:
                out.append(np.array([h.GetPoint_2D(i) for i in range(h.GetPointCount())], float))
    return out


def write_scene(path, grid, dem, flood, wse, depth, evidence=None, other=None, backdrop=None, style=None,
                extent_wkts=None, extent_crs=None, lines=None, title="Flood depth - 3D scene",
                max_vertices=MAX_VERTICES, vertical_exaggeration=None):
    """Write the 3D scene HTML. `lines`: [(name, kind, wkts, crs_wkt, '#rrggbb')] drawn on the
    terrain (connectivity features, area of interest). Returns a dict with what was written."""
    st = RP.MapStyle.from_dict(style)
    flood = np.asarray(flood, bool)
    if not flood.any():
        raise ValueError("No flood cells: nothing to show in 3D.")
    rows, cols = flood.shape
    rr, cc = np.nonzero(flood)
    pad = max(30, int(0.3 * max(np.ptp(rr), np.ptp(cc))))
    r0, r1 = max(0, rr.min() - pad), min(rows, rr.max() + pad + 1)
    c0, c1 = max(0, cc.min() - pad), min(cols, cc.max() + pad + 1)
    win = (slice(r0, r1), slice(c0, c1))
    step = max(1, int(math.ceil(math.sqrt((r1 - r0) * (c1 - c0) / float(max_vertices)))))
    cell = grid.cell_x * step
    x0 = grid.gt[0] + c0 * grid.gt[1]
    y0 = grid.gt[3] + r0 * grid.gt[5]
    z = _block_nanmean(np.asarray(dem, "float64")[win], step)
    fl = _block_frac(flood[win], step) >= 0.5
    ws = _block_nanmean(np.where(flood, wse, np.nan)[win], step)
    dp = _block_nanmean(np.where(flood, depth, np.nan)[win], step)

    # terrain colours: the image (or hillshade + tint), permanent water and doubtful cells painted on
    sub = R.Grid(z.shape[1], z.shape[0], (x0, cell, 0.0, y0, 0.0, -cell), grid.wkt)
    rgb, drape = None, "hillshade and elevation tint"
    if backdrop is not None:
        try:
            bd = backdrop(sub)
            if bd is not None:
                rgb, drape = bd[0], bd[1]
        except Exception:  # noqa: BLE001 - the scene must not fail on the drape
            rgb = None
    if rgb is not None:
        ok = np.all(np.isfinite(rgb), axis=-1)
        base = np.where(ok[..., None], 255.0 * np.nan_to_num(rgb), RP._base_image(z, cell, None))
    else:
        base = RP._base_image(z, cell, None)
    if other is not None:
        ow = _block_frac(np.asarray(other, bool)[win], step) >= 0.5
        base[ow] = 0.15 * base[ow] + 0.85 * np.array(st.permanent, float)
    mud = np.zeros(z.shape, bool)
    if evidence is not None:
        ev = np.asarray(evidence)[win]
        dbt = _block_frac(ev == 5, step) >= 0.5
        base[dbt] = 0.25 * base[dbt] + 0.75 * np.array(st.doubtful, float)
        # erosion / deposition: the mud and debris a flash flood leaves after the water drained
        mud = (_block_frac(ev == 3, step) >= 0.5) & fl
        base[mud] = 0.2 * base[mud] + 0.8 * np.array(st.mud, float)
    base = np.clip(base, 0, 255).astype("uint8")

    t_xyz, t_faces, t_idx = grid_mesh(z, np.isfinite(z), x0, y0, cell)
    t_rgb = base[np.isfinite(z)]
    w_keep = fl & np.isfinite(ws) & ~mud
    w_xyz, w_faces, w_idx = grid_mesh(ws, w_keep, x0, y0, cell)
    w_depth = np.nan_to_num(dp[w_keep], nan=0.0)
    stops = st.depth_stops()
    w_rgb = np.clip(RP.ramp(w_depth, stops), 0, 255).astype("uint8")

    m_keep = mud & np.isfinite(z)
    m_xyz, m_faces, _ = grid_mesh(np.where(m_keep, z + 0.3, np.nan), m_keep, x0, y0, cell)
    m_depth = np.nan_to_num(dp[m_keep], nan=0.0)
    shade = np.clip(0.75 + 0.25 * (m_xyz[:, 2] - np.nanmin(z)) / max(float(np.nanmax(z) - np.nanmin(z)), 1.0), 0.6, 1.0) \
        if len(m_xyz) else np.zeros(0)
    m_rgb = np.clip(np.array(st.mud, float)[None, :] * shade[:, None], 0, 255).astype("uint8")

    allz = np.concatenate([t_xyz[:, 2], w_xyz[:, 2]]) if len(w_xyz) else t_xyz[:, 2]
    origin = np.array([float(t_xyz[:, 0].mean()), float(t_xyz[:, 1].mean()), float(np.nanmin(allz))])
    lo, hi = np.nanpercentile(t_xyz[:, 2], [5, 95])
    size = math.hypot((c1 - c0) * grid.cell_x, (r1 - r0) * grid.cell_y)
    ve = vertical_exaggeration or float(np.clip(round(0.06 * size / max(hi - lo, 1.0), 1), 1.0, 8.0))

    def rel(p):
        return (np.asarray(p, float) - origin).astype("float32")

    srs = R.srs_from_wkt(grid.wkt)
    try:
        srs.AutoIdentifyEPSG()
    except RuntimeError:
        pass
    code = srs.GetAuthorityCode(None)
    surfaces = [
        {"id": "terrain", "name": "Terrain", "role": "existing", "mesh": "terrain", "style": None,
         "opacity": 1.0, "lift": 0.0},
        {"id": "water", "name": "Flood water (by depth)", "role": "water", "mesh": "water",
         "style": {"attribute": "depth", "unit": "m", "classes": _classes(stops)}, "opacity": 0.9,
         "lift": 0.03, "depth_label": "Depth"},
    ]
    meshes = {"terrain": {"positions": _pack(rel(t_xyz), "<f4"), "colors": _pack(t_rgb, "u1"),
                          "indices": _pack(t_faces.ravel(), "<u4", delta=True), "depth": None}}
    if len(w_faces):
        meshes["water"] = {"positions": _pack(rel(w_xyz), "<f4"), "colors": _pack(w_rgb, "u1"),
                           "indices": _pack(w_faces.ravel(), "<u4", delta=True),
                           "depth": _pack(w_depth.astype("float32"), "<f4")}
    else:
        surfaces.pop()
    if len(m_faces):
        hx = "#{:02x}{:02x}{:02x}".format
        surfaces.append({"id": "mud", "name": "Mud and debris (erosion / deposition)", "role": "mud", "mesh": "mud",
                         "style": {"attribute": "evidence", "unit": "",
                                   "classes": [{"low": None, "high": None, "color": hx(*st.mud),
                                                "label": "Mud, debris, scour (flood trace)"}]},
                         "opacity": 1.0, "lift": 0.0, "depth_label": "Flood depth when it passed"})
        meshes["mud"] = {"positions": _pack(rel(m_xyz), "<f4"), "colors": _pack(m_rgb, "u1"),
                         "indices": _pack(m_faces.ravel(), "<u4", delta=True),
                         "depth": _pack(m_depth.astype("float32"), "<f4")}

    def on_terrain(xy, lift=0.0):
        rr_ = np.clip(((y0 - xy[:, 1]) / cell).astype(int), 0, z.shape[0] - 1)
        cc_ = np.clip(((xy[:, 0] - x0) / cell).astype(int), 0, z.shape[1] - 1)
        zz = np.where(np.isfinite(ws[rr_, cc_]), ws[rr_, cc_], z[rr_, cc_])
        ok = np.isfinite(zz)
        return np.column_stack([xy[ok], zz[ok] + lift])

    out_lines = []
    budget = 80_000

    def add_lines(name, kind, wkts, crs, color, width=2.5):
        nonlocal budget
        for ln in _ring_lines(wkts, crs, grid.wkt, simplify=cell * 0.75):
            if budget <= 0:
                return
            if len(ln) < 2:
                continue
            p = on_terrain(ln, 0.05)
            if len(p) < 2:
                continue
            budget -= len(p)
            out_lines.append({"id": f"{kind}-{len(out_lines) + 1}", "name": name, "kind": kind,
                              "points": [[round(float(v), 3) for v in q] for q in rel(p)],
                              "color": color, "width": width, "lift": 0.1, "properties": {}})

    hx = "#{:02x}{:02x}{:02x}".format
    if extent_wkts:
        big = [w for w in extent_wkts]
        add_lines("Flood edge", "Flood edge", big, extent_crs, hx(*st.edge), 2.0)
    for name, kind, wkts, crs, color in lines or []:
        add_lines(name, kind, wkts, crs, color, 3.0)

    labels = []
    if len(m_xyz):
        k = int(np.argmax(m_depth))
        labels.append({"text": f"Mud and debris: {mud.sum() * cell * cell / 1e6:.2f} km²",
                       "position": [round(float(v), 3) for v in rel(m_xyz[k])], "lift": 8.0, "kind": "label"})
    if len(w_xyz):
        k = int(np.argmax(w_depth))
        labels.append({"text": f"Deepest: {w_depth[k]:.1f} m", "position": [round(float(v), 3) for v in rel(w_xyz[k])],
                       "lift": 5.0, "kind": "label"})
    ext_min = rel(np.nanmin(np.vstack([t_xyz, w_xyz]) if len(w_xyz) else t_xyz, axis=0))
    ext_max = rel(np.nanmax(np.vstack([t_xyz, w_xyz]) if len(w_xyz) else t_xyz, axis=0))
    centre = ((ext_min + ext_max) / 2.0).tolist()
    views = [{"name": n, "direction": d, "distance": 0.95 * size, "target": [round(v, 3) for v in centre],
              "view_angle": 30.0} for n, d in (("South-west", [-0.6, -0.8, 0.55]), ("South-east", [0.8, -0.5, 0.6]),
                                               ("From above", [0.0, -0.05, 1.0]))]
    data = {"format": "landxmltools-scene", "version": 2, "generator": "FloodDepth",
            "created": time.strftime("%Y-%m-%dT%H:%M:%S"), "title": title,
            "source": {"file": title}, "crs": {"authid": f"EPSG:{code}" if code else "", "wkt": ""},
            "units": {"horizontal": "meter", "vertical": "meter"}, "origin": [float(v) for v in origin],
            "extent": {"min": [float(v) for v in ext_min], "max": [float(v) for v in ext_max]},
            "vertical_exaggeration": ve, "background": "#ffffff",
            "ground_mode": "exported", "contours": False,
            "surfaces": surfaces, "meshes": meshes, "lines": out_lines, "labels": labels, "views": views}
    with open(TEMPLATE, encoding="utf-8") as f:
        template = f.read()
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    text = template.replace("__TITLE__", html.escape(title)).replace("__SCENE__", payload)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    return {"path": str(path), "step": step, "cell_m": cell, "terrain_triangles": int(len(t_faces)),
            "water_triangles": int(len(w_faces)), "mud_triangles": int(len(m_faces)), "drape": drape, "vertical_exaggeration": ve,
            "lines": len(out_lines), "size_mb": round(os.path.getsize(path) / 1e6, 2)}


def scene_from_folder(out_dir, path=None, max_vertices=MAX_VERTICES, style=None, extra_lines=None,
                      vertical_exaggeration=None):
    """Rebuild the arrays of a finished Quick flood depth run (run.json and its rasters) and write
    the scene. The DEM is read again from the path recorded in run.json."""
    from . import dem as DM
    from . import imagery as IMG
    with open(os.path.join(out_dir, "run.json"), encoding="utf-8") as f:
        info = json.load(f)
    o, outs = info.get("options", {}), info.get("outputs", {})
    depth, grid, _ = R.read_band(outs["depth"], 1)
    flood = np.isfinite(depth)
    src, _ = DM.dem_source([o["dem_path"]] + list(o.get("dem_tiles") or []), grid.bounds, grid.wkt)
    dem, _ = R.read_on_grid(src, grid, 1, "bilinear", zero_nodata_auto=True)
    evid = None
    if outs.get("evidence") and os.path.exists(outs["evidence"]):
        evid, _, _ = R.read_band(outs["evidence"], 1)
    wse = np.where(flood, dem + depth, np.nan)      # water level (also for HAND runs, whose wse.tif is a height)
    other = None
    if outs.get("water_mask") and os.path.exists(outs["water_mask"]):
        wm, _ = R.read_on_grid(outs["water_mask"], grid, 1, "near")
        other = wm == 2
    if o.get("permanent_raster"):
        from .pipeline import permanent_from_raster
        pr, _ = R.read_on_grid(o["permanent_raster"], grid, 1, "near")
        pw, _ = permanent_from_raster(pr, o.get("permanent_values") or "auto")
        other = pw if other is None else (other | pw)
    sc = None
    if o.get("image_paths"):
        try:
            sc = IMG.scene_from_paths(o["image_paths"])
        except ValueError:
            sc = None

    def bd(g):
        return IMG.backdrop(sc, g) if sc is not None else None
    ext = []
    gp = outs.get("extent")
    if gp and os.path.exists(gp):
        from osgeo import ogr
        ds = ogr.Open(gp)
        lay = ds.GetLayerByName("flood_extent")
        ext = [f.GetGeometryRef().ExportToWkt() for f in lay] if lay is not None else []
        ds = None
    return write_scene(path or os.path.join(out_dir, "scene3d.html"), grid, dem, flood, wse, depth, evidence=evid,
                       other=other, backdrop=bd, style=style or o.get("map_style"), extent_wkts=ext,
                       extent_crs=grid.wkt, lines=extra_lines, title=f"Flood depth - {os.path.basename(out_dir)}",
                       max_vertices=max_vertices, vertical_exaggeration=vertical_exaggeration)
