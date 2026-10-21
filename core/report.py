# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Self-contained HTML run report: maps (inline PNG inside SVG), SVG chart, tables.

Maps are framed on the flood (with an overview inset when the frame is a
small part of the analysis area) and drawn over the satellite image when one
was used, else over a relief-stretched hillshade with an elevation tint so
that flat floodplains still show their relief.
"""

import base64
import html
import math
import os
import re

import numpy as np

from . import raster as R

# sequential blue ramp for depth (light = shallow, dark = deep)
DEPTH_RAMP = [(0.0, (198, 219, 239)), (0.5, (107, 174, 214)), (1.0, (49, 130, 189)),
              (2.0, (8, 81, 156)), (3.0, (8, 48, 107))]
# flood probability (share of ensemble members flooded)
PROB_RAMP = [(0.0, (255, 247, 188)), (0.5, (254, 196, 79)), (0.8, (217, 95, 14)), (1.0, (153, 52, 4))]
# elevation tint for the relief backdrop (low -> high, all light)
TINT_RAMP = [(0.0, (222, 233, 214)), (0.5, (236, 228, 204)), (1.0, (219, 200, 172))]
PERMANENT_RGB = (120, 80, 160)       # same as the water_mask colour table
OUTLINE_RGB = (8, 37, 78)
NODATA_RGB = (205, 205, 205)
DOUBTFUL_RGB = (231, 138, 195)      # evidence class 5 (same as evidence.tif)

MAP_PX = 1000
INSET_PX = 210


# ---------------------------------------------------------------- raster helpers

def hillshade(dem, cell, az=315.0, alt=45.0, z=1.0):
    z0 = np.where(np.isfinite(dem), dem, np.nanmean(dem) if np.isfinite(dem).any() else 0.0)
    gy, gx = np.gradient(z0 * z, cell, cell)
    slope = np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, gy)
    a, b = math.radians(az), math.radians(alt)
    hs = np.sin(b) * np.cos(slope) + np.cos(b) * np.sin(slope) * np.cos(a - aspect)
    return np.clip(hs, 0, 1)


def relief(dem, cell):
    """Multidirectional hillshade with a vertical exaggeration chosen from the
    terrain itself (the 90th-percentile slope is shown as ~20 degrees), so
    flat floodplains are not a uniform grey."""
    fin = np.isfinite(dem)
    if fin.sum() < 16:
        return np.ones(dem.shape)
    z0 = np.where(fin, dem, np.nanmean(dem))
    gy, gx = np.gradient(z0, cell, cell)
    g = np.hypot(gx, gy)[fin]
    p90 = float(np.percentile(g, 90)) if g.size else 0.0
    zf = float(np.clip(math.tan(math.radians(20)) / p90, 1.0, 60.0)) if p90 > 0 else 1.0
    hs = sum(w * hillshade(dem, cell, az, 40.0, zf) for az, w in ((270, 0.2), (315, 0.4), (0, 0.2), (225, 0.2)))
    return hs


def ramp(values, stops):
    v = np.asarray(values, dtype="float64")
    xs = [s[0] for s in stops]
    out = np.zeros(v.shape + (3,), dtype="float64")
    for ch in range(3):
        out[..., ch] = np.interp(v, xs, [s[1][ch] for s in stops])
    return out


def _block(a, step, how):
    """Reduce a 2-D array by `step` x `step` blocks: 'mean' / 'max' (NaN-aware) or 'any' / 'frac' (bool)."""
    if step == 1:
        return a.astype("float64") if how in ("mean", "max", "frac") else a.astype(bool)
    r, c = a.shape
    R2, C2 = -(-r // step) * step, -(-c // step) * step
    if how in ("any", "frac"):
        p = np.zeros((R2, C2), dtype="float32")
        p[:r, :c] = a
        b = p.reshape(R2 // step, step, C2 // step, step)
        return b.max(axis=(1, 3)) > 0 if how == "any" else b.mean(axis=(1, 3))
    p = np.full((R2, C2), np.nan, dtype="float64")
    p[:r, :c] = a
    b = p.reshape(R2 // step, step, C2 // step, step)
    with np.errstate(invalid="ignore"), _nanwarn():
        return np.nanmean(b, axis=(1, 3)) if how == "mean" else np.nanmax(b, axis=(1, 3))


class _nanwarn(object):
    def __enter__(self):
        import warnings
        self._w = warnings.catch_warnings()
        self._w.__enter__()
        warnings.simplefilter("ignore", RuntimeWarning)

    def __exit__(self, *a):
        self._w.__exit__(*a)


def _png_b64(img, fmt="JPEG"):
    """uint8 HxWx3 -> base64 image via GDAL (no extra dependency). JPEG by
    default: a satellite backdrop is photographic, and PNG made each map
    ~1.7 MB; JPEG at quality 88 keeps the depth colours and edges legible."""
    gdal = R._gdal()
    h, w, _ = img.shape
    mem = gdal.GetDriverByName("MEM").Create("", w, h, 3, gdal.GDT_Byte)
    for i in range(3):
        mem.GetRasterBand(i + 1).WriteArray(img[..., i])
    drv = gdal.GetDriverByName(fmt) or gdal.GetDriverByName("PNG")
    ext = "jpg" if drv.ShortName == "JPEG" else "png"
    path = "/vsimem/flooddepth_map_%d.%s" % (id(img), ext)
    with R.quiet():
        drv.CreateCopy(path, mem, options=["QUALITY=88"] if ext == "jpg" else [])
    fp = gdal.VSIFOpenL(path, "rb")
    gdal.VSIFSeekL(fp, 0, 2)
    size = gdal.VSIFTellL(fp)
    gdal.VSIFSeekL(fp, 0, 0)
    data = gdal.VSIFReadL(1, size, fp)
    gdal.VSIFCloseL(fp)
    for f in (path, path + ".aux.xml"):
        try:
            if gdal.VSIStatL(f) is not None:
                gdal.Unlink(f)
        except RuntimeError:
            pass
    mime = "image/jpeg" if ext == "jpg" else "image/png"
    return mime + ";base64," + base64.b64encode(data).decode("ascii")


# ---------------------------------------------------------------- framing

def flood_frame(flood, share=0.8, pad_frac=0.12, min_cells=40, max_share=0.6):
    """(r0, r1, c0, c1, zoomed, held): the smallest window holding `share` of the
    flood cells (searched on a coarse grid), padded and widened to a sensible
    aspect. `held` is the share of flood cells inside the frame. The whole grid
    is returned (zoomed False) when the window would cover most of it."""
    rows, cols = flood.shape
    n = int(flood.sum())
    if n == 0:
        return 0, rows, 0, cols, False, 1.0
    f = max(1, int(math.ceil(max(rows, cols) / 64.0)))
    cnt = _block(flood, f, "frac") * (f * f)
    gr, gc = cnt.shape
    P = np.zeros((gr + 1, gc + 1))
    P[1:, 1:] = cnt.cumsum(0).cumsum(1)
    need = share * P[-1, -1]
    best = None
    for a in range(gr):
        for b in range(a + 1, gr + 1):
            col = P[b, :] - P[a, :]            # cumulative column sums of rows a..b-1
            if col[-1] < need:
                continue
            c1 = np.searchsorted(col, col[:-1] + need - 1e-9)   # first c1 with sum(c0..c1) >= need
            ok = c1 <= gc
            if not ok.any():
                continue
            c0s = np.flatnonzero(ok)
            area = (b - a) * (c1[ok] - c0s)
            k = int(np.argmin(area))
            if best is None or area[k] < best[0]:
                best = (area[k], a, b, c0s[k], c1[ok][k])
    _, a, b, c0, c1 = best
    r0, r1, c0, c1 = a * f, min(rows, b * f), c0 * f, min(cols, c1 * f)
    h, w = r1 - r0, c1 - c0
    ph = max(int(h * pad_frac), min_cells // 2)
    pw = max(int(w * pad_frac), min_cells // 2)
    r0, r1, c0, c1 = r0 - ph, r1 + ph, c0 - pw, c1 + pw
    # keep the aspect between 1:1.6 and 1.6:1 so a thin flood is shown in context
    h, w = r1 - r0, c1 - c0
    if w > 1.6 * h:
        g = int((w / 1.6 - h) / 2) + 1
        r0, r1 = r0 - g, r1 + g
    elif h > 1.6 * w:
        g = int((h / 1.6 - w) / 2) + 1
        c0, c1 = c0 - g, c1 + g
    r0, c0 = max(0, r0), max(0, c0)
    r1, c1 = min(rows, r1), min(cols, c1)
    if (r1 - r0) * (c1 - c0) > max_share * rows * cols:
        return 0, rows, 0, cols, False, 1.0
    held = float(flood[r0:r1, c0:c1].sum()) / n
    return r0, r1, c0, c1, True, held


def _sub_grid(grid, r0, c0, rows, cols, step):
    gt = grid.gt
    return R.Grid(cols=cols, rows=rows,
                  gt=(gt[0] + c0 * gt[1], gt[1] * step, 0.0, gt[3] + r0 * gt[5], 0.0, gt[5] * step),
                  wkt=grid.wkt)


def _base_image(dem_d, cell, backdrop_rgb):
    """Background RGB (0-255 float) for a display window."""
    if backdrop_rgb is not None:
        rgb = backdrop_rgb.astype("float64")
        ok = np.all(np.isfinite(rgb), axis=-1)
        gray = rgb.mean(axis=-1, keepdims=True)
        rgb = 0.7 * rgb + 0.3 * gray                      # muted so the depth colours stand out
        rgb = 255.0 * (0.12 + 0.82 * rgb)
        if dem_d is not None and np.isfinite(dem_d).any():
            fill = _base_image(dem_d, cell, None)
            rgb = np.where(ok[..., None], np.nan_to_num(rgb), fill)
        else:
            rgb = np.where(ok[..., None], np.nan_to_num(rgb), 245.0)
        return rgb
    fin = np.isfinite(dem_d)
    if not fin.any():
        return np.full(dem_d.shape + (3,), 245.0)
    lo, hi = np.percentile(dem_d[fin], [2, 98])
    t = np.clip((np.nan_to_num(dem_d, nan=lo) - lo) / max(hi - lo, 1e-6), 0, 1)
    tint = ramp(t, TINT_RAMP)
    hs = relief(dem_d, cell)
    img = tint * (0.5 + 0.5 * hs)[..., None]
    img[~fin] = 250.0
    return img


def _hatch(shape, period=7):
    r, c = np.indices(shape)
    return ((r + c) % period) < 2


def _compose(base, layers):
    img = base.copy()
    for mask, rgb, alpha in layers:
        if mask is None or not mask.any():
            continue
        col = rgb if np.ndim(rgb) == 3 else np.broadcast_to(np.asarray(rgb, dtype="float64"), img.shape)
        img[mask] = (1 - alpha) * img[mask] + alpha * col[mask]
    return np.clip(img, 0, 255).astype("uint8")


def _outline(mask):
    m = mask.astype(bool)
    inner = m.copy()
    inner[1:, :] &= m[:-1, :]
    inner[:-1, :] &= m[1:, :]
    inner[:, 1:] &= m[:, :-1]
    inner[:, :-1] &= m[:, 1:]
    return m & ~inner


def _nice_length(m):
    if m <= 0:
        return 0
    e = 10 ** math.floor(math.log10(m))
    for k in (5, 2, 1):
        if k * e <= m:
            return k * e
    return e


def _svg_map(b64, w, h, metres_per_px, inset=None, label=""):
    """Wrap the PNG in an SVG with a scale bar, a north arrow and an optional inset."""
    fs = max(w, h) / 36.0
    pad = fs * 0.9
    sb_m = _nice_length(metres_per_px * w * 0.22)
    sb_px = sb_m / metres_per_px if metres_per_px else 0
    txt = f"{sb_m / 1000:g} km" if sb_m >= 1000 else f"{sb_m:g} m"
    p = [f"<svg class='map' viewBox='0 0 {w} {h}' role='img' aria-label='{html.escape(label)}' "
         f"xmlns='http://www.w3.org/2000/svg'>",
         f"<image width='{w}' height='{h}' href='data:{b64}'/>"]
    if sb_px > 0:
        x0, y0 = pad, h - pad
        bw = sb_px + 2 * fs * 0.6
        p.append(f"<rect x='{x0 - fs * 0.6:.1f}' y='{y0 - fs * 2.3:.1f}' width='{bw:.1f}' height='{fs * 2.9:.1f}' "
                 "rx='4' fill='#fff' fill-opacity='.85'/>")
        p.append(f"<rect x='{x0:.1f}' y='{y0 - fs * 0.55:.1f}' width='{sb_px / 2:.1f}' height='{fs * 0.4:.1f}' fill='#222'/>"
                 f"<rect x='{x0 + sb_px / 2:.1f}' y='{y0 - fs * 0.55:.1f}' width='{sb_px / 2:.1f}' height='{fs * 0.4:.1f}' "
                 "fill='#fff' stroke='#222' stroke-width='1'/>")
        p.append(f"<text x='{x0:.1f}' y='{y0 - fs * 0.95:.1f}' font-size='{fs:.1f}' fill='#222' "
                 f"font-family='system-ui,sans-serif'>{txt}</text>")
    ax, ay, s = w - pad - fs * 0.9, pad + fs * 0.4, fs
    if inset:
        ax = pad + s * 0.9
    p.append(f"<g transform='translate({ax:.1f},{ay:.1f})'><circle cx='0' cy='{s * 1.25:.1f}' r='{s * 1.35:.1f}' "
             "fill='#fff' fill-opacity='.85'/>"
             f"<path d='M0 {s * 0.2:.1f} L{s * 0.55:.1f} {s * 1.9:.1f} L0 {s * 1.5:.1f} L{-s * 0.55:.1f} {s * 1.9:.1f}Z' "
             "fill='#222'/>"
             f"<text x='0' y='{s * 2.35 + s * 0.3:.1f}' font-size='{s * 0.8:.1f}' text-anchor='middle' fill='#222' "
             "font-family='system-ui,sans-serif' font-weight='600'>N</text></g>")
    if inset:
        ib64, iw, ih, (fx, fy, fw, fh) = inset
        k = 0.3 * max(w, h) / max(iw, ih)          # inset at 30 % of the map, whatever the map's pixel size
        iw, ih, fx, fy, fw, fh = iw * k, ih * k, fx * k, fy * k, fw * k, fh * k
        ix, iy = w - iw - pad, pad
        bw_ = max(1.5, fs * 0.2)
        p.append(f"<rect x='{ix - bw_:.1f}' y='{iy - bw_:.1f}' width='{iw + 2 * bw_:.1f}' height='{ih + 2 * bw_:.1f}' "
                 "fill='#fff'/>"
                 f"<image x='{ix:.1f}' y='{iy:.1f}' width='{iw:.1f}' height='{ih:.1f}' href='data:{ib64}'/>"
                 f"<rect x='{ix + fx:.1f}' y='{iy + fy:.1f}' width='{max(fw, 3):.1f}' height='{max(fh, 3):.1f}' "
                 f"fill='none' stroke='#d7301f' stroke-width='{bw_:.1f}'/>")
    p.append("</svg>")
    return "".join(p)


def maps(depth, flood, dem, grid, valid=None, other=None, prob=None, backdrop=None, max_px=MAP_PX,
         doubtful=None):
    """Build the report maps. Returns dict with 'depth' (svg), optional 'prob' (svg),
    'backdrop' (label) and 'frame' (description)."""
    rows, cols = flood.shape
    r0, r1, c0, c1, zoomed, held = flood_frame(flood)
    h, w = r1 - r0, c1 - c0
    step = max(1, int(math.ceil(max(h, w) / float(max_px))))
    sl = (slice(r0, r1), slice(c0, c1))
    fl = _block(flood[sl], step, "frac")
    fany = fl > 0
    fmaj = fl >= 0.25
    dd = _block(np.where(flood[sl], depth[sl], np.nan), step, "mean")
    dem_d = _block(dem[sl], step, "mean") if dem is not None else None
    vd = _block(valid[sl], step, "frac") >= 0.5 if valid is not None else np.ones(fl.shape, bool)
    od = _block(other[sl], step, "frac") >= 0.5 if other is not None else None
    dh, dw = fl.shape
    sub = _sub_grid(grid, r0, c0, dh, dw, step)
    rgb, blabel = None, "relief (DEM hillshade, elevation tint)"
    if backdrop is not None:
        try:
            bd = backdrop(sub)
            if bd is not None:
                rgb, blabel = bd
        except Exception:  # a backdrop must never fail the report
            rgb = None
    base = _base_image(dem_d, grid.cell_x * step, rgb)
    nod = ~vd & ~fany
    hatch = _hatch(nod.shape) & nod
    common = [(nod, NODATA_RGB, 0.35), (hatch, (150, 150, 150), 0.8)]
    if od is not None:
        common.append((od & ~fany, PERMANENT_RGB, 0.85))
    dcol = ramp(np.nan_to_num(dd), DEPTH_RAMP)
    dbt = None
    if doubtful is not None and np.any(doubtful):
        dbt = (_block(doubtful[sl], step, "frac") >= 0.5) & fany
    img = _compose(base, common + [(fany & np.isfinite(dd), dcol, 0.9), (dbt, DOUBTFUL_RGB, 0.85),
                                   (_outline(fmaj) if step > 1 else _outline(fany), OUTLINE_RGB, 0.9)])
    mpp = grid.cell_x * step
    inset = None
    if zoomed:
        istep = max(1, int(math.ceil(max(rows, cols) / float(INSET_PX))))
        idem = _block(dem, istep, "mean") if dem is not None else np.zeros((1, 1))
        ifl = _block(flood, istep, "any")
        irgb = None
        if backdrop is not None:
            try:
                bd = backdrop(_sub_grid(grid, 0, 0, ifl.shape[0], ifl.shape[1], istep))
                irgb = bd[0] if bd is not None else None
            except Exception:
                irgb = None
        ibase = _base_image(idem, grid.cell_x * istep, irgb)
        iimg = _compose(ibase, [(ifl, OUTLINE_RGB, 0.95)])
        ih, iw = iimg.shape[:2]
        sc = float(INSET_PX) / max(ih, iw)
        dispw, disph = iw * sc, ih * sc
        fx, fy = c0 / istep * sc, r0 / istep * sc
        fw, fh = w / istep * sc, h / istep * sc
        inset = (_png_b64(iimg), round(dispw), round(disph), (fx, fy, fw, fh))
    out = {"backdrop": blabel, "zoomed": zoomed, "held": held, "doubtful": bool(dbt is not None and dbt.any()),
           "frame_km": (w * grid.cell_x / 1000.0, h * grid.cell_y / 1000.0),
           "depth": _svg_map(_png_b64(img), dw, dh, mpp, inset, "Flood depth map")}
    if prob is not None:
        pd = _block(prob[sl], step, "max")
        pm = np.isfinite(pd) & (pd >= 0.05)
        pcol = ramp(np.nan_to_num(pd), PROB_RAMP)
        pimg = _compose(base, common + [(pm, pcol, 0.9), (_outline(fmaj) if step > 1 else _outline(fany),
                                                           OUTLINE_RGB, 0.6)])
        out["prob"] = _svg_map(_png_b64(pimg), dw, dh, mpp, None, "Flood probability map")
    return out


def quicklook_png(depth, flood, dem, grid, max_px=900):
    """Kept for compatibility: the depth map as a bare base64 PNG (no frame furniture)."""
    m = maps(depth, flood, dem, grid, max_px=max_px)
    found = re.search(r"base64,([A-Za-z0-9+/=]+)", m["depth"])  # JPEG since 0.4.2
    return found.group(1) if found else ""


# ---------------------------------------------------------------- charts and tables

def classes_svg(classes):
    if not classes:
        return ""
    w, bh, gap, lw = 520, 22, 8, 90
    mx = max(c["share"] for c in classes) or 1.0
    h = len(classes) * (bh + gap)
    parts = [f'<svg viewBox="0 0 {w} {h}" width="100%" style="max-width:{w}px" role="img" '
             'aria-label="Flooded area by depth class">']
    for i, c in enumerate(classes):
        y = i * (bh + gap)
        bw = (w - lw - 140) * c["share"] / mx
        mid = (DEPTH_RAMP[min(i, len(DEPTH_RAMP) - 1)][1])
        parts.append(f'<text x="{lw - 8}" y="{y + bh * 0.7}" text-anchor="end" class="ax">{html.escape(c["label"])}</text>')
        parts.append(f'<rect x="{lw}" y="{y}" width="{max(bw, 1):.1f}" height="{bh}" rx="3" '
                     f'fill="rgb{mid}"/>')
        parts.append(f'<text x="{lw + bw + 6:.1f}" y="{y + bh * 0.7}" class="ax">'
                     f'{100 * c["share"]:.1f} % · {c["area_km2"]:,.2f} km²</text>')
    parts.append("</svg>")
    return "".join(parts)


def _gradient_legend(stops, unit, title, top_label=None):
    vmax = stops[-1][0]
    css = ", ".join(f"rgb{c} {100 * v / vmax:.0f}%" for v, c in stops)
    ticks = "".join(f"<span style='left:{100 * v / vmax:.0f}%'>{v:g}</span>" for v, _ in stops)
    return (f"<div class='lg'><div class='lgt'>{html.escape(title)}</div>"
            f"<div class='bar' style='background:linear-gradient(90deg, {css})'></div>"
            f"<div class='ticks'>{ticks}<span class='u'>{html.escape(top_label or unit)}</span></div></div>")


def _swatches(items):
    return "<div class='sw'>" + "".join(
        f"<span><i style='{st}'></i>{html.escape(t)}</span>" for t, st in items) + "</div>"


def _row(k, v):
    return f"<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>"


def _fmt(v, nd=3):
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:,.{nd}f}"
    return str(v)


def _duration(s):
    try:
        s = float(s)
    except (TypeError, ValueError):
        return ""
    if s < 90:
        return f"{s:.0f} s"
    m, sec = divmod(int(round(s)), 60)
    if m < 60:
        return f"{m} min {sec:02d} s"
    h, m = divmod(m, 60)
    return f"{h} h {m:02d} min"


_RENDERED = re.compile(r"^(pre-event: )?(.+?): rendered '([^']+)' product")


def group_notes(warnings, input_notes):
    """Split warnings into (results, inputs): result warnings first; input notes
    deduplicated, with skipped rendered products collapsed into one line each
    for the post- and pre-event images."""
    seen, res, inp = set(), [], []
    skipped = {"": [], "pre-event: ": []}
    inset = set(input_notes or [])
    for w in warnings or []:
        if w in seen:
            continue
        seen.add(w)
        m = _RENDERED.match(w)
        if m:
            skipped[m.group(1) or ""].append(m.group(2))
            continue
        if w.startswith("Pre-event image is a rendered"):
            res.append((w, None))
            continue
        if w.startswith("Pre-event image not used"):
            res.append(("Pre-event image not used: its folder holds no raw water-index bands (green, NIR, "
                        "SWIR1), so permanent water was not taken from it.", w))
            continue
        (inp if w in inset else res).append((w, None))
    for pre, files in skipped.items():
        if files:
            lab = "pre-event image" if pre else "image"
            inp.append((f"{len(files)} rendered product(s) in the {lab} files skipped - they hold display "
                        "colours, not reflectance. Only the Raw bands are used.", "\n".join(files)))
    return res, inp


def _note_html(text, detail, cls="warn"):
    if detail:
        return (f"<div class='{cls}'>{html.escape(text)}<details><summary>details</summary>"
                f"<pre>{html.escape(detail)}</pre></details></div>")
    return f"<div class='{cls}'>{html.escape(text)}</div>"


def _files_html(paths):
    if not paths:
        return "-"
    if isinstance(paths, str):
        paths = [paths]
    dirs = {}
    for p in paths:
        p = str(p).replace("\\", "/")
        dirs.setdefault(os.path.dirname(p), []).append(os.path.basename(p))
    out = []
    for d, names in dirs.items():
        if len(names) == 1:
            out.append(f"<div><b>{html.escape(names[0])}</b><div class='path'>{html.escape(d)}</div></div>")
        else:
            out.append(f"<details class='compact'><summary>{len(names)} files in <span class='path'>{html.escape(d)}</span>"
                       f"</summary><ul class='files'>" + "".join(f"<li>{html.escape(n)}</li>" for n in names)
                       + "</ul></details>")
    return "".join(out)


CSS = """
:root{--page:#eef0f3;--card:#ffffff;--fg:#1d2329;--mute:#5d6873;--line:#e2e5e9;--acc:#08519c;
--warn:#6f4a00;--warnbg:#fff6e0;--warnline:#f1d79b;--info:#33414f;--infobg:#f3f5f8;--infoline:#dde2e8;--shadow:0 1px 2px rgba(16,24,40,.06)}
@media (prefers-color-scheme: dark){:root{--page:#0f1316;--card:#171c21;--fg:#e6e8ea;--mute:#9aa4ad;--line:#2a3138;
--acc:#6baed6;--warn:#f0c674;--warnbg:#2b2414;--warnline:#4d3f1c;--info:#c6ced6;--infobg:#1d2329;--infoline:#2c343c;--shadow:none}}
*{box-sizing:border-box}
html{background:var(--page)}
body{background:var(--page);color:var(--fg);font:15px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;margin:0}
header{background:var(--card);border-bottom:1px solid var(--line)}
.wrap{max-width:1200px;margin:0 auto;padding:0 16px}
header .wrap{padding-top:22px;padding-bottom:18px}
h1{font-size:23px;margin:0 0 4px;letter-spacing:-.01em}
.meta{color:var(--mute);margin:0;font-size:14px}.meta b{color:var(--fg);font-weight:600}
main.wrap{padding-top:18px;padding-bottom:56px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:10px;margin-top:16px}
.kpi{background:var(--page);border-radius:10px;padding:10px 14px}
.kpi b{display:block;font-size:22px;font-variant-numeric:tabular-nums;line-height:1.25}
.kpi span{color:var(--mute);font-size:13px}.kpi small{display:block;color:var(--mute);font-size:12px}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:16px 18px;box-shadow:var(--shadow);min-width:0}
.card+.card,.cols+.card,.card+.cols,.cols+.cols{margin-top:16px}
.cols>.card{margin-top:0}
.cols{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,460px),1fr));gap:16px;align-items:start}
h2{font-size:16px;margin:0 0 12px}
h3{font-size:14px;margin:16px 0 6px;color:var(--mute);font-weight:600}
.sub{color:var(--mute);font-size:13.5px;margin:8px 0 0}
table{border-collapse:collapse;width:100%;font-size:14px}
th,td{text-align:left;padding:5px 8px;border-bottom:1px solid var(--line);vertical-align:top}
tr:last-child th,tr:last-child td{border-bottom:0}
th{color:var(--mute);font-weight:500;width:42%}
td{font-variant-numeric:tabular-nums;overflow-wrap:anywhere}
.warn,.note{border-radius:8px;padding:8px 12px;margin:6px 0;font-size:14px;border:1px solid}
.warn{background:var(--warnbg);color:var(--warn);border-color:var(--warnline)}
.note{background:var(--infobg);color:var(--info);border-color:var(--infoline)}
details summary{cursor:pointer;color:var(--mute);font-size:13px}
pre{white-space:pre-wrap;font-size:12.5px;margin:6px 0 0}
.maps{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,440px),1fr));gap:16px}
figure{margin:0}
svg.map{width:100%;height:auto;display:block;border-radius:8px;border:1px solid var(--line);background:#f6f6f4}
figcaption{font-size:13px;color:var(--mute);margin-top:8px}
.ax{font-size:12px;fill:var(--fg)}
.lg{margin-top:10px;max-width:360px}.lgt{font-size:12.5px;color:var(--mute);margin-bottom:3px}
.bar{height:10px;border-radius:3px}
.ticks{position:relative;height:16px;font-size:11.5px;color:var(--mute)}
.ticks span{position:absolute;transform:translateX(-50%);top:1px}.ticks span:first-child{transform:none}
.ticks span.u{display:none}
.sw{display:flex;gap:14px;flex-wrap:wrap;font-size:12.5px;color:var(--mute);margin-top:8px}
.sw i{display:inline-block;width:14px;height:10px;border-radius:2px;margin-right:5px;vertical-align:-1px}
.path{color:var(--mute);font-size:12.5px;overflow-wrap:anywhere}
ul.files{margin:6px 0 0;padding-left:18px;font-size:13px}
ol,ul{padding-left:20px;margin:0}li{margin:3px 0}
code{font-size:13px}
@page{size:A4 portrait;margin:12mm 11mm 14mm}
@media print{
:root{--page:#fff;--card:#fff;--fg:#111;--mute:#555;--line:#d4d7db;--acc:#08519c;--warn:#5a3c00;--warnbg:#fff6e0;
--warnline:#e6c77d;--info:#333;--infobg:#f4f5f7;--infoline:#d8dce1;--shadow:none}
*{-webkit-print-color-adjust:exact;print-color-adjust:exact}
html,body{background:#fff;font-size:9.5pt;line-height:1.4}
.wrap{max-width:none;padding:0}
header{border:0;background:#fff}header .wrap{padding:0 0 3mm}
h1{font-size:15pt}h2{font-size:11pt;margin-bottom:2mm}h3{font-size:9.5pt}.meta{font-size:8.5pt}
.kpis{grid-template-columns:repeat(3,1fr);gap:2mm;margin-top:3mm}
.kpi{padding:2mm 3mm;border:1px solid var(--line)}.kpi b{font-size:13pt}.kpi span,.kpi small{font-size:7.5pt}
main.wrap{padding:0}
.card{border:1px solid var(--line);border-radius:2mm;padding:3mm 4mm;box-shadow:none}
.card+.card,.cols+.card,.card+.cols,.cols+.cols{margin-top:3mm}
.cols{grid-template-columns:1fr 1fr;gap:3mm}
h2,h3{break-after:avoid}
figure.prob,.pg{break-before:page}
figure.prob svg.map{max-height:105mm}
.cols>.card,.card.keep,.warn,.note,figure,tr{break-inside:avoid}
.maps{grid-template-columns:1fr;gap:4mm}
section.mapcard{break-inside:auto}
svg.map{max-height:118mm;width:auto;max-width:100%;margin:0 auto}
table{font-size:8.5pt}th,td{padding:.6mm 1.5mm}
.sub,figcaption,.path,ul.files,.sw,.lgt{font-size:7.5pt}
.ticks{font-size:7pt;height:11px}
details>summary{list-style:none;color:var(--fg)}details>summary::-webkit-details-marker{display:none}
.printonly{display:block!important}
.warn,.note{font-size:8.5pt;padding:1.5mm 2.5mm}
}
.printonly{display:none}
"""

LIMITS = [
    "Depths are terrain-based estimates from the flood edge, not a hydraulic simulation: no flow, "
    "velocity, timing or backwater is modelled.",
    "The water surface is interpolated from the terrain at the flood edge. An extent that is too large "
    "(wet soil, shadow, vegetation classed as water) raises the surface; one that is too small lowers it.",
    "Depth accuracy is bounded by the DEM: global 30 m DEMs (Copernicus, ALOS, TanDEM-X) are surface "
    "models with metre-level vertical error and smoothed channels; channel depths are underestimated.",
    "The image shows the water at the overpass time, which may not be the peak.",
    "Permanent water is excluded only where a pre-event image or a permanent-water layer was given.",
    "A flood TRACE (vegetation loss between two dates) is not a water extent: it mixes inundation, scour "
    "and deposition. Where the ground was eroded the pre-event DEM no longer describes it, and the trace "
    "edge is an erosion edge. Depths there are indicative; the evidence and uncertainty layers show where.",
    "The uncertainty bands cover DEM error and extent / threshold choice. They do not cover errors in the "
    "method itself (e.g. a water surface that is not level across the extent) or the image timing.",
]


# ---------------------------------------------------------------- the report

def write_report(path, info, depth=None, flood=None, dem=None, grid=None, valid=None, other=None,
                 prob=None, backdrop=None, doubtful=None):
    d = info.get("depth", {})
    cls = info.get("classification")
    s = info.get("surface", {})
    o = info.get("options", {})
    di = info.get("depth_interior", {})
    u = info.get("uncertainty")
    nan = float("nan")

    def rng(key, scale=1.0, nd=2, unit=""):
        if not u or key not in u:
            return ""
        a = u[key]
        return f"P10-P90 {a[0] * scale:,.{nd}f} - {a[2] * scale:,.{nd}f}{unit}"

    k = [("Flooded area", f"{d.get('area_km2', 0):,.2f} km²", rng("area_km2_p10_p50_p90", 1, 1, " km²")),
         ("Mean depth, interior", f"{di.get('mean_m', nan):.2f} m", "flood without its edge ring"),
         ("Mean depth, all cells", f"{d.get('mean_m', nan):.2f} m", rng("mean_depth_m_p10_p50_p90", 1, 2, " m")),
         ("90th percentile depth", f"{d.get('p90_m', nan):.2f} m", ""),
         ("Maximum depth", f"{d.get('max_m', nan):.2f} m", ""),
         ("Volume", f"{d.get('volume_m3', 0) / 1e6:,.2f} Mm³", rng("volume_m3_p10_p50_p90", 1e-6, 1, " Mm³"))]
    method = "HAND" if o.get("method") == "hand" else "boundary elevation (DEM)"
    when = (cls or {}).get("date")
    sensor = (cls or {}).get("sensor")
    meta = []
    if when:
        meta.append(f"Image <b>{html.escape(str(when))}</b>" + (f" ({html.escape(str(sensor))})" if sensor else ""))
    meta += [f"{html.escape(method)} method", f"run {html.escape(info.get('started', ''))}",
             f"{_duration(info.get('seconds'))}", f"FloodDepth {html.escape(info.get('version', ''))}"]
    h = ["<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<meta name='viewport' content='width=device-width,initial-scale=1'>",
         "<title>FloodDepth run report</title><style>", CSS, "</style></head><body>",
         "<header><div class='wrap'><h1>Flood depth estimate</h1>",
         "<p class='meta'>" + " · ".join(x for x in meta if x) + "</p><div class='kpis'>"]
    h += [f"<div class='kpi'><b>{html.escape(v)}</b><span>{html.escape(n)}</span>"
          + (f"<small>{html.escape(sm)}</small>" if sm else "") + "</div>" for n, v, sm in k]
    h.append("</div></div></header><main class='wrap'>")

    res, inp = group_notes(info.get("warnings"), (cls or {}).get("notes"))
    if res:
        h.append("<section class='card'><h2>Check these</h2>")
        h += [_note_html(t, dt, "warn") for t, dt in res]
        h.append("</section>")

    # maps
    if depth is not None and grid is not None and flood is not None:
        try:
            m = maps(depth, flood, dem, grid, valid=valid, other=other, prob=prob, backdrop=backdrop,
                     doubtful=doubtful)
            fw, fhh = m["frame_km"]
            sw = [("flood edge", f"background:rgb{OUTLINE_RGB}")]
            if other is not None and other.any():
                sw.append(("permanent water", f"background:rgb{PERMANENT_RGB}"))
            if m.get("doubtful"):
                sw.append(("doubtful: ground above the water surface", f"background:rgb{DOUBTFUL_RGB}"))
            sw.append(("no data / cloud", "background:repeating-linear-gradient(45deg,#999 0 2px,#ddd 2px 6px)"))
            where = (f"{fw:,.1f} x {fhh:,.1f} km around the main flood area ({100 * m['held']:.0f} % of the "
                     "flooded cells); the inset shows the whole analysis area (flood in dark blue, frame in red)"
                     if m["zoomed"] else f"the whole analysis area, {fw:,.1f} x {fhh:,.1f} km")
            h.append("<section class='card mapcard'><h2>Maps</h2><div class='maps'>")
            h.append("<figure>" + m["depth"] + _gradient_legend(DEPTH_RAMP, "m", "Depth (m)", "m")
                     + _swatches(sw) + f"<p class='sub'>Background: {html.escape(m['backdrop'])}. "
                     f"Frame: {html.escape(where)}.</p></figure>")
            if "prob" in m:
                h.append("<figure class='prob'>" + m["prob"]
                         + _gradient_legend(PROB_RAMP, "", "Flood probability (share of ensemble members)", " ")
                         + "<p class='sub'>Where the probability is low the extent itself is uncertain (threshold, "
                         "edge position); depths there matter less than whether it flooded at all.</p></figure>")
            h.append("</div></section>")
        except Exception as exc:  # the report must not fail the run
            h.append(f"<section class='card'><p class='sub'>Maps not drawn: {html.escape(str(exc))}</p></section>")

    # depth + uncertainty
    h.append("<div class='cols'><section class='card'><h2>Depth</h2>")
    if d.get("classes"):
        h.append(classes_svg(d["classes"]))
    h.append("<h3>Statistics</h3><table>")
    h.append("<tr><th></th><td><b>All flood cells</b></td><td><b>Interior</b></td></tr>")
    for key, lab, nd in (("cells", "Cells", 0), ("mean_m", "Mean (m)", 2), ("median_m", "Median (m)", 2),
                         ("p90_m", "90th percentile (m)", 2), ("p99_m", "99th percentile (m)", 2),
                         ("max_m", "Maximum (m)", 2), ("volume_m3", "Volume (m³)", 0)):
        if key in d:
            h.append(f"<tr><th>{lab}</th><td>{_fmt(d.get(key), nd)}</td><td>{_fmt(di.get(key), nd)}</td></tr>")
    if "zero_depth_share" in d:
        h.append(f"<tr><th>Share at 0 m (ground above the surface)</th><td>{100 * d['zero_depth_share']:.1f} %</td>"
                 "<td></td></tr>")
    h.append("</table>")
    if di.get("cells"):
        h.append(f"<p class='sub'>Interior = flood cells not on the waterline ({100 * di.get('share_of_flood', 0):.0f} % "
                 "of the flood). Depth on the waterline is about 0 by construction, so in a fragmented extent the "
                 "all-cell mean is dominated by edges.</p>")
    h.append("</section>")
    if u:
        h.append("<section class='card'><h2>Uncertainty</h2><table>")
        a, md, vv = u["area_km2_p10_p50_p90"], u["mean_depth_m_p10_p50_p90"], u["volume_m3_p10_p50_p90"]
        h.append("<tr><th></th><td><b>P10</b></td><td><b>P50</b></td><td><b>P90</b></td></tr>")
        h.append(f"<tr><th>Flooded area (km²)</th><td>{a[0]:,.2f}</td><td>{a[1]:,.2f}</td><td>{a[2]:,.2f}</td></tr>")
        h.append(f"<tr><th>Mean depth (m)</th><td>{md[0]:.2f}</td><td>{md[1]:.2f}</td><td>{md[2]:.2f}</td></tr>")
        h.append(f"<tr><th>Volume (Mm³)</th><td>{vv[0] / 1e6:,.2f}</td><td>{vv[1] / 1e6:,.2f}</td>"
                 f"<td>{vv[2] / 1e6:,.2f}</td></tr></table><h3>Ensemble</h3><table>")
        for lab, v in (("Members (succeeded)", f"{u['members']} ({u['succeeded']})"),
                       ("DEM error (1 sigma)", f"{u['sigma_dem_m']:g} m - {u.get('sigma_source', '')}"),
                       ("DEM error correlation length", f"{u['corr_length_m']:g} m"),
                       ("Method error (edge cross-validation)",
                        "-" if u.get("method_sd_m") is None else f"{u['method_sd_m']:.2f} m"),
                       ("Extra error in erosion / deposition cells", f"{u.get('sigma_erosion_m', 0):g} m"),
                       ("Extent perturbation", u.get("extent_perturbation")),
                       ("Computed at", "full resolution" if not u.get("coarsening") or u.get("coarsening") == 1
                        else f"{u['coarsening']} x coarser grid"),
                       ("Random seed", u.get("seed"))):
            h.append(_row(lab, v))
        h.append("</table><p class='sub'>Monte Carlo ensemble: correlated DEM error, threshold / extent "
                 "perturbation" + (", extra terrain error in erosion / deposition cells" if u.get("sigma_erosion_m")
                                   else "")
                 + ". P10-P90 is the range 80 % of members fall within. Rasters: flood_probability.tif, "
                 "depth_p10.tif, depth_p50.tif, depth_p90.tif.</p></section>")
    h.append("</div>")

    # evidence + surface
    ev = info.get("evidence")
    h.append("<div class='cols pg'>")
    if ev:
        h.append("<section class='card'><h2>What the extent is based on</h2><table>")
        h.append("<tr><th>Evidence</th><td><b>km²</b></td><td><b>Mean depth (m)</b></td>"
                 "<td><b>P90-P10 (m)</b></td><td><b>P(flood)</b></td></tr>")
        for r in ev:
            h.append(f"<tr><th>{html.escape(r['label'])}</th><td>{r['area_km2']:,.2f}</td>"
                     f"<td>{r['mean_depth_m']:.2f}</td><td>{_fmt(r.get('mean_band_m'), 2)}</td>"
                     f"<td>{_fmt(r.get('mean_prob'), 2)}</td></tr>")
        h.append("</table></section>")
    h.append("<section class='card'><h2>Water surface</h2><table>")
    for lab, key in (("Flood cells", "flood_cells"), ("Waterline cells", "waterline_cells"),
                     ("Waterline used", "waterline_used"), ("Dropped as outliers", "waterline_dropped"),
                     ("Edge cells not used (NoData / grid edge)", "edge_cells_excluded"),
                     ("Outlier limit (m)", "outlier_limit"), ("Mean edge lift (m)", "edge_lift_mean_m"),
                     ("Cells clipped to 0 m", "zero_clipped_cells"),
                     ("Cells not reached by interpolation", "unfilled_cells")):
        if key in s:
            v = s[key]
            h.append(_row(lab, _fmt(v, 2) if isinstance(v, float) else f"{v:,}" if isinstance(v, int) else v))
    h.append("</table></section></div>")

    # classification
    if cls:
        t = cls.get("threshold", {})
        h.append("<section class='card'><h2>Classification</h2><div class='cols'><div><table>")
        rows1 = (("Mode", cls.get("mode")), ("Sensor", cls.get("sensor")), ("Date", cls.get("date")),
                 ("Index", cls.get("index")),
                 ("Threshold", f"{t.get('threshold', nan):.4f} ({t.get('method')})"),
                 ("Separability eta", "-" if t.get("eta") is None else f"{t['eta']:.2f} ({cls.get('confidence')})"),
                 ("Near-optimal threshold band", "-" if not cls.get("threshold_band") else
                  f"{cls['threshold_band'][0]:.4f} to {cls['threshold_band'][1]:.4f}"))
        rows2 = (("Water share (before / after cleanup)",
                  f"{100 * cls.get('water_share_raw', 0):.2f} % / {100 * cls.get('water_share', 0):.2f} %"),
                 ("Cloud masked", f"{100 * cls.get('cloud_share', 0):.1f} %"),
                 ("Cleanup", ", ".join(f"{a}={b}" for a, b in cls.get("cleanup", {}).items())),
                 ("Pre-event image", cls.get("pre_bands") or "-"),
                 ("Permanent water share (pre-event image)", _fmt(cls.get("permanent_share"), 4)))
        h += [_row(a, b) for a, b in rows1 if b not in (None, "")]
        h.append("</table></div><div><table>")
        h += [_row(a, b) for a, b in rows2 if b not in (None, "")]
        h.append("</table></div></div>")
        if cls.get("bands"):
            h.append("<details><summary>Bands used</summary><ul class='files'>" + "".join(
                f"<li>{html.escape(b.strip())}</li>" for b in str(cls['bands']).split(";")) + "</ul></details>")
        if inp:
            h.append("<h3>Input notes</h3>")
            h += [_note_html(tx, dt, "note") for tx, dt in inp]
        h.append("</section>")

    # inputs + log
    g = info.get("grid", {})
    h.append("<div class='cols pg'><section class='card'><h2>Inputs and analysis grid</h2><table>")
    h.append(_row("Analysis grid", f"{g.get('cols', 0):,} x {g.get('rows', 0):,} cells of {g.get('cell_x', 0):g} m "
                  f"({(g.get('cols') or 0) * (g.get('rows') or 0) / 1e6:.1f} M cells), {g.get('crs')}"))
    for key, lab in (("dem_path", "DEM"), ("hand_path", "HAND raster"), ("image_paths", "Image files"),
                     ("pre_paths", "Pre-event image files"), ("permanent_raster", "Permanent water raster")):
        if o.get(key):
            h.append(f"<tr><th>{lab}</th><td>{_files_html(o[key])}</td></tr>")
    for key, lab in (("extent_wkts", "Flood extent polygons"), ("aoi_wkts", "Area of interest"),
                     ("seed_wkts", "Connectivity features"), ("permanent_wkts", "Permanent water polygons")):
        if o.get(key):
            h.append(_row(lab, o[key]))
    params = [(kk, o.get(kk)) for kk in ("cell", "min_water_cells", "max_hole_cells", "outlier_radius", "outlier_k",
                                         "outlier_min", "smooth_radius", "smooth_tolerance", "edge_step",
                                         "edge_step_cap", "interpolation", "min_fraction", "uncertainty_members")
              if o.get(kk) not in (None, "")]
    if params:
        h.append("<tr><th>Parameters</th><td><details class='compact'><summary>" + f"{len(params)} settings</summary><table>"
                 + "".join(_row(a, b) for a, b in params) + "</table></details></td></tr>")
    h.append("</table></section><section class='card'><h2>Processing log</h2><ol>")
    h += [f"<li>{html.escape(x)}</li>" for x in info.get("steps", [])]
    h.append("</ol>")
    tm = info.get("timings")
    if tm:
        tot = sum(tm.values()) or 1.0
        h.append("<h3>Time by stage</h3><table>" + "".join(
            f"<tr><th>{html.escape(k)}</th><td>{_duration(v)} ({100 * v / tot:.0f} %)</td></tr>"
            for k, v in tm.items()) + "</table>")
    h.append("</section></div>")
    h.append("<section class='card keep'><details><summary><b>Limitations</b> - read before using the depths</summary><ul>")
    h += [f"<li>{html.escape(x)}</li>" for x in LIMITS]
    h.append("</ul></details></section>")
    h.append("<p class='printonly sub'>FloodDepth run report - output folder: "
             + html.escape(str(o.get("out_dir") or "")) + "</p></main>")
    # print: every collapsed section opens, then closes again afterwards
    h.append("<script>(function(){var o=[];addEventListener('beforeprint',function(){o=[];"
             "document.querySelectorAll('details:not([open]):not(.compact)').forEach(function(d){d.open=true;o.push(d)})});"
             "addEventListener('afterprint',function(){o.forEach(function(d){d.open=false})})})();</script>")
    h.append("</body></html>")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(h))
    return path
