# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Self-contained HTML run report (inline PNG quicklook and SVG chart)."""

import base64
import html
import math

import numpy as np

from . import raster as R

# sequential blue ramp for depth (light = shallow, dark = deep)
DEPTH_RAMP = [(0.0, (198, 219, 239)), (0.5, (107, 174, 214)), (1.0, (49, 130, 189)),
              (2.0, (8, 81, 156)), (3.0, (8, 48, 107))]


def hillshade(dem, cell, az=315.0, alt=45.0):
    z = np.where(np.isfinite(dem), dem, np.nanmean(dem) if np.isfinite(dem).any() else 0.0)
    gy, gx = np.gradient(z, cell, cell)
    slope = np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, gy)
    a, b = math.radians(az), math.radians(alt)
    hs = np.sin(b) * np.cos(slope) + np.cos(b) * np.sin(slope) * np.cos(a - aspect)
    return np.clip(hs, 0, 1)


def ramp(values, stops):
    v = np.asarray(values, dtype="float64")
    xs = [s[0] for s in stops]
    out = np.zeros(v.shape + (3,), dtype="float64")
    for ch in range(3):
        out[..., ch] = np.interp(v, xs, [s[1][ch] for s in stops])
    return out


def quicklook_png(depth, flood, dem, grid, max_px=900):
    gdal = R._gdal()
    rows, cols = depth.shape
    step = max(1, int(math.ceil(max(rows, cols) / float(max_px))))
    d = depth[::step, ::step]
    f = flood[::step, ::step]
    hs = hillshade(dem[::step, ::step], grid.cell_x * step)
    base = (0.35 + 0.65 * hs)[..., None] * np.array([235, 232, 225], dtype="float64")
    col = ramp(np.nan_to_num(d), DEPTH_RAMP)
    img = np.where((f & np.isfinite(d))[..., None], 0.15 * base + 0.85 * col, base)
    img[~np.isfinite(dem[::step, ::step])] = (255, 255, 255)
    img = np.clip(img, 0, 255).astype("uint8")
    h, w, _ = img.shape
    mem = gdal.GetDriverByName("MEM").Create("", w, h, 3, gdal.GDT_Byte)
    for i in range(3):
        mem.GetRasterBand(i + 1).WriteArray(img[..., i])
    path = "/vsimem/flooddepth_quicklook.png"
    gdal.GetDriverByName("PNG").CreateCopy(path, mem)
    fp = gdal.VSIFOpenL(path, "rb")
    gdal.VSIFSeekL(fp, 0, 2)
    size = gdal.VSIFTellL(fp)
    gdal.VSIFSeekL(fp, 0, 0)
    data = gdal.VSIFReadL(1, size, fp)
    gdal.VSIFCloseL(fp)
    gdal.Unlink(path)
    return base64.b64encode(data).decode("ascii")


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
        bw = (w - lw - 120) * c["share"] / mx
        mid = (DEPTH_RAMP[min(i, len(DEPTH_RAMP) - 1)][1])
        parts.append(f'<text x="{lw - 8}" y="{y + bh * 0.7}" text-anchor="end" class="ax">{html.escape(c["label"])}</text>')
        parts.append(f'<rect x="{lw}" y="{y}" width="{max(bw, 1):.1f}" height="{bh}" rx="2" '
                     f'fill="rgb{mid}"/>')
        parts.append(f'<text x="{lw + bw + 6:.1f}" y="{y + bh * 0.7}" class="ax">'
                     f'{100 * c["share"]:.1f} % · {c["area_km2"]:.3f} km²</text>')
    parts.append("</svg>")
    return "".join(parts)


def _row(k, v):
    return f"<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>"


def _fmt(v, nd=3):
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:,.{nd}f}"
    return str(v)


CSS = """
:root{--bg:#fbfaf7;--fg:#1d2329;--mute:#5d6873;--line:#e2ded5;--acc:#08519c;--warn:#8a5a00;--warnbg:#fff5df}
@media (prefers-color-scheme: dark){:root{--bg:#15191d;--fg:#e6e8ea;--mute:#9aa4ad;--line:#2c333a;--acc:#6baed6;--warn:#f0c674;--warnbg:#2e2614}}
body{background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;margin:0}
main{max-width:980px;margin:0 auto;padding:28px 16px 60px}
h1{font-size:24px;margin:0 0 4px}h2{font-size:17px;margin:32px 0 10px;border-bottom:1px solid var(--line);padding-bottom:6px}
.sub{color:var(--mute);margin:0 0 20px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px}
.kpi{border:1px solid var(--line);border-radius:8px;padding:10px 12px}
.kpi b{display:block;font-size:22px;font-variant-numeric:tabular-nums}.kpi span{color:var(--mute);font-size:13px}
table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;padding:5px 8px;border-bottom:1px solid var(--line);vertical-align:top}
th{color:var(--mute);font-weight:500;width:38%}
.warn{background:var(--warnbg);color:var(--warn);border-radius:8px;padding:10px 14px;margin:6px 0}
img.ql{width:100%;border-radius:8px;border:1px solid var(--line)}
.ax{font-size:12px;fill:var(--fg)}
.legend{display:flex;gap:10px;flex-wrap:wrap;font-size:13px;color:var(--mute);margin-top:6px}
.legend i{display:inline-block;width:14px;height:10px;border-radius:2px;margin-right:4px;vertical-align:middle}
code{font-size:13px}
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


def write_report(path, info, depth=None, flood=None, dem=None, grid=None):
    d = info.get("depth", {})
    cls = info.get("classification")
    s = info.get("surface", {})
    o = info.get("options", {})
    k = [("Flooded area", f"{d.get('area_km2', 0):,.3f} km²"),
         ("Mean depth", f"{d.get('mean_m', float('nan')):.2f} m"),
         ("90th pct depth", f"{d.get('p90_m', float('nan')):.2f} m"),
         ("Max depth", f"{d.get('max_m', float('nan')):.2f} m"),
         ("Volume", f"{d.get('volume_m3', 0) / 1e6:,.2f} Mm³")]
    h = ["<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<meta name='viewport' content='width=device-width,initial-scale=1'>",
         "<title>FloodDepth run report</title><style>", CSS, "</style></head><body><main>",
         "<h1>FloodDepth run report</h1>",
         f"<p class='sub'>Version {html.escape(info.get('version', ''))} · {html.escape(info.get('started', ''))} · "
         f"{'HAND' if o.get('method') == 'hand' else 'boundary-elevation (DEM)'} method · "
         f"{info.get('seconds', '')} s</p>",
         "<div class='kpis'>"]
    h += [f"<div class='kpi'><b>{html.escape(v)}</b><span>{html.escape(n)}</span></div>" for n, v in k]
    h.append("</div>")
    if info.get("warnings"):
        h.append("<h2>Check these</h2>")
        h += [f"<div class='warn'>{html.escape(w)}</div>" for w in info["warnings"]]
    if depth is not None and grid is not None:
        try:
            b64 = quicklook_png(depth, flood, dem, grid)
            h.append("<h2>Depth map</h2>")
            h.append(f"<img class='ql' alt='Flood depth over hillshade' src='data:image/png;base64,{b64}'>")
            h.append("<div class='legend'>" + "".join(
                f"<span><i style='background:rgb{c}'></i>{v:g} m</span>" for v, c in DEPTH_RAMP) + "</div>")
        except Exception as exc:  # the report must not fail the run
            h.append(f"<p class='sub'>Quicklook not drawn: {html.escape(str(exc))}</p>")
    if d.get("classes"):
        h.append("<h2>Area by depth class</h2>")
        h.append(classes_svg(d["classes"]))
    h.append("<h2>Depth statistics</h2><table>")
    for key, lab in (("cells", "Flood cells"), ("mean_m", "Mean (m)"), ("median_m", "Median (m)"),
                     ("p90_m", "90th percentile (m)"), ("p99_m", "99th percentile (m)"), ("max_m", "Maximum (m)"),
                     ("volume_m3", "Volume (m³)"), ("zero_depth_share", "Share at 0 m (ground above surface)")):
        if key in d:
            h.append(_row(lab, _fmt(d[key], 0 if key in ("volume_m3", "cells") else 3)))
    h.append("</table>")
    u = info.get("uncertainty")
    if u:
        h.append("<h2>Uncertainty</h2>")
        h.append("<p class='sub'>Monte Carlo ensemble: correlated DEM error, threshold / extent perturbation"
                 + (", extra terrain error in erosion / deposition cells" if u.get("sigma_erosion_m") else "")
                 + ". P10-P90 is the range that 80 % of members fall within. Rasters: flood_probability.tif, "
                 "depth_p10.tif, depth_p50.tif, depth_p90.tif.</p><table>")
        a, md, vv = u["area_km2_p10_p50_p90"], u["mean_depth_m_p10_p50_p90"], u["volume_m3_p10_p50_p90"]
        h.append("<tr><th></th><td><b>P10</b></td><td><b>P50</b></td><td><b>P90</b></td></tr>")
        h.append(f"<tr><th>Flooded area (km²)</th><td>{a[0]:.3f}</td><td>{a[1]:.3f}</td><td>{a[2]:.3f}</td></tr>")
        h.append(f"<tr><th>Mean depth (m)</th><td>{md[0]:.2f}</td><td>{md[1]:.2f}</td><td>{md[2]:.2f}</td></tr>")
        h.append(f"<tr><th>Volume (Mm³)</th><td>{vv[0] / 1e6:.2f}</td><td>{vv[1] / 1e6:.2f}</td><td>{vv[2] / 1e6:.2f}</td></tr>")
        h.append("</table><table>")
        for lab, v in (("Members (succeeded)", f"{u['members']} ({u['succeeded']})"),
                       ("DEM error (1 sigma)", f"{u['sigma_dem_m']:g} m - {u.get('sigma_source', '')}"),
                       ("DEM error correlation length", f"{u['corr_length_m']:g} m"),
                       ("Extra error in erosion / deposition cells", f"{u.get('sigma_erosion_m', 0):g} m"),
                       ("Extent perturbation", u.get("extent_perturbation")), ("Random seed", u.get("seed"))):
            h.append(_row(lab, v))
        h.append("</table>")
    ev = info.get("evidence")
    if ev:
        h.append("<h2>What the extent is based on</h2><table>")
        h.append("<tr><th>Evidence</th><td><b>Area km²</b></td><td><b>Mean depth m</b></td>"
                 "<td><b>Mean P90-P10 m</b></td><td><b>Mean probability</b></td></tr>")
        for r in ev:
            h.append(f"<tr><th>{html.escape(r['label'])}</th><td>{r['area_km2']:.3f}</td>"
                     f"<td>{r['mean_depth_m']:.2f}</td><td>{_fmt(r.get('mean_band_m'), 2)}</td>"
                     f"<td>{_fmt(r.get('mean_prob'), 2)}</td></tr>")
        h.append("</table>")
    if cls:
        t = cls.get("threshold", {})
        h.append("<h2>Classification</h2><table>")
        for lab, v in (("Mode", cls.get("mode")), ("Sensor", cls.get("sensor")), ("Date", cls.get("date")),
                       ("Bands", cls.get("bands")), ("Index", cls.get("index")),
                       ("Threshold", f"{t.get('threshold', float('nan')):.4f} ({t.get('method')})"),
                       ("Separability eta", "-" if t.get("eta") is None else f"{t['eta']:.2f} ({cls.get('confidence')})"),
                       ("Tiles pooled", t.get("tiles_used")),
                       ("Global Otsu threshold", _fmt(t.get("global_threshold"), 4)),
                       ("Near-optimal threshold band", "-" if not cls.get("threshold_band") else
                        f"{cls['threshold_band'][0]:.4f} to {cls['threshold_band'][1]:.4f}"),
                       ("Water share (before / after cleanup)",
                        f"{100 * cls.get('water_share_raw', 0):.2f} % / {100 * cls.get('water_share', 0):.2f} %"),
                       ("Cloud masked", f"{100 * cls.get('cloud_share', 0):.2f} %"),
                       ("Cleanup", ", ".join(f"{a}={b}" for a, b in cls.get("cleanup", {}).items())),
                       ("Pre-event image", cls.get("pre_bands") or "-"),
                       ("Permanent water share", _fmt(cls.get("permanent_share"), 4))):
            h.append(_row(lab, v))
        h.append("</table>")
    h.append("<h2>Water surface</h2><table>")
    for lab, key in (("Flood cells", "flood_cells"), ("Waterline cells", "waterline_cells"),
                     ("Waterline used", "waterline_used"), ("Dropped as outliers", "waterline_dropped"),
                     ("Edge cells not used (NoData / grid edge)", "edge_cells_excluded"),
                     ("Outlier limit (m)", "outlier_limit"), ("Cells clipped to 0 m", "zero_clipped_cells"),
                     ("Cells not reached by interpolation", "unfilled_cells")):
        if key in s:
            h.append(_row(lab, _fmt(s[key], 3)))
    h.append("</table>")
    g = info.get("grid", {})
    h.append("<h2>Analysis grid and inputs</h2><table>")
    h.append(_row("Grid", f"{g.get('cols')} x {g.get('rows')} cells of {g.get('cell_x', 0):g} m, {g.get('crs')}"))
    for key in ("dem_path", "hand_path", "image_paths", "pre_paths", "extent_wkts", "aoi_wkts", "seed_wkts",
                "permanent_wkts", "cell", "min_water_cells", "max_hole_cells", "outlier_radius", "outlier_k",
                "outlier_min", "smooth_radius", "min_fraction"):
        v = o.get(key)
        if v not in (None, "", [], 0) or key in ("outlier_k",):
            h.append(_row(key, ", ".join(map(str, v)) if isinstance(v, list) else v))
    h.append("</table><h2>Processing log</h2><ol>")
    h += [f"<li>{html.escape(x)}</li>" for x in info.get("steps", [])]
    h.append("</ol><h2>Limitations</h2><ul>")
    h += [f"<li>{html.escape(x)}</li>" for x in LIMITS]
    h.append("</ul></main></body></html>")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(h))
    return path
