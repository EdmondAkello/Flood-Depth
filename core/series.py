# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Flood depth for several image dates: one quick run per date, then the
maximum observed depth, the number of dates flooded and the first date
flooded on one common grid, with a summary page (series.html) and table
(series.csv).

Image files are grouped by the acquisition date in their names
(YYYY-MM-DD, YYYY_MM_DD or YYYYMMDD, as in Copernicus Browser, SAFE and
Landsat names). Each date runs exactly like *Quick flood depth* in its own
sub-folder, with the same DEM, permanent water, area of interest and
settings, so the dates are comparable. A date that fails (no water left,
refused threshold) is reported and skipped; the others still run.

The common grid is the union of the per-date analysis grids. They share
the DEM lattice and cell size, so the per-date rasters drop onto it
without resampling ("near").
"""

import copy
import csv
import html
import math
import os
import time

import numpy as np

from . import raster as R
from . import sensors as S


def group_by_date(paths):
    """{date: [paths]} sorted by date; raises ValueError naming undated files."""
    groups, undated = {}, []
    for p in paths:
        d = S.acquisition_date(p)
        if d is None:
            undated.append(os.path.basename(str(p)))
        else:
            groups.setdefault(d, []).append(p)
    if undated:
        raise ValueError("No acquisition date in the file name of: " + ", ".join(undated[:8])
                         + (" ..." if len(undated) > 8 else "")
                         + ". Name the files with their date (YYYY-MM-DD), as the Copernicus Browser does.")
    if not groups:
        raise ValueError("No image files were given.")
    return dict(sorted(groups.items()))


def _union_grid(grids):
    g0 = grids[0]
    x0 = min(g.gt[0] for g in grids)
    y0 = max(g.gt[3] for g in grids)
    x1 = max(g.gt[0] + g.cols * g.gt[1] for g in grids)
    y1 = min(g.gt[3] + g.rows * g.gt[5] for g in grids)
    cx, cy = g0.gt[1], g0.gt[5]
    cols = int(round((x1 - x0) / cx))
    rows = int(round((y1 - y0) / cy))
    return R.Grid(cols, rows, (x0, cx, 0.0, y0, 0.0, cy), g0.wkt)


def _line_svg(dates, series, unit, title, w=560, h=190):
    """Small SVG line chart: series = [(label, values, lows, highs, colour)]."""
    n = len(dates)
    vals = [v for _, vs, lo, hi, _ in series for v in list(vs) + list(lo or []) + list(hi or [])
            if v is not None and np.isfinite(v)]
    if not vals or n == 0:
        return ""
    vmax = max(vals) * 1.1 or 1.0
    L, Rm, T, B = 52, 30, 22, 34
    pw, ph = w - L - Rm, h - T - B
    import datetime as _dt
    days = [_dt.date.fromisoformat(d).toordinal() for d in dates]
    span = float(days[-1] - days[0]) or 1.0
    X = (lambda i: L + (pw * (days[i] - days[0]) / span if n > 1 else pw / 2))   # spaced by time
    Y = (lambda v: T + ph * (1 - v / vmax))
    p = [f"<svg viewBox='0 0 {w} {h}' width='100%' style='max-width:{w}px' role='img' "
         f"aria-label='{html.escape(title)}'><text x='{L}' y='14' class='ax' font-weight='600'>"
         f"{html.escape(title)}</text>"]
    for k in range(5):
        v = vmax * k / 4
        p.append(f"<line x1='{L}' x2='{w - Rm}' y1='{Y(v):.1f}' y2='{Y(v):.1f}' stroke='currentColor' "
                 "stroke-opacity='.12'/>"
                 f"<text x='{L - 6}' y='{Y(v) + 4:.1f}' text-anchor='end' class='ax'>{v:,.3g}</text>")
    for i, d in enumerate(dates):
        p.append(f"<text x='{X(i):.1f}' y='{h - 12}' text-anchor='middle' class='ax'>{html.escape(d[5:])}</text>")
    for lab, vs, lo, hi, col in series:
        if lo and hi:
            for i in range(n):
                if lo[i] is not None and hi[i] is not None and np.isfinite(lo[i]) and np.isfinite(hi[i]):
                    p.append(f"<line x1='{X(i):.1f}' x2='{X(i):.1f}' y1='{Y(lo[i]):.1f}' y2='{Y(hi[i]):.1f}' "
                             f"stroke='{col}' stroke-width='6' stroke-opacity='.25' stroke-linecap='round'/>")
        pts = [(X(i), Y(v)) for i, v in enumerate(vs) if v is not None and np.isfinite(v)]
        if len(pts) > 1:
            p.append("<polyline fill='none' stroke='" + col + "' stroke-width='2' points='"
                     + " ".join(f"{x:.1f},{y:.1f}" for x, y in pts) + "'/>")
        p += [f"<circle cx='{x:.1f}' cy='{y:.1f}' r='3.5' fill='{col}'/>" for x, y in pts]
    p.append(f"<text x='{w - Rm}' y='14' text-anchor='end' class='ax'>{html.escape(unit)}</text></svg>")
    return "".join(p)


def run_series(opt, progress=None, log=print):
    """Run the quick workflow for every image date in opt.image_paths.

    Returns a dict with the per-date rows, the common-grid outputs and the
    summary paths. opt.out_dir receives one sub-folder per date."""
    from .pipeline import run_quick, VERSION
    t0 = time.time()
    groups = group_by_date(opt.image_paths)
    dates = list(groups)
    os.makedirs(opt.out_dir, exist_ok=True)

    def say(m):
        if log:
            log(m)

    say(f"{len(dates)} image date(s): " + ", ".join(dates))
    rows, runs = [], []
    for k, d in enumerate(dates):
        o = copy.copy(opt)
        o.image_paths = groups[d]
        o.out_dir = os.path.join(opt.out_dir, d)

        def prog(f, m="", _k=k):
            if progress:
                progress((_k + f) / len(dates) * 0.92, f"{d}: {m}")
        say(f"--- {d} ({len(groups[d])} file(s)) ---")
        row = {"date": d, "files": len(groups[d]), "status": "ok"}
        try:
            info = run_quick(o, progress=prog, log=log)
        except (ValueError, RuntimeError, R.RasterError) as exc:
            row["status"] = "skipped: " + str(exc).split(". ")[0][:200]
            say(f"{d} skipped: {exc}")
            rows.append(row)
            continue
        dep, di, u = info.get("depth", {}), info.get("depth_interior", {}), info.get("uncertainty") or {}
        cl = info.get("classification") or {}
        row.update({"area_km2": dep.get("area_km2"), "mean_depth_m": dep.get("mean_m"),
                    "interior_mean_m": di.get("mean_m"), "p90_depth_m": dep.get("p90_m"),
                    "max_depth_m": dep.get("max_m"), "volume_m3": dep.get("volume_m3"),
                    "area_p10_km2": (u.get("area_km2_p10_p50_p90") or [None])[0],
                    "area_p90_km2": (u.get("area_km2_p10_p50_p90") or [None, None, None])[2],
                    "volume_p10_m3": (u.get("volume_m3_p10_p50_p90") or [None])[0],
                    "volume_p90_m3": (u.get("volume_m3_p10_p50_p90") or [None, None, None])[2],
                    "sensor": cl.get("sensor"), "index": cl.get("index_key"),
                    "cloud_share": cl.get("cloud_share"), "separability": (cl.get("threshold") or {}).get("eta"),
                    "warnings": len(info.get("warnings", [])), "seconds": info.get("seconds"),
                    "report": os.path.join(d, "report.html")})
        rows.append(row)
        runs.append((d, info))
    if not runs:
        raise ValueError("No date produced a flood map: " + "; ".join(f"{r['date']} {r['status']}" for r in rows))

    # common grid: union of the per-date grids (same lattice)
    if progress:
        progress(0.93, "Combining the dates")
    grids = [R.grid_of(R.open_dataset(i["outputs"]["depth"])) for _, i in runs]
    ug = _union_grid(grids)
    maxd = np.full(ug.shape, np.nan, dtype="float32")
    count = np.zeros(ug.shape, dtype="uint8")
    first = np.full(ug.shape, 255, dtype="uint8")
    for k, (d, info) in enumerate(runs):
        a, g, _ = R.read_band(info["outputs"]["depth"], 1)
        a = R.warp_array(a, g, ug, "near")
        f = np.isfinite(a)
        maxd = np.where(f & ~(np.nan_to_num(maxd, nan=-1.0) >= a), a, maxd)
        count += f.astype("uint8")
        first[(first == 255) & f] = k + 1
    tags = {"FLOODDEPTH": VERSION, "DATES": ", ".join(d for d, _ in runs)}
    out_max = os.path.join(opt.out_dir, "max_depth.tif")
    out_cnt = os.path.join(opt.out_dir, "dates_flooded.tif")
    out_first = os.path.join(opt.out_dir, "first_date_flooded.tif")
    R.write_raster(out_max, maxd, ug, -9999.0, "Float32", dict(tags, CONTENT="maximum observed flood depth (m)"))
    R.write_raster(out_cnt, np.where(count > 0, count, 0).astype("uint8"), ug, 0, "Byte",
                   dict(tags, CONTENT="number of image dates on which the cell was flooded"))
    R.write_raster(out_first, first, ug, 255, "Byte",
                   dict(tags, CONTENT="first image date flooded (1 = first date)",
                        VALUES="; ".join(f"{k + 1} {d}" for k, (d, _) in enumerate(runs))))
    ca = ug.cell_area()
    env = {"area_km2": float(np.isfinite(maxd).sum() * ca / 1e6),
           "volume_m3": float(np.nansum(maxd) * ca),
           "mean_depth_m": float(np.nanmean(maxd)) if np.isfinite(maxd).any() else float("nan"),
           "flooded_every_date_km2": float((count == len(runs)).sum() * ca / 1e6)}

    # summary table
    cpath = os.path.join(opt.out_dir, "series.csv")
    cols = ["date", "status", "files", "sensor", "index", "area_km2", "area_p10_km2", "area_p90_km2", "mean_depth_m", "interior_mean_m",
            "p90_depth_m", "max_depth_m", "volume_m3", "volume_p10_m3", "volume_p90_m3", "cloud_share",
            "separability", "warnings", "seconds"]
    with open(cpath, "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(cols)
        for r in rows:
            wr.writerow(["" if r.get(c) is None else (round(r[c], 4) if isinstance(r[c], float) else r[c])
                         for c in cols])

    # summary page
    hpath = os.path.join(opt.out_dir, "series.html")
    dem, _ = R.read_on_grid(opt.dem_path, ug, 1, "bilinear", zero_nodata_auto=True)
    _write_html(hpath, rows, runs, env, maxd, count, dem, ug, opt, VERSION, round(time.time() - t0, 1))
    if progress:
        progress(1.0, "Done")
    say(f"Series done: {len(runs)} of {len(dates)} date(s) mapped; envelope {env['area_km2']:.2f} km2.")
    return {"rows": rows, "envelope": env, "dates": dates,
            "outputs": {"max_depth": out_max, "dates_flooded": out_cnt, "first_date_flooded": out_first,
                        "csv": cpath, "html": hpath},
            "seconds": round(time.time() - t0, 1)}


def _f(v, nd=2, scale=1.0):
    return "-" if v is None or (isinstance(v, float) and not math.isfinite(v)) else f"{v * scale:,.{nd}f}"


def _write_html(path, rows, runs, env, maxd, count, dem, grid, opt, version, seconds):
    from . import report as RP
    ok = [r for r in rows if r["status"] == "ok"]
    dates = [r["date"] for r in ok]
    chart_a = _line_svg(dates, [("area", [r["area_km2"] for r in ok], [r.get("area_p10_km2") for r in ok],
                                 [r.get("area_p90_km2") for r in ok], "#2171b5")], "km²", "Flooded area")
    chart_v = _line_svg(dates, [("volume", [(r["volume_m3"] or 0) / 1e6 for r in ok],
                                 [(r["volume_p10_m3"] or 0) / 1e6 if r.get("volume_p10_m3") is not None else None
                                  for r in ok],
                                 [(r["volume_p90_m3"] or 0) / 1e6 if r.get("volume_p90_m3") is not None else None
                                  for r in ok], "#08519c")], "Mm³", "Flood volume")
    chart_d = _line_svg(dates, [("interior", [r["interior_mean_m"] for r in ok], None, None, "#6baed6"),
                                ("p90", [r["p90_depth_m"] for r in ok], None, None, "#08306b")], "m",
                        "Mean depth, interior (light) and 90th percentile (dark)")
    h = ["<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<meta name='viewport' content='width=device-width,initial-scale=1'>",
         "<title>FloodDepth time series</title><style>", RP.CSS, "</style></head><body>",
         "<header><div class='wrap'><h1>Flood depth over time</h1>",
         f"<p class='meta'><b>{len(ok)}</b> of {len(rows)} image date(s) mapped · {html.escape(dates[0])} to "
         f"{html.escape(dates[-1])} · {RP._duration(seconds)} · FloodDepth {html.escape(version)}</p>",
         "<div class='kpis'>",
         f"<div class='kpi'><b>{env['area_km2']:,.2f} km²</b><span>Flooded on any date</span></div>",
         f"<div class='kpi'><b>{env['flooded_every_date_km2']:,.2f} km²</b><span>Flooded on every mapped date</span></div>",
         f"<div class='kpi'><b>{_f(env['mean_depth_m'])} m</b><span>Mean of the maximum depth</span></div>",
         f"<div class='kpi'><b>{env['volume_m3'] / 1e6:,.2f} Mm³</b><span>Volume of the maximum depth</span></div>",
         "</div></div></header><main class='wrap'>"]
    warn = [f"{r['date']}: {r['status']}" for r in rows if r["status"] != "ok"]
    kinds = sorted({f"{r.get('sensor') or '?'} / {r.get('index') or '?'}" for r in ok})
    if len(kinds) > 1:
        warn.append("The dates use different sensors or water indices (" + ", ".join(kinds) + "): part of the "
                    "change in area is the sensor, not the flood. Sentinel-1 SAR misses shallow water among crops "
                    "and vegetation and maps less than Sentinel-2 of the same flood; a rendered 'False color' "
                    "image allows NDWI only.")
    if warn:
        h.append("<section class='card'><h2>Check these</h2>" + "".join(
            f"<div class='warn'>{html.escape(w)}</div>" for w in warn) + "</section>")
    h.append("<div class='cols'><section class='card'>" + chart_a + chart_v + "</section>"
             "<section class='card'>" + chart_d
             + "<p class='sub'>Bars: P10-P90 of the uncertainty ensemble for that date. Dates are only comparable "
               "when they cover the same ground: give an area of interest, and watch the cloud share.</p></section></div>")
    try:
        m = RP.maps(maxd, np.isfinite(maxd), dem, grid, valid=np.isfinite(dem))
        h.append("<section class='card mapcard'><h2>Maximum observed depth</h2><div class='maps'><figure>"
                 + m["depth"] + RP._gradient_legend(RP.DEPTH_RAMP, "m", "Maximum depth over the dates (m)", "m")
                 + "<p class='sub'>The deepest water seen at each cell on any mapped date (max_depth.tif). "
                   "dates_flooded.tif counts the dates each cell was flooded; first_date_flooded.tif gives the "
                   "first one.</p></figure></div></section>")
    except Exception as exc:  # the summary must not fail the run
        h.append(f"<section class='card'><p class='sub'>Map not drawn: {html.escape(str(exc))}</p></section>")
    h.append("<section class='card pg'><h2>By date</h2><table><tr><th>Date</th><td><b>Sensor / index</b></td>"
             "<td><b>Area km²</b></td>"
             "<td><b>P10-P90 km²</b></td><td><b>Mean / interior m</b></td><td><b>P90 m</b></td>"
             "<td><b>Volume Mm³</b></td><td><b>Cloud</b></td><td><b>Report</b></td></tr>")
    for r in rows:
        if r["status"] != "ok":
            h.append(f"<tr><th>{html.escape(r['date'])}</th><td colspan='8'>{html.escape(r['status'])}</td></tr>")
            continue
        h.append(f"<tr><th>{html.escape(r['date'])}</th>"
                 f"<td>{html.escape(str(r.get('sensor') or '-'))} / {html.escape(str(r.get('index') or '-'))}</td>"
                 f"<td>{_f(r['area_km2'])}</td>"
                 f"<td>{_f(r.get('area_p10_km2'), 1)} - {_f(r.get('area_p90_km2'), 1)}</td>"
                 f"<td>{_f(r['mean_depth_m'])} / {_f(r['interior_mean_m'])}</td><td>{_f(r['p90_depth_m'])}</td>"
                 f"<td>{_f(r['volume_m3'], 2, 1e-6)}</td><td>{_f(r.get('cloud_share'), 1, 100)} %</td>"
                 f"<td><a href='{html.escape(r['report'].replace(os.sep, '/'))}'>report</a></td></tr>")
    h.append("</table><p class='sub'>Also in series.csv. Each date's folder holds the usual outputs (depth.tif, "
             "evidence.tif, flood_extent.gpkg, report.html ...).</p></section>")
    h.append("<section class='card keep'><details><summary><b>Limitations</b> - read before comparing dates</summary><ul>"
             "<li>Each date is an independent estimate from that day's image: differences smaller than the P10-P90 "
             "bands are not meaningful.</li><li>Cloud, image footprint and permanent-water masks differ between "
             "dates; an area of interest keeps the compared ground the same.</li><li>The maximum depth is the deepest "
             "water seen on the mapped dates, not the flood peak, which may fall between overpasses.</li>"
             + "".join(f"<li>{html.escape(x)}</li>" for x in RP.LIMITS[:3]) + "</ul></details></section>")
    h.append("</main></body></html>")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(h))
