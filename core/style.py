# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Report look shared by the run report and the time-series page.

The palette, type and print layout follow the reports of the LandXML Road
& Terrain GIS Tools (same author): a title band with the report facts,
numbered sections, indigo / teal / amber cards, tabular figures and an A4
print layout with page numbers. No QGIS imports.
"""

import html
import re

CSS = """
:root{--fg:#1f2433;--muted:#5d6377;--line:#dcdcec;--bg:#fff;--soft:#f6f6fc;--navy:#2b2766;--indigo:#4a45a6;
--indigo-soft:#ecebf8;--teal:#1f8f80;--teal-soft:#e6f4f1;--amber:#e39a26;--amber-soft:#fdf1dc;--warn:#7a4b06}
*{-webkit-print-color-adjust:exact;print-color-adjust:exact;box-sizing:border-box}
html{background:#f1f1f7}
body{font:14px/1.5 "Segoe UI",Inter,system-ui,-apple-system,Roboto,"Liberation Sans",sans-serif;color:var(--fg);
background:var(--bg);margin:24px auto;max-width:1180px;padding:0 0 28px;box-shadow:0 2px 18px #2b276614;border-radius:8px}
.band{background:linear-gradient(90deg,var(--indigo-soft),#f8f8fd 70%);border-left:8px solid var(--indigo);
border-bottom:3px solid var(--amber);padding:20px 32px 16px;border-radius:8px 8px 0 0}
.band .kicker{font-size:11px;letter-spacing:.16em;text-transform:uppercase;color:#b06e0c;font-weight:600;margin-bottom:6px}
.band h1{color:var(--navy);font-size:26px;margin:0 0 6px;letter-spacing:-.01em}
.band .facts{display:flex;flex-wrap:wrap;gap:4px 22px;font-size:12.5px;color:var(--muted)}
.band .facts b{color:var(--indigo);font-weight:600;margin-right:4px}
.wrap{padding:0 32px}
main.wrap{padding-top:8px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:10px;margin:16px 0 4px}
.kpi{border:1px solid var(--line);border-left:4px solid var(--indigo);border-radius:6px;padding:9px 12px;background:#fff}
.kpi:nth-child(3n+2){border-left-color:var(--teal)}.kpi:nth-child(3n+3){border-left-color:var(--amber)}
.kpi span{display:block;color:var(--muted);font-size:11.5px;text-transform:uppercase;letter-spacing:.04em;order:-1}
.kpi{display:flex;flex-direction:column}
.kpi b{font-size:19px;font-weight:650;font-variant-numeric:tabular-nums;color:var(--navy);line-height:1.3}
.kpi small{color:var(--muted);font-size:11.5px}
.card{background:#fff;min-width:0;margin-top:22px}
.cols{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,460px),1fr));gap:8px 28px;align-items:start}
.cols>.card{margin-top:22px}
h2{font-size:18px;color:var(--navy);margin:0 0 10px;padding-bottom:5px;border-bottom:2px solid var(--indigo)}
h2 .no{display:inline-block;min-width:26px;padding:0 6px;margin-right:10px;border-radius:4px;background:var(--teal);
color:#fff;font-size:14px;text-align:center;font-variant-numeric:tabular-nums}
h3{font-size:15px;color:var(--indigo);margin:18px 0 6px}
.sub,.meta{color:var(--muted);font-size:12.5px;margin:8px 0 0}
table{border-collapse:collapse;width:100%;font-size:13px;font-variant-numeric:tabular-nums;margin:6px 0 10px}
th,td{border-bottom:1px solid var(--line);padding:5px 9px;vertical-align:top;text-align:left}
th{color:var(--muted);font-weight:600;width:42%}
tr:first-child td b{color:var(--navy)}
tbody tr:nth-child(even) td,tbody tr:nth-child(even) th{background:#fafafe}
td{overflow-wrap:anywhere}
table.grid th{width:auto}table.grid td{white-space:nowrap}
table.grid tr:first-child th{background:var(--indigo-soft)}
thead th,table.grid tr:first-child th,table.grid tr:first-child td{background:var(--indigo-soft);color:var(--navy);
font-weight:650;border-bottom:2px solid var(--indigo)}
.warn,.note{padding:8px 12px;margin:8px 0;font-size:13px;border-radius:0 4px 4px 0}
.warn{background:var(--amber-soft);border-left:3px solid var(--amber);color:var(--warn)}
.note{background:var(--teal-soft);border-left:3px solid var(--teal);color:#2c3a3c}
details summary{cursor:pointer;color:var(--indigo);font-size:13px}
pre{white-space:pre-wrap;font-size:12.5px;margin:6px 0 0}
.maps{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,460px),1fr));gap:18px}
figure{margin:0}
svg.map{width:100%;height:auto;display:block;border-radius:6px;border:1px solid var(--line);background:#f6f6f4}
svg.chart{width:100%;height:auto;border:1px solid var(--line);border-radius:6px;background:#fff}
figcaption{font-size:12.5px;color:var(--muted);margin-top:8px}
.ax{font-size:12px;fill:var(--fg)}
.lg{margin-top:10px;max-width:360px}.lgt{font-size:12px;color:var(--muted);margin-bottom:3px}
.bar{height:10px;border-radius:3px}
.ticks{position:relative;height:16px;font-size:11.5px;color:var(--muted)}
.ticks span{position:absolute;transform:translateX(-50%);top:1px}.ticks span:first-child{transform:none}
.ticks span.u{display:none}
.sw{display:flex;gap:14px;flex-wrap:wrap;font-size:12px;color:var(--muted);margin-top:8px}
.sw i{display:inline-block;width:14px;height:10px;border-radius:2px;margin-right:5px;vertical-align:-1px}
.path{color:var(--muted);font-size:12.5px;overflow-wrap:anywhere}
ul.files{margin:6px 0 0;padding-left:18px;font-size:13px}
ol,ul{padding-left:20px;margin:0}li{margin:3px 0}
code{font-size:13px}
a{color:var(--indigo)}
.toc{border:1px solid var(--line);border-radius:6px;padding:10px 16px;margin:18px 0 6px;background:var(--soft)}
.toc b{display:block;font-size:11px;letter-spacing:.12em;text-transform:uppercase;color:var(--indigo);margin-bottom:4px}
.toc ol{margin:0;padding-left:20px;columns:2;column-gap:28px;font-size:13px}.toc a{color:var(--fg);text-decoration:none}
.print-hint{font-size:12px;color:var(--muted);border:1px dashed var(--line);border-radius:6px;padding:6px 10px;margin:14px 0 0}
.foot{margin:34px 32px 0;padding-top:10px;border-top:1px solid var(--line);font-size:11.5px;color:var(--muted)}
.printonly{display:none}
@media print{
  html{background:#fff}
  body{font-size:9.6pt;margin:0;max-width:none;padding:0;box-shadow:none;border-radius:0}
  .band{border-radius:0;padding:10px 14px 8px}.band h1{font-size:17pt}.wrap{padding:0 2px}.foot{margin:18px 2px 0}
  h2{font-size:12.5pt}h3{font-size:10.5pt;margin:10px 0 4px}h2,h3{break-after:avoid-page}
  .print-hint,.toc,.noprint{display:none}
  .kpis{grid-template-columns:repeat(3,1fr);gap:6px}.kpi{padding:5px 8px}.kpi b{font-size:12pt}
  .kpi span,.kpi small{font-size:7.5pt}
  .card,.cols>.card{margin-top:4mm}.cols{grid-template-columns:1fr 1fr;gap:2mm 5mm}
  figure.prob,.pg{break-before:page}
  figure.prob svg.map{max-height:105mm}
  .cols>.card,.card.keep,.warn,.note,figure,tr,.kpi{break-inside:avoid}
  .maps{grid-template-columns:1fr;gap:4mm}
  section.mapcard{break-inside:auto}
  svg.map{max-height:118mm;width:auto;max-width:100%;margin:0 auto}
  table{font-size:8.3pt}th,td{padding:.6mm 1.5mm}
  .sub,figcaption,.path,ul.files,.sw,.lgt{font-size:7.5pt}
  .ticks{font-size:7pt;height:11px}
  details>summary{list-style:none;color:var(--fg)}details>summary::-webkit-details-marker{display:none}
  .printonly{display:block!important}
  .warn,.note{font-size:8.3pt;padding:1.5mm 2.5mm}
}
"""


def page_rule(title):
    label = title.replace("\\", "").replace('"', "'")
    return ("@page{size:A4 portrait;margin:14mm 12mm 16mm 12mm;"
            f'@bottom-left{{content:"{label}";font:8pt system-ui,sans-serif;color:#5b6575}}'
            '@bottom-right{content:"Page " counter(page) " of " counter(pages);font:8pt system-ui,sans-serif;'
            'color:#5b6575}}')


def head(title):
    return ("<!doctype html><html lang='en'><head><meta charset='utf-8'>"
            "<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{html.escape(title)}</title><style>{CSS}{page_rule(title)}</style></head><body>")


def band(title, facts, kicker="FloodDepth - flood extent, depth and uncertainty"):
    """The title band: kicker, title and [(label, value)] facts (plain text)."""
    return (f"<header class='band'><div class='kicker'>{html.escape(kicker)}</div>"
            f"<h1>{html.escape(title)}</h1><div class='facts'>"
            + "".join(f"<span><b>{html.escape(str(k))}</b>{html.escape(str(v))}</span>" for k, v in facts if v)
            + "</div></header>")


PRINT_HINT = ("<div class='print-hint'>To save as PDF, print this page from a browser (Chrome or Edge: Save as PDF, "
              "A4). Pages are numbered and maps keep their colours; untick <i>Headers and footers</i> to hide the "
              "browser's file path.</div>")


def foot(version, when, extra=""):
    return (f"<footer class='foot'>Prepared with FloodDepth {html.escape(str(version))} on {html.escape(when)}. "
            "Terrain-based estimates from the files named in the report, not a hydraulic model; read the warnings "
            f"and limitations before use.{(' ' + extra) if extra else ''}</footer>")


_H2 = re.compile(r"<h2(?P<attrs>[^>]*)>(?P<text>.*?)</h2>", re.S)


def number_sections(body, toc_min=4):
    """Number the h2 headings, give them anchors and put a contents box
    before the first one when there are at least `toc_min`."""
    found = []

    def repl(m):
        text = m.group("text")
        n = len(found) + 1
        found.append((f"s{n}", re.sub(r"<[^>]+>", "", text).strip()))
        return f"<h2{m.group('attrs')} id='s{n}'><span class='no'>{n}</span>{text}</h2>"

    body = _H2.sub(repl, body)
    if len(found) >= toc_min:
        toc = ("<nav class='toc'><b>Contents</b><ol>"
               + "".join(f"<li><a href='#{a}'>{html.escape(t)}</a></li>" for a, t in found) + "</ol></nav>")
        i = body.find("<section")
        body = body[:i] + toc + body[i:] if i >= 0 else toc + body
    return body
