# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Identify which file and band holds green, NIR, SWIR, VV ... without asking.

Satellite downloads come in three shapes, and FloodDepth accepts all of them
from the same "image file(s)" input:

1. One file per band - Copernicus Browser / EO Browser "Raw" exports
   (``..._Sentinel-2_L2A_B03_(Raw).tiff``), Sentinel-2 SAFE / HLS
   (``..._B03_10m.jp2``, ``....B03.tif``), Landsat Collection 2
   (``..._SR_B3.TIF``). The band is read from the file name.
2. One multiband stack whose bands carry descriptions or wavelengths
   (GDAL's SENTINEL2 driver, most EO exports): read from band metadata.
3. A rendered composite: the Browser "False color" product is
   2.5 x (B08, B04, B03) on all three channels, so NDWI is recoverable from
   channels 3 and 1 (same gain; only saturated pixels are affected).
   Sentinel-1 products are recognised by name or by VV/VH bands.

Anything that cannot be identified with certainty is reported as such and
the user is asked for band numbers; FloodDepth never guesses band order
from the band count alone, except for a 12- or 13-band file explicitly
named as Sentinel-2.
"""

import os
import re

ROLES = ("blue", "green", "red", "nir", "swir1", "swir2", "vv", "vh")

S2 = {"B02": "blue", "B03": "green", "B04": "red", "B08": "nir", "B8A": "nir_narrow",
      "B11": "swir1", "B12": "swir2"}
L89 = {"B2": "blue", "B3": "green", "B4": "red", "B5": "nir", "B6": "swir1", "B7": "swir2"}
L457 = {"B1": "blue", "B2": "green", "B3": "red", "B4": "nir", "B5": "swir1", "B7": "swir2"}
S2_L2A_ORDER = ["B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B09", "B11", "B12"]
S2_L1C_ORDER = ["B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B09", "B10", "B11", "B12"]

WAVELENGTH_NM = (("blue", 440, 520), ("green", 520, 600), ("red", 620, 700),
                 ("nir", 760, 900), ("swir1", 1520, 1760), ("swir2", 2050, 2360))

DESCRIPTION_WORDS = (("swir1", r"\bswir[\s_-]*1\b|\bswir16\b|\bmir\b"),
                     ("swir2", r"\bswir[\s_-]*2\b|\bswir22\b"),
                     ("nir", r"\bnir\b|near[\s_-]*infra|\bnir08\b"),
                     ("green", r"\bgreen\b"), ("red", r"\bred\b"), ("blue", r"\bblue\b"),
                     ("vv", r"\bvv\b"), ("vh", r"\bvh\b"))

_BAND_TOKEN = re.compile(r"(?:^|[_\-\s.(])(?:SR_|band)?B(8A|\d{1,2})(?=[_\-\s.)]|$)", re.I)
_DATE = re.compile(r"(20\d{2})[-_]?([01]\d)[-_]?([0-3]\d)")


def sensor_family(name):
    n = os.path.basename(name).upper()
    if re.search(r"\bS1[ABC]_|SENTINEL[-_ ]?1|_GRD|\bIW\b|[_-]VV[_.-]|[_-]VH[_.-]", n):
        return "s1"
    if re.search(r"LC0?[89]|LANDSAT[-_ ]?[89]|_OLI", n):
        return "landsat89"
    if re.search(r"L[TE]0?[457]|LANDSAT[-_ ]?[457]", n):
        return "landsat457"
    if re.search(r"SENTINEL[-_ ]?2|\bS2[ABC]\b|S2[ABC]_|MSIL[12][AC]|HLS\.S30|_L2A_|_L1C_", n):
        return "s2"
    return None


def band_token(name):
    """'B03', 'B8A', 'B3' ... found in a file name, else None."""
    m = list(_BAND_TOKEN.finditer(os.path.basename(name)))
    if not m:
        return None
    t = m[-1].group(1).upper()
    return "B" + t


def acquisition_date(name):
    m = _DATE.search(os.path.basename(name))
    if not m:
        return None
    y, mo, d = m.groups()
    return f"{y}-{mo}-{d}"


def composite_kind(name):
    n = os.path.basename(name).lower()
    if "false_color" in n or "false color" in n or "false-colour" in n or "false_colour" in n:
        if "urban" in n:
            return "s2_false_color_urban"
        return "s2_false_color"
    if "true_color" in n or "natural_color" in n or "true color" in n:
        return "rgb"
    if "scene_classification" in n:
        return "scl_visual"
    if re.search(r"_ndwi|_ndvi|_ndsi|moisture_index|_swir\.|custom_script", n):
        return "rendered_index"
    return None


def _role_from_token(token, family):
    if token is None:
        return None
    t = token.upper()
    if family == "landsat89":
        return L89.get("B" + str(int(t[1:])) if t[1:].isdigit() else t)
    if family == "landsat457":
        return L457.get("B" + str(int(t[1:])) if t[1:].isdigit() else t)
    if t[1:].isdigit():
        t = "B%02d" % int(t[1:])
    return S2.get(t)


def _role_from_band(b):
    md = {k.upper(): str(v) for k, v in (b.get("metadata") or {}).items()}
    for key in ("WAVELENGTH", "CENTRAL_WAVELENGTH", "CENTRAL_WAVELENGTH_NM"):
        if key in md:
            try:
                wl = float(re.findall(r"[\d.]+", md[key])[0])
                unit = md.get("WAVELENGTH_UNITS", "").lower()
                if wl < 30 or unit.startswith("micro") or unit == "um":
                    wl *= 1000.0
                for role, lo, hi in WAVELENGTH_NM:
                    if lo <= wl <= hi:
                        return role
            except (IndexError, ValueError):
                pass
    for key in ("BANDNAME", "BAND_NAME"):
        if key in md:
            r = _role_from_token(md[key] if md[key].upper().startswith("B") else "B" + md[key], "s2")
            if r:
                return r
    desc = (b.get("description") or "").lower()
    tok = band_token(desc)
    if tok and re.fullmatch(r"b\d{1,2}|b8a", desc.strip()):
        return _role_from_token(tok, "s2")
    for role, pat in DESCRIPTION_WORDS:
        if re.search(pat, desc):
            return role
    return None


class Scene(object):
    """What was found: role -> (path, band); product kind; notes."""

    def __init__(self):
        self.sources = {}          # role -> (path, band index)
        self.kind = None           # 'bands', 's2_false_color', 'sar', 'rgb', 'unknown'
        self.family = None
        self.dtype = {}
        self.offset = {}           # role -> additive DN offset applied before scaling
        self.notes = []
        self.date = None
        self.rendered_sar = None   # name of a Sentinel-1 visualisation that was refused
        self.offset_candidate = False

    def has(self, *roles):
        return all(r in self.sources for r in roles)

    def describe(self):
        parts = [f"{r}: {os.path.basename(p)} band {b}" for r, (p, b) in sorted(self.sources.items())]
        return "; ".join(parts) if parts else "no bands identified"


def identify(paths, band_info_fn, overrides=None):
    """Build a Scene from one or more rasters.

    `band_info_fn(path)` returns raster.band_info(path) (injected so the
    logic can be tested without files). `overrides` maps role -> 1-based
    band number in the FIRST file and always wins.
    """
    sc = Scene()
    infos = [band_info_fn(p) for p in paths]
    if not infos:
        sc.kind = "unknown"
        return sc
    sc.family = next((sensor_family(p) for p in paths if sensor_family(p)), None)
    sc.date = next((acquisition_date(p) for p in paths if acquisition_date(p)), None)
    overrides = {k: int(v) for k, v in (overrides or {}).items() if v}

    first = infos[0]
    for role, band in overrides.items():
        if band > first["count"]:
            raise ValueError(f"{os.path.basename(first['path'])} has {first['count']} band(s); "
                             f"the {role} band number {band} does not exist.")
        sc.sources[role] = (first["path"], band)
    if overrides:
        sc.notes.append("Band numbers set by the user: " + ", ".join(f"{k}={v}" for k, v in overrides.items()))

    # 1. one file per band, named
    single = [i for i in infos if (i["count"] - (1 if i["alpha"] else 0)) == 1]
    if single:
        fam = sc.family
        if fam is None and any(band_token(i["path"]) in ("B8A", "B11", "B12", "B03", "B08")
                               for i in single):
            fam = "s2"
        for i in single:
            tok = band_token(i["path"])
            role = _role_from_token(tok, fam) if fam in ("s2", "landsat89", "landsat457", None) else None
            if role is None and fam == "s1":
                n = os.path.basename(i["path"]).upper()
                role = "vh" if re.search(r"[_-]VH[_.-]", n) else ("vv" if re.search(r"[_-]VV[_.-]", n) else None)
            if role is None:
                role = _role_from_band(i["bands"][0])
            if role == "nir_narrow":
                role = "nir" if "nir" not in sc.sources else None
            if role and role not in sc.sources:
                sc.sources[role] = (i["path"], 1)
        if fam:
            sc.family = fam

    # 2. multiband stacks with metadata, composites, or named S2 stacks
    for i in infos:
        n_data = i["count"] - (1 if i["alpha"] else 0)
        if n_data < 2:
            continue
        kind = composite_kind(i["path"])
        if kind == "s2_false_color":
            for role, b in (("nir", 1), ("red", 2), ("green", 3)):
                sc.sources.setdefault(role, (i["path"], b))
            sc.kind = sc.kind or "s2_false_color"
            sc.notes.append("Sentinel-2 'False color' composite recognised (channels = 2.5 x B08, B04, B03); "
                            "NDWI uses channels 3 and 1. Saturated pixels are less reliable.")
            continue
        if sensor_family(i["path"]) == "s1" or kind is None and any(
                re.search(r"\bv[vh]\b", (b["description"] or "").lower()) for b in i["bands"]):
            found = False
            for b in i["bands"]:
                r = _role_from_band(b)
                if r in ("vv", "vh"):
                    sc.sources.setdefault(r, (i["path"], b["index"]))
                    found = True
            if not found:
                sc.rendered_sar = os.path.basename(i["path"])
            sc.family = "s1"
            continue
        if kind in ("rgb", "s2_false_color_urban", "scl_visual", "rendered_index"):
            sc.notes.append(f"{os.path.basename(i['path'])}: rendered '{kind}' product - no usable "
                            "water-index bands. Download the 'Raw' bands (B03, B08, B11) instead.")
            if kind == "rgb":
                sc.kind = sc.kind or "rgb"
            continue
        named = 0
        for b in i["bands"]:
            if b["alpha"]:
                continue
            r = _role_from_band(b)
            if r == "nir_narrow":
                r = "nir"
            if r and r not in sc.sources:
                sc.sources[r] = (i["path"], b["index"])
                named += 1
        if named == 0 and sensor_family(i["path"]) == "s2" and n_data in (12, 13):
            order = S2_L2A_ORDER if n_data == 12 else S2_L1C_ORDER
            for idx, code in enumerate(order, 1):
                r = S2.get(code)
                if r == "nir_narrow":
                    continue
                if r and r not in sc.sources:
                    sc.sources[r] = (i["path"], idx)
            sc.notes.append(f"{n_data}-band file named as Sentinel-2: band order assumed "
                            f"{'L2A (no B10)' if n_data == 12 else 'L1C'}. Check the band list.")

    if sc.kind is None:
        if sc.has("vv") or sc.has("vh"):
            sc.kind = "sar"
        elif sc.has("green") and (sc.has("nir") or sc.has("swir1")):
            sc.kind = "bands"
        else:
            sc.kind = "unknown"
    if sc.kind == "sar" and sc.family is None:
        sc.family = "s1"

    # 3. Sentinel-2 L2A DN offset (processing baseline 04.00, from 25 Jan 2022)
    for role, (p, b) in sc.sources.items():
        info = next(x for x in infos if x["path"] == p)
        dt = info["bands"][b - 1]["dtype"] if b >= 1 else info["bands"][0]["dtype"]
        sc.dtype[role] = dt
        off = 0.0
        md = {k.upper(): v for k, v in info["metadata"].items()}
        keys = [k for k in md if k.startswith("BOA_ADD_OFFSET")]
        if keys:
            try:
                off = float(md[keys[0]])
            except ValueError:
                off = 0.0
        elif (sc.family == "s2" and dt in ("UInt16", "Int16") and sc.date and sc.date >= "2022-01-25"
              and re.search(r"L2A", os.path.basename(p), re.I)):
            sc.offset_candidate = True   # decided from the data (imagery.decide_offset)
        sc.offset[role] = off
    if any(v for v in sc.offset.values()):
        sc.notes.append("Sentinel-2 L2A BOA_ADD_OFFSET from the file metadata is applied before computing indices.")
    return sc


def reflectance_scale(dtype, sample_p99):
    """Factor that turns stored values into approximate reflectance (0-1).

    Byte    : Browser exports store reflectance x 255 (approximate).
    UInt16  : 0-10000 digital numbers.
    Float   : already reflectance if the 99th percentile is below 2.
    Indices are scale-free; only the optional bright-cloud test uses this.
    """
    if dtype == "Byte":
        return 1.0 / 255.0
    if dtype in ("UInt16", "Int16", "UInt32", "Int32"):
        if sample_p99 is not None and sample_p99 > 12000.0:
            return 1.0 / 65535.0
        return 1.0 / 10000.0
    if sample_p99 is not None and sample_p99 <= 2.0:
        return 1.0
    if sample_p99 is not None and sample_p99 > 12000.0:
        return 1.0 / 65535.0      # Browser uint16 exports: reflectance x 65535
    return 1.0 / 10000.0
