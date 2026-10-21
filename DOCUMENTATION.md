# FloodDepth 0.4.1 — technical documentation

This document describes what FloodDepth computes, every decision it takes
automatically, the uncertainty model, how it was validated, and where it
should not be trusted. The README is the short version.

Contents

1. Workflow
2. Inputs and band recognition
3. Water classification
4. Permanent water, cloud and connectivity
5. Flood-trace mode (receded floods, erosion and deposition)
6. Water surface and depth
7. Uncertainty
8. Outputs
9. Validation
10. Limitations and good practice
11. Implementation notes
12. References

---

## 1. Workflow

```
image file(s) ──► band recognition ──► water index ──► threshold ──► cleanup ─┐
pre-event image / WorldCover / JRC / polygons ──► permanent water ────────────┤
                                                                             ▼
DEM ──► analysis grid (projected) ◄── flood extent (from the image or polygons)
                     │
                     ├─► waterline cells ─► outlier screen ─► water surface ─► depth
                     ├─► evidence classes (observed water / trace / erosion / steep)
                     └─► Monte Carlo ensemble ─► depth P10 / P50 / P90, flood probability
                                                                             ▼
                                                     report.html, run.json, rasters, GeoPackage
```

*Quick flood depth* runs all of it. The step-by-step tools expose each part for
the Graphical Modeler: (1) classify, (1b) supervised classification, (2) water
surface and depth from an extent, (3) depth from a known surface or level, and
the extent comparison.

## 2. Inputs and band recognition

**Files as downloaded.** The image input accepts one or more rasters. The core
(`core/sensors.py`) assigns roles (blue, green, red, NIR, SWIR1, SWIR2, VV, VH)
in this order of evidence:

1. Band numbers given by the user (advanced parameters) — always win.
2. One file per band: the band code in the file name — `…_B03_(Raw).tiff`
   (Copernicus / EO Browser), `…_B03_10m.jp2` (SAFE), `….B03.tif` (HLS),
   `…_SR_B3.TIF` (Landsat). The sensor family decides the table: Landsat 8/9
   B3 = green, B5 = NIR, B6 = SWIR1; Landsat 4–7 B2 = green, B4 = NIR, B5 =
   SWIR1; Sentinel-2 B03, B08 (B8A if B08 is absent), B11.
3. Band metadata: wavelength (nm or µm), `BANDNAME`, or descriptions such as
   "Green", "NIR", "SWIR1", "B3", "VV".
4. A 12- or 13-band file *named* as Sentinel-2 (L2A without B10 / L1C).
5. The Browser **False color** composite: channels are 2.5 × (B08, B04, B03),
   so NDWI is computed from channels 3 and 1 (same gain; saturated pixels are
   less reliable).

FloodDepth never guesses a band order from the band count alone. True-colour
images and rendered Sentinel-1 visualisations ("Enhanced visualization", "RGB
ratio") are refused with an explanation: their colours are display choices,
not reflectance or backscatter, and on the Nyando 6 May 2024 scene a
brightness rule classified 67 % of the image (including the black no-data
margin) as water.

**NoData.** Every reader masks: the declared NoData value, non-finite values,
undeclared float sentinels (|v| > 1e30, e.g. −3.4028e38 in GLO-30 HAND clips),
an alpha / data-mask band (declared, or the last band of a 2- or 4-band integer
raster holding only 0 and 255 — the Browser "Raw" layout), and for DEMs, exact
zeros when no NoData is declared, they cover > 0.5 % of the grid and every
other value is above 50 (a clip with 0 fill).

**Sentinel-2 offset.** Since processing baseline 04.00 (25 January 2022) L2A
digital numbers are reflectance × 10 000 + 1 000. The −1 000 offset is applied
when the file metadata declares `BOA_ADD_OFFSET`, or when the data show it:
every water-index band has its 1st percentile ≥ 1 000 and its 99th ≤ 12 000.
Browser exports are harmonised (float or byte) or rescaled (uint16 =
reflectance × 65 535, as in the Buulo Gaduud–Afmadow 3 Dec 2023 scene) and are
left alone; an offset applied to them would bias MNDWI.

**Analysis grid.** Imagery is classified on its own grid. Depth is computed on
one projected grid: a projected DEM's own lattice (cells are never shifted),
cut to the image or extent; a geographic DEM (Copernicus, HAND) on the WGS 84 /
UTM zone of the area at the DEM's resolution in metres. Masks move between
grids by area fraction (≥ 0.5 = water by default).

## 3. Water classification

**Index.** NDWI = (green − NIR) / (green + NIR) (McFeeters 1996) and MNDWI =
(green − SWIR1) / (green + SWIR1) (Xu 2006). With green, NIR and SWIR1 the
automatic choice is **NDWI + confident MNDWI**: water where NDWI is above its
threshold, plus where MNDWI is above max(its threshold, 0.3). MNDWI catches
turbid and shallow flood water that NDWI misses, but alone it over-maps wet
soil and bright built-up surfaces (area bias 1.49 on Sen1Floods11); only its
confident part is added (area bias 1.01, § 9.2). With green and NIR only:
NDWI; with green and SWIR1 only: MNDWI. Sentinel-1 VV (else VH) in dB, water =
dark, after a 5 × 5 boxcar mean in linear power (speckle).

**Threshold — edge-based Otsu (default).** Following Donchyts et al. (2016):
classify with the physical threshold (0 for NDWI / MNDWI, −16 dB for VV,
−22 dB for VH), take the pixels within two cells of the resulting water / land
edges — a sample in which both classes are present in similar numbers — and
apply Otsu (1979) to that sample only. When one class is absent or the edge is
shorter than 200 pixels, the physical threshold is kept and the report says
so. Tile-based Otsu (Martinis et al. 2009; mixed tiles pooled) and global Otsu
remain as options. The choice of default is empirical (§ 9.2): on the 446
Sen1Floods11 chips edge-based Otsu scored the highest CSI; tile-based Otsu
split land classes on scenes with little water.

**Separability.** Otsu's η = σ²between / σ²total (0–1) is reported. A single
Gaussian population gives η = 2/π ≈ 0.64 (verified in `tests/test_classify.py`),
so FloodDepth calls η ≥ 0.75 *good*, 0.65–0.75 *moderate* and warns below 0.65.
Otsu's maximum can be flat when bins between the classes are empty (8-bit
data); the middle of the flat maximum is used, not its first edge.

**Threshold band.** The thresholds whose between-class variance is ≥ 95 % of
the maximum form the *near-optimal band*: the range over which the split is
almost equally good. It is reported and drives the extent part of the
uncertainty.

**Plausibility guard.** An automatic threshold outside the range where open
water separates — NDWI / MNDWI −0.3 to 0.6 (reflectance data), VV −26 to
−10 dB, VH −32 to −16 dB — is refused, as is a scene where less than 0.05 % of
pixels fall on the water side. The message explains the likely cause (no water
at that date; cloud) and points to flood-trace mode. The 12 May 2024 Mai Mahiu
image is the reference case: the 0.3 Otsu reported 7.5 % "water" in a scene
with none.

**Cleanup.** Water patches smaller than *N* image cells (default 10) are
removed. Dry holes are filled only when asked (default 0 — islands, roads and
buildings are kept). Holes touching NoData or the grid edge are never filled.
Both counts are reported. Connected components are labelled with
`gdal.Polygonize` + `gdal.RasterizeLayer` (C++), so no SciPy is needed.

**Supervised (tool 1b).** For RGB, drone and aerial images: a standardised
two-class nearest-centroid model from water and non-water training polygons.

**Mask values.** 1 flood water, 2 permanent water, 0 dry, 255 NoData.

## 4. Permanent water, cloud and connectivity

* **Pre-event image** (any supported product): classified the same way; its
  water is permanent water.
* **Permanent-water raster**: *auto* recognises ESA WorldCover (all values are
  WorldCover codes → class 80) and JRC Global Surface Water occurrence (values
  0–100 with many levels → ≥ 50 %); otherwise any value > 0; or give class
  values (`80,90`) or a threshold (`>=75`).
* **Permanent-water polygons.**

Permanent water is not flood (no depth is reported there), and the edge
between flood and permanent water is **not** a waterline (§ 6).

**Cloud.** Pixels whose blue reflectance exceeds 0.25 (approximate scale from
the data type) are masked and the mask is grown by two cells; cloud edges are
not waterlines. Shadows are not detected — check the image.

**Connectivity.** *Keep only flooding connected to these features* keeps the
flood regions that touch a river, reach or road layer (within a buffer). On
Mai Mahiu it removes crop-change and cloud patches unrelated to the flood path.

## 5. Flood-trace mode (receded floods, erosion and deposition)

Flash floods and debris flows drain within hours, so the optical image after
the event shows the damage, not the water. Flood-trace mode classifies the
NDVI loss between a pre- and a post-event image (dNDVI = NDVI_post − NDVI_pre,
tile-based Otsu, water side = loss). With the SWIR1, red, NIR and blue bands
the trace is split by the bare-soil index BSI = ((SWIR1 + red) − (NIR + blue))
/ ((SWIR1 + red) + (NIR + blue)) (Rikimaru et al. 2002) of the post-event image:

* BSI ≤ 0 — vegetation damaged but present: **inundation**;
* BSI > 0 — bare soil or fresh sediment: **erosion / deposition**.

A trace is not a water extent. Where the ground was scoured the pre-event DEM
no longer describes it, and the trace edge is the erosion width, not the water
width; depths there are indicative and carry extra terrain uncertainty (§ 7).

## 6. Water surface and depth

FloodDepth estimates the water-surface elevation (WSE) from the terrain along
the flood edge — the boundary method of FwDET (Cohen et al. 2018, 2019) — with
these refinements:

1. **Waterline cells** are flood cells with at least one 4-neighbour that is
   valid dry ground. Edges against NoData, the grid edge, cloud, or permanent
   water are not used: the extent stops there for reasons unrelated to the
   water level.
2. **Edge outliers.** Each waterline elevation is compared with the mean of the
   other waterline cells within 10 cells; values more than 3 robust SD (1.4826 ×
   MAD) and more than 0.5 m away are dropped (embankments, buildings, canopy
   in a DSM). Waterline cells on ground steeper than 10 % are never dropped: on
   an escarpment edge elevations change by metres per cell for real. (Without
   this protection the Mai Mahiu trace produced 130 m "depths".)
3. **Edge water level.** The waterline lies somewhere between the last wet
   cell and the first dry one, so the water level at a waterline cell is
   placed halfway up the step from its ground to the lowest dry 4-neighbour,
   lifted by at most 1 m. Taking the wet cell's own ground (depth 0 at the
   edge, as in FwDET) always underestimates; the halfway step cuts the
   flat-floodplain underestimate by 40–80 % (§ 9.3). The cap keeps a cliff,
   bank or embankment from lifting the surface. Set the share to 0 for the
   FwDET convention.
4. **Surface.** The kept edge water levels are interpolated over the flood by
   inverse distance (`gdal.FillNodata`, C++) — the default — or allocated from
   the nearest waterline cell (FwDET v1; exact Euclidean allocation by jump
   flooding in NumPy, schedule bounded with `gdal.ComputeProximity`). The
   default is empirical (§ 9.3): inverse distance was as good or better in
   every benchmark scenario.
5. **Smoothing.** A focal mean (radius 2 cells) over flood cells, applied only
   where it changes the surface by less than 0.5 m, so it never averages across
   an escarpment or gully.
6. **Depth** = WSE − DEM inside the extent. Cells where the ground is above the
   surface get depth 0 and are counted (reported share).

**HAND method.** The same procedure applied to a HAND raster gives the water
height above the drainage at the edge, and depth = edge HAND − HAND. HAND
values inside the extent are heights above the channel, not water depths;
reading them as depths inverts the result (highest "depth" at the flood edge).

**Constant level / known surface.** depth = level − DEM (bathtub) or WSE raster
− DEM, on the DEM grid, optionally limited to an extent.

All steps are whole-array NumPy or GDAL calls over the flood's bounding box.
The 0.3 per-pixel IDW took 19 s for 0.36 Mcells; the Nyando 30 m scene (2
Mcells, 71 km² of flood) now takes seconds for the surface and under a minute
including a 20-member ensemble.

## 7. Uncertainty

### 7.1 Evidence classes (`evidence.tif`)

| Code | Evidence | Meaning |
|---|---|---|
| 1 | observed water | water in the image at the overpass (ponding, standing flood) |
| 2 | trace — inundation | vegetation lost, still vegetated: water passed |
| 3 | trace — erosion / deposition | bare soil or sediment after the event: scour, deposition; DEM outdated |
| 4 | steep flow path | slope > 10 %: water flowed, did not pond |
| 5 | doubtful — ground above the water surface | a mapped flood cell whose ground is above the estimated surface (depth 0): wet soil, crops or vegetation classed as water, a DSM artefact, or a surface too low there. Replaces classes 0–2 (not 3 or 4); pink on the report's depth map |
| 0 | extent from polygons | the user's layer |

### 7.2 Monte Carlo ensemble

Each member (default 20, seed 42 — reproducible) re-runs the surface and the
depth with:

* **DEM error**: a Gaussian random field, SD σ_DEM, Gaussian correlation length
  150 m (FFT filter). σ_DEM is taken from the DEM file name when it names a
  product — LiDAR / drone 0.15 m, FABDEM 1.5 m, Copernicus GLO-30 and TanDEM-X
  2.0 m, ALOS AW3D30 2.5 m, SRTM 4.0 m — otherwise 2.0 m with a warning; it can
  be set.
* **Extra terrain error** (default 1 m, independent) in erosion / deposition
  cells.
* **Water-surface (method) error**: estimated from the data by a spatially
  blocked cross-validation of the waterline. Edge cells are held out in square
  blocks twice the median flood-cell distance to the edge (so a held-out cell
  is predicted from about as far as an interior cell), predicted from the rest,
  and the robust SD of the residuals is applied as a correlated field. It
  measures how far the real surface departs from what the edge can tell, which
  the DEM term does not cover; without it the P10–P90 band covered only 33–65 %
  of the true depths in the DEM benchmark (§ 9.3).
* **Extent**: for image extents the classification threshold is drawn across
  the near-optimal band and the area-fraction cut-off (image → analysis grid)
  in [0.3, 0.7]; for polygons the edge is shifted by −1, 0 or +1 cell (±2 in
  trace mode). Connectivity to the reference features is re-applied.

Members are independent (one random stream each, drawn from the seed up front)
and run in a thread pool — GDAL, NumPy and SciPy release Python's lock in the
heavy calls — so the result does not depend on the number of threads
(CPU count − 1, at most 8; `FLOODDEPTH_THREADS` overrides).

Outputs: `flood_probability.tif` (share of members flooded), `depth_p10.tif`,
`depth_p50.tif`, `depth_p90.tif` (a dry member counts as 0 m) and the P10 / P50
/ P90 of area, mean depth and volume. Windows above 2 Mcells are computed on a
coarser grid (factor reported); the bands do not need full resolution.

The bands cover DEM error, extent and threshold choice, and the surface error
the edge can reveal. They do not cover a systematic bias of the method (§ 9.3)
or the image timing.

## 8. Outputs

| File | Content |
|---|---|
| `depth.tif` | flood depth (m), NoData outside the flood |
| `wse.tif` | water-surface elevation (DEM method) or water height above drainage (HAND) |
| `evidence.tif` | evidence class (§ 7.1), palette |
| `depth_p10/p50/p90.tif`, `flood_probability.tif` | ensemble results |
| `water_mask.tif` | classification on the image grid, palette |
| `flood_extent.gpkg` | `flood_extent` polygons (area_m2), `waterline` points (used / dropped, z, surface) |
| `report.html` | decisions, warnings, statistics, maps, limitations |
| `run.json` | everything in the report, machine-readable, incl. run time by stage (`timings`) |

The report prints on A4 portrait (about four pages).

### 8.1 Several dates (*Flood depth over time*)

The image files of several dates are grouped by the acquisition date in their
names (YYYY-MM-DD, YYYY_MM_DD or YYYYMMDD) and each date runs the quick
workflow in its own sub-folder with the same DEM, permanent water, area of
interest and settings. A date that fails (no water, refused threshold) is
reported and skipped. The per-date grids share the DEM lattice, so they are
combined without resampling on their union:

| File | Content |
|---|---|
| `max_depth.tif` | deepest water seen at each cell on any mapped date (m) |
| `dates_flooded.tif` | number of dates each cell was flooded |
| `first_date_flooded.tif` | first date flooded (1 = first date; dates in the file tags) |
| `series.csv`, `series.html` | area, depth, volume (with P10–P90), sensor / index and cloud per date; charts spaced by time; the maximum-depth map |

Dates are comparable only when they cover the same ground (give an area of
interest) and use the same kind of sensor: Sentinel-1 maps less of a shallow
flood among crops than Sentinel-2, and a rendered 'False color' image allows
NDWI only; the summary warns when the dates mix them. The maximum depth is the
deepest water *seen*, not the peak, which may fall between overpasses.

## 9. Validation

Three kinds of evidence, each reproducible from the repository (`tests/`,
`benchmarks/`; results in `docs/validation/`). No published report was used as
ground truth.

### 9.1 Synthetic site with a known answer (`tests/`)

A 9 km valley on a 30 m DEM (UTM 36S) with a sloping water surface 2 m above
the channel, a permanent lake, a NoData strip, and Sentinel-2 "Raw" bands
written the way the Copernicus Browser exports them (one 8-bit file per band
with an alpha band, EPSG:4326). Through real QGIS Processing
(`tests/qgis_smoke.py`, QGIS 3.34.4): depth bias −0.03 m, MAE 0.03 m;
classification CSI 0.94 (automatic and supervised); the lake gets no depth;
P10 ≤ P50 ≤ P90 everywhere. The same scene drained of water and stripped of
vegetation is refused in water mode and recovered in flood-trace mode (MAE
0.02 m, bare sediment classed as erosion / deposition). 7 pure suites, 120
checks; QGIS smoke test, 32 checks.

### 9.2 Classification — Sen1Floods11 (446 hand-labelled chips)

Bonafilia et al. (2020): 512 × 512 chips at 10 m from 11 flood events
(Bolivia, Ghana, India, Mekong, Nigeria, Pakistan, Paraguay, Somalia, Spain,
Sri Lanka, USA), Sentinel-2 L1C and Sentinel-1 (VV, VH dB), labels 1 water /
0 not water / −1 no data. Every chip goes through `imagery.classify_image` as
in the plugin; the automatic configurations give no band numbers (bands are
recognised from the chips' descriptions). A refusal is scored as "no water
mapped". Scores pooled over all chips (`benchmarks/sen1floods11.py`):

| Configuration | CSI | F1 | Hit rate | False-alarm ratio | Area bias | Refused (with > 5 % water) |
|---|---|---|---|---|---|---|
| **S2 automatic: NDWI + confident MNDWI, edge Otsu** | **0.737** | **0.849** | 0.854 | 0.157 | **1.01** | 90 (3) |
| S2 NDWI, edge Otsu | 0.716 | 0.835 | 0.816 | 0.145 | 0.95 | 94 (4) |
| S2 MNDWI, edge Otsu | 0.613 | 0.760 | 0.948 | 0.366 | 1.49 | 36 (0) |
| S2 MNDWI, tile Otsu | 0.464 | 0.634 | 0.902 | 0.512 | 1.85 | 27 (0) |
| S2 MNDWI, global Otsu | 0.410 | 0.581 | 0.918 | 0.575 | 2.16 | 47 (0) |
| **S1 automatic: VV, speckle 5 × 5, edge Otsu** | **0.447** | 0.618 | 0.509 | 0.215 | 0.65 | 72 (7) |
| S1 VV, tile Otsu | 0.332 | 0.499 | 0.613 | 0.579 | 1.46 | 176 (26) |
| S1 VV, global Otsu | 0.301 | 0.462 | 0.641 | 0.639 | 1.77 | 188 (32) |

By event (automatic configurations, CSI):

| Event | S2 CSI | S2 hit | S2 FAR | S1 CSI |
|---|---|---|---|---|
| Bolivia | 0.75 | 0.93 | 0.20 | 0.50 |
| Ghana | 0.63 | 0.66 | 0.08 | 0.33 |
| India | 0.70 | 0.89 | 0.24 | 0.47 |
| Mekong | 0.93 | 0.97 | 0.04 | 0.74 |
| Nigeria | 0.89 | 0.93 | 0.04 | 0.53 |
| Pakistan | 0.42 | 0.70 | 0.49 | 0.26 |
| Paraguay | 0.57 | 0.63 | 0.15 | 0.40 |
| Somalia | 0.56 | 0.79 | 0.34 | 0.15 |
| Spain | 0.79 | 0.90 | 0.13 | 0.63 |
| Sri Lanka | 0.86 | 0.93 | 0.08 | 0.19 |
| USA | 0.72 | 0.85 | 0.17 | 0.36 |

What this decided:

* **Edge-based Otsu** is the default: per-chip Otsu on scenes with little
  water splits land classes (tile Otsu false-alarm ratio 0.51 for MNDWI,
  0.58 for VV).
* **NDWI + confident MNDWI** is the default index: MNDWI alone over-maps
  (area bias 1.49 — wet soil, bright built-up); NDWI alone under-maps slightly
  (0.95). Adding only the MNDWI pixels above max(threshold, 0.3) gives the
  highest CSI with an unbiased area (1.01), which matters for depth: a wide
  extent raises the water surface, a narrow one lowers it. The 0.3 level is
  the one calibrated parameter (446 chips, 11 events).
* The **plausibility guard** refused 90 S2 chips of which only 3 had more than
  5 % water: it refuses water-free scenes, not flood scenes.
* **Sentinel-1** is the weaker path (CSI 0.45; hit rate 0.51 — flooded
  vegetation and urban flooding are not dark in VV). For context, an Otsu
  baseline of IoU 0.546 is reported for the Sen1Floods11 *test split*
  (Bonafilia et al. 2020, as tabulated by Garg et al. 2023); that figure uses a
  different subset and averaging and is not directly comparable with the
  pooled score here. Prefer optical imagery when it is cloud-free.

### 9.3 Depth — known water surfaces on real DEMs (`benchmarks/dem_recovery.py`)

The engine is given only an extent and a DEM; the true depth is known because
the water surface is imposed:

* **valley trend** — the DEM smoothed over 1.5 km plus an offset flooding 10 %
  or 25 % of the area (a stress test with deep water);
* **HAND** — water h = 1 m or 3 m above the drainage on GLO-30 HAND (Copernicus
  DEM, from the AWS open-data buckets): wet = HAND < h connected to the channel,
  depth = h − HAND.

Sites: Nyando / Kano plains (flat floodplain) and Mai Mahiu (escarpment,
gorge); DEMs: FABDEM, ALOS AW3D30 and TanDEM-X (from the project's national
mosaics) and Copernicus GLO-30. *Same DEM* isolates the method; *other DEM*
adds the DEM-choice error and tests the uncertainty band (P10–P90, 12 members).

Same DEM — mean absolute error and bias (m):

| Site | Scenario | True mean | Default (IDW, edge halfway, cap 1 m) | Edge at the wet cell (FwDET) | Nearest edge |
|---|---|---|---|---|---|
| Nyando | valley trend 10 % | 3.88 | 1.07 (+0.37) | 1.10 (+0.05) | 2.99 |
| Nyando | valley trend 25 % | 2.84 | 1.10 (+0.74) | 1.07 (+0.51) | 2.40 |
| Nyando | HAND 1 m | 0.67 | **0.24 (−0.06)** | 0.38 (−0.29) | 0.40 |
| Nyando | HAND 3 m | 2.26 | **0.70 (−0.52)** | 0.97 (−0.91) | 1.08 |
| Mai Mahiu | valley trend 10 % | 20.7 | 3.39 (−0.05) | 3.52 (−0.64) | 13.27 |
| Mai Mahiu | valley trend 25 % | 22.4 | 5.38 (+0.39) | 5.47 (−0.26) | 13.72 |
| Mai Mahiu | HAND 1 m | 0.59 | **0.51 (+0.09)** | 0.63 (−0.31) | 0.63 |
| Mai Mahiu | HAND 3 m | 1.73 | **1.10 (−0.81)** | 1.43 (−1.29) | 1.40 |

Other DEM (default settings) — share of true depths inside the ensemble
P10–P90:

| Site | Scenario | ALOS | TanDEM-X | other |
|---|---|---|---|---|
| Nyando | valley trend 10 % / 25 % | 87 % / 84 % | 82 % / 80 % | Copernicus 82 % / 81 % |
| Nyando | HAND 1 m / 3 m | 85 % / 80 % | 91 % / 85 % | FABDEM 66 % / 65 % |
| Mai Mahiu | valley trend 10 % / 25 % | 81 % / 85 % | 81 % / 84 % | Copernicus 96 % / 86 % |
| Mai Mahiu | HAND 1 m / 3 m | 86 % / 74 % | 89 % / 75 % | FABDEM 79 % / 60 % |

What this decided and shows:

* **Inverse distance** is the default: two to four times better than
  nearest-edge allocation where the water surface varies across a wide extent
  (valley trend), and as good or better in the HAND scenarios.
* **Edge water level halfway up the step (cap 1 m)** is the default: it lowers
  the error in 7 of 8 scenarios and cuts the flat-floodplain underestimate by
  40–80 % (Nyando HAND 3 m: bias −0.91 → −0.52 m; HAND 1 m: −0.29 → −0.06 m).
  The cost is a small overestimate in the Nyando valley-trend cases (+0.3 m).
  Uncapped, cliffs at Mai Mahiu inflated depths by 3 m.
* **Tried and dropped — a HAND-shaped surface.** HAND computed from the DEM
  (QEHT routing, streams ≥ 0.1 km²; r = 0.91–0.94 against ASF GLO-30 HAND) was
  used to interpolate the edge height above drainage instead of the edge
  elevation, chosen per flood region by cross-validation. It did not reduce the
  underestimate (Nyando HAND 3 m bias −0.91 m with either surface) and the
  cross-validation kept the elevation surface almost everywhere: the cause of
  the underestimate is where the waterline sits within the edge cells, not the
  shape of the surface.
* The **method-error term** (blocked cross-validation of the waterline) brings
  the P10–P90 coverage to about 77–95 % for valley-type floods (nominal 80 %).
* **Remaining bias — wide flat floodplains.** With water 3 m above the drainage
  over hundreds of km², depths are still about 0.5–0.8 m low and the band
  covers about 60–85 % of true depths. Treat depths of broad river-fed sheets as a
  lower estimate.
* **Coverage varies between ensembles.** The DEM errors are spatially
  correlated, so one 12-member ensemble is effectively a few independent
  samples: at Mai Mahiu the coverage of a case moved by up to ±12 points
  between random seeds (four seeds, 0.4.1 and the 0.4.2 ensemble alike). The
  ranges above are means over seeds; the depth error itself does not depend on
  the seed. In 0.4.2 each member has its own random stream and the noise field
  no longer wraps around the window edges; mean coverage was unchanged for
  valley floods and 5–8 points higher for HAND-type floods.

### 9.4 Case studies

All three were run with *Quick flood depth* at default settings (20-member
ensemble); inputs are the project data as downloaded (Copernicus Browser
"Raw" bands, national DEM mosaics, WorldCover). They are demonstrations of
behaviour, not accuracy assessments: there is no independent depth
measurement for these events.

| | Nyando, Kano plains | Mai Mahiu | Buulo Gaduud–Afmadow (Juba) |
|---|---|---|---|
| Event / image | 5 May 2024, S2 L2A Raw (8-bit + alpha) | 29 Apr 2024 flash flood; S2 7 Apr / 12 May (float) | Dec 2023 Deyr floods; S2 3 Dec (uint16 × 65 535) |
| Mode, index | water; NDWI + confident MNDWI | flood trace (dNDVI), kept connected to the Main Reach | water; NDWI + confident MNDWI |
| Threshold (η) | −0.239 (0.75, good) | −0.080 tile Otsu (0.51, low) | 0.122 (0.73, moderate) |
| Permanent water | WorldCover class 80: 279 km² removed (Winam Gulf) | — | — |
| DEM (σ) | FABDEM 30 m (1.5 m) | TanDEM-X 30 m (2.0 m) | project clip (2.0 m assumed) |
| Flooded area | 71.4 km² (P10–P90 62–116) | 2.19 km² (1.1–24.6) | 16.7 km² (11.9–18.6) |
| Mean depth, all / interior | 0.17 / 0.18 m | 1.33 / 2.41 m | 1.10 / 1.22 m |
| Interior median, P90 | 0.06, 0.38 m | 1.43, 5.22 m | 0.89, 3.06 m |
| Mean depth P10–P90 (ensemble) | 0.21–0.25 m | 1.59–2.32 m | 1.03–2.05 m |
| Method error (cross-validation) | 0.08 m | 2.52 m | 1.01 m |
| Evidence | observed water 70.9 km², steep 0.5 km² | erosion / deposition 1.32, inundation 0.14, steep flow path 0.73 km² | observed water 16.4 km², steep 0.3 km² |
| Mean edge lift (§ 6 step 3) | 0.02 m | 0.45 m | 0.13 m |
| Mean P90 − P10 by evidence | 0.53 m | 3.8 m / 3.0 m / 5.9 m | 3.1 m |
| Run time | 223 s | 16 s | 57 s |

*Nyando.* The mapped water is mostly a shallow sheet over the Ahero and West
Kano rice schemes and the fringes of the plains (interior median 0.04 m); the
coherent bodies along the Nyando reach and the south-west plains show
0.5–1.5 m. 41 % of the flood cells are waterline cells (depth ≈ 0 by
construction) because the 30 m extent is fragmented; the report shows the
interior statistics beside the all-cell ones for this reason. The edge lift
is small here (0.02 m): the Kano plains have almost no step at the flood edge. Without the
permanent-water layer the 279 km² of Winam Gulf were counted as flood with
zero depth — a permanent-water source is essential on lake shores. Given § 9.3,
these depths are more likely low than high.

*Mai Mahiu.* The post-event image shows no water (0.07 % NDWI > 0) and water
mode refuses it. Flood-trace mode maps the debris path from the culvert down
through the town; connectivity to the Main Reach drops 1 697 unrelated change
patches (crops, cloud). 61 % of the trace is bare soil or fresh sediment and
33 % lies on ground steeper than 10 %. The ensemble flags it plainly: the area
varies 22-fold across the threshold band (the NDVI-loss threshold is weakly
constrained, η 0.51), the method error is 2.5 m, and the bands are widest on
the steep flow path (5.3 m). The deepest cells (up to ~30 m) are in the
incised gorge, where the trace edge is the gully rim: an upper bound, not a
water depth. This is the intended reading of a mixed erosion-and-flooding
event — where the evidence is weak the output says so.

*Juba (Buulo Gaduud–Afmadow).* The Browser uint16 export is rescaled
(reflectance × 65 535), so no Sentinel-2 offset is applied (the report says
why). 16.7 km² of water on the floodplain within the scene; interior mean
1.22 m, P90 3.06 m; DEM error is the largest term (σ assumed 2 m because the
file name does not name the product — set it when known).

## 10. Limitations and good practice

* **Not a hydraulic model.** No flow, velocity, timing, backwater, or
  structures. For design, calibrate a hydraulic model; use FloodDepth for rapid
  assessment, damage screening and to check or seed a model.
* **The extent drives everything.** Wet soil, shadow, dark roofs or flooded
  vegetation misclassified as water raise the surface; flooded vegetation and
  urban flooding missed by dark-water SAR lower it. Read the mask before the
  depths.
* **DEMs.** Prefer a bare-earth DEM (LiDAR, FABDEM) to a surface model (GLO-30,
  ALOS, TanDEM-X): in a DSM the edge sits on canopy and roofs. Global 30 m DEMs
  smooth channels — channel depths are underestimated.
* **Wide flat floodplains.** The edge method still underestimates depth (about
  0.5–0.8 m with water 3 m above the drainage) when the water surface follows
  the drainage across a broad plain (§ 9.3).
* **Timing.** The image shows the water at the overpass, which may be before or
  after the peak.
* **Flood traces** mix inundation, scour and deposition; treat depths in
  erosion / deposition cells as indicative.

## 11. Implementation notes

* `core/` never imports QGIS; `processing_provider/` never implements the
  method. Vectors cross the boundary as WKT + CRS WKT, rasters as file paths.
* No subprocess, no pip dependency: NumPy and the GDAL / OGR bindings shipped
  with QGIS. GDAL C++ routines used: Warp, FillNodata, ComputeProximity,
  Polygonize, RasterizeLayer. FillNodata's search stops at twice the farthest
  flood cell's distance to its waterline (+10 cells): its cost grows with the
  search distance for every cell of the window.
* **Optional SciPy** (`core/accel.py`): when SciPy is importable (most QGIS
  installs ship it) connected components, box filters, the nearest-edge
  transform and the FFT use `scipy.ndimage` / `scipy.fft`; otherwise the
  NumPy / GDAL paths run. Both are tested against each other
  (`tests/test_accel.py`) and CI runs every suite both ways
  (`FLOODDEPTH_NO_SCIPY=1` forces the fallbacks). The report says which ran.
* Algorithm ids of 0.3 are kept (`classify_inundation`,
  `estimate_stage_dem_boundary`, `compute_flood_depth`,
  `compute_constant_stage_depth`).
* Copernicus Browser 'Raw' Sentinel-1 in UINT16 is linear backscatter × 65535
  and is read as such (it was read as dB before 0.4.2 and no water was found).
* Tests: `tools/dev/ci_pure.sh` (pure suites), `tests/qgis_smoke.py` (every tool
  through Processing on a synthetic site with a known answer). Benchmarks
  (`benchmarks/`) are in the repository but not in the plugin zip.

## 12. References

* Bonafilia, D., Tellman, B., Anderson, T., Issenberg, E. (2020). Sen1Floods11:
  a georeferenced dataset to train and test deep learning flood algorithms for
  Sentinel-1. CVPR Workshops, 210–211.
* Cohen, S., Brakenridge, G. R., Kettner, A., et al. (2018). Estimating
  floodwater depths from flood inundation maps and topography. JAWRA 54(4),
  847–858.
* Cohen, S., Raney, A., Munasinghe, D., et al. (2019). The Floodwater Depth
  Estimation Tool (FwDET v2.0) for improved remote sensing analysis of coastal
  flooding. NHESS 19, 2053–2065.
* Donchyts, G., Schellekens, J., Winsemius, H., Eisemann, E., van de Giesen, N.
  (2016). A 30 m resolution surface water mask including estimation of
  positional and thematic differences using Landsat 8, SRTM and OpenStreetMap.
  Remote Sensing 8(5), 386.
* Garg, S., Feinstein, B., Timnat, S., et al. (2023). Cross modal distillation for flood
  extent mapping. arXiv:2302.08180.
* Hawker, L., Uhe, P., Paulo, L., et al. (2022). A 30 m global map of
  elevation with forests and buildings removed (FABDEM). Environ. Res. Lett.
  17, 024016.
* Martinis, S., Twele, A., Voigt, S. (2009). Towards operational near real-time
  flood detection using a split-based automatic thresholding procedure on high
  resolution TerraSAR-X data. NHESS 9, 303–314.
* McFeeters, S. K. (1996). The use of the Normalized Difference Water Index
  (NDWI) in the delineation of open water features. IJRS 17(7), 1425–1432.
* Nobre, A. D., Cuartas, L. A., Hodnett, M., et al. (2011). Height Above the
  Nearest Drainage – a hydrologically relevant new terrain model. J. Hydrol.
  404, 13–29.
* Otsu, N. (1979). A threshold selection method from gray-level histograms.
  IEEE Trans. SMC 9(1), 62–66.
* Pekel, J.-F., Cottam, A., Gorelick, N., Belward, A. S. (2016). High-resolution
  mapping of global surface water and its long-term changes. Nature 540,
  418–422.
* Rikimaru, A., Roy, P. S., Miyatake, S. (2002). Tropical forest cover density
  mapping. Tropical Ecology 43(1), 39–47.
* Rong, G., Tan, T.-S. (2006). Jump flooding in GPU with applications to
  Voronoi diagram and distance transform. Proc. I3D 2006, 109–116.
* Xu, H. (2006). Modification of normalised difference water index (NDWI) to
  enhance open water features in remotely sensed imagery. IJRS 27(14),
  3025–3033.
* Zanaga, D., et al. (2022). ESA WorldCover 10 m 2021 v200. Zenodo.
