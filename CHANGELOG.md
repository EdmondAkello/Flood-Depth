# Changelog

All notable changes to FloodDepth are recorded here.

## [Unreleased]

### Added
- **Interactive 3D scene** (`scene3d.html`, written by Quick flood depth; on by default): the terrain draped
  with the satellite image (or a hillshade), permanent water and doubtful cells painted on it, the water
  surface coloured by depth, **mud and debris** (flood-trace erosion / deposition) as their own surface,
  the flood edge, connectivity features and area of interest. One offline HTML file (WebGL 2, the
  LandXML Road & Terrain viewer): rotate, layers and opacity, vertical exaggeration, click to read level
  and depth, views, PNG and A3 export. **Export 3D scene** rebuilds it from a run folder with extra
  layers, another level of detail or exaggeration. The report links to it.
- **Evidence map in the report** (mud and debris in brown) and a **mud and debris colour** among the
  map colours; `evidence.tif` uses the chosen colours.
- **Raster menu and toolbar.** *Raster > FloodDepth* lists Quick Flood Depth, Flood Depth over
  Time, Prepare DEM, Classify Water, Compare Flood Extents and the User Guide; one button on the
  Raster toolbar opens Quick flood depth on a click and lists the other tools on its arrow.
- **User guide installed with the plugin** (`help/user-guide.html`, offline): opened from the menu
  and from every tool's Help button at the tool's section. PDF in `docs/user-guide/`.
- **DEM tiles.** Quick flood depth and the time-series tool take *More DEM tiles*; tiles are
  combined into a virtual mosaic (`core/dem.py`), and when they come in several coordinate systems
  (the Google GEM 15 m DTM ships Kenya in four UTM zones) the one covering the area best is used.
- **Prepare DEM** tool: merge tiles, clip with a buffer, reproject to the local UTM zone with
  bilinear (or cubic) resampling, float32; a DEM already in the target CRS keeps its cells. Nearest
  neighbour is not offered.
- **DEM check** in every run and in Prepare DEM: repeated rows / columns (nearest-neighbour
  resampling upstream) and whole-metre storage on flat ground are reported (after QEHT's DEM audit).
- **Google GEM 15 m DTM** DEM-error preset (2.0 m) from the file name.
- **Map colours**: depth ramp (Blues, Viridis, Yellow-green-blue, Magma), the depth at its dark
  end, and the colours of permanent water, the flood edge and doubtful cells, for the report maps
  and the loaded depth layer (advanced parameters; saved in run.json).
- **Memory check** before the DEM is read: a warning above half the free memory, a refusal with
  advice above 90 % (about 64 bytes per analysis cell, measured).
- **Flood depth over time (several image dates)**: a new tool. Image files of several dates are
  grouped by the date in their names; each date runs Quick flood depth in its own folder with the
  same DEM and settings, then `max_depth.tif` (deepest water seen on any date), `dates_flooded.tif`,
  `first_date_flooded.tif`, `series.csv` and `series.html` (area, depth and volume over time with
  P10-P90, charts spaced by time, the maximum-depth map). Failed dates are reported and skipped;
  dates mixing sensors / indices are flagged. Nyando 30 Apr - 6 May 2024 (S2 False color, S2 Raw,
  S1 Raw): 3 dates in 48 s.
- **Evidence class 5, "doubtful - ground above the water surface"**: mapped flood cells whose
  ground lies above the estimated surface (depth 0; wet soil, crops, vegetation or DSM artefacts),
  in `evidence.tif` and in pink on the report's depth map.
- **Optional SciPy acceleration** (`core/accel.py`): connected components, box filters, the
  nearest-edge transform and the FFT use SciPy when it is importable; the NumPy / GDAL paths are
  kept and tested against it, and CI runs every suite both ways (`FLOODDEPTH_NO_SCIPY=1`).
- **Uncertainty members run in parallel threads** (CPU count - 1, at most 8; `FLOODDEPTH_THREADS`).
  Each member has its own random stream, so results do not depend on the thread count. Correlated
  DEM noise uses FFT-friendly sizes. Nyando quick run 52 -> 24 s on 2 cores.
- Note when the pre-event image is a rendered colour composite (the Raw bands give a cleaner
  permanent-water mask).

### Changed
- **Large DEMs and images are read in windows.** The DEM (and HAND, permanent-water rasters) is
  read only over the analysis grid, not whole: a 16 384 x 16 384 tile went from 5.3 s and 1.1 GB to
  0.3 s for a 2 000 x 2 000 area; a whole 32 768 x 32 768 GEM tile (4.3 GB as float32) no longer
  has to fit in memory. With an area of interest the image is read only around it (the area plus
  half its size each side, for the threshold); Nyando with a 30 % x 30 % AOI: 37 s -> 14 s.
- **Water surface per flood cluster.** FillNodata runs on the boxes around groups of flood
  patches instead of one box over all of them (identical result; 14x faster on two distant
  patches, the Juba layout).
- **Report and series page in the LandXML Road & Terrain look**: title band with the run facts,
  numbered sections with contents, indigo / teal / amber cards, tabular figures, A4 print layout
  with page numbers and a footer.
- **QGIS 3.44 to 4.x**: minimum version 3.44; Processing enums in their `Qgis.*` form
  (`ProcessingSourceType`, `ProcessingParameterFlag`, `ProcessingNumberParameterType`,
  `CrsWktVariant`, `ShaderInterpolationMethod`).
- Quick flood depth dialog: pre-event image and mode next to the image, area of interest earlier,
  connectivity features last; shorter help.
- CI: QGIS smoke test on QGIS 3.44.15 and 4.2.3 (Docker), flake8 (`.flake8`, 130 columns) and
  Bandit; `tools/dev/run_tests.sh` runs the same locally.

### Fixed
- **Flood polygons are valid geometries.** 8-connected regions touched themselves at cell corners, and
  QGIS zonal statistics, clip and dissolve refused them; they are now split at the corner.
- Step-by-step depth tools on a whole large DEM now stop with advice (area of interest) instead of
  trying to read it all.
- Unused imports, ambiguous names and lambdas flagged by flake8.
- **Copernicus Browser Sentinel-1 'Raw' in UINT16** (linear backscatter x 65535) was read as dB,
  so every pixel looked bright and no water was found (Nyando 6 May 2024). Now scaled to linear.
- **Bright-cloud test scaled band by band.** In a Copernicus Browser uint16 export (reflectance
  x 65535) a dark blue band (dry soil) sits below 12000 and was read as x 10000, 6.5 times too
  bright, so clear land was masked as cloud. The scale is now decided once per scene from the
  brightest optical band. Juba (3 Dec 2023, clear): cloud mask 60.7 % -> 0.2 %, flooded area
  16.7 -> 20.7 km2.
- The report said a pre-event image was used when its bands were not found.
- **Pre-event image without SWIR1** (a Browser 'False color' export) was rejected because the
  post-event index (NDWI + MNDWI) needs SWIR1; it now falls back to NDWI for permanent water.

### Performance
- **Water-surface interpolation up to 50x faster on large grids.** `gdal.FillNodata` was asked to
  search the whole window for every cell; the search now stops at twice the farthest flood
  cell's distance to its waterline (+10 cells). Surface unchanged to millimetres. Quick run:
  Nyando 220 -> 52 s, Juba (7.8 M cells) 288 -> 33 s. Same cap in the method cross-validation.
- Per-stage run times in `run.json` (`timings`) and in the report.

### Changed
- **Report redesigned.** Maps are framed on the main flood area (the smallest window holding 80 %
  of the flooded cells) with an overview inset, scale bar and north arrow, drawn over the
  satellite image itself (true or false colour, read from overviews) or, without one, over a
  relief-stretched hillshade with an elevation tint so flat floodplains are not a uniform grey.
  Permanent water and NoData / cloud have their own symbols. A flood-probability map sits beside
  the depth map. Page: section cards, KPIs with their P10-P90 ranges, result warnings separated
  from input notes, skipped rendered products collapsed into one line, inputs listed by file
  name per folder.
- **Report prints on A4 portrait** (about 4 pages): summary and depth map, probability map with depth
  and uncertainty, evidence / surface / classification, inputs / log / limitations. Collapsed
  sections open for printing; colours and legends print.
- Report maps are JPEG instead of PNG and embedded once: the Juba report fell from 6.5 MB to
  well under 1 MB.

## [0.4.1] — 2026-10-20

### Changed
- **Edge water level halfway up the step** to the first dry cell (capped at 1 m) instead of at
  the wet cell's ground. Depth-recovery benchmark: lower error in 7 of 8 scenarios; the
  flat-floodplain underestimate falls by 40–80 % (Nyando, water 3 m above drainage: bias −0.91 →
  −0.52 m). Advanced parameters *Edge water level* (0 = FwDET convention) and *maximum lift*.
- CI byte-compiles every file before the suites (the Processing layer is not imported by them).

### Removed
- The HAND-shaped surface trial merged from `floodplain-hand` (`core/hand.py`, `core/hydro/`,
  per-region surface choice): it did not reduce the underestimate (DOCUMENTATION.md § 9.3).

## [0.4.0] — 2026-10-20 (redesign)

The 0.3 package was audited in QGIS 3.34 before this release
([docs/AUDIT_0.3.0.md](docs/AUDIT_0.3.0.md)). Its stage tool failed on every
input, so the image-to-depth workflow had never run end to end; classification
required a single multiband file and could not read the usual downloads. 0.4.0
is a redesign on the QEHT architecture (pure core, thin Processing layer).

### Added
- **Quick flood depth** — one tool from image file(s) or a flood-extent layer and a DEM to
  depth, uncertainty and an HTML report, in a fixed output folder.
- **Band recognition** from file names and metadata: Copernicus / EO Browser Raw exports,
  Sentinel-2 SAFE / HLS, Landsat 4–9 Collection 2, multiband stacks with descriptions or
  wavelengths, 12/13-band Sentinel-2 stacks named as such, the Browser False-color composite,
  Sentinel-1 VV/VH. True-colour images and rendered SAR visualisations are refused with
  guidance. Band numbers can still be set for unnamed stacks.
- **NoData conventions**: alpha / dataMask bands, undeclared float sentinels (−3.4e38), DEMs
  with 0 fill and no declared NoData.
- **Sentinel-2 −1000 DN offset** applied from metadata, or only when every band's values lie in
  the offset DN window (Browser exports are harmonised or rescaled and are left alone).
- **Edge-based Otsu** (Donchyts et al. 2016) seeded at the physical threshold (0; −16 dB VV;
  −22 dB VH) as the default — chosen on the 446 Sen1Floods11 chips (CSI 0.61 S2, 0.45 S1 vs
  0.46 / 0.31 for tile-based); tile-based (Martinis et al. 2009) and global Otsu as options;
  a 5 × 5 speckle filter for SAR. Separability η (single Gaussian = 2/π),
  the near-optimal threshold band, a flat-maximum fix (middle of an empty-bin gap) and a
  **plausibility guard** that refuses an automatic threshold outside the range where water
  separates (NDWI/MNDWI −0.3…0.6; VV −26…−10 dB; VH −32…−16 dB).
- **NDWI + confident MNDWI** as the automatic Sentinel-2 / Landsat index (NDWI water, or MNDWI >
  max(threshold, 0.3)): Sen1Floods11 CSI 0.74 with an unbiased area (1.01), against 0.61 / bias
  1.49 for MNDWI alone.
- **Interior depth statistics** (waterline ring excluded) beside the all-cell ones in the report.
- **Permanent water** from a pre-event image, a raster (ESA WorldCover class 80, JRC occurrence
  ≥ 50 %, or given values) or polygons.
- **Flood-trace mode** (pre/post NDVI loss) for floods that drained before the overpass, with
  bare-soil index to split inundation from erosion / deposition evidence.
- **Connectivity filter**: keep only flooding connected to a river, reach or road layer.
- **Water-surface engine** (FwDET boundary method): waterline only against valid dry ground,
  robust edge-outlier screening (steep ground protected), inverse-distance surface from the
  edge (default, chosen by benchmark) or nearest-edge allocation by jump flooding in NumPy
  (exact; schedule bounded with GDAL ComputeProximity), slope-aware smoothing; HAND method
  (depth = edge HAND − HAND).
- **Uncertainty**: evidence classes (observed water, trace-inundation, trace-erosion/deposition,
  steep flow path); Monte Carlo ensemble with correlated DEM error (presets from the DEM file
  name), a water-surface error from a blocked cross-validation of the waterline,
  threshold-band and edge perturbation, extra terrain error in eroded cells → depth
  P10/P50/P90 and flood probability; coarsened automatically above 2 Mcells.
- **Compare flood extents** (CSI, hit rate, false-alarm ratio, bias, F1, agreement raster).
- **HTML run report** and `run.json`.
- Automatic grid alignment: geographic DEMs run on the WGS 84 / UTM grid; images in any CRS are
  warped; masks move between grids by area fraction.
- Tests: 7 pure suites (120 checks) and a QGIS smoke test (32 checks) on a synthetic site with a
  known answer; GitHub Actions; Bandit clean.
- Benchmarks (`benchmarks/`, not in the zip): Sen1Floods11 classification; depth recovery on real
  DEMs with imposed water surfaces.

### Changed
- Stage estimation: the per-pixel Python IDW over the whole DTM (minutes to hours) is replaced
  by whole-array operations over the flood only (seconds).
- Cleanup no longer fills dry holes by default (islands, roads and buildings are kept); water
  patch removal and hole filling are separate settings and are reported.
- Licence headers: GPL-3.0-or-later throughout (as in 0.3).

### Fixed
- `estimate_stage_dem_boundary` crashed on every input (`QgsGeometry.boundary()` does not exist).
- Image, NoData and AOI edges were sampled as the waterline.
- `compute_flood_depth` refused memory polygon layers.
- The About dialog reported version 0.1.0; metadata had no icon, tracker or provider flag.

## [0.3.0] — ChatGPT prototype
Automatic NDWI/MNDWI Otsu classification, sieve cleanup, supervised centroid classification,
DEM-boundary IDW stage, depth from a water-surface raster or a constant stage.
