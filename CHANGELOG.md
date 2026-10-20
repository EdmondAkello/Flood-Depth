# Changelog

All notable changes to FloodDepth are recorded here.

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
