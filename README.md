# FloodDepth — flood extent, depth and uncertainty for QGIS

[![tests](https://github.com/EdmondAkello/Flood-Depth/actions/workflows/tests.yml/badge.svg)](https://github.com/EdmondAkello/Flood-Depth/actions/workflows/tests.yml)

Quick, defensible flood-depth estimates from a satellite image of a flood (or a
flood-extent layer) and a DEM, with uncertainty bands, computed entirely inside
QGIS.

Author: Edmond Akello · License: GPL-3.0-or-later · QGIS 3.28 – 4.x · version 0.4.0 (experimental)

## What it does

Select the image files you downloaded, select a DEM, click Run:

1. **Recognises the bands** from file names and metadata — Copernicus / EO
   Browser "Raw" exports (`…_B03_(Raw).tiff`), Sentinel-2 SAFE / HLS, Landsat
   Collection 2 (`…_SR_B3.TIF`), multiband stacks with band names, the Browser
   "False color" product, Sentinel-1 VV/VH. Alpha / data-mask bands,
   undeclared NoData sentinels and the Sentinel-2 −1000 DN offset are handled.
2. **Classifies water** with MNDWI, NDWI or SAR backscatter and a
   **tile-based Otsu** threshold, reporting its separability (η). A threshold
   outside the physically usual range is **refused with an explanation**, not
   used — a post-event image with no water left does not become a flood map.
3. **Removes permanent water** using a pre-event image, ESA WorldCover or JRC
   occurrence, so lakes and rivers are not counted as flood.
4. **Flood-trace mode** for floods that drained before the overpass (flash
   floods, debris flows): maps vegetation lost between two dates and splits it
   into *inundation* and *erosion / deposition* evidence.
5. **Estimates the water surface** from the terrain along the flood edge
   (FwDET boundary method) — edges against NoData, cloud or permanent water are
   not used, edge outliers (embankments, canopy) are screened, steep ground is
   protected — and computes **depth = surface − DEM**. A HAND raster can be
   used instead of the DEM.
6. **Quantifies uncertainty** with a Monte Carlo ensemble (spatially
   correlated DEM error, threshold band, edge shifts, extra terrain error where
   the ground was eroded): **depth P10 / P50 / P90** and **flood probability**,
   plus an **evidence layer** saying what each flooded cell is based on.
7. Writes an **HTML report** with every automatic decision, statistics,
   warnings and limitations.

Every computation uses NumPy and the GDAL bindings QGIS ships. Nothing is
installed with pip and no external process is started, so it runs on managed
workstations where SAGA / GRASS / WhiteboxTools calls are blocked.

## Tools (Processing Toolbox → FloodDepth)

| Tool | Use |
|---|---|
| **Quick flood depth (image or extent + DEM)** | The one-click workflow: classification → extent → water surface → depth → uncertainty → report, in a fixed output folder |
| 1. Classify water from imagery (automatic) | Water mask + polygons from image files (water or flood-trace mode) |
| 1b. Classify water from training polygons | RGB, drone or aerial images; scenes where indices fail |
| 2. Water surface and depth from a flood extent | Polygons or a mask raster + DEM (or HAND) → WSE and depth rasters |
| 3. Flood depth from a water-surface raster | A known surface (hydraulic model, surveyed marks) − DEM |
| Flood depth below a constant water level | Lake / reservoir level, flood mark, what-if level |
| Compare flood extents (CSI, hit rate, false alarms) | Validate an extent against a reference |

All tools work in batch mode, the Graphical Modeler and `processing.run()`.
The 0.3 algorithm ids are kept.

## Quick start

1. In the Copernicus Browser, open the post-event Sentinel-2 L2A scene, *Download → Analytical*,
   tick the **Raw** bands B02, B03, B04, B08, B11 (GeoTIFF). Optionally do the same for a clear
   pre-event date.
2. Get a DEM covering the area (FABDEM or Copernicus GLO-30 are good global choices; a LiDAR DTM
   is better). A geographic DEM is fine: the analysis runs on the WGS 84 / UTM grid of the area.
3. Run **Quick flood depth**: *Satellite image file(s)* = the band files, *DEM* = the DEM,
   *Pre-event image file(s)* = the pre-event bands (or a WorldCover tile as *Permanent water
   raster*), and an output folder.
4. Open `report.html` in the output folder and read the *Check these* box first.

Outputs (output folder): `depth.tif`, `wse.tif`, `evidence.tif`, `depth_p10.tif`,
`depth_p50.tif`, `depth_p90.tif`, `flood_probability.tif`, `water_mask.tif`,
`flood_extent.gpkg` (layers `flood_extent`, `waterline`), `report.html`, `run.json`.

## Design

    flooddepth/
      core/                    ← NO QGIS IMPORTS. NumPy + GDAL/OGR only.
        raster.py              GDAL I/O, NoData conventions, analysis grid, warping
        sensors.py             band recognition from names / metadata
        classify.py            indices, Otsu (tile-based, separability, band), components, cleanup
        imagery.py             files → water / trace mask (permanent water, cloud, offsets)
        surface.py             waterline, outliers, nearest-edge allocation (jump flooding), depth
        uncertainty.py         evidence classes, correlated DEM error, Monte Carlo ensemble
        pipeline.py            the quick workflow (testable without QGIS)
        compare.py, stats.py   validation scores, depth statistics
        report.py, vector.py   HTML report, polygons ↔ grids
      processing_provider/     ← the only place QGIS and core meet
      tests/                   ← pure suites (no QGIS) + qgis_smoke.py
      tools/dev/               ← test runners, zip builder

The core is importable without QGIS:

    python -m flooddepth.tests.test_surface     # analytic surfaces with known depths
    python -m flooddepth.tests.test_pipeline    # the whole workflow on a synthetic site

`tools/dev/ci_pure.sh` runs the 7 pure suites (118 checks); `tools/dev/run_tests.sh` adds the
QGIS smoke test (`tests/qgis_smoke.py`, 32 checks: every tool through real Processing).
GitHub Actions runs the pure suites and a Bandit scan on every push.

## Method, assumptions and limits

See [DOCUMENTATION.md](DOCUMENTATION.md) for the method, every automatic decision, the
uncertainty model, validation results and references. In short:

* Depths are **terrain-based estimates from the flood edge, not a hydraulic simulation**.
* An extent that is too large (wet soil, shadow, vegetation classed as water) raises the
  surface; one too small lowers it. The uncertainty bands and the evidence layer show where.
* Accuracy is bounded by the DEM. Global 30 m DEMs are surface models with metre-level error
  and smoothed channels: channel depths are underestimated.
* A **flood trace is not a water extent**: it mixes inundation, scour and deposition; where the
  ground was eroded the pre-event DEM no longer describes it.

## Install

From a release zip: *Plugins → Manage and Install Plugins → Install from ZIP*. To build the zip
from a clone: `bash tools/dev/make_zip.sh` → `dist/FloodDepth_<version>.zip`.

## Citation

Akello, E. (2026). *FloodDepth: flood extent, depth and uncertainty from satellite imagery and
a DEM for QGIS* (version 0.4.0). https://github.com/EdmondAkello/Flood-Depth
