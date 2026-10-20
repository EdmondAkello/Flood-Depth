# Audit of FloodDepth 0.3.0 (ChatGPT prototype)

Date: 2026-10-20. Environment: QGIS 3.34.4, GDAL 3.8.4, Python 3.12 (Ubuntu), offscreen.
Method: the 0.3.0 zip was installed as a Processing provider and every algorithm was run through
`processing.run()` on a synthetic 200 × 200 and 600 × 600 site, then on the case-study data
(Nyando 2024, Mai Mahiu 2024). The inherited claims ("6 unit tests pass, compilation passes")
were true but did not cover the stage or depth tools or any QGIS execution.

| # | Severity | Module | Finding | Reproduction | Resolution in 0.4.0 |
|---|---|---|---|---|---|
| 1 | Critical | `algorithms/estimate_dem_stage.py` | `QgsGeometry.boundary()` does not exist in the QGIS API: the stage tool raises `AttributeError` on every input. The image → depth chain has never run. | Any polygon + DEM | Replaced by a grid-based waterline (no vector boundary needed) |
| 2 | Critical | same | With #1 patched, IDW is evaluated per pixel in a Python double loop over the **whole DTM**, not the flood: 2.1 s for 0.04 Mcells, 19 s for 0.36 Mcells (with SciPy). A 10 m scene of 10 Mcells ≈ 9 min; without SciPy the fallback sorts every sample per pixel (hours). | 600 × 600 synthetic | Whole-array NumPy / GDAL over the flood's bounding box only |
| 3 | Critical | `algorithms/classify_inundation.py` | Requires one multiband file with ≥ 3 bands. Browser "Raw" exports, Sentinel-2 SAFE and Landsat are one file per band; single-band SAR is rejected. Default band numbers (G = 2, NIR = 4, SWIR = 5) are wrong for a Sentinel-2 stack (B03 = 3, B08 = 8, B11 = 11/12). | Case-study imagery | Band recognition from names / metadata; multiple files; overrides only for unnamed stacks |
| 4 | High | `estimate_dem_stage.py` | Every polygon boundary is treated as the waterline, including image edges, AOI edges, NoData and cloud edges. | Flood touching the image edge | Waterline only against valid dry ground |
| 5 | High | all | No grid alignment: classification is written in the image CRS (EPSG:4326 for Browser exports), `compute_flood_depth` requires identical geotransforms, the stage tool rejects geographic DEMs (Copernicus, HAND). | Case-study data | One projected analysis grid; automatic warping |
| 6 | High | `classify_inundation.py` | Sieve filter also fills dry holes (islands, roads, buildings) and the original mask is overwritten. | Flood with islands | Separate, reported water-patch removal and hole filling; holes off by default |
| 7 | High | `core/image_classification.py` | Otsu without a separability check: on a scene without water (Mai Mahiu, 12 May 2024) it still classifies 7.5 % "water" at MNDWI −0.41. | Mai Mahiu post-event bands | η reported; plausibility guard refuses the threshold and points to flood-trace mode |
| 8 | High | `classify_inundation.py` | Sentinel-2 L2A −1000 DN offset (baseline 04.00, 2022+) not handled; alpha / dataMask bands ignored (NoData classified). | Browser Raw exports | Offset from metadata or data; alpha masks read |
| 9 | Medium | `compute_depth.py` | Opens the mask with `ogr.Open(layer.source())` and fails for memory layers although the features are then read through QGIS. | Memory polygon layer | Geometries passed as WKT |
| 10 | Medium | `plugin.py`, `metadata.txt` | About dialog says 0.1.0 and lists no classification; empty icon, no tracker, no `hasProcessingProvider`. | — | Rewritten |
| 11 | Low | `core/image_classification.py` | Otsu returns the first bin of a flat maximum: with an empty-bin gap (8-bit data) the threshold sits on the edge of the land mode. | Two separated populations | Middle of the flat maximum |
| 12 | Low | tests | 6 numerical tests, none for stage, depth or QGIS. | — | 7 pure suites (118 checks) + QGIS smoke (32 checks) |

The supervised nearest-centroid classifier and the constant-stage and depth-from-surface
calculations were correct and are kept (ported to the new core).
