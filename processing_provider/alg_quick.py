# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Quick flood depth: satellite image(s) or a flood extent + DEM -> depth."""

from qgis.core import (
    Qgis,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterColor,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterMultipleLayers,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
    QgsProcessingOutputFile,
    QgsProcessingOutputRasterLayer,
    QgsProcessingOutputVectorLayer,
)

from qgis.PyQt.QtGui import QColor

from .base import FloodDepthAlgorithm, RasterStyle

MODES = ["water", "trace"]
INDEXES = ["auto", "consensus", "ndwi", "mndwi", "sar"]
METHODS = ["dem", "hand"]
INTERP = ["idw", "nearest"]
RAMPS = ["blues", "viridis", "ylgnbu", "magma"]


class QuickFloodDepthAlgorithm(FloodDepthAlgorithm):

    def name(self):
        return "quickflooddepth"

    def displayName(self):
        return "Quick flood depth (image or extent + DEM)"

    def group(self):
        return "Flood depth"

    def groupId(self):
        return "flooddepth"

    def shortHelpString(self):
        return (
            "From a satellite image (or a flood-extent layer) and a DEM to flood depth with uncertainty bands "
            "and a report, in one run.\n\n"
            "<b>1. Image.</b> Select the downloaded files as they are: Copernicus / EO Browser 'Raw' bands "
            "(B03, B08, B11 ...), Sentinel-2 SAFE, Landsat SR_B3 ..., a multiband stack, or Sentinel-1 VV / VH. "
            "Bands are recognised from file names and metadata. True-colour images and rendered SAR pictures "
            "are refused with the reason. Or give a flood-extent polygon layer instead.\n\n"
            "<b>2. DEM.</b> Any DEM GDAL reads: FABDEM, Copernicus GLO-30, ALOS, TanDEM-X, the Google GEM 15 m "
            "DTM or LiDAR. Tiled DEMs: pick one tile as the DEM and the rest under 'More DEM tiles' (only the "
            "cells under the image are read). Geographic DEMs are analysed in the local UTM zone.\n\n"
            "<b>3. Optional.</b> A pre-event image removes permanent water (and is required for flood-trace "
            "mode, which maps vegetation stripped by a flood that drained before the overpass). An area of "
            "interest limits the run and makes it faster. Permanent-water layers (ESA WorldCover, JRC "
            "occurrence, lakes) are excluded from the flood.\n\n"
            "<b>Method.</b> Water from NDWI + confident MNDWI (or SAR backscatter) with an edge-based Otsu "
            "threshold; implausible thresholds are refused. Depth = water surface - DEM, the surface taken from "
            "the terrain along the flood edge (FwDET), or from a HAND raster. A Monte Carlo ensemble with "
            "correlated DEM error and extent perturbation gives depth P10 / P50 / P90 and flood probability.\n\n"
            "<b>Outputs</b> (output folder): depth.tif, flood_extent.gpkg, evidence.tif, wse.tif, "
            "depth_p10/p50/p90.tif, flood_probability.tif, water_mask.tif, report.html (print it to PDF), "
            "run.json.\n\n"
            "Terrain-based estimates, not a hydraulic model: read the report's warnings.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterMultipleLayers(
            "IMAGES", "Satellite image file(s) (post-event)", Qgis.ProcessingSourceType.Raster, optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "EXTENT", "...or a flood-extent polygon layer instead of an image",
            [Qgis.ProcessingSourceType.VectorPolygon], optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", "DEM (terrain elevation)"))
        self.addParameter(QgsProcessingParameterMultipleLayers(
            "DEM_TILES", "More DEM tiles, mosaicked with the DEM (optional; e.g. Google GEM 15 m tiles)",
            Qgis.ProcessingSourceType.Raster, optional=True))
        self.addParameter(QgsProcessingParameterMultipleLayers(
            "PRE_IMAGES", "Pre-event image file(s) (removes permanent water; required for flood trace)",
            Qgis.ProcessingSourceType.Raster, optional=True))
        self.addParameter(QgsProcessingParameterEnum(
            "MODE", "What the image shows",
            ["Water at the overpass (ponding, standing flood)",
             "Flood trace after the water drained (needs a pre-event image)"], defaultValue=0))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "AOI", "Area of interest (faster: only this area is read)", [Qgis.ProcessingSourceType.VectorPolygon], optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "PERM_RASTER", "Permanent water raster (ESA WorldCover, JRC occurrence or a water mask)",
            optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "PERM_POLYS", "Permanent water polygons (lakes, rivers)", [Qgis.ProcessingSourceType.VectorPolygon],
            optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "SEEDS", "Keep only flooding connected to these features (river, road, reach)",
            [Qgis.ProcessingSourceType.VectorAnyGeometry], optional=True))
        self.addParameter(QgsProcessingParameterEnum(
            "METHOD", "Depth method", ["Flood-edge elevation on the DEM (FwDET)",
                                       "Flood-edge height on a HAND raster"], defaultValue=0))
        self.addParameter(QgsProcessingParameterRasterLayer("HAND", "HAND raster (HAND method)", optional=True))
        self.number("MEMBERS", "Uncertainty ensemble members (0 = off)", "int", 20, 0, 200)
        self.number("DEM_SIGMA", "DEM vertical error, 1 sigma (m; 0 = from the DEM name)", "double", 0.0, 0.0, 50.0)
        # advanced
        self.addParameter(self.advanced(QgsProcessingParameterEnum(
            "INDEX", "Water index", ["Automatic", "NDWI + confident MNDWI (green, NIR, SWIR1)",
                                     "NDWI (green, NIR)", "MNDWI (green, SWIR1)", "SAR backscatter (VV / VH)"], defaultValue=0)))
        self.number("THRESHOLD", "Manual threshold (empty = automatic)", "double", None, optional=True, adv=True)
        self.addParameter(self.advanced(QgsProcessingParameterEnum(
            "TMETHOD", "Automatic threshold", ["Edge-based Otsu (recommended)", "Tile-based Otsu", "Global Otsu"],
            defaultValue=0)))
        for key, label in (("B_GREEN", "green"), ("B_NIR", "NIR"), ("B_SWIR1", "SWIR1"), ("B_RED", "red"),
                           ("B_BLUE", "blue")):
            self.number(key, f"Band number of {label} in a multiband file (0 = automatic)", "int", 0, 0, 500, adv=True)
        self.number("MIN_WATER", "Remove water patches smaller than (image cells)", "int", 10, 0, adv=True)
        self.number("MAX_HOLE", "Fill dry holes up to (image cells; 0 = keep islands)", "int", 0, 0, adv=True)
        self.number("CLOUD", "Bright-cloud blue reflectance limit (0 = off)", "double", 0.25, 0.0, 1.0, adv=True)
        self.number("SEED_BUFFER", "Connection distance to the features (m)", "double", 30.0, 0.0, adv=True)
        self.addParameter(self.advanced(QgsProcessingParameterString(
            "PERM_VALUES", "Permanent-water raster values (auto, 80, 80,90 or >=50)", defaultValue="auto")))
        self.number("CELL", "Analysis cell size (m; 0 = DEM resolution)", "double", 0.0, 0.0, adv=True)
        self.addParameter(self.advanced(QgsProcessingParameterEnum(
            "INTERP", "Surface interpolation", ["Inverse distance from the flood edge (GDAL)", "Nearest flood edge (FwDET v1)"],
            defaultValue=0)))
        self.number("SMOOTH", "Smoothing radius (cells)", "int", 2, 0, 50, adv=True)
        self.number("SMOOTH_TOL", "Smoothing tolerance (m)", "double", 0.5, 0.0, adv=True)
        self.number("EDGE_STEP", "Edge water level: share of the step to the first dry cell (0 = FwDET, 0.5 = halfway)",
                    "double", 0.5, 0.0, 1.0, adv=True)
        self.number("EDGE_CAP", "Edge water level: maximum lift (m)", "double", 1.0, 0.0, adv=True)
        self.number("OUT_RADIUS", "Edge outlier window radius (cells)", "int", 10, 0, 200, adv=True)
        self.number("OUT_K", "Edge outlier limit (robust SD)", "double", 3.0, 0.0, adv=True)
        self.number("OUT_MIN", "Edge outlier minimum (m)", "double", 0.5, 0.0, adv=True)
        self.number("STEEP", "Steep ground, flow path not ponding (slope %)", "double", 10.0, 0.0, adv=True)
        self.number("CORR", "DEM error correlation length (m)", "double", 150.0, 0.0, adv=True)
        self.number("ERO_SIGMA", "Extra terrain error in erosion / deposition cells (m)", "double", 1.0, 0.0, adv=True)
        self.number("RNG_SEED", "Random seed", "int", 42, 0, adv=True)
        self.addParameter(self.advanced(QgsProcessingParameterEnum(
            "DEPTH_RAMP", "Map colours: depth ramp", ["Blues", "Viridis", "Yellow-green-blue", "Magma"],
            defaultValue=0)))
        self.number("DEPTH_MAX", "Map colours: depth at the dark end of the ramp (m)", "double", 3.0, 0.1, 100.0,
                    adv=True)
        for key, label, default in (("COL_PERM", "permanent water", "#7850a0"),
                                    ("COL_EDGE", "flood edge", "#08254e"),
                                    ("COL_DOUBT", "doubtful cells (ground above the water surface)", "#e78ac3")):
            self.addParameter(self.advanced(QgsProcessingParameterColor(
                key, f"Map colours: {label}", defaultValue=QColor(default), opacityEnabled=False)))
        self.addParameter(self.advanced(QgsProcessingParameterBoolean(
            "LOAD_ALL", "Load every output (otherwise depth, extent and evidence)", defaultValue=False)))
        self.addParameter(QgsProcessingParameterFolderDestination("OUT_DIR", "Output folder"))
        self.addOutput(QgsProcessingOutputRasterLayer("DEPTH", "Flood depth"))
        self.addOutput(QgsProcessingOutputVectorLayer("FLOOD_EXTENT", "Flood extent", Qgis.ProcessingSourceType.VectorPolygon))
        self.addOutput(QgsProcessingOutputFile("REPORT", "Run report (HTML)"))

    def options(self, parameters, context, feedback):
        """QuickOptions from the dialog (shared with the time-series tool)."""
        from ..core.pipeline import QuickOptions
        images = self.raster_paths(parameters, "IMAGES", context)
        ext_wkts, ext_crs = self.vector_wkts(parameters, "EXTENT", context, feedback)
        if not images and not ext_wkts:
            raise QgsProcessingException("Give the satellite image file(s) or a flood-extent polygon layer.")
        if images and ext_wkts:
            feedback.pushWarning("Both an image and a flood-extent layer were given: the image is used.")
            ext_wkts = []
        mode = MODES[self.parameterAsEnum(parameters, "MODE", context)]
        method = METHODS[self.parameterAsEnum(parameters, "METHOD", context)]
        thr = parameters.get("THRESHOLD")
        threshold = None if thr in (None, "") else self.parameterAsDouble(parameters, "THRESHOLD", context)
        overrides = {}
        for key, role in (("B_GREEN", "green"), ("B_NIR", "nir"), ("B_SWIR1", "swir1"), ("B_RED", "red"),
                          ("B_BLUE", "blue")):
            v = self.parameterAsInt(parameters, key, context)
            if v:
                overrides[role] = v
        aoi, aoi_crs = self.vector_wkts(parameters, "AOI", context, feedback)
        seeds, seed_crs = self.vector_wkts(parameters, "SEEDS", context, feedback)
        perm, perm_crs = self.vector_wkts(parameters, "PERM_POLYS", context, feedback)
        out_dir = self.parameterAsString(parameters, "OUT_DIR", context)
        opt = QuickOptions(
            dem_path=self.raster_path(parameters, "DEM", context), out_dir=out_dir,
            dem_tiles=self.raster_paths(parameters, "DEM_TILES", context),
            image_paths=images, pre_paths=self.raster_paths(parameters, "PRE_IMAGES", context), mode=mode,
            index=INDEXES[self.parameterAsEnum(parameters, "INDEX", context)], threshold=threshold,
            threshold_method=["edge", "tile", "global"][self.parameterAsEnum(parameters, "TMETHOD", context)],
            band_overrides=overrides, extent_wkts=ext_wkts, extent_crs=ext_crs, aoi_wkts=aoi, aoi_crs=aoi_crs,
            seed_wkts=seeds, seed_crs=seed_crs,
            seed_buffer_m=self.parameterAsDouble(parameters, "SEED_BUFFER", context),
            permanent_wkts=perm, permanent_crs=perm_crs,
            permanent_raster=self.raster_path(parameters, "PERM_RASTER", context, required=False) or "",
            permanent_values=self.parameterAsString(parameters, "PERM_VALUES", context) or "auto",
            method=method, hand_path=self.raster_path(parameters, "HAND", context, required=False) or "",
            cell=self.parameterAsDouble(parameters, "CELL", context),
            min_water_cells=self.parameterAsInt(parameters, "MIN_WATER", context),
            max_hole_cells=self.parameterAsInt(parameters, "MAX_HOLE", context),
            cloud_threshold=self.parameterAsDouble(parameters, "CLOUD", context),
            outlier_radius=self.parameterAsInt(parameters, "OUT_RADIUS", context),
            outlier_k=self.parameterAsDouble(parameters, "OUT_K", context),
            outlier_min=self.parameterAsDouble(parameters, "OUT_MIN", context),
            smooth_radius=self.parameterAsInt(parameters, "SMOOTH", context),
            smooth_tolerance=self.parameterAsDouble(parameters, "SMOOTH_TOL", context),
            edge_step=self.parameterAsDouble(parameters, "EDGE_STEP", context),
            edge_step_cap=self.parameterAsDouble(parameters, "EDGE_CAP", context),
            steep_slope_pct=self.parameterAsDouble(parameters, "STEEP", context),
            interpolation=INTERP[self.parameterAsEnum(parameters, "INTERP", context)],
            uncertainty_members=self.parameterAsInt(parameters, "MEMBERS", context),
            dem_sigma=self.parameterAsDouble(parameters, "DEM_SIGMA", context),
            dem_corr_m=self.parameterAsDouble(parameters, "CORR", context),
            erosion_sigma=self.parameterAsDouble(parameters, "ERO_SIGMA", context),
            seed=self.parameterAsInt(parameters, "RNG_SEED", context),
            map_style=self.map_style(parameters, context))
        if method == "hand" and not opt.hand_path:
            raise QgsProcessingException("The HAND method needs a HAND raster.")
        if mode == "trace" and not opt.pre_paths:
            raise QgsProcessingException("Flood-trace mode needs the pre-event image file(s).")
        return opt

    def map_style(self, parameters, context):
        st = {"depth_ramp": RAMPS[self.parameterAsEnum(parameters, "DEPTH_RAMP", context)],
              "depth_max": self.parameterAsDouble(parameters, "DEPTH_MAX", context)}
        for key, name in (("COL_PERM", "permanent"), ("COL_EDGE", "edge"), ("COL_DOUBT", "doubtful")):
            c = self.parameterAsColor(parameters, key, context)
            if c is not None and c.isValid():
                st[name] = c.name()
        return st

    def depth_style(self, opt):
        from ..core.report import MapStyle
        return RasterStyle("depth", MapStyle.from_dict(opt.map_style).depth_stops())

    def processAlgorithm(self, parameters, context, feedback):
        from ..core.pipeline import run_quick
        from ..core.raster import RasterError
        opt = self.options(parameters, context, feedback)
        out_dir = opt.out_dir
        try:
            info = run_quick(opt, progress=self.progress(feedback), log=self.log(feedback))
        except QgsProcessingException:
            raise
        except (ValueError, RuntimeError, RasterError, MemoryError) as exc:
            raise QgsProcessingException(str(exc) or type(exc).__name__)
        for w in info.get("warnings", []):
            feedback.pushWarning(w)
        o = info["outputs"]
        load_all = self.parameterAsBoolean(parameters, "LOAD_ALL", context)
        self.load_raster(context, o["depth"], "Flood depth (m)", self.depth_style(opt))
        self.load_raster(context, o["evidence"], "Flood evidence", RasterStyle("palette"))
        self.load_vector(context, o["extent"] + "|layername=flood_extent", "Flood extent")
        if load_all:
            for key, name, style in (("wse", "Water surface", None), ("p10", "Depth P10 (m)", "depth"),
                                     ("p50", "Depth P50 (m)", "depth"), ("p90", "Depth P90 (m)", "depth"),
                                     ("prob", "Flood probability", "prob"), ("water_mask", "Water mask", "palette")):
                if o.get(key):
                    rs = self.depth_style(opt) if style == "depth" else RasterStyle(style) if style else None
                    self.load_raster(context, o[key], name, rs)
        d = info.get("depth", {})
        feedback.pushInfo(f"Flooded {d.get('area_km2', 0):.3f} km2, mean depth {d.get('mean_m', float('nan')):.2f} m, "
                          f"P90 {d.get('p90_m', float('nan')):.2f} m. Report: {o['report']}")
        return {"OUT_DIR": out_dir, "DEPTH": o["depth"], "FLOOD_EXTENT": o["extent"] + "|layername=flood_extent",
                "REPORT": o["report"]}
