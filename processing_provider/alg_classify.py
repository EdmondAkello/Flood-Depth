# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Step 1: classify water (or flood trace) from satellite imagery."""

import os

from qgis.core import (
    QgsProcessing,
    QgsProcessingException,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterMultipleLayers,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterVectorDestination,
)

from .base import FloodDepthAlgorithm

MODES = ["water", "trace"]
INDEXES = ["auto", "mndwi", "ndwi", "sar"]


class ClassifyWaterAlgorithm(FloodDepthAlgorithm):

    def name(self):
        # id kept from 0.3 so existing models keep resolving
        return "classify_inundation"

    def displayName(self):
        return "1. Classify water from imagery (automatic)"

    def group(self):
        return "Step by step"

    def groupId(self):
        return "steps"

    def shortHelpString(self):
        return (
            "Water mask from satellite image file(s) with no training data. Bands are recognised from file "
            "names and metadata (Browser 'Raw' exports, Sentinel-2 SAFE, Landsat, multiband stacks with band "
            "names); set band numbers in the advanced parameters only for unnamed stacks.\n\n"
            "Index: MNDWI when green and SWIR1 are present, else NDWI; Sentinel-1 VV/VH in dB or linear. "
            "Threshold: tile-based Otsu (Martinis et al. 2009) with the separability eta reported; an automatic "
            "threshold outside the physically usual range is refused with an explanation rather than used.\n\n"
            "Mask values: 1 flood water, 2 permanent water (seen in the pre-event image), 0 dry, 255 NoData. "
            "Cleanup removes water patches below a size and, only if asked, fills small dry holes - islands, "
            "roads and buildings are kept by default.\n\n"
            "Flood-trace mode maps vegetation lost between a pre- and a post-event image (NDVI drop) for "
            "floods that drained before the overpass.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterMultipleLayers(
            "IMAGES", "Satellite image file(s)", QgsProcessing.TypeRaster))
        self.addParameter(QgsProcessingParameterEnum(
            "MODE", "What the image shows", ["Water at the overpass", "Flood trace (needs a pre-event image)"],
            defaultValue=0))
        self.addParameter(QgsProcessingParameterMultipleLayers(
            "PRE_IMAGES", "Pre-event image file(s)", QgsProcessing.TypeRaster, optional=True))
        self.addParameter(QgsProcessingParameterEnum(
            "INDEX", "Water index", ["Automatic", "MNDWI (green, SWIR1)", "NDWI (green, NIR)",
                                     "SAR backscatter (VV / VH)"], defaultValue=0))
        self.number("THRESHOLD", "Manual threshold (empty = automatic)", "double", None, optional=True)
        self.addParameter(QgsProcessingParameterEnum(
            "TMETHOD", "Automatic threshold", ["Tile-based Otsu", "Global Otsu"], defaultValue=0))
        self.number("MIN_WATER", "Remove water patches smaller than (cells)", "int", 10, 0)
        self.number("MAX_HOLE", "Fill dry holes up to (cells; 0 = keep islands)", "int", 0, 0)
        self.number("CLOUD", "Bright-cloud blue reflectance limit (0 = off)", "double", 0.25, 0.0, 1.0, adv=True)
        for key, label in (("B_GREEN", "green"), ("B_NIR", "NIR"), ("B_SWIR1", "SWIR1"), ("B_RED", "red"),
                           ("B_BLUE", "blue")):
            self.number(key, f"Band number of {label} (0 = automatic)", "int", 0, 0, 500, adv=True)
        self.addParameter(QgsProcessingParameterFeatureSource(
            "SEEDS", "Keep only water connected to these features", [QgsProcessing.TypeVectorAnyGeometry],
            optional=True))
        self.number("SEED_BUFFER", "Connection distance (map units of the image)", "double", 0.0, 0.0, adv=True)
        self.addParameter(QgsProcessingParameterRasterDestination("MASK", "Water mask"))
        self.addParameter(QgsProcessingParameterVectorDestination(
            "POLYGONS", "Water polygons", QgsProcessing.TypeVectorPolygon))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core import imagery as I, raster as R, vector as V, classify as C
        images = self.raster_paths(parameters, "IMAGES", context)
        if not images:
            raise QgsProcessingException("Select at least one image file.")
        thr = parameters.get("THRESHOLD")
        threshold = None if thr in (None, "") else self.parameterAsDouble(parameters, "THRESHOLD", context)
        overrides = {role: self.parameterAsInt(parameters, key, context)
                     for key, role in (("B_GREEN", "green"), ("B_NIR", "nir"), ("B_SWIR1", "swir1"),
                                       ("B_RED", "red"), ("B_BLUE", "blue"))}
        overrides = {k: v for k, v in overrides.items() if v}
        mode = MODES[self.parameterAsEnum(parameters, "MODE", context)]
        pre = self.raster_paths(parameters, "PRE_IMAGES", context)
        if mode == "trace" and not pre:
            raise QgsProcessingException("Flood-trace mode needs the pre-event image file(s).")
        try:
            res = I.classify_image(
                images, overrides, INDEXES[self.parameterAsEnum(parameters, "INDEX", context)], threshold,
                ["tile", "global"][self.parameterAsEnum(parameters, "TMETHOD", context)],
                self.parameterAsInt(parameters, "MIN_WATER", context),
                self.parameterAsInt(parameters, "MAX_HOLE", context),
                self.parameterAsDouble(parameters, "CLOUD", context), pre or None, None, mode,
                progress=self.progress(feedback))
        except (ValueError, R.RasterError) as exc:
            raise QgsProcessingException(str(exc))
        info = res["info"]
        water = res["water"]
        seeds, seed_crs = self.vector_wkts(parameters, "SEEDS", context, feedback)
        if seeds:
            g = res["grid"]
            buf = self.parameterAsDouble(parameters, "SEED_BUFFER", context)
            sm = V.rasterize_wkts(seeds, seed_crs, g, all_touched=True)
            sm = C.dilate(sm, int(round(buf / g.cell_x)) if buf > 0 else 1)
            water, kept, dropped = C.keep_connected(water, sm)
            feedback.pushInfo(f"Connectivity: kept {kept} region(s), dropped {dropped}.")
        t = info["threshold"]
        feedback.pushInfo(f"Bands: {info['bands']}")
        feedback.pushInfo(f"Index: {info['index']}; threshold {t['threshold']:.4f} ({t['method']}"
                          + (f", eta {t['eta']:.2f} - {info.get('confidence')}" if t.get("eta") is not None else "")
                          + (f"; near-optimal band {info['threshold_band'][0]:.4f} to {info['threshold_band'][1]:.4f}"
                             if info.get("threshold_band") else "") + ")")
        feedback.pushInfo(f"Water: {100 * info['water_share']:.2f} % of valid pixels "
                          f"(before cleanup {100 * info['water_share_raw']:.2f} %); cleanup {info['cleanup']}")
        for n in info.get("notes", []):
            feedback.pushWarning(n)
        out = self.parameterAsOutputLayer(parameters, "MASK", context)
        R.write_raster(out, I.mask_raster_values(water, res["valid"], res["permanent"]), res["grid"], 255, "Byte",
                       tags={"VALUES": "1 flood water, 2 permanent water, 0 dry, 255 NoData",
                             "INDEX": info["index"], "THRESHOLD": t["threshold"]},
                       color_table=I.MASK_COLORS)
        vout = self.parameterAsOutputLayer(parameters, "POLYGONS", context)
        flood = water if res["permanent"] is None else (water & ~res["permanent"])
        if vout.lower().endswith(".gpkg") and os.path.exists(vout):
            os.remove(vout)
        n = V.polygonize(flood, res["grid"], vout, "water", {"index": info["index_key"],
                                                            "threshold": float(t["threshold"])})
        feedback.pushInfo(f"{n} polygon(s) written.")
        return {"MASK": out, "POLYGONS": vout}
