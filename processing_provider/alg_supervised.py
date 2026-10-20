# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Supervised water classification from training polygons (RGB or any image)."""

import numpy as np
from qgis.core import (
    QgsProcessing,
    QgsProcessingException,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
    QgsProcessingParameterVectorDestination,
)

from .base import FloodDepthAlgorithm


class SupervisedWaterAlgorithm(FloodDepthAlgorithm):

    def name(self):
        return "supervisedwater"

    def displayName(self):
        return "1b. Classify water from training polygons (RGB, drone, aerial)"

    def group(self):
        return "Step by step"

    def groupId(self):
        return "steps"

    def shortHelpString(self):
        return ("For images without NIR or SWIR (true-colour satellite, drone or aerial photos) or scenes "
                "where the automatic indices fail. Draw a few polygons over water and over non-water (include "
                "shadows, dark roofs and wet soil in the non-water set). Each pixel goes to the nearer class "
                "centroid in standardised band space. Review the result: two classes cannot separate every "
                "surface that looks like water.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer("IMAGE", "Image"))
        self.addParameter(QgsProcessingParameterString("BANDS", "Bands to use (comma-separated, empty = all data bands)",
                                                       defaultValue="", optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "WATER_TRAIN", "Water training polygons", [QgsProcessing.TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "LAND_TRAIN", "Non-water training polygons", [QgsProcessing.TypeVectorPolygon]))
        self.number("MIN_WATER", "Remove water patches smaller than (cells)", "int", 10, 0)
        self.addParameter(QgsProcessingParameterRasterDestination("MASK", "Water mask"))
        self.addParameter(QgsProcessingParameterVectorDestination(
            "POLYGONS", "Water polygons", QgsProcessing.TypeVectorPolygon))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core import raster as R, vector as V, classify as C, imagery as I
        path = self.raster_path(parameters, "IMAGE", context)
        ds = R.open_dataset(path)
        alpha = R.alpha_band_index(ds)
        txt = (self.parameterAsString(parameters, "BANDS", context) or "").strip()
        try:
            bands = [int(x) for x in txt.replace(";", ",").split(",") if x.strip()] if txt else \
                [i for i in range(1, ds.RasterCount + 1) if i != alpha]
        except ValueError:
            raise QgsProcessingException("Bands must be comma-separated band numbers, e.g. 1,2,3.")
        if not bands or any(b < 1 or b > ds.RasterCount for b in bands):
            raise QgsProcessingException(f"The image has {ds.RasterCount} band(s); check the band list.")
        stack, grid = [], None
        for b in bands:
            a, g, _ = R.read_band(path, b, ds)
            stack.append(a)
            grid = g
        stack = np.stack(stack)
        ww, wc = self.vector_wkts(parameters, "WATER_TRAIN", context, feedback)
        lw, lc = self.vector_wkts(parameters, "LAND_TRAIN", context, feedback)
        wm = V.rasterize_wkts(ww, wc, grid, all_touched=True)
        lm = V.rasterize_wkts(lw, lc, grid, all_touched=True)
        try:
            model = C.fit_centroids(stack[:, wm].T, stack[:, lm].T)
        except ValueError as exc:
            raise QgsProcessingException(f"{exc} Draw larger training polygons inside the image.")
        water, valid = C.classify_centroids(stack, model)
        water, cs = C.clean_mask(water, valid, self.parameterAsInt(parameters, "MIN_WATER", context), 0)
        feedback.pushInfo(f"Training pixels: water {int(wm.sum()):,}, non-water {int(lm.sum()):,}. "
                          f"Water {100 * water.sum() / max(1, valid.sum()):.2f} % of valid pixels; cleanup {cs}.")
        out = self.parameterAsOutputLayer(parameters, "MASK", context)
        R.write_raster(out, I.mask_raster_values(water, valid), grid, 255, "Byte",
                       {"VALUES": "1 water, 0 dry, 255 NoData", "METHOD": "nearest centroid"},
                       color_table=I.MASK_COLORS)
        vout = self.parameterAsOutputLayer(parameters, "POLYGONS", context)
        import os
        if vout.lower().endswith(".gpkg") and os.path.exists(vout):
            os.remove(vout)
        V.polygonize(water, grid, vout, "water", {"method": "supervised"})
        return {"MASK": out, "POLYGONS": vout}
