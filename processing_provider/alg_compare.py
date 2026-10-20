# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Validation: compare a flood extent with a reference extent."""

from qgis.core import (
    QgsProcessing,
    QgsProcessingException,
    QgsProcessingOutputNumber,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterRasterLayer,
)

from .base import FloodDepthAlgorithm


class CompareExtentsAlgorithm(FloodDepthAlgorithm):

    def name(self):
        return "compareextents"

    def displayName(self):
        return "Compare flood extents (CSI, hit rate, false alarms)"

    def group(self):
        return "Validation"

    def groupId(self):
        return "validation"

    def shortHelpString(self):
        return ("Scores a flood extent against a reference (field survey, Copernicus EMS / UNOSAT map, a "
                "hydraulic model, another image). Each may be polygons or a mask raster (cells equal to the "
                "water value are flooded; 255, negative values and NoData are not counted).\n\n"
                "Scores over cells valid in both: CSI = hits / (hits + false alarms + misses); hit rate = "
                "hits / (hits + misses); false-alarm ratio = false alarms / (hits + false alarms); bias = "
                "flooded area / reference area; F1. The agreement raster shows 1 hit, 2 false alarm, 3 miss, "
                "0 dry in both.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(
            "PRED_POLYS", "Flood extent to test (polygons)", [QgsProcessing.TypeVectorPolygon], optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer("PRED_RASTER", "...or as a mask raster", optional=True))
        self.number("PRED_VALUE", "Water value in that raster", "double", 1.0)
        self.addParameter(QgsProcessingParameterFeatureSource(
            "REF_POLYS", "Reference extent (polygons)", [QgsProcessing.TypeVectorPolygon], optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer("REF_RASTER", "...or as a mask raster", optional=True))
        self.number("REF_VALUE", "Water value in that raster", "double", 1.0)
        self.addParameter(QgsProcessingParameterFeatureSource(
            "AOI", "Area of interest", [QgsProcessing.TypeVectorPolygon], optional=True))
        self.number("CELL", "Comparison cell (m; 0 = a raster's own grid, else 10 m)", "double", 0.0, 0.0)
        self.addParameter(QgsProcessingParameterRasterDestination("AGREEMENT", "Agreement raster"))
        for k in ("CSI", "HIT_RATE", "FALSE_ALARM_RATIO", "BIAS", "F1"):
            self.addOutput(QgsProcessingOutputNumber(k, k.replace("_", " ").title()))

    def _src(self, parameters, context, feedback, poly, ras, val):
        path = self.raster_path(parameters, ras, context, required=False)
        if path:
            return {"raster": path, "value": self.parameterAsDouble(parameters, val, context)}
        w, c = self.vector_wkts(parameters, poly, context, feedback)
        return {"wkts": w, "crs": c} if w else None

    def processAlgorithm(self, parameters, context, feedback):
        from ..core.compare import compare_extents, AGREEMENT, AGREEMENT_COLORS
        from ..core import raster as R
        pred = self._src(parameters, context, feedback, "PRED_POLYS", "PRED_RASTER", "PRED_VALUE")
        ref = self._src(parameters, context, feedback, "REF_POLYS", "REF_RASTER", "REF_VALUE")
        if pred is None or ref is None:
            raise QgsProcessingException("Give the extent to test and the reference (polygons or raster each).")
        aoi, ac = self.vector_wkts(parameters, "AOI", context, feedback)
        try:
            s, ag, grid = compare_extents(pred, ref, aoi, ac, self.parameterAsDouble(parameters, "CELL", context))
        except (ValueError, R.RasterError) as exc:
            raise QgsProcessingException(str(exc))
        p = self.parameterAsOutputLayer(parameters, "AGREEMENT", context)
        R.write_raster(p, ag, grid, 255, "Byte", {"VALUES": "; ".join(f"{k} {v}" for k, v in AGREEMENT.items())},
                       color_table=AGREEMENT_COLORS)
        feedback.pushInfo(f"Hits {s['tp']:,}, false alarms {s['fp']:,}, misses {s['fn']:,}, dry in both {s['tn']:,}")
        feedback.pushInfo(f"CSI {s['csi']:.3f} | hit rate {s['hit_rate']:.3f} | false-alarm ratio "
                          f"{s['false_alarm_ratio']:.3f} | bias {s['bias']:.3f} | F1 {s['f1']:.3f}")
        if "pred_area_km2" in s:
            feedback.pushInfo(f"Area tested {s['pred_area_km2']:.3f} km2, reference {s['ref_area_km2']:.3f} km2")
        return {"AGREEMENT": p, "CSI": s["csi"], "HIT_RATE": s["hit_rate"],
                "FALSE_ALARM_RATIO": s["false_alarm_ratio"], "BIAS": s["bias"], "F1": s["f1"]}
