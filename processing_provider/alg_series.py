# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Flood depth for several image dates (time series)."""

from qgis.core import (
    QgsProcessingException,
    QgsProcessingOutputFile,
    QgsProcessingOutputRasterLayer,
)

from .alg_quick import QuickFloodDepthAlgorithm
from .base import RasterStyle


class FloodSeriesAlgorithm(QuickFloodDepthAlgorithm):
    """Same inputs and settings as Quick flood depth; the image files may
    hold several dates, which are grouped by the date in their names."""

    def name(self):
        return "flooddepthseries"

    def displayName(self):
        return "Flood depth over time (several image dates)"

    def shortHelpString(self):
        return (
            "Runs <b>Quick flood depth</b> once per image date and combines the results: how the flooded "
            "area, depth and volume changed between overpasses, the <b>maximum depth seen on any date</b>, "
            "how many dates each cell was flooded, and the first date it flooded.\n\n"
            "<b>Image files.</b> Select the band files of all dates together. They are grouped by the date in "
            "their names (YYYY-MM-DD as in Copernicus Browser exports, or SAFE / Landsat names). Every date "
            "uses the same DEM, permanent water, area of interest and settings.\n\n"
            "<b>Comparable dates.</b> Give an area of interest so every date covers the same ground; cloud and "
            "image footprints otherwise change the mapped area. A date with no flood left (or a refused "
            "threshold) is reported and skipped; the others still run.\n\n"
            "Outputs (output folder): one sub-folder per date with the usual outputs and report, plus "
            "max_depth.tif, dates_flooded.tif, first_date_flooded.tif, series.csv and series.html.\n\n"
            "Tip: 10 uncertainty members per date keep long series quick.")

    def initAlgorithm(self, config=None):
        super().initAlgorithm(config)
        self.parameterDefinition("IMAGES").setDescription("Satellite image files of all dates (post-event)")
        self.parameterDefinition("MEMBERS").setDefaultValue(10)
        self.removeParameter("EXTENT")
        self.addOutput(QgsProcessingOutputRasterLayer("MAX_DEPTH", "Maximum observed depth"))
        self.addOutput(QgsProcessingOutputFile("SUMMARY", "Series summary (HTML)"))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core.series import run_series
        from ..core import raster as R
        if not self.raster_paths(parameters, "IMAGES", context):
            raise QgsProcessingException("Select the image files of the dates to map.")
        opt = self.options(parameters, context, feedback)
        try:
            res = run_series(opt, progress=self.progress(feedback), log=self.log(feedback))
        except QgsProcessingException:
            raise
        except (ValueError, RuntimeError, R.RasterError) as exc:
            raise QgsProcessingException(str(exc))
        for r in res["rows"]:
            if r["status"] != "ok":
                feedback.pushWarning(f"{r['date']}: {r['status']}")
        o = res["outputs"]
        self.load_raster(context, o["max_depth"], "Maximum observed depth (m)", RasterStyle("depth"))
        self.load_raster(context, o["dates_flooded"], "Dates flooded", None)
        e = res["envelope"]
        feedback.pushInfo(f"Flooded on any date {e['area_km2']:.3f} km2; on every date "
                          f"{e['flooded_every_date_km2']:.3f} km2. Summary: {o['html']}")
        return {"OUT_DIR": opt.out_dir, "MAX_DEPTH": o["max_depth"], "SUMMARY": o["html"]}
