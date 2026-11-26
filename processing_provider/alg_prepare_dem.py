# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Prepare DEM: merge tiles, clip, reproject to UTM (bilinear) and audit."""

from qgis.core import (
    Qgis,
    QgsProcessingException,
    QgsProcessingParameterCrs,
    QgsProcessingParameterEnum,
    QgsProcessingParameterExtent,
    QgsProcessingParameterMultipleLayers,
    QgsProcessingParameterRasterDestination,
)

from .base import FloodDepthAlgorithm, _crs_wkt

RESAMPLING = ["bilinear", "cubic"]


class PrepareDemAlgorithm(FloodDepthAlgorithm):

    def name(self):
        return "preparedem"

    def displayName(self):
        return "Prepare DEM (merge tiles, clip, reproject)"

    def group(self):
        return "DEM preparation"

    def groupId(self):
        return "dempreparation"

    def shortHelpString(self):
        return (
            "Makes one analysis-ready DEM from downloaded tiles: merged, clipped to the area with a buffer, "
            "in the WGS 84 / UTM zone of the area (or a CRS you choose), float32, resampled with "
            "<b>bilinear</b> interpolation (nearest neighbour repeats rows and columns and puts steps in the "
            "terrain, so it is not offered).\n\n"
            "<b>Tiles.</b> Select every tile; tiles in other coordinate systems are ignored (the Google GEM 15 m "
            "DTM for Kenya ships the country in four UTM zones: the zone best covering the area is used).\n\n"
            "<b>Check.</b> The output is audited for repeated rows / columns (nearest-neighbour resampling "
            "upstream) and for elevations stored in whole metres on flat ground; findings are listed in the "
            "log and in a .json file next to the output.\n\n"
            "No depression filling is done: filling raises the ground in the sinks where flood water stands, and "
            "the depth is water surface minus the real ground.\n\n"
            "Quick flood depth reads tiles directly too (DEM tiles parameter): use this tool when the same area "
            "is analysed many times, or to share a prepared DEM.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterMultipleLayers(
            "INPUT", "DEM file(s) or tiles", Qgis.ProcessingSourceType.Raster))
        self.addParameter(QgsProcessingParameterExtent("CLIP", "Clip to this extent", optional=True))
        self.number("BUFFER", "Buffer around the extent (m)", "double", 1000.0, 0.0, 100000.0)
        self.number("CELL", "Cell size (m; 0 = the DEM's native resolution)", "double", 0.0, 0.0, 1000.0)
        self.addParameter(self.advanced(QgsProcessingParameterCrs(
            "TARGET_CRS", "Target CRS (empty = WGS 84 / UTM zone of the area)", optional=True)))
        self.addParameter(self.advanced(QgsProcessingParameterEnum(
            "RESAMPLING", "Resampling", ["Bilinear", "Cubic"], defaultValue=0)))
        self.addParameter(QgsProcessingParameterRasterDestination("OUTPUT", "Prepared DEM"))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core import dem as DM, raster as R
        paths = self.raster_paths(parameters, "INPUT", context)
        if not paths:
            raise QgsProcessingException("Select the DEM file(s).")
        clip, clip_wkt = None, None
        if parameters.get("CLIP"):
            rect = self.parameterAsExtent(parameters, "CLIP", context)
            if rect is not None and not rect.isEmpty():
                crs = self.parameterAsExtentCrs(parameters, "CLIP", context)
                clip = (rect.xMinimum(), rect.yMinimum(), rect.xMaximum(), rect.yMaximum())
                clip_wkt = _crs_wkt(crs) if crs.isValid() else None
        target = None
        tcrs = self.parameterAsCrs(parameters, "TARGET_CRS", context)
        if tcrs is not None and tcrs.isValid():
            target = _crs_wkt(tcrs)
        out = self.parameterAsOutputLayer(parameters, "OUTPUT", context)
        if not out.lower().endswith((".tif", ".tiff")):
            out += ".tif"
        try:
            info = DM.prepare_dem(paths, out, clip, clip_wkt, self.parameterAsDouble(parameters, "BUFFER", context),
                                  self.parameterAsDouble(parameters, "CELL", context), target,
                                  RESAMPLING[self.parameterAsEnum(parameters, "RESAMPLING", context)],
                                  progress=self.progress(feedback))
        except (ValueError, RuntimeError, R.RasterError) as exc:
            raise QgsProcessingException(str(exc))
        for n in info["notes"]:
            feedback.pushInfo(n)
        for w in info["warnings"]:
            feedback.pushWarning(w)
        feedback.pushInfo(f"Prepared DEM: {info['cols']:,} x {info['rows']:,} cells of {info['cell_m']:g} m, "
                          f"{info['crs']}; {100 * info['valid_share']:.1f} % valid.")
        return {"OUTPUT": out}
