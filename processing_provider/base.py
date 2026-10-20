# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Shared base for FloodDepth Processing algorithms.

This is the ONLY place where QGIS and `flooddepth.core` meet. The core
never imports QGIS; the algorithms never implement the method. Layers are
handed to the core as file paths (rasters) or WKT + CRS WKT (vectors).
"""

import os

from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterNumber,
)


def _number_type(kind):
    """Integer / Double parameter type on QGIS 3.28 - 4."""
    try:
        from qgis.core import Qgis
        t = Qgis.ProcessingNumberParameterType
        return t.Integer if kind == "int" else t.Double
    except AttributeError:
        t = QgsProcessingParameterNumber
        return t.Integer if kind == "int" else t.Double


def _crs_wkt(crs):
    try:
        from qgis.core import QgsCoordinateReferenceSystem
        return crs.toWkt(QgsCoordinateReferenceSystem.WKT_PREFERRED_GDAL)
    except (AttributeError, TypeError):
        return crs.toWkt()


class FloodDepthAlgorithm(QgsProcessingAlgorithm):
    """Common metadata and helpers."""

    def group(self):
        return "Flood depth"

    def groupId(self):
        return "flooddepth"

    def createInstance(self):
        return type(self)()

    def helpUrl(self):
        return "https://github.com/EdmondAkello/Flood-Depth/blob/main/DOCUMENTATION.md"

    # -- parameters ---------------------------------------------------

    def advanced(self, param):
        """Mark a parameter as advanced (QGIS 3.28 - 4)."""
        try:
            from qgis.core import Qgis
            flag = Qgis.ProcessingParameterFlag.Advanced
        except AttributeError:
            from qgis.core import QgsProcessingParameterDefinition
            flag = QgsProcessingParameterDefinition.FlagAdvanced
        param.setFlags(param.flags() | flag)
        return param

    def number(self, key, label, kind="double", default=0.0, minimum=None, maximum=None,
               optional=False, adv=False):
        p = QgsProcessingParameterNumber(key, label, type=_number_type(kind), defaultValue=default,
                                         optional=optional)
        if minimum is not None:
            p.setMinimum(minimum)
        if maximum is not None:
            p.setMaximum(maximum)
        self.addParameter(self.advanced(p) if adv else p)
        return p

    # -- inputs ---------------------------------------------------------

    def raster_path(self, parameters, name, context, required=True):
        """A raster parameter as a file path GDAL can open (None if absent)."""
        layer = self.parameterAsRasterLayer(parameters, name, context)
        if layer is None:
            if required:
                raise QgsProcessingException(f"Could not resolve raster parameter '{name}'.")
            return None
        path = layer.source().split("|")[0]
        if not os.path.exists(path):
            raise QgsProcessingException(
                f"Layer '{layer.name()}' is not a file GDAL can open directly ({path}). "
                "Save it to GeoTIFF first.")
        return path

    def raster_paths(self, parameters, name, context):
        """A multiple-layer parameter as a list of file paths."""
        out = []
        for layer in self.parameterAsLayerList(parameters, name, context) or []:
            path = layer.source().split("|")[0]
            if not os.path.exists(path):
                raise QgsProcessingException(
                    f"Layer '{layer.name()}' is not a file on disk ({path}). Save it to GeoTIFF first.")
            out.append(path)
        return out

    def vector_wkts(self, parameters, name, context, feedback=None):
        """(list of WKT, CRS WKT) of a vector parameter; ([], '') if absent."""
        source = self.parameterAsSource(parameters, name, context)
        if source is None:
            return [], ""
        wkts = []
        for f in source.getFeatures():
            g = f.geometry()
            if g is None or g.isEmpty():
                continue
            wkts.append(g.asWkt())
        if not wkts and feedback is not None:
            feedback.pushWarning(f"Layer '{name}' has no geometry; ignored.")
        return wkts, _crs_wkt(source.sourceCrs())

    def progress(self, feedback):
        """Adapt core's progress(fraction, message) to QGIS feedback."""
        def _p(fraction, message=""):
            if feedback.isCanceled():
                raise QgsProcessingException("Cancelled by user.")
            feedback.setProgress(int(max(0.0, min(1.0, fraction)) * 100))
            if message:
                feedback.setProgressText(message)
        return _p

    def log(self, feedback):
        return lambda msg: feedback.pushInfo(msg)

    # -- outputs ------------------------------------------------------

    def load_raster(self, context, path, name, style=None):
        """Add a raster to the project when the algorithm finishes, styled."""
        from qgis.core import QgsProcessingContext
        details = QgsProcessingContext.LayerDetails(name, context.project(), name)
        if style is not None:
            details.setPostProcessor(style.instance)
        context.addLayerToLoadOnCompletion(path, details)

    def load_vector(self, context, uri, name):
        from qgis.core import QgsProcessingContext
        context.addLayerToLoadOnCompletion(uri, QgsProcessingContext.LayerDetails(name, context.project(), name))


class RasterStyle(object):
    """Post-processor that styles a raster on load: depth ramp, probability
    ramp or a palette (classes). Cosmetic: failures never fail a run."""

    DEPTH = [(0.0, "#c6dbef", "0"), (0.5, "#6baed6", "0.5"), (1.0, "#3182bd", "1"),
             (2.0, "#08519c", "2"), (3.0, "#08306b", "3+")]
    PROB = [(0.0, "#fff5eb", "0"), (0.5, "#fd8d3c", "0.5"), (1.0, "#7f2704", "1")]

    _keep = []   # post-processors must outlive the algorithm call

    def __init__(self, kind):
        from qgis.core import QgsProcessingLayerPostProcessorInterface
        stops = {"depth": self.DEPTH, "prob": self.PROB}.get(kind)

        class _P(QgsProcessingLayerPostProcessorInterface):
            def postProcessLayer(self, layer, context, feedback):
                try:
                    if stops is None:
                        return  # palette rasters carry their own colour table
                    from qgis.core import (QgsColorRampShader, QgsRasterShader,
                                           QgsSingleBandPseudoColorRenderer)
                    from qgis.PyQt.QtGui import QColor
                    shader = QgsColorRampShader(stops[0][0], stops[-1][0])
                    try:
                        shader.setColorRampType(QgsColorRampShader.Interpolated)
                    except AttributeError:
                        from qgis.core import Qgis
                        shader.setColorRampType(Qgis.ShaderInterpolationMethod.Linear)
                    shader.setColorRampItemList([QgsColorRampShader.ColorRampItem(v, QColor(c), lab)
                                                 for v, c, lab in stops])
                    rs = QgsRasterShader()
                    rs.setRasterShaderFunction(shader)
                    r = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, rs)
                    r.setClassificationMin(stops[0][0])
                    r.setClassificationMax(stops[-1][0])
                    layer.setRenderer(r)
                    layer.triggerRepaint()
                except Exception as exc:   # cosmetic only
                    feedback.pushInfo(f"(style not applied: {exc})")
        self.instance = _P()
        RasterStyle._keep.append(self)
