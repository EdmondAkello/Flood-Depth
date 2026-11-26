# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Export 3D scene: an interactive WebGL HTML scene from a finished run."""

import os

from qgis.core import (
    Qgis,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterFile,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterMultipleLayers,
)

from .base import FloodDepthAlgorithm, _crs_wkt


class Scene3dAlgorithm(FloodDepthAlgorithm):

    def name(self):
        return "scene3d"

    def displayName(self):
        return "Export 3D scene (interactive HTML)"

    def shortHelpString(self):
        return (
            "An interactive 3D scene of a finished <b>Quick flood depth</b> run, as one HTML file that opens "
            "offline in any current browser (WebGL 2).\n\n"
            "<b>Shows</b> the terrain draped with the satellite image (or a hillshade), permanent water and "
            "doubtful cells painted on it, the water surface coloured by depth (the run's map colours), the flood "
            "edge, and any line or polygon layers you add (roads, rivers, the area of interest).\n\n"
            "<b>In the viewer</b>: drag to rotate, right-drag to pan, wheel to zoom; switch layers and their "
            "opacity; change the vertical exaggeration; click the water to read its level and depth; saved views; "
            "PNG and A3 sheet export.\n\n"
            "Quick flood depth writes scene3d.html itself (on by default); use this tool to add layers, change "
            "the level of detail or the exaggeration without running again. The DEM is read again from the path "
            "recorded in run.json.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(
            "RUN_DIR", "Output folder of a Quick flood depth run (holds run.json)",
            behavior=Qgis.ProcessingFileParameterBehavior.Folder))
        self.addParameter(QgsProcessingParameterMultipleLayers(
            "LINES", "Line or polygon layers to draw on the terrain (roads, rivers, boundaries)",
            Qgis.ProcessingSourceType.Vector, optional=True))
        self.number("DETAIL", "Level of detail (vertices per surface, thousands)", "int", 600, 50, 4000)
        self.number("VE", "Vertical exaggeration (0 = automatic)", "double", 0.0, 0.0, 50.0)
        self.addParameter(self.advanced(QgsProcessingParameterBoolean(
            "OPEN", "Open the scene in the browser when done", defaultValue=False)))
        self.addParameter(QgsProcessingParameterFileDestination(
            "OUTPUT", "3D scene", "HTML files (*.html)", optional=True))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core import raster as R
        from ..core.scene3d import scene_from_folder
        run_dir = self.parameterAsFile(parameters, "RUN_DIR", context)
        if not os.path.exists(os.path.join(run_dir, "run.json")):
            raise QgsProcessingException("The folder has no run.json: choose the output folder of a Quick flood "
                                         "depth run.")
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context) or os.path.join(run_dir, "scene3d.html")
        if not out.lower().endswith(".html"):
            out += ".html"
        palette = ["#c2185b", "#e39a26", "#1f8f80", "#4a45a6", "#555555"]
        lines = []
        for i, layer in enumerate(self.parameterAsLayerList(parameters, "LINES", context) or []):
            if layer.type() != Qgis.LayerType.Vector:
                continue
            wkts = [f.geometry().asWkt() for f in layer.getFeatures() if not f.geometry().isEmpty()]
            if wkts:
                lines.append((layer.name(), layer.name(), wkts, _crs_wkt(layer.crs()), palette[i % len(palette)]))
        feedback.setProgressText("Building the 3D scene")
        try:
            info = scene_from_folder(run_dir, out, max_vertices=1000 * self.parameterAsInt(parameters, "DETAIL", context),
                                     extra_lines=lines,
                                     vertical_exaggeration=self.parameterAsDouble(parameters, "VE", context) or None)
        except (ValueError, RuntimeError, KeyError, OSError, R.RasterError) as exc:
            raise QgsProcessingException(f"The 3D scene could not be written: {exc}")
        feedback.pushInfo(f"3D scene: {info['path']} ({info['size_mb']} MB; terrain at {info['cell_m']:g} m, "
                          f"{info['terrain_triangles']:,} + {info['water_triangles']:,} triangles; "
                          f"{info['drape']}; vertical exaggeration {info['vertical_exaggeration']:g}).")
        if self.parameterAsBoolean(parameters, "OPEN", context):
            from qgis.PyQt.QtCore import QUrl
            from qgis.PyQt.QtGui import QDesktopServices
            QDesktopServices.openUrl(QUrl.fromLocalFile(info["path"]))
        return {"OUTPUT": info["path"]}
