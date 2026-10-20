# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Step 3: depth from a known water surface (raster or one level)."""

import numpy as np
from qgis.core import (
    QgsProcessing,
    QgsProcessingException,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterRasterLayer,
)

from .base import FloodDepthAlgorithm


def _dem_grid_and_extent(alg, parameters, context, feedback, dem_path):
    from ..core import raster as R
    wkts, wcrs = alg.vector_wkts(parameters, "EXTENT", context, feedback)
    if wkts:
        from ..core.pipeline import _wkts_bounds
        ext, ew = _wkts_bounds(wkts, wcrs, pad_frac=0.05)
        grid, notes = R.analysis_grid(dem_path, ext, ew)
    else:
        grid, notes = R.analysis_grid(dem_path)
    for n in notes:
        feedback.pushInfo(n)
    return grid, wkts, wcrs


class DepthFromSurfaceAlgorithm(FloodDepthAlgorithm):

    def name(self):
        return "compute_flood_depth"   # id kept from 0.3

    def displayName(self):
        return "3. Flood depth from a water-surface raster"

    def group(self):
        return "Step by step"

    def groupId(self):
        return "steps"

    def shortHelpString(self):
        return ("depth = water-surface elevation - DEM, on the DEM grid (a geographic DEM is processed on the "
                "WGS 84 / UTM grid). The water-surface raster may be on any grid and CRS: it is warped "
                "(bilinear). Optionally limited to a flood-extent layer. Cells where the ground is above the "
                "surface are NoData. The two rasters must share a vertical datum and units.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", "DEM"))
        self.addParameter(QgsProcessingParameterRasterLayer("STAGE", "Water-surface elevation raster"))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "EXTENT", "Flood extent (optional)", [QgsProcessing.TypeVectorPolygon], optional=True))
        self.addParameter(QgsProcessingParameterRasterDestination("OUTPUT", "Flood depth"))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core import raster as R, vector as V
        dem_path = self.raster_path(parameters, "DEM", context)
        try:
            grid, wkts, wcrs = _dem_grid_and_extent(self, parameters, context, feedback, dem_path)
            dem, dn = R.read_on_grid(dem_path, grid, 1, "bilinear", zero_nodata_auto=True)
            wse, _ = R.read_on_grid(self.raster_path(parameters, "STAGE", context), grid, 1, "bilinear")
        except R.RasterError as exc:
            raise QgsProcessingException(str(exc))
        d = wse - dem
        ok = np.isfinite(d) & (d > 0)
        if wkts:
            ok &= V.rasterize_wkts(wkts, wcrs, grid)
        out = np.where(ok, d, np.nan)
        p = self.parameterAsOutputLayer(parameters, "OUTPUT", context)
        R.write_raster(p, out, grid, -9999.0, "Float32", {"CONTENT": "flood depth (m)"})
        if ok.any():
            feedback.pushInfo(f"{int(ok.sum()):,} wet cells; mean depth {float(np.nanmean(out)):.2f} m, "
                              f"max {float(np.nanmax(out)):.2f} m.")
        else:
            feedback.pushWarning("No cell has the water surface above the ground.")
        return {"OUTPUT": p}


class DepthConstantLevelAlgorithm(FloodDepthAlgorithm):

    def name(self):
        return "compute_constant_stage_depth"   # id kept from 0.3

    def displayName(self):
        return "Flood depth below a constant water level"

    def group(self):
        return "Step by step"

    def groupId(self):
        return "steps"

    def shortHelpString(self):
        return ("Bathtub depth: depth = level - DEM where positive, optionally limited to a flood extent. "
                "Use for a lake or reservoir level, a surveyed flood mark, or a what-if level. The level is in "
                "the DEM's vertical datum and units. Without an extent every cell below the level is wet, "
                "connected or not.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", "DEM"))
        self.number("LEVEL", "Water level (DEM vertical units)", "double", 0.0)
        self.addParameter(QgsProcessingParameterFeatureSource(
            "EXTENT", "Flood extent (optional)", [QgsProcessing.TypeVectorPolygon], optional=True))
        self.addParameter(QgsProcessingParameterRasterDestination("OUTPUT", "Flood depth"))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core import raster as R, vector as V, surface as SF
        dem_path = self.raster_path(parameters, "DEM", context)
        try:
            grid, wkts, wcrs = _dem_grid_and_extent(self, parameters, context, feedback, dem_path)
            dem, _ = R.read_on_grid(dem_path, grid, 1, "bilinear", zero_nodata_auto=True)
        except R.RasterError as exc:
            raise QgsProcessingException(str(exc))
        ext = V.rasterize_wkts(wkts, wcrs, grid) if wkts else None
        out = SF.depth_constant_stage(self.parameterAsDouble(parameters, "LEVEL", context), dem, ext)
        p = self.parameterAsOutputLayer(parameters, "OUTPUT", context)
        R.write_raster(p, out, grid, -9999.0, "Float32", {"CONTENT": "flood depth (m)"})
        n = int(np.isfinite(out).sum())
        feedback.pushInfo(f"{n:,} wet cells" + (f"; mean depth {float(np.nanmean(out)):.2f} m" if n else ""))
        return {"OUTPUT": p}
