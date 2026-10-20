# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Step 2: water surface and depth from a flood extent (polygons or mask)."""

import numpy as np
from qgis.core import (
    QgsProcessing,
    QgsProcessingException,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterVectorDestination,
)

from .base import FloodDepthAlgorithm

INTERP = ["idw", "nearest"]


class SurfaceFromExtentAlgorithm(FloodDepthAlgorithm):

    def name(self):
        # id kept from 0.3 (was an IDW interpolation over the whole DTM)
        return "estimate_stage_dem_boundary"

    def displayName(self):
        return "2. Water surface and depth from a flood extent"

    def group(self):
        return "Step by step"

    def groupId(self):
        return "steps"

    def shortHelpString(self):
        return (
            "Estimates the water-surface elevation over a flood extent from the terrain along its edge, and "
            "the depth (FwDET boundary method, Cohen et al. 2018, 2019).\n\n"
            "Waterline cells are flood cells next to valid dry ground; edges against NoData, the grid edge or "
            "permanent water are not used. Edge elevations far from their neighbours (embankments, canopy in a "
            "DSM) are dropped unless the ground there is steep. Each flood cell takes the elevation of its "
            "nearest waterline cell, then a focal mean smooths it where that changes it by less than the "
            "tolerance. depth = surface - DEM; ground above the surface gets depth 0.\n\n"
            "Extent: polygons, or a mask raster (value 1 = flood, 2 = permanent water, as written by step 1). "
            "With a HAND raster, the same is done on HAND (depth = edge HAND - HAND). "
            "A geographic DEM is processed on the WGS 84 / UTM grid of the area.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(
            "EXTENT", "Flood extent polygons", [QgsProcessing.TypeVectorPolygon], optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer(
            "MASK", "...or a flood mask raster (1 = flood, 2 = permanent water)", optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", "DEM"))
        self.addParameter(QgsProcessingParameterRasterLayer("HAND", "HAND raster (optional: HAND method)",
                                                            optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "PERM_POLYS", "Permanent water polygons", [QgsProcessing.TypeVectorPolygon], optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource(
            "AOI", "Area of interest", [QgsProcessing.TypeVectorPolygon], optional=True))
        self.number("CELL", "Analysis cell size (m; 0 = DEM resolution)", "double", 0.0, 0.0, adv=True)
        self.addParameter(self.advanced(QgsProcessingParameterEnum(
            "INTERP", "Surface interpolation", ["Inverse distance from the flood edge (GDAL)", "Nearest flood edge (FwDET v1)"],
            defaultValue=0)))
        self.number("SMOOTH", "Smoothing radius (cells)", "int", 2, 0, 50, adv=True)
        self.number("SMOOTH_TOL", "Smoothing tolerance (m)", "double", 0.5, 0.0, adv=True)
        self.number("OUT_RADIUS", "Edge outlier window radius (cells)", "int", 10, 0, 200, adv=True)
        self.number("OUT_K", "Edge outlier limit (robust SD)", "double", 3.0, 0.0, adv=True)
        self.number("OUT_MIN", "Edge outlier minimum (m)", "double", 0.5, 0.0, adv=True)
        self.number("STEEP", "Steep ground (slope %)", "double", 10.0, 0.0, adv=True)
        self.addParameter(QgsProcessingParameterRasterDestination("WSE", "Water surface"))
        self.addParameter(QgsProcessingParameterRasterDestination("DEPTH", "Flood depth"))
        self.addParameter(QgsProcessingParameterVectorDestination(
            "WATERLINE", "Waterline cells (used / dropped)", QgsProcessing.TypeVectorPoint, optional=True,
            createByDefault=False))

    def processAlgorithm(self, parameters, context, feedback):
        from ..core import raster as R, vector as V, surface as SF, stats as ST
        dem_path = self.raster_path(parameters, "DEM", context)
        hand_path = self.raster_path(parameters, "HAND", context, required=False)
        mask_path = self.raster_path(parameters, "MASK", context, required=False)
        wkts, wcrs = self.vector_wkts(parameters, "EXTENT", context, feedback)
        if not wkts and not mask_path:
            raise QgsProcessingException("Give flood-extent polygons or a flood mask raster.")
        aoi, aoi_crs = self.vector_wkts(parameters, "AOI", context, feedback)
        try:
            if mask_path:
                mg = R.grid_of(R.open_dataset(mask_path))
                ext, ext_wkt = mg.bounds, mg.wkt
            else:
                from ..core.pipeline import _wkts_bounds
                ext, ext_wkt = _wkts_bounds(wkts, wcrs, pad_frac=0.1)
            if aoi:
                from ..core.pipeline import _wkts_bounds
                ab, aw = _wkts_bounds(aoi, aoi_crs)
                ext = R.intersect(ext, R.transform_bounds(ab, aw, ext_wkt))
                if ext is None:
                    raise QgsProcessingException("The area of interest does not overlap the extent.")
            grid, notes = R.analysis_grid(dem_path, ext, ext_wkt,
                                          self.parameterAsDouble(parameters, "CELL", context))
            for n in notes:
                feedback.pushInfo(n)
            dem, dn = R.read_on_grid(dem_path, grid, 1, "bilinear", zero_nodata_auto=True)
            for n in dn:
                feedback.pushWarning("DEM: " + n)
            valid = np.isfinite(dem)
            other = None
            if mask_path:
                mk, _ = R.read_on_grid(mask_path, grid, 1, "near")
                flood = mk == 1
                other = mk == 2
                valid &= np.isfinite(mk)
            else:
                flood = V.rasterize_wkts(wkts, wcrs, grid)
            perm, pcrs = self.vector_wkts(parameters, "PERM_POLYS", context, feedback)
            if perm:
                pw = V.rasterize_wkts(perm, pcrs, grid)
                other = pw if other is None else (other | pw)
            if other is not None:
                flood &= ~other
            if aoi:
                valid &= V.rasterize_wkts(aoi, aoi_crs, grid)
            flood &= valid
            elev = dem
            if hand_path:
                elev, hn = R.read_on_grid(hand_path, grid, 1, "bilinear")
                valid &= np.isfinite(elev)
                feedback.pushInfo("HAND method: depth = edge HAND - HAND.")
            surf, info, used, dropped = SF.estimate_surface(
                elev, flood, valid, grid, other_water=other,
                outlier_radius=self.parameterAsInt(parameters, "OUT_RADIUS", context),
                outlier_k=self.parameterAsDouble(parameters, "OUT_K", context),
                outlier_min=self.parameterAsDouble(parameters, "OUT_MIN", context),
                smooth_radius=self.parameterAsInt(parameters, "SMOOTH", context),
                smooth_tolerance=self.parameterAsDouble(parameters, "SMOOTH_TOL", context),
                interpolation=INTERP[self.parameterAsEnum(parameters, "INTERP", context)],
                steep_slope_pct=self.parameterAsDouble(parameters, "STEEP", context),
                progress=self.progress(feedback))
            depth, neg = SF.depth_from_surface(surf, elev, flood & valid)
        except (ValueError, R.RasterError) as exc:
            raise QgsProcessingException(str(exc))
        feedback.pushInfo(f"Waterline {info['waterline_cells']:,} cells: {info['waterline_used']:,} used, "
                          f"{info['waterline_dropped']:,} dropped; {info['edge_cells_excluded']:,} edge cells at "
                          "NoData / grid edge not used.")
        st = ST.depth_stats(depth, grid.cell_area())
        feedback.pushInfo(f"Flooded {st.get('area_km2', 0):.3f} km2; depth mean {st.get('mean_m', float('nan')):.2f} m, "
                          f"P90 {st.get('p90_m', float('nan')):.2f} m, max {st.get('max_m', float('nan')):.2f} m; "
                          f"{neg:,} cells above the surface set to 0.")
        wp = self.parameterAsOutputLayer(parameters, "WSE", context)
        dp = self.parameterAsOutputLayer(parameters, "DEPTH", context)
        R.write_raster(wp, surf, grid, -9999.0, "Float32", {"CONTENT": "water surface" if not hand_path
                                                             else "water height above drainage"})
        R.write_raster(dp, depth, grid, -9999.0, "Float32", {"CONTENT": "flood depth (m)"})
        out = {"WSE": wp, "DEPTH": dp}
        if parameters.get("WATERLINE"):
            lp = self.parameterAsOutputLayer(parameters, "WATERLINE", context)
            rr, cc = np.nonzero(used | dropped)
            V.write_points(lp, "waterline", grid, rr, cc, {"used": used[rr, cc].astype("int32"),
                                                          "z": elev[rr, cc].astype("float64")})
            out["WATERLINE"] = lp
        return out
