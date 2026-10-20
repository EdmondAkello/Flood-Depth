# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Processing provider registration.

Everything is a QgsProcessingAlgorithm: parameter widgets, batch mode,
the Graphical Modeler, the history log and processing.run() come free.
"""

import os

from qgis.core import QgsProcessingProvider
from qgis.PyQt.QtGui import QIcon

from .alg_quick import QuickFloodDepthAlgorithm
from .alg_classify import ClassifyWaterAlgorithm
from .alg_supervised import SupervisedWaterAlgorithm
from .alg_surface import SurfaceFromExtentAlgorithm
from .alg_depth import DepthFromSurfaceAlgorithm, DepthConstantLevelAlgorithm
from .alg_compare import CompareExtentsAlgorithm


class FloodDepthProvider(QgsProcessingProvider):

    def loadAlgorithms(self):
        for alg in (QuickFloodDepthAlgorithm(), ClassifyWaterAlgorithm(), SupervisedWaterAlgorithm(),
                    SurfaceFromExtentAlgorithm(), DepthFromSurfaceAlgorithm(), DepthConstantLevelAlgorithm(),
                    CompareExtentsAlgorithm()):
            self.addAlgorithm(alg)

    def id(self):
        return "flooddepth"

    def name(self):
        return "FloodDepth"

    def longName(self):
        return "FloodDepth - flood extent, depth and uncertainty from imagery and a DEM"

    def icon(self):
        path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "icon.svg")
        return QIcon(path) if os.path.exists(path) else QgsProcessingProvider.icon(self)
