# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Plugin shell. Deliberately thin: it registers the Processing provider
and does nothing else. All functionality lives in the algorithms."""

from qgis.core import QgsApplication

from .processing_provider.provider import FloodDepthProvider


class FloodDepthPlugin(object):

    def __init__(self, iface):
        self.iface = iface
        self.provider = None

    def initGui(self):
        self.initProcessing()

    def initProcessing(self):
        self.provider = FloodDepthProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def unload(self):
        if self.provider is not None:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None
