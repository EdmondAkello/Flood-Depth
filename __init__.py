# -*- coding: utf-8 -*-
# FloodDepth - flood extent and depth from satellite imagery and a DEM.
# Licensed under the GNU General Public License v3 or later.


def classFactory(iface):
    from .plugin import FloodDepthPlugin
    return FloodDepthPlugin(iface)
