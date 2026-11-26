# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Plugin shell: the Processing provider, a Raster menu and one toolbar
button. Everything else lives in the algorithms; nothing heavy is imported
until a tool is opened, so QGIS starts as fast with the plugin as without."""

import os

from qgis.core import QgsApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QMenu

try:  # Qt 6 moved QAction to QtGui
    from qgis.PyQt.QtGui import QAction
except ImportError:  # Qt 5
    from qgis.PyQt.QtWidgets import QAction

from .processing_provider.provider import FloodDepthProvider

MENU = "&FloodDepth"
HERE = os.path.dirname(__file__)
GUIDE = os.path.join(HERE, "help", "user-guide.html")
GUIDE_URL = "https://github.com/EdmondAkello/Flood-Depth/blob/main/docs/USER_GUIDE.md"

# (algorithm id, menu text, tooltip, theme icon or None for the plugin icon)
TOOLS = (
    ("flooddepth:quickflooddepth", "Quick Flood Depth…",
     "Satellite image (or flood extent) + DEM -> flood depth, uncertainty and a report", None),
    ("flooddepth:flooddepthseries", "Flood Depth over Time…",
     "Several image dates -> depth per date, maximum depth and dates flooded", "/mIconTemporal.svg"),
    ("flooddepth:preparedem", "Prepare DEM…",
     "Merge DEM tiles, clip to the area and reproject to UTM (bilinear), with a DEM check", "/mActionAddRasterLayer.svg"),
    (None, None, None, None),
    ("flooddepth:classify_inundation", "Classify Water from Imagery…",
     "Step 1 on its own: a water mask and polygons from the image", "/mActionCalculateRaster.svg"),
    ("flooddepth:compareextents", "Compare Flood Extents…",
     "CSI, hit rate and false alarms of an extent against a reference", "/mActionMeasureArea.svg"),
)


class FloodDepthPlugin(object):

    def __init__(self, iface):
        self.iface = iface
        self.provider = None
        self.actions = []
        self.toolbar_action = None
        self.toolbar_menu = None

    def initProcessing(self):
        if self.provider is None:
            self.provider = FloodDepthProvider()
            QgsApplication.processingRegistry().addProvider(self.provider)

    def initGui(self):
        self.initProcessing()
        if self.iface is None:
            return
        window = self.iface.mainWindow()
        icon = QIcon(os.path.join(HERE, "icon.svg"))
        menu_actions = []
        for alg_id, text, tip, theme in TOOLS:
            if alg_id is None:
                menu_actions.append(None)
                continue
            if QgsApplication.processingRegistry().algorithmById(alg_id) is None:
                continue
            a = QAction(QgsApplication.getThemeIcon(theme) if theme else icon, text, window)
            a.setToolTip(tip)
            a.setStatusTip(tip)
            a.triggered.connect(lambda _checked=False, i=alg_id: self.open_tool(i))
            menu_actions.append(a)
        guide = QAction(QgsApplication.getThemeIcon("/mActionHelpContents.svg"), "User Guide", window)
        guide.setToolTip("Open the FloodDepth user guide (installed with the plugin; works offline)")
        guide.triggered.connect(self.open_guide)
        menu_actions += [None, guide]

        self.toolbar_menu = QMenu("FloodDepth", window)
        for a in menu_actions:
            if a is None:
                self.toolbar_menu.addSeparator()
                continue
            self.iface.addPluginToRasterMenu(MENU, a)
            self.actions.append(a)
            self.toolbar_menu.addAction(a)
        # one toolbar button: a click opens Quick flood depth, the arrow lists every tool
        self.toolbar_action = QAction(icon, "FloodDepth", window)
        self.toolbar_action.setToolTip("FloodDepth: Quick flood depth (click) - more tools on the arrow")
        self.toolbar_action.setMenu(self.toolbar_menu)
        self.toolbar_action.triggered.connect(lambda _checked=False: self.open_tool("flooddepth:quickflooddepth"))
        self.iface.addRasterToolBarIcon(self.toolbar_action)
        try:
            from qgis.PyQt.QtWidgets import QToolButton
            button = self.iface.rasterToolBar().widgetForAction(self.toolbar_action)
            if isinstance(button, QToolButton):
                button.setPopupMode(QToolButton.ToolButtonPopupMode.MenuButtonPopup)
        except (AttributeError, RuntimeError, TypeError):
            pass

    def open_tool(self, alg_id):
        """The algorithm's Processing dialog (batch, history and the modeller still work)."""
        try:
            from processing import execAlgorithmDialog
            execAlgorithmDialog(alg_id, {})
        except Exception as exc:  # noqa: BLE001 - shown to the user, never raised into QGIS
            from qgis.core import Qgis
            self.iface.messageBar().pushMessage("FloodDepth", f"Could not open {alg_id}: {exc}",
                                                level=Qgis.MessageLevel.Warning, duration=8)

    def open_guide(self, *_args):
        from qgis.PyQt.QtCore import QUrl
        from qgis.PyQt.QtGui import QDesktopServices
        url = QUrl.fromLocalFile(GUIDE) if os.path.exists(GUIDE) else QUrl(GUIDE_URL)
        QDesktopServices.openUrl(url)

    def unload(self):
        for a in self.actions:
            if self.iface is not None:
                self.iface.removePluginRasterMenu(MENU, a)
        self.actions = []
        if self.toolbar_action is not None and self.iface is not None:
            self.iface.removeRasterToolBarIcon(self.toolbar_action)
        self.toolbar_action = None
        self.toolbar_menu = None
        if self.provider is not None:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None
