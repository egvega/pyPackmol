"""
pyPackmol - PyMOL plugin for Packmol
Described at PyMOL wiki: http://www.pymolwiki.org/index.php/
Authors: Esteban Gabriel Vega-Hissi
email: egvega@gmail.com, egvega@unsl.edu.ar
Date: October 2026
License: GNU General Public License
Version 0.56 - Qt (PyMOL 2.0+) port of the Tkinter/Pmw plugin
"""

__version__ = '0.56'

dialog = None


def __init_plugin__(app=None):
    """Adds this plugin to the PyMOL Plugin menu."""
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('pyPackmol', run_plugin_gui)


def run_plugin_gui():
    """Opens the plugin window, reusing it if it was already created."""
    global dialog
    if dialog is None:
        from .gui.dialog import PackmolDialog
        dialog = PackmolDialog()  # applies the startup path
    elif not dialog.isVisible():
        dialog.apply_startup_path()
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()
