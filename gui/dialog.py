"""Main plugin window: status label + one tab per tool + working directory row."""
import os

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from .. import __version__
from .tab_settings import SettingsTab
from .tab_setup import SetupTab
from .tab_pack import PackTab
from .tab_sandbox import SandboxTab
from .tab_volume import VolumeTab
from .tab_restart import RestartTab
from .widgets import ClickableLabel, ask_directory, ask_save, button, error, hbox, save_widget_image

STATUS = {'stdby': 'Packmol : stand-by', 'check': 'Packmol : Calculating Volume',
          'run': 'Packmol :  p a c k i n g ', 'restart': 'Packmol : loading restart file'}

ABOUT = """pyPackmol was developed by Esteban Gabriel Vega-Hissi.
For instructions on how to use this plugin, please
read the user guide (docs/USER-GUIDE.md in the
pyPackmol distribution).
v %s

pyPackmol is free software, distributed under the
GNU General Public License v3 (see LICENSE).

MoloVol volume method: pure-Python port of MoloVol
(https://github.com/jmaglic/MoloVol), Copyright (c) 2020-2025
Jasmin B. Maglic, Roy Lavendomme. MIT License (see LICENSE-MoloVol.txt).
"""


class PackmolDialog(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent, QtCore.Qt.Window)
        self.setWindowTitle('pyPackmol')
        self.resize(820, 680)
        self.vol_cgo = {}  # CGO name -> (volume, max length)
        self.vol_mol = {}  # molecule name -> (volume, length)

        main = QtWidgets.QVBoxLayout(self)
        main.setContentsMargins(4, 4, 4, 4)
        self.statusLabel = ClickableLabel()
        self.statusLabel.on_double_click = self.save_screenshot
        self.statusLabel.setToolTip('Double-click to save an image of this window')
        self.statusLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.statusLabel.setMinimumHeight(40)
        self.statusLabel.setStyleSheet('background: black; color: darkorange; '
                                       'font-family: Courier, monospace; font-size: 10pt; font-weight: bold;')
        main.addWidget(self.statusLabel)
        self.set_status('stdby')

        # settings first, the other tabs read it
        self.settings = SettingsTab(self)
        self.setup = SetupTab(self)
        self.pack = PackTab(self)
        self.sandbox = SandboxTab(self)
        self.volume = VolumeTab(self)
        self.restart = RestartTab(self)

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(self.setup, 'Set Up Input')
        self.tabs.addTab(self.pack, 'Pack')
        self.tabs.addTab(self.sandbox, 'Sandbox')
        self.tabs.addTab(self.volume, 'Volume Calculator')
        self.tabs.addTab(self.restart, 'Restart Files')
        self.tabs.addTab(self.settings, 'Plugin Settings')
        about = QtWidgets.QPlainTextEdit(ABOUT % __version__)
        about.setReadOnly(True)
        self.tabs.addTab(about, 'About')
        self.tabs.currentChanged.connect(self._tab_changed)
        main.addWidget(self.tabs, 1)

        # PyMOL working directory, used by every tab (relative paths, Packmol runs)
        self.cwdEdit = QtWidgets.QLineEdit()
        self.cwdEdit.setReadOnly(True)
        main.addLayout(hbox(button('Actual Path:', self.change_cwd), self.cwdEdit))
        self.apply_startup_path()

        self.setup.input_file.textChanged.connect(self.pack.inputDisplay.setText)

    def update_cwd(self):
        self.cwdEdit.setText(os.path.abspath(os.getcwd()))

    def apply_startup_path(self):
        """Change to the startup folder of Plugin Settings, if enabled (window opened)."""
        path = self.settings['startPath']
        if self.settings['startPathOn'] and path:
            if os.path.isdir(path):
                cmd.cd(path)
            else:
                print('pyPackmol: startup path not found: ' + path)
        self.update_cwd()

    def change_cwd(self):
        fn = ask_directory(self, 'Select New Working Directory')
        if fn:
            cmd.cd(fn)
        self.update_cwd()

    def set_status(self, key):
        self.statusLabel.setText(STATUS.get(key, key))
        QtWidgets.QApplication.processEvents()

    def save_screenshot(self):
        """Save the plugin window as a 600 dpi image (figures of the GUI)."""
        fn = ask_save(self, 'Save window image', 'image')
        if fn and not save_widget_image(self, fn):
            error(self, 'Error', 'Could not save ' + fn)

    def _tab_changed(self, index):
        if self.tabs.widget(index) is self.volume:
            self.volume.refresh()

    def changeEvent(self, event):
        if event.type() == QtCore.QEvent.ActivationChange and self.isActiveWindow():
            self.update_cwd()  # PyMOL 'cd' may have changed it
        QtWidgets.QWidget.changeEvent(self, event)
