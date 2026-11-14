"""'Plugin Settings' tab. Also the settings store read by the other tabs."""
import pymol.plugins
from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

from ..core.pymolutil import show_axis
from ..volume import calcVolume3
from ..molovol import calcVolumeMoloVol
from .widgets import NumEdit, FileField, button, color_rgb, group, hbox, style_color_button

CHECKS = [('fullPath', 'Use full path for files'),
          ('loadStruct', 'Load structures into pymol'),
          ('loadVisual', 'Visualize constraints'),
          ('loadPartial', 'Load Partial Solved structures'),
          ('loadOutput', 'Load Output Final Solved structure'),
          ('showAxis', 'Show 3D axis object'),
          ('showBaricenter', 'Show Baricenter for Restart Files')]

DEFAULTS = {'fullPath': True, 'loadStruct': False, 'loadVisual': True, 'loadPartial': False,
            'loadOutput': True, 'showAxis': True, 'showBaricenter': False,
            'solventRadius': 0.0, 'lenTolerance': 0.10, 'packWarningN': 5,
            'Npoints': 100, 'cgoColor': '#ffffff', 'cgoLineWidth': 2.0,
            'startPathOn': False, 'startPath': '',
            'volMethod': 'ConvexHull', 'mvGrid': 0.2, 'mvProbe2': 0.0,
            'mvRotate': True, 'mvBlocks': True}
VOL_METHODS = ['ConvexHull', 'MoloVol']
METHOD_SOLVENT = {'ConvexHull': 0.0, 'MoloVol': 1.2}  # Solvent Radius set when the user picks a method


class SettingsTab(QtWidgets.QWidget):
    """settings['key'] returns the current value (bool, float or int)."""

    def __init__(self, dlg, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.dlg = dlg
        main = QtWidgets.QVBoxLayout(self)
        g, lay = group('Settings')
        main.addWidget(g)
        main.addStretch(1)

        # folder set as Actual Path whenever the plugin window is opened
        self.startPathOn = QtWidgets.QCheckBox('Startup path')
        self.startPath = FileField('...', 'Select startup folder', 'dir', button_after=True)
        self.startPathOn.toggled.connect(self.startPath.setEnabled)
        self.startPath.setEnabled(False)
        lay.addLayout(hbox(self.startPathOn, self.startPath))

        self.checks = {}
        for key, text in CHECKS:
            cb = QtWidgets.QCheckBox(text)
            self.checks[key] = cb
            lay.addWidget(cb)
        self.checks['showAxis'].toggled.connect(lambda on: show_axis(on))

        self.solventRadius = NumEdit(0.0, float, 10)
        lay.addLayout(hbox('Solvent Radius :', self.solventRadius, None))

        self.lenTolLabel = QtWidgets.QLabel()
        self.lenTolSlider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.lenTolSlider.setRange(5, 100)  # 0.05 - 1.00 A
        self.lenTolSlider.setMaximumWidth(200)
        self.lenTolSlider.valueChanged.connect(
            lambda v: self.lenTolLabel.setText('%.2f' % (v / 100.)))
        lay.addLayout(hbox('Length Tolerance (A): ', self.lenTolLabel, self.lenTolSlider, None))

        self.packWarningN = NumEdit(5, int, 10)
        lay.addLayout(hbox('Packing Warnings Limit :', self.packWarningN, None))

        self.volMethod = QtWidgets.QComboBox()
        self.volMethod.addItems(VOL_METHODS)
        self.volMethod.currentTextChanged.connect(self._method_changed)
        # activated is emitted only for user picks, so loading saved settings keeps their Solvent Radius
        self.volMethod.activated.connect(
            lambda i: self.solventRadius.setValue(METHOD_SOLVENT[self.volMethod.currentText()]))
        lay.addLayout(hbox('Volume Calculation method :', self.volMethod, None))

        self.Npoints = NumEdit(100, int, 4)
        lay.addLayout(hbox('N. points for Volume Calculation (ConvexHull) :', self.Npoints, None))
        self.mvGrid = NumEdit(0.2, float, 6)
        lay.addLayout(hbox('MoloVol grid resolution (A) :', self.mvGrid, None))
        self.mvProbe2 = NumEdit(0.0, float, 6)
        lay.addLayout(hbox('MoloVol large probe (A, 0 = off) :', self.mvProbe2, None))
        self.mvRotate = QtWidgets.QCheckBox('MoloVol: rotate to principal axes')
        self.mvBlocks = QtWidgets.QCheckBox('MoloVol: skip empty grid blocks')
        lay.addLayout(hbox(self.mvRotate, self.mvBlocks, None))

        self._cgoColor = DEFAULTS['cgoColor']
        self.colorButton = button(self._cgoColor, self.pick_color)
        lay.addLayout(hbox('CGO objects color :', self.colorButton, None))
        self.cgoLineWidth = NumEdit(2.0, float, 6)
        self.cgoLineWidth.editingFinished.connect(self._apply_line_width)
        lay.addLayout(hbox('CGO line width :', self.cgoLineWidth, None))

        lay.addLayout(hbox(None, button('Defaults', self.defaults), button('Save', self.save)))

        self.defaults()
        self.load()

    def _method_changed(self, method):
        molovol = method == 'MoloVol'
        self.Npoints.setEnabled(not molovol)
        self.mvGrid.setEnabled(molovol)
        self.mvProbe2.setEnabled(molovol)
        self.mvRotate.setEnabled(molovol)
        self.mvBlocks.setEnabled(molovol)
        self.dlg.vol_mol.clear()  # cached molecular volumes were computed with the other method

    def mol_volume(self, name):
        """Molecular volume of a PyMOL object with the selected method and its settings.

        MoloVol raises molovol.InputError on wrong settings and MemoryError if the grid is too fine.
        """
        if self['volMethod'] == 'MoloVol':
            return calcVolumeMoloVol(name, self['solventRadius'], self['mvProbe2'] or None, self['mvGrid'],
                                     rotate=self['mvRotate'], blocks=self['mvBlocks'])
        return calcVolume3(name, self['solventRadius'], self['Npoints'])

    def _apply_line_width(self):
        """PyMOL's global cgo_line_width (all CGO lines: constraints and Sandbox objects)."""
        width = self.cgoLineWidth.value()
        if width <= 0:
            width = DEFAULTS['cgoLineWidth']
            self.cgoLineWidth.setValue(width)
        cmd.set('cgo_line_width', width)

    def cgo_color(self):
        """CGO color as (r, g, b) floats in 0-1."""
        return color_rgb(self._cgoColor)

    def _set_color(self, name):
        name = QtGui.QColor(name).name()  # normalized '#rrggbb'
        changed = name != self._cgoColor
        self._cgoColor = name
        style_color_button(self.colorButton, name)
        # redraw constraints with the new color (setup tab does not exist during __init__)
        if changed and getattr(self.dlg, 'setup', None) is not None:
            self.dlg.setup.tree_op('clearVis')
            self.dlg.setup.tree_op('visAll')

    def pick_color(self):
        color = QtWidgets.QColorDialog.getColor(QtGui.QColor(self._cgoColor), self,
                                                'Choose color for CGO objects')
        if color.isValid():  # not cancelled
            self._set_color(color.name())

    def __getitem__(self, key):
        if key in self.checks:
            return self.checks[key].isChecked()
        if key == 'solventRadius':
            return self.solventRadius.value()
        if key == 'lenTolerance':
            return self.lenTolSlider.value() / 100.
        if key == 'packWarningN':
            return self.packWarningN.value()
        if key == 'Npoints':
            return max(1, self.Npoints.value())
        if key == 'cgoColor':
            return self._cgoColor
        if key == 'cgoLineWidth':
            return self.cgoLineWidth.value()
        if key == 'startPathOn':
            return self.startPathOn.isChecked()
        if key == 'startPath':
            return self.startPath.text()
        if key == 'volMethod':
            return self.volMethod.currentText()
        if key == 'mvGrid':
            return self.mvGrid.value()
        if key == 'mvProbe2':
            return self.mvProbe2.value()
        if key == 'mvRotate':
            return self.mvRotate.isChecked()
        if key == 'mvBlocks':
            return self.mvBlocks.isChecked()
        raise KeyError(key)

    def __setitem__(self, key, value):
        if key in self.checks:
            self.checks[key].setChecked(bool(value))
        elif key == 'solventRadius':
            self.solventRadius.setValue(float(value))
        elif key == 'lenTolerance':
            self.lenTolSlider.setValue(int(round(float(value) * 100)))
        elif key == 'packWarningN':
            self.packWarningN.setValue(int(value))
        elif key == 'Npoints':
            self.Npoints.setValue(int(value))
        elif key == 'startPathOn':
            self.startPathOn.setChecked(bool(value))
        elif key == 'startPath':
            self.startPath.setText(str(value))
        elif key == 'cgoLineWidth':
            if float(value) <= 0:
                raise ValueError(value)
            self.cgoLineWidth.setValue(float(value))
            self._apply_line_width()
        elif key == 'cgoColor':
            if not QtGui.QColor(str(value)).isValid():
                raise ValueError(value)
            self._set_color(str(value))
        elif key == 'volMethod':
            self.volMethod.setCurrentText(value if value in VOL_METHODS else DEFAULTS['volMethod'])
            self._method_changed(self.volMethod.currentText())
        elif key == 'mvGrid':
            if float(value) <= 0:
                raise ValueError(value)
            self.mvGrid.setValue(float(value))
        elif key == 'mvProbe2':
            self.mvProbe2.setValue(max(0.0, float(value)))
        elif key == 'mvRotate':
            self.mvRotate.setChecked(bool(value))
        elif key == 'mvBlocks':
            self.mvBlocks.setChecked(bool(value))
        else:
            raise KeyError(key)

    def defaults(self):
        for key, value in DEFAULTS.items():
            self[key] = value

    def load(self):
        """Override defaults with the values stored in PyMOL preferences."""
        for key in DEFAULTS:
            value = pymol.plugins.pref_get(key, None)
            if value is not None:
                try:
                    self[key] = value
                except (TypeError, ValueError):
                    pass

    def save(self):
        for key in DEFAULTS:
            pymol.plugins.pref_set(key, self[key])
        if not pymol.plugins.pref_get('instantsave', True):
            pymol.plugins.pref_save()
