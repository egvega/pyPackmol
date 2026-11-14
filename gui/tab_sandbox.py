"""'Sandbox' tab: get coordinates from objects, create and move CGO helper objects."""
import re

import numpy as np
from pymol import cmd, stored
from pymol.Qt import QtCore, QtGui, QtWidgets

from ..buildCGO import CGO_DICT, crisscross, cgo_prisma_lines
from ..volume import volCGO
from ..core.rules import REGION_NPARAMS, PARAM_NAMES
from ..core.pymolutil import is_label, exclude_labels
from .widgets import NumEdit, RefreshComboBox, ColorPicker, button, color_rgb, group, hbox

Qt = QtCore.Qt
COLUMNS = ['suffix', 'a', 'b', 'c', 'd', 'e', 'f', 'color']
COLOR_COL = len(COLUMNS) - 1  # own color '#rrggbb' in Qt.UserRole, None = default CGO color
SANDBOX_KINDS = ['fixed', 'box', 'sphere']


def _elem(name):
    """CGO type from an object name like 'box_mymol'."""
    return name[:name.find('_')]


class SandboxTab(QtWidgets.QWidget):
    def __init__(self, dlg, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.dlg = dlg
        main = QtWidgets.QVBoxLayout(self)

        # Get Coordinates
        g, gl = group('Get Coordinates')
        main.addWidget(g)
        self.objPymol = RefreshComboBox(
            lambda: ['(all)'] + [n for n in cmd.get_names_of_type('object:molecule') if not is_label(n)]
            + cmd.get_names('selections'))
        gl.addLayout(hbox('Object/Selection', self.objPymol, None))
        self.optFixed = QtWidgets.QCheckBox('Baricenter')
        self.optBox = QtWidgets.QCheckBox('Box')
        self.optSphere = QtWidgets.QCheckBox('Sphere')
        gl.addLayout(hbox(self.optFixed, self.optBox, self.optSphere, None, button('Get', self.get_coords)))

        # Create Objects
        g, gl = group('Create Objects')
        main.addWidget(g, 1)
        self.objCGO = QtWidgets.QComboBox()
        self.objCGO.addItems(SANDBOX_KINDS)
        self.objCGO.currentTextChanged.connect(self._kind_changed)
        self.nameCGO = QtWidgets.QLineEdit()
        self.nameCGO.setMaximumWidth(120)
        grid = QtWidgets.QGridLayout()
        grid.addWidget(QtWidgets.QLabel('Object'), 0, 0)
        grid.addWidget(self.objCGO, 1, 0)
        grid.addWidget(QtWidgets.QLabel('suffix'), 0, 1)
        grid.addWidget(self.nameCGO, 1, 1)
        self.cgo_labels, self.cgo_params = [], []
        for i, name in enumerate(['x', 'y', 'z', 'a', 'b', 'c']):
            label = QtWidgets.QLabel(name)
            edit = NumEdit(0.0, float, 5)
            self.cgo_labels.append(label)
            self.cgo_params.append(edit)
            grid.addWidget(label, 0, i + 2)
            grid.addWidget(edit, 1, i + 2)
        self.cgoColor = ColorPicker(self.dlg.settings['cgoColor'])
        grid.addWidget(self.cgoColor, 1, 8)
        grid.addWidget(button('Add', self.add_obj), 1, 9)
        grid.setColumnStretch(10, 1)
        gl.addLayout(grid)
        self._kind_changed(self.objCGO.currentText())

        self.table = QtWidgets.QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.table.cellChanged.connect(self._cell_changed)
        self.table.cellDoubleClicked.connect(self._cell_double_clicked)
        side = QtWidgets.QVBoxLayout()
        side.addWidget(button('Remove', lambda: self.remove_obj('remove')))
        side.addWidget(button('Clear', lambda: self.remove_obj('clear')))
        side.addWidget(button('Show', self.vis_obj))
        side.addWidget(button('Use', self.use_coords))
        side.addStretch(1)
        gl.addLayout(hbox(self.table, side), 1)

        # Translate (Move) Objects
        g, gl = group('Translate (Move) Objects')
        main.addWidget(g)
        self.objMove = RefreshComboBox(
            lambda: [n for n in cmd.get_names_of_type('object:molecule') if not is_label(n)]
            + cmd.get_names_of_type('object:cgo'))
        row = [ 'Object', self.objMove]
        self.steps = []
        for axis, name in enumerate('XYZ'):
            step = NumEdit(1.0, float, 5)
            self.steps.append(step)
            up = button(u'▲', lambda a=axis: self.change_loc(a, 1))
            down = button(u'▼', lambda a=axis: self.change_loc(a, -1))
            for b in (up, down):
                b.setMaximumWidth(28)
            row += [name + ' step', step, up, down]
        gl.addLayout(hbox(*(row + [None])))

    # -------------------------------------------------------------------
    def _kind_changed(self, kind):
        n = REGION_NPARAMS[kind]
        for i, (label, edit) in enumerate(zip(self.cgo_labels, self.cgo_params)):
            edit.setEnabled(i < n)
            if i < n:
                label.setText(PARAM_NAMES[kind][i])

    def _add_row(self, name, values, color=None):
        self.table.blockSignals(True)
        row = self.table.rowCount()
        self.table.insertRow(row)
        item = QtWidgets.QTableWidgetItem(name)
        item.setData(Qt.UserRole, name)
        self.table.setItem(row, 0, item)
        for col, v in enumerate(values[:COLOR_COL - 1]):
            self.table.setItem(row, col + 1, QtWidgets.QTableWidgetItem('%.3f' % float(v)))
        item = QtWidgets.QTableWidgetItem()
        item.setFlags(item.flags() & ~Qt.ItemIsEditable)  # double-click opens a color dialog
        item.setToolTip('Double-click to change the color')
        self.table.setItem(row, COLOR_COL, item)
        self._set_row_color(row, color)
        self.table.blockSignals(False)

    def _set_row_color(self, row, color):
        item = self.table.item(row, COLOR_COL)
        item.setData(Qt.UserRole, color)
        item.setText(color or 'default')
        item.setBackground(QtGui.QColor(color) if color else QtGui.QBrush())

    def _row_color(self, row):
        return self.table.item(row, COLOR_COL).data(Qt.UserRole)

    def _row_name(self, row):
        return self.table.item(row, 0).data(Qt.UserRole)

    def _row_values(self, row):
        name = self._row_name(row)
        values = []
        for col in range(1, REGION_NPARAMS[_elem(name)] + 1):
            it = self.table.item(row, col)
            try:
                values.append(float(it.text()) if it else 0.)
            except ValueError:
                values.append(0.)
        return np.array(values)

    def _current_row(self):
        rows = self.table.selectionModel().selectedRows()
        return rows[0].row() if rows else -1

    def _load(self, elem, values, objName, color=None):
        color = color_rgb(color) if color else self.dlg.settings.cgo_color()
        if elem == 'fixed':
            crisscross(values[0], values[1], values[2], 0.5, objName, color)
            return
        if elem == 'box':
            obj = cgo_prisma_lines(*[float(v) for v in values], color)
        else:  # sphere
            obj = CGO_DICT['sphere']([float(v) for v in values], color)
        currentView = cmd.get_view()
        cmd.load_cgo(obj, objName)
        cmd.set_view(currentView)

    # -------------------------------------------------------------------
    def get_coords(self):
        sele = self.objPymol.currentText() or '(all)'
        base = re.sub(r'\W', '', sele) or 'all'
        sele = exclude_labels(sele)  # constraint label pseudoatoms are not part of the molecules
        boundingBox = np.array(cmd.get_extent(sele))  # [[x0,y0,z0],[x1,y1,z1]]
        color = self.cgoColor.color()
        if self.optFixed.isChecked():
            f = boundingBox[0] + (boundingBox[1] - boundingBox[0]) / 2.
            objName = cmd.get_unused_name('fixed_' + base)
            self._load('fixed', f, objName, color)
            self._add_row(objName, list(f) + [0., 0., 0.], color)
            self.dlg.vol_cgo[objName] = (0, 0)
        if self.optBox.isChecked() or self.optSphere.isChecked():
            stored.vdwList = [0.]
            cmd.iterate_state(-1, sele, 'stored.vdwList.append(vdw)')
            bb = np.array(boundingBox)
            bb[0] -= max(stored.vdwList)
            bb[1] += max(stored.vdwList)
        if self.optBox.isChecked():
            c = bb.flatten()
            objName = cmd.get_unused_name('box_' + base)
            self._load('box', c, objName, color)
            self._add_row(objName, c, color)
            self.dlg.vol_cgo[objName] = volCGO(['box'] + list(c))
        if self.optSphere.isChecked():
            ce = bb[0] + (bb[1] - bb[0]) / 2.
            r = np.sqrt(((bb[1] - ce) ** 2).sum())
            s = np.append(ce, r)
            objName = cmd.get_unused_name('sphere_' + base)
            self._load('sphere', s, objName, color)
            self._add_row(objName, s, color)
            self.dlg.vol_cgo[objName] = volCGO(['sphere'] + list(s))

    def add_obj(self):
        kind = self.objCGO.currentText()
        objName = cmd.get_unused_name(kind + '_' + self.nameCGO.text().strip())
        values = [self.cgo_params[i].value() for i in range(REGION_NPARAMS[kind])]
        color = self.cgoColor.color()
        self._add_row(objName, values, color)
        self.dlg.vol_cgo[objName] = volCGO([kind] + values)
        self._load(kind, values, objName, color)

    def use_coords(self):
        row = self._current_row()
        if row < 0:
            return
        self.dlg.setup.set_constraint(_elem(self._row_name(row)), self._row_values(row),
                                      self._row_color(row))

    def remove_obj(self, op):
        if op == 'clear':
            for row in range(self.table.rowCount()):
                cmd.delete(self._row_name(row))
            self.table.setRowCount(0)
        elif op == 'remove':
            row = self._current_row()
            if row >= 0:
                cmd.delete(self._row_name(row))
                self.table.removeRow(row)

    def vis_obj(self, row=None):
        if row is None:
            row = self._current_row()
        if row < 0:
            return
        objName = self._row_name(row)
        elem = _elem(objName)
        values = self._row_values(row)
        self.dlg.vol_cgo[objName] = volCGO([elem] + list(values))
        cmd.delete(objName)
        self._load(elem, values, objName, self._row_color(row))

    def _cell_double_clicked(self, row, col):
        if col != COLOR_COL:
            return
        current = self._row_color(row) or self.dlg.settings['cgoColor']
        color = QtWidgets.QColorDialog.getColor(QtGui.QColor(current), self, 'Choose CGO color')
        if color.isValid():  # not cancelled
            self.table.blockSignals(True)
            self._set_row_color(row, color.name())
            self.table.blockSignals(False)
            self.vis_obj(row)

    def _cell_changed(self, row, col):
        if col == 0:  # rename, the type prefix is kept
            item = self.table.item(row, 0)
            old = item.data(Qt.UserRole)
            prefix = old[:old.find('_') + 1]
            new = item.text().strip()
            if not new.startswith(prefix):
                new = prefix + new
            self.table.blockSignals(True)
            if new == prefix or (new != old and new in cmd.get_names('all')):
                item.setText(old)  # empty or already used name
                self.table.blockSignals(False)
                return
            item.setText(new)
            item.setData(Qt.UserRole, new)
            self.table.blockSignals(False)
            cmd.delete(old)
            self.dlg.vol_cgo.pop(old, None)
        self.table.selectRow(row)
        self.vis_obj(row)

    def change_loc(self, axis, d):
        name = self.objMove.currentText()
        if not name:
            return
        delta = self.steps[axis].value() * float(d)
        for row in range(self.table.rowCount()):
            if self._row_name(row).lower() == name.lower():
                self.table.selectRow(row)
                elem = _elem(name)
                values = self._row_values(row)
                values[axis] += delta
                if elem == 'box':
                    values[axis + 3] += delta
                self.table.blockSignals(True)
                for col, v in enumerate(values):
                    self.table.item(row, col + 1).setText('%.3f' % v)
                self.table.blockSignals(False)
                self.dlg.vol_cgo[name] = volCGO([elem] + list(values))
                break
        vector = [0., 0., 0.]
        vector[axis] = delta
        cmd.translate(vector, object=name)
