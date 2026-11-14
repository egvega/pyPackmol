"""'Set Up Input' tab: build the Packmol rules tree, visualize and check constraints."""
import os
import re

import numpy as np
from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

from ..buildCGO import CGO_DICT
from ..volume import volCGO
from ..molovol import InputError
from ..core import rules
from ..core.rules import head
from ..core.checks import check_volumes
from ..core.pymolutil import unique_name, add_cgo_label, is_label
from .widgets import (NumEdit, FileField, LogView, ColorPicker, add_button, button, group, hbox,
                      ask_open, ask_save, info, error, ask_yes_no, color_rgb, random_color)

Qt = QtCore.Qt
SURFACES = ('plane', 'xygauss')  # constraints without volume
COLOR_ROLE = Qt.UserRole + 1  # own CGO color of a constraint rule ('#rrggbb'), unset = default


def _head(item):
    return head(item.text(0)) if item is not None else ''


def _is_constraint(text):
    """Constraint (or pbc) rule that is drawn as a CGO and can have its own color."""
    sub = text.split()
    return len(sub) > 1 and (sub[1].lower() in CGO_DICT or sub[0].lower() in rules.GLOBAL_DIRECTIONS)


def _is_volume(sub):
    """Split rule is a volume constraint with valid parameters (it can give a plane its size)."""
    if len(sub) < 2 or sub[1].lower() not in CGO_DICT or sub[1].lower() in SURFACES:
        return False
    return rules.valid_region(sub[1], sub[2:])


def _pbc_box(subunits):
    """Box parameters of a 'pbc' line: 3 values mean 0..max, 6 values min..max (None if wrong)."""
    try:
        values = [float(par) for par in subunits[1:]]
    except ValueError:
        return None
    if len(values) == 3:
        return [0., 0., 0.] + values
    if len(values) == 6:
        return values
    return None


class SetupTab(QtWidgets.QWidget):
    def __init__(self, dlg, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.dlg = dlg
        self.cgo_loaded = {}  # item key -> CGO object names loaded into PyMOL
        self._next_key = 0  # QTreeWidgetItem is unhashable, so each item gets an int key
        self.constraint_list = []
        self.new_const = False
        self._edit_vis = False  # re-visualizing after a manual edit: don't ask to reload structures
        self._old_names = []  # CGO names drawn before for the item being redrawn (reused)

        main = QtWidgets.QVBoxLayout(self)
        # Packmol input file LOAD/SAVE
        self.input_file = QtWidgets.QLineEdit()
        main.addLayout(hbox('Packmol input:', button('...', self.load_input), self.input_file,
                            button('Load', self.load_typed), button('Save', self.save_input)))

        cols = QtWidgets.QHBoxLayout()
        main.addLayout(cols, 1)
        sub = QtWidgets.QTabWidget()
        sub.addTab(self._build_basic(), 'Basic')
        sub.addTab(self._build_additional(), 'Additional')
        cols.addWidget(sub, 1)
        cols.addWidget(self._build_rules(), 1)

    # =========================================================== widgets ==
    def _build_basic(self):
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)

        # General Input
        g, gl = group('General Input', QtWidgets.QHBoxLayout)
        lay.addWidget(g)
        col = QtWidgets.QVBoxLayout()
        self.tolerance = NumEdit(2.0, float, 4)
        self.seed = NumEdit(-1, int, 7)
        col.addLayout(hbox('Tolerance:', self.tolerance, 'Seed:', self.seed, None))
        self.filetype = QtWidgets.QComboBox()
        self.filetype.addItems(['pdb', 'xyz', 'tinker'])
        col.addLayout(hbox('File Type: ', self.filetype, None))
        self.output = FileField('Out:', 'Select Output file', 'save', 'output')
        col.addWidget(self.output)
        gl.addLayout(col, 1)
        gl.addWidget(add_button('>', lambda: self.add_rule('basic')))

        # Structures
        g, gl = group('Structures')
        lay.addWidget(g)
        self.structure = FileField('...', 'Open Structure file', 'open', 'structure')
        gl.addLayout(hbox(self.structure, add_button('>', lambda: self.add_rule('structure'))))
        self.number = NumEdit(0, int, 10)
        gl.addLayout(hbox('Number:', self.number, None, add_button('>', lambda: self.add_rule('number'))))
        self.atoms = QtWidgets.QLineEdit()
        gl.addLayout(hbox('Atom numbers:', self.atoms, add_button('>', lambda: self.add_rule('atoms'))))
        self.atomRadius = NumEdit(0.0, float, 4)
        gl.addLayout(hbox('Atom radius:', self.atomRadius, None, add_button('>', lambda: self.add_rule('radius'))))
        self.resNumber = QtWidgets.QComboBox()
        self.resNumber.addItems(['0', '1', '2', '3'])
        self.changeChains = QtWidgets.QCheckBox('Change chains')
        gl.addLayout(hbox('ResId', self.resNumber, self.changeChains, None,
                          add_button('>', lambda: self.add_rule('relabel'))))

        # Constraints
        g, gl = group('Constraints')
        lay.addWidget(g)
        self.direction = QtWidgets.QComboBox()
        self.direction.addItems(rules.DIRECTIONS)
        self.region = QtWidgets.QComboBox()
        self.region.addItems(rules.REGIONS)
        self.region.currentTextChanged.connect(self._region_changed)
        self.direction.currentTextChanged.connect(self._direction_changed)
        grid = QtWidgets.QGridLayout()
        self.param_labels, self.params = [], []
        for i, name in enumerate('abcdefgh'):
            label = QtWidgets.QLabel(name)
            edit = NumEdit(0.0, float, 5)
            self.param_labels.append(label)
            self.params.append(edit)
            grid.addWidget(label, i // 3, 2 * (i % 3))
            grid.addWidget(edit, i // 3, 2 * (i % 3) + 1)
        col = QtWidgets.QVBoxLayout()
        col.addLayout(hbox('Direction', self.direction, 'Region', self.region, None))
        col.addLayout(grid)
        self.cgoColor = ColorPicker(self.dlg.settings['cgoColor'])
        col.addLayout(hbox(self.cgoColor, None))
        gl.addLayout(hbox(col, add_button('>', lambda: self.add_rule('constraint'))))
        self.centerKind = QtWidgets.QComboBox()
        self.centerKind.addItems(['center', 'centerofmass'])
        gl.addLayout(hbox(self.centerKind, None, add_button('>', lambda: self.add_rule('center'))))
        self.rotAxis = QtWidgets.QComboBox()
        self.rotAxis.addItems(['x', 'y', 'z'])
        self.rotAngle = NumEdit(0.0, float, 3)
        self.rotRange = NumEdit(0.0, float, 3)
        col = QtWidgets.QVBoxLayout()
        col.addWidget(QtWidgets.QLabel('Constrain Rotation'))
        col.addLayout(hbox('Axis', self.rotAxis, 'Angle:', self.rotAngle, 'Range:', self.rotRange, None))
        gl.addLayout(hbox(col, add_button('>', lambda: self.add_rule('constrain_rotation'))))
        lay.addStretch(1)

        self._region_changed(self.region.currentText())
        return w

    def _build_additional(self):
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)

        g, gl = group('Additional Input')
        lay.addWidget(g)
        self.additional_values = {}
        for a in rules.ADDITIONAL_SETTINGS:
            if a in rules.INT_SETTINGS:
                self.additional_values[a] = '0'
            elif a in rules.FLAG_SETTINGS:
                self.additional_values[a] = ''
            else:
                self.additional_values[a] = '0.0'
        self.additional_values['sidemax'] = '2000.0'  # Packmol default
        self.additional_values.update({'movefrac': '0.05', 'nloop0': '1000', 'maxmove': '10',
                                       'writecrd': 'output.crd', 'chain': 'A', 'segid': 'A',
                                       'connect': 'no'})
        self.additionalMenu = QtWidgets.QComboBox()
        self.additionalMenu.addItems(rules.ADDITIONAL_SETTINGS)
        self.additionalEntry = QtWidgets.QLineEdit()
        self.additionalEntry.textEdited.connect(
            lambda t: self.additional_values.__setitem__(self.additionalMenu.currentText(), t))
        self.additionalMenu.currentTextChanged.connect(self._additional_changed)
        col = QtWidgets.QVBoxLayout()
        col.addWidget(self.additionalMenu)
        col.addWidget(self.additionalEntry)
        gl.addLayout(hbox(col, add_button('>', lambda: self.add_rule('additional'))))
        self._additional_changed(self.additionalMenu.currentText())

        self.useShortTol = QtWidgets.QCheckBox('use_short_tol')
        self.shortTolDist = NumEdit(0.0, float, 7)
        self.shortTolScale = NumEdit(0.0, float, 7)
        self.useShortTol.toggled.connect(self._short_tol_toggled)
        self._short_tol_toggled(False)
        col = QtWidgets.QVBoxLayout()
        col.addWidget(self.useShortTol)
        col.addLayout(hbox('short_tol_dist:', self.shortTolDist, None))
        col.addLayout(hbox('short_tol_scale:', self.shortTolScale, None))
        gl.addLayout(hbox(col, add_button('>', lambda: self.add_rule('additionalUST'))))

        g, gl = group('Restarts')
        lay.addWidget(g)
        self.restartTo = FileField('To...', 'Select Restart file', 'save', 'restart')
        gl.addLayout(hbox(self.restartTo, add_button('>', lambda: self.add_rule('restartTo'))))
        self.restartFrom = FileField('From...', 'Open Restart file', 'open', 'restart')
        gl.addLayout(hbox(self.restartFrom, add_button('>', lambda: self.add_rule('restartFrom'))))
        gl.addLayout(hbox('check', None, add_button('>', lambda: self.add_rule('check'))))
        lay.addStretch(1)
        return w

    def _build_rules(self):
        g, gl = group('Packmol Rules')
        self.tree = QtWidgets.QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        # double-click is handled by _item_double_clicked (color swatch or text edit)
        self.tree.setEditTriggers(QtWidgets.QAbstractItemView.EditKeyPressed)
        self.tree.itemDoubleClicked.connect(self._item_double_clicked)
        self.tree.itemChanged.connect(self._item_edited)
        gl.addWidget(self.tree, 2)
        gl.addLayout(hbox(button('Add Blank', lambda: self.tree_op('blank')),
                          button(u'↑', lambda: self.tree_op('moveUp')),
                          button(u'↓', lambda: self.tree_op('moveDown')),
                          button('Remove', lambda: self.tree_op('remove')),
                          button('Clear Rules', lambda: self.tree_op('clear'))))
        gl.addLayout(hbox(button('Visualize', lambda: self.tree_op('visSel')),
                          button('Visualize All', lambda: self.tree_op('visAll')),
                          button('Clear Visualization', lambda: self.tree_op('clearVis'))))
        gl.addLayout(hbox(button('Check Volume Constraints', lambda: self.tree_op('check')),
                          button('Clear Check', lambda: self.tree_op('clearCheck')), None))
        self.checkText = LogView()
        gl.addWidget(self.checkText, 1)
        return g

    # ===================================================== field helpers ==
    def _region_changed(self, choice):
        if choice == 'fixed':
            self.direction.setCurrentIndex(-1)
        elif self.direction.currentIndex() < 0:
            self.direction.setCurrentIndex(0)
        n = rules.REGION_NPARAMS[choice]
        names = rules.PARAM_NAMES[choice]
        for i, (label, edit) in enumerate(zip(self.param_labels, self.params)):
            edit.setEnabled(i < n)
            label.setText(names[i] if i < n else 'abcdefgh'[i])

    def _direction_changed(self, choice):
        # pbc takes the box parameters, so the region is locked to box
        if choice in rules.GLOBAL_DIRECTIONS:
            self.region.setCurrentText('box')
            self.region.setEnabled(False)
        else:
            self.region.setEnabled(True)

    def set_constraint(self, region, values, color=None):
        """Fill the constraint fields (used by the Sandbox 'Use' button)."""
        if self.direction.currentText() in rules.GLOBAL_DIRECTIONS and region != 'box':
            self.direction.setCurrentText('inside')
        self.region.setCurrentText(region)
        self._region_changed(region)
        for edit, v in zip(self.params, values):
            edit.setValue(float(v))
        self.cgoColor.setColor(color)

    def _rule_path(self, fn, what):
        """
        File name as written in a structure, output or restart rule: absolute with 'Use full path
        for files', otherwise relative to the Actual Path (absolute if on another
        drive). Warns about spaces, which end the file name for Packmol.
        """
        if not fn:
            return fn
        if self.dlg.settings['fullPath']:
            fn = os.path.abspath(fn)
        else:
            try:
                fn = os.path.relpath(os.path.abspath(fn))
            except ValueError:  # Windows: another drive than the working directory
                fn = os.path.abspath(fn)
        if ' ' in fn:
            info(self, 'Warning', 'The %s path contains spaces: %s\nPackmol reads only up to the first '
                 'space. Rule has been added anyway.' % (what, fn))
        return fn

    def _additional_changed(self, key):
        self.additionalEntry.setText(self.additional_values[key])
        self.additionalEntry.setEnabled(key not in rules.FLAG_SETTINGS)

    def _additional_value(self, key):
        text = self.additional_values[key].strip()
        if key in rules.FLAG_SETTINGS:
            return ''
        if key == 'connect':
            return 'no' if text.lower() in ('no', '0', 'false') else 'yes'
        if key in rules.STRING_SETTINGS:
            return text.split()[0] if text else ''
        try:
            value = float(text)
        except ValueError:
            value = 0.0
        if key in rules.INT_SETTINGS:
            return int(value)
        if key == 'avoid_overlap':
            return 'yes' if value == 1 else 'no'
        return value

    def sidemax(self):
        return float(self._additional_value('sidemax'))

    def _short_tol_toggled(self, on):
        self.shortTolDist.setEnabled(on)
        self.shortTolScale.setEnabled(on)

    # ====================================================== tree helpers ==
    def _selected(self):
        items = self.tree.selectedItems()
        return items[0] if items else None

    def _new_item(self, parent, text, index=None, color=None):
        item = QtWidgets.QTreeWidgetItem([text])
        item.setFlags(item.flags() | Qt.ItemIsEditable)
        # set before insertion so no itemChanged signal is emitted
        self._next_key += 1
        item.setData(0, Qt.UserRole, self._next_key)
        if color:
            item.setData(0, COLOR_ROLE, color)
        self._set_swatch(item)
        if parent is None:
            if index is None:
                self.tree.addTopLevelItem(item)
            else:
                self.tree.insertTopLevelItem(index, item)
        else:
            if index is None:
                parent.addChild(item)
            else:
                parent.insertChild(index, item)
            parent.setExpanded(True)
        return item

    def _index(self, item):
        parent = item.parent()
        if parent is None:
            return self.tree.indexOfTopLevelItem(item)
        return parent.indexOfChild(item)

    def _take(self, item):
        parent = item.parent()
        if parent is None:
            return self.tree.takeTopLevelItem(self.tree.indexOfTopLevelItem(item))
        return parent.takeChild(parent.indexOfChild(item))

    def _children(self, parent):
        if parent is None:
            return [self.tree.topLevelItem(i) for i in range(self.tree.topLevelItemCount())]
        return [parent.child(i) for i in range(parent.childCount())]

    def _descendants(self, parent):
        out = []
        for ch in self._children(parent):
            out.append(ch)
            out.extend(self._descendants(ch))
        return out

    def _structure_parent(self, item):
        """Parent for keywords that belong to a structure block (None = top level)."""
        if _head(item) == 'structure':
            return item
        parent = item.parent()
        if parent is not None and _head(parent) == 'atoms':
            return parent.parent()
        return parent

    def _block_parent(self, item):
        """Parent for keywords that belong to a structure or atoms block."""
        if _head(item) in rules.BLOCKS:
            return item
        return item.parent()

    def to_nodes(self, parent=None, colors=False):
        """Rule tree of the items; colors=True appends the own constraint colors as comments."""
        out = []
        for it in self._children(parent):
            text = it.text(0)
            own = it.data(0, COLOR_ROLE)
            if colors and own:
                text = text.rstrip() + '  ' + rules.color_comment(own)
            out.append([text, self.to_nodes(it, colors)])
        return out

    def _from_nodes(self, nodes, parent=None, random_colors=False):
        """
        Build tree items from nodes. A '# color #rrggbb' comment gives the rule its color;
        random_colors gives the other constraints a random one.
        """
        for text, children in nodes:
            text, color = rules.pop_color(text)
            if color is None and random_colors and _is_constraint(text):
                color = random_color()
            item = self._new_item(parent, text, color=color)
            self._from_nodes(children, item, random_colors)

    def _set_swatch(self, item):
        """Color swatch in the tree: own color, the default CGO color for constraints, else none."""
        own = item.data(0, COLOR_ROLE)
        if own:
            color = QtGui.QColor(own)
        elif _is_constraint(item.text(0)):
            color = QtGui.QColor(self.dlg.settings['cgoColor'])
        else:
            color = None
        tree = item.treeWidget()
        blocked = tree.blockSignals(True) if tree is not None else False
        try:
            item.setData(0, Qt.DecorationRole, color)
            item.setToolTip(0, 'Double-click the color swatch to change the CGO color'
                            if color is not None else '')
        finally:
            if tree is not None:
                tree.blockSignals(blocked)

    def _on_swatch(self, item):
        """True if the mouse is on the color swatch of the item."""
        if item.data(0, Qt.DecorationRole) is None:
            return False
        rect = self.tree.visualItemRect(item)
        size = self.tree.iconSize()
        width = size.width() if size.isValid() and size.width() > 0 else \
            self.tree.style().pixelMetric(QtWidgets.QStyle.PM_SmallIconSize)
        pos = self.tree.viewport().mapFromGlobal(QtGui.QCursor.pos())
        return rect.left() <= pos.x() <= rect.left() + width + 6

    def _item_double_clicked(self, item, column):
        if not self._on_swatch(item):
            self.tree.editItem(item, 0)
            return
        start = item.data(0, COLOR_ROLE) or self.dlg.settings['cgoColor']
        color = QtWidgets.QColorDialog.getColor(QtGui.QColor(start), self, 'Choose CGO color')
        if not color.isValid():  # cancelled
            return
        self._set_item_color(item, color.name())

    def _set_item_color(self, item, color):
        blocked = self.tree.blockSignals(True)  # setData would emit itemChanged
        try:
            item.setData(0, COLOR_ROLE, color)
        finally:
            self.tree.blockSignals(blocked)
        self._set_swatch(item)
        self.tree.setCurrentItem(item)
        self._edit_vis = True
        try:
            self.tree_op('visSel')
        finally:
            self._edit_vis = False

    # ======================================================== add rules ==
    def add_rule(self, element):
        settings = self.dlg.settings
        top = None  # insert at top level
        if element == 'basic':
            self._new_item(top, 'tolerance ' + str(self.tolerance.value()))
            self._new_item(top, 'filetype ' + self.filetype.currentText())
            self._new_item(top, 'seed ' + str(self.seed.value()))
            if not self.output.text():  # if empty, set some default output file name
                self.output.setText('output.' + self.filetype.currentText())
            self._new_item(top, 'output ' + self._rule_path(self.output.text(), 'output'))
            return
        if element == 'check':
            self._new_item(top, 'check')
            return
        if element == 'structure':
            if not self.structure.text():
                self.structure.setText('input.pdb')
            fn = self._rule_path(self.structure.text(), 'structure')
            self._new_item(top, 'structure ' + fn)
            if os.path.isfile(fn):
                if settings['loadStruct']:
                    cmd.load(fn, unique_name(os.path.splitext(os.path.basename(fn))[0]))
            else:
                info(self, 'Warning', 'Warning, file to read structure from: ' + fn +
                     ' not found! \n Rule has been added anyway.')
            return
        if element == 'additionalUST':
            if self.useShortTol.isChecked():
                self._new_item(top, 'use_short_tol')
                self._new_item(top, 'short_tol_dist ' + str(self.shortTolDist.value()))
                self._new_item(top, 'short_tol_scale ' + str(self.shortTolScale.value()))
            return

        if element == 'constraint' and self.direction.currentText() in rules.GLOBAL_DIRECTIONS:
            self._add_pbc()
            return

        key = self.additionalMenu.currentText()
        if element == 'additional' and key in rules.STRING_SETTINGS and not self._additional_value(key):
            info(self, 'Warning', key + ' needs a value. Rule not added.')
            return
        if element == 'additional' and key not in rules.STRUCTURE_SETTINGS + rules.ATOMS_SETTINGS:
            self._new_item(top, (key + ' ' + str(self._additional_value(key))).strip())
            return

        sel = self._selected()
        if element in ('restartTo', 'restartFrom'):
            parent = self._structure_parent(sel) if sel is not None else None
            if element == 'restartTo':
                self._new_item(parent, 'restart_to ' + self._rule_path(self.restartTo.text(), 'restart'))
            else:
                fn = self._rule_path(self.restartFrom.text(), 'restart')
                if os.path.isfile(fn):
                    with open(fn, 'rt') as f:
                        msg = 'Restart File has ' + str(len(f.readlines())) + ' molecular elements.'
                else:
                    msg = 'Warning, file to restart from: ' + fn + ' not found! \n Rule has been added anyway.'
                info(self, 'Information', msg)
                self._new_item(parent, 'restart_from ' + fn)
            return

        # elements that need a selection to be inserted
        if sel is None:
            print('Nothing selected')
            return
        if element in ('number', 'atoms', 'relabel', 'constrain_rotation') or \
                (element == 'additional' and key in rules.STRUCTURE_SETTINGS):
            parent = self._structure_parent(sel)
        else:  # constraint, center, radius, short_radius*
            parent = self._block_parent(sel)
        if parent is None:
            print('Select a structure (or atoms) block first')
            return

        if element == 'number':
            self._new_item(parent, 'number ' + str(self.number.value()))
        elif element == 'atoms':
            self._new_item(parent, 'atoms ' + self.atoms.text().strip())
        elif element == 'constraint':
            region = self.region.currentText()
            lin = self.direction.currentText() + ' ' + region
            for i in range(rules.REGION_NPARAMS[region]):
                lin += ' ' + str(self.params[i].value())
            newItem = self._new_item(parent, lin.strip(), color=self.cgoColor.color())
            if settings['loadVisual']:  # visualize while inserting
                self.tree.setCurrentItem(newItem)
                self.tree_op('visSel')
        elif element == 'center':
            self._new_item(parent, self.centerKind.currentText())
        elif element == 'radius':
            self._new_item(parent, 'radius ' + str(self.atomRadius.value()))
        elif element == 'relabel':  # resID and changeChains
            # Packmol turns changechains on whenever the keyword is present, whatever its value
            if self.changeChains.isChecked():
                self._new_item(parent, 'changechains')
            self._new_item(parent, 'resnumbers ' + self.resNumber.currentText())
        elif element == 'constrain_rotation':
            self._new_item(parent, 'constrain_rotation %s %6.3f %6.3f' %
                           (self.rotAxis.currentText(), self.rotAngle.value(), self.rotRange.value()))
        elif element == 'additional':
            self._new_item(parent, key + ' ' + str(self._additional_value(key)))

    def _add_pbc(self):
        """pbc is global: always inserted at top level, with 3 (max only) or 6 values."""
        v = [p.value() for p in self.params[:6]]
        if v[3:] == [0., 0., 0.]:
            values = v[:3]
        elif all(v[i + 3] > v[i] for i in range(3)):
            values = v
        else:
            error(self, 'Error', 'pbc: x1, y1, z1 must be greater than x0, y0, z0 (or all 0 to give\n'
                  'only the box maximum x0, y0, z0).')
            return
        newItem = self._new_item(None, 'pbc ' + ' '.join(str(x) for x in values),
                                 color=self.cgoColor.color())
        if self.dlg.settings['loadVisual']:
            self.tree.setCurrentItem(newItem)
            self.tree_op('visSel')

    # =================================================== tree operations ==
    def _item_edited(self, item, column):
        # a typed '# color #rrggbb' comment sets the rule color
        text, color = rules.pop_color(item.text(0))
        blocked = self.tree.blockSignals(True)
        try:
            if color:
                item.setText(0, text)
                item.setData(0, COLOR_ROLE, color)
        finally:
            self.tree.blockSignals(blocked)
        self._set_swatch(item)
        if self.tree.currentItem() is item:
            self._edit_vis = True
            try:
                self.tree_op('visSel')
                self._redraw_dependent_planes(item)
            finally:
                self._edit_vis = False

    @staticmethod
    def _key(item):
        return item.data(0, Qt.UserRole)

    def _delete_cgo(self, key):
        for o in self.cgo_loaded.pop(key, []):
            cmd.delete(o)

    def tree_op(self, op=''):
        sel = self._selected()
        if op == 'blank':
            parent = sel.parent() if sel is not None else None
            pos = self._index(sel) if sel is not None else 0
            newItem = self._new_item(parent, '', pos)
            self.tree.setCurrentItem(newItem)
            self.tree.editItem(newItem, 0)
        elif op in ('moveUp', 'moveDown'):
            if sel is None:
                return
            parent = sel.parent()
            idx = self._index(sel)
            count = len(self._children(parent))
            new = idx - 1 if op == 'moveUp' else idx + 1
            if 0 <= new < count:
                expanded = sel.isExpanded()
                self._take(sel)
                if parent is None:
                    self.tree.insertTopLevelItem(new, sel)
                else:
                    parent.insertChild(new, sel)
                sel.setExpanded(expanded)
                self.tree.setCurrentItem(sel)
        elif op == 'clear':
            self.tree_op('clearVis')
            self.tree.clear()
        elif op == 'remove':
            if self.tree.topLevelItemCount() == 0:
                print('Nothing to remove.')
            elif sel is None:
                print('No item selected.')
            else:
                for it in [sel] + self._descendants(sel):
                    self._delete_cgo(self._key(it))
                self._take(sel)
        elif op == 'check':
            self.dlg.set_status('check')
            settings = self.dlg.settings
            try:
                self.checkText.write('Volume method: %s\n' % settings['volMethod'])
                check_volumes(self.to_nodes(), self.sidemax(), settings['lenTolerance'],
                              self.checkText.write, settings.mol_volume)
            except InputError as e:
                self.checkText.write('MoloVol error: %s\nCan not check further!\n\n' % e)
                error(self, 'Error', 'MoloVol: %s' % e)
            except MemoryError:
                msg = 'MoloVol: out of memory, increase the grid resolution in Plugin Settings.'
                self.checkText.write(msg + '\nCan not check further!\n\n')
                error(self, 'Error', msg)
            finally:
                self.dlg.set_status('stdby')
        elif op == 'visSel':
            if sel is None:
                print('No item selected.')
                return
            if _head(sel) in rules.GLOBAL_DIRECTIONS:
                self.vis([sel])
                return
            # climb up to the structure block
            unit = _head(sel)
            parent = sel.parent()
            selitems = [sel]
            while unit != 'structure' and parent is not None:
                unit = _head(parent)
                selitems = [parent]
                parent = parent.parent()
            self.constraint_list = self._constraints_before(sel)
            self.new_const = True
            extraitems = self._descendants(sel) or [sel]
            self.vis(selitems + extraitems)
            self.new_const = False
        elif op == 'visAll':
            self.vis(self._descendants(None))
        elif op == 'clearVis':
            for key in list(self.cgo_loaded):
                self._delete_cgo(key)
        elif op == 'clearCheck':
            self.checkText.clear()

    def _constraints_before(self, item):
        """Volume constraints before item in its structure, in tree order, as 'visAll' sees them."""
        found = []
        for it in self._descendants(None):
            sub = it.text(0).split()
            is_struct = sub[:1] == ['structure']
            if is_struct:
                found = []
            if self._key(it) == self._key(item):
                break
            if not is_struct and _is_volume(sub):
                found.append(sub)
        return found

    def _dependent_planes(self, item):
        """Drawn planes sized by item: those after it, up to the next volume constraint or structure."""
        items = self._descendants(None)
        keys = [self._key(it) for it in items]
        planes = []
        for it in items[keys.index(self._key(item)) + 1:]:
            sub = it.text(0).split()
            if sub[:1] == ['structure'] or _is_volume(sub):
                break
            if len(sub) > 1 and sub[1].lower() == 'plane' and self._key(it) in self.cgo_loaded:
                planes.append(it)
        return planes

    def _redraw_dependent_planes(self, item):
        """Redraw the planes whose size comes from the edited rule."""
        if _head(item) == 'structure':  # redrawing a structure already redraws its rules
            return
        planes = self._dependent_planes(item)
        if planes:
            # only volume constraints change the size, so all these planes share the same list
            self.constraint_list = self._constraints_before(planes[0])
            self.vis(planes, self._structure_name(planes[0]))

    @staticmethod
    def _structure_name(item):
        """Prefix of the constraint objects drawn for item ('NN' without structure)."""
        while item is not None:
            sub = item.text(0).split()
            if len(sub) > 1 and sub[0] == 'structure':
                return os.path.splitext(os.path.basename(sub[1].strip()))[0]
            item = item.parent()
        return 'NN'

    def vis(self, items, structName='NN'):
        # structName: prefix for constraint objects without structure
        for it in items:
            # redraw in place: remove what was drawn for this rule, keep its names for reuse
            self._old_names = self.cgo_loaded.pop(self._key(it), [])
            for o in self._old_names:
                cmd.delete(o)
                self.dlg.vol_cgo.pop(o, None)
            self._set_swatch(it)  # follows the default CGO color
            subunits = it.text(0).split()
            if len(subunits) > 1:
                if subunits[0].lower() in rules.GLOBAL_DIRECTIONS:
                    self._vis_pbc(it, subunits)
                elif subunits[0] == 'structure':
                    if not self.new_const:
                        self.constraint_list = []
                    fn = subunits[1].strip()
                    structName = os.path.splitext(os.path.basename(fn))[0]
                    if self.dlg.settings['loadStruct']:
                        if not os.path.isfile(fn):
                            info(self, 'Warning', 'Structure file ' + fn + ' not found !')
                        elif structName not in cmd.get_names('objects'):
                            cmd.load(fn, structName)
                        elif not self._edit_vis and ask_yes_no(self, 'Warning', 'Structure object with that name already loaded '
                                        'into PyMOL\n Do you want to reload?'):
                            cmd.delete(structName)
                            cmd.load(fn, structName)
                elif subunits[1].lower() in CGO_DICT:
                    self._vis_constraint(it, subunits, structName)
            elif subunits and subunits[0] == 'structure':  # kw without argument
                error(self, 'Error', '"structure" keyword do not have argument !')
                if not self.new_const:
                    self.constraint_list = []
                structName = 'NN'

    def _item_color(self, it):
        own = it.data(0, COLOR_ROLE)
        return color_rgb(own) if own else self.dlg.settings.cgo_color()

    def _load_cgo(self, it, cgo_pymol, base, element):
        # same structure and region as before: keep the old name (water_box, water_box1, ...)
        same = [o for o in self._old_names if not is_label(o) and re.match(re.escape(base) + r'\d*$', o)]
        cgoName = same[0] if same and same[0] not in cmd.get_names('objects') else unique_name(base, sep='')
        cmd.load_cgo(cgo_pymol, cgoName)
        loaded = self.cgo_loaded.setdefault(self._key(it), [])
        loaded.append(cgoName)
        self.dlg.vol_cgo[cgoName] = volCGO(element)
        # hidden pseudoatom labelled with the CGO name, 5 A away from the other labels if possible
        label = add_cgo_label(cgoName, cgo_pymol, self._item_color(it), self._label_coords())
        if label:
            loaded.append(label)

    def _label_coords(self):
        """Positions of the constraint labels drawn by this tab."""
        names = set(cmd.get_names('objects'))
        coords = [cmd.get_coords(o) for objs in self.cgo_loaded.values() for o in objs
                  if is_label(o) and o in names]
        coords = [c for c in coords if c is not None]
        return np.concatenate(coords) if coords else np.zeros((0, 3))

    def _vis_pbc(self, it, subunits):
        """pbc is drawn as a box (not added to constraint_list: it is not a structure constraint)."""
        params = _pbc_box(subunits)
        if params is None:
            print('Wrong parameters in rule: ' + it.text(0))
            return
        self._load_cgo(it, CGO_DICT['box'](params, self._item_color(it)), 'pbc_box', ['box'] + params)

    def _vis_constraint(self, it, subunits, structName):
        kind = subunits[1].lower()
        color = self._item_color(it)
        if not rules.valid_region(kind, subunits[2:]):
            print('Wrong parameters in rule: ' + it.text(0))
            return
        cgo_packmol = [float(par) for par in subunits[2:]]
        if kind == 'plane':
            # limit plane size to the previous constraint, default 0-50 A
            prismPlane = []
            if self.constraint_list:
                last = self.constraint_list[-1]
                if last[1].lower() == 'cube':
                    prismPlane = np.float64(last[2:5])
                    prismPlane = np.append(prismPlane, np.float64(last[2:5]) + np.float64(last[-1]))
                elif last[1].lower() == 'box':
                    prismPlane = np.float64(last[2:])
                elif last[1].lower() == 'sphere':
                    sphere = np.float64(last[2:])
                    prismPlane = np.append(sphere[:3] - sphere[-1], sphere[:3] + sphere[-1])
            min1, max1, min2, max2 = 0., 50., 0., 50.
            if len(prismPlane) != 0:
                vec = cgo_packmol[:3]
                if vec == [1., 0., 0.]:
                    min1, max1, min2, max2 = prismPlane[2], prismPlane[5], prismPlane[1], prismPlane[4]
                elif vec == [0., 1., 0.]:
                    min1, max1, min2, max2 = prismPlane[0], prismPlane[3], prismPlane[2], prismPlane[5]
                elif vec == [0., 0., 1.]:
                    min1, max1, min2, max2 = prismPlane[0], prismPlane[3], prismPlane[1], prismPlane[4]
            cgo_pymol = CGO_DICT[kind](cgo_packmol, min1, max1, min2, max2, color)
        else:
            if kind not in SURFACES:
                self.constraint_list.append(subunits)
            cgo_pymol = CGO_DICT[kind](cgo_packmol, color)
        self._load_cgo(it, cgo_pymol, structName + '_' + subunits[1], subunits[1:])

    # ============================================================ files ==
    def load_input(self):
        fn = ask_open(self, 'Open Input file', 'input')
        if not fn:
            return
        self.input_file.setText(fn)
        self.read_input(fn)

    def load_typed(self):
        """Read the input file typed in the entry, or ask for one if it does not exist."""
        fn = self.input_file.text().strip()
        if fn and os.path.isfile(fn):
            self.read_input(fn)
        else:
            self.load_input()

    def save_input(self):
        fn = ask_save(self, 'Save Input file', 'input')
        if not fn:
            return
        self.input_file.setText(fn)
        rules.write_inp(fn, self.to_nodes(colors=True))

    def read_input(self, fn):
        settings = self.dlg.settings
        nodes = rules.read_inp(fn)
        self.tree_op('clearVis')
        self.tree.clear()
        self._from_nodes(nodes, random_colors=True)
        for it in self._descendants(None):
            sub = it.text(0).split()
            if len(sub) < 2:
                continue
            if sub[0].lower() == 'output':
                self.output.setText(sub[1])
            # with visualization on, vis() loads the structures
            elif sub[0].lower() == 'structure' and settings['loadStruct'] and not settings['loadVisual']:
                if os.path.isfile(sub[1]):
                    cmd.load(sub[1], unique_name(os.path.splitext(os.path.basename(sub[1]))[0]))
                else:
                    info(self, 'Warning', 'Structure file ' + sub[1] + ' not found !')
        if settings['loadVisual']:
            self.tree_op('visAll')
