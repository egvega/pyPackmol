"""'Volume Calculator' tab: molecular and constraint volumes, packing estimates."""
import numpy as np
from pymol import cmd
from pymol.Qt import QtWidgets

from ..core.pymolutil import AXIS_NAME, is_label
from ..molovol import InputError
from .widgets import NumEdit, LogView, HoverListWidget, button, group, hbox, ask_save, error


class VolumeTab(QtWidgets.QWidget):
    def __init__(self, dlg, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.dlg = dlg
        main = QtWidgets.QHBoxLayout(self)
        left = QtWidgets.QVBoxLayout()
        main.addLayout(left, 1)

        g, gl = group('Molecules')
        left.addWidget(g)
        self.listMol = HoverListWidget()
        self.listMol.on_enter = self.update_mol_list
        self.listMol.currentTextChanged.connect(lambda t: t and self.molToPack.setText(t))
        self.listMol.itemDoubleClicked.connect(lambda it: self.check_mol(it.text()))
        gl.addWidget(self.listMol)
        self.Nmol = NumEdit(1, int, 5)
        gl.addLayout(hbox('Nmol', self.Nmol, button('Sphere', lambda: self.guess('sphere')),
                          button('Cube', lambda: self.guess('cube')), None))

        g, gl = group('Constraints (CGO elements)')
        left.addWidget(g)
        self.listCGO = HoverListWidget()
        self.listCGO.on_enter = self.update_cgo_list
        self.listCGO.itemDoubleClicked.connect(lambda it: self.check_cgo(it.text()))
        gl.addWidget(self.listCGO)
        self.molToPack = QtWidgets.QLineEdit()
        self.molToPack.setReadOnly(True)
        gl.addLayout(hbox(button('try', self.guess_nmol), self.molToPack))

        g, gl = group('History')
        main.addWidget(g, 1)
        self.hisText = LogView()
        gl.addWidget(self.hisText)
        gl.addLayout(hbox(None, button('Clear', self.hisText.clear), button('Save', self.save_log)))

    # -------------------------------------------------------------------
    @staticmethod
    def _set_list(listw, names):
        current = listw.currentItem().text() if listw.currentItem() else None
        listw.blockSignals(True)
        listw.clear()
        listw.addItems(names)
        if current in names:
            listw.setCurrentRow(names.index(current))
        listw.blockSignals(False)

    def update_mol_list(self):
        self._set_list(self.listMol, [n for n in cmd.get_names_of_type('object:molecule') if not is_label(n)])

    def update_cgo_list(self):
        self._set_list(self.listCGO, cmd.get_names_of_type('object:cgo'))

    def refresh(self):
        self.update_mol_list()
        self.update_cgo_list()

    def _selected(self, listw):
        it = listw.currentItem()
        return it.text() if it else ''

    # -------------------------------------------------------------------
    def check_mol(self, name, quiet=False):
        """Volume and length of a molecule into dlg.vol_mol; False if it could not be computed."""
        if not name:
            print('No molecule selected')
            return False
        self.dlg.set_status('check')
        try:
            volMol = self.dlg.settings.mol_volume(name)
            bBox = np.array(cmd.get_extent(name))  # [[x0,y0,z0],[x1,y1,z1]]
            dx, dy, dz = np.abs(bBox[0] - bBox[1])
            lenMol = np.sqrt(dx * dx + dy * dy + dz * dz)
            self.dlg.vol_mol[name] = (volMol, lenMol)
        except InputError as e:
            error(self, 'Error', 'MoloVol: %s' % e)
            return False
        except MemoryError:
            error(self, 'Error', 'MoloVol: out of memory, increase the grid resolution in Plugin Settings.')
            return False
        finally:
            self.dlg.set_status('stdby')
        if not quiet:
            self.hisText.write("Molecule: %s \t volume (%s)= %.3f A3 \t length (1D)= %.3f A\n"
                               % (name, self.dlg.settings['volMethod'], volMol, lenMol))
            if self.Nmol.value() > 1:
                self.hisText.write("\t total volume for %d molecules= %.3f A3 \n"
                                   % (self.Nmol.value(), volMol * self.Nmol.value()))
        return True

    def _cgo_volume(self, name):
        if name in self.dlg.vol_cgo and name != AXIS_NAME:
            return self.dlg.vol_cgo[name]
        error(self, 'Error', 'Nothing to do with selected CGO element.\nCGO: ' + name +
              ' is either an axis or not present into plugin instance dictionary !')
        return None

    def check_cgo(self, name):
        vol = self._cgo_volume(name)
        if vol is not None:
            self.hisText.write("CGO: %s \t volume= %.3f A3 \t length (1D)= %.3f A\n" % (name, vol[0], vol[1]))

    def guess(self, shape):
        name = self._selected(self.listMol)
        if not name:
            print('No molecule selection')
            return
        if name not in self.dlg.vol_mol and not self.check_mol(name, quiet=True):
            return
        volMol, lenMol = self.dlg.vol_mol[name]
        n = self.Nmol.value()
        lenTol = self.dlg.settings['lenTolerance']
        if shape == 'sphere':
            r = ((n * volMol) * 3 / 4. / np.pi) ** (1 / 3.)
            r = max(r * 2, lenMol) / 2. + lenTol  # diameter or molLength, then divide to get radius
            line = ("To pack %d %s molecules, a sphere of radius %.3f A (length 1D of %.3f A)  should be needed. \n"
                    % (n, name, r, r * 2))
        else:
            r = (n * volMol) ** (1 / 3.)
            r = max(np.sqrt(3) * r, lenMol) / np.sqrt(3) + lenTol  # body diag from side or molLength
            line = ("To pack %d %s molecules, a cube with a side length of %.3f A (body diagonal %.3f A) "
                    "should be needed. \n" % (n, name, r, np.sqrt(3) * r))
        self.hisText.write(line)

    def guess_nmol(self):
        cgoName = self._selected(self.listCGO)
        mol = self.molToPack.text()
        if not cgoName:
            print('No CGO selection')
            return
        if not mol:
            print('No Mol selection')
            return
        vol = self._cgo_volume(cgoName)
        if vol is None:
            return
        volObj, lenObj = vol
        if mol not in self.dlg.vol_mol and not self.check_mol(mol, quiet=True):
            return
        volMol, lenMol = self.dlg.vol_mol[mol]
        lenTol = self.dlg.settings['lenTolerance']
        self.hisText.write("According to volume, in %s approximately %d %s molecules can be packed. \n"
                           % (cgoName, int(volObj / volMol), mol))
        if lenObj < lenMol + lenTol:
            self.hisText.write("\t...but, molecular length is larger than constraint length !\n"
                               "A minimum constraint length of %.3f A is required.\n" % (lenMol + lenTol))

    def save_log(self):
        fn = ask_save(self, 'Select Log file', 'log')
        if fn:
            with open(fn, 'wt') as f:
                f.write(self.hisText.toPlainText() + '\n')
