"""'Restart Files' tab: rebuild a packed system from templates and a Packmol restart file."""
import os

from pymol import cmd
from pymol.Qt import QtWidgets

from ..core import rules
from ..core.restart import iter_restart, count_positions, TEMPLATE_EXPR
from .widgets import FileField, LogView, button, group, hbox


def _template_atoms(fn):
    """Atoms of the first state of a structure file as TEMPLATE_EXPR rows."""
    templateName = cmd.get_unused_name('_packmol_template')
    cmd.load(fn, templateName)
    template = []
    cmd.iterate_state(1, templateName, 'template.append(%s)' % TEMPLATE_EXPR,
                      space={'template': template})
    cmd.delete(templateName)
    return template


class RestartTab(QtWidgets.QWidget):
    def __init__(self, dlg, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.dlg = dlg
        main = QtWidgets.QVBoxLayout(self)

        g, gl = group('Restart File')
        main.addWidget(g)
        # whole system: structures and numbers from a Packmol input file
        # partial: one structure template for every restart position
        self.whole = QtWidgets.QRadioButton('Whole')
        self.partial = QtWidgets.QRadioButton('Partial')
        self.partial.setChecked(True)
        self.inputFile = FileField('Packmol Input', 'Open Input file', 'open', 'input')
        self.template = FileField('Structure Template', 'Open Structure file', 'open', 'structure')
        self.restart = FileField('Packmol Restart File', 'Open Restart file', 'open', 'restart')
        grid = QtWidgets.QGridLayout()
        grid.addWidget(self.whole, 0, 0)
        grid.addWidget(self.inputFile, 0, 1)
        grid.addWidget(self.partial, 1, 0)
        grid.addWidget(self.template, 1, 1)
        grid.addWidget(self.restart, 2, 0, 1, 2)
        gl.addLayout(hbox(grid, button('Load', self.load), button('Clear', self.clear)))
        self.whole.toggled.connect(self._mode_changed)
        self._mode_changed()

        self.progress = QtWidgets.QProgressBar()
        self.progress.setValue(0)
        main.addWidget(self.progress)

        g, gl = group('Restart Log')
        main.addWidget(g, 1)
        self.log = LogView()
        gl.addWidget(self.log)

    def _mode_changed(self, *args):
        # only the file field of the chosen mode can be used
        self.inputFile.setEnabled(self.whole.isChecked())
        self.template.setEnabled(self.partial.isChecked())

    def clear(self):
        self.progress.setValue(0)
        self.log.clear()
        self.inputFile.setText('')
        self.template.setText('')
        self.restart.setText('')
        self.partial.setChecked(True)  # default mode; toggled updates the fields

    def _structures(self, npos):
        """[(structure file, number), ...] to rebuild, or None (error logged)."""
        if self.whole.isChecked():
            inputFn = self.inputFile.text()
            if not os.path.isfile(inputFn):
                self.log.write('ERROR: Input File to whole system restart not found \n')
                return None
            structs = rules.structure_counts(rules.read_inp(inputFn))
            for fn, number in structs:
                self.log.write('%s %d\n' % (fn, number))
            missing = [fn for fn, _ in structs if not os.path.isfile(fn)]
            if missing:
                self.log.write('ERROR: Structure file %s not found \n' % missing[0])
                return None
            if sum(number for _, number in structs) != npos:
                self.log.write('ERROR: Number of Molecules to build Restart is different '
                               'from the number of positions \n')
                return None
            return structs
        templateFn = self.template.text()
        if not os.path.isfile(templateFn):
            self.log.write('ERROR: Structure Template file not found \n')
            return None
        return [(templateFn, npos)]

    def load(self):
        restartFn = self.restart.text()
        if not os.path.isfile(restartFn):
            self.log.write('ERROR: Restart file not found \n')
            return
        with open(restartFn, 'rt') as f:
            lines = f.readlines()
        npos = count_positions(lines)
        self.log.write('Restart file has %d coordinates to be processed\n' % npos)
        structs = self._structures(npos)
        if structs is None:
            return
        self.dlg.set_status('restart')
        try:
            templates = [(_template_atoms(fn), number) for fn, number in structs]
            for (fn, _), (atoms, _) in zip(structs, templates):
                if not atoms:
                    self.log.write('ERROR: no atoms read from %s \n' % fn)
                    return
            self.progress.setMaximum(max(1, npos))
            objRestartName = cmd.get_unused_name(os.path.splitext(os.path.basename(restartFn))[0])
            showBaricenter = self.dlg.settings['showBaricenter']
            pdb = []
            done = 0
            for n, transformation, trmat, shift, text in iter_restart(templates, lines):
                self.log.write('mol %d\n' % n)
                self.log.write('packmol rst: ' + ('%.3f,' * 6) % tuple(transformation) + '\n')
                if showBaricenter:
                    cmd.pseudoatom('baricenter_' + objRestartName, pos=list(transformation[:3]), elem='H', name='X')
                self.log.write('rotation matrix\n')
                for r in trmat:
                    self.log.write(('%.3f,' * 3) % tuple(r) + '\n')
                self.log.write('final translation to fit baricenter ' + ('%.3f,' * 3) % tuple(shift) + '\n')
                pdb.append(text)
                done += 1
                self.progress.setValue(done)
                QtWidgets.QApplication.processEvents()
            self.log.write('EOF \n')
            cmd.read_pdbstr(''.join(pdb), objRestartName)
            cmd.zoom()
        finally:
            self.dlg.set_status('stdby')
