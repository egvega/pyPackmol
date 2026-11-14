"""'Pack' tab: run Packmol (non-blocking) with repetitions."""
import os
import shutil

from pymol.Qt import QtWidgets

from ..core.rules import output_from_inp, seed_from_inp
from ..core.runner import PackmolRunner
from .widgets import NumEdit, FileField, LogView, button, group, hbox, error, ask_yes_no


class PackTab(QtWidgets.QWidget):
    def __init__(self, dlg, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.dlg = dlg
        main = QtWidgets.QVBoxLayout(self)

        g, gl = group('Running Settings', QtWidgets.QHBoxLayout)
        main.addWidget(g)
        col = QtWidgets.QVBoxLayout()
        gl.addLayout(col, 1)
        self.packmolPath = FileField('...', 'Search for Packmol executable', 'open', 'any', button_after=True)
        location = shutil.which('packmol')
        if location:
            print('Packmol found in Path !')
            self.packmolPath.setText(location)
        col.addLayout(hbox('Packmol Path:', self.packmolPath))
        self.inputDisplay = QtWidgets.QLineEdit()
        self.inputDisplay.setReadOnly(True)
        col.addLayout(hbox('Input File:', self.inputDisplay))
        self.nconfig = NumEdit(1, int, 10)
        col.addLayout(hbox('Repetitions:', self.nconfig, None))

        self.buttonRun = button('Start \n Packing', self.run)
        self.buttonKill = button('Kill \n Packing', self.kill)
        self.buttonKill.setEnabled(False)
        gl.addWidget(self.buttonRun)
        gl.addWidget(self.buttonKill)
        gl.addWidget(button('Clear \n', self.reset))

        self.progress = QtWidgets.QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        main.addWidget(self.progress)

        g, gl = group('Running Status')
        main.addWidget(g, 1)
        self.statusText = LogView()
        gl.addWidget(self.statusText)

        self.runner = PackmolRunner(dlg.settings,
                                    lambda title, msg: ask_yes_no(self, title, msg), self)
        self.runner.log.connect(self.statusText.write)
        self.runner.progress.connect(self.progress.setValue)
        self.runner.status.connect(dlg.set_status)
        self.runner.runningChanged.connect(self._running_changed)

    def _running_changed(self, running):
        self.buttonRun.setEnabled(not running)
        self.buttonKill.setEnabled(running)

    def run(self):
        packmol = self.packmolPath.text()
        inp = self.inputDisplay.text().strip()
        if not os.path.isfile(packmol):
            errMsg = 'Packmol not found !'
        elif not os.path.isfile(inp):
            errMsg = 'Input file not found !'
        else:
            output = output_from_inp(inp)
            if output:
                # Packmol's default seed (no 'seed' line) is fixed too; only -1 randomizes each run
                if self.nconfig.value() > 1 and seed_from_inp(inp) != -1:
                    ans = QtWidgets.QMessageBox.question(
                        self, 'Warning',
                        'Warning: a fixed seed was detected and multiple runs were requested. '
                        'Only one run will be performed. Proceed?',
                        QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.Cancel)
                    if ans != QtWidgets.QMessageBox.Yes:
                        return
                    self.nconfig.setValue(1)
                self.runner.start(packmol, inp, output, self.nconfig.value())
                return
            errMsg = '"output" keyword Not or Wrong defined in input file !'
        error(self, 'Error', errMsg)

    def kill(self):
        self.runner.kill()

    def reset(self):
        if self.runner.is_running():
            self.runner.kill()
        self.statusText.clear()
        self.progress.setValue(0)
