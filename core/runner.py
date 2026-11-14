"""Non-blocking Packmol execution with QProcess (repetitions, warnings, output loading)."""
import os

from pymol import cmd
from pymol.Qt import QtCore

from .pymolutil import unique_name

Signal = getattr(QtCore, 'Signal', None) or QtCore.pyqtSignal


class PackmolRunner(QtCore.QObject):
    """
    Runs Packmol nconfig times. settings is a mapping with the plugin settings
    ('loadPartial', 'loadOutput', 'packWarningN'). ask_continue(title, message) -> bool
    is called when a packing stops and the user may skip to the next repetition.
    """
    log = Signal(str)
    progress = Signal(int)  # 0-100
    status = Signal(str)  # status key: 'run' or 'stdby'
    runningChanged = Signal(bool)

    def __init__(self, settings, ask_continue, parent=None):
        QtCore.QObject.__init__(self, parent)
        self.settings = settings
        self.ask_continue = ask_continue
        self.process = None
        self.nconfig = 1
        self._reset()

    def _reset(self):
        self.killed = False
        self.n = 0
        self.warningN = 0
        self.normalTermination = False
        self.override = False  # flag to pack although stopped packing detected

    def is_running(self):
        return self.process is not None

    def start(self, packmol, input_file, output_file, nconfig):
        if self.is_running():
            return
        self.packmol = packmol
        self.input_file = input_file
        self.output_file = output_file
        self.nconfig = max(1, int(nconfig))
        self._reset()
        self.progress.emit(0)
        self.status.emit('run')
        self.runningChanged.emit(True)
        self._start_one()

    def kill(self):
        if self.is_running():
            self.killed = True
            self.process.kill()  # finished() follows and closes the run
        else:
            self._done()

    # -- process handling ---------------------------------------------------
    def _start_one(self):
        self._buffer = b''
        p = QtCore.QProcess(self)
        p.setWorkingDirectory(os.getcwd())
        p.setStandardInputFile(self.input_file)
        p.readyReadStandardOutput.connect(self._read)
        p.finished.connect(self._finished)
        if hasattr(p, 'errorOccurred'):
            p.errorOccurred.connect(self._error)
        self.process = p
        p.start(self.packmol, [])

    def _error(self, err):
        if err == QtCore.QProcess.FailedToStart:  # finished() is not emitted in this case
            self.log.emit('Could not start Packmol: %s\n' % self.packmol)
            self.process = None
            self._done()

    def _read(self):
        self._buffer += bytes(self.process.readAllStandardOutput())
        lines = self._buffer.split(b'\n')
        self._buffer = lines.pop()
        for raw in lines:
            self._parse(raw.decode('utf-8', 'replace') + '\n')

    def _parse(self, line):
        if line.find('  Packing ') != -1 or line.find('  Current solution') != -1:  # partial text status
            self.log.emit(line)
        if self.settings['loadPartial'] and line.find('  Current solution') != -1:  # load partial solution
            partialFn = line[line.find(':') + 1:].strip()
            if os.path.isfile(partialFn):
                cmd.load(partialFn, 'partial')
        if line.find('Success!') != -1:
            self.normalTermination = True  # only change flag if succeed !
            self.log.emit(line)
        if line.find('ERROR') != -1:
            self.log.emit(line)
            if line.find('Could not open file.') != -1:  # structure file not found
                self.log.emit('Structure kw refers to an inexistent file! \n\n')
                self.n = self.nconfig
                self.warningN = self.settings['packWarningN']
            elif line.find('Packmol was unable to put the molecules') != -1:
                self.log.emit('Initial constraints too tight ! \n\n')
        if line.find('STOP') != -1:
            self.log.emit(line)

    def _finished(self, *args):
        if self._buffer:
            self._parse(self._buffer.decode('utf-8', 'replace') + '\n')
            self._buffer = b''
        self.process = None
        self.log.emit('Finished output capture: \n  Normal Termination: ' + str(self.normalTermination) +
                      '\n  Killed: ' + str(self.killed) +
                      '\n  Rep(N) & Warn(N): ' + str(self.n + 1) + ' - ' + str(self.warningN) + '\n')
        if self.killed:
            self.log.emit('\nKilled\n')
            self._done()
            return
        self.log.emit('Finished\n')

        warnLimit = self.settings['packWarningN']
        if not self.normalTermination:
            self.log.emit('Packing was NOT successfully completed.\n')
            if self.nconfig > 1 and (self.n + 1) < self.nconfig and not self.override and self.warningN >= warnLimit:
                ans = self.ask_continue('Packing STOPPED',
                                        'Packing was NOT successfully completed. Proceed to next packing?')
                self.warningN = 0
                if not ans:
                    self.progress.emit(0)
                    self._done()
                    return
                self.override = True  # continue packing
                self.normalTermination = True
                self.log.emit('User override.\nCurrent solution could(not) be loaded.\n'
                              'The next packing attempt will be made.\n')
            elif self.warningN >= warnLimit:
                self._done()
                return
            else:
                self.warningN += 1  # retry this repetition

        if self.normalTermination:
            self._collect_output()
            if self.n >= self.nconfig:
                self._done()
                return
            self.status.emit('run')
            self.log.emit('Running next packing !\n')
        self._start_one()

    def _collect_output(self):
        out = self.output_file
        structName = os.path.splitext(os.path.split(out)[1])[0]
        self.n += 1
        self.warningN = 0
        self.normalTermination = False
        self.override = False
        if os.path.isfile(out):
            if self.settings['loadOutput']:
                if out[-4:].find('xyz') != -1:  # PyMOL fails to read first line of xyz when not properly formatted
                    with open(out, 'r') as fi:
                        li = fi.readlines()
                    if li:
                        li[0] = li[0].strip(' ')
                        with open(out, 'w') as fo:
                            fo.writelines(li)
                structureName = unique_name(structName)
                cmd.load(out, structureName)
                self.log.emit('Output file loaded as %s.\n\n' % structureName)
            if self.nconfig > 1:  # only rename output file for multiple repetitions
                folder, fname = os.path.split(out)
                os.replace(out, os.path.join(folder, fname.replace(structName, structName + '_' + str(self.n), 1)))
        else:
            self.log.emit('Output file not found. Could not be loaded.\n\n')
        self.progress.emit(int(self.n * 100 / self.nconfig))

    def _done(self):
        self._reset()
        self.status.emit('stdby')
        self.runningChanged.emit(False)
