"""Reusable Qt widgets and helpers (replace the old Tk makeentry/callFileDialog code)."""
import os
import random

from pymol.Qt import QtCore, QtGui, QtWidgets

Qt = QtCore.Qt


def char_width(widget, n):
    fm = widget.fontMetrics()
    w = fm.horizontalAdvance('0') if hasattr(fm, 'horizontalAdvance') else fm.width('0')
    return w * (n + 2) + 8


class NumEdit(QtWidgets.QLineEdit):
    """Line edit holding a number. kind is float or int."""

    def __init__(self, value=0, kind=float, width=None, parent=None):
        QtWidgets.QLineEdit.__init__(self, parent)
        self.kind = kind
        if kind is int:
            validator = QtGui.QIntValidator(self)
        else:
            validator = QtGui.QDoubleValidator(self)
            validator.setNotation(QtGui.QDoubleValidator.ScientificNotation)
        validator.setLocale(QtCore.QLocale.c())  # always '.' as decimal separator
        self.setValidator(validator)
        if width:
            self.setMaximumWidth(char_width(self, width))
        self.setValue(value)

    def value(self):
        text = self.text().strip()
        try:
            return self.kind(text)
        except ValueError:
            try:
                return self.kind(float(text))
            except ValueError:
                return self.kind(0)

    def setValue(self, v):
        self.setText(str(self.kind(v)))


class LogView(QtWidgets.QPlainTextEdit):
    """Read-only text output; write() appends raw text (keeps the caller's newlines)."""

    def __init__(self, parent=None):
        QtWidgets.QPlainTextEdit.__init__(self, parent)
        self.setReadOnly(True)
        self.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)

    def write(self, text):
        self.moveCursor(QtGui.QTextCursor.End)
        self.insertPlainText(text)
        self.moveCursor(QtGui.QTextCursor.End)
        self.ensureCursorVisible()


class RefreshComboBox(QtWidgets.QComboBox):
    """Combo box that asks provider() for its items every time the popup opens."""

    def __init__(self, provider, parent=None):
        QtWidgets.QComboBox.__init__(self, parent)
        self.provider = provider
        self.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToContents)
        self.refresh()

    def refresh(self):
        current = self.currentText()
        items = list(self.provider())
        self.blockSignals(True)
        self.clear()
        self.addItems(items)
        if current in items:
            self.setCurrentIndex(items.index(current))
        self.blockSignals(False)

    def showPopup(self):
        self.refresh()
        QtWidgets.QComboBox.showPopup(self)


class HoverListWidget(QtWidgets.QListWidget):
    """List widget calling on_enter() when the mouse enters it."""
    on_enter = None

    def enterEvent(self, event):
        if self.on_enter:
            self.on_enter()
        QtWidgets.QListWidget.enterEvent(self, event)


# -- file dialogs -------------------------------------------------------------
FILTERS = {
    'input': ('Packmol input file (*.inp);;All files (*)', '.inp'),
    'output': ('PDB (*.pdb);;XYZ (*.xyz);;All files (*)', '.pdb'),
    'structure': ('PDB (*.pdb);;XYZ (*.xyz);;All files (*)', '.pdb'),
    'restart': ('PACK (*.pack);;All files (*)', '.pack'),
    'log': ('Text file (*.txt);;All files (*)', '.txt'),
    'image': ('PNG (*.png);;TIFF (*.tif);;JPEG (*.jpg);;All files (*)', '.png'),
    'any': ('All files (*)', ''),
}


def ask_open(parent, title, kind='any'):
    """Returns the chosen file name or '' when cancelled."""
    filt = FILTERS[kind][0]
    fn = QtWidgets.QFileDialog.getOpenFileName(parent, title, os.getcwd(), filt)
    return fn[0] if isinstance(fn, tuple) else fn


def ask_save(parent, title, kind='any'):
    filt, ext = FILTERS[kind]
    fn = QtWidgets.QFileDialog.getSaveFileName(parent, title, os.getcwd(), filt)
    fn = fn[0] if isinstance(fn, tuple) else fn
    if fn and ext and not os.path.splitext(fn)[1]:
        fn += ext
    return fn


def ask_directory(parent, title):
    return QtWidgets.QFileDialog.getExistingDirectory(parent, title, os.getcwd())


class FileField(QtWidgets.QWidget):
    """Button + line edit. The button opens a file dialog and fills the edit."""

    def __init__(self, button_text, title, mode='open', kind='any', button_after=False,
                 on_chosen=None, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.title, self.mode, self.kind, self.on_chosen = title, mode, kind, on_chosen
        lay = QtWidgets.QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.button = QtWidgets.QPushButton(button_text)
        self.button.clicked.connect(self.choose)
        self.edit = QtWidgets.QLineEdit()
        if button_after:
            lay.addWidget(self.edit)
            lay.addWidget(self.button)
        else:
            lay.addWidget(self.button)
            lay.addWidget(self.edit)

    def choose(self):
        if self.mode == 'save':
            fn = ask_save(self, self.title, self.kind)
        elif self.mode == 'dir':
            fn = ask_directory(self, self.title)
        else:
            fn = ask_open(self, self.title, self.kind)
        if fn:  # cancel keeps old value
            self.edit.setText(fn)
            if self.on_chosen:
                self.on_chosen(fn)

    def text(self):
        return self.edit.text().strip()

    def setText(self, text):
        self.edit.setText(text)


def save_widget_image(widget, fn, dpi=600):
    """
    Render widget (with its children) into an image file at the given dpi.
    The widget is repainted at the higher resolution (sharp text, not an upscaled
    screen grab) and covered windows do not get into the image.
    Returns True on success.
    """
    scale = dpi / float(widget.logicalDpiX())
    w, h = widget.width(), widget.height()
    img = QtGui.QImage(QtCore.QSize(int(round(w * scale)), int(round(h * scale))),
                       QtGui.QImage.Format_ARGB32)
    img.fill(Qt.white)
    dpm = int(round(dpi / 0.0254))  # dots per metre, stored in the file
    img.setDotsPerMeterX(dpm)
    img.setDotsPerMeterY(dpm)
    painter = QtGui.QPainter(img)
    painter.setRenderHints(QtGui.QPainter.Antialiasing | QtGui.QPainter.TextAntialiasing |
                           QtGui.QPainter.SmoothPixmapTransform)
    painter.scale(scale, scale)
    widget.render(painter, QtCore.QPoint(), QtGui.QRegion(), QtWidgets.QWidget.DrawChildren)
    painter.end()
    return img.save(fn)


class ClickableLabel(QtWidgets.QLabel):
    """Label calling on_double_click() when double-clicked."""
    on_double_click = None

    def mouseDoubleClickEvent(self, event):
        if self.on_double_click:
            self.on_double_click()
        QtWidgets.QLabel.mouseDoubleClickEvent(self, event)


# -- colors ---------------------------------------------------------------------
def color_rgb(name):
    """Color name ('#rrggbb') as (r, g, b) floats in 0-1."""
    c = QtGui.QColor(name)
    return (round(c.redF(), 3), round(c.greenF(), 3), round(c.blueF(), 3))


def random_color():
    """Random bright color as '#rrggbb'."""
    return QtGui.QColor.fromHsvF(random.random(), 0.75, 0.95).name()


def style_color_button(b, name):
    """Show the color as the button background, with readable text."""
    b.setText(name)
    b.setStyleSheet('background-color: %s; color: %s;' % (
        name, 'black' if QtGui.QColor(name).lightness() > 127 else 'white'))


class ColorPicker(QtWidgets.QWidget):
    """'Color' check box plus a swatch button. color() is None (use the default) when unchecked."""

    def __init__(self, name='#ffffff', parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self._name = QtGui.QColor(name).name()
        self.check = QtWidgets.QCheckBox('Color')
        self.check.setToolTip('Use this color instead of the default CGO color')
        self.swatch = button(self._name, self.pick)
        style_color_button(self.swatch, self._name)
        lay = hbox(self.check, self.swatch)
        self.setLayout(lay)

    def pick(self):
        color = QtWidgets.QColorDialog.getColor(QtGui.QColor(self._name), self, 'Choose CGO color')
        if color.isValid():  # not cancelled
            self.setColor(color.name())

    def color(self):
        return self._name if self.check.isChecked() else None

    def setColor(self, name):
        """Set and check a color; None unchecks (default color)."""
        if name is None:
            self.check.setChecked(False)
            return
        self._name = QtGui.QColor(name).name()
        style_color_button(self.swatch, self._name)
        self.check.setChecked(True)


# -- layout helpers -------------------------------------------------------------
def button(text, slot, tooltip=None):
    b = QtWidgets.QPushButton(text)
    b.clicked.connect(lambda *a: slot())
    if tooltip:
        b.setToolTip(tooltip)
    return b


def add_button(text, slot):
    """The '>' button that adds a rule to the rules tree."""
    b = button('>', slot, 'Add to Packmol rules')
    b.setMaximumWidth(char_width(b, 1))
    b.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Expanding)
    return b


def hbox(*items, margins=0):
    """Horizontal layout; items are widgets, layouts, a str (label) or None (stretch)."""
    lay = QtWidgets.QHBoxLayout()
    lay.setContentsMargins(margins, margins, margins, margins)
    for it in items:
        if it is None:
            lay.addStretch(1)
        elif isinstance(it, str):
            lay.addWidget(QtWidgets.QLabel(it))
        elif isinstance(it, QtWidgets.QLayout):
            lay.addLayout(it)
        else:
            lay.addWidget(it)
    return lay


def group(title, layout_cls=QtWidgets.QVBoxLayout):
    g = QtWidgets.QGroupBox(title)
    lay = layout_cls(g)
    return g, lay


def info(parent, title, msg):
    QtWidgets.QMessageBox.information(parent, title, msg)


def error(parent, title, msg):
    QtWidgets.QMessageBox.critical(parent, title, msg)


def ask_yes_no(parent, title, msg):
    ans = QtWidgets.QMessageBox.question(parent, title, msg,
                                         QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
    return ans == QtWidgets.QMessageBox.Yes
