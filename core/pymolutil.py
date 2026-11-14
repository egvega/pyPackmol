"""Small PyMOL helpers shared by the tabs (no GUI code)."""
from threading import Thread

import numpy as np

from pymol import cmd, cgo
from chempy import cpv

AXIS_NAME = 'axis3D'
LABEL_SUFFIX = '_lbl'  # label pseudoatoms of Set Up Input constraints: <cgo name>_lbl


class PutCenterCallback(object):
    prev_v = None

    def __init__(self, name, corner=0):
        self.name = name
        self.corner = corner
        self.cb_name = cmd.get_unused_name('_cb')

    def load(self):
        cmd.load_callback(self, self.cb_name)

    def __call__(self):
        if self.name not in cmd.get_names('objects'):
            Thread(None, cmd.delete, args=(self.cb_name,)).start()
            return

        v = cmd.get_view()
        if v == self.prev_v:
            return
        self.prev_v = v

        t = v[12:15]

        if self.corner:
            vp = cmd.get_viewport()
            R_mc = [v[0:3], v[3:6], v[6:9]]
            off_c = [0.15 * v[11] * vp[0] / vp[1], 0.15 * v[11], 0.0]
            if self.corner in [2, 3]:
                off_c[0] *= -1
            if self.corner in [3, 4]:
                off_c[1] *= -1
            off_m = cpv.transform(R_mc, off_c)
            t = cpv.add(t, off_m)

        z = -v[11] / 30.0
        m = [z, 0, 0, 0, 0, z, 0, 0, 0, 0, z, 0, t[0] / z, t[1] / z, t[2] / z, 1]
        cmd.set_object_ttt(self.name, m)


def show_axis(enabled, name=AXIS_NAME):
    """Puts coordinate axes to the lower left corner of the viewport (or removes them)."""
    if not enabled:
        cmd.delete(name)
        return
    cmd.set('auto_zoom', 0)

    w = 0.06  # cylinder width
    l = 0.75  # cylinder length
    h = 0.25  # cone hight
    d = w * 1.618  # cone base diameter

    obj = [cgo.CYLINDER, 0.0, 0.0, 0.0,   l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
           cgo.CYLINDER, 0.0, 0.0, 0.0, 0.0,   l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
           cgo.CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
           cgo.CONE,   l, 0.0, 0.0, h + l, 0.0, 0.0, d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0,
           cgo.CONE, 0.0,   l, 0.0, 0.0, h + l, 0.0, d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0,
           cgo.CONE, 0.0, 0.0,   l, 0.0, 0.0, h + l, d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    PutCenterCallback(name, 1).load()
    cmd.load_cgo(obj, name)


def unique_name(base, sep='_'):
    """base, base_1, base_2, ... first name not used by a PyMOL object."""
    name = base
    num = 1
    while name in cmd.get_names('objects'):
        name = base + sep + str(num)
        num += 1
    return name


def load_cgo_keep_view(obj, name):
    view = cmd.get_view()
    cmd.load_cgo(obj, name)
    cmd.set_view(view)


# ------------------------------------------------------ constraint labels --
def is_label(name):
    """True for the label pseudoatom objects of constraints (<cgo name>_lbl, _lbl1, ...)."""
    return LABEL_SUFFIX in name


def label_objects():
    return [n for n in cmd.get_names('objects') if is_label(n)]


def exclude_labels(sele):
    """Selection without the constraint label pseudoatoms."""
    names = label_objects()
    return '(%s) and not (%s)' % (sele, ' '.join(names)) if names else sele


def cgo_vertices(obj):
    """(n, 3) VERTEX coordinates of a CGO list made of COLOR, BEGIN, VERTEX and END."""
    out = []
    i, n = 0, len(obj)
    while i < n:
        op = obj[i]
        if op == cgo.COLOR:
            i += 4
        elif op == cgo.BEGIN:
            i += 2
        elif op == cgo.VERTEX:
            out.append(obj[i + 1:i + 4])
            i += 4
        elif op == cgo.END:
            i += 1
        else:  # other opcodes are not used by the constraint builders
            break
    out = np.array(out, dtype=np.float64).reshape(-1, 3)
    return out[np.isfinite(out).all(axis=1)]


def label_point(candidates, others, tol=5.0):
    """First candidate at least tol from every point in others; else the one farthest from them."""
    if len(others) == 0:
        return candidates[0]
    d2 = ((candidates[:, None, :] - np.asarray(others)[None, :, :]) ** 2).sum(axis=2).min(axis=1)
    ok = np.nonzero(d2 >= tol * tol)[0]
    return candidates[ok[0]] if len(ok) else candidates[np.argmax(d2)]


def add_cgo_label(cgo_name, obj, color, others, tol=5.0):
    """
    Hidden pseudoatom on a vertex of the CGO list obj, labelled cgo_name in color (r, g, b),
    at least tol from the points in others when possible. Returns its object name, or None.
    """
    points = cgo_vertices(obj)
    if not len(points):
        return None
    name = unique_name(cgo_name + LABEL_SUFFIX, sep='')
    cmd.pseudoatom(name, pos=[float(x) for x in label_point(points, others, tol)], label=cgo_name)
    cmd.hide('everything', name)
    cmd.show('labels', name)
    cmd.set_color(name + '_color', [float(c) for c in color])
    cmd.set('label_color', name + '_color', name)
    return name
