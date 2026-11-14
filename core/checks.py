"""Volume checks of Packmol constraints against molecular volumes (no GUI code)."""
import os

import numpy as np
from pymol import cmd

from .rules import GLOBAL_DIRECTIONS, head, valid_region
from ..volume import (volCGO, interVolPrisms, interVolSpheres, interVolSP,
                      totalInsidePrism, totalInsideSphere)


def check_volumes(nodes, sidemax, len_tolerance, write, mol_volume):
    """
    Check every structure block of a rule tree. write(text) receives the report.
    sidemax is Packmol's maximum box side (system limits).
    mol_volume(object name) returns the molecular volume (settings.mol_volume: method chosen in
    Plugin Settings); its exceptions are passed to the caller.
    """
    lenStdBox = sidemax
    volStdBox = lenStdBox ** 3
    for text, children in nodes:
        subunits = text.split()
        if not subunits:  # blank line outside structure/atom kw
            continue
        if subunits[0] != 'structure':
            continue
        nested = []  # new nested constraint list for new structure
        if len(subunits) <= 1:
            write("Structure keyword do not have any argument \nCan not check further!\n\n")
            break
        for _ in _global_keywords(children):
            write("NOTE: pbc found inside structure %s. pbc is a global keyword; it works anyway, "
                  "but consider moving it out of the structure block.\n" % subunits[1])
        heads = [head(t) for t, _ in children]
        if 'chain' in heads and 'changechains' in heads:
            write("NOTE: structure %s has both chain and changechains. Packmol does not accept both "
                  "in one structure and stops with an error.\n" % subunits[1])
        if not os.path.isfile(subunits[1]):
            write("Structure %s not found in actual path!\n\tactual path: %s\nCan not check further!\n\n"
                  % (subunits[1], os.path.abspath(os.getcwd())))
            break
        currentView = cmd.get_view()  # save camera position, loading object changes
        tmp = cmd.get_unused_name('_packmol_tmp')
        cmd.load(subunits[1], tmp)
        try:
            volMol = mol_volume(tmp)
            bBox = np.array(cmd.get_extent(tmp))  # [[x0,y0,z0],[x1,y1,z1]]
        finally:
            cmd.delete(tmp)
            cmd.set_view(currentView)
        dx, dy, dz = np.abs(bBox[0] - bBox[1])
        lenMol = np.sqrt(dx * dx + dy * dy + dz * dz) + len_tolerance
        for text2, _ in children:
            subunits2 = text2.split()
            if len(subunits2) <= 1:
                continue
            if subunits2[0] == 'number':
                volMol *= int(subunits2[1])  # compute for all molecules in structure
            if subunits2[0] in ['inside', 'outside']:
                if not valid_region(subunits2[1], subunits2[2:]):
                    write("Wrong parameters in rule: %s\nSkipped.\n\n" % text2)
                    continue
                volConstraint, lenConstraint = volCGO(subunits2[1:])  # (volume, max length)
                nested.append(subunits2)
                if subunits2[0] == 'outside':
                    volConstraint = volStdBox - volConstraint
                    lenConstraint = lenStdBox - lenConstraint
                line = "volMol= %.3f \t volConst=%.3f \n" % (volMol, volConstraint)
                line += "lenMol= %.3f \t lenConst=%.3f \n" % (lenMol, lenConstraint)
                write(line)
                if volMol > volConstraint:
                    write("WARNING volume of constraint not enough to pack molecular volume.\n")
                elif lenMol > lenConstraint:  # if volume is OK, then check length
                    write("WARNING length of constraint not enough to pack molecular length.\n")
                else:
                    write("\n OK.\n")
        if len(nested) > 1:
            write("WARNING nested constraints. Consider revising.\n")
            if len(nested) == 2:
                check_nested_volume(volMol, nested, sidemax, write)
            else:
                write("Still not available volume check for constraints nested level larger than 2.\n")
        else:
            write("\n End.\n")


def _global_keywords(children):
    """Global keyword lines (pbc) inside a structure block, including its atoms blocks."""
    out = []
    for text, sub in children:
        if head(text) in GLOBAL_DIRECTIONS:
            out.append(text)
        elif head(text) == 'atoms':
            out.extend(_global_keywords(sub))
    return out


def _cube_to_box(c):
    return [c[0]] + ['box'] + c[2:5] + list(np.float64(c[2:5]) + np.float64(c[-1]))


def check_nested_volume(volMol, nested, sidemax, write):
    volConstraint = 0
    volInter = 0
    volStdBox = sidemax ** 3

    nested = [_cube_to_box(n) if n[1] == 'cube' else n for n in nested]

    # calculate INTERSECTION Volume
    if nested[0][1] == nested[1][1]:  # same geometric elements
        if nested[0][1] == 'box':
            volInter = interVolPrisms(nested[0][2:], nested[1][2:])
        elif nested[0][1] == 'sphere':
            volInter = interVolSpheres(nested[0][2:], nested[1][2:])
    else:
        kinds = (nested[0][1], nested[1][1])
        if kinds in (('sphere', 'box'), ('box', 'sphere')):
            if nested[0][1] == 'sphere':
                s, p = nested[0], nested[1]
            else:
                s, p = nested[1], nested[0]
            if totalInsidePrism(s[2:], p[2:]):
                volInter = volCGO(s[1:])[0]
            elif totalInsideSphere(s[2:], p[2:]):
                volInter = volCGO(p[1:])[0]
            else:
                volInter = interVolSP(.75, s[2:], p[2:])

    d0, d1 = nested[0][0], nested[1][0]
    if d0 == 'inside' and d1 == 'inside':
        volConstraint = volInter
    elif d0 == 'inside' and d1 == 'outside':
        volConstraint = volCGO(nested[0][1:])[0] - volInter
    elif d0 == 'outside' and d1 == 'inside':
        volConstraint = volCGO(nested[1][1:])[0] - volInter
    elif d0 == 'outside' and d1 == 'outside':
        v1 = volCGO(nested[0][1:])[0]
        v2 = volCGO(nested[1][1:])[0]
        volConstraint = volStdBox - (v1 + v2 - volInter)

    write("Intersecting volume:  %.3f.\n" % volInter)
    write("Final constraint volume:  %.3f.\n" % volConstraint)
    if volMol > volConstraint:
        write("WARNING volume of nested constraint %.3f not enough to pack molecular volume.\n" % volConstraint)
    else:
        write("\nOK.\n")
