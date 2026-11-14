"""Rebuild molecules from a Packmol restart file (no GUI or PyMOL code)."""
import numpy as np

# type, ID, name, alt, resn, chain, resi, alt, x, y, z, q, b, elem
PDB_FMT = '{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n'
# no builtins here: iterate_state with a custom space has none (ID and resv are already ints)
TEMPLATE_EXPR = '[type, ID, name, alt, resn, chain, resv, alt, x, y, z, q, b, elem]'


def euler_matrix(beta, gamma, theta):
    """Packmol rotation matrix from its three euler angles (radians)."""
    cb, sb = np.cos(beta), np.sin(beta)
    cg, sg = np.cos(gamma), np.sin(gamma)
    ct, st = np.cos(theta), np.sin(theta)
    trmat = np.zeros((3, 3))
    trmat[0, 0] = -sb * sg * ct + cb * cg
    trmat[0, 1] = -sb * cg * ct - cb * sg
    trmat[0, 2] = sb * st
    trmat[1, 0] = cb * sg * ct + sb * cg
    trmat[1, 1] = cb * cg * ct - sb * sg
    trmat[1, 2] = -cb * st
    trmat[2, 0] = sg * st
    trmat[2, 1] = cg * st
    trmat[2, 2] = ct
    return trmat


def place_molecule(template, transformation, resi):
    """
    Rotate and translate a copy of template atoms (TEMPLATE_EXPR rows).
    Returns (atoms, trmat, shift) where shift is the final translation that puts
    the barycenter on the restart position.
    """
    trmat = euler_matrix(*transformation[3:6])
    xyz = np.array([a[8:11] for a in template], dtype=float)
    new = transformation[:3] + xyz.dot(trmat)
    shift = transformation[:3] - new.mean(axis=0)
    new += shift
    atoms = []
    for a, (x, y, z) in zip(template, new):
        at = list(a)
        at[6] = resi
        at[8], at[9], at[10] = x, y, z
        atoms.append(at)
    return atoms, trmat, shift


def count_positions(lines):
    """Number of restart lines before the first empty one."""
    n = 0
    for line in lines:
        if not line.strip():
            break
        n += 1
    return n


def iter_restart(templates, lines):
    """
    templates is a list of (template_atoms, count): the first count restart lines
    place the first template, the next ones the second template, and so on.
    Yields (n, transformation, trmat, shift, pdb_text) per restart line, n counting
    from 0 within each template (residue numbers restart at 1 for each one).
    Stops at the first empty line, like Packmol does, or when templates run out.
    """
    lines = iter(lines)
    for template, count in templates:
        for n in range(count):
            line = next(lines, '').strip()
            if not line:
                return
            transformation = np.array(line.split(), dtype=float)
            atoms, trmat, shift = place_molecule(template, transformation, n + 1)
            pdb = ''.join(PDB_FMT.format(*at) for at in atoms)
            yield n, transformation, trmat, shift, pdb
