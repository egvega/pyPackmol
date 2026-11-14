"""
Packmol input rules: keyword tables and .inp read/write.
No GUI and no PyMOL code here. A rule tree is a list of nodes, each node
being [text, children] where children is again a list of nodes.
"""
import re

REGIONS = ['cube', 'box', 'sphere', 'ellipsoid', 'cylinder', 'plane', 'xygauss', 'fixed']
DIRECTIONS = ['inside', 'outside', 'over', 'above', 'below', 'pbc']
GLOBAL_DIRECTIONS = ['pbc']  # written at top level, without a region word

# number of parameters for each region keyword
REGION_NPARAMS = {'fixed': 6, 'cube': 4, 'box': 6, 'sphere': 4, 'ellipsoid': 7, 'plane': 4,
                  'xygauss': 6, 'cylinder': 8, 'None': 0}
PARAM_NAMES = {'fixed': ['x', 'y', 'z', 'a', 'b', 'g'], 'cube': ['x0', 'y0', 'z0', 'd'],
               'box': ['x0', 'y0', 'z0', 'x1', 'y1', 'z1'], 'sphere': ['a', 'b', 'c', 'd'],
               'ellipsoid': ['a1', 'b1', 'c1', 'a2', 'b2', 'c2', 'd'], 'plane': ['a', 'b', 'c', 'd'],
               'cylinder': ['a1', 'b1', 'c1', 'a2', 'b2', 'c2', 'd', 'l'],
               'xygauss': ['a1', 'b1', 'a2', 'b2', 'c', 'h']}

ADDITIONAL_SETTINGS = ['discale', 'maxit', 'add_box_sides', 'sidemax', 'avoid_overlap', 'nloop', 'precision',
                       'writeout', 'iprint1', 'iprint2', 'fbins', 'fscale', 'movebadrandom', 'add_amber_ter',
                       'randominitialpoint', 'writebad', 'chkgrad', 'short_radius', 'short_radius_scale',
                       'movefrac', 'disable_movebad', 'packall', 'amber_ter_preserve', 'writecrd',
                       'hexadecimal_indices', 'ignore_conect', 'non_standard_conect', 'nloop0',
                       'maxmove', 'chain', 'segid', 'connect']
INT_SETTINGS = ['maxit', 'nloop', 'writeout', 'iprint1', 'iprint2', 'nloop0', 'maxmove']
FLAG_SETTINGS = ['movebadrandom', 'add_amber_ter', 'randominitialpoint', 'writebad', 'chkgrad',
                 'disable_movebad', 'packall', 'amber_ter_preserve', 'hexadecimal_indices',
                 'ignore_conect', 'non_standard_conect']
STRING_SETTINGS = ['writecrd', 'chain', 'segid', 'connect']  # written as typed
STRUCTURE_SETTINGS = ['nloop', 'fscale', 'maxmove', 'chain', 'segid', 'connect']  # go inside a structure block
ATOMS_SETTINGS = ['short_radius', 'short_radius_scale']  # go inside a structure or atoms block

BLOCKS = ('structure', 'atoms')

# own CGO color of a constraint, saved as a comment at the end of its line
COLOR_RE = re.compile(r'\s*#\s*color\s+(#[0-9a-fA-F]{6})\b')


def valid_region(region, params):
    """True if params (strings) are the right number of numeric values for the region keyword."""
    if len(params) != REGION_NPARAMS.get(region.lower()):
        return False
    try:
        [float(par) for par in params]
    except ValueError:
        return False
    return True


def head(text):
    """First keyword of a rule line, lower case ('' for blank lines)."""
    parts = text.split()
    return parts[0].lower() if parts else ''


def pop_color(text):
    """(text without its color comment, '#rrggbb' or None)."""
    m = COLOR_RE.search(text)
    if not m:
        return text, None
    rest = (text[:m.start()] + text[m.end():]).rstrip()
    if rest.endswith('#'):  # nothing else was in the comment
        rest = rest[:-1].rstrip()
    return rest, m.group(1).lower()


def color_comment(color):
    return '# color ' + color


def _strip_comment(lin):
    """Line without its comment (everything from the first '#')."""
    return lin.split('#', 1)[0].strip()


def is_constraint(text):
    """Constraint (region) or pbc line."""
    sub = text.split()
    return len(sub) > 1 and (sub[1].lower() in REGIONS or sub[0].lower() in GLOBAL_DIRECTIONS)


def parse_inp(lines):
    """
    Parse .inp lines into a rule tree. Comments (full-line and end-of-line) and empty lines are
    skipped. Only the '# color #rrggbb' comment of a constraint line is kept, for the GUI to read.
    """
    root = []
    stack = [root]
    for lin in lines:
        rest, color = pop_color(lin)
        lin = _strip_comment(rest)
        if not lin:
            continue
        if color and is_constraint(lin):
            lin += '  ' + color_comment(color)
        key = head(lin)
        if key in BLOCKS:
            node = [lin, []]
            stack[-1].append(node)
            stack.append(node[1])
        elif key == 'end':
            if len(stack) > 1:
                stack.pop()
        else:
            stack[-1].append([lin, []])
    return root


def format_inp(nodes, depth=0):
    """Rule tree to .inp lines (without newline), indented two spaces per level."""
    out = []
    for text, children in nodes:
        out.append('  ' * depth + text)
        key = head(text)
        if key in BLOCKS:
            out.extend(format_inp(children, depth + 1))
            out.append('  ' * depth + 'end ' + key)
    return out


def read_inp(fn):
    with open(fn, 'rt') as f:
        return parse_inp(f.readlines())


def write_inp(fn, nodes):
    with open(fn, 'wt') as f:
        for line in format_inp(nodes):
            f.write(line + '\n')


def structure_counts(nodes):
    """[(structure file, number), ...] for each structure block that has a 'number' keyword."""
    out = []
    for text, children in nodes:
        parts = text.split()
        if head(text) != 'structure' or len(parts) < 2:
            continue
        for ctext, _ in children:
            cparts = ctext.split()
            if head(ctext) == 'number' and len(cparts) > 1:
                try:
                    out.append((parts[1], int(cparts[1])))
                except ValueError:
                    pass
                break
    return out


def output_from_inp(fn):
    """Value of the 'output' keyword of an input file, or None if missing/wrong."""
    out = None
    with open(fn, 'rt') as f:
        for lin in f:
            sl = _strip_comment(lin).split()
            if len(sl) == 2 and sl[0].lower() == 'output':
                out = sl[1]
    return out


def seed_from_inp(fn):
    """Integer value of the last 'seed' keyword of an input file, or None if missing/wrong."""
    seed = None
    with open(fn, 'rt') as f:
        for lin in f:
            sl = _strip_comment(lin).split()
            if len(sl) == 2 and sl[0].lower() == 'seed':
                try:
                    seed = int(sl[1])
                except ValueError:
                    seed = None
    return seed
