"""
Molecular volume Vmol = Vvdw + Vvoid (A^3) with the MoloVol algorithm, in pure Python (numpy only).

Merged from the pure-Python port of MoloVol (standalone_cli/pure_python/molovol_vmol.py and
molovol_pymol.py). MoloVol: https://github.com/jmaglic/MoloVol
Copyright (c) 2020-2025 Jasmin B. Maglic, Roy Lavendomme (MoloVol).
Distributed under the MIT License; the full text is in LICENSE-MoloVol.txt,
which must be kept with this file.

Coordinates and van der Waals radii come from PyMOL (each atom's `vdw` property); no elements file
is used. Algorithm (equivalent to the bottom level of the MoloVol octree):
  1. Cubic grid of voxels with edge GRID around all atoms.
  2. A voxel whose centre lies inside an atom (d < r) is an atom voxel. A voxel whose centre is
     at least r + probe away from every atom is a probe core voxel (the probe centre fits there).
  3. A remaining voxel is probe shell if a probe core voxel lies within the probe radius
     (MoloVol search limit: squared index distance <= floor((probe/grid + sqrt(2)/4)^2)),
     otherwise it is excluded void.
  4. Vmol = (atom voxels + void voxels) * grid^3.
  In two-probe mode the large probe runs first and only the voxels it cannot reach are probed again
  with the small probe.

Speedups (merged from molovol_vmol_pca.py and molovol_vmol_blocks.py), both used by
calcVolumeMoloVol:
  - The molecule is rotated to its principal axes when that gives fewer voxels (orient). Vmol is
    rotation invariant; only the grid discretisation changes (< 0.1 %).
  - Block skipping: only tiles of DEFAULT_TILE voxels reached by the r + probe box of an atom are
    processed, each with a halo of isqrt(search limit) voxels (two-probe: sum of both). Voxel
    counts are identical to the dense molecular_volume on the same grid.
"""
import math
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from pymol import cmd

# bit flags of the voxel grid
ATOM = 1      # centre inside an atom
NEAR_S = 2    # centre closer than r + small probe radius to an atom
NEAR_L = 4    # centre closer than r + large probe radius to an atom

# large numpy operations release the GIL, so the distance transform runs on all cores
N_THREADS = os.cpu_count() or 1


class InputError(Exception):
    pass


def make_grid(coords, radii, pad, grid, depth):
    """Origin and number of voxels per axis, as in the MoloVol Space constructor."""
    add = pad + radii.max()
    cart_min = coords.min(axis=0) - add
    cart_max = coords.max(axis=0) + add
    for k in range(3):
        mod = cart_min[k] - math.trunc(cart_min[k] / grid) * grid
        if mod != 0:
            cart_min[k] -= mod
    block = 2 ** depth
    shape = tuple(int(math.ceil(math.ceil((cart_max[k] - cart_min[k]) / grid) / block)) * block
                  for k in range(3))
    return cart_min, shape


def _index_box(coords, reach, origin, grid, shape):
    """Per-atom [lo, hi) voxel index ranges covering centre +- reach."""
    lo = np.floor((coords - reach[:, None] - origin) / grid - 0.5).astype(np.int64)
    hi = np.ceil((coords + reach[:, None] - origin) / grid - 0.5).astype(np.int64) + 1
    return np.maximum(lo, 0), np.minimum(hi, np.asarray(shape))


def mark_atoms(coords, radii, origin, shape, grid, probe_s, probe_l):
    """Grid of bit flags ATOM / NEAR_S / NEAR_L.

    Squared distances are compared with squared thresholds (no sqrt). In two-probe mode the atom and
    small probe tests only run on the inner part of each atom's large probe box.
    """
    flags = np.zeros(shape, dtype=np.uint8)
    centres = [origin[k] + grid * (np.arange(shape[k]) + 0.5) for k in range(3)]
    reach_s = radii + probe_s
    lo_s, hi_s = _index_box(coords, reach_s, origin, grid, shape)
    t_atom = radii ** 2
    t_s = reach_s ** 2
    if probe_l:
        reach_l = radii + probe_l
        lo, hi = _index_box(coords, reach_l, origin, grid, shape)
        t_l = reach_l ** 2
    else:
        reach_l, lo, hi = reach_s, lo_s, hi_s
    u_atom, u_s, u_l = np.uint8(ATOM), np.uint8(NEAR_S), np.uint8(NEAR_L)

    for i in range(len(radii)):
        if reach_l[i] <= 0:
            continue
        c = coords[i]
        (x0, y0, z0), (x1, y1, z1) = lo[i], hi[i]
        d2 = (((centres[0][x0:x1] - c[0]) ** 2)[:, None, None]
              + ((centres[1][y0:y1] - c[1]) ** 2)[None, :, None]
              + ((centres[2][z0:z1] - c[2]) ** 2)[None, None, :])
        if probe_l:
            local = (d2 < t_l[i]).view(np.uint8) * u_l
            # inner box of the small probe, relative to the large probe box
            inner = tuple(slice(lo_s[i, k] - lo[i, k], hi_s[i, k] - lo[i, k]) for k in range(3))
            d2_in, local_in = d2[inner], local[inner]
        else:
            local = np.zeros(d2.shape, dtype=np.uint8)
            d2_in, local_in = d2, local
        local_in |= (d2_in < t_s[i]).view(np.uint8) * u_s
        if radii[i] > 0:
            local_in |= (d2_in < t_atom[i]).view(np.uint8) * u_atom
        flags[x0:x1, y0:y1, z0:z1] |= local
    return flags


def search_limit(probe, grid):
    """Squared index distance within which a probe core voxel makes a voxel probe shell."""
    return int((probe / grid + math.sqrt(2) / 4) ** 2)


def _min_plus_pass(d, axis, h, cap):
    """In place: d[i] = min(cap, min over |s| <= h of d[i+s] + s^2) along one axis."""
    src = d.copy()
    tmp = np.empty_like(d)
    n = d.shape[axis]
    for s in range(1, min(h, n - 1) + 1):
        w = d.dtype.type(s * s)
        fwd = [slice(None)] * 3
        bwd = [slice(None)] * 3
        fwd[axis] = slice(s, n)
        bwd[axis] = slice(0, n - s)
        fwd, bwd = tuple(fwd), tuple(bwd)
        t = tmp[fwd]
        np.add(src[bwd], w, out=t)
        np.minimum(d[fwd], t, out=d[fwd])
        np.add(src[fwd], w, out=t)
        np.minimum(d[bwd], t, out=d[bwd])
    np.minimum(d, d.dtype.type(cap), out=d)


def _isqrt(n):
    """math.isqrt for Python < 3.8 (PyMOL bundles 3.7)."""
    h = int(math.sqrt(n))
    while h * h > n:
        h -= 1
    while (h + 1) * (h + 1) <= n:
        h += 1
    return h


def near_core(core, limit):
    """True where a core voxel lies within squared index distance <= limit (includes core voxels).

    Truncated separable squared distance transform: d^2 = dx^2 + dy^2 + dz^2 is minimised one axis at
    a time, and offsets beyond sqrt(limit) cannot contribute, so each pass only needs 2*isqrt(limit)
    shifted minimum operations on small integers. Exact for the <= limit test.
    """
    h = _isqrt(limit)
    cap = limit + 1
    dtype = np.uint8 if 2 * cap <= 255 else np.uint16
    d = (~core).view(np.uint8).astype(dtype)
    d *= dtype(cap)
    if h == 0:
        return d <= limit

    nx, ny, _ = d.shape
    n_chunks = 4 * N_THREADS
    x_edges = np.linspace(0, nx, min(nx, n_chunks) + 1).astype(int)
    y_edges = np.linspace(0, ny, min(ny, n_chunks) + 1).astype(int)

    def yz_passes(k):
        slab = d[x_edges[k]:x_edges[k + 1]]
        _min_plus_pass(slab, 1, h, cap)
        _min_plus_pass(slab, 2, h, cap)

    def x_pass(k):
        _min_plus_pass(d[:, y_edges[k]:y_edges[k + 1]], 0, h, cap)

    # y and z passes never cross x slabs, the x pass never crosses y chunks: no halos needed
    with ThreadPoolExecutor(N_THREADS) as ex:
        list(ex.map(yz_passes, range(len(x_edges) - 1)))
        list(ex.map(x_pass, range(len(y_edges) - 1)))
    return d <= limit


def molecular_volume(coords, radii, probe, probe2, grid, depth):
    two_probe = probe2 is not None and probe2 > 0
    pad = (probe2 if two_probe else probe) + 2 * grid
    origin, shape = make_grid(coords, radii, pad, grid, depth)
    flags = mark_atoms(coords, radii, origin, shape, grid, probe, probe2 if two_probe else None)

    n_atom = int(np.count_nonzero(flags & ATOM))
    if two_probe:
        # large probe first: everything it reaches (core or shell) is outside the molecular volume
        reached = near_core((flags & (ATOM | NEAR_L)) == 0, search_limit(probe2, grid))
        rest = ~reached
        del reached
        rest &= (flags & ATOM) == 0
        # small probe only inside the region the large probe could not reach
        core_s = rest & ((flags & NEAR_S) == 0)
        void = rest
        void &= (flags & NEAR_S) != 0
        void &= ~near_core(core_s, search_limit(probe, grid))
    else:
        reached = near_core((flags & (ATOM | NEAR_S)) == 0, search_limit(probe, grid))
        void = (flags & (ATOM | NEAR_S)) == NEAR_S
        void &= ~reached
    return (n_atom + int(np.count_nonzero(void))) * grid ** 3


# ---- principal-axes rotation (from molovol_vmol_pca.py) ----

def principal_axes_rotation(coords):
    """Coordinates centred and rotated so the principal axes become x, y, z (proper rotation)."""
    centred = coords - coords.mean(axis=0)
    _, vecs = np.linalg.eigh(np.cov(centred.T))  # eigenvalues ascending, vectors in columns
    if np.linalg.det(vecs) < 0:
        vecs[:, 0] = -vecs[:, 0]
    return centred @ vecs


def _n_voxels(coords, radii, pad, grid, depth):
    _, shape = make_grid(coords, radii, pad, grid, depth)
    return int(np.prod(shape, dtype=np.int64)), shape


def orient(coords, radii, pad, grid, depth):
    """Original or principal-axes coordinates, whichever gives fewer voxels.

    Returns (coords, rotated, shape_original, shape_used).
    """
    n0, shape0 = _n_voxels(coords, radii, pad, grid, depth)
    if len(coords) < 2:
        return coords, False, shape0, shape0
    rotated = principal_axes_rotation(coords)
    n1, shape1 = _n_voxels(rotated, radii, pad, grid, depth)
    if n1 < n0:
        return rotated, True, shape0, shape1
    return coords, False, shape0, shape0


# ---- block skipping (from molovol_vmol_blocks.py) ----

DEFAULT_TILE = 64


def near_core_local(core, limit):
    """As near_core, single-threaded (tiles are parallelised instead)."""
    h = _isqrt(limit)
    cap = limit + 1
    dtype = np.uint8 if 2 * cap <= 255 else np.uint16
    d = (~core).view(np.uint8).astype(dtype)
    d *= dtype(cap)
    if h:
        for axis in (1, 2, 0):
            _min_plus_pass(d, axis, h, cap)
    return d <= limit


def active_tiles(coords, radii, probe, origin, shape, grid, tile):
    """Indices (n, 3) of tiles intersected by the r + probe index box of at least one atom."""
    reach = radii + probe
    lo, hi = _index_box(coords, reach, origin, grid, shape)
    ok = (reach > 0) & np.all(hi > lo, axis=1)
    lo, hi = lo[ok] // tile, (hi[ok] - 1) // tile
    n_tiles = tuple(int(math.ceil(s / tile)) for s in shape)
    active = np.zeros(n_tiles, dtype=bool)
    for (x0, y0, z0), (x1, y1, z1) in zip(lo, hi):
        active[x0:x1 + 1, y0:y1 + 1, z0:z1 + 1] = True
    return np.argwhere(active)


def molecular_volume_blocks(coords, radii, probe, probe2, grid, depth, tile=DEFAULT_TILE, info=None):
    two_probe = probe2 is not None and probe2 > 0
    pad = (probe2 if two_probe else probe) + 2 * grid
    origin, shape = make_grid(coords, radii, pad, grid, depth)
    flags = mark_atoms(coords, radii, origin, shape, grid, probe, probe2 if two_probe else None)

    lim_s = search_limit(probe, grid)
    halo = _isqrt(lim_s)
    if two_probe:
        lim_l = search_limit(probe2, grid)
        halo += _isqrt(lim_l)
    tiles = active_tiles(coords, radii, probe, origin, shape, grid, tile)
    shape_a = np.asarray(shape)

    def work(t):
        a = t * tile
        b = np.minimum(a + tile, shape_a)
        ra = np.maximum(a - halo, 0)
        rb = np.minimum(b + halo, shape_a)
        f = flags[ra[0]:rb[0], ra[1]:rb[1], ra[2]:rb[2]]
        inner = tuple(slice(a[k] - ra[k], b[k] - ra[k]) for k in range(3))
        if two_probe:
            # large probe first; region edges are wrong by at most h_large voxels, which the small
            # probe (reach h_small) never sees from inside the tile
            rest = ~near_core_local((f & (ATOM | NEAR_L)) == 0, lim_l)
            rest &= (f & ATOM) == 0
            core_s = rest & ((f & NEAR_S) == 0)
            void = rest
            void &= (f & NEAR_S) != 0
            void &= ~near_core_local(core_s, lim_s)
        else:
            void = (f & (ATOM | NEAR_S)) == NEAR_S
            void &= ~near_core_local((f & (ATOM | NEAR_S)) == 0, lim_s)
        return (int(np.count_nonzero(f[inner] & ATOM)) + int(np.count_nonzero(void[inner])),
                int(np.prod(b - a)))

    with ThreadPoolExecutor(N_THREADS) as ex:
        results = list(ex.map(work, tiles))
    n_in = sum(r[0] for r in results)
    if info is not None:
        info.update(shape=shape, voxels=int(np.prod(shape, dtype=np.int64)),
                    tiles_active=len(tiles),
                    tiles_total=int(np.prod([math.ceil(s / tile) for s in shape])),
                    voxels_active=sum(r[1] for r in results))
    return n_in * grid ** 3


def get_atoms(mol='(all)', state=-1):
    """Coordinates (n, 3) and vdW radii (n,) of a PyMOL selection, in one consistent pass."""
    rows = []
    cmd.iterate_state(int(state), '(%s)' % mol, 'rows.append((x, y, z, vdw))', space={'rows': rows})
    if not rows:
        raise InputError('Empty selection: %s' % mol)
    data = np.array(rows, dtype=np.float64)
    return data[:, :3], data[:, 3]


def calcVolumeMoloVol(mol='(all)', probe=1.2, probe2=None, grid=0.2, depth=4, state=-1,
                      tile=DEFAULT_TILE, rotate=True, blocks=True):
    """
    Returns molecular volume Vmol (A^3) of a PyMOL selection with the MoloVol algorithm.

    Parameters
    ----------
    mol    : string. PyMOL object name or selection
    probe  : float. Probe radius (0 gives the van der Waals volume only)
    probe2 : float or None. Large probe radius; > 0 enables two-probe mode
    grid   : float. Grid resolution (voxel edge), smaller is more accurate, slower and uses more memory
    depth  : int. Octree depth, only sets the grid padding
    state  : int. Object state, -1 for current
    tile   : int. Tile edge in voxels for block skipping
    rotate : bool. Grid in the principal-axes frame when that gives fewer voxels
    blocks : bool. Skip empty tiles (same result as the dense grid, faster)
    Raises InputError on invalid input or an empty selection.
    """
    if probe2 is not None and probe2 <= 0:
        probe2 = None
    if grid <= 0 or probe < 0 or not 0 <= depth <= 12:
        raise InputError('Invalid numeric input (grid must be > 0, probe >= 0).')
    if probe2 is not None and probe2 < probe:
        raise InputError('The large probe must have a larger radius than the small probe '
                         '(Solvent Radius).')
    coords, radii = get_atoms(mol, state)
    radii = np.maximum(radii, 0.0)
    if not (radii > 0).any():
        raise InputError('No atoms with non-zero radius in %s.' % mol)
    probe, grid, depth = float(probe), float(grid), int(depth)
    if rotate:
        pad = (probe2 if probe2 else probe) + 2 * grid  # same padding as molecular_volume
        coords, _, _, _ = orient(coords, radii, pad, grid, depth)
    if blocks:
        return molecular_volume_blocks(coords, radii, probe, probe2, grid, depth, int(tile))
    return molecular_volume(coords, radii, probe, probe2, grid, depth)
