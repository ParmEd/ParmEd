"""
A module for efficiently building a pairlist with memory and CPU cost O(N)
"""
from collections import defaultdict
import numpy as np
from parmed.utils.six.moves import zip

def find_atom_pairs(struct, dist):
    """ Finds all pairs of atoms in a structure within a requested distance

    Parameters
    ----------
    struct : :class:`parmed.structure.Structure`
        The structure for which to find all atom pairs
    dist : float
        The maximum distance between identified pairs. No pairs separated by
        less than dist will be omitted, although some pairs greater than dist
        will be identified.

    Returns
    -------
    pairs : list of set of Atom
        The list of set of atoms that may be pairs

    Raises
    ------
    ValueError if ``struct`` does not have any coordinates.

    Notes
    -----
    This function does not do anything with periodic boundary conditions.
    """
    coords = struct.coordinates
    if coords is None:
        raise ValueError('Cannot find pairlist without coordinates')
    # Find voxel assignments for each atom -- need to move enclosing box so
    # origin is at lower left
    voxels = np.array((coords - coords.min(axis=0)) / dist, dtype=int)
    atom_voxel_map = defaultdict(set)
    for a, v in zip(struct.atoms, voxels):
        atom_voxel_map[tuple(v)].add(a)
    pairs = [set() for a in struct.atoms]
    for i, v in enumerate(voxels):
        for ii in range(-1, 2):
            for jj in range(-1, 2):
                for kk in range(-1, 2):
                    tup = tuple(v + [ii, jj, kk])
                    if tup in atom_voxel_map:
                        pairs[i] |= atom_voxel_map[tup]
    # Remove each atom from its own pairlist
    for i, a in enumerate(struct.atoms):
        pairs[i] -= {a}
    return pairs
