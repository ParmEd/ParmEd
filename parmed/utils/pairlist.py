"""
A module for efficiently building a pairlist with memory and CPU cost O(N)
"""

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
    if struct.coordinates is None:
        raise ValueError('Cannot find pairlist without coordinates')
