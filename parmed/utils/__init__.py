""" Various utilities used by ParmEd that don't really fit elsewhere """
from parmed.exceptions import MoleculeError as _MoleculeError

__all__ = ['six', 'io', 'timer', 'which', 'tag_molecules']

def which(prog):
    """ Returns the full path of a program if it exists in PATH

    Parameters
    ----------
    prog : str
        The name of a program to try and locate in PATH

    Returns
    -------
    path : str or None
        The full path of the program. If it cannot be found, None
    """
    import os
    def is_exe(fpath):
        if os.path.isdir(fpath): return False
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    fpath, fprog = os.path.split(prog)
    if fpath:
        if is_exe(prog):
            return prog
        return None
    for fpath in os.environ['PATH'].split(os.pathsep):
        trial = os.path.join(fpath, prog)
        if is_exe(trial):
            return trial
    return None

def tag_molecules(struct):
    """
    Sets the ``marked`` attribute of every Atom in struct to the molecule number
    it is a part of. If no bonds are present, every atom is its own molecule.

    Parameters
    ----------
    struct : :class:`parmed.Structure`
        Input structure to tag the molecules for
    """
    # Make sure our recursion limit is large enough, but never shrink it
    from sys import setrecursionlimit, getrecursionlimit
    setrecursionlimit(max(len(struct.atoms), getrecursionlimit()))

    if not struct.bonds:
        for i, atom in enumerate(struct.atoms):
            atom.marked = i + 1
        return
    # We do have bonds, this is the interesting part
    struct.atoms.unmark()
    mol_id = 1
    for atom in struct.atoms:
        if atom.marked: continue
        atom.marked = mol_id
        _set_owner(atom, mol_id)
        mol_id += 1

def _set_owner(atm, mol_id):
    """ Recursively sets ownership of given atom and all bonded partners """
    for partner in atm.bond_partners:
        if not partner.marked:
            partner.marked = mol_id
            _set_owner(partner, mol_id)
        elif partner.marked != mol_id:
            raise _MoleculeError('Atom %d in multiple molecules' %
                                 partner.idx)
