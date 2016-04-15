""" Various utilities used by ParmEd that don't really fit elsewhere """
from parmed.exceptions import MoleculeError as _MoleculeError
from parmed.utils.pairlist import find_atom_pairs
import sys

__all__ = ['six', 'io', 'timer', 'which', 'tag_molecules', 'PYPY',
           'canonical_improper_order', 'find_atom_pairs']

PYPY = '__pypy__' in sys.builtin_module_names

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

def canonical_improper_order(atom1, atom2, atom3, atom4, center=1):
    """
    Controls how improper torsion keys are generated from Structures.
    Different programs have different conventions as far as where the
    "central" atom is placed.

    Note, different programs use different conventions for where the "central"
    atom comes in the overall torsion ordering. CHARMM puts the central atom
    first, whereas AMBER puts the central atom third. A central atom is defined
    as the one that is bonded to the other 3.

    Parameters
    ----------
    atom1 : :class:`parmed.topologyobjects.Atom`
        The first atom in the improper
    atom2 : :class:`parmed.topologyobjects.Atom`
        The second atom in the improper
    atom3 : :class:`parmed.topologyobjects.Atom`
        The third atom in the improper
    atom4 : :class:`parmed.topologyobjects.Atom`
        The fourth atom in the improper
    center : int, optional
        Which location represents the *center* atom. Default is 1 (first
        location).

    Returns
    -------
    a1, a2, a3, a4 : tuple of :class:`parmed.topologyobjects.Atom`
        The atoms in the necessary parameter order
    """
    if center not in (1, 2, 3, 4):
        raise ValueError('center must be 1, 2, 3, or 4')
    all_atoms = set([atom1, atom2, atom3, atom4])
    for atom in all_atoms:
        for atom2 in all_atoms:
            if atom2 is atom: continue
            if not atom2 in atom.bond_partners:
                break
        else:
            # We found our central atom
            others = sorted(all_atoms - set([atom]))
            central = atom
            break
    else:
        # TODO use "center" keyword to determine which atom is center
        # No atom identified as "central". Just assume that the third is
        central = atom3
        others = sorted([atom1, atom2, atom4])
    if center == 1:
        return central, others[0], others[1], others[2]
    elif center == 2:
        return others[0], central, others[1], others[2]
    elif center == 3:
        return others[0], others[1], central, others[2]
    elif center == 4:
        return others[0], others[1], others[2], central
    assert False, 'Should not be here'
