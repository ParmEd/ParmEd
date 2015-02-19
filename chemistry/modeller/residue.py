"""
This contains the basic residue template and residue building libraries
typically used in modelling applications
"""

import copy
from chemistry.topologyobjects import Atom, Bond, AtomList, TrackedList
try:
    import numpy as np
except ImportError:
    np = None
import warnings

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ResidueTemplate(object):
    """
    This is a residue template, which contains a listing of the atoms in the
    residue template as well as a mapping of which atoms are bonded to other
    atoms.

    Parameters
    ----------
    name : str, optional
        If provided, this is the name of the residue

    Attributes
    ----------
    atoms : :class:`AtomList`
        List of atoms in this residue
    bonds : :class:`TrackedList`
        List of the bonds between the atoms in this residue
    coordinates : numpy array or list
        The partial atomic coordinates in a linear iterable. If numpy is
        available, the coordinates are a numpy ndarray. Otherwise it is a list.
    connections : list of :class:`Atom`
        A list of all atoms that should form connections with atoms of another
        residue *besides* the head and tail atoms
    head : :class:`Atom` or None
        The atom that is connected to the residue that comes before this one
    tail : :class:`Atom` or None
        The atom that is connected to the *next* residue after this one
    """

    def __init__(self, name=''):
        self.atoms = AtomList()
        self.bonds = TrackedList()
        self.name = name
        self.head = None
        self.tail = None
        self.connections = []
        self._atomnames = set()

    def add_atom(self, atom):
        """ Adds an atom to this residue template

        Parameters
        ----------
        atom : :class:`Atom`
            The atom to add to this residue

        Raises
        ------
        ValueError if ``atom`` has the same name as another atom in this
        residue already
        """
        if atom.name in self._atomnames:
            raise ValueError('Residue already has atom named %s' % atom.name)
        atom.residue = self
        self.atoms.append(atom)
        self._atomnames.add(atom.name)

    def add_bond(self, atom1, atom2):
        """ Adds a bond between the two provided atoms in the residue

        Parameters
        ----------
        atom1 : :class:`Atom` or int
            One of the atoms in the bond. It must be in the ``atoms`` list of
            this ResidueTemplate. It can also be the atom index (index from 0)
            of the atom in the bond.
        atom2 : :class:`Atom` or int
            The other atom in the bond. It must be in the ``atoms`` list of this
            ResidueTemplate. It can also be the atom index (index from 0) of the
            atom in the bond.

        Raises
        ------
        IndexError if atom1 or atom2 are integers that are out of range of the
        number of atoms already in this template

        RuntimeError if atom1 or atom2 are :class:`Atom` instances but they are
        *not* in the atoms list of this ResidueTemplate

        Notes
        -----
        If atom1 and atom2 are already bonded, this routine does nothing
        """
        if isinstance(atom1, int):
            atom1 = self.atoms[atom1]
        if isinstance(atom2, int):
            atom2 = self.atoms[atom2]
        if atom1.list is not self.atoms or atom2.list is not self.atoms:
            raise RuntimeError('Both atoms must belong to template.atoms')
        # Do not add the same bond twice
        if atom1 not in atom2.bond_partners:
            self.bonds.append(Bond(atom1, atom2))

    @classmethod
    def from_residue(cls, residue):
        """
        This constructor creates a ResidueTemplate from a particular Residue
        object

        Parameters
        ----------
        residue : :class:`Residue`
            The residue from which to create a template
        """
        inst = cls(name=residue.name)
        for atom in residue:
            inst.add_atom(copy.copy(atom))
        for atom in residue:
            for bond in atom.bonds:
                try:
                    i1 = residue.atoms.index(bond.atom1)
                    i2 = residue.atoms.index(bond.atom2)
                except ValueError:
                    if bond.atom1 in residue:
                        oatom = bond.atom2
                        idx = residue.atoms.index(bond.atom1)
                    else:
                        oatom = bond.atom1
                        idx = residue.atoms.index(bond.atom2)
                    if oatom.residue.idx == residue.idx - 1:
                        inst.head = inst.atoms[idx]
                    elif oatom.residue.idx == residue.idx + 1:
                        inst.tail = inst.atoms[idx]
                    elif oatom.residue.idx == residue.idx:
                        # Don't know WHAT to do with it
                        warnings.warn('Cannot determine head/tail for '
                                      'unordered residues.')
                    else:
                        # Disulfide or something... not head or tail
                        inst.connections.append(inst.atoms[idx])
                else:
                    inst.add_bond(i1, i2)
        return inst

    @property
    def coordinates(self):
        """ Atomic coordinates, in Angstroms, of all atoms in the template """
        try:
            return self._crd
        except AttributeError:
            if np is None:
                raise ImportError('numpy is required to get ResidueTemplate '
                                  'coordinates')
            self._crd = np.array([[a.xx, a.xy, a.xz] for a in self]).flatten()
        return self._crd

    # Make ResidueTemplate look like a container of atoms
    def __len__(self):
        return len(self.atoms)
    def __getitem__(self, idx):
        return self.atoms[idx]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ResidueTemplateContainer(list):
    """
    A container of ResidueTemplate objects representing a unit with multiple
    residues
    """
    def __init__(self):
        self.box = None
