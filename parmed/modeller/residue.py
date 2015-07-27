"""
This contains the basic residue template and residue building libraries
typically used in modelling applications
"""

import copy
try:
    import pandas as pd
except ImportError:
    pd = None
from parmed.residue import AminoAcidResidue, RNAResidue, DNAResidue
from parmed.topologyobjects import Atom, Bond, AtomList, TrackedList
try:
    import numpy as np
except ImportError:
    np = None
import warnings

__all__ = ['PROTEIN', 'NUCLEIC', 'SOLVENT', 'UNKNOWN', 'ResidueTemplate',
           'ResidueTemplateContainer', 'PatchTemplate']

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _ResidueType(object):
    """ Singleton for various types of residues """
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return '<ResidueType %s>' % self.name

    def __str__(self):
        return self.name

PROTEIN = _ResidueType('PROTEIN')
NUCLEIC = _ResidueType('NUCLEIC')
SOLVENT = _ResidueType('SOLVENT')
UNKNOWN = _ResidueType('UNKNOWN')

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
    coordinates : np.ndarray(natom, 3)
        The partial atomic coordinates
    connections : list of :class:`Atom`
        A list of all atoms that should form connections with atoms of another
        residue *besides* the head and tail atoms
    head : :class:`Atom` or None
        The atom that is connected to the residue that comes before this one
    tail : :class:`Atom` or None
        The atom that is connected to the *next* residue after this one
    first_patch : :class:`ResidueTemplate` or None
        If it is not None, this is the patch whose tail is added to the head
        atom of this residue when this residue is the first in a chain
    last_patch : :class:`ResidueTemplate` or None
        If it is not None, this is the patch whose head is added to the tail
        atom of this residue when this residue is the last in a chain
    groups : list of list(:class:`Atom`)
        If set, each group is a list of Atom instances making up each group
    """

    def __init__(self, name=''):
        self.atoms = AtomList()
        self.bonds = TrackedList()
        self.name = name
        self.head = None
        self.tail = None
        self.connections = []
        self._atomnames = set()
        self.type = UNKNOWN
        self.first_patch = None
        self.last_patch = None
        self.groups = []
        self.cmaps = []

    def __repr__(self):
        if self.head is not None:
            head = self.head.name
        else:
            head = 'None'
        if self.tail is not None:
            tail = self.tail.name
        else:
            tail = 'None'
        return '<%s %s: %d atoms; %d bonds; head=%s; tail=%s>' % (
                    type(self).__name__, self.name, len(self.atoms),
                    len(self.bonds), head, tail)

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
        atom1 : :class:`Atom` or int or str
            One of the atoms in the bond. It must be in the ``atoms`` list of
            this ResidueTemplate. It can also be the atom index (index from 0)
            of the atom in the bond.
        atom2 : :class:`Atom` or int or str
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
        If atom1 and atom2 are already bonded, this routine does nothing. If
        atom1 or atom2 are strings, then they will match the first instance of
        the atom name that is the same as the atom name passed.
        """
        if not isinstance(atom1, Atom):
            atom1 = self[atom1]
        if not isinstance(atom2, Atom):
            atom2 = self[atom2]
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
            self._crd = np.array([[a.xx, a.xy, a.xz] for a in self])
        return self._crd

    # Make ResidueTemplate look like a container of atoms, also indexable by the
    # atom name
    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)

    def __contains__(self, atom):
        if isinstance(atom, Atom):
            return atom in self.atoms
        if isinstance(atom, str):
            for a in self.atoms:
                if a.name == atom:
                    return True
        return False

    def __getitem__(self, idx):
        if isinstance(idx, str):
            for atom in self.atoms:
                if atom.name == idx:
                    return atom
            raise IndexError('Atom %s not found in %s' % (idx, self.name))
        return self.atoms[idx]

    def to_dataframe(self):
        """ Create a pandas dataframe from the atom information

        Returns
        -------
        df : :class:`pandas.DataFrame`
            The pandas DataFrame with all of the atomic properties

        Notes
        -----
        The DataFrame will be over all atoms. The columns will be the attributes
        of the atom (as well as its containing residue). Some columns will
        *always* exist. Others will only exist if those attributes have been set
        on the Atom instances (see the :class:`Atom` docs for possible
        attributes and their meaning). The columns that will always be present
        are:

            - number : int
            - name : str
            - type : str
            - atomic_number : int
            - charge : float
            - mass : float
            - nb_idx : int
            - radii : float
            - screen : float
            - occupancy : float
            - bfactor : float
            - altloc : str
            - tree : str
            - join : int
            - irotat : int
            - rmin : float
            - epsilon : float
            - rmin_14 : float
            - epsilon_14 : float

        The following attributes are optionally present if they were present in
        the original file defining the structure:

            - xx : float (x-coordinate position)
            - xy : float (y-coordinate position)
            - xz : float (z-coordinate position)
            - vx : float (x-coordinate velocity)
            - vy : float (y-coordinate velocity)
            - vz : float (z-coordinate velocity)
        """
        if pd is None:
            raise ImportError('pandas is not available; cannot create a pandas '
                              'DataFrame from this Structure')
        ret = pd.DataFrame()

        ret['number'] = [atom.number for atom in self.atoms]
        ret['name'] = [atom.name for atom in self.atoms]
        ret['type'] = [atom.type for atom in self.atoms]
        ret['atomic_number'] = [atom.atomic_number for atom in self.atoms]
        ret['charge'] = [atom.charge for atom in self.atoms]
        ret['mass'] = [atom.mass for atom in self.atoms]
        ret['nb_idx'] = [atom.nb_idx for atom in self.atoms]
        ret['radii'] = [atom.radii for atom in self.atoms]
        ret['screen'] = [atom.screen for atom in self.atoms]
        ret['occupancy'] = [atom.occupancy for atom in self.atoms]
        ret['bfactor'] = [atom.bfactor for atom in self.atoms]
        ret['altloc'] = [atom.altloc for atom in self.atoms]
        ret['tree'] = [atom.tree for atom in self.atoms]
        ret['join'] = [atom.join for atom in self.atoms]
        ret['irotat'] = [atom.irotat for atom in self.atoms]
        ret['rmin'] = [atom.rmin for atom in self.atoms]
        ret['epsilon'] = [atom.epsilon for atom in self.atoms]
        ret['rmin_14'] = [atom.rmin_14 for atom in self.atoms]
        ret['epsilon_14'] = [atom.epsilon_14 for atom in self.atoms]
        ret['resname'] = [atom.residue.name for atom in self.atoms]

        # Now for optional attributes
        # Coordinates
        try:
            coords = pd.DataFrame(
                    [[atom.xx, atom.xy, atom.xz] for atom in self.atoms],
                    columns=['xx', 'xy', 'xz']
            )
        except AttributeError:
            pass
        else:
            ret = ret.join(coords)
        # Velocities
        try:
            vels = pd.DataFrame(
                    [[atom.vx, atom.vy, atom.vz] for atom in self.atoms],
                    columns=['vx', 'vy', 'vz']
            )
        except AttributeError:
            pass
        else:
            ret = ret.join(vels)
        return ret

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PatchTemplate(ResidueTemplate):
    """
    A residue patch (typically used for CHARMM) that is used to modify existing
    residues in some way (e.g., terminal patches, disulfide bridges, etc.)

    Parameters
    ----------
    name : str, optional
        If provided, this is the name of the residue

    Attributes
    ----------
    delete : list of str
        List of atoms that need to be deleted in applying the patch

    See Also
    --------
    :class:`ResidueTemplate`

    Notes
    -----
    This class basically just provides an additional list of atoms that need to
    be deleted when applying this patch -- something that does not apply to
    standard Residues
    """
    def __init__(self, name=''):
        super(PatchTemplate, self).__init__(name)
        self.delete = []

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ResidueTemplateContainer(list):
    """
    A container of ResidueTemplate objects representing a unit with multiple
    residues

    Parameters
    ----------
    name : str, optional
        The name of the residue container
    """
    def __init__(self, name=''):
        self.box = None
        self.name = name

    @classmethod
    def from_structure(cls, struct, term_decorate=True):
        """
        Instantiates a ResidueTemplateContainer from a Structure instance filled
        with residues

        Parameters
        ----------
        struct : :class:`parmed.structure.Structure`
            The structure from which to generate the ResidueTemplateContainer
            from
        term_decorate : bool, optional
            If True, terminal amino and nucleic acid residues will be adorned as
            follows:

                * N-prepended if it is an N-terminal amino acid
                * C-prepended if it is a C-terminal amino acid
                * 5-appended if it is a 5'-terminal nucleic acid
                * 3-appended if it is a 3'-terminal nucleic acid

            For example, an N-terminal GLY will become NGLY, while a 5'-terminal
            DA will become DA5. Default is True
        """
        inst = cls()
        for res in struct.residues:
            rt = ResidueTemplate.from_residue(res)
            # See if we need to decorate the termini names
            if rt.head is None and rt.tail is not None and term_decorate:
                if AminoAcidResidue.has(rt.name):
                    rt.name = 'N%s' % rt.name
                elif RNAResidue.has(rt.name) or DNAResidue.has(rt.name):
                    rt.name = '%s5' % rt.name
            elif rt.tail is None and rt.head is not None and term_decorate:
                if AminoAcidResidue.has(rt.name):
                    rt.name = 'C%s' % rt.name
                elif RNAResidue.has(rt.name) or DNAResidue.has(rt.name):
                    rt.name = '%s3' % rt.name
            inst.append(rt)
        return inst

    def __getitem__(self, value):
        if isinstance(value, str):
            # Behave like a dict here... albeit a slow one
            for res in self:
                if res.name == value: return res
        return list.__getitem__(self, value)

    def to_library(self):
        """
        Converts the ResidueTemplateContainer instance to a library of unique
        :class:`ResidueTemplate` instances. The first of each kind of residue is
        taken

        Returns
        -------
        residues : dict {str : :class:`ResidueTemplate}
            The residue library with all residues from this residue collection
        """
        ret = dict()
        for res in self:
            if res.name in ret: continue
            ret[res.name] = res
        return ret
