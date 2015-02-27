"""
This module contains objects that deal with system topology and parameters, such
as atoms, residues, bonds, angles, etc.

by Jason Swails
"""
from __future__ import division

from chemistry.exceptions import (BondError, DihedralError, CmapError,
                                  AmoebaError, MissingParameter)
from chemistry.constants import TINY, DEG_TO_RAD, RAD_TO_DEG
from chemistry.periodic_table import Mass, Element as _Element
import chemistry.unit as u
from compat24 import all, property
import copy
import math
import warnings
try:
    from itertools import izip as zip
except ImportError:
    pass # Must be Python 3... zip _is_ izip

__all__ = ['Angle', 'AngleType', 'Atom', 'AtomList', 'Bond', 'BondType',
           'ChiralFrame', 'Cmap', 'CmapType', 'Dihedral', 'DihedralType',
           'DihedralTypeList', 'Improper', 'ImproperType', 'MultipoleFrame',
           'OutOfPlaneBend', 'PiTorsion', 'Residue', 'ResidueList',
           'StretchBend', 'StretchBendType', 'TorsionTorsion',
           'TorsionTorsionType', 'TrigonalAngle', 'TrackedList', 'UreyBradley',
           'OutOfPlaneBendType', 'NonbondedException', 'NonbondedExceptionType',
           'AcceptorDonor', 'Group', 'AtomType', 'NoUreyBradley', 'ExtraPoint',
           'TwoParticleExtraPointFrame', 'ThreeParticleExtraPointFrame',
           'OutOfPlaneExtraPointFrame']

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Private classes and methods

class _ListItem(object):
    """
    Helpful methods for items that appear in some kind of tracked list.
    Subclasses of this method interact with their `list` attribute, if they have
    one, in order to determine _where_ in the list they occur.

    Attributes
    ----------
    idx : ``int``
        This is intended to be a read-only variable that determines where in the
        list this particular object is. If there is no `list` attribute for this
        object, or the item is not in the list at all, `idx` is -1
        
    Notes
    -----
    For lists that support indexing its members and tracking when that list is
    changed (so we know when to update indexes), this is a fast lookup, as
    indexing only needs to be done once, at most, for the entire list (until
    changes are made that could change the indexing).

    The ``idx`` lookup in a TrackedList is therefore an O(1) operation, while
    it is O(N) for standard containers (O(N^2) for looking up *all* indexes).
    """

    @property
    def idx(self):
        try:
            mylist = self.list
        except AttributeError:
            return -1

        if mylist is None:
            return -1

        try:
            needs_indexing = mylist.needs_indexing
        except AttributeError:
            # This isn't a tracked list, so just look through the list
            for i, item in enumerate(mylist):
                if item is self: return i
            return -1
        else:
            if needs_indexing or self._idx == -1:
                try:
                    self._idx = -1
                    mylist.index_members()
                    return self._idx
                except AttributeError:
                    for i, item in enumerate(mylist):
                        if item is self: return i
                    return -1
            else:
                return self._idx

class _FourAtomTerm(object):
    """
    A base class for a parameter that spans 4 atoms

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom in the term
    atom2 : :class:`Atom`
        The second atom in the term
    atom3 : :class:`Atom`
        The third atom in the term
    atom4 : :class:`Atom`
        The fourth atom in the term

    Notes
    -----
    This is a base class that should be overridden by terms consisting of 4
    atoms (like torsions). There are a lot of such terms in the Amoeba force
    field, so this functionality is broken out into a new form
    """
    def __init__(self, atom1, atom2, atom3, atom4):
        if (atom1 is atom2 or atom1 is atom3 or atom1 is atom4 or
            atom2 is atom3 or atom2 is atom4 or atom3 is atom4):
            raise BondError('4-atom term cannot have duplicate atoms')
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

    def __contains__(self, thing):
        return (self.atom1 is thing or self.atom2 is thing or
                self.atom3 is thing or self.atom4 is thing)

    def delete(self):
        """ Sets all atoms in this term to None, and its type if it exists """
        self.atom1 = self.atom2 = self.atom3 = self.atom4 = None
        if hasattr(self, 'type'): self.type = None

class _ParameterType(object):
    """
    A parameter type that defines the nature of a particular molecular
    interaction. This is a base class that simply indicates whether a particular
    parameter type is "used", which in turn is used to determine whether or not
    this parameter type needs to be printed

    Attributes
    ----------
    used : ``bool``
        If ``True``, then this parameter type should be considered *used*. If
        ``False``, it is not being used and does not need to be printed.
    """

    def __init__(self):
        self.used = False

def _delete_from_list(list, item):
    """
    Deletes a requested item from a list. If the item does not exist in the
    list, a ValueError is raised

    Parameters
    ----------
    list : ``list``
        The list from which an item will be deleted
    item : ``object``
        The object to delete from the list
    """
    list.pop(list.index(item))

def _safe_assigns(dest, source, attrs):
    """
    Shallow-copies all requested attributes from `source` to `dest` if they
    are present. If not present, nothing is done
    """
    for attr in attrs:
        if not hasattr(source, attr): continue
        myattr = getattr(source, attr)
        setattr(dest, attr, myattr)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Atom(_ListItem):
    """ 
    An atom. Only use these as elements in AtomList instances, since AtomList
    will keep track of when indexes and other stuff needs to be updated. All
    parameters are optional.

    Parameters
    ----------
    atomic_number : ``int``
        The atomic number of this atom
    name : ``str``
        The name of this atom
    type : ``str``
        The type name of this atom
    charge : ``float``
        The partial atomic charge of this atom in fractions of an electron
    mass : ``float``
        The atomic mass of this atom in daltons
    nb_idx : ``int``
        The nonbonded index. This is a pointer that is relevant in the context
        of an Amber topology file and identifies its Lennard-Jones atom type
    radii : ``float``
        The intrinsic solvation radius of this atom.
    screen : ``float``
        The Generalized Born screening factor for this atom.
    occupancy : ``float``
        The occupancy of the atom (see PDB file)
    bfactor : ``float``
        The B-factor of the atom (see PDB file)
    altloc : ``str``
        Alternate location indicator (see PDB file)

    Other Parameters
    ----------------
    list : :class:`AtomList`
        The AtomList that this atom belongs to. If None, this atom does not
        belong to any list. This can be any iterable, but should typically be an
        AtomList or None
    tree : ``str``
        The tree chain identifier assigned to this atom. Relevant in the context
        of an Amber topology file, and not used for very much.
    join : ``int``
        The 'join` property of atoms stored in the Amber topology file. At the
        time of writing this class, `join` is unused, but still oddly required.
        Add support for future-proofing
    irotat : ``int``
        The `irotat` property of atoms stored in the Amber topology file.
        Unused, but included for future-proofing.
    number : ``int``
        The serial number given to the atom (see PDB file)
    rmin : ``float``
        The Rmin/2 Lennard-Jones parameter for this atom. Default evaluates to 0
    epsilon : ``float``
        The epsilon (well depth) Lennard-Jones parameter for this atom. Default
        evaluates to 0
    rmin14 : ``float``
        The Rmin/2 Lennard-Jones parameter for this atom in 1-4 interactions.
        Default evaluates to 0
    epsilon14 : ``float``
        The epsilon (well depth) Lennard-Jones parameter for this atom in 1-4
        interactions. Default evaluates to 0

    Other Attributes
    ----------------
    element : ``int``
        This is an alias for atomic_number
    atom_type : :class:`AtomType`
        In some cases, "type" is an ambiguous choice of an integer serial number
        or a string descriptor. In this case, atom_type is an AtomType instance
        that disambiguates this discrepancy.
    anisou : ``numpy.ndarray(float64) (or list of floats)``
        Anisotropic temperature scaling factors. This is a 6-element numpy array
        (if numpy is not available, it is a 6-element list) of floating point
        numbers. They are the 3x3 symmetric matrix elements U(1,1), U(2,2),
        U(3,3), U(1,2), U(1,3), U(2,3). If no factors available, it is None.
    idx : ``int``
        The index of this atom in the list. Set to -1 if this atom is not part
        of a list or the index cannot otherwise be determined (i.e., if the
        containing list does not support indexing its members)
    residue : :class:`Residue`
        The Residue that this atom belongs to. This is assigned when this atom
        is passed to `Residue.add_atom` -- see below for more information. Until
        it is set there, it is None
    other_locations : ``dict`` of :class:`Atom`
        A dict of Atom instances that represent alternate conformers of this
        atom. The keys are the `altloc` characters for those Atoms.
    bonds : ``list`` of :class:`Bond`
        list of Bond objects in which this atom is a member. This attribute
        should not be modified.
    angles : ``list`` of :class:`Angle`
        list of Angle objects in which this atom is a member. This attribute
        should not be modified.
    dihedrals : ``list`` of :class:`Dihedral`
        list of Dihedral objects in which this atom is a member. This attribute
        should not be modified.
    urey_bradleys : ``list`` of :class:`UreyBradley`
        list of UreyBradley objects in which this atom is a member (CHARMM,
        AMOEBA). This attribute should not be modified.
    impropers : ``list`` of :class:`Improper`
        list of Improper objects in which the atom is a member (CHARMM). This
        attribute should not be modified.
    cmaps : ``list`` of :class:`Cmap`
        list of Cmap objects in which the atom is a member (CHARMM, AMOEBA).
        This attribute should not be modified.
    tortors : ``list`` of :class:`TorsionTorsion`
        list of TorsionTorsion objects in which the atom is a member (AMOEBA).
        This attribute should not be modified.
    bond_partners : ``list`` of :class:`Atom`
        list of Atoms to which this atom is bonded. Do not modify this
        attribute -- it will almost certainly not do what you think it will
    angle_partners : ``list`` of :class:`Atom`
        list of Atoms to which this atom forms an angle, but not a bond. Do not
        modify this attribute -- it will almost certainly not do what you think
        it will
    dihedral_partners : ``list`` of :class:`Atom`
        list of Atoms to which this atom forms an dihedral, but not a bond or
        angle. Do not modify this attribute -- it will almost certainly not do
        what you think it will
    tortor_partners : ``list`` of :class:`Atom`
        list of Atoms to which this atom forms a coupled Torsion-Torsion, but
        not a bond or angle (AMOEBA). Do not modify this attribute -- it will
        almost certainly not do what you think it will
    exclusion_partners : ``list`` of :class:`Atom`
        list of Atoms with which this atom is excluded, but not bonded, angled,
        or dihedraled to. Do not modify this attribute -- it will almost
        certainly not do what you think it will
    marked : ``int``
        Mainly for internal use, it is used to indicate when certain atoms have
        been "marked" when traversing the bond network identifying topological
        features (like molecules and rings)
    children : ``list`` of :class:`ExtraPoint`
        This is the list of "child" ExtraPoint objects bonded to this atom
    number : ``int``
        The serial number of the atom in the input structure (e.g., PDB file).
        If not present in the original structure, a default value of -1 is used
    rmin : ``float``
        The Rmin/2 Lennard-Jones parameter for this atom. Default value is 0.
        If not set, it is taken from the `atom_type` attribute (if available)
    epsilon : ``float``
        The epsilon (well depth) Lennard-Jones parameter for this atom. Default
        value is 0. If not set, it is taken from the `atom_type` attribute (if
        available)
    rmin_14 : ``float``
        The Rmin/2 L-J parameter for 1-4 pairs. Default value is `rmin` (see
        above). If not set, it is taken from the `atom_type` attribute (if
        available).
    epsilon_14 : ``float``
        The epsilon L-J parameter for 1-4 pairs. Default value is `epsilon` (see
        above). If not set, it is taken from the `atom_type` attribute (if
        available).

    Possible Attributes
    -------------------
    xx : ``float``
        The X-component of the position of this atom. Only present if
        coordinates have been loaded. Otherwise raises AttributeError
    xy : ``float``
        The Y-component of the position of this atom. Only present if
        coordinates have been loaded. Otherwise raises AttributeError
    xz : ``float``
        The Z-component of the position of this atom. Only present if
        coordinates have been loaded. Otherwise raises AttributeError
    vx : ``float``
        The X-component of the velocity of this atom. Only present if the
        velocities have been loaded. Otherwise raises AttributeError
    vy : ``float``
        The Y-component of the velocity of this atom. Only present if the
        velocities have been loaded. Otherwise raises AttributeError
    vz : ``float``
        The Z-component of the velocity of this atom. Only present if the
        velocities have been loaded. Otherwise raises AttributeError
    type_idx : ``int``
        The AMOEBA atom type index. Only present if initialized with the AMOEBA
        force field. Otherwise raises AttributeError.
    class_idx : ``int``
        The AMOEBA class type index Only present if initialized with the AMOEBA
        force field. Otherwise raises AttributeError.
    multipoles : ``list(float)``
        The list of the 10 multipole moments up through quadrupoles Only present
        if initialized with the AMOEBA force field. Otherwise raises
        AttributeError.
    polarizability : ``list(float)``
        The polarizability of the atom. Only present if initialized with the
        AMOEBA force field. Otherwise raises AttributeError.
    vdw_parent : :class:`Atom`
        In the AMOEBA force field, this is the parent atom for the van der Waals
        term
    vdw_weight : ``float``
        In the AMOEBA force field, this is the weight of the van der Waals
        interaction on the parent atom
    segid : ``str``
        In CHARMM PDB and PSF files, the SEGID behaves similarly to the residue
        chain ID and is used to separate the total system into representative
        parts. This will only be set if read in from the input structure.

    Notes
    -----
    The bond_partners, angle_partners, dihedral_partners, and exclusion_partners
    arrays are actually generated as properties by taking differences of sets
    and sorting them. As a result, accessing this attribute constantly can be
    less efficient than you would expect. Iterating over them in a loop requires
    minimal overhead. But if frequent access is needed and these sets are
    guaranteed not to change, you should save a reference to the object and use
    that instead.

    Binary comparisons are done by atom index and are used primarily for sorting
    sublists of atoms. The == operator is not defined, so Atom equality should
    be done using the `is` operator. This allows Atom instances to be hashable
    (and so used as `dict` keys and put in `set`s)

    Examples
    --------
    >>> a1 = Atom(name='CO', type='C', charge=0.5, mass=12.01)
    >>> a2 = Atom(name='OC', type='O', charge=-0.5, mass=12.01)
    >>> a1.bond_to(a2)
    >>> a1 in a2.bond_partners and a2 in a1.bond_partners
    True
    >>> a1.idx # Not part of a container
    -1

    This object also supports automatic indexing when it is part of a container

    >>> atom_list = []
    >>> atom_list.append(Atom(list=atom_list, name='CO', charge=0.5))
    >>> atom_list.append(Atom(list=atom_list, name='OC', charge=-0.5))
    >>> atom_list[0].idx
    0
    >>> atom_list[1].idx
    1
    """
    #===================================================

    def __init__(self, list=None, atomic_number=0, name='', type='',
                 charge=0.0, mass=0.0, nb_idx=0, radii=0.0, screen=0.0,
                 tree='BLA', join=0.0, irotat=0.0, occupancy=0.0,
                 bfactor=0.0, altloc='', number=-1, rmin=None, epsilon=None,
                 rmin14=None, epsilon14=None):
        self.list = list
        self._idx = -1
        self.atomic_number = atomic_number
        self.name = name.strip()
        try:
            self.type = type.strip()
        except AttributeError:
            self.type = type
        self.charge = charge
        self.mass = mass
        self.nb_idx = nb_idx
        self.radii = radii
        self.screen = screen
        self.tree = tree
        self.join = join
        self.irotat = irotat
        self.bfactor = bfactor
        self.altloc = altloc
        self.occupancy = occupancy
        self._bond_partners = []
        self._angle_partners = []
        self._dihedral_partners = []
        self._tortor_partners = []
        self._exclusion_partners = [] # For arbitrary/other exclusions
        self.residue = None
        self.marked = 0 # For setting molecules
        self.bonds, self.angles, self.dihedrals = [], [], []
        self.urey_bradleys, self.impropers, self.cmaps = [], [], []
        self.tortors = []
        self.other_locations = {} # A dict of Atom instances
        self.atom_type = _UnassignedAtomType
        self.number = number
        self.anisou = None
        self._rmin = rmin
        self._epsilon = epsilon
        self._rmin14 = rmin14
        self._epsilon14 = epsilon14
        self.children = []
   
    #===================================================

    @classmethod
    def _copy(cls, item):
        new = cls(atomic_number=item.atomic_number, name=item.name,
                  type=item.type, charge=item.charge, mass=item.mass,
                  nb_idx=item.nb_idx, radii=item.radii,
                  screen=item.screen, tree=item.tree, join=item.join,
                  irotat=item.irotat, occupancy=item.occupancy,
                  bfactor=item.bfactor, altloc=item.altloc)
        new.atom_type = item.atom_type
        for key in item.other_locations:
            new.other_locations[key] = copy.copy(item.other_locations[key])
        _safe_assigns(new, item, ('xx', 'xy', 'xz', 'vx', 'vy', 'vz',
                      'type_idx', 'class_idx', 'multipoles', 'polarizability',
                      'vdw_parent', 'vdw_weight'))
        return new

    def __copy__(self):
        """ Returns a deep copy of this atom, but not attached to any list """
        return type(self)._copy(self)

    #===================================================

    # To alert people to changes in the API
    @property
    def starting_index(self):
        warnings.warn('starting_index has been replaced by idx',
                      DeprecationWarning)
        return self.idx

    @property
    def atname(self):
        warnings.warn('atname has been replaced by name', DeprecationWarning)
        return self.name
    @atname.setter
    def atname(self, thing):
        warnings.warn('atname has been replaced by name', DeprecationWarning)
        self.name = thing

    @property
    def attype(self):
        warnings.warn('attype has been replaced by type', DeprecationWarning)
        return self.type
    @attype.setter
    def attype(self, thing):
        warnings.warn('attype has been replaced by type', DeprecationWarning)
        self.type = thing

    #===================================================

    @property
    def bond_partners(self):
        """ Go through all bonded partners """
        bp = set(self._bond_partners)
        for p in self._bond_partners:
            for c in p.children:
                bp.add(c)
        return sorted(list(bp))

    @property
    def angle_partners(self):
        """ List of all angle partners that are NOT bond partners """
        bp = set(self._bond_partners)
        ap = set(self._angle_partners) - bp
        for p in ap:
            for c in p.children:
                ap.add(c)
        return sorted(list(ap))

    @property
    def dihedral_partners(self):
        " List of all dihedral partners that are NOT angle or bond partners "
        bp = set(self._bond_partners)
        ap = set(self._angle_partners)
        dp = set(self._dihedral_partners) - ap - bp
        for p in dp:
            for c in p.children:
                dp.add(c)
        return sorted(list(dp))

    @property
    def tortor_partners(self):
        """
        List of all 1-5 partners that are NOT in angle or bond partners. This is
        *only* used in the Amoeba force field
        """
        bp = set(self._bond_partners)
        ap = set(self._angle_partners)
        dp = set(self._dihedral_partners)
        tp = set(self._tortor_partners) - dp - ap - bp
        for p in tp:
            for c in p.children:
                tp.add(c)
        return sorted(list(tp))

    @property
    def exclusion_partners(self):
        """
        List of all exclusions not otherwise excluded by bonds/angles/torsions
        """
        bp = set(self._bond_partners)
        ap = set(self._angle_partners)
        dp = set(self._dihedral_partners)
        tp = set(self._tortor_partners)
        ep = set(self._exclusion_partners) - tp - dp - ap - bp
        for p in ep:
            for c in p.children:
                ep.add(c)
        return sorted(list(ep))

    #===================================================

    # Lennard-Jones parameters... can be taken from atom_type if it is set.
    # Otherwise take it from _rmin and _epsilon attributes

    @property
    def rmin(self):
        """ Lennard-Jones Rmin/2 parameter (the to Lennard-Jones radius) """
        if self._rmin is None:
            if (self.atom_type is _UnassignedAtomType or
                    self.atom_type.rmin is None):
                return 0.0
            return self.atom_type.rmin
        return self._rmin

    @rmin.setter
    def rmin(self, value):
        """ Lennard-Jones Rmin/2 parameter (the Lennard-Jones radius) """
        self._rmin = value

    @property
    def epsilon(self):
        """ Lennard-Jones epsilon parameter (the Lennard-Jones well depth) """
        if self._epsilon is None:
            if (self.atom_type is _UnassignedAtomType or
                    self.atom_type.epsilon is None):
                return 0.0
            return self.atom_type.epsilon
        return self._epsilon

    @epsilon.setter
    def epsilon(self, value):
        """ Lennard-Jones epsilon parameter (the Lennard-Jones well depth) """
        self._epsilon = value

    @property
    def rmin_14(self):
        """ The 1-4 Lennard-Jones Rmin/2 parameter """
        if self._rmin14 is None:
            if (self.atom_type is _UnassignedAtomType or
                    self.atom_type.rmin_14 is None):
                return self.rmin
            return self.atom_type.rmin_14
        return self._rmin14

    @rmin_14.setter
    def rmin_14(self, value):
        """ The 1-4 Lennard-Jones Rmin/2 parameter """
        self._rmin14 = value

    @property
    def epsilon_14(self):
        """ The 1-4 Lennard-Jones epsilon parameter """
        if self._epsilon14 is None:
            if (self.atom_type is _UnassignedAtomType or
                    self.atom_type.epsilon_14 is None):
                return self.epsilon
            return self.atom_type.epsilon_14
        return self._epsilon14

    @epsilon_14.setter
    def epsilon_14(self, value):
        """ The 1-4 Lennard-Jones Rmin/2 parameter """
        self._rmin14 = value

    #===================================================

    def nonbonded_exclusions(self, only_greater=True, index_from=0):
        """
        Returns the total number of nonbonded atom exclusions for this atom. The
        general rules for building the exclusion list is to include both
        exceptions AND exclusions (i.e., the Amber scaling of 1-4 interactions
        means that the 1-4 terms are excluded and a special pairlist is built to
        handle those exceptions).

        All atoms in the `_partners` arrays are nonbonded exclusions.

        Parameters
        ----------
        only_greater : ``bool``
            If True, only atoms whose `idx` value is greater than this `Atom`s
            `idx` will be counted as an exclusion (to avoid double-counting
            exclusions). If False, all exclusions will be counted.
        index_from : ``int``
            This is the index of the first atom, and is intended to be 0 (for C-
            and Python-style numbering) or 1 (for Fortran-style numbering, such
            as that used in the Amber and CHARMM topology files)

        Returns
        -------
        ``list of int``
            The returned list will be the atom indexes of the exclusion partners
            for this atom (indexing starts from ``index_from``)

        Notes
        -----
        If this instance's `idx` attribute evaluates to -1 -- meaning it is not
        in an AtomList -- a IndexError will be raised. If you have two extra
        points (i.e., those with atomic numbers of 0) bonded to each other, this
        routine may raise a ``RuntimeError`` if the recursion depth is exceeded.
        """
        if self.idx < 0:
            raise IndexError('Cannot find exclusions of an unindexed Atom')
        if only_greater:
            baseline = self.idx
        else:
            baseline = -1
        excl = []
        for atm in self.bond_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        for atm in self.angle_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        for atm in self.dihedral_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        for atm in self.tortor_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        for atm in self.exclusion_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        return sorted(excl)

    #===================================================

    # Make 'element' an alias for 'atomic_number'

    @property
    def element(self):
        return self.atomic_number
    @element.setter
    def element(self, value):
        self.atomic_number = value

    #===================================================

    def bond_to(self, other):
        """
        Log this atom as bonded to another atom.
        
        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `bond_partners`

        Notes
        -----
        This action adds `self` to `other.bond_partners`. Raises
        :class:`BondError` if `other is self`
        """
        if isinstance(other, ExtraPoint):
            self.children.append(other)
        elif isinstance(self, ExtraPoint):
            other.children.append(self)
        if self is other:
            raise BondError("Cannot bond atom to itself!")
        self._bond_partners.append(other)
        other._bond_partners.append(self)

    #===================================================
      
    def angle_to(self, other):
        """
        Log this atom as angled to another atom.

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `angle_partners`

        Notes
        -----
        This action adds `self` to `other.angle_partners`. Raises
        :class:`BondError` if `other is self`
        """
        if self is other:
            raise BondError("Cannot angle an atom with itself!")
        self._angle_partners.append(other)
        other._angle_partners.append(self)
   
    #===================================================

    def dihedral_to(self, other):
        """
        Log this atom as dihedral-ed to another atom.
        
        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `dihedral_partners`

        Notes
        -----
        This action adds `self` to `other.dihedral_partners`. Raises
        :class:`BondError` if `other is self`
        """
        if self is other:
            raise BondError("Cannot dihedral an atom with itself!")
        self._dihedral_partners.append(other)
        other._dihedral_partners.append(self)
      
    #===================================================

    def tortor_to(self, other):
        """
        Log this atom as 1-5 partners to another atom

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `tortor_partners`

        Notes
        -----
        This action adds `self` to `other.tortor_partners`. Raises
        :class:`BondError` if `other is self`
        """
        if self is other:
            raise BondError('Cannot coupled-dihedral atom to itself')
        self._tortor_partners.append(self)
        other._tortor_partners.append(self)

    #===================================================

    def exclude(self, other):
        """
        Add one atom to my arbitrary exclusion list

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `exclusion_partners`.

        Notes
        -----
        This action adds `self` to `other.exclusion_partners`. Raises
        :class:`BondError` if `other is self`
        """
        if self is other:
            raise BondError("Cannot exclude an atom from itself")
        self._exclusion_partners.append(other)
        other._exclusion_partners.append(self)

    #===================================================

    def reset_topology(self, keep_exclusions=True):
        """
        Empties all of the bond, angle, and dihedral partners so they can be set
        up again with updated data.

        Parameters
        ----------
        keep_exclusions : ``bool``
            If True, the `exclusion_partners` array is left alone. If False,
            `exclusion_partners` is emptied as well.
        """
        self._bond_partners = []
        self._angle_partners = []
        self._dihedral_partners = []
        self._tortor_partners = []
        if not keep_exclusions:
            self._exclusion_partners = []

    #===================================================

    def prune_exclusions(self):
        """
        For extra points, the exclusion partners may be filled before the bond,
        angle, dihedral, and tortor partners. Since we don't want memory of
        these exclusions if any of those topological features were to break, we
        want to *remove* those from the exclusion list. This function makes sure
        that nothing in the bond, angle, dihedral, and tortor lists appears in
        the exclusion list.
        """
        excludes = (set(self._exclusion_partners) - set(self._tortor_partners) -
                    set(self._dihedral_partners) - set(self._angle_partners) -
                    set(self._bond_partners))
        self._exclusion_partners = sorted(list(excludes))

    #===================================================

    # Comparisons are done by comparing indexes

    def __gt__(self, other):
        return self.idx > other.idx

    def __lt__(self, other):
        return self.idx < other.idx

    def __ge__(self, other):
        return not self < other

    def __le__(self, other):
        return not self > other

    def __repr__(self):
        start = '<Atom %s [%d]' % (self.name, self.idx)
        if self.residue is not None and hasattr(self.residue, 'idx'):
            return start + '; In %s %d>' % (self.residue.name, self.residue.idx)
        elif self.residue is not None:
            return start + '; In %s>' % self.residue.name
        return start + '>'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ExtraPoint(Atom):
    """
    An extra point is a massless, virtual site that is used to extend the
    flexibility of partial atomic charge fitting. They can be used to, for
    instance, give atoms permanent multipoles in a fixed-charge format.

    However, virtual sites are massless particles whose position is fixed with
    respect to a particular frame of reference. As a result, they must be
    treated specially when running dynamics. This class extends the `Atom` class
    with extra functionality specific to these "Extra points" or "virtual
    sites". See the documentation for the :class:`Atom` class for more
    information.
    """
    def __init__(self, *args, **kwargs):
        super(ExtraPoint, self).__init__(*args, **kwargs)
        self._frame_type = None

    #===================================================

    @property
    def parent(self):
        """
        The parent atom of this extra point is the atom it is bonded to (there
        should only be 1 bond partner, but the parent is defined as its first)
        """
        try:
            return self._bond_partners[0]
        except IndexError:
            return None

    #===================================================

    @property
    def bond_partners(self):
        """ List of all atoms bonded to this atom and its parent """
        try:
            return sorted([self.parent] +
                    [x for x in self.parent.bond_partners if x is not self])
        except AttributeError:
            if self.parent is not None:
                return [self.parent]
            return []

    @property
    def angle_partners(self):
        """ List of all atoms angled to this atom and its parent """
        try:
            return self.parent.angle_partners
        except AttributeError:
            return []

    @property
    def dihedral_partners(self):
        try:
            return self.parent.dihedral_partners
        except AttributeError:
            return []

    @property
    def tortor_partners(self):
        try:
            return self.parent.tortor_partners
        except AttributeError:
            return []

    @property
    def exclusion_partners(self):
        try:
            return self.parent.exclusion_partners
        except AttributeError:
            return []

    #===================================================

    @property
    def frame_type(self):
        """
        The type of frame used for this extra point. Whether the EP lies
        beyond (outside) or between (inside) the is determined from the
        positions of each atom in the frame and the extra point. If those are
        not set, the EP is assumed to be between the two points in the frame
        """
        if self._frame_type is not None:
            return self._frame_type
        if len(self.parent.bonds) == 2:
            mybond = None
            otherbond = None
            for bond in self.parent.bonds:
                if self in bond:
                    mybond = bond
                else:
                    otherbond = bond
            if mybond is None or otherbond is None:
                raise RuntimeError('Strange bond pattern detected')
            bonddist = otherbond.type.req
            if otherbond.atom1 is self.parent:
                otheratom = otherbond.atom2
            else:
                otheratom = otherbond.atom1
            try:
                x1, y1, z1 = self.xx, self.xy, self.xz
                x2, y2, z2 = otheratom.xx, otheratom.xy, otheratom.xz
            except AttributeError:
                self._frame_type = TwoParticleExtraPointFrame(self, True)
            else:
                dx, dy, dz = x1-x2, y1-y2, z1-z2
                if dx*dx + dy*dy + dz*dz > bonddist:
                    self._frame_type = TwoParticleExtraPointFrame(self, False)
                else:
                    self._frame_type = TwoParticleExtraPointFrame(self, True)
            return self._frame_type
        elif len(self.parent.bonds) == 3:
            # Just take 1 other atom -- an acute angle is "inside", an obtuse
            # angle is "outside"
            other_atom = None
            for bond in self.parent.bonds:
                if not self in bond:
                    if bond.atom1 is self.parent:
                        other_atom = bond.atom2
                    else:
                        other_atom = bond.atom1
                    break
            if other_atom is None:
                raise RuntimeError('Strange bond pattern detected')
            try:
                x1, y1, z1 = other_atom.xx, other_atom.xy, other_atom.xz
                x2, y2, z2 = self.parent.xx, self.parent.xy, self.parent.xz
                x3, y3, z3 = self.xx, self.xy, self.xz
            except AttributeError:
                # See if both other atoms are hydrogens and the current atom is
                # an oxygen. If so, this is TIP4P and we go inside. Otherwise,
                # we go outside
                other_atom2 = None
                for bond in self.parent.bonds:
                    if not self in bond and not other_atom in bond:
                        if bond.atom1 is self.parent:
                            other_atom2 = bond.atom2
                        else:
                            other_atom2 = bond.atom1
                        break
                if other_atom2 is None:
                    raise RuntimeError('Strange bond pattern detected')
                if (self.parent.element == 8 and other_atom.element == 1 and
                        other_atom2.element == 1): # TIP4P!
                    self._frame_type = ThreeParticleExtraPointFrame(self, True)
                else:
                    self._frame_type = ThreeParticleExtraPointFrame(self, False)
            else:
                vec1 = [x1-x2, y1-y2, z1-z2]
                vec2 = [x3-x2, y3-y2, z3-z2]
                dot11 = sum([x*y for x, y in zip(vec1, vec1)])
                dot22 = sum([x*y for x, y in zip(vec2, vec2)])
                dot12 = sum([x*y for x, y in zip(vec1, vec2)])
                angle = math.acos(dot12 / (math.sqrt(dot11)*math.sqrt(dot22)))
                if angle < math.pi / 2:
                    self._frame_type = ThreeParticleExtraPointFrame(self, True)
                else:
                    self._frame_type = ThreeParticleExtraPointFrame(self, False)
        elif len(self.parent.bonds) == 4:
            # Only supported option here is the OOP one (e.g., TIP5P) -- always
            # assume tetrahedral angle
            self._frame_type = OutOfPlaneExtraPointFrame(self)

        return self._frame_type

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TwoParticleExtraPointFrame(object):
    r"""
    This class defines a frame of reference for a given extra point with a frame
    of reference defined by 2 particles

    Parameters
    ----------
    ep : :class:`ExtraPoint`
        The extra point defined by this frame
    inside : ``bool``
        If True, the extra point is contained inside the curve connecting all
        points in the frame of reference. If False, the point is outside that
        curve. See figures below for example

    Figures
    -------
    In the figures, x marks the real particles, e marks the extra point
    (properly numbered). The numbers refer to the ordering this class expects
    those atoms to be in. The real particle that is the parent of the extra
    point is shown in upper-case

    Inside  : x1-----e--X2

    Outside : x1--------X2--e
    """
    def __init__(self, ep, inside=False):
        # Don't do error-checking now, since we want to give it time for the
        # bonds to be generated -- only error check when trying to get weights
        self.ep = ep
        self.inside = inside

    def get_atoms(self):
        """
        Returns the particles involved in the frame

        Returns
        -------
        a1, a2 : :class:`Atom`, :class:`Atom`
            a1 is the parent atom of the extra point. a2 is the other atom
            bonded to the parent atom
        """
        try:
            mybond, = [bond for bond in self.ep.parent.bonds
                                if self.ep not in bond]
        except TypeError:
            raise RuntimeError("Bad bond pattern in EP frame")

        if mybond.atom1 is self.ep.parent:
            return self.ep.parent, mybond.atom2
        else:
            return self.ep.parent, mybond.atom1

    def get_weights(self):
        """
        Returns the weights for the two particles

        Returns
        -------
        w1, w2 : ``float, float``
            w1 is the weight of particle 1 (the parent atom), whereas w2 is the
            weight of particle 2 (the atom bonded to the parent atom)
        """
        ep = self.ep
        if len(ep.parent.bonds) != 2:
            raise ValueError('EP parent bond pattern inconsistent with a 2-'
                             'point virtual site frame')
        b1, b2 = ep.parent.bonds
        if ep in b1:
            r1 = b2.type.req
            r2 = b1.type.req
        else:
            r1 = b2.type.req
            r2 = b1.type.req

        if self.inside:
            # It is closer to atom 1, but both weights are positive and add to 1
            return ((r1 - r2) / r1), (r2 / r1)
        else:
            return ((r1 + r2) / r1), -(r2 / r1)

class ThreeParticleExtraPointFrame(object):
    r"""
    This class defines a frame of reference for a given extra point with a frame
    of reference defined by 3 particles

    Parameters
    ----------
    ep : :class:`ExtraPoint`
        The extra point defined by this frame
    inside : ``bool``
        If True, the extra point is contained inside the curve connecting all
        points in the frame of reference. If False, the point is outside that
        curve. See figures below for example

    Figures
    -------
    In the figures, x marks the real particles, e marks the extra point
    (properly numbered). The numbers refer to the ordering this class expects
    those atoms to be in. The real particle that is the parent of the extra
    point is shown in upper-case

    Inside :         X
                    /|\
                   / e \
                  x     x
    Outside :
                     e
                     |
                     X
                    / \
                   /   \
                  x     x
    """
    def __init__(self, ep, inside=True):
        # Don't do error-checking now, since we want to give it time for the
        # bonds to be generated -- only error check when trying to get weights
        self.ep = ep
        self.inside = inside

    def get_atoms(self):
        """
        Returns the particles involved in the frame

        Returns
        -------
        a1, a2, a3 : :class:`Atom`, :class:`Atom`, :class:`Atom`
            a1 is the parent atom of the extra point. a2 and a3 are the other
            atoms that are both bonded to the parent atom
        """
        try:
            b1, b2 = [bond for bond in self.ep.parent.bonds
                                if self.ep not in bond]
        except TypeError:
            raise RuntimeError('Unsupported bonding pattern in EP frame')
        if b1.atom1 is self.ep.parent:
            oatom1 = b1.atom2
        else:
            oatom1 = b1.atom1
        if b2.atom1 is self.ep.parent:
            oatom2 = b2.atom2
        else:
            oatom2 = b2.atom1
        return self.ep.parent, oatom1, oatom2

    def get_weights(self):
        """
        Returns the weights for the three particles

        Returns
        -------
        w1, w2, w3 : ``float, float, float``
            w1 is the weight of particle 1 (the parent atom), whereas w2 and w3
            are the weights of particles 2 and 3 (the atoms bonded to the
            parent atom)
        """
        ep = self.ep
        if len(ep.parent.bonds) != 3:
            raise ValueError('EP parent bond pattern inconsistent with a 3-'
                             'point virtual site frame')
        b1, b2, b3 = ep.parent.bonds
        # There are 2 possibilities here -- there is an angle between the 3
        # atoms in the frame OR there is a triangle of bonds (e.g., TIP4P).
        # Compute the 'ideal' distance between the 3 frame-of-ref. atoms
        if ep in b1:
            b1, b3 = b3, b1
        elif ep in b2:
            b2, b3 = b3, b2
        # See if there is an angle with both b1 and b2 in it
        found = False
        for angle in ep.parent.angles:
            if b1 in angle and b2 in angle:
                found = True
                break
        if found:
            # Compute the 2-3 distance from the two bond lengths and the angles
            # using law of cosines
            r1 = b1.type.req
            r2 = b2.type.req
            theta = angle.type.theteq * DEG_TO_RAD
            req23 = math.sqrt(r1*r1 + r2*r2 - 2*r1*r2*math.cos(theta))
        else:
            # See if there is a bond between particles 2 and 3
            if b1.atom1 is ep.parent:
                a1 = b1.atom2
            else:
                a1 = b1.atom1
            if b2.atom1 is ep.parent:
                a2 = b2.atom2
            else:
                a2 = b2.atom1
            if a1 not in a2.bond_partners:
                raise RuntimeError('EP frame definition incomplete for 3-point '
                                   'virtual site... cannot determine distance '
                                   'between particles 2 and 3')
            req23 = None
            for bond in a1.bonds:
                if a2 in bond:
                    req23 = bond.type.req
            if req23 is None:
                raise RuntimeError('Cannot determine 2-3 particle distance in '
                                   'three-site virtual site frame')
        req12 = b1.type.req
        req13 = b2.type.req
        weight = b3.type.req / math.sqrt(req12*req13 - 0.25*req23*req23)

        if self.inside:
            return 1 - weight, weight / 2, weight / 2
        else:
            return 1 + weight, -weight / 2, -weight / 2

class OutOfPlaneExtraPointFrame(object):
    r"""
    This class defines a frame of reference for a given extra point with a frame
    of reference defined by 3 particles, but with the virtual site out of the
    plane of those 3 particles. For example, TIP5P

    Parameters
    ----------
    ep : :class:`ExtraPoint`
        The extra point defined by this frame
    angle : ``float``
        The angle out-of-plane that the extra point is. By default, it is half
        of the angle of a tetrahedron. Given in degrees

    Figures
    -------
    In the figures, x marks the real particles, e marks the extra point
    (properly numbered). The numbers refer to the ordering this class expects
    those atoms to be in. The real particle that is the parent of the extra
    point is shown in upper-case

                   e   e
                    \ /
                     X
                    / \
                   /   \
                  x     x
    """
    def __init__(self, ep, angle=54.735):
        # Don't do error-checking now, since we want to give it time for the
        # bonds to be generated -- only error check when trying to get weights
        self.ep = ep
        self.angle = angle

    def get_atoms(self):
        """
        Returns the particles involved in the frame

        Returns
        -------
        a1, a2, a3 : :class:`Atom`, :class:`Atom`, :class:`Atom`
            a1 is the parent atom of the extra point. a2 and a3 are the other
            atoms that are both bonded to the parent atom but are not EPs
        """
        try:
            b1, b2 = [bond for bond in self.ep.parent.bonds
                                if (not isinstance(bond.atom1, ExtraPoint) and
                                    not isinstance(bond.atom2, ExtraPoint))
            ]
        except TypeError:
            raise RuntimeError('Unsupported bonding pattern in EP frame')
        if b1.atom1 is self.ep.parent:
            oatom1 = b1.atom2
        else:
            oatom1 = b1.atom1
        if b2.atom1 is self.ep.parent:
            oatom2 = b2.atom2
        else:
            oatom2 = b2.atom1
        return self.ep.parent, oatom2, oatom1

    def get_weights(self):
        """
        Returns the weights for the three particles

        Returns
        -------
        w1, w2, w3 : ``float, float, float``
            w1 and w2 are the weights with respect to the second two particles
            in the frame (i.e., NOT the parent atom). w3 is the weight of the
            cross-product
        """
        ep = self.ep
        # Find the two bonds that do not involve any extra points
        regbonds = []
        mybond = None
        for bond in ep.parent.bonds:
            if (not isinstance(bond.atom1, ExtraPoint) and
                    not isinstance(bond.atom2, ExtraPoint)):
                regbonds.append(bond)
            elif ep in bond:
                if mybond is not None:
                    raise ValueError('multiple bonds detected to extra point')
                mybond = bond
        #FIXME -- is this necessary?
        if len(regbonds) != 2:
            raise ValueError('EP parent bond pattern inconsistent with an '
                             'out-of-plane, 3-point virtual site frame')
        if mybond is None:
            raise RuntimeError('No EP bond found... should not be here')
        b1, b2 = regbonds[:2]
        req12 = b1.type.req
        req13 = b2.type.req
        # See if there is an angle with both b1 and b2 in it
        found = False
        for angle in ep.parent.angles:
            if b1 in angle and b2 in angle:
                found = True
                break
        if found:
            # Compute the 2-3 distance from the two bond lengths and the angles
            # using law of cosines
            t213 = angle.theteq
            r1 = b1.type.req
            r2 = b2.type.req
            theta = angle.type.theteq * DEG_TO_RAD
            req23 = math.sqrt(r1*r1 + r2*r2 - 2*r1*r2*math.cos(theta))
        else:
            # See if there is a bond between particles 2 and 3
            if b1.atom1 is ep.parent:
                a1 = b1.atom2
            else:
                a1 = b1.atom1
            if b2.atom1 is ep.parent:
                a2 = b2.atom2
            else:
                a2 = b2.atom1
            if a1 not in a2.bond_partners:
                raise RuntimeError('EP frame definition incomplete for 3-point '
                                   'virtual site... cannot determine distance '
                                   'between particles 2 and 3')
            req23 = None
            for bond in a1.bonds:
                if a2 in bond:
                    req23 = bond.type.req
            if req23 is None:
                raise RuntimeError('Cannot determine 2-3 particle distance in '
                                   'three-site virtual site frame')
            # Now calculate the angle
            t213 = 2 * math.asin(req23 / (req12 + req13)) * RAD_TO_DEG
        # Some necessary constants
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        req12 *= length_conv
        req13 *= length_conv
        req23 *= length_conv
        sinOOP = math.sin(self.angle * DEG_TO_RAD)
        cosOOP = math.cos(self.angle * DEG_TO_RAD)
        sin213 = math.sin(t213 * DEG_TO_RAD)
        # Find how big the cross product is
        lenCross = req12 * req13 * sin213
        # Find our weights (assume symmetric)
        weightCross = sinOOP * mybond.type.req * length_conv / lenCross
        weight = (cosOOP * mybond.type.req * length_conv /
                        math.sqrt(req12*req13 - 0.25*req23*req23))
        # Find out of the cross product weight should be positive or negative.
        a1, a2, a3 = self.get_atoms()
        try:
            v12 = (a2.xx-a1.xx, a2.xy-a1.xy, a2.xz-a1.xz)
            v13 = (a3.xx-a1.xx, a3.xy-a1.xy, a3.xz-a1.xz)
            v1e = (self.ep.xx-a1.xx, self.ep.xy-a1.xy, self.ep.xz-a1.xz)
        except AttributeError:
            # No coordinates... we have to guess. The first EP will have a
            # positive weight, the second will have a negative (this matches
            # what happens with Amber prmtop files...)
            for a in self.ep.parent.bond_partners:
                if a is self.ep:
                    break
                if isinstance(a, ExtraPoint):
                    weightCross = -weightCross
                    break
        else:
            # Take the cross product of v12 and v13, then dot that with the EP
            # vector. An acute angle is a positive weightCross. An obtuse one is
            # negative
            cross = (v12[1]*v13[2] - v12[2]*v13[1],
                     v12[2]*v13[0] - v12[0]*v13[2],
                     v12[0]*v13[1] - v12[1]*v13[0])
            lencross = math.sqrt(sum([cross[i]*cross[i] for i in xrange(3)]))
            lenv1e = math.sqrt(sum([v1e[i]*v1e[i] for i in xrange(3)]))
            v1edotcross = sum([v1e[i]*cross[i] for i in xrange(3)])
            costheta = v1edotcross / (lenv1e*lencross)
            if costheta < 0: weightCross = -weightCross
        return weight / 2, weight / 2, weightCross

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Bond(object):
    """
    A covalent bond connecting two atoms.

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom involved in the bond
    atom2 : :class:`Atom`
        The other atom involved in the bond
    type : :class:`BondType`
        The bond type that defines the parameters for this bond

    Notes
    -----
    You can test whether an :class:`Atom` is contained within the bond using the
    `in` operator. A `BondError` is raised if `atom1` and `atom2` are identical.
    This bond instance is `append`ed to the `bonds` list for both `atom1` and
    `atom2` and is automatically removed from those lists upon garbage
    collection

    Examples
    --------
    >>> a1, a2 = Atom(), Atom()
    >>> bond = Bond(a1, a2)
    >>> a1 in bond and a2 in bond
    True
    """

    def __init__(self, atom1, atom2, type=None):
        """ Bond constructor """
        # Make sure we're not bonding me to myself
        if atom1 is atom2:
            raise BondError('Cannot bond atom to itself!')
        # Order the atoms so the lowest atom # is first
        self.atom1 = atom1
        self.atom2 = atom2
        # Load this struct into the atoms
        self.atom1.bonds.append(self)
        self.atom2.bonds.append(self)
        atom1.bond_to(atom2)
        self.type = type

    @property
    def bond_type(self):
        warnings.warn("bond_type has been replaced by type", DeprecationWarning)
        return self.type

    def __contains__(self, thing):
        """ Quick and easy way to see if an Atom is in this Bond """
        return thing is self.atom1 or thing is self.atom2

    def delete(self):
        """
        Deletes this bond from the atoms that make it up. This method removes
        this bond from the `bonds` list for both atom1 and atom2 as well as
        removing atom1 from atom2.bond_partners and vice versa.
        """
        _delete_from_list(self.atom1.bonds, self)
        _delete_from_list(self.atom2.bonds, self)
        _delete_from_list(self.atom1._bond_partners, self.atom2)
        _delete_from_list(self.atom2._bond_partners, self.atom1)

        self.atom1 = self.atom2 = self.type = None

    def __repr__(self):
        return '<%s %r--%r; type=%r>' % (type(self).__name__,
                self.atom1, self.atom2, self.type)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondType(_ListItem, _ParameterType):
    """
    A bond type with a set of bond parameters

    Parameters
    ----------
    k : ``float``
        Force constant in kcal/mol/Angstrom^2
    req : ``float``
        Equilibrium bond distance in Angstroms
    list : :class:`TrackedList`
        A list of :class:`BondType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : int
        The index of this BondType inside its containing list

    Notes
    -----
    Two `BondType`s are equal if their `k` and `req` attributes are equal

    Examples
    --------
    >>> bt1 = BondType(10.0, 1.0)
    >>> bt2 = BondType(10.0, 1.0)
    >>> bt1 is bt2
    False
    >>> bt1 == bt2
    True
    >>> bt1.idx # Not in a list or container
    -1

    As part of a list, they can be indexed

    >>> bond_list = []
    >>> bond_list.append(BondType(10.0, 1.0, list=bond_list))
    >>> bond_list.append(BondType(10.0, 1.0, list=bond_list))
    >>> bond_list[0].idx
    0
    >>> bond_list[1].idx
    1
    """

    def __init__(self, k, req, list=None):
        _ParameterType.__init__(self)
        self.k = k
        self.req = req
        self.list = list
        self._idx = -1

    def __eq__(self, other):
        return self.k == other.k and self.req == other.req

    def __repr__(self):
        return '<%s; k=%.3f, Req=%.3f>' % (type(self).__name__,
                self.k, self.req)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Angle(object):
    """
    A valence angle between 3 atoms separated by two covalent bonds.

    Parameters
    ----------
    atom1 : :class:`Atom`
        An atom one end of the valence angle
    atom2 : :class:`Atom`
        The atom in the middle of the valence angle bonded to both other atoms
    atom3 : :class:`Atom`
        An atom on the other end of the valence angle to atom1
    type : :class:`AngleType`
        The AngleType object containing the parameters for this angle

    Notes
    -----
    An Angle can contain bonds or atoms. A bond is contained if it exists
    between atoms 1 and 2 or atoms 2 and 3.

    Examples
    --------
    >>> a1, a2, a3 = Atom(), Atom(), Atom()
    >>> angle = Angle(a1, a2, a3)
    >>> Bond(a1, a2) in angle and Bond(a3, a2) in angle
    True
    >>> a1 in angle and a2 in angle and a3 in angle
    True
    >>> Bond(a1, a3) in angle # this is not part of the angle definition
    False
    """
      
    def __init__(self, atom1, atom2, atom3, type=None):
        # Make sure we're not angling me to myself
        if atom1 is atom2 or atom1 is atom3 or atom2 is atom3:
            raise BondError('Cannot angle atom to itself!')
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        # Log these angles in each atom
        atom1.angles.append(self)
        atom2.angles.append(self)
        atom3.angles.append(self)
        # Load the force constant and equilibrium angle
        self.type = type
        atom1.angle_to(atom2)
        atom1.angle_to(atom3)
        atom2.angle_to(atom3)

    @property
    def angle_type(self):
        warnings.warn("angle_type has been replaced with type",
                      DeprecationWarning)
        return self.type

    def __contains__(self, thing):
        """ Quick and easy way to see if an Atom or a Bond is in this Angle """
        if isinstance(thing, Atom):
            return (thing is self.atom1 or thing is self.atom2 or
                    thing is self.atom3)
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing))

    def delete(self):
        """
        Deletes this angle from the atoms that make it up. This method removes
        this angle from the `angles` list for atom1, atom2, and atom3 as well as
        removing each atom form each others' angle partner lists
        """
        _delete_from_list(self.atom1.angles, self)
        _delete_from_list(self.atom2.angles, self)
        _delete_from_list(self.atom3.angles, self)

        _delete_from_list(self.atom1._angle_partners, self.atom2)
        _delete_from_list(self.atom1._angle_partners, self.atom3)
        _delete_from_list(self.atom2._angle_partners, self.atom1)
        _delete_from_list(self.atom2._angle_partners, self.atom3)
        _delete_from_list(self.atom3._angle_partners, self.atom1)
        _delete_from_list(self.atom3._angle_partners, self.atom2)

        self.atom1 = self.atom2 = self.atom3 = self.type = None

    def __repr__(self):
        return '<%s %r--%r--%r; type=%r>' % (type(self).__name__,
                self.atom1, self.atom2, self.atom3, self.type)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleType(_ListItem, _ParameterType):
    """
    An angle type with a set of angle parameters

    Parameters
    ----------
    k : ``float``
        Force constant in kcal/mol/radians^2
    theteq : ``float``
        Equilibrium angle in Degrees
    list : :class:`TrackedList`
        A list of `AngleType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : ``int``
        The index of this AngleType inside its containing list

    Notes
    -----
    Two `AngleType`s are equal if their `k` and `theteq` attributes are equal

    Examples
    --------
    >>> at1 = AngleType(10.0, 180.0)
    >>> at2 = AngleType(10.0, 180.0)
    >>> at1 is at2
    False
    >>> at1 == at2
    True
    >>> at1.idx # not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> angle_list = []
    >>> angle_list.append(AngleType(10.0, 180.0, list=angle_list))
    >>> angle_list.append(AngleType(10.0, 180.0, list=angle_list))
    >>> angle_list[0].idx
    0
    >>> angle_list[1].idx
    1
    """
    def __init__(self, k, theteq, list=None):
        _ParameterType.__init__(self)
        self.k = k
        self.theteq = theteq
        self._idx = -1
        self.list = list

    def __eq__(self, other):
        return self.k == other.k and self.theteq == other.theteq

    def __repr__(self):
        return '<%s; k=%.3f, THETAeq=%.3f>' % (type(self).__name__,
                self.k, self.theteq)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Dihedral(_FourAtomTerm):
    """
    A valence dihedral between 4 atoms separated by three covalent bonds.

    Parameters
    ----------
    atom1 : :class:`Atom`
        An atom on one end of the valence dihedral bonded to atom 2
    atom2 : :class:`Atom`
        An atom in the middle of the valence dihedral bonded to atom1 and atom3
    atom3 : :class:`Atom`
        An atom in the middle of the valence dihedral bonded to atom2 and atom4
    atom4 : :class:`Atom`
        An atom on the other end of the valence dihedral bonded to atom 3
    improper : ``bool``
        If True, this is an Amber-style improper torsion, where atom3 is the
        "central" atom bonded to atoms 1, 2, and 4 (atoms 1, 2, and 4 are *only*
        bonded to atom 3 in this instance)
    ignore_end : ``bool``
        If True, the end-group interactions for this torsion are ignored, either
        because it is involved in a ring where the end-group atoms are excluded
        or because it is one term in a multi-term dihedral expansion
    type : :class:`DihedralType`
        The DihedralType object containing the parameters for this dihedral

    Notes
    -----
    A Dihedral can contain bonds or atoms. A bond is contained if it exists
    between atoms 1 and 2, between atoms 2 and 3, or between atoms 3 and 4.

    Examples
    --------
    >>> a1, a2, a3, a4 = Atom(), Atom(), Atom(), Atom()
    >>> dihed = Dihedral(a1, a2, a3, a4)
    >>> Bond(a1,a2) in dihed and Bond(a3,a2) in dihed and Bond(a3,a4) in dihed
    True
    >>> a1 in dihed and a2 in dihed and a3 in dihed and a4 in dihed
    True
    >>> Bond(a1, a4) in dihed # this is not part of the angle definition
    False

    For improper torsions, the bond pattern is different

    >>> a1, a2, a3, a4 = Atom(), Atom(), Atom(), Atom()
    >>> dihed = Dihedral(a1, a2, a3, a4, improper=True)
    >>> Bond(a1,a3) in dihed and Bond(a2,a3) in dihed and Bond(a3,a4) in dihed
    True
    >>> Bond(a1,a2) in dihed # Not like a normal dihedral!
    False
    >>> a1 in dihed and a2 in dihed and a3 in dihed and a4 in dihed
    True
    """
      
    def __init__(self, atom1, atom2, atom3, atom4, improper=False,
                 ignore_end=False, type=None):
        _FourAtomTerm.__init__(self, atom1, atom2, atom3, atom4)
        # improper _implies_ ignore_end
        ignore_end = improper or ignore_end
        # Log these dihedrals in each atom
        atom1.dihedrals.append(self)
        atom2.dihedrals.append(self)
        atom3.dihedrals.append(self)
        atom4.dihedrals.append(self)
        self.type = type
        self.improper = improper
        self.ignore_end = ignore_end
        self.signs = [1, 1]
        if ignore_end: self.signs[0] = -1
        if improper: self.signs[1] = -1
        if not improper:
            atom1.dihedral_to(atom2)
            atom1.dihedral_to(atom3)
            atom1.dihedral_to(atom4)
            atom2.dihedral_to(atom3)
            atom2.dihedral_to(atom4)
            atom3.dihedral_to(atom4)

    @property
    def dihed_type(self):
        warnings.warn('dihed_type has been replaced by type',
                      DeprecationWarning)
        return self.type

    def __contains__(self, thing):
        """
        Quick and easy way to find out if an Atom or Bond is in this Dihedral
        """
        if isinstance(thing, Atom):
            return _FourAtomTerm.__contains__(self, thing)
        # A dihedral is made up of 3 bonds
        if self.improper:
            # An improper is different... Atom 3 is the central atom
            return ((self.atom1 in thing and self.atom3 in thing) or
                    (self.atom2 in thing and self.atom3 in thing) or
                    (self.atom4 in thing and self.atom3 in thing))
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing) or
                (self.atom3 in thing and self.atom4 in thing))

    def same_atoms(self, thing):
        """
        Determines if this dihedral has the same atoms (or atom indexes) as
        another object

        Parameters
        ----------
        thing : :class:`Dihedral` or other ``iterable``
            A Dihedral or an iterable with 4 indexes that will be used to
            identify whether this torsion has the same atoms

        Returns
        -------
        is_same : ``bool``
            True if ``thing`` and ``self`` have the same atoms. ``False``
            otherwise

        Notes
        -----
        This raises a ``TypeError`` if thing is not a Dihedral and is not
        iterable
        """
        if isinstance(thing, Dihedral):
            # I'm comparing with another Dihedral here
            return ((self.atom1 is thing.atom1 and self.atom2 is thing.atom2 and
                     self.atom3 is thing.atom3 and self.atom4 is thing.atom4) or
                    (self.atom1 is thing.atom4 and self.atom2 is thing.atom3 and
                     self.atom4 is thing.atom1)
            )
        thing = list(thing)
        # Here, atoms are expected to index from 0 (Python standard) if we
        # are comparing with a list or tuple
        if len(thing) != 4:
            raise DihedralError('comparative %s has %d elements! Expect 4.'
                                % (type(thing).__name__, len(thing)))
        # Compare starting_index, since we may not have an index right now
        return ( (self.atom1.idx == thing[0] and 
                self.atom2.idx == thing[1] and
                self.atom3.idx == thing[2] and
                self.atom4.idx == thing[3]) or
                (self.atom1.idx == thing[3] and 
                self.atom2.idx == thing[2] and
                self.atom3.idx == thing[1] and
                self.atom4.idx == thing[0]) )

    def delete(self):
        """
        Deletes this dihedral from the atoms that make it up. This method
        removes this dihedral from the `dihedrals` list for atom1, atom2, atom3,
        and atom4 as well as removing each atom form each others' dihedral
        partner lists
        """
        _delete_from_list(self.atom1.dihedrals, self)
        _delete_from_list(self.atom2.dihedrals, self)
        _delete_from_list(self.atom3.dihedrals, self)
        _delete_from_list(self.atom4.dihedrals, self)

        if not self.improper:
            _delete_from_list(self.atom1._dihedral_partners, self.atom2)
            _delete_from_list(self.atom1._dihedral_partners, self.atom3)
            _delete_from_list(self.atom1._dihedral_partners, self.atom4)
            _delete_from_list(self.atom2._dihedral_partners, self.atom1)
            _delete_from_list(self.atom2._dihedral_partners, self.atom3)
            _delete_from_list(self.atom2._dihedral_partners, self.atom4)
            _delete_from_list(self.atom3._dihedral_partners, self.atom1)
            _delete_from_list(self.atom3._dihedral_partners, self.atom2)
            _delete_from_list(self.atom3._dihedral_partners, self.atom4)
            _delete_from_list(self.atom4._dihedral_partners, self.atom1)
            _delete_from_list(self.atom4._dihedral_partners, self.atom2)
            _delete_from_list(self.atom4._dihedral_partners, self.atom3)

        self.atom1 = self.atom2 = self.atom3 = self.atom4 = self.type = None

    def __repr__(self):
        if self.improper:
            name = '%s [imp]' % (type(self).__name__)
        elif self.ignore_end:
            name = '%s [ign]' % (type(self).__name__)
        else:
            name = type(self).__name__
        return '<%s; %r--%r--%r--%r; type=%r>' % (name, self.atom1,
                self.atom2, self.atom3, self.atom4, self.type)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralType(_ListItem, _ParameterType):
    """
    A dihedral type with a set of dihedral parameters

    Parameters
    ----------
    phi_k : ``float``
        The force constant in kcal/mol
    per : ``int``
        The dihedral periodicity
    phase : ``float``
        The dihedral phase in degrees
    scee : ``float``
        1-4 electrostatic scaling factor
    scnb : ``float``
        1-4 Lennard-Jones scaling factor
    list : :class:`TrackedList`
        A list of `DihedralType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : ``int``
        The index of this DihedralType inside its containing list

    Notes
    -----
    Two :class:`DihedralType`s are equal if their `phi_k`, `per`, `phase`,
    `scee`, and `scnb` attributes are equal

    Examples
    --------
    >>> dt1 = DihedralType(10.0, 2, 180.0, 1.2, 2.0)
    >>> dt2 = DihedralType(10.0, 2, 180.0, 1.2, 2.0)
    >>> dt1 is dt2
    False
    >>> dt1 == dt2
    True
    >>> dt1.idx # not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> dihedral_list = []
    >>> dihedral_list.append(DihedralType(10.0, 2, 180.0, 1.2, 2.0,
    ...                                   list=dihedral_list))
    >>> dihedral_list.append(DihedralType(10.0, 2, 180.0, 1.2, 2.0,
    ...                                   list=dihedral_list))
    >>> dihedral_list[0].idx
    0
    >>> dihedral_list[1].idx
    1
    """

    #===================================================
   
    def __init__(self, phi_k, per, phase, scee=1.0, scnb=1.0, list=None):
        """ DihedralType constructor """
        _ParameterType.__init__(self)
        self.phi_k = phi_k
        self.per = per
        self.phase = phase
        self.scee = scee
        self.scnb = scnb
        self.list = list
        self._idx = -1

   #===================================================

    def write_info(self, parm):
        """ Write out the dihedral parameters """
        # If our idx == -1, we're not being used, so just return
        if self.idx == -1: return
        parm.parm_data['DIHEDRAL_FORCE_CONSTANT'][self.idx] = self.phi_k
        parm.parm_data['DIHEDRAL_PERIODICITY'][self.idx] = self.per
        parm.parm_data['DIHEDRAL_PHASE'][self.idx] = self.phase
        try:
            parm.parm_data['SCEE_SCALE_FACTOR'][self.idx] = self.scee
        except KeyError:
            pass
        try:
            parm.parm_data['SCNB_SCALE_FACTOR'][self.idx] = self.scnb
        except KeyError:
            pass
   
    #===================================================

    def __eq__(self, other):
        return (self.phi_k == other.phi_k and self.per == other.per and 
                self.phase == other.phase and self.scee == other.scee and
                self.scnb == other.scnb)

    def __repr__(self):
        return ('<%s; k=%.3f, periodicity=%d, phase=%.3f, '
                'scee=%.3f, scnb=%.3f>' % (type(self).__name__, self.phi_k,
                    self.per, self.phase, self.scee, self.scnb))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralTypeList(list, _ListItem):
    """
    Dihedral types are a Fourier expansion of terms. In some cases, they are
    stored in a list like this one inside another TrackedList. In other cases,
    each term is a separate entry in the :class:`TrackedList`.

    In cases where `DihedralType`s are stored with every term in the same
    container, this object supports list assignment and indexing like
    :class:`DihedralType`.
    """
    def __init__(self, *args, **kwargs):
        list.__init__(self, *args, **kwargs)
        self._idx = -1
        self.list = None
        self.used = False

    def __eq__(self, other):
        if len(self) != len(other): return False
        for t1, t2 in zip(self, other):
            if not t1 == t2:
                return False
        return True

    def __repr__(self):
        return 'DihedralTypes %s' % (super(DihedralTypeList, self).__repr__())

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradley(object):
    """
    A Urey-Bradley angle type with a set of parameters. It has the same
    functional form as a bond, but it is defined between two atoms forming a
    valence angle separated by two bonds.

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom involved in the Urey-Bradley bond
    atom2 : :class:`Atom`
        The other atom involved in the Urey-Bradley bond
    type : :class:`BondType`
        The Urey-Bradley bond type that defines the parameters for this bond

    Notes
    -----
    You can test whether an :class:`Atom` is contained within the bond using the
    `in` operator. A :class:`BondError` is raised if `atom1` and `atom2` are
    identical.  You can also test that a :class:`Bond` is contained in this
    Urey-Bradley valence angle

    Examples
    --------
    >>> a1, a2, a3 = Atom(), Atom(), Atom()
    >>> b1 = Bond(a1, a2)
    >>> b2 = Bond(a2, a3)
    >>> angle = Angle(a1, a2, a3)
    >>> urey = UreyBradley(a1, a3)
    >>> a1 in urey and a3 in urey
    True
    >>> b1 in urey and b2 in urey
    True
    """
    def __init__(self, atom1, atom2, type=None):
        """ Bond constructor """
        # Make sure we're not bonding me to myself
        if atom1 is atom2:
            raise BondError('Cannot angle atom to itself!')
        # Order the atoms so the lowest atom # is first
        self.atom1 = atom1
        self.atom2 = atom2
        # Log this urey-bradley in the atoms
        atom1.urey_bradleys.append(self)
        atom2.urey_bradleys.append(self)
        # Load the force constant and equilibrium distance
        self.type = type

    @property
    def ub_type(self):
        warnings.warn("ub_type has been replaced by type", DeprecationWarning)
        return self.type

    def __contains__(self, thing):
        " Quick and easy way to see if an Atom or Bond is in this Urey-Bradley "
        if isinstance(thing, Atom):
            return thing is self.atom1 or thing is self.atom2
        # If this is a bond things are a bit more complicated since we don't
        # know the central atom of this Urey-Bradley. We need to make sure that
        # one of the atoms of the bond is either atom1 or atom2 and that the
        # OTHER atom in Bond "thing" has the OTHER atom in this Urey-Bradley in
        # its list of bonded partners.
        if not thing.atom1 in self:
            if not thing.atom2 in self:
                # Neither atom is in this Urey-Bradley...
                return False
            # If we are here, thing.atom2 is in self
            end1 = thing.atom2
            cent = thing.atom1
        else:
            # If we are here, thing.atom1 is in self
            end1 = thing.atom1
            cent = thing.atom2
        # If we are here, end1 and cent are set. Look through the bond
        # partners of cent(er) and see if any of them is in this
        # Urey-Bradley (but ONLY if that atom is not the original end1)
        for atm in cent.bond_partners:
            if atm is end1: continue
            if atm in self:
            # If we got here, we found both atoms in this Urey-Bradley
            # separated by 2 bonds, so "thing" IS in this Urey-Bradley
                return True
        # If we got here, we could not find the other atom connected by 2
        # bonds
        return False

    def delete(self):
        """
        Deletes this Urey-Bradley from the atoms that make it up. This method
        removes this urey-bradley from the `urey_bradleys` list for atom1 and
        atom2.
        """
        _delete_from_list(self.atom1.urey_bradleys, self)
        _delete_from_list(self.atom2.urey_bradleys, self)

        self.atom1 = self.atom2 = self.type = None

    def __repr__(self):
        return '<%s %r--%r; type=%r>' % (type(self).__name__,
                self.atom1, self.atom2, self.type)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Improper(_FourAtomTerm):
    """
    A CHARMM-style improper torsion between 4 atoms. The first atom must be the
    central atom, as shown in the schematic below
      
                      A3
                      |
                      |
             A4 ----- A1 ----- A2

    Parameters
    ----------
    atom1 : :class:`Atom`
        The central atom A1 in the schematic above
    atom2 : :class:`Atom`
        An atom in the improper, A2 in the schematic above
    atom3 : :class:`Atom`
        An atom in the improper, A3 in the schematic above
    atom4 : :class:`Atom`
        An atom in the improper, A4 in the schematic above
    type : :class:`ImproperType`
        The ImproperType object containing the parameters for this improper
        torsion

    Notes
    -----
    An Improper torsion can contain bonds or atoms. A bond is contained if it
    exists between atom 1 and any other atom. Raises :class:`BondError` if any
    of the atoms are duplicates.

    Examples
    --------
    >>> a1, a2, a3, a4 = Atom(), Atom(), Atom(), Atom()
    >>> imp = Improper(a1, a2, a3, a4)
    >>> Bond(a1, a2) in imp and Bond(a1, a3) in imp and Bond(a1, a4) in imp
    True
    >>> Bond(a2, a3) in imp
    False
    """
    def __init__(self, atom1, atom2, atom3, atom4, type=None):
        _FourAtomTerm.__init__(self, atom1, atom2, atom3, atom4)
        # Log these impropers in each atom
        atom1.impropers.append(self)
        atom2.impropers.append(self)
        atom3.impropers.append(self)
        atom4.impropers.append(self)
        # Load the force constant and equilibrium angle
        self.type = type

    @property
    def improp_type(self):
        warnings.warn('improp_type has been replaced by type',
                      DeprecationWarning)
        return self.type

    def __contains__(self, thing):
        """
        Quick and easy way to find out if an Atom or Bond is in this Improper
        """
        if isinstance(thing, Atom):
            return _FourAtomTerm.__contains__(self, thing)
        # Treat it like a Bond
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom1 in thing and self.atom3 in thing) or
                (self.atom1 in thing and self.atom4 in thing))

    def same_atoms(self, thing):
        """
        An improper has the same 4 atoms if atom1 is the same for both and the
        other 3 can be in any order
        """
        if isinstance(thing, Improper):
            # I'm comparing with another Improper here. Central atom must be
            # the same. Others can be in any order
            if self.atom1 is not thing.atom1:
                return False
            # Make a set with the remaining atoms. If they are equal, the
            # impropers are equivalent
            selfset = set([self.atom2, self.atom3, self.atom4])
            otherset = set([thing.atom2, thing.atom3, thing.atom4])
            return selfset == otherset
        thing = list(thing)
        # Here, atoms are expected to index from 0 (Python standard) if we
        # are comparing with a list or tuple
        if len(thing) != 4:
            raise DihedralError('Impropers have 4 atoms, not %s' % len(thing))
        if self.atom1.idx != thing[0]:
            return False
        selfset = set([self.atom2.idx, self.atom3.idx, self.atom4.idx])
        otherset = set([thing[1], thing[2], thing[3]])
        return selfset == otherset

    def delete(self):
        """
        Deletes this Improper from the atoms that make it up. This method
        removes this Improper from the `impropers` list for atom1, atom2, atom3,
        and atom4
        """
        _delete_from_list(self.atom1.impropers, self)
        _delete_from_list(self.atom2.impropers, self)
        _delete_from_list(self.atom3.impropers, self)
        _delete_from_list(self.atom4.impropers, self)

        self.atom1 = self.atom2 = self.atom3 = self.atom4 = self.type = None

    def __repr__(self):
        return '<%s; %r--(%r,%r,%r); type=%r>' % (type(self).__name__,
                self.atom1, self.atom2, self.atom3, self.atom4, self.type)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ImproperType(_ListItem, _ParameterType):
    """
    An improper type with a set of improper torsion parameters

    Parameters
    ----------
    psi_k : ``float``
        Force constant in kcal/mol/radians^2
    psi_eq : ``float``
        Equilibrium torsion angle in Degrees
    list : :class:`TrackedList`
        A list of `ImproperType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : ``int``
        The index of this ImproperType inside its containing list

    Notes
    -----
    Two `ImproperType`s are equal if their `psi_k` and `psi_eq` attributes are
    equal

    Examples
    --------
    >>> it1 = ImproperType(10.0, 180.0)
    >>> it2 = ImproperType(10.0, 180.0)
    >>> it1 is it2
    False
    >>> it1 == it2
    True
    >>> it1.idx # Not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> improper_list = []
    >>> improper_list.append(ImproperType(10.0, 180.0, list=improper_list))
    >>> improper_list.append(ImproperType(10.0, 180.0, list=improper_list))
    >>> improper_list[0].idx
    0
    >>> improper_list[1].idx
    1
    """
    def __init__(self, psi_k, psi_eq, list=None):
        _ParameterType.__init__(self)
        self.psi_k = psi_k
        self.psi_eq = psi_eq
        self.list = list
        self._idx = -1

    def __eq__(self, other):
        return self.psi_k == other.psi_k and self.psi_eq == other.psi_eq

    def __repr__(self):
        return '<%s; k=%.3f, PSIeq=%.3f>' % (type(self).__name__,
                self.psi_k, self.psi_eq)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Cmap(object):
    """
    A coupled-torsion correction map term defined between 5 atoms connected by
    four covalent bonds. This is a coupled-torsion potential in which the
    torsions are consecutive.

    Parameters
    ----------
    atom1 : :class:`Atom`
        An atom on one end of the valence coupled-torsion bonded to atom2
    atom2 : :class:`Atom`
        An atom in the middle of the CMAP bonded to atoms 1 and 3
    atom3 : :class:`Atom`
        An atom in the middle of the CMAP bonded to atoms 2 and 4
    atom4 : :class:`Atom`
        An atom in the middle of the CMAP bonded to atoms 3 and 5
    atom5 : :class:`Atom`
        An atom in the middle of the CMAP bonded to atom 4
    type : :class:`CmapType`
        The CmapType object containing the parameter map for this term

    Notes
    -----
    A CMAP can contain bonds or atoms. A bond is contained if it exists between
    atoms 1 and 2, between atoms 2 and 3, between atoms 3 and 4, or between
    atoms 4 and 5.

    Examples
    --------
    >>> a1, a2, a3, a4, a5 = Atom(), Atom(), Atom(), Atom(), Atom()
    >>> cmap = Cmap(a1, a2, a3, a4, a5)
    >>> Bond(a1, a2) in cmap and Bond(a2, a3) in cmap
    True
    >>> Bond(a1, a3) in cmap
    False
    """
    def __init__(self, atom1, atom2, atom3, atom4, atom5, type=None):
        # Make sure we're not CMAPping me to myself
        atmlist = [atom1, atom2, atom3, atom4, atom5]
        for i in xrange(len(atmlist)):
            for j in xrange(i+1, len(atmlist)):
                if atmlist[i] is atmlist[j]:
                    raise BondError('Cannot cmap atom to itself!')
        # Set up instances
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.atom5 = atom5
        # Log these cmaps in each atom
        atom1.cmaps.append(self)
        atom2.cmaps.append(self)
        atom3.cmaps.append(self)
        atom4.cmaps.append(self)
        atom5.cmaps.append(self)
        # Load the CMAP interpolation table
        self.type = type

    @classmethod
    def extended(cls, atom1, atom2, atom3, atom4,
                 atom5, atom6, atom7, atom8, type=None):
        """
        Alternative constructor for correction maps defined with 8 atoms (each
        torsion being separately specified). Correction maps are, to the best of
        my knowledge, used to parametrized "consecutive" torsions (i.e., atoms
        2, 3, and 4 are common to both torsions). However, the CHARMM definition
        specifies the torsions separately.

        If the torsions are _not_ consecutive (i.e., if atoms 2, 3, and 4 are
        not the same as atoms 5, 6, and 7), NotImplementedError is raised.

        Parameters
        ----------
        atom1 : :class:`Atom`
            The first atom of the first torsion
        atom2 : :class:`Atom`
            The second atom of the first torsion
        atom3 : :class:`Atom`
            The third atom of the first torsion
        atom4 : :class:`Atom`
            The fourth atom of the first torsion
        atom5 : :class:`Atom`
            The first atom of the second torsion
        atom6 : :class:`Atom`
            The second atom of the second torsion
        atom7 : :class:`Atom`
            The third atom of the second torsion
        atom8 : :class:`Atom`
            The fourth atom of the second torsion
        type : :class:`CmapType`
            The CmapType object containing the parameter map for this term
        """
        if atom2 is not atom5 or atom3 is not atom6 or atom7 is not atom4:
            raise NotImplementedError('Only consecutive coupled-torsions are '
                                      'supported by CMAPs currently')
        return cls(atom1, atom2, atom3, atom4, atom8, type=type)

    @property
    def cmap_type(self):
        warnings.warn('cmap_type has been replaced by type', DeprecationWarning)
        return self.type

    def __contains__(self, thing):
        """
        Quick and easy way to find out an atom or bond is in this
        coupled-torsion
        """
        if isinstance(thing, Atom):
            return (thing is self.atom1 or thing is self.atom2 or
                    thing is self.atom3 or thing is self.atom4 or
                    thing is self.atom5)
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing) or
                (self.atom3 in thing and self.atom4 in thing) or
                (self.atom4 in thing and self.atom5 in thing))

    def same_atoms(self, thing):
        """
        A coupled-torsion is equivalent if the 5 atoms are in the same or
        reverse order Allow comparison with another type of cmap or with a
        sequence of 5 indexes
        """
        if isinstance(thing, Cmap):
            # I'm comparing with another Improper here. Central atom must be the
            # same. Others can be in any order
            return ((self.atom1 is thing.atom1 and self.atom2 is thing.atom2 and
                    self.atom3 is thing.atom3 and self.atom4 is thing.atom4 and
                    self.atom5 is thing.atom5) or
                    (self.atom1 is thing.atom5 and self.atom2 is thing.atom4 and
                    self.atom3 is thing.atom3 and self.atom4 is thing.atom2 and
                    self.atom5 is thing.atom1))
        thing = list(thing)
        # Here, atoms are expected to index from 0 (Python standard) if we
        # are comparing with a list or tuple
        if len(thing) != 5:
            raise DihedralError('CMAP can compare to 5 elements, not %d'
                                % (type(thing).__name__, len(thing)))
        return ((self.atom1.idx == thing[0] and self.atom2.idx == thing[1] and
                 self.atom3.idx == thing[2] and self.atom4.idx == thing[3] and
                 self.atom5.idx == thing[4]) or
                (self.atom1.idx == thing[4] and self.atom2.idx == thing[3] and
                 self.atom3.idx == thing[2] and self.atom4.idx == thing[1] and
                 self.atom5.idx == thing[0]))

    def delete(self):
        """
        Deletes this Cmap from the atoms that make it up. This method removes
        the Cmap from the `cmaps` list for atom1, atom2, atom3, atom4, and atom5
        """
        _delete_from_list(self.atom1.cmaps, self)
        _delete_from_list(self.atom2.cmaps, self)
        _delete_from_list(self.atom3.cmaps, self)
        _delete_from_list(self.atom4.cmaps, self)
        _delete_from_list(self.atom5.cmaps, self)

        self.atom1 = self.atom2 = self.atom3 = self.atom4 = self.atom5 = None
        self.type = None

    def __repr__(self):
        return '<%s; %r--%r--%r--%r--%r; type=%r>' % (type(self).__name__,
                self.atom1, self.atom2, self.atom3, self.atom4, self.atom5,
                self.type)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CmapType(_ListItem, _ParameterType):
    """
    A CMAP type with a potential energy interpoloation grid mapping out the 2-D
    potential of coupled torsions.

    Parameters
    ----------
    resolution : ``int``
        The number of grid points in the correction map potential in both
        torsion dimensions
    grid : ``iterable of floats``
        This must be a 1-dimensional list of grid values. The dimension must be
        `resolution*resolution`, and must be row-major (i.e., the second
        dimension changes the fastest) when indexed with 2 indices.
    list : :class:`TrackedList`
        A list of `CmapType`s in which this is a member

    Attributes
    ----------
    resolution : ``int``
        Potential grid resolution (see description in Parameters)
    grid : :class:`_CmapGrid`
        A _CmapGrid object defining the interpolating potential energy grid,
        with each point having the units kcal/mol
    comments : ``list(str)``
        List of strings that represent comments about this parameter type
    list : :class:`TrackedList`
        If not None, this is a list in which this instance _may_ be a member
    idx : ``int``
        The index of this CmapType inside its containing list

    Notes
    -----
    A CmapError is raised if the `grid` does not have the correct number of
    elements for the given `resolution`. Two `CmapType`s are equal if their
    resolution is the same and each grid point is the same to within 1E-8

    See the docs for `_CmapGrid` for information on how to access the
    interpolating grid data.

    Examples
    --------
    >>> ct = CmapType(2, [0, 1, 2, 3])
    >>> ct.grid[0,0], ct.grid[0]
    (0, 0)
    >>> ct.grid[0,1], ct.grid[1]
    (1, 1)
    >>> ct.grid[1,0], ct.grid[2]
    (2, 2)
    >>> ct.grid[1,1], ct.grid[3]
    (3, 3)
    >>> ct.idx # not part of a list or iterable
    -1
    """
    def __init__(self, resolution, grid, comments=None, list=None):
        _ParameterType.__init__(self)
        self.resolution = resolution
        self.grid = _CmapGrid(resolution, grid)
        if len(grid) != self.resolution * self.resolution:
            raise CmapError('CMAP grid does not match expected resolution')
        if comments is None:
            self.comments = []
        else:
            self.comments = comments
        self._idx = -1
        self.list = list

    def __eq__(self, other):
        return (self.resolution == other.resolution and
                all([abs(i - j) < TINY for i, j in zip(self.grid, other.grid)]))

    def __repr__(self):
        return '<%s; res=%d>' % (type(self).__name__, self.resolution)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _CmapGrid(object):
    """
    A grid object for storing Correction map data. Data can be accessed in one
    of two ways; either with 1 or 2 indexes. If 2 indexes [i,j] are given, the
    index into the flattened array is i*resolution+j. Indexing starts from 0.

    The _CmapGrid usually has ranges for the two angles from -180 to 180. Some
    places will expect the range to be 0-360 degrees (e.g., OpenMM). The
    switch_range method returns a _CmapGrid with the "other" range than the
    current object.  i.e., if the current range is from 0 to 360 degrees, the
    _CmapGrid returned by `switch_range` will be from -180 to 180, and vice
    versa.

    Example:
    >>> g = _CmapGrid(2, [0, 1, 2, 3])
    >>> print('%s %s' % (g[0], g[0,0]))
    0 0
    >>> print('%s %s' % (g[1], g[0,1]))
    1 1
    >>> print(g[1,0])
    2
    >>> g[1,1] = 10
    >>> print(g[3])
    10
    >>> print(g.switch_range())
    [10.0000, 2.0000
     1.0000, 0.0000]
    """

    def __init__(self, resolution, data=None):
        self.resolution = resolution
        if data is None:
            self._data = [0 for i in xrange(self.resolution*self.resolution)]
        else:
            self._data = data

    @property
    def transpose(self):
        """ The transpose of the potential grid """
        try:
            return self._transpose
        except AttributeError:
            pass
        _transpose = []
        size = len(self._data)
        for i in xrange(self.resolution):
            piece = [self[j] for j in xrange(i, size, self.resolution)]
            _transpose += piece
        self._transpose = _CmapGrid(self.resolution, _transpose)
        return self._transpose

    T = transpose

    def __getitem__(self, idx):
        if isinstance(idx, tuple):
            if idx[0] >= self.resolution or idx[1] >= self.resolution:
                raise IndexError('_CmapGrid: Index out of range')
            return self._data[self.resolution*idx[0]+idx[1]]
        return self._data[idx]

    def __setitem__(self, idx, val):
        if isinstance(idx, tuple):
            if idx[0] >= self.resolution or idx[1] >= self.resolution:
                raise IndexError('_CmapGrid: Index out of range')
            self._data[self.resolution*idx[0]+idx[1]] = val
        else:
            try:
                indices = xrange(*idx.indices(len(self._data)))
            except AttributeError:
                self._data[idx] = val
            else:
                try:
                    lenval = len(val)
                except TypeError:
                    lenval = 1
                if lenval == 1:
                    for x in indices:
                        self._data[x] = val
                elif lenval != len(indices):
                    raise ValueError('Wrong number of values setting a slice')
                else:
                    for x, y in zip(indices, val):
                        self._data[x] = y

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __repr__(self):
        return '<_CmapGrid: %dx%d>' % (self.resolution, self.resolution)

    def __str__(self):
        retstr = '[%.4f,' % self._data[0]
        fmt = ' %.4f'
        for i, val in enumerate(self):
            if i == 0: continue
            retstr += fmt % val
            if (i+1) % self.resolution == 0 and i != len(self._data) - 1:
                retstr += '\n'
            elif i != len(self) - 1:
                retstr += ','
        return retstr + ']'

    def __eq__(self, other):
        try:
            if self.resolution != other.resolution:
                return False
            for x, y in zip(self, other):
                if abs(x - y) > TINY:
                    return False
            return True
        except AttributeError:
            return TypeError('Bad type comparison with _CmapGrid')

    def switch_range(self):
        """
        Returns a grid object whose range is 0 to 360 degrees in both dimensions
        instead of -180 to 180 degrees (or -180 to 180 degrees if the range is
        already 0 to 360 degrees)
        """
        res = self.resolution
        mid = res // 2
        newgrid = _CmapGrid(res)
        for i in xrange(res):
            ii = (i + mid) % res
            for j in xrange(res):
                jj = (j + mid) % res
                # Start from the middle
                newgrid[i, j] = self[ii, jj]
        return newgrid

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TrigonalAngle(_FourAtomTerm):
    """
    A trigonal-angle term in the AMOEBA force field. It exists in a pattern like
    the one shown below

                                A1
                                |
                                |
                          A4----A2----A3

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom involved in the trigonal angle
    atom2 : :class:`Atom`
        The central atom involved in the trigonal angle
    atom3 : :class:`Atom`
        The third atom involved in the trigonal angle
    atom4 : :class:`Atom`
        The fourth atom involved in the trigonal angle
    type : :class:`AngleType`
        The angle type containing the parameters

    Notes
    -----
    Either `Atom`s or `Bond`s can be contained within this trigonal angle
    """
    def __init__(self, atom1, atom2, atom3, atom4, type=None):
        _FourAtomTerm.__init__(self, atom1, atom2, atom3, atom4)
        self.type = type

    @property
    def trigang_type(self):
        warnings.warn('trigang_type has been replaced with type',
                      DeprecationWarning)
        return type

    def __contains__(self, thing):
        if isinstance(thing, Atom):
            return _FourAtomTerm.__contains(self, thing)
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing) or
                (self.atom2 in thing and self.atom4 in thing))

    def __repr__(self):
        return '<%s; %r--(%r,%r,%r); type=%r>' % (type(self).__name__,
                self.atom2, self.atom1, self.atom3, self.atom4, self.type)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class OutOfPlaneBend(_FourAtomTerm):
    """
    Out-of-plane bending term in the AMOEBA force field. The bond pattern is the
    same as :class:`TrigonalAngle`

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom involved in the trigonal angle
    atom2 : :class:`Atom`
        The central atom involved in the trigonal angle
    atom3 : :class:`Atom`
        The third atom involved in the trigonal angle
    atom4 : :class:`Atom`
        The fourth atom involved in the trigonal angle
    type : :class:`OutOfPlaneBendType`
        The angle type containing the parameters

    Notes
    -----
    Either `Atom`s or `Bond`s can be contained within this trigonal angle
    """
    def __init__(self, atom1, atom2, atom3, atom4, type=None):
        _FourAtomTerm.__init__(self, atom1, atom2, atom3, atom4)
        self.type = type

    @property
    def oopbend_type(self):
        warnings.warn('oopbend_type has been replaced with type',
                      DeprecationWarning)
        return type

    def __contains__(self, thing):
        if isinstance(thing, Atom):
            return _FourAtomTerm.__contains(self, thing)
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing) or
                (self.atom2 in thing and self.atom4 in thing))

    def __repr__(self):
        return '<%s; %r--(%r,%r,%r); type=%r>' % (type(self).__name__,
                self.atom2, self.atom1, self.atom3, self.atom4, self.type)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class OutOfPlaneBendType(_ListItem, _ParameterType):
    """
    An angle type with a set of angle parameters

    Parameters
    ----------
    k : ``float``
        Force constant in kcal/mol/radians^2
    list : :class:`TrackedList`
        A list of `OutOfPlaneBendType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : ``int``
        The index of this OutOfPlaneBendType inside its containing list

    Notes
    -----
    Two `OutOfPlaneBendType`s are equal if their `k` attribute is equal

    Examples
    --------
    >>> ot1 = OutOfPlaneBendType(10.0)
    >>> ot2 = OutOfPlaneBendType(10.0)
    >>> ot1 is ot2
    False
    >>> ot1 == ot2
    True
    >>> ot1.idx # not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> oopbend_list = []
    >>> oopbend_list.append(OutOfPlaneBendType(10.0, list=oopbend_list))
    >>> oopbend_list.append(OutOfPlaneBendType(10.0, list=oopbend_list))
    >>> oopbend_list[0].idx
    0
    >>> oopbend_list[1].idx
    1
    """
    def __init__(self, k, list=None):
        _ParameterType.__init__(self)
        self.k = k
        self._idx = -1
        self.list = list

    def __eq__(self, other):
        return self.k == other.k

    def __repr__(self):
        return '<%s; k=%.3f>' % (type(self).__name__, self.k)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PiTorsion(object):
    r"""
    Defines a pi-torsion term in the AMOEBA force field. The Pi-torsion is
    defined around a sp2-hybridized pi-delocalized orbital (like an amide) by 6
    atoms, as shown in the schematic below.


         A2           A5-AAA
           \         /
            \       /
             A3---A4
            /       \
           /         \
     AAA-A1           A6

    In the above schematic, A3 and A4 are sp2-hybridized, and atoms A2 and A6
    are bonded *only* to A3 and A4, respectively. Atoms A1 and A5 are each
    bonded to 3 other atoms.

    Parameters
    ----------
    atom1 : :class:`Atom`
        atom A1 in the schematic above
    atom2 : :class:`Atom`
        atom A2 in the schematic above
    atom3 : :class:`Atom`
        atom A3 in the schematic above
    atom4 : :class:`Atom`
        atom A4 in the schematic above
    atom5 : :class:`Atom`
        atom A5 in the schematic above
    atom6 : :class:`Atom`
        atom A6 in the schematic above
    type : :class:`DihedralType`
        The parameters for this Pi-torsion

    Notes
    -----
    Both :class:`Bond`s and :class:`Atom`s can be contained in a pi-torsion
    """
    def __init__(self, atom1, atom2, atom3, atom4, atom5, atom6, type=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.atom5 = atom5
        self.atom6 = atom6
        self.type = type

    @property
    def pitor_type(self):
        warnings.warn("pitor_type has been replaced by type",
                      DeprecationWarning)
        return self.type

    def __contains__(self, thing):
        if isinstance(thing, Atom):
            return (thing is self.atom1 or thing is self.atom2 or
                    thing is self.atom3 or thing is self.atom4 or
                    thing is self.atom5 or thing is self.atom6)
        # Assume Bond
        return ((self.atom2 in thing and self.atom3 in thing) or
                (self.atom1 in thing and self.atom3 in thing) or
                (self.atom3 in thing and self.atom4 in thing) or
                (self.atom4 in thing and self.atom5 in thing) or
                (self.atom4 in thing and self.atom6 in thing))

    def delete(self):
        """ Sets all atoms to None and deletes the type """
        self.type = self.atom1 = self.atom2 = self.atom3 = None
        self.atom4 = self.atom5 = self.atom6 = None

    def __repr__(self):
        return '<%s; (%r,%r)--%r--%r--(%r,%r); type=%r>' % (type(self).__name__,
                self.atom1, self.atom2, self.atom3, self.atom4, self.atom5,
                self.atom6, self.type)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class StretchBend(object):
    """
    This term models the stretching and bending of a standard valence angle, and
    is used in the AMOEBA force field

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom on one end of the angle
    atom2 : :class:`Atom`
        The central atom in the angle
    atom3 : :class:`Atom`
        The atom on the other end of the angle
    type : :class:`StretchBendType`
        The type containing the stretch-bend parameters

    Notes
    -----
    Both :class:`Bond`s and :class:`Atom`s can be contained in a stretch-bend term
    """
    def __init__(self, atom1, atom2, atom3, type=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.type = type

    @property
    def strbnd_type(self):
        warnings.warn("strbnd_type has been replaced by type",
                      DeprecationWarning)
        return self.type

    def __contains__(self, thing):
        if isinstance(thing, Atom):
            return (self.atom1 is thing or self.atom2 is thing or
                    self.atom3 is thing)
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing))

    def delete(self):
        """ Sets all of the atoms and parameter type to None """
        self.atom1 = self.atom2 = self.atom3 = self.type = None

    def __repr__(self):
        return '<%s %r--%r--%r; type=%r>' % (type(self).__name__,
                self.atom1, self.atom2, self.atom3, self.type)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class StretchBendType(_ListItem, _ParameterType):
    """
    A stretch-bend type with two distances and an angle in AMOEBA

    Parameters
    ----------
    k1 : ``float``
        First force constant in kcal/mol/(radians*angstroms)
    k2 : ``float``
        Second force constant in kcal/mol/(radians*angstroms)
    req1 : ``float``
        Equilibrium bond distance for bond between the first and second atoms in
        Angstroms
    req2 : ``float``
        Equilibrium bond distance for bond between the second and third atoms in
        Angstroms
    theteq : ``float``
        Equilibrium angle in degrees
    list : :class:`TrackedList`
        A list of `StretchBendType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : ``int``
        The index of this StretchBendType inside its containing list

    Notes
    -----
    Two `StretchBendType`s are equal if their `req1`, `req2`, `theteq`, and `k`
    attributes are equal

    Examples
    --------
    >>> sbt1 = StretchBendType(10.0, 10.0, 1.0, 1.0, 180.0)
    >>> sbt2 = StretchBendType(10.0, 10.0, 1.0, 1.0, 180.0)
    >>> sbt1 is sbt2
    False
    >>> sbt1 == sbt2
    True
    >>> sbt1.idx # Not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> strbnd_list = []
    >>> strbnd_list.append(StretchBendType(10.0, 10.0, 1.0, 1.0, 180.0, strbnd_list))
    >>> strbnd_list.append(StretchBendType(10.0, 10.0, 1.0, 1.0, 180.0, strbnd_list))
    >>> strbnd_list[0].idx
    0
    >>> strbnd_list[1].idx
    1
    """
    def __init__(self, k1, k2, req1, req2, theteq, list=None):
        _ParameterType.__init__(self)
        self.k1 = k1
        self.k2 = k2
        self.req1 = req1
        self.req2 = req2
        self.theteq = theteq
        self._idx = -1
        self.list = list

    def __eq__(self, other):
        return (self.k1 == other.k1 and self.k2 == other.k2 and
                self.req1 == other.req1 and self.req2 == other.req2 and
                self.theteq == other.theteq)

    def __repr__(self):
        return '<%s; Req_1=%.3f, Req_2=%.3f, THETAeq=%.3f, k1=%.3f, k2=%.3f>' \
                % (type(self).__name__, self.req1, self.req2, self.theteq,
                   self.k1, self.k2)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TorsionTorsion(Cmap):
    """
    This is a coupled-torsion map used in the AMOEBA force field similar to the
    correction-map (CMAP) potential used by the CHARMM force field

    Parameters
    ----------
    atom1 : :class:`Atom`
        An atom on one end of the valence torsion-torsion bonded to atom2
    atom2 : :class:`Atom`
        An atom in the middle of the torsion-torsion bonded to atoms 1 and 3
    atom3 : :class:`Atom`
        An atom in the middle of the torsion-torsion bonded to atoms 2 and 4
    atom4 : :class:`Atom`
        An atom in the middle of the torsion-torsion bonded to atoms 3 and 5
    atom5 : :class:`Atom`
        An atom in the middle of the torsion-torsion bonded to atom 4
    type : :class:`TorsionTorsionType`
        The TorsionTorsionType object containing the parameter map for this term

    Notes
    -----
    A TorsionTorsion can contain bonds or atoms. A bond is contained if it
    exists between atoms 1 and 2, between atoms 2 and 3, between atoms 3 and 4,
    or between atoms 4 and 5.

    Examples
    --------
    >>> a1, a2, a3, a4, a5 = Atom(), Atom(), Atom(), Atom(), Atom()
    >>> tortor = TorsionTorsion(a1, a2, a3, a4, a5)
    >>> Bond(a1, a2) in tortor and Bond(a2, a3) in tortor
    True
    >>> Bond(a1, a3) in tortor
    False
    """
    def __init__(self, atom1, atom2, atom3, atom4, atom5, type=None):
        Cmap.__init__(self, atom1, atom2, atom3, atom4, atom5, type)
        atom1.tortors.append(self)
        atom2.tortors.append(self)
        atom3.tortors.append(self)
        atom4.tortors.append(self)
        atom5.tortors.append(self)
        atom1.tortor_to(atom2)
        atom1.tortor_to(atom3)
        atom1.tortor_to(atom4)
        atom1.tortor_to(atom5)
        atom2.tortor_to(atom3)
        atom2.tortor_to(atom4)
        atom2.tortor_to(atom5)
        atom3.tortor_to(atom4)
        atom3.tortor_to(atom5)
        atom4.tortor_to(atom5)

    @property
    def tortor_type(self):
        warnings.warn("tortor_type has been replaced by type",
                      DeprecationWarning)
        return self.type

    def delete(self):
        """
        Deletes this TorsionTorsion from the atoms that make it up. This method
        removes the TorsionTorsion from the `tortors` list for atom1, atom2,
        atom3, atom4, and atom5, and removes each atom from the others'
        tortor_partners list.
        """
        _delete_from_list(self.atom1.tortors, self)
        _delete_from_list(self.atom2.tortors, self)
        _delete_from_list(self.atom3.tortors, self)
        _delete_from_list(self.atom4.tortors, self)
        _delete_from_list(self.atom5.tortors, self)

        _delete_from_list(self.atom1._tortor_partners, self.atom2)
        _delete_from_list(self.atom1._tortor_partners, self.atom3)
        _delete_from_list(self.atom1._tortor_partners, self.atom4)
        _delete_from_list(self.atom1._tortor_partners, self.atom5)
        _delete_from_list(self.atom2._tortor_partners, self.atom1)
        _delete_from_list(self.atom2._tortor_partners, self.atom3)
        _delete_from_list(self.atom2._tortor_partners, self.atom4)
        _delete_from_list(self.atom2._tortor_partners, self.atom5)
        _delete_from_list(self.atom3._tortor_partners, self.atom1)
        _delete_from_list(self.atom3._tortor_partners, self.atom2)
        _delete_from_list(self.atom3._tortor_partners, self.atom4)
        _delete_from_list(self.atom3._tortor_partners, self.atom5)
        _delete_from_list(self.atom4._tortor_partners, self.atom1)
        _delete_from_list(self.atom4._tortor_partners, self.atom2)
        _delete_from_list(self.atom4._tortor_partners, self.atom3)
        _delete_from_list(self.atom4._tortor_partners, self.atom5)
        _delete_from_list(self.atom5._tortor_partners, self.atom1)
        _delete_from_list(self.atom5._tortor_partners, self.atom2)
        _delete_from_list(self.atom5._tortor_partners, self.atom3)
        _delete_from_list(self.atom5._tortor_partners, self.atom4)

        self.atom1 = self.atom2 = self.atom3 = self.atom4 = self.atom5 = None
        self.type = None

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _TorTorTable(object):
    """
    Contains an interpolating potential grid for a coupled-torsion in the AMOEBA
    force field.

    Parameters
    ----------
    ang1 : ``list of floats``
        Angles in the first dimension of the interpolation table
    ang2 : ``list of floats``
        Angles in the second dimension of the interpolation table
    data : ``list of floats``
        Value of the potential grid at each point (ang2 varies fastest)

    Notes
    -----
    Raises `AmoebaError` if the dimension of the data array does not match the
    number of data required by ang1 and ang2

    Elements from the table are obtained and/or set by using angles as indexes.
    If the pair of angles was not one of the original angles passed to this
    table, a KeyError is raised.

    Examples
    --------
    >>> table = _TorTorTable([1.0, 2.0], [1.0, 2.0, 3.0], [1, 2, 3, 4, 5, 6])
    >>> table[1.0,1.0]
    1
    >>> table[1.0,2.0]
    2
    >>> table[2.0,2.0]
    5
    >>> table[2.0,2.0] = 10
    >>> table.data
    [1, 2, 3, 4, 10, 6]
    """
    def __init__(self, ang1, ang2, data):
        if len(data) != len(ang1) * len(ang2):
            raise AmoebaError('Coupled torsion parameter size mismatch. %dx%d '
                              'grid expects %d elements (got %d)' % (len(ang1),
                              len(ang2), len(ang1)*len(ang2), len(data)))
        self.data = data
        self._indexes = dict()
        i = 0
        for a1 in ang1:
            for a2 in ang2:
                self._indexes[(a1, a2)] = i
                i += 1

    def __getitem__(self, idx, idx2=None):
        if idx2 is None:
            return self.data[self._indexes[idx]]
        return self.data[self._indexes[(idx, idx2)]]

    def __setitem__(self, idx, second, third=None):
        if third is not None:
            idx = self._indexes[(idx, second)]
            value = third
        else:
            idx = self._indexes[idx]
            value = second
        self.data[idx] = value

    def __eq__(self, other):
        try:
            for idx in self._indexes.keys():
                if abs(self[idx] - other[idx]) > TINY:
                    return False
        except KeyError:
            return False
        else:
            return True

    def __ne__(self, other):
        return not self.__eq__(other)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TorsionTorsionType(_ListItem, _ParameterType):
    """
    The type containing the parameter maps for the Amoeba torsion-torsion
    potentials. It contains the original potential as well as interpolated first
    and second derivatives for the AMOEBA force field.

    Parameters
    ----------
    dims : ``tuple of 2 ints``
        The table dimensions
    ang1 : ``list of floats``
        The list of angles in the first dimension
    ang2 : ``list of floats``
        The list of angles in the second dimension
    f : ``list of floats``
        The interpolation table for the energy
    dfda1 : ``list of floats``
        The interpolation table of the gradient w.r.t. angle 1
    dfda2 : ``list of floats``
        The interpolation table of the gradient w.r.t. angle 2
    d2fda1da2 : ``list of floats``
        The interpolation table of the 2nd derivative w.r.t. both angles
    list : :class:`TrackedList`
        The list containing this coupled torsion-torsion map

    Attributes
    ----------
    dims : ``tuple of 2 ints``
        The table dimensions
    ang1 : ``list of floats``
        The list of angles in the first dimension
    ang2 : ``list of floats``
        The list of angles in the second dimension
    f : :class:`_TorTorTable`
        The interpolation table for the energy as a _TorTorTable
    dfda1 : :class:`_TorTorTable`
        The interpolation table for the first gradient as a _TorTorTable
    dfda2 : :class:`_TorTorTable`
        The interpolation table for the second gradient as a _TorTorTable
    d2fda1da2 : :class:`_TorTorTable`
        The interpolation table for the second derivative as a _TorTorTable
    list : :class:`TrackedList`
        The list that may, or may not, contain this TorsionTorsionType
    idx : ``int``
        The index of this item in the list or iterable defined by `list`
    """
    def __init__(self, dims, ang1, ang2, f,
                 dfda1=None, dfda2=None, d2fda1da2=None, list=None):
        _ParameterType.__init__(self)
        if len(dims) != 2:
            raise ValueError('dims must be a 2-dimensional iterable')
        if len(ang1) != dims[0] or len(ang2) != dims[1]:
            raise ValueError('dims does match the angle definitions')
        self.dims = tuple(dims)
        self.ang1 = ang1
        self.ang2 = ang2
        self.f = _TorTorTable(ang1, ang2, f)
        if dfda1 is None:
            self.dfda1 = None
        else:
            self.dfda1 = _TorTorTable(ang1, ang2, dfda1)
        if dfda2 is None:
            self.dfda2 = None
        else:
            self.dfda2 = _TorTorTable(ang1, ang2, dfda2)
        if d2fda1da2 is None:
            self.d2fda1da2 = None
        else:
            self.d2fda1da2 = _TorTorTable(ang1, ang2, d2fda1da2)
        self._idx = -1
        self.list = list

    def __eq__(self, other):
        if self.dims != other.dims: return False
        if self.ang1 != other.ang1: return False
        if self.ang2 != other.ang2: return False
        return (self.f == other.f and self.dfda1 == other.dfda1 and
                self.dfda2 == other.dfda2 and self.d2fda1da2 == other.d2fda1da2)

    def __repr__(self):
        return '<%s; %dx%d>' % (type(self).__name__, self.dims[0], self.dims[1])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ChiralFrame(object):
    """
    A chiral frame as defined in the AMOEBA force field. It defines the frame of
    reference for a chiral center

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom defined in the chiral frame
    atom2 : :class:`Atom`
        The second atom defined in the chiral frame
    chirality : ``int``
        Either 1 or -1 to identify directionality. A ValueError is raised if a
        different value is provided

    Notes
    -----
    A chiral frame can only contain atoms.
    """
    def __init__(self, atom1, atom2, chirality):
        self.atom1 = atom1
        self.atom2 = atom2
        if chirality != 1 and chirality != -1:
            raise ValueError('chirality must be 1 or -1')
        self.chirality = chirality

    def __contains__(self, thing):
        return thing is self.atom1 or thing is self.atom2

    def __repr__(self):
        return '<%s; %r--%r, direction=%d>' % (type(self).__name__, self.atom1,
                self.atom2, self.chirality)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class MultipoleFrame(object):
    """
    This defines the frame of reference for computing multipole interactions in
    the AMOEBA force field.

    Parameters
    ----------
    atom : :class:`Atom`
        The atom for which the frame of reference is defined
    frame_pt_num : ``int``
        The frame point number
    vectail : ``int``
        The vector tail index
    vechead : ``int``
        The vector head index
    nvec : ``int``
        The number of vectors

    Examples
    --------
    >>> atom = Atom()
    >>> mf = MultipoleFrame(atom, 0, 1, 2, 3)
    >>> atom in mf
    True
    >>> mf.frame_pt_num
    0
    >>> mf.vectail
    1
    >>> mf.vechead
    2
    >>> mf.nvec
    3
    """
    def __init__(self, atom, frame_pt_num, vectail, vechead, nvec):
        self.atom = atom
        self.frame_pt_num = frame_pt_num
        self.vectail = vectail
        self.vechead = vechead
        self.nvec = nvec

    def __contains__(self, thing):
        return self.atom is thing

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Residue(_ListItem):
    """
    A single residue that is composed of a small number of atoms

    Parameters
    ----------
    name : ``str``
        Name of the residue. Typical convention is to choose a name that is 4
        characters or shorter
    number : ``int``
        Residue number assigned in the input structure
    chain : ``str``
        The 1-letter chain identifier for this residue
    insertion_code : ``str``
        The insertion code (used in PDB files) for this residue
    list : :class:`TrackedList`
        List of residues in which this residue is a member

    Attributes
    ----------
    name : ``str``
        The name of this residue
    number : ``int``
        The number of this residue in the input structure
    idx : ``int``
        The index of this residue inside the container. If this residue has no
        container, or it is not present in the container, idx is -1
    chain : ``str``
        The 1-letter chain identifier for this residue
    insertion_code : ``str``
        The insertion code (used in PDB files) for this residue
    ter : ``bool``
        If True, there is a TER card directly after this residue (i.e., a
        molecule or chain ends). By default, it is False
    list : :class:`TrackedList`
        The container that _may_ have this residue contained inside
    atoms : ``list of`` :`class`Atom` ``instances``
        This is the list of `Atom`s that make up this residue

    Notes
    -----
    - Iterating over a residue will iterate over the atoms. It is exactly
      equivalent to iterating over the `atoms` attribute
    - Supports testing if an Atom instance is contained `in` this residue
    - `len()` returns the number of atoms in this residue
    """

    def __init__(self, name, number=-1, chain='', insertion_code='', list=None):
        self.name = name.strip()
        self.number = number
        self.chain = chain.strip()
        self.insertion_code = insertion_code.strip()
        self.list = list
        self._idx = -1
        self.atoms = []
        self.ter = False

    def add_atom(self, atom):
        """ Adds an atom to this residue

        Parameters
        ----------
        atom : :class:`Atom`
            The atom to add to this residue
        
        Notes
        -----
        This action assigns the `residue` attribute to `atom`
        """
        atom.residue = self
        self.atoms.append(atom)

    def delete_atom(self, atom):
        """
        If an atom is present in this residue, delete it from the list of
        atoms. No change if an atom is not present in this residue.
        """
        for a in self.atoms:
            if atom.residue is self:
                atom.residue = None
                self.atoms = [a for a in self.atoms if a is not atom]

    # Implement some container methods over the list of atoms
    def __contains__(self, thing):
        """ True if an atom is present in this residue """
        return thing in self.atoms

    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)

    def is_empty(self):
        """
        Determines if there are any atoms in this residue

        Returns
        -------
        empty: ``bool``
            ``True`` if there are no atoms left. ``False`` otherwise.
        """
        return len(self) == 0

    def sort(self):
        """ Sorts the atoms in this list by atom index """
        self.atoms.sort()

    def __getitem__(self, idx):
        return self.atoms.__getitem__(idx)

    # Sort by atom indices
    def __lt__(self, other):
        return self.atoms[0].idx < other.atoms[0].idx
    def __gt__(self, other):
        return self.atoms[0].idx > other.atoms[0].idx
    def __le__(self, other):
        return not self.atoms[0].idx > other.atoms[0].idx
    def __ge__(self, other):
        return not self.atoms[0].idx < other.atoms[0].idx

    def __repr__(self):
        if self.number == -1:
            num = self.idx
        else:
            num = self.number
        rep = '<%s %s[%d]' % (type(self).__name__, self.name, num)
        if self.chain:
            rep += '; chain=%s' % self.chain
        if self.insertion_code:
            rep += '; insertion_code=%s' % self.insertion_code
        return rep + '>'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _changes(func):
    """ Decorator to indicate the list has changed """
    def new_func(self, *args, **kwargs):
        self.changed = True
        self.needs_indexing = True
        return func(self, *args, **kwargs)
    return new_func

class TrackedList(list):
    """
    This creates a list type that allows you to see if anything has changed

    Attributes
    ----------
    changed : ``bool``
        Determines if something has been done to fundamentally change the
        underlying topology defined by this list such that the topology needs to
        be rebuilt
    needs_indexing : ``bool``
        A flag to determine whether or not the items in a tracked list need to
        be indexed or not.

    Examples
    --------
    >>> tl = TrackedList()
    >>> tl.append(Atom())
    >>> tl.append(Atom())
    >>> tl.append(Atom())
    >>> tl.needs_indexing, tl.changed
    (True, True)
    >>> tl.index_members()
    >>> tl.needs_indexing, tl.changed
    (False, True)
    >>> tl.changed = False # Must do when changes have been incorporated
    >>> tl.needs_indexing, tl.changed
    (False, False)
    """
    def __init__(self, arg=[]):
        self.changed = False
        self.needs_indexing = False
        return list.__init__(self, arg)

    @_changes
    def __delitem__(self, item):
        """ Deletes items and slices. Make sure all items """
        try:
            indices = xrange(*item.indices(len(self)))
        except AttributeError:
            indices = [item]

        for index in indices:
            try:
                self[index]._idx = -1
            except AttributeError:
                pass
            try:
                self[index].list = None
            except AttributeError:
                pass

        return list.__delitem__(self, item)

    @_changes
    def __delslice__(self, start, stop):
        """ Python 2 still uses __delslice__... """
        if not self: return
        indices = xrange(start, min(stop, len(self)))
        for index in indices:
            try:
                self[index]._idx = -1
            except AttributeError:
                pass
            try:
                self[index].list = None
            except AttributeError:
                pass
        return list.__delslice__(self, start, stop)

    @_changes
    def pop(self, idx=-1):
        item = list.pop(self, idx)
        try:
            item._idx = -1
        except AttributeError:
            # Must be an immutable type, so don't complain
            pass
        return item

    append = _changes(list.append)
    extend = _changes(list.extend)
    insert = _changes(list.insert)
    __setitem__ = _changes(list.__setitem__)
    __iadd__ = _changes(list.__iadd__)
    __imul__ = _changes(list.__imul__)

    # Type-safe methods that return another instance
    def __add__(self, other):
        return TrackedList(list.__add__(self, other))

    def __mul__(self, fac):
        return TrackedList(list.__mul__(self, fac))

    def __getitem__(self, thing):
        retval = list.__getitem__(self, thing)
        if hasattr(thing, 'indices'):
            return TrackedList(retval)
        return retval

    def __getslice__(self, start, end):
        return TrackedList(list.__getslice__(self, start, end))

    def index_members(self):
        """
        Assigns the idx variable for every member of this list to its place in
        the list, if the members of this list permit
        """
        for i, item in enumerate(self):
            try:
                item._idx = i
            except AttributeError:
                # Must be some kind of immutable type... don't worry
                pass
        self.needs_indexing = False

    def claim(self):
        """
        This method causes this list to "claim" all of the items it contains and
        subsequently indexes all of its items.
        """
        for i, item in enumerate(self):
            try:
                item.list = self
            except AttributeError:
                # Must be some kind of immutable type... don't worry
                pass
        self.index_members()

    def prune_unused(self):
        """
        This method inspects the `used` attribute of all of its members, if it
        has one, and deletes any item in which it is set to `False`
        """
        for i in reversed(xrange(len(self))):
            try:
                if not self[i].used:
                    del self[i]
            except AttributeError:
                # Don't worry if we don't have a `used` attribute
                pass

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ResidueList(TrackedList):
    """ Array of `Residue` instances """

    def add_atom(self, atom, resname, resnum, chain='', inscode=''):
        """
        Adds a new atom to the ResidueList, adding a new residue to this list if
        it has a different name or number as the last residue

        Parameters
        ----------
        atom : :class:`Atom`
            The atom to add to this residue list
        resname : ``str``
            The name of the residue this atom belongs to
        resnum : ``int``
            The number of the residue this atom belongs to
        chain : ``str``
            The chain ID character for this residue
        inscode : ``str``
            The insertion code ID character for this residue (it is stripped)

        Notes
        -----
        If the residue name and number differ from the last residue in this
        list, a new residue is added and the atom is added to that residue
        """
        inscode = inscode.strip()
        try:
            last = self[-1]
        except IndexError:
            # Empty list -- add our first residue
            new_res = Residue(resname, resnum, chain, inscode, list=self)
            new_res.add_atom(atom)
            self.append(new_res)
        else:
            if (last.number != resnum or last.name != resname.strip() or
                last.chain != chain.strip() or
                last.insertion_code != inscode.strip()):
                new_res = Residue(resname, resnum, chain, inscode, list=self)
                new_res.add_atom(atom)
                self.append(new_res)
            else:
                last.add_atom(atom)

    def prune(self):
        """
        This function goes through the residue list and removes all empty
        residues from the list. This isn't done automatically when atoms are
        deleted, since it will become very slow. You must remember to do this to
        avoid including empty residues
        """
        # Delete from the back to avoid indexes changing as we iterate
        for i in reversed(xrange(len(self))):
            res = self[i]
            if res.is_empty(): del self[i]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AtomList(TrackedList):
    """
    Array of Atoms

    Notes
    -----
    Deleting an atom from the AtomList also deletes that atom from the residue
    it belongs to.
    """

    @_changes
    def __delitem__(self, idx):
        """ Deleting an atom also needs to delete it from the residue """
        try:
            indices = xrange(*idx.indices(len(self)))
        except AttributeError:
            indices = [idx]

        for index in indices:
            atom = self[index]
            atom._idx = -1
            atom.list = None
            # Make sure we delete this atom from its respective residue
            if atom.residue is not None: atom.residue.delete_atom(atom)

        list.__delitem__(self, idx)

    @_changes
    def __delslice__(self, start, stop):
        """ Python 2 still uses __delslice__... sigh. """
        indices = xrange(start, min(stop, len(self)))
        for index in indices:
            atom = self[index]
            atom._idx = -1
            atom.list = None
            if atom.residue is not None: atom.residue.delete_atom(atom)
        list.__delslice__(self, start, stop)

    @_changes
    def pop(self, idx=-1):
        atom = list.pop(self, idx)
        atom._idx = -1
        atom.list = None
        if atom.residue is not None: atom.residue.delete_atom(atom)
        return atom

    def unmark(self):
        """ Unmark all atoms in this list """
        for atm in self: atm.marked = 0

    @_changes
    def append(self, item):
        """
        Add an Atom to the end of the list and have this list claim ownership of
        the item.

        Parameters
        ----------
        item : :class:`Atom`
            The atom to add to this list

        Notes
        -----
        Only Atom objects should be added here, so if `item` does not have a
        `list` attribute, an AttributeError will be raised. This action assigns
        this list as the `list` attribute to the passed Atom.
        """
        item.list = self
        return list.append(self, item)

    @_changes
    def extend(self, items):
        """
        Add an iterable of `Atom`s to the end of the list and have this list
        claim ownership of the items.

        Parameters
        ----------
        items : iterable of :class:`Atom`s
            The iterable containing Atom instances to add to the end of this
            list

        Notes
        -----
        Any generator passed here will be exhausted. The `list` attribute for
        every object in `items` will have their `list` attribute set to this
        list, and an AttributeError will be raised if this is not possible.
        """
        for item in items:
            item.list = self
        return list.extend(self, items)

    @_changes
    def insert(self, idx, item):
        """
        Insert an Atom into the atom list

        Parameters
        ----------
        idx : ``int``
            The index in front of (i.e., before) which to insert the item
        item : :class:`Atom`
            The atom to insert in the desired index. This atom will be claimed
            by the AtomList
        """
        item.list = self
        return list.insert(self, idx, item)

    def assign_nbidx_from_types(self):
        """
        Assigns the nb_idx attribute of every atom inside here from the
        atom_type definition. If the atom_type is not assigned, RuntimeError is
        raised.

        Returns
        -------
        ``list of dict``
            Each element is a `set` of the `nb_idx` indices for which NBFIX
            alterations are defined for the type with that given that index
            (minus 1, to adjust for indexing from 0 and nb_idx starting from 1)
        """
        idx = 1
        natoms = len(self)
        atom_type_lookups = dict() # For fast lookups
        atom_type_list = []
        try:
            for i, atom in enumerate(self):
                atom.atom_type._idx = -1
        except AttributeError:
            raise RuntimeError('atom types are not assigned')

        for i, atom in enumerate(self):
            type1 = atom.atom_type
            # Skip atom types that have already been assigned
            if type1._idx != -1: continue
            type1._idx = idx
            atom_type_lookups[str(type1)] = type1
            atom_type_list.append(type1)
            for j in xrange(i+1, natoms):
                atom2 = self[j]
                type2 = atom2.atom_type
                # Skip atom types that have already been assigned
                if type2._idx != -1: continue
                if type1 == type2:
                    type2._idx = idx
            idx += 1
        # Now go back through and assign nb_idx from type._idx
        for atom in self:
            atom.nb_idx = atom.atom_type._idx
        # Now collect the nbfixes
        nbfix_list = [set() for i in xrange(idx-1)]
        # Look through all of the atom types and add the nbfixes
        for i, type in enumerate(atom_type_list):
            for key in type.nbfix:
                rmin, eps, rmin14, eps14 = type.nbfix[key]
                try:
                    otype = atom_type_lookups[key]
                except KeyError:
                    continue
                else:
                    obj = (otype._idx, rmin, eps, rmin14, eps14)
                    nbfix_list[i].add(obj)

        return nbfix_list

    def find_original_index(self, idx):
        """
        Finds an atom with the given original index. Cannot assume that the
        original indexes are in order, since reordering may have been necessary.
        As a result, the complexity of this algorithm is O(N)

        Parameters
        ----------
        idx : int
            The integer corresponding to the original index

        Returns
        -------
        atom : :class:`Atom`
            The atom with the original index ``idx``

        Raises
        ------
        IndexError
            If no atom has the original index ``idx``
        """
        for atom in self:
            if atom.number == idx: return atom
        raise IndexError('No atom found with index %d' % idx)

    def __iadd__(self, other):
        return NotImplemented

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class NonbondedException(object):
    """
    The AMOEBA force field has complex exclusion and exception rules (referred
    to as "adjustments" in the Amber-converted files). This class stores
    per-particle exceptions.

    Parameters
    ----------
    atom1 : :class:`Atom`
        One of the atoms in the exclusion pair
    atom2 : :class:`Atom`
        The other atom in the exclusion pair
    type : :class:`NonbondedExceptionType`
        The nonbonded exception type that describes how the various nonbonded
        interactions between these two atoms should work

    Notes
    -----
    NonbondedException objects "contain" two atoms and will return True when
    used with the binary `in` operator
    """
    def __init__(self, atom1, atom2, type=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.type = type

    def __contains__(self, thing):
        return thing is self.atom1 or thing is self.atom2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class NonbondedExceptionType(_ListItem):
    """
    A parameter describing how the various nonbonded interactions between a
    particular pair of atoms is scaled in the AMOEBA force field

    Parameters
    ----------
    vdw_weight : ``float``
        The scaling factor by which van der Waals interactions are multiplied
    multipole_weight : ``float``
        The scaling factor by which multipole interactions are multiplied
    direct_weight : ``float``
        The scaling factor by which direct-space interactions are multiplied
    polar_weight : ``float``
        The scaling factor by which polarization interactions are multiplied
    mutual_weight : ``float``
        The scaling factor by which mutual interactions are multiplied
    list : :class:`TrackedList`
        The list containing this nonbonded exception

    Other Attributes
    ----------------
    idx : ``int``
        The index of this term in the list that contains it
    """
    def __init__(self, vdw_weight, multipole_weight, direct_weight,
                 polar_weight, mutual_weight, list=None):
        self.vdw_weight = vdw_weight
        self.multipole_weight = multipole_weight
        self.direct_weight = direct_weight
        self.polar_weight = polar_weight
        self.mutual_weight = mutual_weight
        self._idx = -1
        self.list = None

    def __eq__(self, other):
        return (self.vdw_weight == other.vdw_weight and
                self.multipole_weight == other.multipole_weight and
                self.direct_weight == other.direct_weight and
                self.polar_weight == other.polar_weight and
                self.mutual_weight == other.mutual_weight)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AtomType(object):
    """
    Atom types can either be compared by indexes or names. Can be assigned with
    a string, integer, (string is not automatically cast) or with both.

    Parameters
    ----------
    name : ``str``
        The name of the atom type
    number : ``int``
        The serial index of the atom type
    mass : ``float``
        The mass of the atom type
    atomic_number : ``int``
        The atomic number of the element of the atom type

    Other Attributes
    ----------------
    epsilon : ``float``
        If set, it is the Lennard-Jones well depth of this atom type
    rmin : ``float``
        If set, it is the Lennard-Jones Rmin/2 parameter of this atom type
    epsilon_14 : ``float``
        If set, it is the Lennard-Jones well depth of this atom type in 1-4
        nonbonded interactions
    rmin_14 : ``float``
        If set, it is the Lennard-Jones Rmin/2 parameter of this atom type in
        1-4 nonbonded interactions
    nbfix : ``dict(str:tuple)``
        A hash that maps atom type names of other atom types with which _this_
        atom type has a defined NBFIX with a tuple containing the terms
        (Rmin, epsilon, Rmin14, Epsilon14)

    Notes
    -----
    This object is primarily used to build parameter databases from parameter
    files

    Examples
    --------
    >>> at = AtomType('HA', 1, 1.008, 1)
    >>> at.name, at.number
    ('HA', 1)
    >>> at2 = AtomType('CA', 2, 12.01, 6)
    >>> at2.name, at2.number
    ('CA', 2)
    >>> print("%s: %d" % (str(at), int(at)))
    HA: 1
    """

    def __init__(self, name, number, mass, atomic_number):
        if number is None and name is not None:
            # If we were given an integer, assign it to number. Otherwise,
            # assign it to the name
            if isinstance(name, int):
                self.number = name
                self.name = None
            else:
                self.name = name
                self.number = None
        else:
            self.name = name
            self.number = int(number)
        self.mass = mass
        self.atomic_number = atomic_number
        # We have no LJ parameters as of yet
        self.epsilon = self.rmin = self.epsilon_14 = self.rmin_14 = None
        # Store each NBFIX term as a dict with the atom type string matching to
        # a 2-element tuple that is rmin, epsilon
        self.nbfix = dict()
        self._idx = -1 # needed for some internal bookkeeping

    def __eq__(self, other):
        """
        Compares based on available properties (name and number, just name,
        or just number)
        """
        if isinstance(other, AtomType):
            if self.name != other.name or self.number != other.number:
                return False
            # Now check if LJ parameters are defined, and make sure those are
            # also equal
            return (abs(self.epsilon - other.epsilon) < TINY and
                    abs(self.rmin - other.rmin) < TINY and
                    abs(self.epsilon_14 - other.epsilon_14) > TINY and
                    abs(self.rmin_14 - other.rmin_14) > TINY and
                    self.nbfix == other.nbfix)
        if isinstance(other, basestring):
            return self.name == other
        if isinstance(other, int):
            return self.number == other
        return other == (self.number, self.name)

    def set_lj_params(self, eps, rmin, eps14=None, rmin14=None):
        """ Sets Lennard-Jones parameters on this atom type """
        if eps14 is None:
            eps14 = eps
        if rmin14 is None:
            rmin14 = rmin
        self.epsilon = eps
        self.rmin = rmin
        self.epsilon_14 = eps14
        self.rmin_14 = rmin14

    def __int__(self):
        """ The integer representation of an AtomType is its index """
        return self.number

    def add_nbfix(self, typename, rmin, epsilon, rmin14, epsilon14):
        """ Adds a new NBFIX exclusion for this atom """
        if rmin14 is None: rmin14 = rmin
        if epsilon14 is None: epsilon14 = epsilon
        self.nbfix[typename] = (rmin, epsilon, rmin14, epsilon14)

    def __str__(self):
        return self.name

    # Comparisons are all based on number
    def __gt__(self, other):
        return self._member_number > other._member_number
    def __lt__(self, other):
        return self._member_number < other._member_number
    def __ge__(self, other):
        return self._member_number > other._member_number or self == other
    def __le__(self, other):
        return self._member_number < other._member_number or self == other

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _UnassignedAtomType(object):
    """
    This raises the appropriate exceptions (MissingParameter) when you try to
    access its properties
    """

    def __int__(self):
        raise MissingParameter('Atom type is not defined')

    def __str__(self):
        raise MissingParameter('Atom type is not defined')

_UnassignedAtomType = _UnassignedAtomType() # Make it a singleton

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AcceptorDonor(object):
    """
    Just a holder for donors and acceptors in CHARMM speak
    
    Parameters
    ----------
    atom1 : :class:`Atom`
        First atom in the donor/acceptor group
    atom2 : :class:`Atom`
        Second atom in the donor/acceptor group
    """
    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2

    def __repr__(self):
        return '<AcceptorDonor; %r %r>' % (self.atom1, self.atom2)

    def __contains__(self, thing):
        """ See if the atom is in this donor/acceptor """
        return thing is self.atom1 or thing is self.atom2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Group(object):
    """
    An 'interacting' group defined by CHARMM PSF files

    Parameters
    ----------
    bs : ``int``
        Not sure
    type : ``int``
        The group type (??)
    move : ``int``
        If the group moves (??)

    Disclaimer
    ----------
    I really don't know what these numbers mean. I'm speculating based on the
    source code of 'chamber', and this section is simply ignored there.
    """
    def __init__(self, bs, type, move):
        self.bs = bs
        self.type = type
        self.move = move

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def Element(mass):
    """ Determines what element the given atom is based on its mass """

    diff = mass
    best_guess = 'EP'

    for element in _Element:
        if abs(Mass[element] - mass) < diff:
            best_guess = element
            diff = abs(Mass[element] - mass)

    return best_guess

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NoUreyBradley = BondType(0.0, 0.0) # singleton representing lack of a U-B term

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == '__main__':
    import doctest
    doctest.testmod()
