"""
This module contains objects that deal with system topology and parameters, such
as atoms, residues, bonds, angles, etc.

by Jason Swails
"""
from __future__ import division

from chemistry.exceptions import (BondError, DihedralError, CmapError,
                                  AmoebaError)
from chemistry.amber.constants import TINY
from chemistry.periodic_table import Mass, Element as _Element
from compat24 import all, property
import warnings
try:
    from itertools import izip as zip
except ImportError:
    pass # Must be Python 3... zip _is_ izip

__all__ = ['Angle', 'AngleType', 'Atom', 'AtomList', 'Bond', 'BondType', 'Cmap',
           'CmapType', 'Dihedral', 'DihedralType', 'Improper', 'ImproperType',
           'OutOfPlaneBend', 'PiTorsion', 'Residue', 'ResidueList',
           'StretchBend', 'StretchBendType', 'TorsionTorsion',
           'TorsionTorsionType', 'TrigonalAngle', 'TrackedList', 'UreyBradley']

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Private classes and methods

class _ListItem(object):
    """
    Helpful methods for items that appear in some kind of tracked list.
    Subclasses of this method interact with their `list` attribute, if they have
    one, in order to determine _where_ in the list they occur.

    Attributes
    ----------
    idx : int
        This is intended to be a read-only variable that determines where in the
        list this particular object is. If there is no `list` attribute for this
        object, or the item is not in the list at all, `idx` is -1
        
    Notes
    -----
    For lists that support indexing its members and tracking when that list is
    changed (so we know when to update indexes), this is a fast lookup, as
    indexing only needs to be done once, at most, for the entire list (until
    changes are made that could change the indexing).

    The `idx` lookup in a TrackedList is therefore an O(1) operation, while it
    is O(N^2) for standard containers.
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
    atom1 : Atom
        The first atom in the term
    atom2 : Atom
        The second atom in the term
    atom3 : Atom
        The third atom in the term
    atom4 : Atom
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

def _delete_from_list(list, item):
    """
    Deletes a requested item from a list. If the item does not exist in the
    list, a ValueError is raised

    Parameters
    ----------
    list : list
        The list from which an item will be deleted
    item : object
        The object to delete from the list
    """
    list.pop(list.index(item))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Atom(_ListItem):
    """ 
    An atom. Only use these as elements in AtomList instances, since AtomList
    will keep track of when indexes and other stuff needs to be updated.

    Parameters
    ----------
    list : AtomList=None
        The AtomList that this atom belongs to. If None, this atom does not
        belong to any list. This can be any iterable, but should typically be an
        AtomList or None
    atomic_number : int=0
        The atomic number of this atom
    name : str=''
        The name of this atom
    type : str=''
        The type name of this atom
    charge : float=0.0
        The partial atomic charge of this atom in fractions of an electron
    mass : float=0.0
        The atomic mass of this atom in daltons
    nb_idx : int=0
        The nonbonded index. This is a pointer that is relevant in the context
        of an Amber topology file and identifies its Lennard-Jones atom type
    radii : float=0.0
        The intrinsic solvation radius of this atom.
    screen : float=0.0
        The Generalized Born screening factor for this atom.
    tree : str='BLA'
        The tree chain identifier assigned to this atom. Relevant in the context
        of an Amber topology file, and not used for very much.

    Attributes
    ----------
    atomic_number : int
        The atomic number of this atom
    list : AtomList (or other iterable)
        The iterable that (possibly) contains this atom
    element : int
        This is an alias for atomic_number
    name : str
        The name of this atom
    type : str
        The name of the atom type assigned to this atom
    mass : float
        The mass, in daltons, of this atom
    charge : float
        The charge, in fractions of an electron, of this atom
    nb_idx : int
        The nonbonded Lennard-Jones index. Required when it is part of an Amber
        topology file instance
    tree : str
        The tree chain classification string. Applies to the Amber topology file
        instance, but is not used for much.
    radii : float
        The intrinsic solvation radius of the atom
    screen : float
        The GB screening factor of the atom
    idx : int
        The index of this atom in the list. Set to -1 if this atom is not part
        of a list or the index cannot otherwise be determined (i.e., if the
        containing list does not support indexing its members)
    residue : Residue
        The Residue that this atom belongs to. This is assigned when this atom
        is passed to `Residue.add_atom` -- see below for more information. Until
        it is set there, it is None
    bonds : list of Bond instances
        list of Bond objects in which this atom is a member
    angles : list of Angle instances
        list of Angle objects in which this atom is a member
    dihedrals : list of Dihedral instances
        list of Dihedral objects in which this atom is a member
    urey_bradleys : list of UreyBradley instances
        list of UreyBradley objects in which this atom is a member (CHARMM,
        AMOEBA)
    impropers : list of Improper instances
        list of Improper objects in which the atom is a member (CHARMM)
    cmaps : list of Cmap instances
        list of Cmap objects in which the atom is a member (CHARMM, AMOEBA)
    tortors : list of TorsionTorsion instances
        list of TorsionTorsion objects in which the atom is a member (AMOEBA)
    bond_partners : list of Atom instances
        list of Atoms to which this atom is bonded.
    angle_partners : list of Atom instances
        list of Atoms to which this atom forms an angle, but not a bond
    dihedral_partners : list of Atom instances
        list of Atoms to which this atom forms an dihedral, but not a bond or
        angle
    tortor_partners : list of Atom instances
        list of Atoms to which this atom forms a coupled Torsion-Torsion, but
        not a bond or angle (AMOEBA)
    exclusion_partners : list of Atom instances
        list of Atoms with which this atom is excluded, but not bonded, angled,
        or dihedraled to
    marked : int
        Mainly for internal use, it is used to indicate when certain atoms have
        been "marked" when traversing the bond network identifying topological
        features (like molecules and rings)

    Possible Attributes
    -------------------
    xx : float
        The X-component of the position of this atom. Only present if
        coordinates have been loaded. Otherwise raises AttributeError
    xy : float
        The Y-component of the position of this atom. Only present if
        coordinates have been loaded. Otherwise raises AttributeError
    xz : float
        The Z-component of the position of this atom. Only present if
        coordinates have been loaded. Otherwise raises AttributeError
    vx : float
        The X-component of the velocity of this atom. Only present if the
        velocities have been loaded. Otherwise raises AttributeError
    vy : float
        The Y-component of the velocity of this atom. Only present if the
        velocities have been loaded. Otherwise raises AttributeError
    vz : float
        The Z-component of the velocity of this atom. Only present if the
        velocities have been loaded. Otherwise raises AttributeError

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
                 tree='BLA'):
        self.list = list
        self._idx = -1
        self.atomic_number = atomic_number
        self.name = name
        self.type = type
        self.charge = charge
        self.mass = mass
        self.nb_idx = nb_idx
        self.radii = radii
        self.screen = screen
        self.tree = tree
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
        return sorted(list(bp))

    @property
    def angle_partners(self):
        """ List of all angle partners that are NOT bond partners """
        bp = set(self._bond_partners)
        ap = set(self._angle_partners)
        return sorted(list(ap - bp))

    @property
    def dihedral_partners(self):
        " List of all dihedral partners that are NOT angle or bond partners "
        dp = set(self._dihedral_partners)
        ap = set(self._angle_partners)
        bp = set(self._bond_partners)
        return sorted(list(dp - ap - bp))

    @property
    def tortor_partners(self):
        """
        List of all 1-5 partners that are NOT in angle or bond partners. This is
        _only_ used in the Amoeba force field
        """
        tp = set(self._tortor_partners)
        dp = set(self._dihedral_partners)
        ap = set(self._angle_partners)
        bp = set(self._bond_partners)
        return sorted(list(tp - dp - ap - bp))

    @property
    def exclusion_partners(self):
        """
        List of all exclusions not otherwise excluded by bonds/angles/torsions
        """
        ep = set(self._exclusion_partners)
        tp = set(self._tortor_partners)
        dp = set(self._dihedral_partners)
        ap = set(self._angle_partners)
        bp = set(self._bond_partners)
        return sorted(list(ep - tp - dp - ap - bp))

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
        other : Atom
            An atom that will be added to `bond_partners`

        Notes
        -----
        This action adds `self` to `other.bond_partners`. Raises `BondError` if
        `other is self`
        """
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
        other : Atom
            An atom that will be added to `angle_partners`

        Notes
        -----
        This action adds `self` to `other.angle_partners`. Raises `BondError` if
        `other is self`
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
        other : Atom
            An atom that will be added to `dihedral_partners`

        Notes
        -----
        This action adds `self` to `other.dihedral_partners`. Raises `BondError`
        if `other is self`
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
        other : Atom
            An atom that will be added to `tortor_partners`

        Notes
        -----
        This action adds `self` to `other.tortor_partners`. Raises `BondError`
        if `other is self`
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
        other : Atom
            An atom that will be added to `exclusion_partners`.

        Notes
        -----
        This action adds `self` to `other.exclusion_partners`. Raises
        `BondError` if `other is self`
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
        keep_exclusions : bool=True
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
        return "<Atom %s [%d]; In %s %d>" % (self.atname, self.starting_index+1,
                self.residue.resname, self.residue.idx)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Bond(object):
    """
    A covalent bond connecting two atoms.

    Parameters (and Attributes)
    ---------------------------
    atom1 : Atom
        The first atom involved in the bond
    atom2 : Atom
        The other atom involved in the bond
    type : BondType=None
        The bond type that defines the parameters for this bond

    Notes
    -----
    You can test whether an Atom is contained within the bond using the `in`
    operator. A `BondError` is raised if `atom1` and `atom2` are identical. This
    bond instance is `append`ed to the `bonds` list for both `atom1` and `atom2`
    and is automatically removed from those lists upon garbage collection

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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondType(_ListItem):
    """
    A bond type with a set of bond parameters

    Parameters (and Attributes)
    ---------------------------
    k : float
        Force constant in kcal/mol/Angstrom^2
    req : float
        Equilibrium bond distance in Angstroms
    list : TrackedList=None
        A list of `BondType`s in which this is a member

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
        self.k = k
        self.req = req
        self.list = list
        self._idx = -1

    def __eq__(self, other):
        return self.k == other.k and self.req == other.req

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Angle(object):
    """
    A valence angle between 3 atoms separated by two covalent bonds.

    Parameters (and Attributes)
    ---------------------------
    atom1 : Atom
        An atom one end of the valence angle
    atom2 : Atom
        The atom in the middle of the valence angle bonded to both other atoms
    atom3 : Atom
        An atom on the other end of the valence angle to atom1
    type : AngleType=None
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleType(_ListItem):
    """
    An angle type with a set of angle parameters

    Parameters (and Attributes)
    ---------------------------
    k : float
        Force constant in kcal/mol/radians^2
    theteq : float
        Equilibrium angle in Degrees
    list : TrackedList=None
        A list of `AngleType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : int
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
        self.k = k
        self.theteq = theteq
        self._idx = -1
        self.list = list

    def __eq__(self, other):
        return self.k == other.k and self.theteq == other.theteq

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Dihedral(_FourAtomTerm):
    """
    A valence dihedral between 4 atoms separated by three covalent bonds.

    Parameters (and Attributes)
    ---------------------------
    atom1 : Atom
        An atom on one end of the valence dihedral bonded to atom 2
    atom2 : Atom
        An atom in the middle of the valence dihedral bonded to atom1 and atom3
    atom3 : Atom
        An atom in the middle of the valence dihedral bonded to atom2 and atom4
    atom4 : Atom
        An atom on the other end of the valence dihedral bonded to atom 3
    improper : bool=False
        If True, this is an Amber-style improper torsion, where atom3 is the
        "central" atom bonded to atoms 1, 2, and 4 (atoms 1, 2, and 4 are _only_
        bonded to atom 3 in this instance)
    ignore_end : bool=False
        If True, the end-group interactions for this torsion are ignored, either
        because it is involved in a ring where the end-group atoms are excluded
        or because it is one term in a multi-term dihedral expansion
    type : DihedralType=None
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
        self._signs = [1, 1]
        if ignore_end: self._signs[0] = -1
        if improper: self._signs[1] = -1
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

    @property
    def signs(self):
        """
        For Amber topology files, the signs of the 3rd and 4th atoms indicate
        whether the end-group interactions (i.e., 1-4 nonbonded terms) are
        ignored or if the torsion is improper, respectively.

        This is a 2-element list with elements:
            If end-groups are ignored, signs[0] = -1
            If the torsion is improper, signs[1] = -1
        """
        return self._signs

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
        thing : Dihedral or other iterable
            A Dihedral or an iterable with 4 indexes that will be used to
            identify whether this torsion has the same atoms

        Returns
        -------
        bool - True if `thing` and `self` have the same atoms, False otherwise

        Notes
        -----
        This raises a `TypeError` if thing is not a Dihedral and is not iterable
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

    def __repr__(self):
        return "<Dihedral %r--%r--%r--%r>" % (self.atom1, self.atom2,
                self.atom3, self.atom4)

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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralType(_ListItem):
    """
    A dihedral type with a set of dihedral parameters

    Parameters (and Attributes)
    ---------------------------
    phi_k : float
        The force constant in kcal/mol
    per : int
        The dihedral periodicity
    phase : float
        The dihedral phase in degrees
    scee : float
        1-4 electrostatic scaling factor
    scnb : float
        1-4 Lennard-Jones scaling factor
    list : TrackedList=None
        A list of `DihedralType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : int
        The index of this DihedralType inside its containing list

    Notes
    -----
    Two `DihedralType`s are equal if their `phi_k`, `per`, `phase`, `scee`, and
    `scnb` attributes are equal

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
   
    def __init__(self, phi_k, per, phase, scee, scnb, list=None):
        """ DihedralType constructor """
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradley(object):
    """
    A Urey-Bradley angle type with a set of parameters. It has the same
    functional form as a bond, but it is defined between two atoms forming a
    valence angle separated by two bonds.

    Parameters (and Attributes)
    ---------------------------
    atom1 : Atom
        The first atom involved in the Urey-Bradley bond
    atom2 : Atom
        The other atom involved in the Urey-Bradley bond
    type : BondType=None
        The Urey-Bradley bond type that defines the parameters for this bond

    Notes
    -----
    You can test whether an Atom is contained within the bond using the `in`
    operator. A `BondError` is raised if `atom1` and `atom2` are identical. You
    can also test that a `Bond` is contained in this Urey-Bradley valence angle

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
    atom1 : Atom
        The central atom A1 in the schematic above
    atom2 : Atom
        An atom in the improper, A2 in the schematic above
    atom3 : Atom
        An atom in the improper, A3 in the schematic above
    atom4 : Atom
        An atom in the improper, A4 in the schematic above
    type : ImproperType=None
        The ImproperType object containing the parameters for this improper
        torsion

    Notes
    -----
    An Improper torsion can contain bonds or atoms. A bond is contained if it
    exists between atom 1 and any other atom. Raises `BondError` if any of the
    atoms are duplicates.

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
            if self.atom1 != thing.atom1:
                return False
            # Make a set with the remaining atoms. If they are equal, the
            # impropers are equivalent
            selfset = set([self.atom2, self.atom3, self.atom4])
            otherset = set([thing.atom2, thing.atom3, thing.atom3])
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ImproperType(_ListItem):
    """
    An improper type with a set of improper torsion parameters

    Parameters (and Attributes)
    ---------------------------
    psi_k : float
        Force constant in kcal/mol/radians^2
    psi_eq : float
        Equilibrium torsion angle in Degrees
    list : TrackedList=None
        A list of `ImproperType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : int
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
        self.psi_k = psi_k
        self.psi_eq = psi_eq
        self.list = list
        self._idx = -1

    def __eq__(self, other):
        return self.psi_k == other.psi_k and self.psi_eq == other.psi_eq

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Cmap(object):
    """
    A coupled-torsion correction map term defined between 5 atoms connected by
    four covalent bonds. This is a coupled-torsion potential in which the
    torsions are consecutive.

    Parameters (and Attributes)
    ---------------------------
    atom1 : Atom
        An atom on one end of the valence coupled-torsion bonded to atom2
    atom2 : Atom
        An atom in the middle of the CMAP bonded to atoms 1 and 3
    atom3 : Atom
        An atom in the middle of the CMAP bonded to atoms 2 and 4
    atom4 : Atom
        An atom in the middle of the CMAP bonded to atoms 3 and 5
    atom5 : Atom
        An atom in the middle of the CMAP bonded to atom 4
    type : CmapType=None
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
        atom4.cmaps.append(self)
        # Load the CMAP interpolation table
        self.type = type

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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CmapType(_ListItem):
    """
    A CMAP type with a potential energy interpoloation grid mapping out the 2-D
    potential of coupled torsions.

    Parameters
    ----------
    resolution : int
        The number of grid points in the correction map potential in both
        torsion dimensions
    grid : iterable of floats
        This must be a 1-dimensional list of grid values. The dimension must be
        `resolution*resolution`, and must be row-major (i.e., the second
        dimension changes the fastest) when indexed with 2 indices.
    list : TrackedList=None
        A list of `CmapType`s in which this is a member

    Attributes
    ----------
    resolution : int
        Potential grid resolution (see description in Parameters)
    grid : _CmapGrid
        A _CmapGrid object defining the interpolating potential energy grid,
        with each point having the units kcal/mol
    list : TrackedList=None
        If not None, this is a list in which this instance _may_ be a member
    idx : int
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
    def __init__(self, resolution, grid, list=None):
        self.resolution = resolution
        self.grid = _CmapGrid(resolution, grid)
        if len(grid) != self.resolution * self.resolution:
            raise CmapError('CMAP grid does not match expected resolution')
        self._idx = -1
        self.list = list

    def __eq__(self, other):
        return (self.resolution == other.resolution and
                all([abs(i - j) < TINY for i, j in zip(self.grid, other.grid)]))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _CmapGrid(object):
    """
    A grid object for storing Correction map data. Data can be accessed in one
    of two ways; either with 1 or 2 indexes. If 2 indexes [i,j] are given, the
    index into the flattened array is i*resolution+j. Indexing starts from 0.

    The _CmapGrid usually has ranges for the two angles from -180 to 180. Some
    places will expect the range to be 0-360 degrees (e.g., OpenMM). The
    switch_range method returns a _CmapGrid with this range. This method will
    also work to go backwards.

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
            for j in xrange(res):
                # Start from the middle
                newgrid[i, j] = self[(i+mid)%res, (j+mid)%res]
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

    Parameters (and Attributes)
    ---------------------------
    atom1 : Atom
        The first atom involved in the trigonal angle
    atom2 : Atom
        The central atom involved in the trigonal angle
    atom3 : Atom
        The third atom involved in the trigonal angle
    atom4 : Atom
        The fourth atom involved in the trigonal angle
    type : AngleType=None
        The angle type containing the parameters

    Notes
    -----
    Either `Atom`s or `Bond`s can be contained within this trigonal angle
    """
    def __init__(self, atom1, atom2, atom3, atom4, type=None):
        _FourAtomTerm.__init__(self)
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class OutOfPlaneBend(_FourAtomTerm):
    """
    Out-of-plane bending term in the AMOEBA force field. The bond pattern is the
    same as `TrigonalAngle`

    Parameters (and Attributes)
    ---------------------------
    atom1 : Atom
        The first atom involved in the trigonal angle
    atom2 : Atom
        The central atom involved in the trigonal angle
    atom3 : Atom
        The third atom involved in the trigonal angle
    atom4 : Atom
        The fourth atom involved in the trigonal angle
    type : AngleType=None
        The angle type containing the parameters

    Notes
    -----
    Either `Atom`s or `Bond`s can be contained within this trigonal angle
    """
    def __init__(self, atom1, atom2, atom3, atom4, type=None):
        _FourAtomTerm.__init__(self)
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
    are bonded _only_ to A3 and A4, respectively. Atoms A1 and A5 are each
    bonded to 3 other atoms.

    Parameters (and Attributes)
    ---------------------------
    atom1 : Atom
        atom A1 in the schematic above
    atom2 : Atom
        atom A2 in the schematic above
    atom3 : Atom
        atom A3 in the schematic above
    atom4 : Atom
        atom A4 in the schematic above
    atom5 : Atom
        atom A5 in the schematic above
    atom6 : Atom
        atom A6 in the schematic above
    type : DihedralType=None
        The parameters for this Pi-torsion

    Notes
    -----
    Both `Bond`s and `Atom`s can be contained in a pi-torsion
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
                    thing is self.atom5 or thing is self.atom5)
        # Assume Bond
        return ((self.atom2 in thing and self.atom3 in thing) or
                (self.atom1 in thing and self.atom3 in thing) or
                (self.atom3 in thing and self.atom4 in thing) or
                (self.atom4 in thing and self.atom5 in thing) or
                (self.atom4 in thing and self.atom6 in thing))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class StretchBend(object):
    """
    This term models the stretching and bending of a standard valence angle, and
    is used in the AMOEBA force field

    Parameters (and Attributes)
    ---------------------------
    atom1 : Atom
        The first atom on one end of the angle
    atom2 : Atom
        The central atom in the angle
    atom3 : Atom
        The atom on the other end of the angle
    type : StretchBendType=None
        The type containing the stretch-bend parameters

    Notes
    -----
    Both `Bond`s and `Atom`s can be contained in a stretch-bend term
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class StretchBendType(_ListItem):
    """
    A stretch-bend type with two distances and an angle in AMOEBA

    Parameters (and Attributes)
    ---------------------------
    k : float
        Force constant in kcal/mol/radians^2
    req1 : float
        Equilibrium bond distance for bond between the first and second atoms in
        Angstroms
    req2 : float
        Equilibrium bond distance for bond between the second and third atoms in
        Angstroms
    theteq : float
        Equilibrium angle in degrees
    list : TrackedList=None
        A list of `StretchBendType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : int
        The index of this StretchBendType inside its containing list

    Notes
    -----
    Two `StretchBendType`s are equal if their `req1`, `req2`, `theteq`, and `k`
    attributes are equal

    Examples
    --------
    >>> sbt1 = StretchBendType(10.0, 1.0, 1.0, 180.0)
    >>> sbt2 = StretchBendType(10.0, 1.0, 1.0, 180.0)
    >>> sbt1 is sbt2
    False
    >>> sbt1 == sbt2
    True
    >>> sbt1.idx # Not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> strbnd_list = []
    >>> strbnd_list.append(StretchBendType(10.0, 1.0, 1.0, 180.0, strbnd_list))
    >>> strbnd_list.append(StretchBendType(10.0, 1.0, 1.0, 180.0, strbnd_list))
    >>> strbnd_list[0].idx
    0
    >>> strbnd_list[1].idx
    1
    """
    def __init__(self, k, req1, req2, theteq, list=None):
        self.k = k
        self.req1 = req1
        self.req2 = req2
        self.theteq = theteq
        self._idx = -1
        self.list = list

    def __eq__(self, other):
        return (self.k == other.k and self.req1 == other.req1 and
                self.req2 == other.req2 and self.theteq == other.theteq)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TorsionTorsion(Cmap):
    """
    This is a coupled-torsion map used in the AMOEBA force field similar to the
    correction-map (CMAP) potential used by the CHARMM force field

    Parameters (and Attributes)
    ---------------------------
    atom1 : Atom
        An atom on one end of the valence torsion-torsion bonded to atom2
    atom2 : Atom
        An atom in the middle of the torsion-torsion bonded to atoms 1 and 3
    atom3 : Atom
        An atom in the middle of the torsion-torsion bonded to atoms 2 and 4
    atom4 : Atom
        An atom in the middle of the torsion-torsion bonded to atoms 3 and 5
    atom5 : Atom
        An atom in the middle of the torsion-torsion bonded to atom 4
    type : TorsionTorsionType=None
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

    Parameters (and Attributes)
    ---------------------------
    ang1 : list of floats
        Angles in the first dimension of the interpolation table
    ang2 : list of floats
        Angles in the second dimension of the interpolation table
    data : list of floats
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
        for x, y in zip(self.data, other.data):
            if x != y: return False
        return True

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TorsionTorsionType(_ListItem):
    """
    The type containing the parameter maps for the Amoeba torsion-torsion
    potentials. It contains the original potential as well as interpolated first
    and second derivatives for the AMOEBA force field.

    Parameters
    ----------
    dims : tuple of 2 ints
        The table dimensions
    ang1 : list of floats
        The list of angles in the first dimension
    ang2 : list of floats
        The list of angles in the second dimension
    f : list of floats
        The interpolation table for the energy
    dfda1 : list of floats
        The interpolation table of the gradient w.r.t. angle 1
    dfda2 : list of floats
        The interpolation table of the gradient w.r.t. angle 2
    d2fda1da2 : list of floats
        The interpolation table of the 2nd derivative w.r.t. both angles
    list : TrackedList=None
        The list containing this coupled torsion-torsion map

    Attributes
    ----------
    dims : tuple of 2 ints
        The table dimensions
    ang1 : list of floats
        The list of angles in the first dimension
    ang2 : list of floats
        The list of angles in the second dimension
    f : _TorTorTable
        The interpolation table for the energy as a _TorTorTable
    dfda1 : _TorTorTable
        The interpolation table for the first gradient as a _TorTorTable
    dfda2 : _TorTorTable
        The interpolation table for the second gradient as a _TorTorTable
    d2fda1da2 : _TorTorTable
        The interpolation table for the second derivative as a _TorTorTable
    list : TrackedList
        The list that may, or may not, contain this TorsionTorsionType
    idx : int
        The index of this item in the list or iterable defined by `list`

    Notes
    -----
    Since the derivatives are uniquely determined by the original potential,
    equality between two coupled-coupled torsions can be uniquely determined by
    comparing the energy table. The other tables are assumed to follow suit.
    """
    def __init__(self, dims, ang1, ang2, f, dfda1, dfda2, d2fda1da2, list=None):
        if len(dims) != 2:
            raise ValueError('dims must be a 2-dimensional iterable')
        if len(ang1) != dims[0] or len(ang2) != dims[1]:
            raise ValueError('dims does match the angle definitions')
        self.dims = tuple(dims)
        self.ang1 = ang1
        self.ang2 = ang2
        self.f = _TorTorTable(ang1, ang2, f)
        self.dfda1 = _TorTorTable(ang1, ang2, dfda1)
        self.dfda2 = _TorTorTable(ang1, ang2, dfda2)
        self.d2fda1da2 = _TorTorTable(ang1, ang2, d2fda1da2)
        self._idx = -1
        self.list = None

    def __eq__(self, other):
        if self.dims != other.dims: return False
        if self.ang1 != other.ang1: return False
        if self.ang2 != other.ang2: return False
        return self.f == other.f

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ChiralFrame(object):
    """
    A chiral frame as defined in the AMOEBA force field. It defines the frame of
    reference for a chiral center

    Parameters (and Attributes)
    ---------------------------
    atom1 : Atom
        The first atom defined in the chiral frame
    atom2 : Atom
        The second atom defined in the chiral frame
    chirality : int
        Either 1 or -1 to identify directionality

    Notes
    -----
    A chiral frame can only contain atoms
    """
    def __init__(self, atom1, atom2, chirality):
        self.atom1 = atom1
        self.atom2 = atom2
        self.chirality = chirality

    def __contains__(self, thing):
        return thing is self.atom1 or thing is self.atom2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class MultipoleFrame(object):
    """
    This defines the frame of reference for computing multipole interactions in
    the AMOEBA force field.

    Parameters (and Attributes)
    ---------------------------
    atom : Atom
        The atom for which the frame of reference is defined
    frame_pt_num : int
        The frame point number
    vectail : int
        The vector tail index
    vechead : int
        The vector head index
    nvec : int
        The number of vectors
    """
    def __init__(self, atom, frame_pt_num, vectail, vechead, nvec):
        self.atom = atom
        self.frame_pt_num = frame_pt_num
        self.vectail = vectail
        self.vechead = vechead
        self.nvec = nvec

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Residue(_ListItem):
    """
    A single residue that is composed of a small number of atoms

    Parameters
    ----------
    name : str
        Name of the residue. Typical convention is to choose a name that is 4
        characters or shorter
    number : int
        Residue number assigned in the input structure
    list : TrackedList=None
        List of residues in which this residue is a member

    Attributes
    ----------
    name : str
        The name of this residue
    number : int=-1
        The number of this residue in the input structure
    idx : int
        The index of this residue inside the container. If this residue has no
        container, or it is not present in the container, idx is -1
    list : TrackedList
        The container that _may_ have this residue contained inside
    atoms : list of Atom instances
        This is the list of `Atom`s that make up this residue

    Notes
    -----
    - Iterating over a residue will iterate over the atoms. It is exactly
      equivalent to iterating over the `atoms` attribute
    - Supports testing if an Atom instance is contained `in` this residue
    - `len()` returns the number of atoms in this residue
    """

    def __init__(self, name, number=-1, list=None):
        self.name = name
        self.number = number
        self.list = list
        self._idx = -1
        self.atoms = []

    def add_atom(self, atom):
        """ Adds an atom to this residue

        Parameters
        ----------
        atom : Atom
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
        True if there are no atoms left. False if this residue still has atoms
        """
        return not bool(self.atoms)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _changes(func):
    """ Decorator to indicate the list has changed """
    def new_func(self, *args):
        self.changed = True
        self.needs_indexing = True
        return func(self, *args)
    return new_func

class TrackedList(list):
    """
    This creates a list type that allows you to see if anything has changed

    Attributes
    ----------
    changed : bool
        Determines if something has been done to fundamentally change the
        underlying topology defined by this list such that the topology needs to
        be rebuilt
    needs_indexing : bool
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
        list.__init__(self, arg)

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
            except IndexError:
                print(index)
                raise
            try:
                self[index].list = None
            except AttributeError:
                pass

        list.__delitem__(self, item)

    @_changes
    def pop(self, idx):
        item = list.pop(self, idx)
        try:
            item._idx = -1
        except IndexError:
            # Must be an immutable type, so don't complain
            pass
        return item

    append = _changes(list.append)
    extend = _changes(list.extend)
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
        try:
            return TrackedList(retval)
        except TypeError:
            return retval

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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ResidueList(TrackedList):
    """ Array of `Residue` instances """

    def add_atom(self, atom, resname, resnum):
        """
        Adds a new atom to the ResidueList, adding a new residue to this list if
        it has a different name or number as the last residue

        Parameters
        ----------
        atom : Atom
            The atom to add to this residue list
        resname : str
            The name of the residue this atom belongs to
        resnum : int
            The number of the residue this atom belongs to

        Notes
        -----
        If the residue name and number differ from the last residue in this
        list, a new residue is added and the atom is added to that residue
        """
        try:
            last = self[-1]
        except IndexError:
            # Empty list -- add our first residue
            new_res = Residue(resname, resnum, list=self)
            new_res.add_atom(atom)
            self.append(new_res)
        else:
            if last.number != resnum or last.name != resname:
                new_res = Residue(resname, resnum, list=self)
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
        for i in xrange(len(self)-1, -1, -1):
            res = self[i]
            if res.empty(): del self[i]

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

    def unmark(self):
        """ Unmark all atoms in this list """
        for atm in self: atm.marked = 0

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

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)
