"""
This module contains objects that deal with system topology and parameters, such
as atoms, residues, bonds, angles, etc.

by Jason Swails
"""
from __future__ import division

from chemistry.exceptions import (BondError, DihedralError, AmberParmError,
                                  CmapError)
from chemistry.amber.constants import NATOM, TINY
from chemistry.periodic_table import AtomicNum, Mass, Element as _Element
from compat24 import all, property
import warnings

__all__ = ['Atom', 'Bond', 'BondType', 'Angle', 'AngleType', 'Dihedral',
           'DihedralType', 'Residue', 'ResidueList', 'AtomList', 'BondTypeList',
           'AngleTypeList', 'DihedralTypeList', 'TrackedList']

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
    changed (so we know when to update indexes), this is a fast lookup.
    Otherwise, it can be slow.
    """

    @property
    def idx(self):
        try:
            mylist = self.list
        except AttributeError:
            return -1

        try:
            list_changed = mylist.changed
        except AttributeError:
            # This isn't a tracked list, so just look through the list
            for i, item in enumerate(mylist):
                if item is self: return i
            return -1
        else:
            if list_changed or self._idx == -1:
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
        list of UreyBradley objects in which this atom is a member
    impropers : list of Improper instances
        list of Improper objects in which the atom is a member
    cmaps : list of Cmap instances
        list of Cmap objects in which the atom is a member
    bond_partners : list of Atom instances
        list of Atoms to which this atom is bonded.
    angle_partners : list of Atom instances
        list of Atoms to which this atom forms an angle, but not a bond
    dihedral_partners : list of Atom instances
        list of Atoms to which this atom forms an dihedral, but not a bond or
        angle
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
    >>> a1 is in a2.bond_partners and a2 is in a1.bond_partners
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
                 charge=0.0, mass=0.0, nb_idx=0, radii=0.0, tree='BLA'):
        self.list = list
        self._idx = -1
        self.atomic_number = atomic_number
        self.name = name
        self.type = type
        self.charge = charge
        self.mass = mass
        self.nb_idx = nb_idx
        self.radii = radii
        self.tree = tree
        self._bond_partners = set()
        self._angle_partners = set()
        self._dihedral_partners = set()
        self._exclusion_partners = set() # For arbitrary/other exclusions
        self.residue = None
        self.marked = 0 # For setting molecules
        self.bonds, self.angles, self.dihedrals = [], [], []
        self.urey_bradleys, self.impropers, self.cmaps = [], [], []
   
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
        return sorted(list(self._bond_partners))

    @property
    def angle_partners(self):
        """ List of all angle partners that are NOT bond partners """
        return sorted(list(self._angle_partners - self._bond_partners))

    @property
    def dihedral_partners(self):
        " List of all dihedral partners that are NOT angle or bond partners "
        return sorted(list(self._dihedral_partners - self._angle_partners -
                           self._bond_partners))

    @property
    def exclusion_partners(self):
        """
        List of all exclusions not otherwise excluded by bonds/angles/torsions
        """
        return sorted(list(self._exclusion_partners - self._dihedral_partners -
                           self._angle_partners - self._bond_partners))

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
        self._bond_partners.add(other)

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
        self._angle_partners.add(other)
   
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
        self._dihedral_partners.add(other)
      
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
        self._exclusion_partners.add(other)
        # If he is excluded from me, then I am excluded from him
        other._exclusion_partners.add(self)

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
        self._bond_partners = set()
        self._angle_partners = set()
        self._dihedral_partners = set()
        if not keep_exclusions:
            self._exclusion_partners = set()

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
        atom1.angle_to(atom3)

    def __contains__(self, thing):
        """ Quick and easy way to see if an Atom or a Bond is in this Angle """
        if isinstance(thing, Atom):
            return (thing is self.atom1 or thing is self.atom2 or
                    thing is self.atom3)
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom1 in thing))

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

class Dihedral(object):
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
    >>> Bond(a1, a4) in angle # this is not part of the angle definition
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
        # Make sure we're not dihedraling me to myself
        atmlist = [atom1, atom2, atom3, atom4]
        for i in xrange(len(atmlist)):
            for j in xrange(i+1, len(atmlist)):
                if atmlist[i] is atmlist[j]:
                    raise BondError('Cannot dihedral atom to itself!')
        # Set up instances
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
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
        atom1.dihedral_to(atom4)

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
            return (thing is self.atom1 or thing is self.atom2 or
                    thing is self.atom3 or thing is self.atom4)
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
    >>> dt.idx # not part of any list or iterable
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
    type : UreyBradleyType=None
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradleyType(_ListItem):
    """
    A Urey-Bradley bond type with a set of bond parameters

    Parameters (and Attributes)
    ---------------------------
    k : float
        Force constant in kcal/mol/Angstrom^2
    req : float
        Equilibrium distance in Angstroms
    list : TrackedList=None
        A list of `UreyBradleyType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : int
        The index of this UreyBradleyType inside its containing list

    Notes
    -----
    Two `UreyBradleyType`s are equal if their `k` and `req` attributes are
    equal

    Examples
    --------
    >>> ubt1 = UreyBradleyType(10.0, 1.0)
    >>> ubt2 = UreyBradleyType(10.0, 1.0)
    >>> ubt1 is ubt2
    False
    >>> ubt1 == ubt2
    True

    As part of a list, they can be indexed

    >>> urey_list = []
    >>> urey_list.append(UreyBradleyType(10.0, 1.0, list=urey_list))
    >>> urey_list.append(UreyBradleyType(10.0, 1.0, list=urey_list))
    >>> urey_list[0].idx
    0
    >>> urey_list[1].idx
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

class Improper(object):
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
        # Make sure we're not dihedraling me to myself
        atmlist = [atom1, atom2, atom3, atom4]
        for i in xrange(len(atmlist)):
            for j in xrange(i+1, len(atmlist)):
                if atmlist[i] is atmlist[j]:
                    raise BondError('Cannot improper atom to itself!')
        # Set up instances
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        # Log these impropers in each atom
        atom1.impropers.append(self)
        atom2.impropers.append(self)
        atom3.impropers.append(self)
        atom4.impropers.append(self)
        # Load the force constant and equilibrium angle
        self.type = type

    def __contains__(self, thing):
        """
        Quick and easy way to find out if an Atom or Bond is in this Improper
        """
        if isinstance(thing, Atom):
            return (thing is self.atom1 or thing is self.atom2 or
                    thing is self.atom3 or thing is self.atom4)
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
    >>> it.idx # Not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> improper_list = []
    >>> improper_list.append(ImproperType(10.0, 180.0, list=improper_type))
    >>> improper_list.append(ImproperType(10.0, 180.0, list=improper_type))
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
    >>> at1, at2, at3, at4, at5 = Atom(), Atom(), Atom(), Atom(), Atom()
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
                indices = idx.indices(len(self._data))
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

class Residue(_ListItem):
    """
    A single residue that is composed of a small number of atoms

    Parameters
    ----------
    name : str
        Name of the residue. Typical convention is to choose a name that is 4
        characters or shorter
    list : TrackedList=None
        List of residues in which this residue is a member

    Attributes
    ----------
    name : str
        The name of this residue
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

    def __init__(self, name, list=None):
        self.name = name
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _tracking(fcn):
    """ Decorator to indicate the list has changed """
    def new_fcn(self, *args):
        self.changed = True
        return fcn(self, *args)
    return new_fcn

class TrackedList(list):
    """
    This creates a list type that allows you to see if anything has changed
    """
    def __init__(self, arg=[]):
        self.changed = True
        list.__init__(self, arg)

    @_tracking
    def __delitem__(self, item):
        """ Deletes items and slices. Make sure all items """
        try:
            indices = item.indices(len(self))
        except AttributeError:
            indices = [item]

        try:
            for index in indices:
                self[index]._idx = -1
        except AttributeError:
            pass

        list.__delitem__(self, item)

    append = _tracking(list.append)
    extend = _tracking(list.extend)
    __setitem__ = _tracking(list.__setitem__)
    __iadd__ = _tracking(list.__iadd__)
    __imul__ = _tracking(list.__imul__)
    pop = _tracking(list.pop)

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
        try:
            for i, item in enumerate(self):
                item._idx = i
            self.changed = False
        except AttributeError:
            # This must be some kind of immutable type, so don't worry about it
            self.changed = False

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ResidueList(TrackedList):
    """ Array of Residues. """

    def __init__(self, parm):
        list.__init__(self, [Residue(parm.parm_data['RESIDUE_LABEL'][i], i+1)
                             for i in xrange(parm.ptr('nres'))])
        for i, val in enumerate(parm.parm_data['RESIDUE_POINTER']):
            start = val - 1
            try:
                end = parm.parm_data['RESIDUE_POINTER'][i+1] - 1
            except IndexError:
                end = parm.parm_data['POINTERS'][NATOM]
            for j in xrange(start, end):
                self[i].add_atom(parm.atom_list[j])
        self.parm = parm
   
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AtomList(list):
    """ Array of Atoms """
    #===================================================

    def __init__(self, parm, fill_from=None):
        self.parm = parm
        if fill_from is None:
            list.__init__(self, [Atom(self.parm, i) for i in
                                xrange(self.parm.ptr('natom'))])
        else:
            list.__init__(self, [0 for i in xrange(self.parm.ptr('natom'))])
            for i, atm in enumerate(fill_from): self[i] = atm
        self.changed = False

    #===================================================

    def __delitem__(self, idx):
        """ Deletes this atom then re-indexes everybody else """
        self[idx].idx = -1
        # Delete this atom from its residue as well
        self[idx].residue.delete_atom(self[idx])
        list.__delitem__(self, idx)
        self.changed = True

    #===================================================
   
    def unmark(self):
        """ Unmark all atoms in this list """
        for atm in self: atm.marked = 0

    #===================================================
   
    def _index_us(self):
        """ We have deleted an atom, so now we have to re-index everybody """
        for i, atom in enumerate(self): atom.idx = atom.starting_index = i

    #===================================================

    def append(self, item):
        """ Don't allow this! """
        raise AmberParmError("Cannot add to an AtomList!")

    #===================================================

    def find_extra_exclusions(self):
        " Load all extra exclusions that may be stored in the topology file "
        first = 0
        for atom in self:
            first += atom.load_exclusions(first)

    #===================================================

    def write_to_parm(self):
        """ Writes all of the atom data to the topology file """
        # Write all of the arrays here
        self.parm.parm_data['POINTERS'][NATOM] = len(self)
        # Array slices are faster than copy() and creating new arrays
        # each time
        zeros = [0 for i in xrange(self.parm.parm_data['POINTERS'][NATOM])]
        self.parm.parm_data['ATOM_NAME'] = zeros[:]
        self.parm.parm_data['CHARGE'] = zeros[:]
        self.parm.parm_data['MASS'] = zeros[:]
        self.parm.parm_data['ATOM_TYPE_INDEX'] = zeros[:]
        self.parm.parm_data['NUMBER_EXCLUDED_ATOMS'] = zeros[:]
        self.parm.parm_data['AMBER_ATOM_TYPE'] = zeros[:]
        self.parm.parm_data['JOIN_ARRAY'] = zeros[:]
        self.parm.parm_data['TREE_CHAIN_CLASSIFICATION'] = zeros[:]
        self.parm.parm_data['IROTAT'] = zeros[:]
        self.parm.parm_data['RADII'] = zeros[:]
        self.parm.parm_data['SCREEN'] = zeros[:]
        self._index_us()
        self._determine_exclusions()
        for atm in self: 
            atm.add_data()
            atm.starting_index = atm.idx # arrays are updated...

    #===================================================

    def _determine_exclusions(self):
        """
        Figures out the EXCLUDED_ATOMS_LIST. Only do this right before you write
        the topology file, since it's expensive
        """
        self.parm.parm_data['EXCLUDED_ATOMS_LIST'] = []
        # We have to do something different for extra points. See the top of
        # extra_pts.f in the sander src/ directory. Effectively, the EP is
        # excluded from every atom that the atom it's attached to is excluded
        # from. So here we go through and exclude every extra point from every
        # other atom my bonded pair is excluded from:
        for atm in self:
            if not atm.attype[:2] in ['EP', 'LP']: continue
            partner = atm.bond_partners[0]
            # Now add all bond partners
            for patm in partner.bond_partners:
                # Don't add myself
                if patm is atm: continue
                atm.exclude(patm)
            # Now add all angle partners
            for patm in partner.angle_partners: atm.exclude(patm)
            # Now add all dihedral partners
            for patm in partner.dihedral_partners: atm.exclude(patm)
            # Now add all other arbitrary exclusions
            for patm in partner.exclusion_partners:
                if patm is atm: continue
                atm.exclude(patm)

        for atm in self:
            vals_to_add = []
            for member in atm.bond_partners:
                if member.idx > atm.idx: vals_to_add.append(member.idx+1)
            for member in atm.angle_partners:
                if member.idx > atm.idx: vals_to_add.append(member.idx+1)
            for member in atm.dihedral_partners:
                if member.idx > atm.idx: vals_to_add.append(member.idx+1)
            # Enable additional (arbitrary) exclusions
            for member in atm.exclusion_partners:
                if member.idx > atm.idx: vals_to_add.append(member.idx+1)
            vals_to_add.sort()
            # See comment above about numex = 0 --> numex = 1
            if not vals_to_add: vals_to_add = [0]
            self.parm.parm_data['EXCLUDED_ATOMS_LIST'].extend(vals_to_add)

    #===================================================

    def refresh_data(self):
        """
        Re-loads all data in parm.parm_data for each atom in case we changed
        any of it
        """
        for atm in self:
            atm.load_from_parm()

    #===================================================

    def __setitem__(self, idx, thing):
        """
        This means we changed things...
        """
        self.changed = True
        list.__setitem__(self, idx, thing)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _TypeList(TrackedList):
    """ Base class for all type lists """

    #===================================================

    def __init__(self, parm):
        """ Constructs a list of bond types from the topology file """
        self.parm = parm
        self._make_array()
        self.changed = False

    #===================================================

    def _make_array(self):
        """ 
        This method fills self with whichever element we need. This MUST be
        overwritten, so I force it here
        """
        raise NotImplemented('Subclasses must implement _make_array')

    #===================================================

    def write_to_parm(self):
        """ Writes the data here to the parm data """
        for item in self: item.write_info(self.parm)

    #===================================================

    def reset(self):
        """ Reset indexes to -1 to allow multiple remake_parm calls """
        for thing in self: thing.idx = -1

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondTypeList(_TypeList):
    """ Bond type list """

    #===================================================

    def _make_array(self):
        kl = self.parm.parm_data['BOND_FORCE_CONSTANT']
        eql = self.parm.parm_data['BOND_EQUIL_VALUE']
        list.__init__(self, [BondType(k, eq, -1) for k, eq in zip(kl, eql)])
      
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleTypeList(_TypeList):
    """ Angle type list """

    #===================================================

    def _make_array(self):
        kl = self.parm.parm_data['ANGLE_FORCE_CONSTANT']
        eql = self.parm.parm_data['ANGLE_EQUIL_VALUE']
        list.__init__(self, [AngleType(k, eq, -1) for k, eq in zip(kl, eql)])
      
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralTypeList(_TypeList):
    """ Dihedral type list """

    #===================================================

    def _make_array(self):
        kl = self.parm.parm_data['DIHEDRAL_FORCE_CONSTANT']
        pel = self.parm.parm_data['DIHEDRAL_PERIODICITY']
        phl = self.parm.parm_data['DIHEDRAL_PHASE']
        if (not 'SCEE_SCALE_FACTOR' in self.parm.parm_data.keys() or
            not 'SCNB_SCALE_FACTOR' in self.parm.parm_data.keys()):
            list.__init__(self,
                    [DihedralType(k, pe, ph, 1.2, 2.0, -1)
                     for k, pe, ph in zip(kl, pel, phl)]
            )
        else:
            scel = self.parm.parm_data['SCEE_SCALE_FACTOR']
            scnl = self.parm.parm_data['SCNB_SCALE_FACTOR']
            list.__init__(self,
                    [DihedralType(k, pe, ph, sce, scn, -1)
                     for k, pe, ph, sce, scn in zip(kl, pel, phl, scel, scnl)]
            )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradleyTypeList(_TypeList):
    """ Urey-Bradley type list """

    #===================================================

    def _make_array(self):
        kl = self.parm.parm_data['CHARMM_UREY_BRADLEY_FORCE_CONSTANT']
        eql = self.parm.parm_data['CHARMM_UREY_BRADLEY_EQUIL_VALUE']
        list.__init__(self,
                [UreyBradleyType(k, eq, -1) for k, eq in zip(kl, eql)]
        )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ImproperTypeList(_TypeList):
    """ CHARMM Improper torsion type list """

    #===================================================

    def _make_array(self):
        kl = self.parm.parm_data['CHARMM_IMPROPER_FORCE_CONSTANT']
        pl = self.parm.parm_data['CHARMM_IMPROPER_PHASE']
        list.__init__(self, [ImproperType(k, p, -1) for k, p in zip(kl, pl)])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CmapTypeList(_TypeList):
    """ CHARMM correction-map type list """

    #===================================================

    def _make_array(self):
        # Need to get all of the CMAP types
        ncmaps = self.parm.parm_data['CHARMM_CMAP_COUNT'][1]
        list.__init__(self)
        for i in xrange(ncmaps):
            res = self.parm.parm_data['CHARMM_CMAP_RESOLUTION'][i]
            grid = self.parm.parm_data['CHARMM_CMAP_PARAMETER_%02d' % (i+1)]
            cmts = self.parm.parm_comments['CHARMM_CMAP_PARAMETER_%02d' % (i+1)]
            list.append(self, CmapType(res, grid, cmts, -1))

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
    doctest.testmod()
