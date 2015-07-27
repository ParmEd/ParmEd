The :mod:`topologyobjects` module
=================================

The :mod:`topologyobjects` module contains a wide array of classes used
extensively in the core :class:`structure.Structure` class.

.. currentmodule:: parmed.topologyobjects
.. autosummary::
    :toctree: topobj/

    Atom
    Residue
    AtomList
    ResidueList
    Bond
    Angle
    UreyBradley
    Dihedral
    Improper
    Cmap
    AtomType
    BondType
    AngleType
    DihedralType
    DihedralTypeList
    ImproperType

Using the :class:`Atom` class
-----------------------------

The :class:`Atom` class behaves much as you would expect it to.  It has various
attributes, including the ``name`` of the atom, the ``type`` of the atom, its
``mass``, ``atomic_number``, and many other properties typically stored in one
of the supported file formats.  For example::

    >>> atom = Atom(name="CA", mass=12.01, type="CT", atomic_number=12)
    >>> atom.name
    'CA'
    >>> atom.type
    'CT'
    >>> atom.atomic_number
    12
    >>> atom.mass
    12.01

* It also stores references to all of the parameter types in which this belongs.
  For instance, the ``bonds``, ``angles``, ``dihedrals``, ``impropers``, and
  ``cmaps`` attributes are *lists* of references to the :class:`Bond`,
  :class:`Angle`, :class:`Dihedral`, :class:`Improper`, and :class:`Cmap`
  instances in which this atom is a member.

* Instances of :class:`Atom` should be considered *singletons* -- that is, their
  equality is based on identity. As such, you should use the ``is`` operator
  rather than ``==`` when comparing the equality of two atoms. In fact,
  ``Atom.__eq__`` is not implemented at all, so it works the same way that
  ``is`` does in this case. This is an important distinction, because it makes
  :class:`Atom` instances *hashable*, meaning that they can serve as keys in a
  ``dict`` or can be added to a ``set``.

* Atom ordering is determined based on positions within a particular
  :class:`AtomList`. When an :class:`Atom` belongs to an :class:`AtomList`, it
  knows its absolute position within that list via the :attr:`idx` attribute
  (which is set to ``-1`` if the Atom is not part of an :class:`AtomList`). The
  :attr:`idx` attribute is automatically updated if the parent :class:`AtomList`
  is changed in any way (e.g., some atoms are deleted or inserted). This allows
  sublists of :class:`Atom` instances to be sorted according to their order in
  their main :class:`AtomList`. (Note, each :class:`Atom` can belong to no more
  than one :class:`AtomList`)

* Atoms store lists of other :class:`Atom` instances that are connected to that
  atom via covalent bonds, angles, torsions, or coupled-torsions maps in the
  ``bond_partners``, ``angle_partners``, ``dihedral_partners``, and
  ``tortor_partners``, respectively. These are primarily used to determine pairs
  of atoms that are excluded from a nonbonded interaction in traditional force
  fields. There is an additional ``exclusion_partners`` array that contains a
  list of arbitrary pairs of atoms that are not connected by any of the valence
  terms mentioned above. Each atom can appear in one, and only one list of
  another atom, with preference given to the "partner" arrays in the order they
  were listed above. Membership in these arrays is always symmetric---that is,
  if one atom is in the ``bond_partners`` array of another, that other atom is
  also in the ``bond_partners`` array of the first. These lists are populated
  automatically when the relevant valence terms are created (see below).

You will typically not have to instantiate an :class:`Atom` or :class:`AtomList`
instance on your own, as they are created as part of a
:class:`structure.Structure` instance created by the various structure file
parsers.

Using the :class:`Residue` class
--------------------------------

The :class:`Residue` class is an object that represents a typical *residue* in a
biomolecular system. Examples include a single amino or nucleic acid in a
peptide or nucleotide chain.  A :class:`Residue` is a collection of atoms, each
of which is connected to the other ones through a network of bonds (and may be
bonded to other residues as well, although that is not required).

Residues have 4 primary attributes of general use:

    * ``name`` : The name of the residue
    * ``chain`` : The chain identifier of the residue (a common identifier in
      the PDB)
    * ``insertion_code`` : Also used as part of the PDB, these are primarily
      used to align homologous sequences when some homologies have an extra
      residue compared to the other(s)
    * ``atoms`` : A ``list`` of :class:`Atom` instances contained within the
      residue

Residues are uniquely identified by their ``name``, ``chain``, and
``insertion_code`` within a single system.

Common things you can do with a :class:`Residue`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the following code sample, assume ``atom`` is an instance of :class:`Atom`
which is one of the atoms inside ``residue`` (an instance of :class:`Residue`).

1. Iterate through the atoms::

    for atom in residue:
        print(atom.name)

2. Get the i'th atom from a residue::

    atomI = residue[i]

3. Determine if an atom is part of a residue::

    if atom in residue:
        print("Atom %s is inside residue %s" % (atom.name, residue.name))

4. Find out how many atoms are in a residue::

    n_atoms_in_res = len(residue)

5. Add an atom to a residue::

    residue.add_atom(atom)
    assert atom in residue

6. Delete an atom from a residue::

    residue.delete_atom(atom)
    assert atom not in residue

The :class:`AtomList` class
---------------------------

:class:`AtomList` is a subclass of the basic Python *list* type, with some added
features:

* The :class:`AtomList` claims ownership of any :class:`Atom` added to it and
  sets the ``idx`` attribute of each of the atoms to their position in the
  array. Atoms should never be part of multiple :class:`AtomList` instances.

* The :class:`AtomList` has a ``changed`` attribute that is set to ``True`` any
  time the list is changed (i.e., if an atom is added, deleted, inserted, moved,
  etc.) -- you must manaully set this ``changed`` attribute to ``False``
  afterwards if you wish to track any future changes.

This class is also typically instantiated by the :class:`structure.Structure`
parsers (e.g., :func:`structure.read_PDB`, :func:`structure.read_CIF`, ...)

The :class:`ResidueList` class
------------------------------

:class:`ResidueList` is also a subclass of *list* with some added features of
its own:

* The :class:`ResidueList` contains a list of :class:`Residue` instances, and
  sets the ``idx`` attribute on them to point to their current location in the
  list.

* Atoms can be added directly to the residue list. This is done because most
  file formats list atoms sequentially, and the beginning and end of each
  residue is detected by either the residue name, number, chain ID, or insertion
  code changing. As long as these residue properties are the same as the last
  :class:`Atom` instance you added, the atom will be added to that last residue.
  If it is a new residue, :class:`ResidueList` will create a new
  :class:`Residue` object, append it to itself, and add that atom to it. For
  example::

    >>> residues = ResidueList()
    >>> residues.add_atom(Atom(name='C1'), 'RE1', 1, 'A')
    >>> residues.add_atom(Atom(name='C2'), 'RE1', 1, 'A')
    >>> residues.add_atom(Atom(name='C3'), 'RE1', 1, 'A')
    >>> residues.add_atom(Atom(name='C4'), 'RE2', 2, 'A') # new residue
    >>> residues.add_atom(Atom(name='C5'), 'RE2', 3, 'A') # new residue
    >>> len(residues)
    3
    >>> len(residues[0])
    3
    >>> len(residues[1])
    1
    >>> len(residues[2])
    1

The :class:`Bond` class
-----------------------

The covalent bond is the simplest of topological features in a molecule. The
:class:`Bond` class represents a covalent bonds between two :class:`Atom`
instances. Creating a :class:`Bond` between two :class:`Atom` instances adds
the resulting :class:`Bond` to the ``bonds`` list of both atoms. It also adds
each atom to the ``bond_partners`` list of the other atom. For example::

    >>> a1 = Atom(name='CA')
    >>> a2 = Atom(name='CB')
    >>> bond = Bond(a1, a2)
    >>> a1.bonds[0] is bond
    True
    >>> a2.bonds[0] is bond
    True
    >>> a2 in a1.bond_partners
    True
    >>> a1 in a2.bond_partners
    True

Bonds can also *contain* :class:`Atom` instances. Continuing from the example
above::

    >>> a1 in bond
    True
    >>> a2 in bond
    True
    >>> a3 = Atom(name='HA')
    >>> a3 in bond
    False

Bonds are unique valence terms since they are the simplest of valence terms.
Other valence terms (described below) can *contain* both :class:`Atom` and
:class:`Bond` instances.

The ordering of the two atoms in the :class:`Bond` constructor does not matter,
as the bond ``a1--a2`` and ``a2--a1`` are equivalent.

The :class:`Angle` class
------------------------

The valence angle is formed from two adjacent bonds that form a common center.
The :class:`Angle` class represents a valence angle between 3 atoms in which two
of the three atoms are bonded to a central atom. In this case, the order of the
*middle* atom in the constructor matters, but the order of the outer atoms does
not. For example::

    >>> a1 = Atom(name='CA')
    >>> a2 = Atom(name='CB')
    >>> a3 = Atom(name='CG')
    >>> b1 = Bond(a1, a2)
    >>> b2 = Bond(a2, a3)
    >>> b3 = Bond(a1, a3) # let's make this a triangle
    >>> angle = Angle(a1, a2, a3)

In the above example, our angle is formed between the bonds ``a1--a2--a3``, so
``angle`` contains the three atoms as well as ``b1`` and ``b2``, but *not*
``b3``::

    >>> a1 in angle
    True
    >>> a2 in angle
    True
    >>> a3 in angle
    True
    >>> b1 in angle
    True
    >>> b2 in angle
    True
    >>> b3 in angle
    False

Like with the :class:`Bond` class, creating an :class:`Angle` adds the created
angle to the ``angles`` list on all three of the atoms.  Furthermore, the
``angle_partners`` attribute of ``a1`` and ``a3`` have ``a3`` and ``a1`` added
to them, respectively. As mentioned above, the same :class:`Atom` instance will
never appear in two ``xxx_partners`` arrays. For example::

    >>> a1 in a3.angle_partners
    True
    >>> a1 in a3.bond_partners
    False
    >>> a1 in a2.bond_partners
    True
    >>> a1 in a3.angle_partners
    False

The :class:`UreyBradley` class
------------------------------

The :class:`UreyBradley` class defines a covalent bond-like interaction between
two atoms separated by a valence angle, adding anharmonicity to purely harmonic
angle terms. The interface to :class:`UreyBradley` terms is analogous to those
for :class:`Bond` and :class:`Angle` above.

The :class:`Dihedral` class
---------------------------

The :class:`Dihedral` class follows the same trend as the :class:`Angle` class
described above, and behaves in an analogous way. There are a couple unique
aspects to Dihedrals, though. A :class:`Dihedral` can be an improper dihedral
(in which all 4 atoms are ideally coplanar, with the third atom being the
"central" atom bonded to all of the others). The :class:`Dihedral` instance has
a boolean attribute ``improper`` that indicates whether or not it is an improper
torsion.

The :class:`Improper` class
---------------------------

The :class:`Improper` class is *strictly* an improper torsion, and is used where
certain force fields utilize a separate functional form for "proper" and
"improper" torsions (e.g., the CHARMM force field). You are forwarded to the API
documentation for :class:`Improper` for more details.

The :class:`Cmap` class
-----------------------

The :class:`Cmap` class is a correction map defined between adjacent torsions
(i.e., 5 atoms separated by 4 bonds). Like the other valence terms described
above, it contains 5 :class:`Atom` instances and may contain 4 :class:`Bond`
instances (as tested via the ``in`` operator). See the API documentation for
:class:`Cmap` for more details.

The ``XyzType`` classes
-----------------------

Every valence term described above contains a ``type`` attribute that is either
``None`` if unassigned or the relevant ``ParameterType`` (e.g.,
:class:`BondType` for :class:`Bond` and :class:`UreyBradley` classes) instance.
Certain functionality (e.g., simulations using OpenMM or writing Amber topology
files) requires all parameters be present and assigned. These types have, where
applicable, force constants and equilibrium values for the various parameter
types.

Amoeba-specific terms
---------------------

The AMOEBA force field has nominal support in ParmEd (with improved support
planned). This force field contains a large number of additional terms and
parameters for use with the AMOEBA force field. For convenience, they are listed
below. Consult the API documentation in the resulting links for more
information.

.. currentmodule:: parmed.topologyobjects
.. autosummary::
    :toctree: topobj/

    TrigonalAngle
    OutOfPlaneBend
    OutOfPlaneBendType
    PiTorsion
    StretchBend
    StretchBendType
    TorsionTorsion
    TorsionTorsionType
    ChiralFrame
    MultipoleFrame
    NonbondedException
    NonbondedExceptionType
