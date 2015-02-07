The :mod:`chemistry.amber` package
==================================

The :mod:`amber` package contains classes that can parse most of the file
formats used by the Amber molecular dynamics package.  In particular are Amber
parameter and topology files as well as the various formats of coordinate and
trajectory files that Amber supports (and creates).

It is out of the scope of this site to discuss the details of the file formats
themselves, and you are instead forwarded to the `official Amber website
<http://ambermd.org>`_ and `file format specifications
<http://ambermd.org/formats.html>`_\ .

The parameter-topology (*prmtop*) file
--------------------------------------

The primary file in Amber defining the system topology and the parameters
defining the force field for that system is called the *prmtop* file, whose
format is detailed on http://ambermd.org/formats.html. There are several
variants of this extensible file format, and these parsers rigorously implement
the format specification.

.. currentmodule:: chemistry.amber
.. autosummary::
    :toctree: amberobj/

    AmberFormat
    AmberParm
    ChamberParm
    AmoebaParm
    LoadParm
    AmberMask

The :class:`AmberFormat` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class simply parses the raw data from a file with the format of an Amber
*prmtop* file and makes the data available in a ``dict`` whose keys are the
``FLAG_NAME`` in the ``%FLAG <FLAG_NAME>`` lines. Some files---like the MDL
files used in 3D-RISM calculations---use this format to store raw, structured
data, and this class is appropriate for use with those files.

In most cases, however, this simply serves as a base class, along with
:class:`Structure <chemistry.structure.Structure>`, for the *prmtop* parsers.

An example using this class on `one of the test files in the ParmEd repostory
<https://github.com/swails/ParmEd/blob/master/test/files/cSPCE.mdl>`_ is shown
below::

    >>> amb = AmberFormat('cSPCE.mdl')
    >>> amb.parm_data['ATMTYP']
    [1, 2]
    >>> amb.parm_data['COORD']
    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.333314, 0.942816, 0.0]

If you want to, you can make modifications to the data and then write a new
file::

    >>> amb.parm_data['COORD'][0] = 10
    >>> amb.write_parm('new.mdl')
    >>> # Let's check that our new file has the changed COORD
    >>> new_amb = AmberFormat('new.mdl')
    >>> new_amb.parm_data['COORD']
    [10.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.333314, 0.942816, 0.0]

The :class:`AmberParm`, :class:`ChamberParm`, and :class:`AmoebaParm` classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These classes are subclasses of both :class:`Structure
<chemistry.structure.Structure>` and :class:`AmberFormat`. (The special classes
:class:`ChamberParm` and :class:`AmoebaParm`, implementing the CHARMM and Amoeba
force fields, respectively, both inherit from :class:`AmberParm`.) You need to
make sure that you use the correct class for the variant of *prmtop* file you
are using.  If you are unsure, the :func:`LoadParm` function will automatically
determine the correct class and return an instance of it (see below).

The :class:`AmberParm` class parses a *prmtop* object and stores both the raw
data (see the :class:`AmberFormat` description, above), as well as the standard
topological data structures used in the :class:`Structure
<chemistry.structure.Structure>` class. Different subclasses of
:class:`AmberParm` will populate different term and parameter lists depending on
the terms present in the force field.

The example of an Amber topology file shown here is ``trx.prmtop``, taken from
the ParmEd unit test directory `available here
<http://github.com/swails/ParmEd/tree/master/test/files/trx.prmtop>`_::

    >>> parm = AmberParm('trx.prmtop')
    >>> parm.atoms[0].name
    'N'
    >>> parm.residues[0].name
    'SER'
    >>> parm.atoms[0].residue is parm.residues[0]
    True
    >>> len(parm.residues[0])
    13
    >>> parm.atoms[12].residue.name
    'SER'
    >>> parm.atoms[13].residue.name # this should be a new residue!
    'ASP'
    >>> parm.bonds[0].atom1
    <Atom C [11]; In SER 0>
    >>> parm.bonds[0].atom2
    <Atom O [12]; In SER 0>
    >>> parm.bonds[0].type.k # Check out our parameter types
    570.0
    >>> parm.bonds[0].type.req
    1.229

WARNING
~~~~~~~

Since the raw data and topological features are inextricably linked, there are
risks associated with modifying the raw data in :attr:`AmberParm.parm_data`
*and* the parameters and properties of the topological data structures in the
lists inherited from :class:`Structure <chemistry.structure.Structure>`.

The required overhead to keep :attr:`AmberParm.parm_data` and the topology
structures synchronized at all times introduces too great a cost, so keeping
them synchronized largely falls on the programmer, although there is
functionality introduced to help:

* Changes to the size or contents of any of the lists in the structure will
  trigger :meth:`AmberParm.remake_parm` to be called inside
  :meth:`AmberParm.write_parm` to make sure that new topology files are
  written with the updated parameter definitions.
* There are a handful of synchronization routines designed to copy information
  from the parameter and topology arrays into :attr:`AmberParm.parm_data` and
  vice-versa. These are summarized in the table below. These are not called
  automatically (with the exception of :meth:`AmberParm.remake_parm` prior to
  writing a new *prmtop* file) for performance reasons.
* To modify an :class:`AmberParm` instance safely, follow these guidelines:

    - Modify the attributes on the parameter and topology objects, and call
      :attr:`AmberParm.remake_parm` if you ever need to access the contents
      of :attr:`AmberParm.parm_data` (but try to avoid using the raw data if
      possible)
    - Use an available class from :mod:`ParmedTools` from the ParmEd API when
      available, as those classes ensure that the raw data and topology objects
      are always synchronized when the action completes.
    - Feel free to use the :class:`Action <ParmedTools.ParmedActions.Action>`
      classes as examples for working with an :class:`AmberParm`.

+----------------------------------+---------------------------------------------------------------------------------------------+
| Method name                      | Synchronization functionality                                                               |
+==================================+=============================================================================================+
| :meth:`AmberParm.load_atom_info` | Loads the atom properties from ``parm_data`` to the ``atoms`` array                         |
+----------------------------------+---------------------------------------------------------------------------------------------+
| :meth:`AmberParm.load_structure` | Loads all parameter/topology data structures from the information in ``parm_data``          |
+----------------------------------+---------------------------------------------------------------------------------------------+
| :meth:`AmberParm.remake_parm`    | Flushes all changes to the parameter/topology data structures to the raw ``parm_data`` dict |
+----------------------------------+---------------------------------------------------------------------------------------------+

An Example
~~~~~~~~~~

An example of the warnings described here is illustrated below, showing how the
raw data and parameter/topology data structures can become out-of-sync and lead
to potential errors (like above, we use ``trx.prmtop`` from the ParmEd test
suite, `available here
<http://github.com/swails/ParmEd/tree/master/test/files/trx.prmtop>`_::

    >>> parm = AmberParm('trx.prmtop')
    >>> parm.ptr('NATOM')
    1654
    >>> len(parm.atoms)
    1654
    >>> parm.atoms[0].name
    'N'
    >>> parm.parm_data['ATOM_NAME'][0]
    'N'

Notice how the number of atoms according to the raw data ``POINTERS`` array
matches the total number of atoms in the structure, and that the names of the
atoms match between the raw data and :class:`Atom <topologyobjects.Atom>`
instances.  Now let's see what happens if we modify our atom name and pop our
last atom off the structure::

    >>> parm.atoms.pop()
    <Atom OXT [-1]>
    >>> len(parm.atoms)
    1653
    >>> parm.ptr('natom')
    1654
    >>> parm.atoms[0].name = 'New'
    >>> parm.parm_data['ATOM_NAME'][0]
    'N'

Here we can clearly see that the raw data is not being updated with the changes
we are making to our structure.  Now let's use :meth:`AmberParm.remake_parm` to
synchronize the raw data with the parameter/topology changes we are making::

    >>> parm.remake_parm()
    >>> len(parm.atoms)
    1653
    >>> parm.ptr('natom')
    1653
    >>> parm.atoms[0].name
    'New'
    >>> parm.parm_data['ATOM_NAME'][0]
    'New'

They match again. Whew! This goes both ways, though. Let's modify some of the
raw data, and see how it does not propagate up to the parameter/topology data
structures::

    >>> parm.parm_data['ATOM_NAME'][1]
    'H1'
    >>> parm.parm_data['ATOM_NAME'][1] = 'Hey'
    >>> parm.atoms[1].name
    'H1'

Just like with modifying the parameter/topology objects, the changes to the raw
data do not propagate. You need to use :meth:`AmberParm.load_atom_info` (if you
made modifications to arrays that were *not* atomic property arrays, you need to
use :meth:`AmberParm.load_structure`, but this is more expensive and arrays that
are not atomic properties are *much* harder to directly modify correctly)::

    >>> parm.load_atom_info()
    >>> parm.parm_data['ATOM_NAME'][1]
    'Hey'
    >>> parm.atoms[1].name
    'Hey'

Hopefully this little demonstration showed you the pitfalls of modifying an
:class:`AmberParm` instance (or one of its subclasses), and shows you how to do
it the *right* way. But be careful! These classes give you plenty of rope to
hang yourself with.

Generalizing with :func:`LoadParm`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you use the :class:`AmberParm`, :class:`ChamberParm`, or :class:`AmoebaParm`
classes directly, you need to know *which* flavor of topology file you have, or
you will get numerous errors.

In cases where you *don't* know a priori what kind of topology file to expect,
you can use :func:`LoadParm` (with exactly the same semantics as the class
constructors), but the return type will always be the right kind (based on the
flags present in the prmtop file).

For example, let's try it on ``trx.prmtop`` (an Amber prmtop file) and
``dhfr_cmap_pbc.parm7`` (a Chamber prmtop file), both available in the ParmEd
test suite::

    >>> trxparm = LoadParm('trx.prmtop')
    >>> dhfrparm = LoadParm('dhfr_cmap_pbc.parm7')
    >>> type(trxparm) # see that this is an AmberParm
    <class 'chemistry.amber._amberparm.AmberParm'>
    >>> type(dhfrparm) # see that this is a ChamberParm
    <class 'chemistry.amber._chamberparm.ChamberParm'>

You may ask "why don't I always use :func:`LoadParm`?"  And the answer is, "you
could."  But perhaps you only want to accept a certain flavor of prmtop
file---the default constructors provide a more reliable error mechanism.
:func:`LoadParm` is also marginally slower (but not much!), because it first
instantiates an :class:`AmberFormat`, checks the available flags, and then takes
a view of that class in the format of the correct :class:`AmberParm` variant.
Note that no data is explicitly *copied* during this conversion, though.

The coordinate and trajectory files
-----------------------------------

.. currentmodule:: chemistry.amber
.. autosummary::
    :toctree: amberobj/

    Rst7
    AmberAsciiRestart
    NetCDFRestart
    NetCDFTraj
    AmberMdcrd
    netcdffiles.use
