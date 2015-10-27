The :mod:`parmed.amber` package
===============================

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

.. currentmodule:: parmed.amber
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
:class:`Structure <parmed.structure.Structure>`, for the *prmtop* parsers.

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
<parmed.structure.Structure>` and :class:`AmberFormat`. (The special classes
:class:`ChamberParm` and :class:`AmoebaParm`, implementing the CHARMM and Amoeba
force fields, respectively, both inherit from :class:`AmberParm`.) You need to
make sure that you use the correct class for the variant of *prmtop* file you
are using.  If you are unsure, the :func:`LoadParm` function will automatically
determine the correct class and return an instance of it (see below).

The :class:`AmberParm` class parses a *prmtop* object and stores both the raw
data (see the :class:`AmberFormat` description, above), as well as the standard
topological data structures used in the :class:`Structure
<parmed.structure.Structure>` class. Different subclasses of
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
lists inherited from :class:`Structure <parmed.structure.Structure>`.

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
    - Use an available class from :mod:`parmed.tools` from the ParmEd API when
      available, as those classes ensure that the raw data and topology objects
      are always synchronized when the action completes.
    - Feel free to use the :class:`Action <parmed.tools.actions.Action>`
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
    <class 'parmed.amber._amberparm.AmberParm'>
    >>> type(dhfrparm) # see that this is a ChamberParm
    <class 'parmed.amber._chamberparm.ChamberParm'>

You may ask "why don't I always use :func:`LoadParm`?"  And the answer is, "you
could."  But perhaps you only want to accept a certain flavor of prmtop
file---the default constructors provide a more reliable error mechanism.
:func:`LoadParm` is also marginally slower (but not much!), because it first
instantiates an :class:`AmberFormat`, checks the available flags, and then takes
a view of that class in the format of the correct :class:`AmberParm` variant.
Note that no data is explicitly *copied* during this conversion, though.

The coordinate and trajectory files
-----------------------------------

ParmEd also supplies classes for parsing all of the various formats of Amber
coordinate files, including both restarts and trajectories (in both ASCII and
NetCDF format). If ``numpy`` is available, all coordinates, velocities, and/or
forces are added to ``numpy`` arrays. If not, ``array.array`` is the return
object. These classes have *very* different APIs, so the onus is on you to
only use these arrays as ``numpy`` arrays if you are sure ``numpy`` is
available. (Since all NetCDF-capable libraries are built on ``numpy``, this
warning really only applies to ASCII files... see below).

Note that the :mod:`parmed.amber` package needs some help from third-party
packages to read NetCDF-format files.  For more information on NetCDF support,
see the last section of this page.

The corresponding classes are:

.. currentmodule:: parmed.amber
.. autosummary::
    :toctree: amberobj/

    Rst7
    AmberAsciiRestart
    NetCDFRestart
    NetCDFTraj
    AmberMdcrd
    netcdffiles.use

Inpcrd and Restart files
~~~~~~~~~~~~~~~~~~~~~~~~

The main class that deals with this is the :class:`Rst7` class.  This class will
automatically detect the file format of any inpcrd/restart file, and whether
velocities and/or unit cell information is present.  There are 2 ways you can
instantiate a Rst7.

To read a restart file, provide a filename to the :class:`Rst7` constructor. If
you plan on setting the data in order to *write* a restart file, you should
provide the argument ``natom`` to the :class:`Rst7` constructor. Examples using
the NetCDF restart file ``ncinpcrd.rst7`` from the ParmEd unit tests are shown
below (you are, of course, free to use an ASCII restart file instead)::

    >>> # Read an Amber NetCDF restart file
    >>> ncrst = Rst7("ncinpcrd.rst7")
    >>> ncrst.natom
    2101
    >>> ncrst.coordinates[:3]
    array([ 6.82122493,  6.62762507, -8.51669   ])
    >>> ncrst.vels[:3]
    array([-2.87565556, -3.02095312,  3.73882771])
    >>> ncrst.box
    array([  30.2642725,   30.2642725,   30.2642725,  109.471219 ,
            109.471219 ,  109.471219 ])

For OpenMM, coordinates and velocities need to have a specific shape as well as
units (since the OpenMM core unit system is different than the one used by
Amber). To facilitate this, the ``positions``, ``velocities``, and
``box_vectors`` attributes of the :class:`Rst7` class satisfy this requirement::

    >>> ncrst.positions[0]
    Quantity(value=array([ 6.82122493,  6.62762507, -8.51669   ]), unit=angstrom)
    >>> ncrst.velocities[0]
    Quantity(value=array([-2.87565556, -3.02095312,  3.73882771]), unit=angstrom/picosecond)
    >>> print('\n'.join(repr(x) for x in ncrst.box_vectors)) # Prettier printing
    Quantity(value=[30.264272500000001, 0.0, 0.0], unit=angstrom)
    Quantity(value=[-10.088090019353215, 28.533430037688813, 0.0], unit=angstrom)
    Quantity(value=[-10.088090019353215, -13.45078642114427, 25.164120774794483], unit=angstrom)

The :class:`Rst7` class uses the :class:`AmberAsciiRestart` and
:class:`NetCDFRestart` classes internally with a much simpler interface. Unless
you need the flexibility afforded by the raw classes, you are encouraged to use
:class:`Rst7` instead.

Writing a restart file is easy as well.  If you provide velocities, velocities
will be written.  If you supply unit cell dimensions, unit cell dimensions will
be written.  For example::

    >>> import numpy as np
    >>> new_rst7 = Rst7(natom=10)
    >>> new_rst7.coordinates = np.random.random(30)
    >>> new_rst7.box = [1, 1, 1, 90, 90, 90]
    >>> new_rst7.write('new_rst7.ncrst', netcdf=True) # Write NetCDF format
    >>> new_rst7.write('new_rst7.rst7')               # Write ASCII format

After this short script, you will have two coordinate files, one in NetCDF
format (``new_rst7.ncrst``) and one in ASCII format (``new_rst7.rst7``) with the
same (fake) coordinates. You can use ``ncdump`` to print the contents of the
NetCDF file to compare, if you wish.

Trajectory files
~~~~~~~~~~~~~~~~

Amber generates two different formats of trajectory files: *NetCDF* trajectory
files and a standardized ASCII-format (raw text) trajectory files. The
:mod:`parmed` package contains classes that support both formats, described
below:

.. currentmodule:: parmed.amber
.. autosummary::
    :toctree: amberobj/

    AmberMdcrd
    NetCDFTraj

First we will discuss the ASCII trajectory file (:class:`AmberMdcrd`). We will
first demonstrate how to *parse* an existing mdcrd file, then we will
demonstrate how to *write* a new one.

You must specify---for both an existing mdcrd and one that you wish to
create---the total number of atoms specified by the system and whether or not
the trajectory will contain unit cell dimensions.  The :meth:`coordinates
<AmberMdcrd.coordinates>` method, which takes as an optional argument a frame
number, returns either all or just the requested coordinates. The analogous
:meth:`box <AmberMdcrd.box>` method returns the unit cell lengths. The example
below uses the file ``tz2.truncoct.crd`` from the ParmEd unit tests, which has
5827 atoms and box lengths::

    >>> from parmed.amber import AmberMdcrd
    >>> crd = AmberMdcrd('tz2.truncoct.crd', natom=5827, hasbox=True, mode='r')
    >>> crd.coordinates(0)
    array([ 0.078,  3.174, -8.844, ...,  2.178, -5.224,  8.991])
    >>> len(crd.coordinates()) # number of total frames
    10
    >>> crd.coordinates(0).shape
    (17481,)
    >>> crd.coordinates(0).reshape((5827, 3)) # It is a numpy array
    array([[ 0.078,  3.174, -8.844],
           [ 0.999,  3.431, -8.521],
           [-0.183,  2.306, -8.398],
           ..., 
           [ 1.761, -6.031,  8.689],
           [ 2.047, -6.123,  7.781],
           [ 2.178, -5.224,  8.991]])
    >>> crd.box(0)
    array([ 42.439,  42.439,  42.439])
    >>> len(crd.box()) # Should also be the number of frames
    10

When you are writing an :class:`AmberMdcrd` file, you need to use the methods
:meth:`add_coordinates <AmberMdcrd.add_coordinates>` and :meth:`add_box
<AmberMdcrd.add_box>` instead of :meth:`coordinates <AmberMdcrd.coordinates>`
and :meth:`box <AmberMdcrd.box>`, as shown in the example below::

    import numpy as np
    from parmed.amber import AmberMdcrd
    crd = AmberMdcrd('new.mdcrd', natom=30, hasbox=True, mode='w')
    # Add the first frame
    crd.add_coordinates(np.random.random((30,3)))
    crd.add_box([1, 1, 1])
    # Add the second frame
    crd.add_coordinates(np.random.random((30,3)))
    crd.add_box([1.1, 1.1, 1.1])
    # Close the file
    crd.close()

You should be able to open the file ``new.mdcrd`` and see 90 numbers spread over
9 lines (10 numbers per line), then a line with the three box dimensions,
followed by 9 more lines with 90 total numbers and another line with box
dimensions.

Now I will demonstrate the use of the :class:`NetCDFTraj` class, which is
designed to be very similar to the :class:`AmberMdcrd` class. Creating an
instance is different, but querying and adding coordinates (and unit cell
dimensions, velocities, and forces) is done in the same way.

To parse an existing NetCDF trajectory, use the :meth:`open_old
<NetCDFTraj.open_old>` constructor. To write a new one, use :meth:`open_new
<NetCDFTraj.open_new>`. An example parsing an existing NetCDF trajectory (using
the ``tz2.truncoct.nc`` file in the ParmEd unit tests) is shown below::

    >>> from parmed.amber import NetCDFTraj
    >>> traj = NetCDFTraj.open_old('tz2.truncoct.nc')
    >>> traj.atom  # number of atoms in the file
    5827
    >>> traj.frame # number of frames currently in the file
    10
    >>> traj.coordinates(0)
    array([ 0.07762779,  3.1744082 , -8.84385777, ...,  2.17817903,
           -5.22395372,  8.99102497], dtype=float32)
    >>> traj.coordinates(9)
    array([ 0.33575553,  4.19061375, -9.66586971, ...,  2.48203635,
           -3.30589318,  8.45326519], dtype=float32)
    >>> traj.box(0)
    array([  42.43884853,   42.43884853,   42.43884853,  109.471219  ,
            109.471219  ,  109.471219  ])
    >>> traj.box(9)
    array([  42.42843203,   42.42843203,   42.42843203,  109.471219  ,
            109.471219  ,  109.471219  ])

One thing you'll notice is that the :class:`NetCDFTraj` class has the
:attr:`atom <NetCDFTraj.atom>` and :attr:`frame <NetCDFTraj.frame>` attributes
that indicate how many atoms are in the trajectory and how many frames are
present. Furthermore, unlike the :class:`AmberMdcrd` format---which stores only
the unit cell lengths---the :class:`NetCDFTraj` stores both unit cell lengths
*and* angles.

Furthermore, velocities and forces can be saved in NetCDF trajectory files,
allowing you to extract them with the analogous :meth:`velocities
<NetCDFTraj.velocities>` and :meth:`forces <NetCDFTraj.forces>` methods.
Likewise, you can set them with :meth:`add_velocities
<NetCDFTraj.add_velocities>` and :meth:`add_forces <NetCDFTraj.add_forces>`.

An example of writing a NetCDF trajectory (with velocities and forces) is shown
below::

    from parmed.amber import NetCDFTraj
    import numpy as np
    # Open the NetCDF trajectory file for writing
    crd = NetCDFTraj.open_new('new_crd.nc', natom=30, box=True,
                              crds=True, vels=True, frcs=True)
    # Add random coordinates, velocities, and forces for first frame
    crd.add_coordinates(np.random.random((30,3)))
    crd.add_velocities(np.random.random((30,3)) * 10)
    crd.add_forces(np.random.random((30,3)) / 10)
    crd.add_box([1, 1, 1, 90, 90, 90])
    # Add random second frame
    crd.add_coordinates(np.random.random((30,3)))
    crd.add_velocities(np.random.random((30,3)) * 10)
    crd.add_forces(np.random.random((30,3)) / 10)
    crd.add_box([1.1, 1.1, 1.1, 90, 90, 90])
    # Close the file
    crd.close()

As an exercise, try using :meth:`NetCDFTraj.open_old` to check the contents of
``new_crd.nc`` and verify that it contains the information you expect.

Reading and writing NetCDF files
--------------------------------

NetCDF---standing for *Net*\ work *C*\ ommon *D*\ ata *F*\ ormat---is a binary
file format that is flexible, portable (i.e., between `big endian and little
endian <http://en.wikipedia.org/wiki/Endianness>`_), and extensible (i.e., you
can extend a file convention to include more data *without* breaking existing
parsers).

Python cannot read NetCDF files without the support of an external package. In
particular, the following packages provide Python bindings to the NetCDF API:

    1. `scipy <http://www.scipy.org/>`_
    2. `netCDF4 <https://github.com/Unidata/netcdf4-python>`_
    3. `ScientificPython <http://dirac.cnrs-orleans.fr/plone/software/scientificpython/>`_

All three of these files provide similar, although not identical interfaces. The
ParmEd API supports all three NetCDF backends. You can swap out backends at any
time using the :func:`use <parmed.amber.netcdffiles.use>` function. If you do
not choose one, one will be picked automatically based on the available
packages.  If several of those packages are available, the first one in the
above list that is available will be picked.

For most users, simply ensuring that at least one of the above packages is
available and working and using the default will suffice. For the other small
percentage that may want more fine-grained control, the :func:`use
<parmed.amber.netcdffiles.use>` function will let you specify a particular
implementation to use.  You can change the back-end at any time, but note that
the behavior of a NetCDF object *created* with one implementation is undefined
if you switch to a different one. So you are cautioned against switching NetCDF
implementations when open NetCDF file handles are still open.

An example of switching between NetCDF implementations is shown below::

    from parmed.amber import netcdffiles
    # Use the default NetCDF package (if one is already selected,
    # this does nothing)
    netcdffiles.use() # You can also pass "None"
    # Use scipy for NetCDF files
    netcdffiles.use('scipy')
    # Use ScientificPython for NetCDF files
    netcdffiles.use('Scientific') # Can also be ScientificPython
    # Use netCDF4
    netcdffiles.use('netCDF4')

Amber mask syntax
-----------------

The Amber programs use a selection syntax very similar to Chimera's selection
syntax, referred to as an *Amber mask*. The character ``:`` is used for residue
selections and ``@`` is used for atom selections.  Both numbers and names are
supported in the atom and residue selections (where the serial numbers range
from 1 to *N*, where *N* is either the total number of atoms or residues,
depending on what is being selected).

The atom/residue name or number list can be specified as a comma-delimited list
of entries (*inclusive* number ranges are also supported).  You can also use
logical operators ``&`` (and), ``|`` (or), and ``!`` (not) and group expressions
with parenthesis to denote order-of-operations.

Residue selections
~~~~~~~~~~~~~~~~~~

Examples of residue selections, accompanied by comments annotating what is being
selected, is shown below::

    :1,3,6-10    Select residues 1, 3, 6, 7, 8, 9, and 10
    :ALA,LYS,10  Select all residues with names ALA or LYS, and residue 10
    :1-3,6-10    Select residues 1, 2, 3, 6, 7, 8, 9, and 10

Atom selections
~~~~~~~~~~~~~~~

Examples of residue selections, accompanied by comments annotating what is being
selected, is shown below::

    @1-100,150-160   Select the first 100 atoms as well as the 11 atoms 150-160
    @CA,C,10-15      Select all atoms whose names are CA or C, and 6 atoms 10-15

Atom element selections
~~~~~~~~~~~~~~~~~~~~~~~

You can select all residues that are a specific *element* if you want to.
Examples are shown below, accompanied by comments annotating what is being
selected::

    @/C          Select all carbon atoms
    @/H          Select all hydrogen atoms
    @/C,N,O,H,S  Select all carbon, nitrogen, oxygen, hydrogen, and sulfur atoms
    @/Fe         Select all iron atoms

Be careful, though!  These element names are case-sensitive, and bad element
names will cause the mask parser to raise an exception.

Atom type selections
~~~~~~~~~~~~~~~~~~~~

In some cases, it is more succinct to select by atom *type* name, rather than
atom name. Be careful though! Atom type names can change from force field to
force field, and are considered an implementation detail.  While atom *names*
should never change (and are based on the PDB standard, typically), atom *type*
names can be set and/or changed to whatever the force field developer wants.

You have been warned. If you wish to push on, though, you can use the ``@%``
characters to select from type names rather than atom names. Note, no indices
(or index ranges) can be used, as atom types have no defined order.  Examples
are shown below::

    @%CT,CX         Select all atoms whose types are CT or CX

Logical operators
~~~~~~~~~~~~~~~~~

You often want to select certain atoms from certain residues. You can use the
``&`` binary operator to do this, as shown in examples below::

    (:1-100)&(@CA)      Selects atoms named CA in the first 100 residues
    (:1-100)&(@CA,CB)   Selects atoms named CA or CB in the first 100 residues
    (:ASP,GLU)&(@CA)    Selects atoms named CA in all ASP or GLU residues

This is so common, you can use a simple shortcut, demonstrated below, giving
equivalent expressions to the examples above::

    :1-100@CA
    :1-100@CA,CB
    :ASP,GLU@CA

You can also expand your selection using the ``|`` (or) operator. For example::

    :1-100|@CA      Selects all atoms in the first 100 residues AND all atoms named CA

Don't get confused between what the ``&`` (and) and ``|`` (or) operators do.
The above example selects all atoms if it satisfies the criteria of being in the
first 100 residues *or* if it satisfies the criteria of having the name ``CA``,
meaning it selects all atoms in the first 100 residues *and* all atoms named
``CA``.

You can also select all atoms *except* a particular selection::

    !(:LIG)         Selects all atoms that are NOT in a residue named LIG

Wild-cards
~~~~~~~~~~

Let's suppose you want to select all hydrogen atoms. These atoms usually have
many names, but they almost always start with "H", so we would like a way to
select all atoms that start with the letter "H".  We can do this with
*wild-cards*, which is the ``=`` character in Amber masks, and works for both
atoms and residues. For example::

    :AS=@H=   Select all atoms whose residue starts with AS AND atom name starts with H

Distance-based selections
~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes you also want to be able to select all atoms within a certain distance
(or all atoms within a *residue* that is a certain distance) from a specific set
of atoms.  For this, the ``<`` and ``>`` operators are used to denote distances
that are less than or greater than a particular value, respectively.

This is probably the most complex part of the Amber mask syntax, but it is also
quite powerful.  In addition to the ``<`` and ``>`` operators, you must use
either ``@`` or ``:`` to indicate whether you want to select individual *atoms*
or whole *residues* satisfying the distance.  When selecting residues, the whole
residue is selected if *any atom* satisfies the distance criteria.  The general
syntax is ``SELECTION<:DISTANCE``, where ``SELECTION`` is an Amber mask (``:``
can also be ``@`` to select only atoms), and ``DISTANCE`` is a floating point
distance that defines the criteria.  Examples are shown below with annotations::

    :8<@5.0     All atoms that are LESS THAN 5 Angstroms away from any atom in residue 8
    :1>:10.0    All atoms in any residue in which one atom is MORE THAN 10 angstroms away from residue 1
