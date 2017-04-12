The :mod:`parmed.gromacs` package
==================================

The :mod:`gromacs <parmed.gromacs>` package contains classes that can parse the
Gromacs topology and coordinate files. Like `grompp`, ParmEd pre-processes the
topology file, automatically finding and parsing the include topology files
(ITP) referenced by your topology file.  It also recognizes ``#define`` tokens,
which can be used in ParmEd the same way you can with Gromacs.

Note, like with the standard C and Gromacs preprocessors, include files will
*first* be looked for in the directory in which the processed file resides. If
it is not found there, ParmEd will look in the Gromacs installation directory.

Setting up Gromacs
------------------

Most Gromacs topology files (particularly those *written* by Gromacs) include
so-called "include topology files" which are located in the
``share/gromacs/top`` directory of the Gromacs installation. As a result, you
need to have `Gromacs <http://www.gromacs.org>`_ installed and visible to
ParmEd in many cases.

ParmEd will look for Gromacs in the following places (in this order):

1. The ``GMXDATA`` environment variable, which should point to ``share/gromacs``
   in the installation directory.
2. The ``GMXBIN`` environment variable, which should point to the location where
   the Gromacs program(s) are installed.
3. It will look for the ``share/gromacs`` directory in the following locations:
   ``/usr/local``, ``/usr``, ``/opt/local``, and ``/opt``.
4. It will look for the ``gmx`` program (the main program for Gromacs versions 5
   and higher), and assume that it is installed in the ``bin`` directory of
   the Gromacs install location.
5. If none of the above work, the default install location of
   ``/usr/local/gromacs/share/gromacs`` is used.

If you wish to specify the location, you can do so by modifying the
``GROMACS_TOPDIR`` variable in the ``parmed.gromacs`` package in the
following way::

    from parmed import gromacs
    gromacs.GROMACS_TOPDIR = "/path/to/gromacs/installation/share/gromacs/top"

Note, because strings are immutable in Python, the following will *not* work::

    from parmed.gromacs import GROMACS_TOPDIR
    GROMACS_TOPDIR = "/path/to/gromacs/installation/share/gromacs/top"

In the above example, ``parmed.gromacs.GROMACS_TOPDIR`` will remain
unchanged.

The topology file
-----------------

The primary file in Gromacs defining the system topology and the parameters
defining the force field for that system is called the *topology* file, whose
format is detailed in the Gromacs User's manual.

.. currentmodule:: parmed.gromacs
.. autosummary::
    :toctree: gromacsobj/

    GromacsTopologyFile
    GromacsGroFile

The :class:`GromacsTopologyFile` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class inherits from the base class :class:`Structure
<parmed.structure.Structure>`, adding the attributes ``parameterset`` (a
:class:`ParameterSet <parmed.parameters.ParameterSet>` object with all
parameters parsed from the include topology files) and ``itps`` (a list of all
include topology files that were included and parsed while processing the
topology file).

An example using this class on `one of the test files in the ParmEd repostory
<https://github.com/swails/ParmEd/blob/master/test/files/1aki.charmm27.top>`_
is shown below::

    >>> top = GromacsTopologyFile('1aki.charmm27.top')
    >>> top
    <GromacsTopologyFile 1960 atoms; 129 residues; 1984 bonds; parametrized>
    
You can manipulate the ``top`` instance any way you can a :class:`Structure
<parmed.structure.Structure>` instance.

Writing GROMACS Topology Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Molecule type definitions**

The GROMACS topology is constructed by piecing together *fragments* of a system,
where each fragment is, by standard convention, an individual *molecule*
(meaning that you can "reach" every atom in a molecule from every other atom in
that molecule by traversing some number of bonds -- i.e., it is a connected
graph).

These "molecules" are defined in separate ``[ moleculetype ]`` sections of the
topology file, and then the ``[ molecules ]`` section gives the order, and
number, of each molecule *type* in the system.

GROMACS can make use of the definitions of a molecule *type*, and this
distinction is particularly important for free energy calculations. For this
reason, the ``combine`` keyword in :meth:`GromacsTopologyFile.write` allows
fine-grained control over how ``[ moleculetype ]`` sections are defined. By
default, (when ``combine`` is assigned ``None``), the structure is split into
molecules and written as separate ``[ moleculetype ]`` sections. When
``combine`` is assigned the string ``all``, the entire structure is considered
to be a single molecule and is stored as a single ``[ moleculetype ]`` section.

The last allowable mode is to assign lists of lists of molecule indices for
``combine``, telling :class:`GromacsTopologyFile` exactly which molecules to
"combine" into larger ``[ moleculetype ]`` sections. There is one major
restriction here: combined molecules *must* be contiguously ordered in the
original atom sequence of the structure! The reason for this is that if the
original order of the atoms is changed, then the topology atom order will no
longer match that from the coordinate file! An annotated example may help
demonstrate this. Consider the following topology file (it is actually the
``topol3.top`` file from the ``12.DPPC/`` directory in the test suite):

::

    #include "DPPC_2.itp"
    
    ; System specifications
    [ system ]
    DPPC Bilayer
    
    [ molecules ]
    ; molecule name nr.
    DPPC 4
    SOL	122
    DPPC 4
    SOL 122
    
The DPPC and SOL molecule types are defined in the ``DPPC_2.itp`` include
topology file. When this topology is read in, there will be 252 residues -- 8
DPPC residues and 244 SOL residues::

    >>> import parmed as pmd
    >>> # You can safely ignore the warning printed here
    ... system = pmd.load_file('topol3.top')
    >>> system
    <GromacsTopologyFile 1132 atoms; 252 residues; 880 bonds; parametrized>

Now suppose that we want to do something special with a DPPC-solvent pair, and
then something special with all of the DPPC residues in the second set (possibly
a second "leaflet" in a bilayer system) again with a single solvent.  So what we
want to do is combine the 4th and 5th molecules (since indices start from 0,
this would be molecules 3 and 4), as well as the five residues 127 through 131
(indices 126, 127, 128, 129, and 130). Our ``combine`` would then look like
``[[3,4],[126,127,128,129,130]]``::

    >>> system.write('test.top', combine=[[3, 4], [126, 127, 128, 129, 130]])

If we look at the bottom of the ``test.top`` file we generated, we should see
this:

::

    [ molecules ]
    ; Compound       #mols
    DPPC                 3
    system1              1
    SOL                121
    system2              1
    SOL                121
    
If we think about it, this is exactly what we want.  We left our first 3 DPPC
residues alone, then combined DPPC and SOL into a single molecule (since the
molecule is defined by more than 1 residue, it is given a generic name
``system#``, where ``#`` increments for each multi-residue molecule that is
encountered. We left the last 121 of 122 SOL residues alone, then combined all 4
DPPC residues and the first of the second set of 122 SOL residues into a single
molecule type. Leaving us with ``system2`` and 121 ``SOL`` left to define our
original system.

Note, if you read a GROMACS topology file that combines multiple molecules into
individual molecule types, writing a copy back to a different topology file will
*not* preserve the original ``[ moleculetype ]`` definitions (since the ParmEd
``Structure`` object carries no memory of such decompositions and always
flattens the list of atoms and residues).  For the most part, you should not
bother with this added complication unless you *know* it is necessary.

**Parameters**

There are 3 ways that parameters can be specified in GROMACS topologies --
either through include topology files, as separate ``[ <parameter>type ]``
sections in the main topology file, or in-line with the actual topologies
themselves.

You can specify where the parameters are placed using the ``parameters`` keyword
in :meth:`GromacsTopologyFile.write`. The allowed values are ``inline``, a
string representing the name of a file, or an open file object. If the file name
is the same as the file name of the topology file itself, then the parameter
sections will be written to the beginning of the topology file. If it is a
different file, the name of the file will be included in the generated topology
file. The default value is ``inline``.

The GRO coordinate
------------------

ParmEd also supplies classes for parsing the various variants of the GROMACS GRO
file. Because the GRO file resembles a PDB file, the :class:`GromacsGroFile
<parmed.gromacs.GromacsGroFile>` has a ``parse`` method that will parse a GRO
file and instantiate a :class:`Structure <parmed.structure.Structure>`.

The coordinates can be extracted from the ``xx``, ``xy``, and ``xz`` attributes
on each :class:`Atom <parmed.topologyobjects.Atom>` instance in the
Structure.

The GRO file is automatically detected by :meth:`parmed.load_file`, and this
is the recommended way to parse these files.

The corresponding classes are:

.. currentmodule:: parmed.gromacs
.. autosummary::
    :toctree: amberobj/

    GromacsTopologyFile
    GromacsGroFile

Example usage
-------------

Many exciting possibilities are available with :class:`GromacsTopologyFile` and
:class:`GromacsGroFile`.  For example, ParmEd is the first software package to
enable conversion of Gromacs simulation input files to AMBER format in just a
few lines of code::

    >>> from parmed import gromacs, amber, unit as u
    >>> gmx_top = GromacsTopologyFile('topol.top')
    >>> gmx_gro = GromacsGroFile.parse('conf.gro')
    >>> gmx_top.box = gmx_gro.box # Needed because .prmtop contains box info
    >>> gmx_top.positions = gmx_gro.positions
    >>> amb_prm = AmberParm.from_structure(gmx_top)
    >>> amb_prm.write("prmtop")
    >>> amb_inpcrd = amber.AmberAsciiRestart("inpcrd", mode="w")
    >>> amb_inpcrd.coordinates = gmx_top.coordinates
    >>> amb_inpcrd.box = gmx_top.box
    >>> amb_inpcrd.close()

Furthermore, you may check the correctness of topology loading by using the
:class:`GromacsTopologyFile` object to calculate an OpenMM potential energy and
force, then comparing that result with your own output from Gromacs::

    >>> import simtk.openmm as mm
    >>> import simtk.openmm.app as app
    >>> system = gmx_top.createSystem() # OpenMM System creation
    >>> integ = mm.VerletIntegrator(1.0*u.femtosecond)
    >>> plat = mm.Platform.getPlatformByName('Reference')
    >>> simul = app.Simulation(gmx_top.topology, system, integ, plat)
    >>> simul.context.setPositions(gmx_top.positions)
    >>> simul.context.applyConstraints(1e-12)
    >>> state = simul.context.getState(getEnergy=True, getForces=True)
    >>> print(state.getPotentialEnergy())

A typical protein/water system with 23,000 atoms at ambient conditions with
periodic boundary conditions and PME electrostatics has average forces on the
order of 20 kcal/mol/Angstrom. ParmEd allows us to run this same simulation in
OpenMM or AMBER with a RMS force difference of 0.002 kcal/mol/Angstrom, i.e.
the forces between the software packages are accurate to 1 part in 10,000. The
remaining differences are due to how the different software packages treat
nonbonded interactions in the cut-off region, use of single precision in the
computation, and other small factors that are not expected to affect the
simulation results.
