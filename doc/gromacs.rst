The :mod:`chemistry.gromacs` package
==================================

The :mod:`gromacs` package contains classes that can parse the Gromacs topology
and coordinate files. Like `grompp`, ParmEd pre-processes the the topology file,
automatically finding and parsing the include topology files (ITP) referenced by
your topology file.  It also recognizes ``#define`` tokens, which can be used in
ParmEd the same way you can with Gromacs.

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
``GROMACS_TOPDIR`` variable in the ``chemistry.gromacs`` package in the
following way::

    from chemistry import gromacs
    gromacs.GROMACS_TOPDIR = "/path/to/gromacs/installation/share/gromacs/top"

Note, because strings are immutable in Python, the following will *not* work::

    from chemistry.gromacs import GROMACS_TOPDIR
    GROMACS_TOPDIR = "/path/to/gromacs/installation/share/gromacs/top"

In the above example, ``chemistry.gromacs.GROMACS_TOPDIR`` will remain
unchanged.

The topology file
-----------------

The primary file in Gromacs defining the system topology and the parameters
defining the force field for that system is called the *topology* file, whose
format is detailed in the Gromacs User's manual.

.. currentmodule:: chemistry.gromacs
.. autosummary::
    :toctree: gromacsobj/

    GromacsTopologyFile
    GromacsGroFile

The :class:`GromacsTopologyFile` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class inherits from the base class :class:`Structure
<chemistry.structure.Structure>`, adding the attributes ``parameterset`` (a
:class:`ParameterSet <chemistry.parameters.ParameterSet>` object with all
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
<chemistry.structure.Structure>` instance.

The GRO coordinate
------------------

ParmEd also supplies classes for parsing the various variants of the GROMACS GRO
file. Because the GRO file resembles a PDB file, the :class:`GromacsGroFile
<chemistry.gromacs.GromacsGroFile>` has a ``parse`` method that will parse a GRO
file and instantiate a :class:`Structure <chemistry.structure.Structure>`.

The coordinates can be extracted from the ``xx``, ``xy``, and ``xz`` attributes
on each :class:`Atom <chemistry.topologyobjects.Atom>` instance in the
Structure.

The GRO file is automatically detected by :meth:`chemistry.load_file`, and this
is the recommended way to parse these files.

The corresponding classes are:

.. currentmodule:: chemistry.gromacs
.. autosummary::
    :toctree: amberobj/

    GromacsTopologyFile
    GromacsGroFile

Example usage
-------------

Many exciting possibilities are available with `GromacsTopologyFile` and `GromacsGroFile`.
For example, ParmEd is the first software package to enable conversion of Gromacs simulation
input files to AMBER format in just a few lines of code::

    >>> from chemistry import gromacs, amber
    >>> gmx_top = GromacsTopologyFile('topol.top')
    >>> gmx_gro = GromacsGroFile('conf.gro')
    >>> gmx_top.box = gmx_gro.box # Needed because .prmtop contains box info
    >>> gmx_top.positions = gmx_gro.positions
    >>> amb_prm = AmberParm.load_from_structure(gmx_top)
    >>> amb_prm.write("prmtop")
    >>> amb_inpcrd = amber.AmberAsciiRestart("inpcrd", mode="w")
    >>> amb_inpcrd.coordinates = np.array(gmx_top.positions.value_in_unit(u.angstrom)).reshape(-1,3)
    >>> amb_inpcrd.box = gmx_top.box
    >>> amb_inpcrd.close()

Furthermore, you may check the correctness of topology loading by using the `GromacsTopologyFile`
object to calculate an OpenMM potential energy and force, then comparing that result with your own 
output from Gromacs. 

A typical protein/water system with 23,000 atoms at ambient conditions with periodic boundary conditions
and PME electrostatics has average forces on the order of 20 kcal/mol/Angstrom. ParmEd allows us to run 
this same simulation in OpenMM or AMBER with a RMS force difference of 0.002 kcal/mol/Angstrom, i.e.
the forces between the software packages are accurate to 1 part in 10,000. The remaining differences 
are due to how the different software packages treat nonbonded interactions in the cut-off region, use
of single precision in the computation, and other small factors that are not expected to affect the 
simulation results.

    >>> import simtk.openmm as mm
    >>> import simtk.openmm.app as app
    >>> system = gmx_top.createSystem() # OpenMM system creation; make sure to pass keyword arguments consistent with your .mdp file
    >>> integ = mm.VerletIntegrator(1.0*u.femtosecond)
    >>> plat = mm.Platform.getPlatformByName('Reference')
    >>> simul = app.Simulation(gmx_top.topology, system, integ, plat)
    >>> simul.context.setPositions(gmx_top.positions)
    >>> simul.context.applyConstraints(1e-12)
    >>> state = simul.context.getState(getEnergy=True, getForces=True)
    >>> print(state.getPotentialEnergy())

