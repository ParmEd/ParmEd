OpenMM Functionality
====================

Since ParmEd is a library for reading, converting, and modifying full molecular
mechanical descriptions of chemical systems in a wide variety of different
families of force fields, supporting molecular simulation directly using the
fantastic `OpenMM Python API <http://simtk.org/home/openmm>`_ was a natural
extension.

This page is not meant as an exhaustive description of the OpenMM library and
its usage. Instead, you should visit the `OpenMM website
<http://simtk.org/home/openmm>`_ for that. However, this page *will* provide a
brief description of OpenMM, what I consider to be its strengths that can make
it an invaluable tool in the field of molecular mechanics. It will also present
an introduction to using OpenMM with the tools provided by ParmEd through a
handful of examples.

What is OpenMM
--------------

OpenMM is not a *program* in the traditional sense. There is no *OpenMM* program
that you can run from the terminal like you can with AMBER, CHARMM, or GROMACS,
for example. Instead, it is a library of routines written in C++ that allow you
to program your own molecular models in C, C++, Fortran, or, as we will
demonstrate here, Python. Its basic features include:

    * Basic valence-term forces, like bonds, angles, torsions, and
      coupled-torsion correction maps.
    * Nonbonded potential terms, like the typical electrostatic and
      Lennard-Jones potentials, Generalized Born implicit solvent models, and
      more, with accurate long-range potentials.
    * Various integrators to carry out molecular dynamics, like the Verlet
      integrator for *traditional* MD, or stochastic integrators like those
      for Langevin and Brownian dynamics.
    * Thermostats and barostats for sampling from various statistical ensembles.

What makes OpenMM so awesome?
-----------------------------

There are two key features that not only make *OpenMM* awesome, but make it
unlike any other molecular dynamics software in existence.

Stellar performance on GPUs
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Among the primarily advertised selling-points of OpenMM is the computational
performance that is possible by running OpenMM on NVidia and AMD/ATI graphics
processing units (GPUs). OpenMM utilizes both CUDA and OpenCL to program
molecular models for GPUs, and runs the entire calculation on the accelerator
device to eliminate the GPU-to-CPU communications that can significantly limit
performance.

In general, the CUDA platform is faster than the OpenCL platform, but OpenCL
works on a much wider array of hardware. I will not expound on the details of
the GPU performance, as you are directed to look at the various GPU-related
publications for such comparisons.  Suffice it to say that, currently (February,
2015), GPU performance is ca. one, possibly up to two orders of magnitude
*faster* than running on every core of a standard multicore server processor (8
to 16 cores) with a multithreaded MD application.

Custom Forces
~~~~~~~~~~~~~

In my opinion, what makes OpenMM truly awesome is its custom force capabilities.
OpenMM allows you to implement a new potential for bonds, angles, torsions, and
even nonbonded pairwise and multi-particle interactions simply by providing an
equation for the potential energy.

OpenMM will then analytically differentiate this analytical function to get an
analytical, closed-form expression for the gradients (forces), and build a
modestly optimized GPU implementation of that term!

For example, suppose you want to model a few bonds in your biomolecule using a
typical Morse potential, the OpenMM Python code that will do this is shown
below::

    force = CustomBondForce("D*(1-exp(-alpha*(r-req)))^2")

And that's it! Once you add the bonds to that force, with their parameters,
OpenMM will build an optimized kernel for the calculation. Not only does this
eliminate the development cost of building and optimizing a GPU implementation
of a new potential, it also provides a pre-tested, and pre-debugged version of
that kernel.

How does ParmEd enhance OpenMM?
-------------------------------

.. currentmodule:: parmed.openmm
.. autosummary::
    :toctree: openmmobj/

    StateDataReporter
    NetCDFReporter
    MdcrdReporter
    RestartReporter
    ProgressReporter
    XmlFile
    load_topology
    energy_decomposition

ParmEd provides a common framework for constructing and representing fully
parametrized force field models for various systems instantiated from a wide
range of file formats. In particular, it supports the following file formats:

    * Standard Amber topology file
    * Amber Chamber-style topology file*
    * Amber AMOEBA-style topology file*
    * CHARMM PSF file
    * Gromacs topology file

The *-marked options are *not* available in the OpenMM Python application layer.
While the OpenMM Python application layer *does* support standard Amber topology
files, it does not support either old-style topology files or those that define
a 10-12 nonbonded potentials between certain pairs of atoms. ParmEd supports
both of these. With Gromacs topology files, the OpenMM Python application layer
does not correctly handle virtual sites, while the Gromacs topology file parser
in ParmEd does.

In addition to these file formats, ParmEd also supports several new reporter
classes in addition to the small number provided by ParmEd:

    * :class:`StateDataReporter` -- This takes extra arguments specifying the
      units of each of the types of data (like the energy, volume, and time
      units). The defaults correspond to the AKMA unit system, which is more
      familiar to Amber and CHARMM users.
    * :class:`NetCDFReporter` -- This allows you to write a trajectory file in
      the Amber NetCDF format.
    * :class:`MdcrdReporter` -- This allows you to write a trajectory file in
      the Amber ASCII mdcrd format.
    * :class:`RestartReporter` -- This allows you to periodically write NetCDF
      or ASCII restart files during the course of the calculation
    * :class:`ProgressReporter` -- This prints a file during the course of the
      simulation tracking the runtime speed of the calculation and predicting
      the amount of time remaining.

The :func:`energy_decomposition` function takes as input a :class:`Structure
<parmed.structure.Structure>` instance, OpenMM ``Context``, and an optional
energy unit (``nrg``) and returns a dictionary of all energy components for the
different force groups. This permits an form of energy decomposition that allows
energy components to be compared between programs more effectively. For
example::

    >>> import parmed as pmd
    >>> from simtk.openmm import app
    >>> from simtk import openmm as mm
    >>> # Instantiate the parm and create the system
    ... parm = pmd.load_file('tip4p.parm7', 'tip4p.rst7')
    >>> system = parm.createSystem(nonbondedMethod=app.PME,
    ...                            nonbondedCutoff=8*pmd.unit.angstrom)
    >>> # Make the context and set the positions
    ... context = mm.Context(system, mm.VerletIntegrator(0.001))
    >>> context.setPositions(parm.positions)
    >>> # Find the energy decomposition
    ... pmd.openmm.energy_decomposition(parm, context)
    {'total': -2133.295388974015, 'nonbonded': -2133.2953890231834, 'bond': 1.0126508125518531e-07}

Loading OpenMM Objects
----------------------

The OpenMM ``Topology`` object and ``System`` object contain the information
stored in :class:`Structure <parmed.structure.Structure>`. You can use
:func:`load_topology` to load an OpenMM ``Topology`` and create a
:class:`Structure <parmed.structure.Structure>` instance from it. If you provide
either a file containing a serialized ``System`` in XML format or a ``System``
object directly, parameters will be extracted from the various forces and added
to the generated :class:`Structure <parmed.structure.Structure>`. You can also
pass coordinates (or any coordinate file, including an OpenMM XML-serialized
State file) with the ``xyz`` argument and unit cell dimensions with the ``box``
argument (which will override any unit cell information contained in the input
``Topology``, ``System`` or coordinate file if applicable).

The :class:`XmlFile` class can parse and return a deserialized object from an
OpenMM-generated XML file. Supported XML files are:

    - Serialized ``System`` (returns an OpenMM :class:`System` instance)
    - Serialized ``State`` (returns a container object with attributes
      ``coordinates``, ``velocities``, ``forces``, ``energy``, and ``time``)
    - Serialized ``Integrator`` (returns an OpenMM :class:`Integrator` subclass)
    - XML ``ForceField`` file (returns an OpenMM :class:`ForceField` instance)

Examples
--------

Right now, the two main pathways to run an OpenMM simulation starting from a
fully parametrized molecular mechanical model (click on either option to view
the annotated and explained example):

.. toctree::
    :maxdepth: 2

    Starting from AMBER prmtop and inpcrd files <omm_amber>
    Starting from CHARMM PSF and coordinate files <omm_charmm>
    Starting from Gromacs TOP and GRO files <omm_gromacs>

Taking OpenMM ``Topology`` and ``System`` to a ParmEd ``Structure``
-------------------------------------------------------------------

While the above sections described how you would generate an OpenMM ``System``
and ``Topology`` instance from any of a number of file formats (e.g., Amber
topology file, Gromacs topology file, or CHARMM PSF file), it is also possible
to go in the reverse direction. That is, to start with ``Topology`` and
``System`` instances and convert those to a :class:`Structure
<parmed.structure.Structure>` instance.

The function for this is :func:`load_openmm
<parmed.openmm.topsystem.load_openmm>`, and takes a mandatory ``Topology``
instance, along with either an optional ``System`` instance or name of a
serialized ``System`` XML file defining an OpenMM ``System``, and returns a
populated :class:`Structure <parmed.structure.Structure>` from it. This is
particularly useful when you parametrize a system using the OpenMM modelling
capabilities, but want to use that parametrized system in another program, like
Amber or Gromacs.

An example is shown below, using the OpenMM functionality to parametrize a PDB
file with the ff99SB force field::

    import parmed as chem
    import parmed.unit as u

    from simtk.openmm import app
    from simtk import openmm as mm

    pdb = app.PDBFile('input.pdb')
    forcefield = app.ForceField('amber99sb.xml', 'tip3p.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,
                                     nonbondedCutoff=1*u.nanometer)

    struct = chem.openmm.load_topology(pdb.topology, system)

There are some limitations to this functionality, itemized below:

- ParmEd does not recognize all of the different OpenMM Systems it can generate,
  such as any of those features implemented using a ``CustomNonbondedForce``
  (e.g., NBFIX, 10-12 potential, 12-6-4 potential, etc.).
- You should make sure to create the ``System`` with no constraints, since
  OpenMM may be missing bond or angle terms associated with the constrained
  degrees of freedom.
- With the exception of certain ``CustomTorsionForce`` which are recognized as
  quadratic improper torsions, the presence of any ``CustomForce`` instances
  prevents ParmEd from recognizing the potential.
