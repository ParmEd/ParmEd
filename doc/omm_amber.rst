Running OpenMM Simulation with AMBER Files
==========================================

There are popular ways to model biomolecular systems.  Because solvation effects
are often (always?) critically important to biological function, we need some
way to model the solvent.  The two popular approaches are to employ a continuum
model with a dielectric constant equal to that of the bulk solvent or to model
the solvent models directly in your system.  These two approaches are termed
_implicit_ and _explicit_ solvation, respectively.

The next 2 sections present examples using a Generalized Born implicit solvent
model and explicit solvent based on the TIP3P water model. All of the files and
examples here are included in the ``examples/amber`` directory of the ParmEd
release.

Generalized Born
----------------

For the purposes of this example, we are using an alanine pentapeptide. You can
find the following files that you will need for this demonstration in the
``examples/amber`` directory of the ParmEd distribution:

    * ``ala5_gas.parm7``
    * ``ala5_gas.rst7``

The following sample script (``simulate_amber_gb.py`` in the ParmEd
distribution) will set up and run the simulation using OpenMM::

    #!/usr/bin/env python
    from __future__ import division, print_function
    
    import sys
    
    # OpenMM Imports
    import simtk.openmm as mm
    import simtk.openmm.app as app
    
    # ParmEd Imports
    from parmed import load_file, unit as u
    from parmed.openmm import StateDataReporter, NetCDFReporter
    
    # Load the Amber files
    print('Loading AMBER files...')
    ala5_gas = load_file('ala5_gas.parm7', 'ala5_gas.rst7')
    
    # Create the OpenMM system
    print('Creating OpenMM System')
    system = ala5_gas.createSystem(nonbondedMethod=app.NoCutoff,
                                   constraints=app.HBonds, implicitSolvent=app.GBn2,
                                   implicitSolventSaltConc=0.1*u.moles/u.liter,
    )
    
    # Create the integrator to do Langevin dynamics
    integrator = mm.LangevinIntegrator(
                            300*u.kelvin,       # Temperature of heat bath
                            1.0/u.picoseconds,  # Friction coefficient
                            2.0*u.femtoseconds, # Time step
    )
    
    # Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
    # the platform to use the default (fastest) platform
    platform = mm.Platform.getPlatformByName('CUDA')
    prop = dict(CudaPrecision='mixed') # Use mixed single/double precision
    
    # Create the Simulation object
    sim = app.Simulation(ala5_gas.topology, system, integrator, platform, prop)
    
    # Set the particle positions
    sim.context.setPositions(ala5_gas.positions)
    
    # Minimize the energy
    print('Minimizing energy')
    sim.minimizeEnergy(maxIterations=500)
    
    # Set up the reporters to report energies and coordinates every 100 steps
    sim.reporters.append(
            StateDataReporter(sys.stdout, 100, step=True, potentialEnergy=True,
                                   kineticEnergy=True, temperature=True)
    )
    sim.reporters.append(
            NetCDFReporter('ala5_gb.nc', 100, crds=True)
    )
    
    # Run dynamics
    print('Running dynamics')
    sim.step(10000)

Now I'll dissect the script to help you understand what is happening at each
step. We will divide the script into the sections following the ``print``
statements that announce when each stage begins.

Loading Amber files
~~~~~~~~~~~~~~~~~~~

In this stage, we simply instantiate the :class:`AmberParm
<parmed.amber.AmberParm>` object from the input topology and coordinate
files. After this command, ``ala5_gas`` will contain a full description of every
particle, the parameters defining their interactions, and their positions.

Create the OpenMM System
~~~~~~~~~~~~~~~~~~~~~~~~

This command creates an OpenMM ``System`` object from the information stored in
``ala5_gas``. It contains multiple ``Force`` instances for the bonds, angles,
periodic torsions, and nonbonded (electrostatic and van der Waals) interactions.
It is in this function that we define the potential parameters we want to use.
In this example, we have chosen the default values for each parameter except the
ones specified. In particular:

    * ``nonbondedMethod=app.NoCutoff`` indicates we do not want to use a cutoff
      for nonbonded interactions. If you wanted to use a cutoff, you could use
      ``app.CutoffNonPeriodic`` instead (since this system does *not* use
      periodic boundary conditions)
    * ``constraints=app.HBonds`` indicates we want to constrain all bonds in
      which at least one atom is a Hydrogen (i.e., SHAKE or SETTLE for water).
      Other options are ``None`` (no constraints), ``app.AllBonds``, or
      ``app.HAngles``. For the most part, these are self-explanatory, but it is
      worth noting that ``app.HAngles`` constrains all bonds and the distance
      between the 1-3 pairs of angles in which one of those atoms is a hydrogen.
    * ``implicitSolvent=app.GBn2`` indicates we want to use the second GBneck
      model described in Nguyen et al., J. Chem. Theory Comput., 2014 9(4) p.
      2020-2034. Other options are ``app.HCT``, ``app.OBC1``, ``app.OBC2``,
      and ``app.GBn``. These correspond to values of 1, 2, 5, 7 (and 8 for
      ``app.GBn2``) to the ``igb`` variable in AMBER input files.
    * ``implicitSolventSaltConc=0.1*u.liters/u.mole`` indicates we want to model
      a ca. 0.1 molar solution of monovalent ions using a Debye screening model.

Create the integrator to do Langevin Dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this stage we specify an integrator. Common choices are
``LangevinIntegrator`` (as we've chosen here) to do simulations in the NVT
ensemble and ``VerletIntegrator`` that allows us to do simulations either at
constant energy or temperature if using the ``AndersenThermostat``.  In this
example, we've chosen the Langevin integrator with a target temperature of
300 K, a friction coefficient of 1/ps and a time step of 2 fs.

Define the platform
~~~~~~~~~~~~~~~~~~~

In this stage, we define the platform we want to use. In this example, we have
chosen the ``CUDA`` platform, but this may not be available on every machine
since it only runs on NVidia GPU hardware. Other choices are ``OpenCL`` (which
will run on a variety of GPUs including those made by AMD/ATI and CPUs), ``CPU``
(which is an optimized version that runs natively on CPUs), and ``Reference``
(often quite slow).

The properties can be set for each platform. In this case, we specified that we
wanted to use a ``mixed`` precision model (a good compromise between speed and
precision).

Create the ``Simulation`` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step creates a ``Simulation`` object that will be used to run the actual
simulations.  If we wanted OpenMM to simply pick the fastest platform for us
(rather than specify one directly), we could omit the ``platform`` and ``prop``
arguments.

Set the particle positions
~~~~~~~~~~~~~~~~~~~~~~~~~~

This stage is very important.  In this step, we set the particle positions
stored in the ``ala5_gas`` object to our object. If you omit this step, you can
get strange results or other errors like segmentation violations. These particle
positions have been parsed from the input coordinate file, although if you had a
PDB file you could use the OpenMM ``PDBFile`` object as a source of coordinates
instead.

Minimize the energy
~~~~~~~~~~~~~~~~~~~

This stage performs a basic energy minimization to relax particle positions.
This particular invocation will perform at most 500 iterations.

Set up the reporters
~~~~~~~~~~~~~~~~~~~~

This stage defines reporters that will "report" on the status of the simulation
periodically throughout the simulation. The first is an ``StateDataReporter``
which will print out a summary of energies and temperatures every 100 steps.
Unlike the ``StateDataReporter`` in OpenMM, this reporter prints values in the
AKMA unit system (Angstrom, Kilocalorie per mole, and atomic mass units). This
reporter directs the printout to standard output (the screen), ``sys.stdout``
can be replaced with a different file-like object or a file name.

The second reporter is a NetCDF trajectory reporter, which is written in the
Amber NetCDF format.  You can also use the native ``DCDReporter`` reporter in
OpenMM to print DCD-format trajectories.

Running dynamics
~~~~~~~~~~~~~~~~

This is the stage that actually runs the MD. In this case, we are running 10,000
steps of MD.  The wiki page with "Common recipes" has information regarding
running a long simulation in chunks.

Explicit Solvent
----------------

For the purposes of this example, we are using an alanine dipeptide solvated in
a box of water. You can find the following files that you will need for this
demonstration in the ``examples/amber`` directory of the ParmEd distribution:

    * ``ala2_solv.parm7``
    * ``ala2_solv.rst7``

The following sample script (``simulate_amber_pme.py`` in the ParmEd
distribution) will set up and run the simulation using OpenMM::

    #!/usr/bin/env python
    from __future__ import division, print_function
    
    import sys
    
    # OpenMM Imports
    import simtk.openmm as mm
    import simtk.openmm.app as app
    
    # ParmEd Imports
    from parmed import load_file, unit as u
    from parmed.openmm import StateDataReporter, NetCDFReporter
    
    # Load the Amber files
    print('Loading AMBER files...')
    ala2_solv = load_file('ala2_solv.parm7', 'ala2_solv.rst7')
    
    # Create the OpenMM system
    print('Creating OpenMM System')
    system = ala2_solv.createSystem(nonbondedMethod=app.PME,
                                    nonbondedCutoff=8.0*u.angstroms,
                                    constraints=app.HBonds,
    )
    
    # Create the integrator to do Langevin dynamics
    integrator = mm.LangevinIntegrator(
                            300*u.kelvin,       # Temperature of heat bath
                            1.0/u.picoseconds,  # Friction coefficient
                            2.0*u.femtoseconds, # Time step
    )
    
    # Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
    # the platform to use the default (fastest) platform
    platform = mm.Platform.getPlatformByName('CUDA')
    prop = dict(CudaPrecision='mixed') # Use mixed single/double precision
    
    # Create the Simulation object
    sim = app.Simulation(ala2_solv.topology, system, integrator, platform, prop)
    
    # Set the particle positions
    sim.context.setPositions(ala2_solv.positions)
    
    # Minimize the energy
    print('Minimizing energy')
    sim.minimizeEnergy(maxIterations=500)
    
    # Set up the reporters to report energies and coordinates every 100 steps
    sim.reporters.append(
            StateDataReporter(sys.stdout, 100, step=True, potentialEnergy=True,
                              kineticEnergy=True, temperature=True, volume=True,
                              density=True)
    )
    sim.reporters.append(NetCDFReporter('ala2_solv.nc', 100, crds=True))
    
    # Run dynamics
    print('Running dynamics')
    sim.step(10000)

Now we'll dissect the script to help you understand what is happening at each
step. We will divide the script into the sections following the ``print``
statements that announce when each stage begins.

Loading Amber files
~~~~~~~~~~~~~~~~~~~

In this stage, we simply instantiate the :class:`AmberParm
<parmed.amber.AmberParm>` object from the input topology and coordinate
files. After this command, ``ala2_solv`` will contain a full description of
every particle, the parameters defining their interactions, and their positions.

Create the OpenMM system
~~~~~~~~~~~~~~~~~~~~~~~~

This command creates an OpenMM ``System`` object from the information stored in
``ala5_gas``. It contains multiple ``Force`` instances for the bonds, angles,
periodic torsions, and nonbonded (electrostatic and van der Waals) interactions.
It is in this function that we define the potential parameters we want to use.
In this example, we have chosen the default values for each parameter except the
ones specified. In particular:

    * ``nonbondedMethod=app.PME`` indicates we want to use the Particle Mesh
      Ewald method to compute the full-range electrostatics.
    * ``nonbondedCutoff=8.0*u.angstroms`` indicates we want to use an 8 Angstrom
      cutoff for the Lennard-Jones interaction (as well as the direct-space part
      of the Ewald sum).
    * ``constraints=app.HBonds`` indicates that we want to constrain all bonds
      in which at least one atom is hydrogen

If there are any other force objects you want to define, they can be added to
the system after this step (like, for instance, positional restraints to a
reference structure).

Create the integrator to do Langevin Dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this stage we specify an integrator. Common choices are
``LangevinIntegrator`` (as we've chosen here) to do simulations in the NVT
ensemble and ``VerletIntegrator`` that allows us to do simulations either at
constant energy or temperature if using the ``AndersenThermostat``.  In this
example, we've chosen the Langevin integrator with a target temperature of
300 K, a friction coefficient of 1/ps and a time step of 2 fs.

Define the platform
~~~~~~~~~~~~~~~~~~~

In this stage, we define the platform we want to use.  In this example we have
chosen the CUDA platform, but this may not be available on every machine since
it only runs on NVidia GPU hardware.  Other choices are OpenCL (which will run
on a variety of GPUs including those made by AMD/ATI and CPUs), CPU (which is
an optimized version that runs natively on CPUs), and Reference (often quite
slow).

The properties can be set for each platform. In this case, we specified that we
wanted to use a mixed precision model (a good compromise between speed and
precision).

Create the Simulation object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step creates a ``Simulation`` object that will be used to run the actual
simulations.  If we wanted OpenMM to simply pick the fastest platform for us
(rather than specify one directly), we could omit the ``platform`` and ``prop``
arguments.

Set the particle positions
~~~~~~~~~~~~~~~~~~~~~~~~~~

This stage is very important.  In this step, we set the particle positions
stored in the ``ala5_gas`` object to our object.  If you omit this step, you can
get strange results or other errors like segmentation violations. These particle
positions have been parsed from the input coordinate file, although if you had a
PDB file you could use the OpenMM ``PDBFile`` object as a source of coordinates
instead.

Minimize the energy
~~~~~~~~~~~~~~~~~~~

This stage performs a basic energy minimization to relax particle positions.
This particular invocation will perform at most 500 iterations.

Set up the reporters
~~~~~~~~~~~~~~~~~~~~

This stage defines reporters that will "report" on the status of the simulation
periodically throughout the simulation. The first is an
:class:`StateDataReporter` which will print out a summary of energies and
temperatures every 100 steps.  Unlike the ``StateDataReporter`` in OpenMM, this
reporter prints values in the AKMA unit system (Angstrom, Kilocalorie per mole,
and atomic mass units).

The second reporter is a NetCDF trajectory reporter, which is written in the
Amber NetCDF format. You can also use the native ``DCDReporter`` reporter in
OpenMM to print DCD-format trajectories.

Running dynamics
~~~~~~~~~~~~~~~~

This is the stage that actually runs the MD.  In this case, we are running
10,000 steps of MD.  The wiki page with "Common recipes" has information
regarding running a long simulation in chunks.
