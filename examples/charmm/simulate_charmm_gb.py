#!/usr/bin/env python
from __future__ import division, print_function

import sys

# OpenMM Imports
import simtk.openmm as mm
import simtk.openmm.app as app

# ParmEd Imports
from parmed import load_file, unit as u
from parmed.charmm import CharmmParameterSet
from parmed.openmm import StateDataReporter

# Load the CHARMM files
print('Loading CHARMM files...')
params = CharmmParameterSet('toppar/par_all36_prot.prm')
ala5_gas = load_file('ala5_autopsf.psf')
ala5_crds = load_file('ala5_autopsf.pdb')

# NOTE NOTE
# The parameter set we used here is the CHARMM 36 force field, but this is
# strictly an example. It is important that you use the most accurate (typically
# most up-to-date) force fields for your own simulation. See the CHARMM
# parameter web page for updates:
# http://mackerell.umaryland.edu/CHARMM_ff_params.html
# END NOTE

# Create the OpenMM system
print('Creating OpenMM System')
system = ala5_gas.createSystem(params, nonbondedMethod=app.NoCutoff,
                               constraints=app.HBonds, implicitSolvent=app.HCT,
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
sim.context.setPositions(ala5_crds.positions)

# Minimize the energy
print('Minimizing energy')
sim.minimizeEnergy(maxIterations=500)

# Set up the reporters to report energies and coordinates every 100 steps
sim.reporters.append(
        StateDataReporter(sys.stdout, 100, step=True, potentialEnergy=True,
                          kineticEnergy=True, temperature=True,
                          volume=True, density=True)
)
sim.reporters.append(app.DCDReporter('ala5_gb.dcd', 100))

# Run dynamics
print('Running dynamics')
sim.step(10000)
