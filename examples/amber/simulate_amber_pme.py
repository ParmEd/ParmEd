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
