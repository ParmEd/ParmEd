#! /usr/bin/env python

from __future__ import print_function

# OpenMM imports
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

# ParmEd imports
from parmed import load_rosetta

# PyRosetta imports
from toolbox import mutate_residue
from rosetta import init, pose_from_sequence

# Initialize PyRosetta
init()

# Create Ala12
print('Creating Ala12...')
seq = 12*'A'
pose = pose_from_sequence(seq)

# Mutate residue 5 to K
print('Mutating Fifth Ala to Lys...')
mutant = mutate_residue(pose, 5, 'K')

# Load mutant into ParmEd
print('Loading into ParmEd...')
struct = load_rosetta(mutant)

# Load AMBER-99SBildn and TIP3P parameters
print('Loading AMBER parameters...')
ff = ForceField('amber99sbildn.xml', 'tip3p.xml')

# Solvate Structure
print('Solvating...')
mod = Modeller(struct.topology, struct.positions)
mod.addSolvent(ff, model='tip3p', boxSize=Vec3(4, 4, 4)*nanometers,
               positiveIon='Na+', negativeIon='Cl-',
               ionicStrength=0.1*molar)

# Create OpenMM System
print('Creating OpenMM System...')
system = ff.createSystem(mod.topology, nonbondedMethod=PME,
                         nonbondedCutoff=1*nanometers,
                         constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picoseconds,
                                2*femtoseconds)

system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

simulation = Simulation(mod.topology, system, integrator)
simulation.context.setPositions(mod.positions)

# Minimize System
print('Minimizing...')
simulation.minimizeEnergy(maxIterations=1000)

# Write to PDB
simulation.reporters.append(
    PDBReporter('dodecaalanine.solv.pdb', 50000))

# Equilibrate System
print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(300)

simulation.step(50000)

print('Done!')
print('Results saved to dodecaalanine.solv.pdb')
