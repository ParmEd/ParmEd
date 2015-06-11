#!/usr/bin/env python

"""
Cross-compare energies for two sets of input files to validate a conversion tool. One set must be AMBER prmtop/prmcrd files and the other GROMACS topology/coordinate files. Both will be used to create ParmEd Amber/Gromacs topology files and then used to set up systems in OpenMM for energy evaluation.

Written (based on the below) May 2015 by David Mobley.

Originally written as test.py by Lee-Ping Wang and provided in GmxTests of ParmEd. Test.py did the following:

"Test ParmEd's ability to process a Gromacs position/topology file
by comparing Gromacs energy/force to OpenMM-via-ParmEd energy/force.

This script uses ForceBalance to obtain the Gromacs energy/force
and also reads parts of the Gromacs .mdp file to set up the system.

There are also some OpenMM imports for calculating the OpenMM energy/force..

To run this script, provide a gro, top and mdp file.

Author: Lee-Ping Wang
"

"""

# General import
from collections import OrderedDict
import numpy as np
import os, sys, re, subprocess

# ForceBalance convenience functions
from nifty import printcool, printcool_dictionary, _exec, which, wopen, isint, logger
# Only needed for writing constrained .gro files
# from molecule import Molecule

# ParmEd import
from chemistry import gromacs
from chemistry import amber

# OpenMM import
import simtk.unit as u
import simtk.openmm as mm
import simtk.openmm.app as app


def energy_components(Sim, verbose=False):
    # Before using EnergyComponents, make sure each Force is set to a different group.
    EnergyTerms = OrderedDict()
    if type(Sim.integrator) in [mm.LangevinIntegrator, mm.VerletIntegrator]:
        for i in range(Sim.system.getNumForces()):
            EnergyTerms[Sim.system.getForce(i).__class__.__name__] = Sim.context.getState(getEnergy=True,groups=2**i).getPotentialEnergy() / u.kilojoules_per_mole
    EnergyTerms['Potential'] = Sim.context.getState(getEnergy=True).getPotentialEnergy() / u.kilojoules_per_mole
    return EnergyTerms

def Calculate_ParmEd(gro_file, top_file, sysargs, defines):
    #===============================#
    #|   ParmEd object creation    |#
    #===============================#
    # Make sure the proper defines from the .mdp file are passed into the GromacsTopologyFile() :)
    ParmEd_GmxTop = gromacs.GromacsTopologyFile(top_file )
    ParmEd_GmxGro = gromacs.GromacsGroFile.parse(gro_file)
    ParmEd_GmxTop.box = ParmEd_GmxGro.box
    ParmEd_GmxTop.positions = ParmEd_GmxGro.positions
    #===============================#
    #|   OpenMM simulation setup   |#
    #===============================#
    # ParmEd creates System object
    system = ParmEd_GmxTop.createSystem(**sysargs)
    # Keep a record of which atoms are real (not virtual sites)
    isAtom = []
    for i in range(system.getNumParticles()):
        isAtom.append(system.getParticleMass(i).value_in_unit(u.dalton) > 0.0)
    # Setting force groups enables energy components analysis
    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)
        if isinstance(f, mm.NonbondedForce):
            f.setUseDispersionCorrection(True)
    integ = mm.VerletIntegrator(1.0*u.femtosecond)
    plat = mm.Platform.getPlatformByName('Reference')
    # Create Simulation object
    simul = app.Simulation(ParmEd_GmxTop.topology, system, integ, plat)
    simul.context.setPositions(ParmEd_GmxGro.positions)
    simul.context.applyConstraints(1e-12)
    # Obtain OpenMM potential energy
    state = simul.context.getState(getPositions=True,getEnergy=True,getForces=True)
    parmed_energy = state.getPotentialEnergy()
    parmed_forces = state.getForces()
    pos = np.array(state.getPositions().value_in_unit(u.angstrom)).reshape(-1,3)
    # Obtain and save constrained positions
    # M = Molecule(gro_file)
    # M.xyzs[0] = pos
    # M.write('constrained.gro')
    # Print OpenMM-via-ParmEd energy components
    Ecomps_OMM = energy_components(simul)
    printcool_dictionary(Ecomps_OMM, title="OpenMM energy components via ParmEd")
    parmed_forces = np.array([f for i, f in enumerate(parmed_forces.value_in_unit(u.kilojoule_per_mole/u.nanometer)) if isAtom[i]])
    return parmed_energy, parmed_forces, Ecomps_OMM

def Calculate_ParmEd_Amber(crd_file, prmtop_file, sysargs):
    #===============================#
    #|   ParmEd object creation    |#
    #===============================#
    ParmEd_ambertop = amber.AmberParm( prmtop_file, crd_file)    

    system = ParmEd_ambertop.createSystem(**sysargs)

    #===============================#
    #|   OpenMM simulation setup   |#
    #===============================#
    # Keep a record of which atoms are real (not virtual sites)
    isAtom = []
    for i in range(system.getNumParticles()):
        isAtom.append(system.getParticleMass(i).value_in_unit(u.dalton) > 0.0)
    # Setting force groups enables energy components analysis
    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)
        if isinstance(f, mm.NonbondedForce):
            f.setUseDispersionCorrection(True)
    integ = mm.VerletIntegrator(1.0*u.femtosecond)
    plat = mm.Platform.getPlatformByName('Reference')
    # Create Simulation object
    simul = app.Simulation(ParmEd_ambertop.topology, system, integ, plat)
    simul.context.setPositions(ParmEd_ambertop.positions)
    simul.context.applyConstraints(1e-12)
    # Obtain OpenMM potential energy
    state = simul.context.getState(getPositions=True,getEnergy=True,getForces=True)
    parmed_energy = state.getPotentialEnergy()
    parmed_forces = state.getForces()
    pos = np.array(state.getPositions().value_in_unit(u.angstrom)).reshape(-1,3)
    # Obtain and save constrained positions
    # M = Molecule(gro_file)
    # M.xyzs[0] = pos
    # M.write('constrained.gro')
    # Print OpenMM-via-ParmEd energy components
    Ecomps_OMM = energy_components(simul)
    printcool_dictionary(Ecomps_OMM, title="OpenMM energy components via ParmEd")
    parmed_forces = np.array([f for i, f in enumerate(parmed_forces.value_in_unit(u.kilojoule_per_mole/u.nanometer)) if isAtom[i]])
    return parmed_energy, parmed_forces, Ecomps_OMM

def main():
    # Command line arguments
    crd_file = sys.argv[1]
    prmtop_file = sys.argv[2]
    gro_file = sys.argv[3]
    top_file = sys.argv[4]

    #Define sysargs for testing
    sysargs = {}
    sysargs['rigidWater'] = False
    sysargs['constraints'] = app.HBonds
    #sysargs['switchDistance'] = 0.9 * u.nanometer
    #sysargs['nonbondedMethod'] = app.PME
    sysargs['nonbondedMethod'] = app.NoCutoff
    #sysargs['ewaldErrorTolerance'] = 1e-5
    #sysargs['nonbondedCutoff'] = 0.9 * u.nanometer

    # ParmEd-OpenMM calculation with GROMACS inputs
    print "Evaluating energy from GROMACS inputs"
    PED_Energy, PED_Force, Ecomps_PED = Calculate_ParmEd(gro_file, top_file, sysargs, defines = [])
    # ParmEd-OpenMM calculation with AMBER inputs
    print "Evaluating energy from AMBER inputs"
    PED_Energy_Amber, PED_Force_Amber, Ecomps_PED_Amber = Calculate_ParmEd_Amber(crd_file, prmtop_file, sysargs)
    print "Finished evaluating energies" 

    # Analyze force differences
    D_Force = PED_Force - PED_Force_Amber

    # Final printout
    print "Energy Difference (kJ/mol):"
    #print Ecomps_PED, Ecomps_PED_Amber
    #print Ecomps_PED['Potential'], Ecomps_PED_Amber['Potential']
    #print Ecomps_PED['Potential']-Ecomps_PED_Amber['Potential']
    print (Ecomps_PED['Potential']-Ecomps_PED_Amber['Potential'])
    print "RMS / Max Force Difference (kJ/mol/nm):"
    print np.sqrt(np.mean([sum(i**2) for i in D_Force])), np.sqrt(np.max(np.array([sum(i**2) for i in D_Force])))

if __name__ == "__main__":
    main()
