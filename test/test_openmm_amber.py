"""
Contains unittests for running OpenMM calculations using the Amber file parsers
"""
from __future__ import division, print_function, absolute_import
import utils

try:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    has_openmm = True
except ImportError:
    from parmed.amber.readparm import AmberParm, ChamberParm, Rst7
    has_openmm = False

from parmed.amber import AmberParm, ChamberParm, Rst7
from parmed.openmm import load_topology
import parmed.unit as u
from parmed.utils.six.moves import range, zip
from copy import copy
from math import sqrt
import os
import parmed.tools as PT
import sys
import unittest

get_fn = utils.get_fn

if has_openmm:
    CPU = mm.Platform.getPlatformByName('CPU')


# OpenMM NonbondedForce methods are enumerated values. From NonbondedForce.h,
# they are:
#   0 - NoCutoff
#   1 - CutoffNonPeriodic
#   2 - CutoffPeriodic
#   3 - Ewald
#   4 - PME

def decomposed_energy(context, parm, NRG_UNIT=u.kilocalories_per_mole):
    """ Gets a decomposed energy for a given system """
    energies = {}
    # Get energy components
    s = context.getState(getEnergy=True,
                         enforcePeriodicBox=parm.ptr('ifbox')>0,
                         groups=2**parm.BOND_FORCE_GROUP)
    energies['bond'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
    s = context.getState(getEnergy=True,
                         enforcePeriodicBox=parm.ptr('ifbox')>0,
                         groups=2**parm.ANGLE_FORCE_GROUP)
    energies['angle'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
    s = context.getState(getEnergy=True,
                         enforcePeriodicBox=parm.ptr('ifbox')>0,
                         groups=2**parm.DIHEDRAL_FORCE_GROUP)
    energies['dihedral'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
    s = context.getState(getEnergy=True,
                         enforcePeriodicBox=parm.ptr('ifbox')>0,
                         groups=2**parm.NONBONDED_FORCE_GROUP)
    energies['nonbond'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
    # Extra energy terms for chamber systems
    if isinstance(parm, ChamberParm):
        s = context.getState(getEnergy=True,
                             enforcePeriodicBox=parm.ptr('ifbox')>0,
                             groups=2**parm.UREY_BRADLEY_FORCE_GROUP)
        energies['urey'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
        s = context.getState(getEnergy=True,
                             enforcePeriodicBox=parm.ptr('ifbox')>0,
                             groups=2**parm.IMPROPER_FORCE_GROUP)
        energies['improper']=s.getPotentialEnergy().value_in_unit(NRG_UNIT)
        s = context.getState(getEnergy=True,
                             enforcePeriodicBox=parm.ptr('ifbox')>0,
                             groups=2**parm.CMAP_FORCE_GROUP)
        energies['cmap'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
    return energies

class TestAmberParm(utils.FileIOTestCase, utils.TestCaseRelative):

    def testEPEnergy(self):
        """ Tests AmberParm handling of extra points in TIP4P water """
        parm = AmberParm(get_fn('tip4p.parm7'), get_fn('tip4p.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   constraints=app.HBonds,
                                   rigidWater=True,
                                   flexibleConstraints=False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
# Etot   =     -1756.2018  EKtot   =       376.7454  EPtot      =     -2132.9472
# BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
# 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =       378.4039
# EELEC  =     -2511.3511  EHBOND  =         0.0000  RESTRAINT  =         0.0000
# EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =      6514.3661
        self.assertAlmostEqual(energies['bond'], 0)
        self.assertAlmostEqual(energies['angle'], 0)
        self.assertAlmostEqual(energies['dihedral'], 0)
        self.assertRelativeEqual(energies['nonbond'], -2133.2963, places=4)
        # Check that we have the correct number of virtual sites
        nvirt = 0
        for i in range(system.getNumParticles()):
            nvirt += system.isVirtualSite(i)
        self.assertEqual(parm.ptr('NUMEXTRA'), nvirt)
        # Now test the forces to make sure that they are computed correctly in
        # the presence of extra points
        pstate = sim.context.getState(getForces=True)
        sstate = mm.XmlSerializer.deserialize(
                open(utils.get_saved_fn('tip4pforces.xml'), 'r').read()
        )
        pf = pstate.getForces().value_in_unit(u.kilocalorie_per_mole/u.angstrom)
        sf = sstate.getForces().value_in_unit(u.kilocalorie_per_mole/u.angstrom)

        for p, s in zip(pf, sf):
            for x1, x2 in zip(p, s):
                # Compare large forces relatively and small ones absolutely
                if abs(x1) > 1 or abs(x2) > 1:
                    self.assertRelativeEqual(x1, x2, places=2)
                else:
                    self.assertAlmostEqual(x1, x2, delta=2e-2)

    @unittest.skipIf(utils.skip_big_tests(), "Skipping long tests")
    def testRoundTripEP(self):
        """ Test ParmEd -> OpenMM round trip with Amber EPs and PME """
        parm = AmberParm(get_fn('tip4p.parm7'), get_fn('tip4p.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   constraints=app.HBonds,
                                   rigidWater=True,
                                   flexibleConstraints=True)
        system2 = load_topology(parm.topology, system).createSystem(
                                    nonbondedMethod=app.PME,
                                    nonbondedCutoff=8*u.angstroms,
                                    constraints=app.HBonds,
                                    rigidWater=True,
                                    flexibleConstraints=True)
        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system2, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(parm.positions)
        con2.setPositions(parm.positions)
        e1 = decomposed_energy(con1, parm)
        e2 = decomposed_energy(con2, parm)
        self.assertAlmostEqual(e1['bond'], e2['bond'])
        self.assertAlmostEqual(e1['angle'], e2['angle'])
        self.assertAlmostEqual(e1['dihedral'], e2['dihedral'])
        self.assertAlmostEqual(e1['nonbond'], e2['nonbond'], places=5)
        # Check that we have the correct number of virtual sites
        nvirt1 = nvirt2 = 0
        for i in range(system.getNumParticles()):
            nvirt1 += system.isVirtualSite(i)
            nvirt2 += system2.isVirtualSite(i)
        self.assertEqual(nvirt1, nvirt2)
        # Now test the forces to make sure that they are computed correctly in
        # the presence of extra points
        state1 = con1.getState(getForces=True)
        state2 = con2.getState(getForces=True)
        f1 = state1.getForces().value_in_unit(u.kilocalorie_per_mole/u.angstrom)
        f2 = state2.getForces().value_in_unit(u.kilocalorie_per_mole/u.angstrom)

        for p, s in zip(f1, f2):
            for x1, x2 in zip(p, s):
                self.assertAlmostEqual(x1, x2, places=3)

    def testEPEnergy2(self):
        """ Tests AmberParm handling of extra points in TIP5P water """
        parm = AmberParm(get_fn('tip5p.parm7'), get_fn('tip5p.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   constraints=app.HBonds,
                                   rigidWater=True,
                                   flexibleConstraints=False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 0)
        self.assertAlmostEqual(energies['angle'], 0)
        self.assertAlmostEqual(energies['dihedral'], 0)
        self.assertRelativeEqual(energies['nonbond'], -2142.418956, places=4)
        # Check that we have the correct number of virtual sites
        nvirt = 0
        for i in range(system.getNumParticles()):
            nvirt += system.isVirtualSite(i)
        self.assertEqual(parm.ptr('NUMEXTRA'), nvirt)
        # Now test the forces to make sure that they are computed correctly in
        # the presence of extra points
        pstate = sim.context.getState(getForces=True)
        sstate = mm.XmlSerializer.deserialize(
                open(utils.get_saved_fn('tip5pforces.xml'), 'r').read()
        )
        pf = pstate.getForces().value_in_unit(u.kilocalorie_per_mole/u.angstrom)
        sf = sstate.getForces().value_in_unit(u.kilocalorie_per_mole/u.angstrom)

        i = 0
        for p, s in zip(pf, sf):
            if system.isVirtualSite(i):
                i += 1
                continue # Skip forces on virtual sites
            for x1, x2 in zip(p, s):
                # Compare large forces relatively and small ones absolutely
                if abs(x1) > 1 or abs(x2) > 1:
                    self.assertRelativeEqual(x1, x2, delta=5e-3)
                else:
                    self.assertAlmostEqual(x1, x2, delta=5e-3)
            i += 1

    def testGasEnergy(self):
        """ Compare Amber and OpenMM gas phase energies """
        parm = AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem() # Default, no cutoff
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =         2.4544  EKtot   =         0.0000  EPtot      =         2.4544
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =         0.0000  RESTRAINT  =         0.0000
        # Compare OpenMM energies with the Amber energies (above)
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbond'], -30.2355, places=3)

    def testRoundTrip(self):
        """ Test ParmEd -> OpenMM round trip with Amber gas phase """
        parm = AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem()
        system2 = load_topology(parm.topology, system).createSystem()
        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(parm.positions)
        con2.setPositions(parm.positions)
        e1 = decomposed_energy(con1, parm)
        e2 = decomposed_energy(con2, parm)
        self.assertAlmostEqual(e1['bond'], e2['bond'])
        self.assertAlmostEqual(e1['angle'], e2['angle'])
        self.assertAlmostEqual(e1['dihedral'], e2['dihedral'])
        self.assertAlmostEqual(e1['nonbond'], e2['nonbond'])

    def testRoundTripXML(self):
        """ Test ParmEd -> OpenMM round trip with Amber gas phase via XML """
        parm = AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem()
        fname = get_fn('ash.xml', written=True)
        with open(fname, 'w') as f:
            f.write(mm.XmlSerializer.serialize(system))
        system2 = load_topology(parm.topology, fname).createSystem()
        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(parm.positions)
        con2.setPositions(parm.positions)
        e1 = decomposed_energy(con1, parm)
        e2 = decomposed_energy(con2, parm)
        self.assertAlmostEqual(e1['bond'], e2['bond'])
        self.assertAlmostEqual(e1['angle'], e2['angle'])
        self.assertAlmostEqual(e1['dihedral'], e2['dihedral'])
        self.assertAlmostEqual(e1['nonbond'], e2['nonbond'])

    def testGB1Energy(self): # HCT (igb=1)
        """ Compare Amber and OpenMM GB (igb=1) energies (w/ and w/out salt) """
        parm = AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(implicitSolvent=app.HCT)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0306  EKtot   =         0.0000  EPtot      =       -24.0306
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -26.4850  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbond'], -56.7205, places=3)
        system = parm.createSystem(implicitSolvent=app.HCT,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0508  EKtot   =         0.0000  EPtot      =       -24.0508
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -26.5051  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbond'], -56.7406, places=3)

    def testGB2Energy(self): # OBC1 (igb=2)
        """ Compare Amber and OpenMM GB (igb=2) energies (w/ and w/out salt) """
        parm = AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(implicitSolvent=app.OBC1)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0306  EKtot   =         0.0000  EPtot      =       -24.0306
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -28.1689  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbond'], -58.4044, places=3)
        system = parm.createSystem(implicitSolvent=app.OBC1,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0508  EKtot   =         0.0000  EPtot      =       -24.0508
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -28.1895  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbond'], -58.4250, places=3)

    def testGB5Energy(self): # OBC2 (igb=5)
        """ Compare Amber and OpenMM GB (igb=5) energies (w/ and w/out salt) """
        parm = AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(implicitSolvent=app.OBC2)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0306  EKtot   =         0.0000  EPtot      =       -24.0306
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -28.4639  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbond'], -55.6994, places=3)
        system = parm.createSystem(implicitSolvent=app.OBC2,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0508  EKtot   =         0.0000  EPtot      =       -24.0508
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -25.4835  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbond'], -55.7190, places=3)

    def testGB7Energy(self): # GBn (igb=7)
        """ Compare Amber and OpenMM GB (igb=7) energies (w/ and w/out salt) """
        parm = AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        PT.changeRadii(parm, 'mbondi3').execute() # Need new radius set
        PT.loadRestrt(parm, get_fn('ash.rst7')).execute() # Load crds into copy
        system = parm.createSystem(implicitSolvent=app.GBn)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0306  EKtot   =         0.0000  EPtot      =       -24.0306
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -21.1012  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbond'], -51.3367, places=3)
        system = parm.createSystem(implicitSolvent=app.GBn,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0508  EKtot   =         0.0000  EPtot      =       -24.0508
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -21.1194  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbond'], -51.3549, places=3)

    def testGB8Energy(self): # GBn2 (igb=8)
        """ Compare Amber and OpenMM GB (igb=8) energies (w/ and w/out salt) """
        parm = AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        PT.changeRadii(parm, 'mbondi3').execute() # Need new radius set
        PT.loadRestrt(parm, get_fn('ash.rst7')).execute() # Load crds into copy
        system = parm.createSystem(implicitSolvent=app.GBn2)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0306  EKtot   =         0.0000  EPtot      =       -24.0306
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -23.4639  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbond'], -53.6994, places=3)
        system = parm.createSystem(implicitSolvent=app.GBn2,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0508  EKtot   =         0.0000  EPtot      =       -24.0508
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -23.4832  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbond'], -53.7187, places=3)

    def testRst7(self):
        """ Test loading coordinates via the OpenMMRst7 class """
        parm = AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem() # Default, no cutoff
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(Rst7.open(get_fn('ash.rst7')).positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =         2.4544  EKtot   =         0.0000  EPtot      =         2.4544
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =         0.0000  RESTRAINT  =         0.0000
        # Compare OpenMM energies with the Amber energies (above)
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbond'], -30.2355, places=3)

    @unittest.skipIf(utils.skip_big_tests(), "Skipping long tests")
    def testEwald(self):
        """ Compare Amber and OpenMM Ewald energies """
        parm = AmberParm(get_fn('solv.prmtop'), get_fn('solv.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.Ewald,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =     250.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =    -91560.0534  EKtot   =         0.0000  EPtot      =    -91560.0534
#BOND   =       495.0414  ANGLE   =      1268.9447  DIHED      =      1764.7201
#1-4 NB =       610.9632  1-4 EEL =      6264.1100  VDWAALS    =     11213.7649
#EELEC  =   -113177.5977  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#Ewald error estimate:   0.2347E-02
        self.assertRelativeEqual(energies['bond'], 495.0414, 4)
        self.assertRelativeEqual(energies['angle'], 1268.9447, 5)
        self.assertRelativeEqual(energies['dihedral'], 1764.7201, 5)
        self.assertRelativeEqual(energies['nonbond'], -95078.7346, 4)

    def testPME(self):
        """ Compare Amber and OpenMM PME energies """
        parm = AmberParm(get_fn('solv.prmtop'), get_fn('solv.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =     250.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =    -91544.8269  EKtot   =         0.0000  EPtot      =    -91544.8269
#BOND   =       495.0414  ANGLE   =      1268.9447  DIHED      =      1764.7201
#1-4 NB =       610.9632  1-4 EEL =      6264.1100  VDWAALS    =     11213.7649
#EELEC  =   -113162.3712  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#Ewald error estimate:   0.8352E-05
        self.assertRelativeEqual(energies['bond'], 495.0414, places=4)
        self.assertRelativeEqual(energies['angle'], 1268.9447, places=4)
        self.assertRelativeEqual(energies['dihedral'], 1764.7201, places=4)
        self.assertRelativeEqual(energies['nonbond'], -95073.5331, places=4)

    def testDispersionCorrection(self):
        """ Compare Amber and OpenMM PME energies w/out vdW correction """
        parm = AmberParm(get_fn('solv.prmtop'), get_fn('solv.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        sys.stdout.flush()
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5)
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setUseDispersionCorrection(False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =     250.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =    -90623.3323  EKtot   =         0.0000  EPtot      =    -90623.3323
#BOND   =       495.0414  ANGLE   =      1268.9447  DIHED      =      1764.7201
#1-4 NB =       610.9632  1-4 EEL =      6264.1100  VDWAALS    =     12135.2596
#EELEC  =   -113162.3712  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#Ewald error estimate:   0.8352E-05
        self.assertRelativeEqual(energies['bond'], 495.0414, places=4)
        self.assertRelativeEqual(energies['angle'], 1268.9447, places=4)
        self.assertRelativeEqual(energies['dihedral'], 1764.7201, places=4)
        self.assertRelativeEqual(energies['nonbond'], -94152.0384, places=4)

    def testSHAKE(self):
        """ Compare Amber and OpenMM PME energies excluding SHAKEn bonds """
        parm = AmberParm(get_fn('solv.prmtop'), get_fn('solv.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5,
                                   flexibleConstraints=False,
                                   constraints=app.HBonds)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        # The only thing that changes here compared to the other periodic tests
        # is the bond energy, which should be slightly smaller than before
        state = sim.context.getState(getEnergy=True, enforcePeriodicBox=True,
                                     groups=2**parm.BOND_FORCE_GROUP)
        bond = state.getPotentialEnergy().value_in_unit(u.kilocalories_per_mole)
        self.assertRelativeEqual(bond, 494.5578, places=4)

    def testNBFIX(self):
        """ Compare Amber and OpenMM PME energies with NBFIX modifications """
        # For now, long-range correction is not available
        parm = AmberParm(get_fn('ff14ipq.parm7'), get_fn('ff14ipq.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        PT.change(parm, 'CHARGE', ':*', 0).execute() # only check LJ energies
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms)
        # Test without long-range correction, but check that the default value
        # is to have the long-range correction turned on (at time of writing
        # this test, the long-range correction segfaulted, but its cause is
        # known and should be fixed soon -- but I don't want to update the test
        # :-P).
        for force in system.getForces():
            if isinstance(force, mm.CustomNonbondedForce):
                self.assertTrue(force.getUseLongRangeCorrection())
                force.setUseLongRangeCorrection(False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =   193.6
#Etot   =     -7042.2475  EKtot   =         0.0000  EPtot      =     -7042.2475
#BOND   =         0.0654  ANGLE   =         0.9616  DIHED      =        -5.4917
#1-4 NB =        12.4186  1-4 EEL =       258.8388  VDWAALS    =      1243.9393
#EELEC  =     -8552.9795  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#EKCMT  =         0.0000  VIRIAL  =      -178.4985  VOLUME     =     42712.9055
#                                                   Density    =         0.6634
        self.assertAlmostEqual(energies['bond'], 0.0654, delta=0.0002)
        self.assertAlmostEqual(energies['angle'], 0.9616, 4)
        self.assertAlmostEqual(energies['dihedral'], -5.4917, 4)
        self.assertAlmostEqual(energies['nonbond'], 1256.3579, 3)

    def test1264(self):
        """ Testing the 12-6-4 LJ potential in OpenMM """
        parm = AmberParm(get_fn('znf_1264.prmtop'), get_fn('znf.rst'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.NoCutoff)
        for force in system.getForces():
            if isinstance(force, mm.CustomNonbondedForce):
                self.assertTrue(force.getUseLongRangeCorrection())
                force.setUseLongRangeCorrection(False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 26.3947079, places=3)
        self.assertAlmostEqual(energies['angle'], 122.8243431, places=3)
        self.assertAlmostEqual(energies['dihedral'], 319.0419347, places=3)
        self.assertAlmostEqual(energies['nonbond'], -2133.6170786, delta=2e-3)

    def test1012(self):
        """ Testing the 10-12 LJ H-bond potential in OpenMM """
        parm = AmberParm(get_fn('ff91.parm7'), get_fn('ff91.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms)
        for force in system.getForces():
            if isinstance(force, mm.CustomNonbondedForce):
                self.assertTrue(force.getUseLongRangeCorrection())
        integrator = mm.VerletIntegrator(1*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 0.9675961, places=3)
        self.assertAlmostEqual(energies['angle'], 82.5853211, places=3)
        self.assertAlmostEqual(energies['dihedral'], 1.2476012, places=3)
        self.assertRelativeEqual(energies['nonbond'], -11191.2563220, delta=3e-5)

    def testHangleConstraints(self):
        """ Tests that HAngle constraints get applied correctly """
        # This used to be a bug
        # Constrain just bonds. Make sure 1000 angles remain, and we have 2000
        # constraints
        parm = AmberParm(get_fn('gaffwat.parm7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=10.0*u.angstroms,
                                   constraints=app.HBonds, rigidWater=False,
                                   flexibleConstraints=False)
        nbonds = 0
        for force in system.getForces():
            if isinstance(force, mm.HarmonicBondForce):
                nbonds += force.getNumBonds()
        self.assertEqual(nbonds, 0)
        nangles = 0
        for force in system.getForces():
            if isinstance(force, mm.HarmonicAngleForce):
                nangles += force.getNumAngles()
        self.assertEqual(nangles, 1000)
        self.assertEqual(system.getNumConstraints(), 2000)
        # Constrain H-angles now, too
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=10.0*u.angstroms,
                                   constraints=app.HAngles, rigidWater=False,
                                   flexibleConstraints=False)
        nbonds = 0
        for force in system.getForces():
            if isinstance(force, mm.HarmonicBondForce):
                nbonds += force.getNumBonds()
        self.assertEqual(nbonds, 0)
        nangles = 0
        for force in system.getForces():
            if isinstance(force, mm.HarmonicAngleForce):
                nangles += force.getNumAngles()
        self.assertEqual(nangles, 0)
        self.assertEqual(system.getNumConstraints(), 3000)

    def testInterfacePBC(self):
        """ Testing all AmberParm.createSystem options (periodic) """
        parm = AmberParm(get_fn('solv.prmtop'), get_fn('solv.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=10.0*u.angstroms,
                                   constraints=None, rigidWater=False,
                                   removeCMMotion=True,
                                   ewaldErrorTolerance=1e-5)
        has_cmmotion = False
        for f in system.getForces():
            has_cmmotion = has_cmmotion or isinstance(f, mm.CMMotionRemover)
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getCutoffDistance(), 10.0*u.angstroms)
                self.assertEqual(f.getNonbondedMethod(), 4)
                self.assertEqual(f.getEwaldErrorTolerance(), 1e-5)
            elif isinstance(f, mm.HarmonicBondForce):
                self.assertEqual(f.getNumBonds(),
                                 parm.ptr('nbona')+parm.ptr('nbonh'))
            elif isinstance(f, mm.HarmonicAngleForce):
                self.assertEqual(f.getNumAngles(),
                                 parm.ptr('ntheth')+parm.ptr('ntheta'))
            elif isinstance(f, mm.PeriodicTorsionForce):
                self.assertEqual(f.getNumTorsions(),
                                 parm.ptr('nphih')+parm.ptr('nphia'))
        self.assertTrue(has_cmmotion)
        self.assertEqual(system.getNumConstraints(), 0)
        system = parm.createSystem(nonbondedMethod=app.Ewald,
                                   nonbondedCutoff=4.0*u.angstroms,
                                   constraints=app.HBonds, rigidWater=True,
                                   removeCMMotion=False)
        has_cmmotion = False
        for f in system.getForces():
            has_cmmotion = has_cmmotion or isinstance(f, mm.CMMotionRemover)
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getCutoffDistance(), 4.0*u.angstroms)
                self.assertEqual(f.getNonbondedMethod(), 3)
                self.assertEqual(f.getEwaldErrorTolerance(), 5e-4)
            elif isinstance(f, mm.HarmonicBondForce):
                self.assertEqual(f.getNumBonds(),
                                 parm.ptr('nbona')+parm.ptr('nbonh'))
            elif isinstance(f, mm.HarmonicAngleForce):
                self.assertEqual(f.getNumAngles(),
                                 parm.ptr('ntheta')+parm.ptr('ntheth'))
            elif isinstance(f, mm.PeriodicTorsionForce):
                self.assertEqual(f.getNumTorsions(),
                                 parm.ptr('nphih')+parm.ptr('nphia'))
        self.assertFalse(has_cmmotion)
        self.assertEqual(system.getNumConstraints(), parm.ptr('nbonh'))
        system = parm.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                   nonbondedCutoff=20.0*u.angstroms,
                                   constraints=app.AllBonds, rigidWater=True,
                                   flexibleConstraints=False)
        has_cmmotion = False
        for f in system.getForces():
            has_cmmotion = has_cmmotion or isinstance(f, mm.CMMotionRemover)
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getCutoffDistance(), 20.0*u.angstroms)
                self.assertEqual(f.getNonbondedMethod(), 2)
            elif isinstance(f, mm.HarmonicBondForce):
                self.assertEqual(f.getNumBonds(), 0)
            elif isinstance(f, mm.HarmonicAngleForce):
                self.assertEqual(f.getNumAngles(),
                                 parm.ptr('ntheta')+parm.ptr('ntheth'))
            elif isinstance(f, mm.PeriodicTorsionForce):
                self.assertEqual(f.getNumTorsions(),
                                 parm.ptr('nphih')+parm.ptr('nphia'))
        self.assertTrue(has_cmmotion)
        self.assertEqual(system.getNumConstraints(),
                         parm.ptr('nbonh')+parm.ptr('nbona'))
        system = parm.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                   nonbondedCutoff=20.0*u.angstroms,
                                   constraints=app.HBonds, rigidWater=True,
                                   flexibleConstraints=False,
                                   hydrogenMass=3.00*u.daltons)
        has_cmmotion = False
        for f in system.getForces():
            has_cmmotion = has_cmmotion or isinstance(f, mm.CMMotionRemover)
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getCutoffDistance(), 20.0*u.angstroms)
                self.assertEqual(f.getNonbondedMethod(), 2)
            elif isinstance(f, mm.HarmonicBondForce):
                self.assertEqual(f.getNumBonds(), parm.ptr('nbona'))
            elif isinstance(f, mm.HarmonicAngleForce):
                self.assertEqual(f.getNumAngles(),
                                 parm.ptr('ntheta')+parm.ptr('ntheth'))
            elif isinstance(f, mm.PeriodicTorsionForce):
                self.assertEqual(f.getNumTorsions(),
                                 parm.ptr('nphih')+parm.ptr('nphia'))
        self.assertTrue(has_cmmotion)
        self.assertEqual(system.getNumConstraints(), parm.ptr('nbonh'))
        totmass = 0
        for i, atom in enumerate(parm.atoms):
            mass = system.getParticleMass(i).value_in_unit(u.daltons)
            totmass += mass
            if atom.element == 1:
                self.assertAlmostEqual(mass, 3)
        self.assertAlmostEqual(totmass, sum(parm.parm_data['MASS']))
        # Trap some illegal options
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod=0))
        self.assertRaises(ValueError, lambda: parm.createSystem(constraints=0))

    def testInterfaceNoPBC(self):
        """ Testing all AmberParm.createSystem options (non-periodic) """
        parm = AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   constraints=app.HBonds,
                                   implicitSolvent=app.HCT,
                                   soluteDielectric=2.0,
                                   solventDielectric=80.0,
                                   flexibleConstraints=False)
        has_cmmotion = False
        for f in system.getForces():
            has_cmmotion = has_cmmotion or isinstance(f, mm.CMMotionRemover)
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)
        self.assertTrue(has_cmmotion)
        system = parm.createSystem(nonbondedMethod=app.CutoffNonPeriodic,
                                   nonbondedCutoff=100*u.angstroms,
                                   implicitSolvent=app.HCT,
                                   implicitSolventSaltConc=0.1*u.molar)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 1)
                self.assertEqual(f.getCutoffDistance(), 100*u.angstroms)
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.OBC1,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.OBC1,
                                   implicitSolventKappa=0.8)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.OBC2,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.OBC2,
                                   implicitSolventKappa=0.8)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn,
                                   implicitSolventSaltConc=10.0,
                                   implicitSolventKappa=0.8) # kappa is priority
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn2,
                                   implicitSolventSaltConc=10.0,
                                   implicitSolventKappa=0.8) # kappa is priority
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn2,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn2,
                                   solventDielectric=80, soluteDielectric=4,
                                   temperature=400.0*u.kelvin,
                                   implicitSolventSaltConc=0.1*u.molar)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)
        # Test some illegal options
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod=app.PME))
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod=app.CutoffPeriodic))
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod=app.Ewald))

class TestChamberParm(utils.TestCaseRelative):
    
    def testGasEnergy(self):
        """ Compare OpenMM and CHAMBER gas phase energies """
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'),
                           get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem() # Default, no cutoff
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =        39.1266  EKtot   =         0.0000  EPtot      =        39.1266
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =         0.0000  RESTRAINT  =         0.0000
        # Compare OpenMM energies with the Amber energies (above)
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbond'], 9.2210, places=4)

    def testGB1Energy(self): # HCT (igb=1)
        """Compare OpenMM and CHAMBER GB (igb=1) energies (w/ and w/out salt)"""
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'),
                           get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(implicitSolvent=app.HCT)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -74.0748  EKtot   =         0.0000  EPtot      =       -74.0748
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -113.2014  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbond'], -103.9804, places=4)
        system = parm.createSystem(implicitSolvent=app.HCT,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -74.4173  EKtot   =         0.0000  EPtot      =       -74.4173
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -113.5439  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbond'], -104.3229, places=4)

    def testGB2Energy(self): # OBC1 (igb=2)
        """Compare OpenMM and CHAMBER GB (igb=2) energies (w/ and w/out salt)"""
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'),
                           get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(implicitSolvent=app.OBC1)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -77.9618  EKtot   =         0.0000  EPtot      =       -77.9618
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -117.0885  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbond'], -107.8675, places=4)
        system = parm.createSystem(implicitSolvent=app.OBC1,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -78.3072  EKtot   =         0.0000  EPtot      =       -78.3072
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -117.4339  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbond'], -108.2129, places=4)

    def testGB5Energy(self): # OBC2 (igb=5)
        """Compare OpenMM and CHAMBER GB (igb=5) energies (w/ and w/out salt)"""
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'),
                           get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(implicitSolvent=app.OBC2)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -73.7129  EKtot   =         0.0000  EPtot      =       -73.7129
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -112.8396  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbond'], -103.6186, places=4)
        system = parm.createSystem(implicitSolvent=app.OBC2,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -74.0546  EKtot   =         0.0000  EPtot      =       -74.0546
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -113.1813  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbond'], -103.9603, places=4)

    def testGB7Energy(self): # GBn (igb=7)
        """Compare OpenMM and CHAMBER GB (igb=7) energies (w/ and w/out salt)"""
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'),
                           get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        PT.changeRadii(parm, 'mbondi3').execute() # Need new radius set
        PT.loadRestrt(parm, get_fn('ala_ala_ala.rst7')).execute()
        system = parm.createSystem(implicitSolvent=app.GBn)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -75.0550  EKtot   =         0.0000  EPtot      =       -75.0550
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -114.1816  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbond'], -104.9606, places=4)
        system = parm.createSystem(implicitSolvent=app.GBn,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -75.3985  EKtot   =         0.0000  EPtot      =       -75.3985
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -114.5251  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbond'], -105.3041, places=4)

    def testGB8Energy(self): # GBn2 (igb=8)
        """Compare OpenMM and CHAMBER GB (igb=8) energies (w/ and w/out salt)"""
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'),
                           get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        PT.changeRadii(parm, 'mbondi3').execute() # Need new radius set
        PT.loadRestrt(parm, get_fn('ala_ala_ala.rst7')).execute()
        system = parm.createSystem(implicitSolvent=app.GBn2)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -78.2339  EKtot   =         0.0000  EPtot      =       -78.2339
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -117.3606  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbond'], -108.1396, places=4)
        system = parm.createSystem(implicitSolvent=app.GBn2,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -78.5802  EKtot   =         0.0000  EPtot      =       -78.5802
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -117.7068  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbond'], -108.4858, places=4)

    def testRst7(self):
        """ Test using OpenMMRst7 to provide coordinates (CHAMBER) """
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'),
                           get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem() # Default, no cutoff
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        rst = Rst7.open(get_fn('ala_ala_ala.rst7')).positions
        sim.context.setPositions(rst)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =        39.1266  EKtot   =         0.0000  EPtot      =        39.1266
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =         0.0000  RESTRAINT  =         0.0000
        # Compare OpenMM energies with the Amber energies (above)
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbond'], 9.2210, places=4)

    def testPME(self):
        """ Compare OpenMM and CHAMBER PME energies """
        parm = ChamberParm(get_fn('ala3_solv.parm7'), get_fn('ala3_solv.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#Bond         =            1.1324222     Angle        =            1.0688008
#Dihedral     =            7.8114302     Urey-Bradley =            0.0614241
#Improper     =            0.0000000     CMAP         =            0.1267899
#Nonbond      =         6514.4460318
#TOTAL        =         6524.6468990
        self.assertAlmostEqual(energies['bond'], 1.1324, delta=5e-3)
        self.assertAlmostEqual(energies['angle'], 1.0688, delta=5e-3)
        self.assertAlmostEqual(energies['urey'], 0.06142, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 7.8114, delta=5e-3)
        self.assertAlmostEqual(energies['improper'], 0)
        self.assertRelativeEqual(energies['cmap'], 0.12679, places=3)
        self.assertRelativeEqual(energies['nonbond'], 6514.4460, places=3)

    def testDispersionCorrection(self):
        """ Compare OpenMM and CHAMBER energies without vdW correction """
        parm = ChamberParm(get_fn('ala3_solv.parm7'), get_fn('ala3_solv.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5)
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setUseDispersionCorrection(False)
            elif isinstance(force, mm.CustomNonbondedForce):
                force.setUseLongRangeCorrection(False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#Bond         =            1.1324222     Angle        =            1.0688008
#Dihedral     =            7.8114302     Urey-Bradley =            0.0614241
#Improper     =            0.0000000     CMAP         =            0.1267899
#Nonbond      =         6584.1603528
#TOTAL        =         6594.3612201
        self.assertAlmostEqual(energies['bond'], 1.13242, delta=5e-5)
        self.assertAlmostEqual(energies['angle'], 1.0688, delta=5e-3)
        self.assertAlmostEqual(energies['urey'], 0.06142, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 7.81143, delta=5e-3)
        self.assertAlmostEqual(energies['improper'], 0, delta=5e-4)
        self.assertRelativeEqual(energies['cmap'], 0.12679, places=3)
        self.assertRelativeEqual(energies['nonbond'], 6584.1604, delta=5e-4)

    def testSHAKE(self):
        """ Compare OpenMM and CHAMBER PME energies excluding SHAKEn bonds """
        parm = ChamberParm(get_fn('ala3_solv.parm7'), get_fn('ala3_solv.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5,
                                   flexibleConstraints=False,
                                   constraints=app.HBonds)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        # The only thing that changes here compared to the other periodic tests
        # is the bond energy, which should be slightly smaller than before
        state = sim.context.getState(getEnergy=True, enforcePeriodicBox=True,
                                     groups=2**parm.BOND_FORCE_GROUP)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.13236, delta=5e-5)
        self.assertAlmostEqual(energies['angle'], 1.0688, delta=5e-3)
        self.assertAlmostEqual(energies['urey'], 0.06142, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 7.81143, delta=5e-3)
        self.assertAlmostEqual(energies['improper'], 0, delta=5e-4)
        self.assertRelativeEqual(energies['cmap'], 0.12679, places=3)
        self.assertRelativeEqual(energies['nonbond'], 6514.4460, places=3)

    @unittest.skipIf(utils.skip_big_tests(), "Skipping OMM tests on large systems")
    def testBigPME(self):
        """ Compare OpenMM and CHAMBER PME energies on big system """
        parm = ChamberParm(get_fn('dhfr_cmap_pbc.parm7'),
                           get_fn('dhfr_cmap_pbc.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =   -228093.7288  EKtot   =         0.0000  EPtot      =   -228093.7288
#BOND   =      8582.2932  ANGLE   =      5019.3573  DIHED      =       740.9529
#UB     =        29.6502  IMP     =        14.2581  CMAP       =      -216.2510
#1-4 NB =       345.6860  1-4 EEL =      6475.4350  VDWAALS    =     28377.6578
#EELEC  =   -277462.7684  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#Ewald error estimate:   0.3355E-03
        self.assertAlmostEqual(energies['bond'], 8582.2932, delta=5e-3)
        self.assertAlmostEqual(energies['angle'], 5019.3573, delta=5e-3)
        self.assertAlmostEqual(energies['urey'], 29.6502, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 740.9529, delta=5e-3)
        self.assertAlmostEqual(energies['improper'], 14.2581, delta=5e-4)
        self.assertRelativeEqual(energies['cmap'], -216.2510, places=3)
        self.assertRelativeEqual(energies['nonbond'], -242263.9896, places=3)

    @unittest.skipIf(utils.skip_big_tests(), "Skipping OMM tests on large systems")
    def testBigDispersionCorrection(self):
        """ Compare OpenMM and CHAMBER w/out vdW corr on big system """
        parm = ChamberParm(get_fn('dhfr_cmap_pbc.parm7'),
                           get_fn('dhfr_cmap_pbc.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5)
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setUseDispersionCorrection(False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = decomposed_energy(sim.context, parm)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =   -226511.4095  EKtot   =         0.0000  EPtot      =   -226511.4095
#BOND   =      8582.2932  ANGLE   =      5019.3573  DIHED      =       740.9529
#UB     =        29.6502  IMP     =        14.2581  CMAP       =      -216.2510
#1-4 NB =       345.6860  1-4 EEL =      6475.4350  VDWAALS    =     29959.9772
#EELEC  =   -277462.7684  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#Ewald error estimate:   0.3355E-03
        self.assertAlmostEqual(energies['bond'], 8582.2932, delta=5e-3)
        self.assertAlmostEqual(energies['angle'], 5019.3573, delta=5e-3)
        self.assertAlmostEqual(energies['urey'], 29.6502, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 740.9529, delta=5e-3)
        self.assertAlmostEqual(energies['improper'], 14.2581, delta=5e-4)
        self.assertRelativeEqual(energies['cmap'], -216.2510, places=3)
        self.assertRelativeEqual(energies['nonbond'], -240681.6702, places=4)

    @unittest.skipIf(utils.skip_big_tests(), "Skipping OMM tests on large systems")
    def testBigSHAKE(self):
        """ Compare OpenMM and CHAMBER PME excluding SHAKEn bonds (big) """
        parm = ChamberParm(get_fn('dhfr_cmap_pbc.parm7'),
                           get_fn('dhfr_cmap_pbc.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5,
                                   flexibleConstraints=False,
                                   constraints=app.HBonds)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        # The only thing that changes here compared to the other periodic tests
        # is the bond energy, which should be slightly smaller than before
        state = sim.context.getState(getEnergy=True, enforcePeriodicBox=True,
                                     groups=2**parm.BOND_FORCE_GROUP)
        bond = state.getPotentialEnergy().value_in_unit(u.kilocalories_per_mole)
        self.assertAlmostEqual(bond, 139.2453, delta=5e-4)

    def testInterfacePBC(self):
        """ Testing all OpenMMChamberParm.createSystem options (periodic) """
        parm = ChamberParm(get_fn('ala3_solv.parm7'),
                           get_fn('ala3_solv.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=10.0*u.angstroms,
                                   constraints=None, rigidWater=False,
                                   removeCMMotion=True,
                                   ewaldErrorTolerance=1e-5)
        has_cmmotion = False
        for f in system.getForces():
            has_cmmotion = has_cmmotion or isinstance(f, mm.CMMotionRemover)
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getCutoffDistance(), 10.0*u.angstroms)
                self.assertEqual(f.getNonbondedMethod(), 4)
                self.assertFalse(f.getUseSwitchingFunction())
                self.assertEqual(f.getEwaldErrorTolerance(), 1e-5)
            elif isinstance(f, mm.HarmonicBondForce):
                # HarmonicBondForce is either bond or Urey-Bradley
                self.assertTrue(f.getNumBonds() in
                        [parm.parm_data['CHARMM_UREY_BRADLEY_COUNT'][0],
                         parm.ptr('nbona')+parm.ptr('nbonh')]
                )
            elif isinstance(f, mm.HarmonicAngleForce):
                self.assertEqual(f.getNumAngles(),
                                 parm.ptr('ntheth')+parm.ptr('ntheta'))
            elif isinstance(f, mm.PeriodicTorsionForce):
                self.assertEqual(f.getNumTorsions(),
                                 parm.ptr('nphih')+parm.ptr('nphia'))
        self.assertTrue(has_cmmotion)
        self.assertEqual(system.getNumConstraints(), 0)
        system = parm.createSystem(nonbondedMethod=app.Ewald,
                                   nonbondedCutoff=30.0*u.angstroms,
                                   switchDistance=10.0*u.angstroms,
                                   constraints=app.HBonds, rigidWater=True,
                                   removeCMMotion=False)
        has_cmmotion = False
        for f in system.getForces():
            has_cmmotion = has_cmmotion or isinstance(f, mm.CMMotionRemover)
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getCutoffDistance(), 30.0*u.angstroms)
                self.assertTrue(f.getUseSwitchingFunction())
                self.assertEqual(f.getSwitchingDistance(), 10.0*u.angstroms)
                self.assertEqual(f.getNonbondedMethod(), 3)
                self.assertEqual(f.getEwaldErrorTolerance(), 5e-4)
            elif isinstance(f, mm.HarmonicBondForce):
                # HarmonicBondForce is either bond or Urey-Bradley
                self.assertTrue(f.getNumBonds() in
                        [parm.parm_data['CHARMM_UREY_BRADLEY_COUNT'][0],
                         parm.ptr('nbona')+parm.ptr('nbonh')]
                )
            elif isinstance(f, mm.HarmonicAngleForce):
                self.assertEqual(f.getNumAngles(),
                                 parm.ptr('ntheta')+parm.ptr('ntheth'))
            elif isinstance(f, mm.PeriodicTorsionForce):
                self.assertEqual(f.getNumTorsions(),
                                 parm.ptr('nphih')+parm.ptr('nphia'))
        self.assertFalse(has_cmmotion)
        self.assertEqual(system.getNumConstraints(), parm.ptr('nbonh'))
        system = parm.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                   nonbondedCutoff=20.0*u.angstroms,
                                   constraints=app.AllBonds, rigidWater=True,
                                   flexibleConstraints=False)
        has_cmmotion = False
        for f in system.getForces():
            has_cmmotion = has_cmmotion or isinstance(f, mm.CMMotionRemover)
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getCutoffDistance(), 20.0*u.angstroms)
                self.assertEqual(f.getNonbondedMethod(), 2)
            elif isinstance(f, mm.HarmonicBondForce):
                # HarmonicBondForce is either bond or Urey-Bradley
                self.assertTrue(f.getNumBonds() in
                        [parm.parm_data['CHARMM_UREY_BRADLEY_COUNT'][0], 0]
                )
            elif isinstance(f, mm.HarmonicAngleForce):
                self.assertEqual(f.getNumAngles(),
                                 parm.ptr('ntheta')+parm.ptr('ntheth'))
            elif isinstance(f, mm.PeriodicTorsionForce):
                self.assertEqual(f.getNumTorsions(),
                                 parm.ptr('nphih')+parm.ptr('nphia'))
        self.assertTrue(has_cmmotion)
        self.assertEqual(system.getNumConstraints(),
                         parm.ptr('nbonh')+parm.ptr('nbona'))
        system = parm.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                   nonbondedCutoff=20.0*u.angstroms,
                                   constraints=app.HBonds, rigidWater=True,
                                   flexibleConstraints=False,
                                   hydrogenMass=3.00*u.daltons)
        has_cmmotion = False
        for f in system.getForces():
            has_cmmotion = has_cmmotion or isinstance(f, mm.CMMotionRemover)
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getCutoffDistance(), 20.0*u.angstroms)
                self.assertEqual(f.getNonbondedMethod(), 2)
            elif isinstance(f, mm.HarmonicBondForce):
                # HarmonicBondForce is either bond or Urey-Bradley
                self.assertTrue(f.getNumBonds() in
                        [parm.parm_data['CHARMM_UREY_BRADLEY_COUNT'][0],
                         parm.ptr('nbona')]
                )
            elif isinstance(f, mm.HarmonicAngleForce):
                self.assertEqual(f.getNumAngles(),
                                 parm.ptr('ntheta')+parm.ptr('ntheth'))
            elif isinstance(f, mm.PeriodicTorsionForce):
                self.assertEqual(f.getNumTorsions(),
                                 parm.ptr('nphih')+parm.ptr('nphia'))
        self.assertTrue(has_cmmotion)
        self.assertEqual(system.getNumConstraints(), parm.ptr('nbonh'))
        totmass = 0
        for i, atom in enumerate(parm.atoms):
            mass = system.getParticleMass(i).value_in_unit(u.daltons)
            totmass += mass
            if atom.element == 1:
                self.assertAlmostEqual(mass, 3)
        self.assertAlmostEqual(totmass, sum(parm.parm_data['MASS']), places=6)
        # Trap some illegal options
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod=0))
        self.assertRaises(ValueError, lambda: parm.createSystem(constraints=0))

    def testInterfaceNoPBC(self):
        """Testing all OpenMMChamberParm.createSystem options (non-periodic)"""
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'),
                           get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   constraints=app.HBonds,
                                   implicitSolvent=app.HCT,
                                   soluteDielectric=2.0,
                                   solventDielectric=80.0,
                                   flexibleConstraints=False)
        has_cmmotion = False
        for f in system.getForces():
            has_cmmotion = has_cmmotion or isinstance(f, mm.CMMotionRemover)
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)
        self.assertTrue(has_cmmotion)
        system = parm.createSystem(nonbondedMethod=app.CutoffNonPeriodic,
                                   nonbondedCutoff=100*u.angstroms,
                                   implicitSolvent=app.HCT,
                                   implicitSolventSaltConc=0.1*u.molar)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 1)
                self.assertEqual(f.getCutoffDistance(), 100*u.angstroms)
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.OBC1,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.OBC1,
                                   implicitSolventKappa=0.8)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.OBC2,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.OBC2,
                                   implicitSolventKappa=0.8)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn,
                                   implicitSolventSaltConc=10.0,
                                   implicitSolventKappa=0.8) # kappa is priority
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn2,
                                   implicitSolventSaltConc=10.0,
                                   implicitSolventKappa=0.8) # kappa is priority
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn2,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn2,
                                   solventDielectric=80, soluteDielectric=4,
                                   temperature=400.0*u.kelvin,
                                   implicitSolventSaltConc=0.1*u.molar)
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)
        # Test some illegal options
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod=app.PME))
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod=app.CutoffPeriodic))
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod=app.Ewald))

if not has_openmm:
    # If we don't have OpenMM installed, delete all of the test classes so we
    # don't try to test anything (nor does it count in the final tally, giving
    # us a sense of false achievement)
    del TestAmberParm, TestChamberParm

if __name__ == '__main__':
    unittest.main()
