"""
Contains unittests for running OpenMM calculations using the Amber file parsers
"""
from __future__ import division, print_function, absolute_import

from collections import defaultdict
from copy import copy
from math import sqrt
import os
from parmed import Structure, topologyobjects, load_file
from parmed.exceptions import ParameterError
from parmed.amber import AmberParm, ChamberParm, Rst7
from parmed.openmm import load_topology, energy_decomposition_system
import parmed.unit as u
from parmed.utils.six.moves import range, zip
import parmed.tools as PT
import sys
import unittest
from utils import (get_fn, CPU, mm, app, has_openmm, FileIOTestCase,
        TestCaseRelative, get_saved_fn, run_all_tests, QuantityTestCase)

# OpenMM NonbondedForce methods are enumerated values. From NonbondedForce.h,
# they are:
#   0 - NoCutoff
#   1 - CutoffNonPeriodic
#   2 - CutoffPeriodic
#   3 - Ewald
#   4 - PME

def energy_decomposition(parm, context):
    from parmed.openmm.utils import energy_decomposition
    ret = defaultdict(float)
    for key, val in energy_decomposition(parm, context).items():
        ret[key] = val
    return ret

@unittest.skipUnless(has_openmm, 'Cannot test without OpenMM')
class TestAmberParm(FileIOTestCase, TestCaseRelative, QuantityTestCase):

    def test_ep_energy(self):
        """ Tests AmberParm handling of extra points in TIP4P water """
        parm = AmberParm(self.get_fn('tip4p.parm7'), self.get_fn('tip4p.rst7'))
        repr(parm) # Make sure it doesn't crash
        parm.join_dihedrals() # Make sure join_dihedrals doesn't do anything w/ no dihedrals
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   constraints=app.HBonds,
                                   rigidWater=True,
                                   flexibleConstraints=False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
# Etot   =     -1756.2018  EKtot   =       376.7454  EPtot      =     -2132.9472
# BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
# 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =       378.4039
# EELEC  =     -2511.3511  EHBOND  =         0.0000  RESTRAINT  =         0.0000
# EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =      6514.3661
        self.assertAlmostEqual(energies['bond'], 0)
        self.assertAlmostEqual(energies['angle'], 0)
        self.assertAlmostEqual(energies['dihedral'], 0)
        self.assertRelativeEqual(energies['nonbonded'], -2133.2963, places=4)
        # Check that we have the correct number of virtual sites
        nvirt = 0
        for i in range(system.getNumParticles()):
            nvirt += system.isVirtualSite(i)
        self.assertEqual(parm.ptr('NUMEXTRA'), nvirt)
        # Now test the forces to make sure that they are computed correctly in
        # the presence of extra points
        pstate = sim.context.getState(getForces=True)
        sstate = mm.XmlSerializer.deserialize(
                open(get_saved_fn('tip4pforces.xml'), 'r').read()
        )
        pf = pstate.getForces().value_in_unit(u.kilocalorie_per_mole/u.angstrom)
        sf = sstate.getForces().value_in_unit(u.kilocalorie_per_mole/u.angstrom)
        # Error checking on omm_set_virtual_sites
        self.assertRaises(ValueError, lambda:
                parm.omm_set_virtual_sites(mm.System()))

        for p, s in zip(pf, sf):
            for x1, x2 in zip(p, s):
                # Compare large forces relatively and small ones absolutely
                if abs(x1) > 1 or abs(x2) > 1:
                    self.assertRelativeEqual(x1, x2, places=2)
                else:
                    self.assertAlmostEqual(x1, x2, delta=2e-2)

    @unittest.skipUnless(run_all_tests, "Skipping long tests")
    def test_round_trip_ep(self):
        """ Test ParmEd -> OpenMM round trip with Amber EPs and PME """
        parm = AmberParm(self.get_fn('tip4p.parm7'), self.get_fn('tip4p.rst7'))
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
        e1 = energy_decomposition(parm, con1)
        e2 = energy_decomposition(parm, con2)
        self.assertAlmostEqual(e1['bond'], e2['bond'])
        self.assertAlmostEqual(e1['angle'], e2['angle'])
        self.assertAlmostEqual(e1['dihedral'], e2['dihedral'])
        self.assertAlmostEqual(e1['nonbonded'], e2['nonbonded'], places=5)
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

    def test_ep_energy2(self):
        """ Tests AmberParm handling of extra points in TIP5P water """
        parm = AmberParm(self.get_fn('tip5p.parm7'), self.get_fn('tip5p.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   constraints=app.HBonds,
                                   rigidWater=True,
                                   flexibleConstraints=False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 0)
        self.assertAlmostEqual(energies['angle'], 0)
        self.assertAlmostEqual(energies['dihedral'], 0)
        self.assertRelativeEqual(energies['nonbonded'], -2142.418956, places=4)
        # Check that we have the correct number of virtual sites
        nvirt = 0
        for i in range(system.getNumParticles()):
            nvirt += system.isVirtualSite(i)
        self.assertEqual(parm.ptr('NUMEXTRA'), nvirt)
        # Now test the forces to make sure that they are computed correctly in
        # the presence of extra points
        pstate = sim.context.getState(getForces=True)
        sstate = mm.XmlSerializer.deserialize(
                open(get_saved_fn('tip5pforces.xml'), 'r').read()
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

    def test_gas_energy(self):
        """ Compare Amber and OpenMM gas phase energies """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem() # Default, no cutoff
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =         2.4544  EKtot   =         0.0000  EPtot      =         2.4544
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =         0.0000  RESTRAINT  =         0.0000
        # Compare OpenMM energies with the Amber energies (above)
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -30.2355, places=3)

    def test_improper_dihedrals(self):
        """ Test splitting improper/proper dihedrals in AmberParm """
        to_delete = []
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        for i, d in enumerate(parm.dihedrals):
            if d.improper:
                continue
            to_delete.append(i)
        for i in reversed(to_delete):
            parm.dihedrals[i].delete()
            del parm.dihedrals[i]
        self.assertGreater(len(parm.dihedrals), 0)
        f = parm.omm_dihedral_force(split=True)
        self.assertIsInstance(f, mm.PeriodicTorsionForce)
        self.assertEqual(f.getNumTorsions(), len(parm.dihedrals))

    def test_round_trip(self):
        """ Test ParmEd -> OpenMM round trip with Amber gas phase """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem()
        system2 = load_topology(parm.topology, system).createSystem()
        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system2, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(parm.positions)
        con2.setPositions(parm.positions)
        e1 = energy_decomposition(parm, con1)
        e2 = energy_decomposition(parm, con2)
        self.assertAlmostEqual(e1['bond'], e2['bond'])
        self.assertAlmostEqual(e1['angle'], e2['angle'])
        self.assertAlmostEqual(e1['dihedral'], e2['dihedral'])
        self.assertAlmostEqual(e1['nonbonded'], e2['nonbonded'])

    def test_round_trip_xml(self):
        """ Test ParmEd -> OpenMM round trip with Amber gas phase via XML """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem()
        fname = self.get_fn('ash.xml', written=True)
        with open(fname, 'w') as f:
            f.write(mm.XmlSerializer.serialize(system))
        system2 = load_topology(parm.topology, fname).createSystem()
        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system2, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(parm.positions)
        con2.setPositions(parm.positions)
        e1 = energy_decomposition(parm, con1)
        e2 = energy_decomposition(parm, con2)
        self.assertAlmostEqual(e1['bond'], e2['bond'])
        self.assertAlmostEqual(e1['angle'], e2['angle'])
        self.assertAlmostEqual(e1['dihedral'], e2['dihedral'])
        self.assertAlmostEqual(e1['nonbonded'], e2['nonbonded'])

    def test_gb1_energy(self): # HCT (igb=1)
        """ Compare Amber and OpenMM GB (igb=1) energies (w/ and w/out salt) """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(implicitSolvent=app.HCT)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0306  EKtot   =         0.0000  EPtot      =       -24.0306
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -26.4850  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -56.7205, places=3)
        system = parm.createSystem(implicitSolvent=app.HCT,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0508  EKtot   =         0.0000  EPtot      =       -24.0508
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -26.5051  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -56.7406, places=3)

    def test_gb2_energy(self): # OBC1 (igb=2)
        """ Compare Amber and OpenMM GB (igb=2) energies (w/ and w/out salt) """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(implicitSolvent=app.OBC1)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0306  EKtot   =         0.0000  EPtot      =       -24.0306
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -28.1689  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -58.4044, places=3)
        system = parm.createSystem(implicitSolvent=app.OBC1,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0508  EKtot   =         0.0000  EPtot      =       -24.0508
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -28.1895  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -58.4250, places=3)

    def test_gb5_energy(self): # OBC2 (igb=5)
        """ Compare Amber and OpenMM GB (igb=5) energies (w/ and w/out salt) """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(implicitSolvent=app.OBC2)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0306  EKtot   =         0.0000  EPtot      =       -24.0306
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -28.4639  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -55.6994, places=3)
        system = parm.createSystem(implicitSolvent=app.OBC2,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0508  EKtot   =         0.0000  EPtot      =       -24.0508
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -25.4835  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -55.7190, places=3)

    def test_gb7_energy(self): # GBn (igb=7)
        """ Compare Amber and OpenMM GB (igb=7) energies (w/ and w/out salt) """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        PT.changeRadii(parm, 'mbondi3').execute() # Need new radius set
        PT.loadRestrt(parm, self.get_fn('ash.rst7')).execute() # Load crds into copy
        system = parm.createSystem(implicitSolvent=app.GBn)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0306  EKtot   =         0.0000  EPtot      =       -24.0306
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -21.1012  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -51.3367, places=3)
        system = parm.createSystem(implicitSolvent=app.GBn,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0508  EKtot   =         0.0000  EPtot      =       -24.0508
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -21.1194  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -51.3549, places=3)

    def test_gb8_energy(self): # GBn2 (igb=8)
        """ Compare Amber and OpenMM GB (igb=8) energies (w/ and w/out salt) """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        PT.changeRadii(parm, 'mbondi3').execute() # Need new radius set
        PT.loadRestrt(parm, self.get_fn('ash.rst7')).execute() # Load crds into copy
        system = parm.createSystem(implicitSolvent=app.GBn2)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0306  EKtot   =         0.0000  EPtot      =       -24.0306
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -23.4639  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -53.6994, places=3)
        system = parm.createSystem(implicitSolvent=app.GBn2,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0508  EKtot   =         0.0000  EPtot      =       -24.0508
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -23.4832  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 5.4435, places=4)
        self.assertRelativeEqual(energies['angle'], 2.8766, places=4)
        self.assertRelativeEqual(energies['dihedral'], 24.3697, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -53.7187, places=3)

    def test_energy_decomp_system(self):
        """ Tests the energy_decomposition_system function """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        PT.changeRadii(parm, 'mbondi3').execute() # Need new radius set
        system = parm.createSystem(implicitSolvent=app.GBn2)
        energies = energy_decomposition_system(parm, system)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -24.0306  EKtot   =         0.0000  EPtot      =       -24.0306
#BOND   =         5.4435  ANGLE   =         2.8766  DIHED      =        24.3697
#1-4 NB =         6.1446  1-4 EEL =        20.8049  VDWAALS    =        44.3715
#EELEC  =      -101.5565  EGB     =       -23.4639  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies[0][1], 5.4435, places=4)
        self.assertRelativeEqual(energies[1][1], 2.8766, places=4)
        self.assertRelativeEqual(energies[2][1], 24.3697, places=4)
        self.assertRelativeEqual(energies[3][1], -30.238225, places=3)
        self.assertRelativeEqual(energies[4][1], -23.464687, places=3)
        self.assertEqual(energies[5][1], 0)
        energies = energy_decomposition_system(parm, system, platform='CPU')
        self.assertRelativeEqual(energies[0][1], 5.4435, places=4)
        self.assertRelativeEqual(energies[1][1], 2.8766, places=4)
        self.assertRelativeEqual(energies[2][1], 24.3697, places=4)
        self.assertRelativeEqual(energies[3][1], -30.238225, places=3)
        self.assertRelativeEqual(energies[4][1], -23.464687, places=3)
        self.assertEqual(energies[5][1], 0)
        parm = AmberParm(self.get_fn('solv2.parm7'), self.get_fn('solv2.rst7'))
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms)
        energies = energy_decomposition_system(parm, system, platform='CPU')

    @unittest.skipUnless(run_all_tests, "Skipping OMM tests on large systems")
    def test_ewald(self):
        """ Compare Amber and OpenMM Ewald energies """
        parm = AmberParm(self.get_fn('solv2.parm7'), self.get_fn('solv2.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.Ewald,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertRelativeEqual(energies['bond'], 12.2920213, 4)
        self.assertRelativeEqual(energies['angle'], 32.3453097, 5)
        self.assertRelativeEqual(energies['dihedral'], 96.0811552, 5)
        self.assertRelativeEqual(energies['nonbonded'], -12926.394844, 4)

    def test_pme(self):
        """ Compare Amber and OpenMM PME energies """
        parm = AmberParm(self.get_fn('solv2.parm7'), self.get_fn('solv2.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
# Etot   =    -12785.6764  EKtot   =         0.0000  EPtot      =    -12785.6764
# BOND   =        12.2920  ANGLE   =        32.3453  DIHED      =        96.0812
# 1-4 NB =        39.1460  1-4 EEL =       420.5797  VDWAALS    =      2836.7832
# EELEC  =    -16222.9037  EHBOND  =         0.0000  RESTRAINT  =         0.0000
        self.assertRelativeEqual(energies['bond'], 12.2920213, 4)
        self.assertRelativeEqual(energies['angle'], 32.3453097, 5)
        self.assertRelativeEqual(energies['dihedral'], 96.0811552, 5)
        self.assertRelativeEqual(energies['nonbonded'], -12926.394844, 4)

    def test_dispersion_correction(self):
        """ Compare Amber and OpenMM PME energies w/out vdW correction """
        parm = AmberParm(self.get_fn('solv2.parm7'), self.get_fn('solv2.rst7'))
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
        energies = energy_decomposition(parm, sim.context)
# Etot   =    -12682.5312  EKtot   =         0.0000  EPtot      =    -12682.5312
# BOND   =        12.2920  ANGLE   =        32.3453  DIHED      =        96.0812
# 1-4 NB =        39.1460  1-4 EEL =       420.5797  VDWAALS    =      2939.9284
# EELEC  =    -16222.9037  EHBOND  =         0.0000  RESTRAINT  =         0.0000
        lrc = 2836.7832 - 2939.9284
        self.assertRelativeEqual(energies['bond'], 12.2920213, 4)
        self.assertRelativeEqual(energies['angle'], 32.3453097, 5)
        self.assertRelativeEqual(energies['dihedral'], 96.0811552, 5)
        self.assertRelativeEqual(energies['nonbonded'], -12926.394844-lrc, 4)

    def test_shake(self):
        """ Compare Amber and OpenMM PME energies excluding SHAKEn bonds """
        parm = AmberParm(self.get_fn('solv2.parm7'), self.get_fn('solv2.rst7'))
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
        self.assertAlmostEqual(bond, 12.134943963951512, places=4)

    def test_nbfix(self):
        """ Compare Amber and OpenMM PME energies with NBFIX modifications """
        # For now, long-range correction is not available
        parm = AmberParm(self.get_fn('ff14ipq.parm7'), self.get_fn('ff14ipq.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        self.assertTrue(parm.has_NBFIX())
        PT.change(parm, 'CHARGE', ':*', 0).execute() # only check LJ energies
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
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
        self.assertAlmostEqual(energies['nonbonded'], 1200.11736916, 3)
        # Check 12-6-4 potential with NBFIX turned on. This was validated
        # against the results generated by sander *without* the long-range
        # correction (but that is kept on for this test -- there is no direct
        # comparison in Amber since the L-R correction is not coded for 12-6-4,
        # but the CustomNonbondedForce L-R correction and the 12-6-4+NBFIX
        # implementations have been separately validated, so this is fine)
        parm.add_flag('LENNARD_JONES_CCOEF',
                      str(parm.formats['LENNARD_JONES_ACOEF']),
                      num_items=len(parm.parm_data['LENNARD_JONES_ACOEF']))
        parm.parm_data['LENNARD_JONES_CCOEF'][0] = 100
        parm.parm_data['LENNARD_JONES_CCOEF'][2] = 100
        parm.parm_data['LENNARD_JONES_CCOEF'][4] = 100
        parm.parm_data['LENNARD_JONES_CCOEF'][6] = 100
        parm.parm_data['LENNARD_JONES_CCOEF'][8] = 100
        parm.parm_data['LENNARD_JONES_CCOEF'][9] = 100
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 0.0654, delta=0.0002)
        self.assertAlmostEqual(energies['angle'], 0.9616, 4)
        self.assertAlmostEqual(energies['dihedral'], -5.4917, 4)
        self.assertAlmostEqual(energies['nonbonded'], 1168.06630486, 3)
        # Now make sure that if we use a switching function, that gets
        # propagated to the CustomNonbondedForce
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   switchDistance=6*u.angstroms)
        n = 0
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce) or isinstance(force,
                    mm.CustomNonbondedForce):
                self.assertTrue(force.getUseSwitchingFunction())
                self.assertAlmostEqualQuantities(force.getSwitchingDistance(),
                                                 6*u.angstroms)
                n += 1
        self.assertEqual(n, 2)

    def test_nbfix_exceptions(self):
        """ Tests exception transfers when NBFIX present """
        s = Structure()
        at1 = topologyobjects.AtomType('A1', 1, 12.01, atomic_number=6)
        at2 = topologyobjects.AtomType('A2', 2, 12.01, atomic_number=6)
        at3 = topologyobjects.AtomType('A3', 3, 12.01, atomic_number=6)
        at4 = topologyobjects.AtomType('A4', 4, 12.01, atomic_number=6)
        at1.set_lj_params(0.5, 1.0)
        at2.set_lj_params(0.6, 1.1)
        at3.set_lj_params(0.7, 1.2)
        at4.set_lj_params(0.8, 1.3)

        # Add some NBFIXes
        at1.add_nbfix('A4', 0.0, 0.0)
        at4.add_nbfix('A1', 0.0, 0.0)
        at3.add_nbfix('A2', 0.9, 4.1)
        at2.add_nbfix('A3', 0.9, 4.1)

        # Add a handful of atoms
        s.add_atom(topologyobjects.Atom(name='AA', type='A1', atomic_number=6), 'LJR', 1)
        s.add_atom(topologyobjects.Atom(name='AB', type='A2', atomic_number=6), 'LJR', 2)
        s.add_atom(topologyobjects.Atom(name='AC', type='A3', atomic_number=6), 'LJR', 3)
        s.add_atom(topologyobjects.Atom(name='AD', type='A4', atomic_number=6), 'LJR', 4)

        # Assign the types
        s[0].atom_type = at1
        s[1].atom_type = at2
        s[2].atom_type = at3
        s[3].atom_type = at4

        # Add bonds, angles, and the one dihedral
        s.bond_types.append(topologyobjects.BondType(0, 1.0))
        s.angle_types.append(topologyobjects.AngleType(0, 90))
        s.dihedral_types.append(topologyobjects.DihedralType(0, 1, 0))
        s.bond_types.claim()
        s.angle_types.claim()
        s.dihedral_types.claim()

        s.bonds.append(topologyobjects.Bond(s[0], s[1], type=s.bond_types[0]))
        s.bonds.append(topologyobjects.Bond(s[1], s[2], type=s.bond_types[0]))
        s.bonds.append(topologyobjects.Bond(s[2], s[3], type=s.bond_types[0]))

        s.angles.append(topologyobjects.Angle(s[0], s[1], s[2], type=s.angle_types[0]))
        s.angles.append(topologyobjects.Angle(s[1], s[2], s[3], type=s.angle_types[0]))

        s.dihedrals.append(
                topologyobjects.Dihedral(s[0], s[1], s[2], s[3],
                                         type=s.dihedral_types[0])
        )

        self.assertTrue(s.has_NBFIX())
        s[0].xx, s[0].xy, s[0].xz = 0, 0, 0
        s[1].xx, s[1].xy, s[1].xz = 0, 1, 0
        s[2].xx, s[2].xy, s[2].xz = 1, 1, 0
        s[3].xx, s[3].xy, s[3].xz = 1, 1, 1

        # Convert to Amber topology file
        parm = AmberParm.from_structure(s)
        self.assertTrue(parm.has_NBFIX())
        system = parm.createSystem()

        # The energy should be zero
        context = mm.Context(system, mm.VerletIntegrator(1*u.femtoseconds), CPU)
        context.setPositions(parm.positions)
        e = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
                u.kilocalories_per_mole)
        self.assertAlmostEqual(e, 0)

    def test_1264(self):
        """ Testing the 12-6-4 LJ potential in OpenMM """
        parm = AmberParm(self.get_fn('znf_1264.prmtop'), self.get_fn('znf.rst'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.NoCutoff)
        for force in system.getForces():
            if isinstance(force, mm.CustomNonbondedForce):
                self.assertTrue(force.getUseLongRangeCorrection())
                force.setUseLongRangeCorrection(False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 26.3947079, places=3)
        self.assertAlmostEqual(energies['angle'], 122.8243431, places=3)
        self.assertAlmostEqual(energies['dihedral'], 319.0419347, places=3)
        self.assertAlmostEqual(energies['nonbonded'], -2133.6170786, delta=2e-3)

    def test_1012(self):
        """ Testing the 10-12 LJ H-bond potential in OpenMM """
        parm = AmberParm(self.get_fn('ff91.parm7'), self.get_fn('ff91.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms)
        for force in system.getForces():
            if isinstance(force, mm.CustomNonbondedForce):
                self.assertTrue(force.getUseLongRangeCorrection())
        integrator = mm.VerletIntegrator(1*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 0.9675961, places=3)
        self.assertAlmostEqual(energies['angle'], 82.5853211, places=3)
        self.assertAlmostEqual(energies['dihedral'], 1.2476012, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -11194.4588654, delta=3e-5)
        # Now test with the 12-6-4 potential. sander does not implement this
        # correctly (or perhaps there is no clear "correct" implementation). But
        # what is done here is that a C-coefficient is counted for EVERY
        # interaction, whether or not there is a 10-12 definition for that pair.
        # sander does not do the r^-4 part for a 10-12 pair, but the 10-12 and
        # r^-4 ixns were developed 2+ decades apart, so using them together was
        # never considered. This test is here just to ensure coverage and
        # provide some kind of documentation about intended behavior.
        for a in parm.atoms: a.charge = 0.0
        parm.add_flag('LENNARD_JONES_CCOEF',
                      str(parm.formats['LENNARD_JONES_ACOEF']),
                      num_items=len(parm.parm_data['LENNARD_JONES_ACOEF']))
        parm.parm_data['LENNARD_JONES_CCOEF'][0] = 100
        parm.parm_data['LENNARD_JONES_CCOEF'][2] = 100
        parm.parm_data['LENNARD_JONES_CCOEF'][4] = 100
        parm.parm_data['LENNARD_JONES_CCOEF'][6] = 100
        parm.parm_data['LENNARD_JONES_CCOEF'][8] = 100
        parm.parm_data['LENNARD_JONES_CCOEF'][9] = 100
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms)
        integrator = mm.VerletIntegrator(1*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 0.9675961, places=3)
        self.assertAlmostEqual(energies['angle'], 82.5853211, places=3)
        self.assertAlmostEqual(energies['dihedral'], 1.2476012, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -48243.2239173, places=5)

    def test_hangle_constraints(self):
        """ Tests that HAngle constraints get applied correctly """
        # This used to be a bug
        # Constrain just bonds. Make sure 1000 angles remain, and we have 2000
        # constraints
        parm = AmberParm(self.get_fn('gaffwat.parm7'))
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

    def test_interface_pbc(self):
        """ Testing all AmberParm.createSystem options (periodic) """
        parm = AmberParm(self.get_fn('solv2.parm7'), self.get_fn('solv2.rst7'))
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
        # Now do rigid water and not flexible constraints, but no other constraints
        system = parm.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                   nonbondedCutoff=20*u.angstroms,
                                   constraints=None, rigidWater=True,
                                   flexibleConstraints=False)
        for f in system.getForces():
            if isinstance(f, mm.HarmonicBondForce):
                self.assertEqual(f.getNumBonds(), 119)
                break
        else:
            assert False, 'Should not be here'
        self.assertEqual(system.getNumConstraints(),
                3*sum([r.name == 'WAT' for r in parm.residues]))
        # Trap some illegal options
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod=0))
        self.assertRaises(ValueError, lambda: parm.createSystem(constraints=0))
        self.assertRaises(ValueError, lambda: parm.createSystem(hydrogenMass=-1))

    def test_dihedral_splitting(self):
        """ Tests proper splitting of torsions into proper/improper groups """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        prop, improp = parm.omm_dihedral_force(split=True)
        self.assertEqual(improp.getNumTorsions(), 5)
        self.assertEqual(prop.getNumTorsions(), len(parm.dihedrals)-5)
        self.assertEqual(parm.omm_dihedral_force(split=False).getNumTorsions(),
                         len(parm.dihedrals))

    def test_interface_no_pbc(self):
        """ Testing all AmberParm.createSystem options (non-periodic) """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        PT.changeRadii(parm, 'mbondi3').execute()
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
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                self.assertEqual(f.getNonbondedMethod(), 0)
        self.assertTrue(has_cmmotion)
        system = parm.createSystem(nonbondedMethod=app.CutoffNonPeriodic,
                                   nonbondedCutoff=100*u.angstroms,
                                   implicitSolvent=app.HCT,
                                   implicitSolventSaltConc=0.1*u.molar)
        for f in system.getForces():
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                self.assertEqual(f.getNonbondedMethod(), 1)
                self.assertEqual(f.getCutoffDistance(), 100*u.angstroms)
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.OBC1,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                self.assertEqual(f.getNonbondedMethod(), 0)
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.OBC1,
                                   implicitSolventKappa=0.8)
        for f in system.getForces():
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.OBC2,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                self.assertEqual(f.getNonbondedMethod(), 0)
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.OBC2,
                                   implicitSolventKappa=0.8)
        for f in system.getForces():
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                self.assertEqual(f.getNonbondedMethod(), 0)
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn,
                                   implicitSolventSaltConc=10.0,
                                   implicitSolventKappa=0.8) # kappa is priority
        for f in system.getForces():
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn2,
                                   implicitSolventSaltConc=10.0,
                                   implicitSolventKappa=0.8) # kappa is priority
        for f in system.getForces():
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn2,
                                   solventDielectric=80, soluteDielectric=4)
        for f in system.getForces():
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                self.assertEqual(f.getNonbondedMethod(), 0)

        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn2,
                                   solventDielectric=80, soluteDielectric=4,
                                   temperature=400.0*u.kelvin,
                                   implicitSolventSaltConc=0.1*u.molar)
        for f in system.getForces():
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                self.assertEqual(f.getNonbondedMethod(), 0)
        # Using ACE SASA method
        system = parm.createSystem(nonbondedMethod=app.NoCutoff,
                                   implicitSolvent=app.GBn2,
                                   solventDielectric=80, soluteDielectric=4,
                                   temperature=300.0*u.kelvin,
                                   implicitSolventSaltConc=0.1*u.molar,
                                   useSASA=True)
        frc = parm.omm_gbsa_force(app.GBn2 ,nonbondedCutoff=1.5,
                                  nonbondedMethod=app.CutoffNonPeriodic)
        self.assertEqual(frc.getCutoffDistance(), 1.5*u.nanometers)
        frc = parm.omm_gbsa_force(app.GBn2 ,nonbondedCutoff=1.5,
                                  nonbondedMethod=app.CutoffPeriodic)
        self.assertEqual(frc.getCutoffDistance(), 1.5*u.nanometers)
        self.assertEqual(frc.getNonbondedMethod(), frc.CutoffPeriodic)
        # Check when CustomNonbondedForce is required that the cutoff methods
        # get appropriately transferred.
        parm.parm_data['LENNARD_JONES_BCOEF'][0] = 0.0
        self.assertTrue(parm.has_NBFIX())
        system = parm.createSystem(nonbondedMethod=app.NoCutoff)
        has_custom = has_normal = False
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)
                has_normal = True
            elif isinstance(f, mm.CustomNonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 0)
                has_custom = True
        self.assertTrue(has_custom and has_normal)
        system = parm.createSystem(nonbondedMethod=app.CutoffNonPeriodic)
        has_custom = has_normal = False
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 1)
                has_normal = True
            elif isinstance(f, mm.CustomNonbondedForce):
                self.assertEqual(f.getNonbondedMethod(), 1)
                has_custom = True
        self.assertTrue(has_custom and has_normal)
        # Test some illegal options
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod=app.PME))
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod=app.CutoffPeriodic))
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod=app.Ewald))
        self.assertRaises(ValueError, lambda:
                parm.createSystem(nonbondedMethod='cat'))
        self.assertRaises(ValueError, lambda:
                parm.createSystem(implicitSolvent='No model'))

@unittest.skipUnless(has_openmm, 'Cannot test without OpenMM')
class TestChamberParm(TestCaseRelative):

    def test_gas_energy(self):
        """ Compare OpenMM and CHAMBER gas phase energies """
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'), get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem() # Default, no cutoff
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =        39.1266  EKtot   =         0.0000  EPtot      =        39.1266
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =         0.0000  RESTRAINT  =         0.0000
        # Compare OpenMM energies with the Amber energies (above)
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbonded'], 9.2210, places=4)
        # Now check that the energy functions in topologyobjects are correct
        self.assertAlmostEqual(sum(b.energy() for b in parm.bonds),
                               energies['bond'])
        self.assertAlmostEqual(sum(a.energy() for a in parm.angles),
                               energies['angle'])
        self.assertAlmostEqual(sum(a.energy() for a in parm.urey_bradleys),
                               energies['urey_bradley'])
        self.assertAlmostEqual(sum(d.energy() for d in parm.dihedrals),
                               energies['dihedral'])
        self.assertAlmostEqual(sum(i.energy() for i in parm.impropers),
                               energies['improper'])

    def test_gb1_energy(self): # HCT (igb=1)
        """Compare OpenMM and CHAMBER GB (igb=1) energies (w/ and w/out salt)"""
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'), get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(implicitSolvent=app.HCT)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -74.0748  EKtot   =         0.0000  EPtot      =       -74.0748
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -113.2014  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbonded'], -103.9804, places=4)
        system = parm.createSystem(implicitSolvent=app.HCT,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -74.4173  EKtot   =         0.0000  EPtot      =       -74.4173
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -113.5439  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbonded'], -104.3229, places=4)

    def test_gb2_energy(self): # OBC1 (igb=2)
        """Compare OpenMM and CHAMBER GB (igb=2) energies (w/ and w/out salt)"""
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'), get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(implicitSolvent=app.OBC1)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -77.9618  EKtot   =         0.0000  EPtot      =       -77.9618
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -117.0885  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbonded'], -107.8675, places=4)
        system = parm.createSystem(implicitSolvent=app.OBC1,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -78.3072  EKtot   =         0.0000  EPtot      =       -78.3072
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -117.4339  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbonded'], -108.2129, places=4)

    def test_gb5_energy(self): # OBC2 (igb=5)
        """Compare OpenMM and CHAMBER GB (igb=5) energies (w/ and w/out salt)"""
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'), get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(implicitSolvent=app.OBC2)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -73.7129  EKtot   =         0.0000  EPtot      =       -73.7129
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -112.8396  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbonded'], -103.6186, places=4)
        system = parm.createSystem(implicitSolvent=app.OBC2,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -74.0546  EKtot   =         0.0000  EPtot      =       -74.0546
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -113.1813  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbonded'], -103.9603, places=4)

    def test_gb7_energy(self): # GBn (igb=7)
        """Compare OpenMM and CHAMBER GB (igb=7) energies (w/ and w/out salt)"""
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'), get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        PT.changeRadii(parm, 'mbondi3').execute() # Need new radius set
        PT.loadRestrt(parm, get_fn('ala_ala_ala.rst7')).execute()
        system = parm.createSystem(implicitSolvent=app.GBn)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -75.0550  EKtot   =         0.0000  EPtot      =       -75.0550
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -114.1816  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbonded'], -104.9606, places=4)
        system = parm.createSystem(implicitSolvent=app.GBn,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -75.3985  EKtot   =         0.0000  EPtot      =       -75.3985
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -114.5251  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbonded'], -105.3041, places=4)

    def test_gb8_energy(self): # GBn2 (igb=8)
        """Compare OpenMM and CHAMBER GB (igb=8) energies (w/ and w/out salt)"""
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'), get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        PT.changeRadii(parm, 'mbondi3').execute() # Need new radius set
        PT.loadRestrt(parm, get_fn('ala_ala_ala.rst7')).execute()
        system = parm.createSystem(implicitSolvent=app.GBn2)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -78.2339  EKtot   =         0.0000  EPtot      =       -78.2339
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -117.3606  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbonded'], -108.1396, places=4)
        system = parm.createSystem(implicitSolvent=app.GBn2,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =       -78.5802  EKtot   =         0.0000  EPtot      =       -78.5802
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =      -117.7068  RESTRAINT  =         0.0000
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbonded'], -108.4858, places=4)

    def test_rst7(self):
        """ Test using OpenMMRst7 to provide coordinates (CHAMBER) """
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'), get_fn('ala_ala_ala.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem() # Default, no cutoff
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        rst = Rst7.open(get_fn('ala_ala_ala.rst7')).positions
        sim.context.setPositions(rst)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =        39.1266  EKtot   =         0.0000  EPtot      =        39.1266
#BOND   =         1.3351  ANGLE   =        14.1158  DIHED      =        14.2773
#UB     =         0.3669  IMP     =         0.3344  CMAP       =        -0.5239
#1-4 NB =         2.1041  1-4 EEL =       278.1614  VDWAALS    =        -1.3323
#EELEC  =      -269.7122  EGB     =         0.0000  RESTRAINT  =         0.0000
        # Compare OpenMM energies with the Amber energies (above)
        self.assertAlmostEqual(energies['bond'], 1.3351, delta=5e-4)
        self.assertAlmostEqual(energies['angle'], 14.1158, delta=5e-4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, delta=5e-4)
        self.assertAlmostEqual(energies['improper'], 0.3344, delta=5e-4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, delta=5e-4)
        self.assertRelativeEqual(energies['nonbonded'], 9.2210, places=4)

    def test_pme(self):
        """ Compare OpenMM and CHAMBER PME energies """
        parm = ChamberParm(get_fn('ala3_solv.parm7'), get_fn('ala3_solv.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#Bond         =            1.1324222     Angle        =            1.0688008
#Dihedral     =            7.8114302     Urey-Bradley =            0.0614241
#Improper     =            0.0000000     CMAP         =            0.1267899
#Nonbond      =         6514.4460318
#TOTAL        =         6524.6468990
        self.assertAlmostEqual(energies['bond'], 1.1324, delta=5e-3)
        self.assertAlmostEqual(energies['angle'], 1.0688, delta=5e-3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.06142, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 7.8114, delta=5e-3)
        self.assertAlmostEqual(energies['improper'], 0)
        self.assertRelativeEqual(energies['cmap'], 0.12679, places=3)
        self.assertRelativeEqual(energies['nonbonded'], 6514.4460, places=3)

    def test_dispersion_correction(self):
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
        energies = energy_decomposition(parm, sim.context)
#Bond         =            1.1324222     Angle        =            1.0688008
#Dihedral     =            7.8114302     Urey-Bradley =            0.0614241
#Improper     =            0.0000000     CMAP         =            0.1267899
#Nonbond      =         6584.1603528
#TOTAL        =         6594.3612201
        self.assertAlmostEqual(energies['bond'], 1.13242, delta=5e-5)
        self.assertAlmostEqual(energies['angle'], 1.0688, delta=5e-3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.06142, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 7.81143, delta=5e-3)
        self.assertAlmostEqual(energies['improper'], 0, delta=5e-4)
        self.assertRelativeEqual(energies['cmap'], 0.12679, places=3)
        self.assertRelativeEqual(energies['nonbonded'], 6584.1604, delta=5e-4)

    def test_shake(self):
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
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.13236, delta=5e-5)
        self.assertAlmostEqual(energies['angle'], 1.0688, delta=5e-3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.06142, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 7.81143, delta=5e-3)
        self.assertAlmostEqual(energies['improper'], 0, delta=5e-4)
        self.assertRelativeEqual(energies['cmap'], 0.12679, places=3)
        self.assertRelativeEqual(energies['nonbonded'], 6514.4460, places=3)

    @unittest.skipUnless(run_all_tests, "Skipping OMM tests on large systems")
    def test_big_pme(self):
        """ Compare OpenMM and CHAMBER PME energies on big system """
        parm = ChamberParm(get_fn('dhfr_cmap_pbc.parm7'), get_fn('dhfr_cmap_pbc.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,
                                   ewaldErrorTolerance=1e-5)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =   -228093.7288  EKtot   =         0.0000  EPtot      =   -228093.7288
#BOND   =      8582.2932  ANGLE   =      5019.3573  DIHED      =       740.9529
#UB     =        29.6502  IMP     =        14.2581  CMAP       =      -216.2510
#1-4 NB =       345.6860  1-4 EEL =      6475.4350  VDWAALS    =     28377.6578
#EELEC  =   -277462.7684  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#Ewald error estimate:   0.3355E-03
        self.assertAlmostEqual(energies['bond'], 8582.2932, delta=5e-3)
        self.assertAlmostEqual(energies['angle'], 5019.3573, delta=5e-3)
        self.assertAlmostEqual(energies['urey_bradley'], 29.6502, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 740.9529, delta=5e-3)
        self.assertAlmostEqual(energies['improper'], 14.2581, delta=5e-4)
        self.assertRelativeEqual(energies['cmap'], -216.2510, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -242263.9896, places=3)

    @unittest.skipUnless(run_all_tests, "Skipping OMM tests on large systems")
    def test_big_dispersion_correction(self):
        """ Compare OpenMM and CHAMBER w/out vdW corr on big system """
        parm = ChamberParm(get_fn('dhfr_cmap_pbc.parm7'), get_fn('dhfr_cmap_pbc.rst7'))
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
        energies = energy_decomposition(parm, sim.context)
#NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
#Etot   =   -226511.4095  EKtot   =         0.0000  EPtot      =   -226511.4095
#BOND   =      8582.2932  ANGLE   =      5019.3573  DIHED      =       740.9529
#UB     =        29.6502  IMP     =        14.2581  CMAP       =      -216.2510
#1-4 NB =       345.6860  1-4 EEL =      6475.4350  VDWAALS    =     29959.9772
#EELEC  =   -277462.7684  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#Ewald error estimate:   0.3355E-03
        self.assertAlmostEqual(energies['bond'], 8582.2932, delta=5e-3)
        self.assertAlmostEqual(energies['angle'], 5019.3573, delta=5e-3)
        self.assertAlmostEqual(energies['urey_bradley'], 29.6502, delta=5e-4)
        self.assertAlmostEqual(energies['dihedral'], 740.9529, delta=5e-3)
        self.assertAlmostEqual(energies['improper'], 14.2581, delta=5e-4)
        self.assertRelativeEqual(energies['cmap'], -216.2510, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -240681.6702, places=4)

    @unittest.skipUnless(run_all_tests, "Skipping OMM tests on large systems")
    def test_big_shake(self):
        """ Compare OpenMM and CHAMBER PME excluding SHAKEn bonds (big) """
        parm = ChamberParm(get_fn('dhfr_cmap_pbc.parm7'), get_fn('dhfr_cmap_pbc.rst7'))
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

    def test_interface_pbc(self):
        """ Testing all ChamberParm.createSystem options (periodic) """
        parm = ChamberParm(get_fn('ala3_solv.parm7'), get_fn('ala3_solv.rst7'))
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

    def test_interface_customforce(self):
        """ Tests AmberParm.createSystem options with CustomNonbondedForce """
        parm = AmberParm(get_fn('solv2.parm7'), get_fn('solv2.rst7'))
        parm.parm_data['LENNARD_JONES_BCOEF'][0] = 0
        self.assertTrue(parm.has_NBFIX())
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=1.2*u.nanometers)
        # Get nonbonded forces
        forces = []
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                forces.append(force)
                self.assertEqual(force.getCutoffDistance(), 1.2*u.nanometers)
                self.assertEqual(force.getNonbondedMethod(), force.PME)
            elif isinstance(force, mm.CustomNonbondedForce):
                forces.append(force)
                self.assertEqual(force.getCutoffDistance(), 1.2*u.nanometers)
                self.assertEqual(force.getNonbondedMethod(), force.CutoffPeriodic)
        self.assertEqual(len(forces), 2)

    def test_interface_no_pbc(self):
        """Testing all ChamberParm.createSystem options (non-periodic)"""
        parm = ChamberParm(get_fn('ala_ala_ala.parm7'), get_fn('ala_ala_ala.rst7'))
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
                                   implicitSolventKappa=0.8*u.nanometers**-1)
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
        for nbmethod in (app.PME, app.CutoffPeriodic, app.Ewald):
            with self.assertRaises(ValueError):
                parm.createSystem(nonbondedMethod=nbmethod)

@unittest.skipUnless(has_openmm, 'Cannot test without OpenMM')
class TestAmoebaParm(TestCaseRelative):
    """ Tests some of the OMM integration with the AmoebaParm classes """

    def test_amoeba_forces(self):
        """ Test creation of some AMOEBA forces """
        parm = load_file(get_fn('amoeba.parm7'))
        self.assertIsInstance(parm.omm_bond_force(rigidWater=False), mm.AmoebaBondForce)
        self.assertIsInstance(parm.omm_angle_force(), mm.AmoebaAngleForce)
        self.assertIsInstance(parm.omm_trigonal_angle_force(), mm.AmoebaInPlaneAngleForce)
        self.assertIsInstance(parm.omm_out_of_plane_bend_force(), mm.AmoebaOutOfPlaneBendForce)
        self.assertIsInstance(parm.omm_pi_torsion_force(), mm.AmoebaPiTorsionForce)
        self.assertIsInstance(parm.omm_stretch_bend_force(), mm.AmoebaStretchBendForce)

        # Now some error handling
        # Trigonal angles
        old_deg = parm.trigonal_angle_types.degree
        del parm.trigonal_angle_types.degree
        self.assertRaises(ParameterError, parm.omm_trigonal_angle_force)
        parm.trigonal_angle_types.degree = old_deg
        parm.trigonal_angles[0].type = None
        self.assertRaises(ParameterError, parm.omm_trigonal_angle_force)
        # Out-of-plane bends
        old_deg = parm.out_of_plane_bend_types.degree
        del parm.out_of_plane_bend_types.degree
        self.assertRaises(ParameterError, parm.omm_out_of_plane_bend_force)
        parm.out_of_plane_bend_types.degree = old_deg
        parm.out_of_plane_bends[0].type = None
        self.assertRaises(ParameterError, parm.omm_out_of_plane_bend_force)
        # Pi-torsion parameters
        parm.pi_torsions[0].type = None
        self.assertRaises(ParameterError, parm.omm_pi_torsion_force)
        # stretch-bend parameters
        parm.stretch_bends[0].type = None
        self.assertRaises(ParameterError, parm.omm_stretch_bend_force)
