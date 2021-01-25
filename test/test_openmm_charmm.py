"""
Contains unittests for running OpenMM calculations using native CHARMM files

The energies here should be the same as those defined in the Chamber tests of
test_openmm_amber.py, with a small caveat.  The energies computed in the chamber
tests all used the same (mbondi2) radii.  The radii for the GB models here is
set internally, so there is no way to reproduce the results from the chamber
tests without changing the tests (too time consuming) or changing the CHARMM
parser code in undesirable ways. So some of the nonbonded energies are close to
the analogous values in the Chamber tests, but have been recomputed using ParmEd
with the appropriate radius set.

The explicit solvent energies are also slightly different, since the coordinates
stored in the CHARMM coordinate file are of higher precision than those stored
in the Amber restart file. The energies are very close, but have been
recalculated using ParmEd for comparison here to improve the precision. The
implementation these energies were computed by has already been validated.
"""
from __future__ import division, print_function, absolute_import

from parmed.amber.readparm import Rst7
from parmed.charmm import CharmmPsfFile, CharmmCrdFile, CharmmRstFile, CharmmParameterSet
from parmed.exceptions import CharmmWarning, ParameterError
from parmed.openmm.utils import energy_decomposition
from parmed import unit as u, openmm, load_file, UreyBradley
from parmed.utils.six.moves import range
from copy import copy
from math import sqrt
import unittest
from utils import get_fn, TestCaseRelative, mm, app, has_openmm, CPU
import warnings

# System
charmm_gas = CharmmPsfFile(get_fn('ala_ala_ala.psf'))
charmm_gas_crds = load_file(get_fn('ala_ala_ala.pdb'))
charmm_nbfix = CharmmPsfFile(get_fn('ala3_solv.psf'))
charmm_nbfix_crds = CharmmCrdFile(get_fn('ala3_solv.crd'))
charmm_nbfix.box = [3.271195e1, 3.299596e1, 3.300715e1, 90, 90, 90]

# Parameter sets
param22 = CharmmParameterSet(get_fn('top_all22_prot.inp'), get_fn('par_all22_prot.inp'))
param36 = CharmmParameterSet(get_fn('par_all36_prot.prm'), get_fn('toppar_water_ions.str'))

@unittest.skipUnless(has_openmm, "Cannot test without OpenMM")
class TestCharmmFiles(TestCaseRelative):

    def test_gas_energy(self):
        """ Compare OpenMM and CHARMM gas phase energies """
        parm = charmm_gas
        system = parm.createSystem(param22)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=4)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, places=4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=4)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=4)
        self.assertRelativeEqual(energies['nonbonded'], 9.2210, places=4)

    def test_round_trip(self):
        """ Test ParmEd -> OpenMM round trip with CHARMM gas phase """
        parm = charmm_gas
        system = parm.createSystem(param22)
        self.assertEqual(parm.combining_rule, 'lorentz')
        system2 = openmm.load_topology(parm.topology, system).createSystem()
        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system2, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(charmm_gas_crds.positions)
        con2.setPositions(charmm_gas_crds.positions)
        energies = energy_decomposition(parm, con1)
        energies2 = energy_decomposition(parm, con2)
        self.assertAlmostEqual(energies['bond'], energies2['bond'])
        self.assertAlmostEqual(energies['angle'], energies2['angle'])
        self.assertAlmostEqual(energies['urey_bradley'], energies2['urey_bradley'])
        self.assertAlmostEqual(energies['dihedral'], energies2['dihedral'])
        self.assertAlmostEqual(energies['improper'], energies2['improper'])
        self.assertAlmostEqual(energies['cmap'], energies2['cmap'])
        self.assertRelativeEqual(energies['nonbonded'], energies2['nonbonded'])

    def test_gb1_energy(self): # HCT (uses mbondi radii internally)
        """ Compare OpenMM and CHARMM GB (igb=1) energies """
        parm = charmm_gas
        system = parm.createSystem(param22, implicitSolvent=app.HCT)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -102.1598379, places=5)
        system = parm.createSystem(param22, implicitSolvent=app.HCT,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -102.5012873, places=5)

    def test_gb2_energy(self): # OBC1 (uses mbondi2 radii internally)
        """ Compare OpenMM and CHARMM GB (igb=2) energies """
        parm = charmm_gas
        system = parm.createSystem(param22, implicitSolvent=app.OBC1)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -107.8675, places=4)
        system = parm.createSystem(param22, implicitSolvent=app.OBC1,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -108.2129, places=4)

    def test_gb5_energy(self): # OBC2 (uses mbondi2 radii internally)
        """ Compare OpenMM and CHARMM GB (igb=5) energies """
        parm = charmm_gas
        system = parm.createSystem(param22, implicitSolvent=app.OBC2)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -103.6186, places=4)
        system = parm.createSystem(param22, implicitSolvent=app.OBC2,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -103.9603, places=4)

    def test_gb7_energy(self): # GBn (uses bondi radii internally)
        """ Compare OpenMM and CHARMM GB (igb=7) energies """
        parm = charmm_gas
        system = parm.createSystem(param22, implicitSolvent=app.GBn)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -109.4987850, places=5)
        system = parm.createSystem(param22, implicitSolvent=app.GBn,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -109.8465917, places=5)

    def test_gb8_energy(self): # GBn2 (uses mbondi3 radii internally)
        """ Compare OpenMM and CHARMM GB (igb=8) energies """
        parm = charmm_gas
        system = parm.createSystem(param22, implicitSolvent=app.GBn2)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -108.1396, places=4)
        system = parm.createSystem(param22, implicitSolvent=app.GBn2,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey_bradley'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbonded'], -108.4858, places=4)

    def test_dispersion_correction(self):
        """ Compare OpenMM and CHARMM PME energies w/out vdW correction """
        parm = charmm_nbfix
        system = parm.createSystem(param36, nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,)
        self.assertEqual(parm.combining_rule, 'lorentz')
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setUseDispersionCorrection(False)
            elif isinstance(force, mm.CustomNonbondedForce):
                force.setUseLongRangeCorrection(False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_nbfix_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.1324212, places=4)
        self.assertAlmostEqual(energies['angle'], 1.06880188, places=4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.06142407, places=4)
        self.assertAlmostEqual(energies['dihedral'], 7.81143025, places=4)
        self.assertAlmostEqual(energies['improper'], 0, places=4)
        self.assertAlmostEqual(energies['cmap'], 0.126790, places=4)
        self.assertRelativeEqual(energies['nonbonded'], 6584.01821319, places=4)

    def test_nbfix(self):
        """ Test PME energies of systems with NBFIX modifications """
        parm = charmm_nbfix
        system = parm.createSystem(param36, nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_nbfix_crds.positions)
        energies = energy_decomposition(parm, sim.context)
        self.assertAlmostEqual(energies['bond'], 1.1324212, places=4)
        self.assertAlmostEqual(energies['angle'], 1.06880188, places=4)
        self.assertAlmostEqual(energies['urey_bradley'], 0.06142407, places=4)
        self.assertAlmostEqual(energies['dihedral'], 7.81143025, places=4)
        self.assertAlmostEqual(energies['improper'], 0, places=4)
        self.assertAlmostEqual(energies['cmap'], 0.126790, places=4)
        self.assertRelativeEqual(energies['nonbonded'], 6514.283116, places=4)

    def test_no_parameters(self):
        """ Test proper error handling when parameters not present """
        parm = CharmmPsfFile(get_fn('ala_ala_ala.psf'))
        parm.rb_torsions.append(parm.dihedrals.pop())
        parm.urey_bradleys.append(UreyBradley(parm.angles[0].atom1, parm.angles[0].atom3))
        self.assertRaises(ParameterError, parm.omm_bond_force)
        self.assertRaises(ParameterError, parm.omm_angle_force)
        self.assertRaises(ParameterError, parm.omm_dihedral_force)
        self.assertRaises(ParameterError, parm.omm_rb_torsion_force)
        self.assertRaises(ParameterError, parm.omm_urey_bradley_force)
        self.assertRaises(ParameterError, parm.omm_improper_force)
        self.assertRaises(ParameterError, parm.omm_cmap_force)

    def test_bad_system(self):
        """ Test error handling of impossible systems """
        parm = copy(charmm_gas)
        parm.load_parameters(param22)
        for d in parm.dihedral_types:
            try:
                for dt in d:
                    dt.scee = 0
            except TypeError:
                d.scee = 0
        self.assertRaises(ValueError, parm.createSystem)
