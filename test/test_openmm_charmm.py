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
import utils

try:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    PDBFile = app.PDBFile
    has_openmm = True
except ImportError:
    # To prevent NameError's
    def PDBFile(*args, **kwargs): return None
    has_openmm = False

from parmed.amber.readparm import Rst7
from parmed.charmm import (CharmmPsfFile, CharmmCrdFile, CharmmRstFile,
                           CharmmParameterSet)
from parmed.exceptions import CharmmWarning
from parmed import unit as u, openmm
from parmed.utils.six.moves import range
from copy import copy
from math import sqrt
import unittest
import warnings
    
get_fn = utils.get_fn

# Suppress warning output from bad psf file... sigh.
warnings.filterwarnings('ignore', category=CharmmWarning)

if has_openmm:
    # System
    charmm_gas = CharmmPsfFile(get_fn('ala_ala_ala.psf'))
    charmm_gas_crds = PDBFile(get_fn('ala_ala_ala.pdb'))
    charmm_solv = CharmmPsfFile(get_fn('dhfr_cmap_pbc.psf'))
    charmm_solv_crds = CharmmCrdFile(get_fn('dhfr_min_charmm.crd'))
    charmm_nbfix = CharmmPsfFile(get_fn('ala3_solv.psf'))
    charmm_nbfix_crds = CharmmCrdFile(get_fn('ala3_solv.crd'))

    # Parameter sets
    param22 = CharmmParameterSet(get_fn('top_all22_prot.inp'),
                                 get_fn('par_all22_prot.inp'))
    param36 = CharmmParameterSet(get_fn('par_all36_prot.prm'),
                                 get_fn('toppar_water_ions.str'))

    CPU = mm.Platform.getPlatformByName('CPU')

def decomposed_energy(context, parm, NRG_UNIT=u.kilocalories_per_mole):
    """ Gets a decomposed energy for a given system """
    energies = {}
    # Get energy components
    s = context.getState(getEnergy=True, groups=1<<parm.BOND_FORCE_GROUP)
    energies['bond'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
    s = context.getState(getEnergy=True, groups=1<<parm.ANGLE_FORCE_GROUP)
    energies['angle'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
    s = context.getState(getEnergy=True, groups=1<<parm.DIHEDRAL_FORCE_GROUP)
    energies['dihedral'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
    s = context.getState(getEnergy=True, groups=1<<parm.NONBONDED_FORCE_GROUP)
    energies['nonbond'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
    s = context.getState(getEnergy=True,
                         groups=1<<parm.UREY_BRADLEY_FORCE_GROUP)
    energies['urey'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
    s = context.getState(getEnergy=True, groups=1<<parm.IMPROPER_FORCE_GROUP)
    energies['improper'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
    s = context.getState(getEnergy=True, groups=1<<parm.CMAP_FORCE_GROUP)
    energies['cmap'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
    return energies

@unittest.skipIf(not has_openmm, "Cannot test without OpenMM")
class TestCharmmFiles(utils.TestCaseRelative):

    def setUp(self):
        if charmm_solv.box is None:
            charmm_solv.box = [95.387, 80.381, 80.225, 90, 90, 90]
        if charmm_nbfix.box is None:
            charmm_nbfix.box = [3.271195e1, 3.299596e1, 3.300715e1, 90, 90, 90]

    def testGasEnergy(self):
        """ Compare OpenMM and CHARMM gas phase energies """
        parm = charmm_gas
        system = parm.createSystem(param22)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=4)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=4)
        self.assertAlmostEqual(energies['urey'], 0.3669, places=4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=4)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=4)
        self.assertRelativeEqual(energies['nonbond'], 9.2210, places=4)

    def testRoundTrip(self):
        """ Test ParmEd -> OpenMM round trip with CHARMM gas phase """
        parm = charmm_gas
        system = parm.createSystem(param22)
        self.assertEqual(parm.combining_rule, 'lorentz')
        system2 = openmm.load_topology(parm.topology, system).createSystem()
        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system2, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(charmm_gas_crds.positions)
        con2.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(con1, parm)
        energies2 = decomposed_energy(con2, parm)
        self.assertAlmostEqual(energies['bond'], energies2['bond'])
        self.assertAlmostEqual(energies['angle'], energies2['angle'])
        self.assertAlmostEqual(energies['urey'], energies2['urey'])
        self.assertAlmostEqual(energies['dihedral'], energies2['dihedral'])
        self.assertAlmostEqual(energies['improper'], energies2['improper'])
        self.assertAlmostEqual(energies['cmap'], energies2['cmap'])
        self.assertRelativeEqual(energies['nonbond'], energies2['nonbond'])

    def testGB1Energy(self): # HCT (uses mbondi radii internally)
        """ Compare OpenMM and CHARMM GB (igb=1) energies """
        parm = charmm_gas
        system = parm.createSystem(param22, implicitSolvent=app.HCT)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbond'], -102.1598379, places=5)
        system = parm.createSystem(param22, implicitSolvent=app.HCT,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbond'], -102.5012873, places=5)

    def testGB2Energy(self): # OBC1 (uses mbondi2 radii internally)
        """ Compare OpenMM and CHARMM GB (igb=2) energies """
        parm = charmm_gas
        system = parm.createSystem(param22, implicitSolvent=app.OBC1)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbond'], -107.8675, places=4)
        system = parm.createSystem(param22, implicitSolvent=app.OBC1,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbond'], -108.2129, places=4)

    def testGB5Energy(self): # OBC2 (uses mbondi2 radii internally)
        """ Compare OpenMM and CHARMM GB (igb=5) energies """
        parm = charmm_gas
        system = parm.createSystem(param22, implicitSolvent=app.OBC2)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbond'], -103.6186, places=4)
        system = parm.createSystem(param22, implicitSolvent=app.OBC2,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbond'], -103.9603, places=4)

    def testGB7Energy(self): # GBn (uses bondi radii internally)
        """ Compare OpenMM and CHARMM GB (igb=7) energies """
        parm = charmm_gas
        system = parm.createSystem(param22, implicitSolvent=app.GBn)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbond'], -109.4987850, places=5)
        system = parm.createSystem(param22, implicitSolvent=app.GBn,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbond'], -109.8465917, places=5)

    def testGB8Energy(self): # GBn2 (uses mbondi3 radii internally)
        """ Compare OpenMM and CHARMM GB (igb=8) energies """
        parm = charmm_gas
        system = parm.createSystem(param22, implicitSolvent=app.GBn2)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbond'], -108.1396, places=4)
        system = parm.createSystem(param22, implicitSolvent=app.GBn2,
                                   implicitSolventSaltConc=1.0*u.molar)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=3)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=3)
        self.assertAlmostEqual(energies['urey'], 0.3669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=3)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=3)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=3)
        self.assertRelativeEqual(energies['nonbond'], -108.4858, places=4)

    def testPME(self):
        """ Compare OpenMM and CHARMM PME energies """
        parm = charmm_solv
        system = parm.createSystem(param22, nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstrom)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_solv_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertRelativeEqual(energies['bond'], 8578.9872739, places=5)
        self.assertRelativeEqual(energies['angle'], 5018.3206306, places=5)
        self.assertRelativeEqual(energies['urey'], 29.6489539, places=5)
        self.assertRelativeEqual(energies['dihedral'], 740.9486106, places=5)
        self.assertRelativeEqual(energies['improper'], 14.2418054, places=5)
        self.assertRelativeEqual(energies['cmap'], -216.1422183, places=5)
        self.assertRelativeEqual(energies['nonbond'], -242262.368372, places=5)

    def testDispersionCorrection(self):
        """ Compare OpenMM and CHARMM PME energies w/out vdW correction """
        parm = charmm_solv
        system = parm.createSystem(param22, nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms,)
        self.assertEqual(parm.combining_rule, 'lorentz')
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setUseDispersionCorrection(False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_solv_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertRelativeEqual(energies['bond'], 8578.9872739, places=5)
        self.assertRelativeEqual(energies['angle'], 5018.3206306, places=5)
        self.assertRelativeEqual(energies['urey'], 29.6489539, places=5)
        self.assertRelativeEqual(energies['dihedral'], 740.9486106, places=5)
        self.assertRelativeEqual(energies['improper'], 14.2418054, places=5)
        self.assertRelativeEqual(energies['cmap'], -216.1422183, places=5)
        self.assertRelativeEqual(energies['nonbond'], -240681.958902, places=5)

    def testNBFIX(self):
        """ Test energies of systems with NBFIX modifications """
        parm = charmm_nbfix
        system = parm.createSystem(param36, nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms)
        self.assertEqual(parm.combining_rule, 'lorentz')
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(charmm_nbfix_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.1324212, places=4)
        self.assertAlmostEqual(energies['angle'], 1.06880188, places=4)
        self.assertAlmostEqual(energies['urey'], 0.06142407, places=4)
        self.assertAlmostEqual(energies['dihedral'], 7.81143025, places=4)
        self.assertAlmostEqual(energies['improper'], 0, places=4)
        self.assertAlmostEqual(energies['cmap'], 0.126790, places=4)
        self.assertRelativeEqual(energies['nonbond'], 6514.283116, places=4)

if __name__ == '__main__':
    unittest.main()
