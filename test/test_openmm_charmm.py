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
from __future__ import division

try:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    import simtk.unit as u
    from chemistry.charmm.openmmloader import (
                OpenMMCharmmPsfFile as CharmmPsfFile,
                OpenMMCharmmCrdFile as CharmmCrdFile,
                OpenMMCharmmRstFile as CharmmRstFile,
    )
    from chemistry.amber.openmmloader import OpenMMRst7 as Rst7
    PDBFile = app.PDBFile
    has_openmm = True
except ImportError:
    from chemistry.charmm.psf import CharmmPsfFile
    from chemistry.charmm.charmmcrds import CharmmCrdFile, CharmmRstFile
    from chemistry.amber.readparm import Rst7
    # To prevent NameError's
    def PDBFile(*args, **kwargs): return None
    has_openmm = False

from chemistry.charmm.parameters import CharmmParameterSet
from copy import copy
from math import sqrt
import unittest
import utils
    
get_fn = utils.get_fn

if has_openmm:
    # System
    charmm_gas = CharmmPsfFile(get_fn('ala_ala_ala.psf'))
    charmm_gas_crds = PDBFile(get_fn('ala_ala_ala.pdb'))
    charmm_solv = CharmmPsfFile(get_fn('dhfr_cmap_pbc.psf'))
    charmm_solv_crds = CharmmCrdFile(get_fn('dhfr_min_charmm.crd'))

    # Parameter sets
    param22 = CharmmParameterSet(get_fn('top_all22_prot.inp'),
                                 get_fn('par_all22_prot.inp'))

    # Make sure all precisions are double
    for i in range(mm.Platform.getNumPlatforms()):
        plat = mm.Platform.getPlatform(i)
        if plat.getName() == 'CUDA':
            plat.setPropertyDefaultValue('CudaPrecision', 'double')
        if plat.getName() == 'OpenCL':
            plat.setPropertyDefaultValue('OpenCLPrecision', 'double')

class TestCharmmFiles(unittest.TestCase):

    def setUp(self):
        if not hasattr(charmm_solv, 'box') or charmm_solv.box is None:
            charmm_solv.setBox(95.387*u.angstrom, 80.381*u.angstrom,
                               80.225*u.angstrom)

    def assertRelativeEqual(self, val1, val2, places=7):
        if val1 == val2: return
        try:
            ratio = val1 / val2
        except ZeroDivisionError:
            return self.assertAlmostEqual(val1, val2, places=places)
        else:
            if abs(round(ratio - 1, places)) == 0:
                return
            raise self.failureException(
                        '%s != %s with relative tolerance %g (%f)' %
                        (val1, val2, 10**-places, ratio)
            )
            return self.assertAlmostEqual(ratio, 1.0, places=places)
    
    def testGasEnergy(self):
        """ Compare OpenMM and CHARMM gas phase energies """
        parm = charmm_gas
        system = parm.createSystem(param22)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator)
        sim.context.setPositions(charmm_gas_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertAlmostEqual(energies['bond'], 1.3351, places=4)
        self.assertAlmostEqual(energies['angle'], 14.1158, places=4)
        self.assertAlmostEqual(energies['urey'], 0.3669, places=4)
        self.assertAlmostEqual(energies['dihedral'], 14.2773, places=4)
        self.assertAlmostEqual(energies['improper'], 0.3344, places=4)
        self.assertAlmostEqual(energies['cmap'], -0.5239, places=4)
        self.assertRelativeEqual(energies['nonbond'], 9.2210, places=4)

    def testGB1Energy(self): # HCT (uses mbondi radii internally)
        """ Compare OpenMM and CHARMM GB (igb=1) energies """
        parm = charmm_gas
        system = parm.createSystem(param22, implicitSolvent=app.HCT)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator)
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
        sim = app.Simulation(parm.topology, system, integrator)
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
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator)
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
        sim = app.Simulation(parm.topology, system, integrator)
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
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator)
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
        sim = app.Simulation(parm.topology, system, integrator)
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
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator)
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
        sim = app.Simulation(parm.topology, system, integrator)
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
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator)
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
        sim = app.Simulation(parm.topology, system, integrator)
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
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator)
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
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setUseDispersionCorrection(False)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator)
        sim.context.setPositions(charmm_solv_crds.positions)
        energies = decomposed_energy(sim.context, parm)
        self.assertRelativeEqual(energies['bond'], 8578.9872739, places=5)
        self.assertRelativeEqual(energies['angle'], 5018.3206306, places=5)
        self.assertRelativeEqual(energies['urey'], 29.6489539, places=5)
        self.assertRelativeEqual(energies['dihedral'], 740.9486106, places=5)
        self.assertRelativeEqual(energies['improper'], 14.2418054, places=5)
        self.assertRelativeEqual(energies['cmap'], -216.1422183, places=5)
        self.assertRelativeEqual(energies['nonbond'], -240681.958902, places=5)

if has_openmm:
    def decomposed_energy(context, parm, NRG_UNIT=u.kilocalories_per_mole):
        """ Gets a decomposed energy for a given system """
        energies = {}
        # Get energy components
        s = context.getState(getEnergy=True,
                             enforcePeriodicBox=parm.boxVectors is not None,
                             groups=2**parm.BOND_FORCE_GROUP)
        energies['bond'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
        s = context.getState(getEnergy=True,
                             enforcePeriodicBox=parm.boxVectors is not None,
                             groups=2**parm.ANGLE_FORCE_GROUP)
        energies['angle'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
        s = context.getState(getEnergy=True,
                             enforcePeriodicBox=parm.boxVectors is not None,
                             groups=2**parm.DIHEDRAL_FORCE_GROUP)
        energies['dihedral'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
        s = context.getState(getEnergy=True,
                             enforcePeriodicBox=parm.boxVectors is not None,
                             groups=2**parm.NONBONDED_FORCE_GROUP)
        energies['nonbond'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
        s = context.getState(getEnergy=True,
                             enforcePeriodicBox=parm.boxVectors is not None,
                             groups=2**parm.UREY_BRADLEY_FORCE_GROUP)
        energies['urey'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
        s = context.getState(getEnergy=True,
                             enforcePeriodicBox=parm.boxVectors is not None,
                             groups=2**parm.IMPROPER_FORCE_GROUP)
        energies['improper'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
        s = context.getState(getEnergy=True,
                             enforcePeriodicBox=parm.boxVectors is not None,
                             groups=2**parm.CMAP_FORCE_GROUP)
        energies['cmap'] = s.getPotentialEnergy().value_in_unit(NRG_UNIT)
        return energies

if not has_openmm:
    del TestCharmmFiles
