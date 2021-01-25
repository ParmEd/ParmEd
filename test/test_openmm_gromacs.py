"""
Contains unittests for running OpenMM calculations using the Amber file parsers
"""
from parmed import load_file, ExtraPoint, openmm, gromacs
from parmed.gromacs import GromacsTopologyFile, GromacsGroFile
from parmed.openmm.utils import energy_decomposition
from parmed.exceptions import GromacsWarning, ParmedError, OpenMMWarning
import parmed.unit as u
from parmed.utils.six.moves import range, zip
from parmed.vec3 import Vec3
import os
import unittest
import warnings
from utils import (
    get_fn, CPU, has_openmm, mm, app, TestCaseRelative, run_all_tests, QuantityTestCase, HAS_GROMACS
)

# OpenMM NonbondedForce methods are enumerated values. From NonbondedForce.h,
# they are:
#   0 - NoCutoff
#   1 - CutoffNonPeriodic
#   2 - CutoffPeriodic
#   3 - Ewald
#   4 - PME

def get_forces_from_xvg(frcfile):
    with open(frcfile, 'r') as f:
        frc = [float(x) for x in f.readline().split()[1:]]
        return [Vec3(frc[i], frc[i+1], frc[i+2]) for i in range(0, len(frc), 3)]

def get_max_diff(f1, f2):
    """ Gets the maximum difference between two force vectors """
    maxfrc = None
    assert len(f1) == len(f2), 'vector size mismatch'
    for v1, v2 in zip(f1, f2):
        norm = u.norm(v1 - v2)
        if maxfrc is None:
            maxfrc = norm
        else:
            maxfrc = max(maxfrc, norm)
    return maxfrc

def zero_ep_frc(frc, struct):
    vec0 = Vec3(0.0, 0.0, 0.0)
    for i, atom in enumerate(struct.atoms):
        if isinstance(atom, ExtraPoint):
            frc[i] = vec0

@unittest.skipUnless(has_openmm, "Cannot test without OpenMM")
@unittest.skipUnless(HAS_GROMACS, "Cannot test without GROMACS")
class TestGromacsTop(TestCaseRelative, QuantityTestCase):
    """ Test ParmEd's energies vs. Gromacs energies as run by Lee-Ping """

    def test_tiny(self):
        """ Test tiny Gromacs system nrg and frc (no PBC) """
        # Load the top and gro files
        top = load_file(os.path.join(get_fn('01.1water'), 'topol.top'),
                        xyz=os.path.join(get_fn('01.1water'), 'conf.gro'))
        self.assertEqual(top.combining_rule, 'lorentz')

        # create the system and context, then calculate the energy decomposition
        system = top.createSystem(constraints=app.HBonds, rigidWater=True)
        context = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        context.setPositions(top.positions)
        energies = energy_decomposition(top, context, nrg=u.kilojoules_per_mole)

        # Make sure they match Lee-Ping's answers (both energies and forces)
        self.assertAlmostEqual(energies['bond'], 1.63805, places=4)
        self.assertAlmostEqual(energies['angle'], 0.014803, places=4)
        gmxfrc = get_forces_from_xvg(os.path.join(get_fn('01.1water'), 'force.xvg'))
        ommfrc = context.getState(getForces=True).getForces().value_in_unit(
                    u.kilojoules_per_mole/u.nanometer)
        max_diff = get_max_diff(gmxfrc, ommfrc)
        self.assertLess(max_diff, 0.01)

    def test_very_small(self):
        """ Test very small Gromacs system nrg and frc (no PBC) """
        # Load the top and gro files
        top = load_file(os.path.join(get_fn('02.6water'), 'topol.top'),
                        xyz=os.path.join(get_fn('02.6water'), 'conf.gro'))
        self.assertEqual(top.combining_rule, 'lorentz')

        # create the system and context, then calculate the energy decomposition
        system = top.createSystem(constraints=app.HBonds, rigidWater=True,
                                  flexibleConstraints=True)
        context = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        context.setPositions(top.positions)
        energies = energy_decomposition(top, context, nrg=u.kilojoules_per_mole)

        # Compare with Lee-Ping's answers. Make sure we zero-out forces for
        # virtual sites, since OMM doesn't do this and Gromacs does.
        self.assertAlmostEqual(energies['bond'], 0, places=4)
        self.assertAlmostEqual(energies['nonbonded'], -108.277803, places=3)
        gmxfrc = get_forces_from_xvg(os.path.join(get_fn('02.6water'), 'force.xvg'))
        ommfrc = context.getState(getForces=True).getForces().value_in_unit(
                    u.kilojoules_per_mole/u.nanometer)
        zero_ep_frc(ommfrc, top)
        max_diff = get_max_diff(gmxfrc, ommfrc)
        self.assertLess(max_diff, 0.05)

    def test_small_peptide(self):
        """ Test alanine dipeptide Gromacs system nrg and frc (no PBC) """
        # Load the top and gro files
        top = load_file(os.path.join(get_fn('04.Ala'), 'topol.top'),
                        xyz=os.path.join(get_fn('04.Ala'), 'conf.gro'))
        self.assertEqual(top.combining_rule, 'lorentz')

        # create the system and context, then calculate the energy decomposition
        system = top.createSystem()
        context = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        context.setPositions(top.positions)
        energies = energy_decomposition(top, context, nrg=u.kilojoules_per_mole)

        # Compare with Lee-Ping's answers. Make sure we zero-out forces for
        # virtual sites, since OMM doesn't do this and Gromacs does.
        self.assertAlmostEqual(energies['bond'], 116.61637407, places=3)
        self.assertAlmostEqual(energies['angle'], 3.00686419356, places=4)
        self.assertAlmostEqual(energies['dihedral'], 46.2010926096, places=4)
        self.assertAlmostEqual(energies['nonbonded'], -103.924745659, places=3)
        gmxfrc = get_forces_from_xvg(os.path.join(get_fn('04.Ala'), 'force.xvg'))
        ommfrc = context.getState(getForces=True).getForces().value_in_unit(
                    u.kilojoules_per_mole/u.nanometer)
        zero_ep_frc(ommfrc, top)
        max_diff = get_max_diff(gmxfrc, ommfrc)
        self.assertLess(max_diff, 0.05)

    def test_small_double_peptide(self):
        """ Test interacting peptides Gromacs system nrg and frc (no PBC) """
        # Load the top and gro files
        top = load_file(os.path.join(get_fn('03.AlaGlu'), 'topol.top'),
                        xyz=os.path.join(get_fn('03.AlaGlu'), 'conf.gro'))
        self.assertEqual(top.combining_rule, 'lorentz')

        # create the system and context, then calculate the energy decomposition
        system = top.createSystem()
        context = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        context.setPositions(top.positions)
        energies = energy_decomposition(top, context, nrg=u.kilojoules_per_mole)

        # Compare with Lee-Ping's answers. Make sure we zero-out forces for
        # virtual sites, since OMM doesn't do this and Gromacs does.
        self.assertAlmostEqual(energies['bond'], 1.5307999, places=4)
        self.assertAlmostEqual(energies['angle'], 6.5804488, places=4)
        self.assertAlmostEqual(energies['dihedral'], 80.379714, places=4)
        self.assertAlmostEqual(energies['nonbonded'], -275.142487, places=3)
        gmxfrc = get_forces_from_xvg(os.path.join(get_fn('03.AlaGlu'), 'force.xvg'))
        ommfrc = context.getState(getForces=True).getForces().value_in_unit(
                    u.kilojoules_per_mole/u.nanometer)
        zero_ep_frc(ommfrc, top)
        max_diff = get_max_diff(gmxfrc, ommfrc)
        self.assertLess(max_diff, 0.05)

    @unittest.skipUnless(run_all_tests, "Skipping long-running tests")
    def test_jac(self):
        """ Tests the JAC benchmark Gromacs system nrg and force (no PBC) """
        # Load the top and gro files
        top = load_file(os.path.join(get_fn('07.DHFR-Liquid-NoPBC'), 'topol.top'),
                        xyz=os.path.join(get_fn('07.DHFR-Liquid-NoPBC'), 'conf.gro'))
        self.assertEqual(top.combining_rule, 'lorentz')

        # Create the system and context, then calculate the energy decomposition
        system = top.createSystem()
        context = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        context.setPositions(top.positions)
        energies = energy_decomposition(top, context, nrg=u.kilojoules_per_mole)

        # Compare with Lee-Ping's answers. Make sure we zero-out forces for
        # virtual sites, since OMM doesn't do this and Gromacs does.
        self.assertAlmostEqual(energies['bond'], 35.142565, places=3)
        self.assertAlmostEqual(energies['angle'], 3735.514669, places=3)
        self.assertAlmostEqual(energies['dihedral'], 7277.741635, delta=0.002)
        self.assertRelativeEqual(energies['nonbonded'], -288718.981405, places=4)
        gmxfrc = get_forces_from_xvg(
                os.path.join(get_fn('07.DHFR-Liquid-NoPBC'), 'force.xvg'))
        ommfrc = context.getState(getForces=True).getForces().value_in_unit(
                    u.kilojoules_per_mole/u.nanometer)
        zero_ep_frc(ommfrc, top)
        max_diff = get_max_diff(gmxfrc, ommfrc)
        self.assertLess(max_diff, 0.5)

    @unittest.skipUnless(run_all_tests, "Skipping long-running tests")
    def test_jac_pme(self):
        """ Tests the JAC benchmark Gromacs system nrg and force (PME) """
        # Load the top and gro files
        top = load_file(os.path.join(get_fn('09.DHFR-PME'), 'topol.top'),
                        xyz=os.path.join(get_fn('09.DHFR-PME'), 'conf.gro'))
        self.assertEqual(top.combining_rule, 'lorentz')

        # Create the system and context, then calculate the energy decomposition
        system = top.createSystem(nonbondedMethod=app.PME,
                                  constraints=app.HBonds,
                                  nonbondedCutoff=0.9*u.nanometers,
                                  ewaldErrorTolerance=1.0e-5)
        context = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        context.setPositions(top.positions)
        energies = energy_decomposition(top, context, nrg=u.kilojoules_per_mole)

        # Compare with Lee-Ping's answers. Make sure we zero-out forces for
        # virtual sites, since OMM doesn't do this and Gromacs does.
        self.assertAlmostEqual(energies['bond'], 1628.54739, places=3)
        self.assertAlmostEqual(energies['angle'], 3604.58751, places=3)
        self.assertAlmostEqual(energies['dihedral'], 6490.00844, delta=0.002)
        self.assertRelativeEqual(energies['nonbonded'], 23616.457584, delta=0.002)
        gmxfrc = get_forces_from_xvg(
                os.path.join(get_fn('09.DHFR-PME'), 'force.xvg'))
        ommfrc = context.getState(getForces=True).getForces().value_in_unit(
                    u.kilojoules_per_mole/u.nanometer)
        zero_ep_frc(ommfrc, top)
        max_diff = get_max_diff(gmxfrc, ommfrc)
        self.assertLess(max_diff, 5)

    def test_pme_switch(self):
        """ Tests the DHFR Gromacs system nrg and force (PME w/ switch) """
        # Load the top and gro files
        top = load_file(get_fn('ildn.solv.top'), xyz=get_fn('ildn.solv.gro'))
        self.assertEqual(top.combining_rule, 'lorentz')

        # Create the system and context, then calculate the energy decomposition
        system = top.createSystem(nonbondedMethod=app.PME,
                                  constraints=app.HBonds,
                                  flexibleConstraints=True,
                                  nonbondedCutoff=0.9*u.nanometers,
                                  ewaldErrorTolerance=1.0e-5,
                                  switchDistance=0.7*u.nanometers)
        context = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        context.setPositions(top.positions)
        energies = energy_decomposition(top, context, nrg=u.kilojoules_per_mole)

        # Compare with Lee-Ping's answers. Make sure we zero-out forces for
        # virtual sites, since OMM doesn't do this and Gromacs does.
        self.assertAlmostEqual(energies['bond'], 399.925189, places=4)
        self.assertAlmostEqual(energies['angle'], 36.18562, places=4)
        self.assertAlmostEqual(energies['dihedral'], 101.92265, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -18587.09715, places=4)

    def test_dppc(self):
        """ Tests non-standard Gromacs force fields and nonbonded exceptions """
        # We know what we're doing
        top = load_file(os.path.join(get_fn('12.DPPC'), 'topol.top'),
                        xyz=os.path.join(get_fn('12.DPPC'), 'conf.gro'))
        self.assertEqual(top.combining_rule, 'lorentz')

        # Create the system and context, then calculate the energy decomposition
        system = top.createSystem()
        context = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        context.setPositions(top.positions)
        energies = energy_decomposition(top, context, nrg=u.kilojoules_per_mole)

        # Compare with Lee-Ping's answers.
        self.assertAlmostEqual(energies['bond'], 0)
        self.assertAlmostEqual(energies['angle'], 1405.7354199, places=4)
        self.assertAlmostEqual(energies['dihedral'], 236.932663255, places=4)
        self.assertAlmostEqual(energies['improper'], 33.201541811, places=4)
        self.assertAlmostEqual(energies['rb_torsion'], 428.0550599, places=4)
        self.assertRelativeEqual(energies['nonbonded'], -16432.8092955, places=4)
        gmxfrc = get_forces_from_xvg(os.path.join(get_fn('12.DPPC'), 'force.xvg'))
        ommfrc = context.getState(getForces=True).getForces().value_in_unit(
                    u.kilojoules_per_mole/u.nanometer)
        zero_ep_frc(ommfrc, top)
        max_diff = get_max_diff(gmxfrc, ommfrc)
        self.assertLess(max_diff, 5)

    def test_opls(self):
        """ Tests geometric combining rules with OPLS/AA in Gromacs topology """
        top = load_file(os.path.join(get_fn('05.OPLS'), 'topol.top'),
                        xyz=os.path.join(get_fn('05.OPLS'), 'conf.gro'))
        self.assertEqual(top.combining_rule, 'geometric')

        # Create the system and context, then calculate the energy decomposition
        system = top.createSystem()
        context = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        context.setPositions(top.positions)
        energies = energy_decomposition(top, context, nrg=u.kilojoules_per_mole)

        # Compare with Gromacs energies
        self.assertAlmostEqual(energies['bond'], 332.178260935, places=4)
        self.assertAlmostEqual(energies['angle'], 29.2231806883, places=4)
        self.assertAlmostEqual(energies['dihedral'], 0.0057758656, places=4)
        self.assertAlmostEqual(energies['rb_torsion'], 55.603030624, places=4)
        self.assertRelativeEqual(energies['nonbonded'], 327.954397827, places=4)
        # Now make sure it can be set with PME
        system = top.createSystem(nonbondedMethod=app.PME,
                                  nonbondedCutoff=8*u.angstroms)
        n = 0
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                self.assertAlmostEqualQuantities(force.getCutoffDistance(),
                                                  8*u.angstroms)
                self.assertEqual(force.getNonbondedMethod(), force.PME)
                n += 1
            elif isinstance(force, mm.CustomNonbondedForce):
                self.assertAlmostEqualQuantities(force.getCutoffDistance(),
                                                  8*u.angstroms)
                self.assertEqual(force.getNonbondedMethod(), force.CutoffPeriodic)
                n += 1
        self.assertEqual(n, 2)
        # Now try with non-periodic cutoff
        system = top.createSystem(nonbondedMethod=app.CutoffNonPeriodic,
                                  nonbondedCutoff=16*u.angstroms,
                                  switchDistance=14*u.angstroms)
        n = 0
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce) or isinstance(force,
                    mm.CustomNonbondedForce):
                self.assertAlmostEqualQuantities(force.getCutoffDistance(),
                                                  16*u.angstroms)
                self.assertEqual(force.getNonbondedMethod(), force.CutoffNonPeriodic)
                self.assertTrue(force.getUseSwitchingFunction())
                self.assertAlmostEqualQuantities(force.getSwitchingDistance(),
                                                 14*u.angstroms)
                n += 1
        self.assertEqual(n, 2)

    def test_round_trip(self):
        """ Test ParmEd -> OpenMM round trip with Gromacs system """
        # Use DPPC to get RB-torsions tested. Also check that it initially fails
        # with the CustomNonbondedForce
        top = load_file(os.path.join(get_fn('12.DPPC'), 'topol.top'),
                        xyz=os.path.join(get_fn('12.DPPC'), 'conf.gro'))
        self.assertEqual(top.combining_rule, 'lorentz')

        system = top.createSystem()
        def bad_system():
            return openmm.load_topology(top.topology, system).createSystem()
        self.assertTrue(openmm.load_topology(top.topology, system).unknown_functional)
        self.assertRaises(ParmedError, bad_system)
        for i in range(len(system.getForces())):
            if isinstance(system.getForce(i), mm.CustomNonbondedForce):
                system.removeForce(i)
                break
        system2 = openmm.load_topology(top.topology, system).createSystem()
        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system2, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(top.positions)
        con2.setPositions(top.positions)
        ene1 = energy_decomposition(top, con1)
        ene2 = energy_decomposition(top, con2)
        self.assertAlmostEqual(ene1['bond'], ene2['bond'])
        self.assertAlmostEqual(ene1['angle'], ene2['angle'])
        self.assertAlmostEqual(ene1['dihedral'], ene2['dihedral'])
        self.assertAlmostEqual(ene1['improper'], ene2['improper'])
        self.assertAlmostEqual(ene1['rb_torsion'], ene2['rb_torsion'])
        self.assertAlmostEqual(ene1['nonbonded'], ene2['nonbonded'])
