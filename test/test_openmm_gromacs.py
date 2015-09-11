"""
Contains unittests for running OpenMM calculations using the Amber file parsers
"""
from __future__ import division, print_function, absolute_import
import utils

try:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    has_openmm = True
    CPU = mm.Platform.getPlatformByName('CPU')
except ImportError:
    from parmed.amber.readparm import AmberParm, ChamberParm, Rst7
    has_openmm = False

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

get_fn = utils.get_fn

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

@unittest.skipIf(not has_openmm, "Cannot test without OpenMM")
class TestGromacsTop(utils.TestCaseRelative):
    """ Test ParmEd's energies vs. Gromacs energies as run by Lee-Ping """

    def setUp(self):
        warnings.filterwarnings('always', category=GromacsWarning)

    def testTiny(self):
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

    def testVerySmall(self):
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

    def testSmallPeptide(self):
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

    def testSmallDoublePeptide(self):
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

    @unittest.skipIf(utils.skip_big_tests(), "Skipping long-running tests")
    def testJAC(self):
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

    @unittest.skipIf(utils.skip_big_tests(), "Skipping long-running tests")
    def testJACPME(self):
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

    def testJACPMESwitch(self):
        """ Tests the DHFR Gromacs system nrg and force (PME w/ switch) """
        # Load the top and gro files
        top = load_file(os.path.join(get_fn('10.DHFR-PME-Switch'), 'topol.top'),
                        xyz=os.path.join(get_fn('10.DHFR-PME-Switch'), 'conf.gro'))
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
        self.assertRelativeEqual(energies['nonbonded'], 23616.457584, places=3)
        gmxfrc = get_forces_from_xvg(
                os.path.join(get_fn('10.DHFR-PME-Switch'), 'force.xvg'))
        ommfrc = context.getState(getForces=True).getForces().value_in_unit(
                    u.kilojoules_per_mole/u.nanometer)
        zero_ep_frc(ommfrc, top)
        max_diff = get_max_diff(gmxfrc, ommfrc)
        self.assertLess(max_diff, 5)

    def testDPPC(self):
        """ Tests non-standard Gromacs force fields and nonbonded exceptions """
        # We know what we're doing
        warnings.filterwarnings('ignore', category=GromacsWarning)
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

    def testOPLS(self):
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

    def testRoundTrip(self):
        """ Test ParmEd -> OpenMM round trip with Gromacs system """
        # Use DPPC to get RB-torsions tested. Also check that it initially fails
        # with the CustomNonbondedForce
        warnings.filterwarnings('ignore', category=GromacsWarning)
        top = load_file(os.path.join(get_fn('12.DPPC'), 'topol.top'),
                        xyz=os.path.join(get_fn('12.DPPC'), 'conf.gro'))
        self.assertEqual(top.combining_rule, 'lorentz')

        system = top.createSystem()
        def bad_system():
            return openmm.load_topology(top.topology, system).createSystem()
        warnings.filterwarnings('ignore', category=OpenMMWarning)
        self.assertTrue(
                openmm.load_topology(top.topology, system).unknown_functional
        )
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
