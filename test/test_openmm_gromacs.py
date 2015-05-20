"""
Contains unittests for running OpenMM calculations using the Amber file parsers
"""
from __future__ import division, print_function, absolute_import

try:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    has_openmm = True
except ImportError:
    from chemistry.amber.readparm import AmberParm, ChamberParm, Rst7
    has_openmm = False

from chemistry import load_file, ExtraPoint
from chemistry.gromacs import GromacsTopologyFile, GromacsGroFile
from chemistry.openmm.utils import energy_decomposition
import chemistry.unit as u
from chemistry.utils.six.moves import range, zip
from chemistry.vec3 import Vec3
import os
import unittest
import utils

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

@unittest.skipIf(not has_openmm(), "Cannot test without OpenMM")
class TestGromacsTop(utils.TestCaseRelative):
    """ Test ParmEd's energies vs. Gromacs energies as run by Lee-Ping """

    def testTiny(self):
        """ Test tiny Gromacs system nrg and frc (no PBC) """
        # Load the top and gro files
        top = load_file(os.path.join(get_fn('01.1water'), 'topol.top'))
        gro = load_file(os.path.join(get_fn('01.1water'), 'conf.gro'))

        # create the system and context, then calculate the energy decomposition
        system = top.createSystem(constraints=app.HBonds, rigidWater=True)
        context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
        context.setPositions(gro.positions)
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
        top = load_file(os.path.join(get_fn('02.6water'), 'topol.top'))
        gro = load_file(os.path.join(get_fn('02.6water'), 'conf.gro'))

        # create the system and context, then calculate the energy decomposition
        system = top.createSystem(constraints=app.HBonds, rigidWater=True,
                                  flexibleConstraints=True)
        context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
        context.setPositions(gro.positions)
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
        top = load_file(os.path.join(get_fn('04.Ala'), 'topol.top'))
        gro = load_file(os.path.join(get_fn('04.Ala'), 'conf.gro'))
        # create the system and context, then calculate the energy decomposition
        system = top.createSystem()
        context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
        context.setPositions(gro.positions)
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
        top = load_file(os.path.join(get_fn('03.AlaGlu'), 'topol.top'))
        gro = load_file(os.path.join(get_fn('03.AlaGlu'), 'conf.gro'))
        # create the system and context, then calculate the energy decomposition
        system = top.createSystem()
        context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
        context.setPositions(gro.positions)
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
        top = load_file(os.path.join(get_fn('07.DHFR-Liquid-NoPBC'), 'topol.top'))
        gro = load_file(os.path.join(get_fn('07.DHFR-Liquid-NoPBC'), 'conf.gro'))

        # Create the system and context, then calculate the energy decomposition
        system = top.createSystem()
        context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
        context.setPositions(gro.positions)
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
        top = load_file(os.path.join(get_fn('09.DHFR-PME'), 'topol.top'))
        gro = load_file(os.path.join(get_fn('09.DHFR-PME'), 'conf.gro'))
        top.box = gro.box[:]

        # Create the system and context, then calculate the energy decomposition
        system = top.createSystem(nonbondedMethod=app.PME,
                                  constraints=app.HBonds,
                                  nonbondedCutoff=0.9*u.nanometers,
                                  ewaldErrorTolerance=1.0e-5)
        context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
        context.setPositions(gro.positions)
        energies = energy_decomposition(top, context, nrg=u.kilojoules_per_mole)

        # Compare with Lee-Ping's answers. Make sure we zero-out forces for
        # virtual sites, since OMM doesn't do this and Gromacs does.
        self.assertAlmostEqual(energies['bond'], 1628.54739, places=3)
        self.assertAlmostEqual(energies['angle'], 3604.58751, places=3)
        self.assertAlmostEqual(energies['dihedral'], 6490.00844, delta=0.002)
        self.assertRelativeEqual(energies['nonbonded'], 23616.457584, places=4)
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
        top = load_file(os.path.join(get_fn('10.DHFR-PME-Switch'), 'topol.top'))
        gro = load_file(os.path.join(get_fn('10.DHFR-PME-Switch'), 'conf.gro'))
        top.box = gro.box[:]

        # Create the system and context, then calculate the energy decomposition
        system = top.createSystem(nonbondedMethod=app.PME,
                                  constraints=app.HBonds,
                                  nonbondedCutoff=0.9*u.nanometers,
                                  ewaldErrorTolerance=1.0e-5)
        context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
        context.setPositions(gro.positions)
        energies = energy_decomposition(top, context, nrg=u.kilojoules_per_mole)

        # Compare with Lee-Ping's answers. Make sure we zero-out forces for
        # virtual sites, since OMM doesn't do this and Gromacs does.
        self.assertAlmostEqual(energies['bond'], 1628.54739, places=3)
        self.assertAlmostEqual(energies['angle'], 3604.58751, places=3)
        self.assertAlmostEqual(energies['dihedral'], 6490.00844, delta=0.002)
        self.assertRelativeEqual(energies['nonbonded'], 23616.457584, places=4)
        gmxfrc = get_forces_from_xvg(
                os.path.join(get_fn('10.DHFR-PME-Switch'), 'force.xvg'))
        ommfrc = context.getState(getForces=True).getForces().value_in_unit(
                    u.kilojoules_per_mole/u.nanometer)
        zero_ep_frc(ommfrc, top)
        max_diff = get_max_diff(gmxfrc, ommfrc)
        self.assertLess(max_diff, 5)
