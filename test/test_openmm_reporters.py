"""
This module tests the various reporters included in the parmed package
"""
import numpy as np
import os
from io import StringIO
from parmed import unit as u, load_file
from parmed.amber import AmberParm, AmberMdcrd, AmberAsciiRestart, NetCDFTraj, NetCDFRestart
from parmed.openmm.reporters import (
    NetCDFReporter, MdcrdReporter, ProgressReporter, RestartReporter, StateDataReporter,
    EnergyMinimizerReporter, _format_time
)
import sys
import unittest
from utils import mm, app, has_openmm, CPU, FileIOTestCase, HAS_GROMACS

@unittest.skipUnless(has_openmm, "Cannot test without OpenMM")
class TestStateDataReporter(FileIOTestCase):

    def setUp(self):
        super().setUp()
        self.amber_gas = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))

    def test_state_data_reporter(self):
        """ Test StateDataReporter with various options """
        system = self.amber_gas.createSystem()
        integrator = mm.LangevinIntegrator(300*u.kelvin, 5.0/u.picoseconds, 1.0*u.femtoseconds)
        sim = app.Simulation(self.amber_gas.topology, system, integrator, platform=CPU)
        sim.context.setPositions(self.amber_gas.positions)
        f = open(self.get_fn('akma5.dat', written=True), 'w')
        sim.reporters.extend([
            StateDataReporter(self.get_fn('akma1.dat', written=True), 10),
            StateDataReporter(self.get_fn('akma2.dat', written=True), 10,
                              time=False, potentialEnergy=False,
                              kineticEnergy=False, totalEnergy=False,
                              temperature=False),
            StateDataReporter(self.get_fn('akma3.dat', written=True), 10, volume=True, density=True),
            StateDataReporter(self.get_fn('akma4.dat', written=True), 10, separator='\t'),
            StateDataReporter(self.get_fn('units.dat', written=True), 10, volume=True, density=True,
                              energyUnit=u.kilojoules_per_mole, volumeUnit=u.nanometers**3),
            StateDataReporter(f, 10)
        ])
        sim.step(500)
        f.close()

        # Now open all of the reporters and check that the information in there
        # is what we expect it to be
        akma1 = open(self.get_fn('akma1.dat', written=True), 'r')
        akma2 = open(self.get_fn('akma2.dat', written=True), 'r')
        akma3 = open(self.get_fn('akma3.dat', written=True), 'r')
        akma4 = open(self.get_fn('akma4.dat', written=True), 'r')
        akma5 = open(self.get_fn('akma5.dat', written=True), 'r')
        units = open(self.get_fn('units.dat', written=True), 'r')
        # AKMA 1 first
        header = akma1.readline().strip()[1:].split(',')
        self.assertEqual(len(header), 6)
        for i, label in enumerate(
            ('Step', 'Time', 'Potential Energy', 'Kinetic Energy', 'Total Energy', 'Temperature')
        ):
            self.assertTrue(label in header[i])
        for i, line in enumerate(akma1):
            words = line.replace('\n', '').split(',')
            self.assertEqual(i*10+10, int(words[0])) # step
        akma1.close()
        # AKMA 2
        header = akma2.readline().strip()[1:].split(',')
        self.assertEqual(len(header), 1)
        self.assertTrue('Step' in header[0])
        for i, line in enumerate(akma2):
            self.assertEqual(int(line.replace('\n', '').split(',')[0]), 10*i+10)
        akma2.close()
        # AKMA 3 -- save energies so we can compare to the file with different
        # units
        header = akma3.readline().strip()[1:].split(',')
        self.assertEqual(len(header), 8)
        for i, label in enumerate(('Step', 'Time', 'Potential Energy',
                                   'Kinetic Energy', 'Total Energy',
                                   'Temperature', 'Box Volume', 'Density')):
            self.assertTrue(label in header[i])
        akma_energies = [[0.0 for i in range(8)] for j in range(50)]
        for i, line in enumerate(akma3):
            words = line.replace('\n', '').split(',')
            akma_energies[i][0] = int(words[0])
            for j in range(1, 8):
                akma_energies[i][j] = float(words[j])
        akma3.close()
        # AKMA 4 -- tab delimiter
        header = akma4.readline().strip()[1:].split('\t')
        self.assertEqual(len(header), 6)
        for i, label in enumerate(('Step', 'Time', 'Potential Energy',
                                   'Kinetic Energy', 'Total Energy',
                                   'Temperature')):
            self.assertTrue(label in header[i])
        for i, line in enumerate(akma4):
            words = line.replace('\n', '').split('\t')
            self.assertEqual(i*10+10, int(words[0])) # step
        akma4.close()
        # AKMA 5 -- write to open file handle
        header = akma5.readline().strip()[1:].split(',')
        self.assertEqual(len(header), 6)
        for i, label in enumerate(('Step', 'Time', 'Potential Energy',
                                   'Kinetic Energy', 'Total Energy',
                                   'Temperature')):
            self.assertTrue(label in header[i])
        for i, line in enumerate(akma5):
            words = line.replace('\n', '').split(',')
            self.assertEqual(i*10+10, int(words[0])) # step
        akma5.close()
        # UNITS -- compare other units
        ene = u.kilojoule_per_mole.conversion_factor_to(u.kilocalorie_per_mole)
        volume = u.nanometers.conversion_factor_to(u.angstroms)**3
        conversions = [1, 1, ene, ene, ene, 1, volume, 1]
        headers = units.readline().strip()[1:].split(',')
        self.assertEqual(len(headers), 8)
        for i, line in enumerate(units):
            words = line.replace('\n', '').split(',')
            self.assertEqual(int(words[0]), akma_energies[i][0])
            for j in range(1, 8):
                self.assertAlmostEqual(float(words[j])*conversions[j],
                                       akma_energies[i][j],
                                       places=5)
        units.close()

    def test_progress_reporter(self):
        """ Test ProgressReporter with various options """
        self.assertRaises(ValueError, lambda: ProgressReporter(sys.stdout, 1, 5))
        system = self.amber_gas.createSystem()
        integrator = mm.LangevinIntegrator(300*u.kelvin, 5.0/u.picoseconds,
                                           1.0*u.femtoseconds)
        sim = app.Simulation(self.amber_gas.topology, system, integrator, platform=CPU)
        sim.context.setPositions(self.amber_gas.positions)
        sim.reporters.append(
            ProgressReporter(self.get_fn('progress_reporter.dat', written=True), 10,
                             500, step=True, time=True, potentialEnergy=True,
                             kineticEnergy=True, totalEnergy=True,
                             temperature=True, volume=True, density=True,
                             systemMass=None)
        )
        sim.step(500)
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 1)
        text = open(self.get_fn('progress_reporter.dat', written=True), 'r').read()
        self.assertTrue('Estimated time to completion' in text)
        self.assertTrue('Total Energy' in text)
        self.assertTrue('Potential Energy' in text)
        self.assertTrue('Kinetic Energy' in text)
        self.assertTrue('Temperature' in text)

@unittest.skipUnless(has_openmm, "Cannot test without OpenMM")
class TestTrajRestartReporter(FileIOTestCase):

    def setUp(self):
        super().setUp()
        self.amber_gas = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))

    def test_reporters(self):
        """ Test NetCDF and ASCII restart and trajectory reporters (no PBC) """
        with self.assertRaises(ValueError):
            NetCDFReporter(self.get_fn('blah', written=True), 1, crds=False)
        with self.assertRaises(ValueError):
            MdcrdReporter(self.get_fn('blah', written=True), 1, crds=False)
        with self.assertRaises(ValueError):
            MdcrdReporter(self.get_fn('blah', written=True), 1, crds=True, vels=True)
        system = self.amber_gas.createSystem()
        integrator = mm.LangevinIntegrator(300*u.kelvin, 5.0/u.picoseconds, 1.0*u.femtoseconds)
        sim = app.Simulation(self.amber_gas.topology, system, integrator, platform=CPU)
        sim.context.setPositions(self.amber_gas.positions)
        sim.reporters.extend([
            NetCDFReporter(self.get_fn('traj1.nc', written=True), 10),
            NetCDFReporter(self.get_fn('traj2.nc', written=True), 10, vels=True),
            NetCDFReporter(self.get_fn('traj3.nc', written=True), 10, frcs=True),
            NetCDFReporter(self.get_fn('traj4.nc', written=True), 10, vels=True, frcs=True),
            NetCDFReporter(self.get_fn('traj5.nc', written=True), 10, crds=False, vels=True),
            NetCDFReporter(self.get_fn('traj6.nc', written=True), 10, crds=False, frcs=True),
            NetCDFReporter(self.get_fn('traj7.nc', written=True), 10, crds=False, vels=True, frcs=True),
            MdcrdReporter(self.get_fn('traj1.mdcrd', written=True), 10),
            MdcrdReporter(self.get_fn('traj2.mdcrd', written=True), 10, crds=False, vels=True),
            MdcrdReporter(self.get_fn('traj3.mdcrd', written=True), 10, crds=False, frcs=True),
            RestartReporter(self.get_fn('restart.ncrst', written=True), 10, write_multiple=True, netcdf=True),
            RestartReporter(self.get_fn('restart.rst7', written=True), 10),
        ])
        sim.step(500)
        for reporter in sim.reporters:
            reporter.finalize()

        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 61)
        ntraj = [NetCDFTraj.open_old(self.get_fn('traj1.nc', written=True)),
                 NetCDFTraj.open_old(self.get_fn('traj2.nc', written=True)),
                 NetCDFTraj.open_old(self.get_fn('traj3.nc', written=True)),
                 NetCDFTraj.open_old(self.get_fn('traj4.nc', written=True)),
                 NetCDFTraj.open_old(self.get_fn('traj5.nc', written=True)),
                 NetCDFTraj.open_old(self.get_fn('traj6.nc', written=True)),
                 NetCDFTraj.open_old(self.get_fn('traj7.nc', written=True))]
        atraj = [
            AmberMdcrd(self.get_fn('traj1.mdcrd', written=True), self.amber_gas.ptr('natom'), hasbox=False, mode='r'),
            AmberMdcrd(self.get_fn('traj2.mdcrd', written=True), self.amber_gas.ptr('natom'), hasbox=False, mode='r'),
            AmberMdcrd(self.get_fn('traj3.mdcrd', written=True), self.amber_gas.ptr('natom'), hasbox=False, mode='r'),
        ]
        for traj in ntraj:
            self.assertEqual(traj.frame, 50)
            self.assertEqual(traj.Conventions, 'AMBER')
            self.assertEqual(traj.ConventionVersion, '1.0')
            self.assertEqual(traj.application, 'AmberTools')
            self.assertEqual(traj.program, 'ParmEd')
            self.assertFalse(traj.hasbox)
        self.assertTrue(ntraj[0].hascrds)
        self.assertFalse(ntraj[0].hasvels)
        self.assertFalse(ntraj[0].hasfrcs)
        self.assertTrue(ntraj[1].hascrds)
        self.assertTrue(ntraj[1].hasvels)
        self.assertFalse(ntraj[1].hasfrcs)
        self.assertTrue(ntraj[2].hascrds)
        self.assertFalse(ntraj[2].hasvels)
        self.assertTrue(ntraj[2].hasfrcs)
        self.assertTrue(ntraj[3].hascrds)
        self.assertTrue(ntraj[3].hasvels)
        self.assertTrue(ntraj[3].hasfrcs)
        self.assertFalse(ntraj[4].hascrds)
        self.assertTrue(ntraj[4].hasvels)
        self.assertFalse(ntraj[4].hasfrcs)
        self.assertFalse(ntraj[5].hascrds)
        self.assertFalse(ntraj[5].hasvels)
        self.assertTrue(ntraj[5].hasfrcs)
        self.assertFalse(ntraj[6].hascrds)
        self.assertTrue(ntraj[6].hasvels)
        self.assertTrue(ntraj[6].hasfrcs)
        for i in (0, 2, 3, 4, 5, 6):
            ntraj[i].close() # still need the 2nd
        for traj in atraj:
            traj.close()
        # Now test the NetCDF restart files
        fn = self.get_fn('restart.ncrst.%d', written=True)
        for i, j in enumerate(range(10, 501, 10)):
            ncrst = NetCDFRestart.open_old(fn % j)
            self.assertEqual(ncrst.coordinates.shape, (1, 25, 3))
            self.assertEqual(ncrst.velocities.shape, (1, 25, 3))
            np.testing.assert_allclose(ncrst.coordinates[0], ntraj[1].coordinates[i])
            np.testing.assert_allclose(ncrst.velocities[0], ntraj[1].velocities[i], rtol=1e-6)
        # Now test the ASCII restart file
        f = AmberAsciiRestart(self.get_fn('restart.rst7', written=True), 'r')
        # Compare to ncrst and make sure it's the same data
        np.testing.assert_allclose(ncrst.coordinates, f.coordinates, atol=1e-3)
        np.testing.assert_allclose(ncrst.velocities, f.velocities, rtol=1e-3)
        # Make sure the EnergyMinimizerReporter does not fail
        f = StringIO()
        rep = EnergyMinimizerReporter(f)
        rep.report(sim, frame=10)
        rep.finalize()

    @unittest.skipUnless(HAS_GROMACS, 'Cannot test without GROMACS')
    def test_reporters_pbc(self):
        """ Test NetCDF and ASCII restart and trajectory reporters (w/ PBC) """
        systemsolv = load_file(self.get_fn('ildn.solv.top'), xyz=self.get_fn('ildn.solv.gro'))
        system = systemsolv.createSystem(nonbondedMethod=app.PME,
                                         nonbondedCutoff=8*u.angstroms)
        integrator = mm.LangevinIntegrator(300*u.kelvin, 5.0/u.picoseconds,
                                           1.0*u.femtoseconds)
        sim = app.Simulation(systemsolv.topology, system, integrator, CPU)
        sim.context.setPositions(systemsolv.positions)
        sim.reporters.extend([
                NetCDFReporter(self.get_fn('traj.nc', written=True), 1, vels=True, frcs=True),
                MdcrdReporter(self.get_fn('traj.mdcrd', written=True), 1),
                RestartReporter(self.get_fn('restart.ncrst', written=True), 1, netcdf=True),
                RestartReporter(self.get_fn('restart.rst7', written=True), 1),
                StateDataReporter(self.get_fn('state.o', written=True), 1, volume=True, density=True, systemMass=1)
        ])
        sim.step(5)
        for reporter in sim.reporters: reporter.finalize()
        ntraj = NetCDFTraj.open_old(self.get_fn('traj.nc', written=True))
        atraj = AmberMdcrd(self.get_fn('traj.mdcrd', written=True), len(systemsolv.atoms), True, mode='r')
        nrst = NetCDFRestart.open_old(self.get_fn('restart.ncrst', written=True))
        arst = AmberAsciiRestart(self.get_fn('restart.rst7', written=True), 'r')
        self.assertEqual(ntraj.frame, 5)
        self.assertEqual(atraj.frame, 5)
        self.assertTrue(ntraj.hasvels)
        self.assertTrue(ntraj.hasfrcs)
        for i in range(ntraj.frame):
            for x1, x2 in zip(ntraj.box[i], atraj.box[i]):
                self.assertAlmostEqual(x1, x2, places=3)
        self.assertEqual(len(nrst.box), 6)
        self.assertEqual(len(arst.box), 6)
        # Make sure the EnergyMinimizerReporter does not fail
        f = StringIO()
        rep = EnergyMinimizerReporter(f, volume=True)
        rep.report(sim)
        rep.finalize()

    def test_private_functions(self):
        """ Tests private helper functions for OMM reporters """
        self.assertEqual(_format_time(7200), (2, 'hr.'))
        self.assertEqual(_format_time(3600), (60, 'min.'))
        self.assertEqual(_format_time(40), (40, 'sec.'))
