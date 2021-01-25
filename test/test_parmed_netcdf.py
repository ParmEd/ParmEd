"""
Tests the NetCDF file parsing capabilities with the different backends
"""

from random import randint
import numpy as np
from parmed import __version__
from parmed.amber.netcdffiles import NetCDFTraj, NetCDFRestart
from parmed import unit as u
from parmed.utils.six.moves import range, zip
from parmed.utils import PYPY
from parmed.utils.netcdf import NetCDFFile
from unittest import skipIf
from utils import get_fn, FileIOTestCase, TestCaseRelative

@skipIf(PYPY, 'NetCDF parsing does not yet work with pypy')
class TestNetCDF(FileIOTestCase, TestCaseRelative):
    """ Test NetCDF Functionality """

    def testNetCDF(self):
        """ Test scipy NetCDF parsing """
        traj = NetCDFTraj.open_old(get_fn('tz2.truncoct.nc'))
        self._check_traj(traj)
        rst = NetCDFRestart.open_old(get_fn('ncinpcrd.rst7'))
        self._check_rst(rst)
        # Now try writing files
        ntraj = NetCDFTraj.open_new(self.get_fn('test.nc', written=True), traj.atom, box=True)
        for crd in traj.coordinates:
            ntraj.add_coordinates(crd)
        for box in traj.box:
            ntraj.add_box(box)
        ntraj.close()
        traj = NetCDFTraj.open_old(self.get_fn('test.nc', written=True))
        self._check_traj(traj, written=True)

        nrst = NetCDFRestart.open_new(
            self.get_fn('test.ncrst', written=True), rst.atom, rst.hasbox, rst.hasvels, title=rst.title
        )
        nrst.coordinates = rst.coordinates
        nrst.velocities = rst.velocities
        nrst.time = rst.time
        nrst.box = rst.box
        nrst.close()
        rst = NetCDFRestart.open_old(self.get_fn('test.ncrst', written=True))
        self._check_rst(rst, written=True)
        # Now write NetCDF trajectory without any optional attributes/variables
        # and check proper reading
        natom = randint(100, 1000)
        nframe = randint(5, 20)
        coords = np.random.rand(nframe, natom, 3) * 20 - 10
        coords = np.array(coords, dtype='f')
        nctraj = NetCDFFile(self.get_fn('test2.nc', written=True), 'w', mmap=False)
        nctraj.Conventions = 'AMBER'
        nctraj.ConventionVersion = '1.0'
        nctraj.program = 'ParmEd'
        nctraj.programVersion = __version__
        nctraj.createDimension('frame', None)
        nctraj.createDimension('spatial', 3)
        nctraj.createDimension('atom', natom)
        v = nctraj.createVariable('spatial', 'c', ('spatial',))
        v[:] = np.asarray(list('xyz'))
        v = nctraj.createVariable('coordinates', 'f', ('frame', 'atom', 'spatial'))
        v[:] = coords
        del v
        nctraj.close()
        traj2 = NetCDFTraj.open_old(self.get_fn('test2.nc', written=True))
        np.testing.assert_allclose(traj2.coordinates, coords)

        # Now write NetCDF restart without any optional attributes/variables and
        # check proper reading
        natom = randint(100, 1000)
        nframe = randint(5, 20)
        coords = np.random.rand(natom, 3) * 20 - 10
        nctraj = NetCDFFile(self.get_fn('test2.ncrst', written=True), 'w', mmap=False)
        nctraj.Conventions = 'AMBERRESTART'
        nctraj.ConventionVersion = '1.0'
        nctraj.program = 'ParmEd'
        nctraj.programVersion = __version__
        nctraj.createDimension('frame', None)
        nctraj.createDimension('spatial', 3)
        nctraj.createDimension('atom', natom)
        v = nctraj.createVariable('spatial', 'c', ('spatial',))
        v[:] = np.asarray(list('xyz'))
        v = nctraj.createVariable('coordinates', 'd', ('atom', 'spatial'))
        v[:] = coords
        del v
        nctraj.close()
        traj2 = NetCDFRestart.open_old(self.get_fn('test2.ncrst', written=True))
        np.testing.assert_allclose(traj2.coordinates[0], coords)

    def testRemdFiles(self):
        """ Test proper reading and writing of NetCDF files with REMD info """
        rstfn = self.get_fn('restart.ncrst', written=True)
        trjfn = self.get_fn('traj.nc', written=True)
        # Do the restart with T-REMD first
        traj = NetCDFRestart.open_new(
            rstfn, 100, True, True, 'Restart w/ REMD', remd='Temperature', temp=300,
        )
        crd = np.random.rand(100, 3)
        vel = np.random.rand(100, 3)
        traj.coordinates = crd
        traj.velocities = vel
        traj.box = [40, 40, 40, 90, 90, 90]
        traj.close()
        traj = NetCDFRestart.open_old(rstfn)
        self.assertEqual(traj.temp0, 300)
        np.testing.assert_allclose(crd, traj.coordinates.squeeze(), atol=1e-5)
        np.testing.assert_allclose(vel, traj.velocities.squeeze(), atol=1e-5)

        traj = NetCDFRestart.open_new(
            rstfn, 100, False, True, 'Restart w/ REMD', remd='Multi', remd_dimtypes=[1, 3, 3, 3],
        )
        crd = np.random.rand(100, 3)
        vel = np.random.rand(100, 3)
        traj.coordinates = crd
        traj.velocities = vel
        remd_indices = np.random.choice(np.arange(20), size=4, replace=True)
        traj.remd_indices = remd_indices
        traj.close()
        traj = NetCDFRestart.open_old(rstfn)
        np.testing.assert_allclose(crd, traj.coordinates.squeeze(), atol=1e-5)
        np.testing.assert_allclose(vel, traj.velocities.squeeze(), atol=1e-5)
        np.testing.assert_equal(remd_indices, traj.remd_indices)
        np.testing.assert_equal([1, 3, 3, 3], traj.remd_dimtype)
        # Do the restart with T-REMD first
        traj = NetCDFTraj.open_new(
            trjfn, 100, True, True, True, title='Traj w/ REMD', remd='Temperature',
        )
        crd = np.random.rand(100, 3)
        vel = np.random.rand(100, 3)
        traj.add_coordinates(crd)
        traj.add_velocities(vel)
        traj.add_temp0(300)
        traj.add_box([40, 40, 40, 90, 90, 90])
        traj.close()
        traj = NetCDFTraj.open_old(trjfn)
        np.testing.assert_equal(traj.temp0, [300])
        np.testing.assert_allclose(crd, traj.coordinates.squeeze(), atol=1e-5)
        np.testing.assert_allclose(vel, traj.velocities.squeeze(), atol=1e-5)

        traj = NetCDFTraj.open_new(
            trjfn, 100, False, True, title='Traj w/ REMD', remd='Multi', remd_dimension=4, vels=True, frcs=True,
        )
        crd = np.random.rand(100, 3)
        vel = np.random.rand(100, 3)
        frc = np.random.rand(100, 3)
        traj.remd_dimtype = [1, 3, 3, 3]
        traj.add_coordinates(crd*u.angstroms)
        traj.add_velocities(vel*u.angstroms/u.picosecond)
        traj.add_forces(frc*u.kilocalories_per_mole/u.angstrom)
        remd_indices = np.random.choice(np.arange(20), size=4, replace=True)
        traj.add_remd_indices(remd_indices)
        traj.close()
        traj = NetCDFTraj.open_old(trjfn)
        np.testing.assert_allclose(crd, traj.coordinates.squeeze(), atol=1e-5)
        np.testing.assert_allclose(vel, traj.velocities.squeeze(), atol=1e-5)
        np.testing.assert_allclose(frc, traj.forces.squeeze(), atol=1e-5)
        np.testing.assert_equal(remd_indices, traj.remd_indices.squeeze())
        np.testing.assert_equal([1, 3, 3, 3], traj.remd_dimtype)

    def _check_traj(self, traj, written=False):
        """ Checks various trajectory properties """
        self.assertEqual(traj.Conventions, 'AMBER')
        if written:
            self.assertEqual(traj.application, 'AmberTools')
            self.assertEqual(traj.program, 'ParmEd')
            self.assertEqual(traj.programVersion, __version__)
        else:
            self.assertEqual(traj.application, 'AMBER')
            self.assertEqual(traj.program, 'sander')
            self.assertEqual(traj.programVersion, '9.0')
        self.assertEqual(traj.ConventionVersion, '1.0')
        runsum = traj.coordinates.sum()
        self.assertRelativeEqual(runsum, 7049.7598, places=4)

        for frame in range(traj.frame):
            box = traj.box[frame]
            self.assertAlmostEqual(box[0], box[1])
            self.assertAlmostEqual(box[0], box[2])
            self.assertAlmostEqual(box[3], 109.471219, places=6)
            self.assertAlmostEqual(box[4], 109.471219, places=6)
            self.assertAlmostEqual(box[5], 109.471219, places=6)

    def testBadNetCDFFiles(self):
        """ Tests error checking for bad usage of NetCDF files """
        with self.assertRaises(ValueError):
            NetCDFRestart.open_new(self.get_fn('test.ncrst', written=True),
                                   natom=10, box=True, vels=True, remd='Temperature')
        with self.assertRaises(ValueError):
            NetCDFRestart.open_new(self.get_fn('test.ncrst', written=True),
                                   natom=10, box=True, vels=True, remd='Multi')
        with self.assertRaises(ValueError):
                NetCDFRestart.open_new(self.get_fn('test.ncrst', written=True),
                                       natom=10, box=True, vels=True,
                                       remd='Multi', remd_dimtypes=[1, 3, 2])
        with self.assertRaises(ValueError):
            NetCDFRestart.open_new(self.get_fn('test.ncrst', written=True),
                                   natom=10, box=True, vels=True, remd='Illegal')
        with self.assertRaises(ValueError):
            NetCDFTraj.open_new(self.get_fn('test.nc', written=True),
                                natom=10, box=True, vels=False, remd='Multidim')
        with self.assertRaises(ValueError):
            NetCDFTraj.open_new(self.get_fn('test.nc', written=True),
                                natom=10, box=True, vels=False, remd='Illegal')

    def _check_rst(self, rst, written=False):
        """ Checks various restart properties """
        self.assertAlmostEqual(rst.time, 30.1)
        self.assertEqual(rst.Conventions, 'AMBERRESTART')
        self.assertEqual(rst.title, 'ACE')
        if written:
            self.assertEqual(rst.application, 'AmberTools')
            self.assertEqual(rst.program, 'ParmEd')
            self.assertEqual(rst.programVersion, __version__)
        else:
            self.assertEqual(rst.application, 'AMBER')
            self.assertEqual(rst.program, 'sander')
            self.assertEqual(rst.programVersion, '11.0')
        self.assertAlmostEqual(rst.cell_lengths[0], 30.2642725)
        self.assertAlmostEqual(rst.cell_lengths[1], 30.2642725)
        self.assertAlmostEqual(rst.cell_lengths[2], 30.2642725)
        self.assertAlmostEqual(rst.cell_angles[0], 109.471219)
        self.assertAlmostEqual(rst.cell_angles[1], 109.471219)
        self.assertAlmostEqual(rst.cell_angles[2], 109.471219)
        self.assertTrue(all([round(x-y, 7) == 0 for x, y in zip(rst.cell_lengths, rst.box[:3])]))
        self.assertTrue(all([round(x-y, 7) == 0 for x, y in zip(rst.cell_angles, rst.box[3:])]))
