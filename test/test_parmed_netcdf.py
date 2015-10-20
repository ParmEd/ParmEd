"""
Tests the NetCDF file parsing capabilities with the different backends
"""

from random import randint
import numpy as np
from parmed import __version__
from parmed.utils.six.moves import range, zip
from parmed.utils import PYPY
from unittest import skipIf
from utils import get_fn, FileIOTestCase, TestCaseRelative

@skipIf(PYPY, 'NetCDF parsing does not yet work with pypy')
class TestNetCDF(FileIOTestCase, TestCaseRelative):
    """ Test NetCDF Functionality """
    
    def testNetCDF(self):
        """ Test scipy NetCDF parsing """
        from parmed.utils.netcdf import NetCDFFile
        from parmed.amber.netcdffiles import NetCDFTraj, NetCDFRestart
        traj = NetCDFTraj.open_old(get_fn('tz2.truncoct.nc'))
        self._check_traj(traj)
        rst = NetCDFRestart.open_old(get_fn('ncinpcrd.rst7'))
        self._check_rst(rst)
        # Now try writing files
        ntraj = NetCDFTraj.open_new(get_fn('test.nc', written=True),
                                    traj.atom, box=True)
        for crd in traj.coordinates:
            ntraj.add_coordinates(crd)
        for box in traj.box:
            ntraj.add_box(box)
        ntraj.close()
        traj = NetCDFTraj.open_old(get_fn('test.nc', written=True))
        self._check_traj(traj, written=True)

        nrst = NetCDFRestart.open_new(get_fn('test.ncrst', written=True),
                                      rst.atom, rst.hasbox, rst.hasvels,
                                      title=rst.title)
        nrst.coordinates = rst.coordinates
        nrst.velocities = rst.velocities
        nrst.time = rst.time
        nrst.box = rst.box
        nrst.close()
        rst = NetCDFRestart.open_old(get_fn('test.ncrst', written=True))
        self._check_rst(rst, written=True)
        # Now write NetCDF trajectory without any optional attributes/variables
        # and check proper reading
        natom = randint(100, 1000)
        nframe = randint(5, 20)
        coords = np.random.rand(nframe, natom, 3) * 20 - 10
        coords = np.array(coords, dtype='f')
        nctraj = NetCDFFile(get_fn('test2.nc', written=True), 'w', mmap=False)
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
        traj2 = NetCDFTraj.open_old(get_fn('test2.nc', written=True))
        np.testing.assert_allclose(traj2.coordinates, coords)

        # Now write NetCDF restart without any optional attributes/variables and
        # check proper reading
        natom = randint(100, 1000)
        nframe = randint(5, 20)
        coords = np.random.rand(natom, 3) * 20 - 10
        nctraj = NetCDFFile(get_fn('test2.ncrst', written=True), 'w', mmap=False)
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
        traj2 = NetCDFRestart.open_old(get_fn('test2.ncrst', written=True))
        np.testing.assert_allclose(traj2.coordinates[0], coords)

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
        self.assertTrue(all([round(x-y, 7) == 0 
                             for x, y in zip(rst.cell_lengths, rst.box[:3])]))
        self.assertTrue(all([round(x-y, 7) == 0 
                             for x, y in zip(rst.cell_angles, rst.box[3:])]))
