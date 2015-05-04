"""
Tests the NetCDF file parsing capabilities with the different backends
"""

from chemistry.utils.six.moves import range, zip
import unittest
import utils

get_fn = utils.get_fn

class TestNetCDF(unittest.TestCase):
    """ Test NetCDF Functionality """
    
    def testScipy(self):
        """ Test scipy NetCDF parsing """
        import chemistry.amber as amber
        if utils.has_scipy():
            amber.use('scipy')
            from chemistry.amber.netcdffiles import NetCDFTraj, NetCDFRestart
            traj = NetCDFTraj.open_old(get_fn('tz2.truncoct.nc'))
            self._check_traj(traj)
            rst = NetCDFRestart.open_old(get_fn('ncinpcrd.rst7'))
            self._check_rst(rst)
        else:
            self.assertRaises(ImportError, lambda: amber.use('scipy'))

    def testNetcdf4(self):
        """ Test netCDF4 parsing """
        import chemistry.amber as amber
        if utils.has_netcdf4():
            amber.use('netCDF4')
            from chemistry.amber.netcdffiles import NetCDFTraj, NetCDFRestart
            traj = NetCDFTraj.open_old(get_fn('tz2.truncoct.nc'))
            self._check_traj(traj)
            rst = NetCDFRestart.open_old(get_fn('ncinpcrd.rst7'))
            self._check_rst(rst)
        else:
            self.assertRaises(ImportError, lambda: amber.use('netCDF4'))

    def testScientificPython(self):
        """ Test ScientificPython parsing """
        import chemistry.amber as amber
        if utils.has_scientific():
            amber.use('Scientific')
            from chemistry.amber.netcdffiles import NetCDFTraj, NetCDFRestart
            traj = NetCDFTraj.open_old(get_fn('tz2.truncoct.nc'))
            self._check_traj(traj)
            rst = NetCDFRestart.open_old(get_fn('ncinpcrd.rst7'))
            self._check_rst(rst)
        else:
            self.assertRaises(ImportError, lambda: amber.use('Scientific'))

    def _check_traj(self, traj):
        """ Checks various trajectory properties """
        self.assertEqual(traj.Conventions, 'AMBER')
        self.assertEqual(traj.application, 'AMBER')
        self.assertEqual(traj.program, 'sander')
        self.assertEqual(traj.programVersion, '9.0')
        self.assertEqual(traj.ConventionVersion, '1.0')
        runsum = 0
        for frame in range(traj.frame):
            for coord in traj.coordinates(frame):
                runsum += coord
        self.assertAlmostEqual(runsum, 7049.7614416316874, places=4)

        for frame in range(traj.frame):
            leng, ang = traj.cell_lengths_angles(frame)
            self.assertAlmostEqual(ang[0], 109.471219, places=6)
            self.assertAlmostEqual(ang[1], 109.471219, places=6)
            self.assertAlmostEqual(ang[2], 109.471219, places=6)
            self.assertAlmostEqual(leng[0], leng[1])
            self.assertAlmostEqual(leng[0], leng[2])

    def _check_rst(self, rst):
        """ Checks various restart properties """
        self.assertAlmostEqual(rst.time, 30.1)
        self.assertEqual(rst.Conventions, 'AMBERRESTART')
        self.assertEqual(rst.title, 'ACE')
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

if __name__ == '__main__':
    unittest.main()
