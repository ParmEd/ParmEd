"""
Tests the functionality in the chemistry.amber package
"""

from array import array
from chemistry.amber import readparm, asciicrd, mask
import os
import unittest
from utils import get_fn, has_numpy

class TestWriteFiles(unittest.TestCase):
    
    def setUp(self):
        try:
            os.makedirs(get_fn('writes'))
        except OSError:
            pass

    def tearDown(self):
        try:
            for f in os.listdir(get_fn('writes')):
                os.unlink(get_fn(f, written=True))
            os.rmdir(get_fn('writes'))
        except OSError:
            pass

    def testWriteAmberParm(self):
        parm = readparm.AmberParm(get_fn('trx.prmtop'))
        parm.writeParm(get_fn('trx.prmtop', written=True))
        f1 = open(get_fn('trx.prmtop'), 'r')
        f2 = open(get_fn('trx.prmtop', written=True), 'r')
        try:
            for line1, line2 in zip(f1, f2):
                if line1.startswith('%VERSION'):
                    self.assertTrue(line2.startswith('%VERSION'))
                    continue
                self.assertEqual(line1.strip(), line2.strip())
        finally:
            f1.close()
            f2.close()

    def testAmberRestart(self):
        Restart = asciicrd.AmberAsciiRestart
        box = [10, 10, 10, 90, 90, 90]
        rst = Restart(get_fn('testc.rst7', written=True), 'w', natom=9)
        rst.coordinates = range(27)
        rst.close()
        rst = Restart(get_fn('testcv.rst7', written=True), 'w', natom=20)
        rst.coordinates = range(60)
        rst.velocities = list(reversed(range(60)))
        rst.close()
        rst = Restart(get_fn('testcb.rst7', written=True), 'w', natom=7)
        rst.coordinates = range(21)
        rst.box = box[:]
        rst.close()
        rst = Restart(get_fn('testcvb.rst7', written=True), 'w', natom=15)
        rst.coordinates = range(45)
        rst.velocities = list(reversed(range(45)))
        rst.box = box[:]
        rst.close()
        self._check_written_restarts(box)
    
    def testAmberRestartNumpy(self):
        import numpy as np
        Restart = asciicrd.AmberAsciiRestart
        box = np.asarray([10, 10, 10, 90, 90, 90])
        rst = Restart(get_fn('testc.rst7', written=True), 'w', natom=9)
        rst.coordinates = np.arange(27).reshape((9,3))
        rst.close()
        rst = Restart(get_fn('testcv.rst7', written=True), 'w', natom=20)
        rst.coordinates = np.arange(60).reshape((20,3))
        rst.velocities = np.asarray(list(reversed(range(60)))).reshape((20,3))
        rst.close()
        rst = Restart(get_fn('testcb.rst7', written=True), 'w', natom=7)
        rst.coordinates = np.arange(21).reshape((7,3))
        rst.box = box
        rst.close()
        rst = Restart(get_fn('testcvb.rst7', written=True), 'w', natom=15)
        rst.coordinates = np.arange(45).reshape((15,3))
        rst.velocities = np.asarray(list(reversed(range(45)))).reshape((15,3))
        rst.box = box
        rst.close()
        self._check_written_restarts(box)

    def testAmberMdcrd(self):
        box = [15, 15, 15]
        Mdcrd = asciicrd.AmberMdcrd
        crd = Mdcrd(get_fn('testc.mdcrd', written=True), natom=15, hasbox=False,
                    mode='w', title='Test file')
        crd.add_coordinates(range(45))
        crd.add_coordinates([x+1 for x in range(45)])
        crd.add_coordinates([x+2 for x in range(45)])
        crd.add_coordinates([x+3 for x in range(45)])
        crd.add_coordinates([x+4 for x in range(45)])
        crd.close()
        crd = Mdcrd(get_fn('testcb.mdcrd', written=True), natom=18, hasbox=True,
                    mode='w', title='Test file')
        crd.add_coordinates(range(54))
        crd.add_box(box)
        crd.add_coordinates([x+1 for x in range(54)])
        crd.add_box(box)                          
        crd.add_coordinates([x+2 for x in range(54)])
        crd.add_box(box)                          
        crd.add_coordinates([x+3 for x in range(54)])
        crd.add_box(box)                          
        crd.add_coordinates([x+4 for x in range(54)])
        crd.add_box(box)
        crd.close()
        self._check_written_mdcrds(box)

    def testAmberMdcrdNumpy(self):
        import numpy as np
        box = np.asarray([15, 15, 15])
        Mdcrd = asciicrd.AmberMdcrd
        crd = Mdcrd(get_fn('testc.mdcrd', written=True), natom=15, hasbox=False,
                    mode='w', title='Test file')
        coorddata = np.arange(45).reshape((15,3))
        crd.add_coordinates(coorddata)
        crd.add_coordinates(coorddata+1)
        crd.add_coordinates(coorddata+2)
        crd.add_coordinates(coorddata+3)
        crd.add_coordinates(coorddata+4)
        crd.close()
        crd = Mdcrd(get_fn('testcb.mdcrd', written=True), natom=18, hasbox=True,
                    mode='w', title='Test file')
        coorddata = np.arange(54).reshape((18,3))
        crd.add_coordinates(coorddata)
        crd.add_box(box)
        crd.add_coordinates(coorddata+1)
        crd.add_box(box)
        crd.add_coordinates(coorddata+2)
        crd.add_box(box)
        crd.add_coordinates(coorddata+3)
        crd.add_box(box)
        crd.add_coordinates(coorddata+4)
        crd.add_box(box)
        crd.close()
        self._check_written_mdcrds(box)

    def testBadFileUsage(self):
        Restart = asciicrd.AmberAsciiRestart
        Mdcrd = asciicrd.AmberMdcrd
        box = [10, 10, 10, 90, 90, 90]
        rst = Restart(get_fn('testc.rst7', written=True), 'w', natom=9,
                      hasbox=True)
        def assign(obj, stmnt):
            rst = crd = obj
            exec(stmnt)
        try:
            self.assertRaises(ValueError,
                              lambda: assign(rst, 'rst.coordinates=range(20)'))
            self.assertRaises(RuntimeError,
                              lambda: assign(rst, 'rst.box=[10]*3+[90]*3'))
            rst.coordinates = range(27)
            rst.box = box
            self.assertRaises(RuntimeError,
                              lambda: assign(rst, 'rst.velocities=range(27)'))
        finally:
            rst.close()
        crd = Mdcrd(get_fn('testc.mdcrd', written=True), natom=15, hasbox=True,
                    mode='w', title='Test file')
        s = 'range(45)'
        s2 = 'range(42)'
        try:
            crd.add_coordinates(eval(s))
            self.assertRaises(RuntimeError,
                              lambda: assign(crd, 'crd.add_coordinates(%s)'%s))
            crd.add_box([10, 10, 10])
            self.assertRaises(ValueError,
                              lambda: assign(crd, 'crd.add_coordinates(%s)'%s2))
        finally:
            crd.close()

    def _check_written_restarts(self, box):
        # Now try to read them and verify the information (keep in mind that the
        # restart velocities are scaled down then back up, so you'll need to use
        # assertAlmostEqual in this case).
        rst = readparm.Rst7.open(get_fn('testc.rst7', written=True))
        self.assertFalse(rst.hasbox)
        self.assertFalse(rst.hasvels)
        for x1, x2 in zip(rst.coordinates, range(27)):
            self.assertEqual(x1, x2)
        rst = readparm.Rst7.open(get_fn('testcb.rst7', written=True))
        self.assertTrue(rst.hasbox)
        self.assertFalse(rst.hasvels)
        for x1, x2 in zip(rst.coordinates, range(21)):
            self.assertEqual(x1, x2)
        for x1, x2 in zip(rst.box, box):
            self.assertEqual(x1, x2)
        rst = readparm.Rst7.open(get_fn('testcv.rst7', written=True))
        self.assertTrue(rst.hasvels)
        self.assertFalse(rst.hasbox)
        for x1, x2 in zip(rst.coordinates, range(60)):
            self.assertEqual(x1, x2)
        for x1, x2 in zip(rst.velocities, reversed(range(60))):
            self.assertAlmostEqual(x1, x2, places=5)
        rst = readparm.Rst7.open(get_fn('testcvb.rst7', written=True))
        self.assertTrue(rst.hasvels)
        self.assertTrue(rst.hasbox)
        for x1, x2 in zip(rst.coordinates, range(45)):
            self.assertEqual(x1, x2)
        for x1, x2 in zip(rst.velocities, reversed(range(45))):
            self.assertAlmostEqual(x1, x2, places=5)
        for x1, x2 in zip(rst.box, box):
            self.assertEqual(x1, x2)

    def _check_written_mdcrds(self, box):
        # Now try to read them and verify the information
        crd = asciicrd.AmberMdcrd(get_fn('testc.mdcrd', written=True),
                                  15, False, 'r')
        self.assertEqual(crd.title, 'Test file')
        self.assertFalse(crd.hasbox)
        for i in range(crd.frame):
            for x1, x2 in zip(crd.coordinates(i), [x+i for x in range(45)]):
                self.assertEqual(x1, x2)
        for i, array in enumerate(crd.coordinates()):
            for x1, x2 in zip(array, [x+i for x in range(45)]):
                self.assertEqual(x1, x2)
        crd.close()

        crd = asciicrd.AmberMdcrd(get_fn('testcb.mdcrd', written=True),
                                  18, True, 'r')
        self.assertEqual(crd.title, 'Test file')
        self.assertTrue(crd.hasbox)
        for i in range(crd.frame):
            for x1, x2 in zip(crd.coordinates(i), [x+i for x in range(54)]):
                self.assertEqual(x1, x2)
            for x1, x2 in zip(crd.box(i), box):
                self.assertEqual(x1, x2)
        for i, (coords, mybox) in enumerate(zip(crd.coordinates(), crd.box())):
            for x1, x2 in zip(coords, [x+i for x in range(45)]):
                self.assertEqual(x1, x2)
            for x1, x2 in zip(mybox, box):
                self.assertEqual(x1, x2)

if not has_numpy():
    del TestWriteFiles.testAmberRestartNumpy, TestWriteFiles.testAmberMdcrdNumpy
