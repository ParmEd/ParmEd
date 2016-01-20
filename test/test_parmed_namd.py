"""
Tests the parmed/namd module
"""
import unittest
import os

from utils import get_fn
import parmed.namd as namd

class TestNamdBin(unittest.TestCase):
    """Test the NamdBinFile classes."""
    def testCoorRead(self):
        """ Test reading NamdBinCoor file """
        coor = namd.NamdBinCoor.read(get_fn('ala_ala_ala.coor'))
        self.assertEqual(coor.natom,33)

    def testVelRead(self):
        """ Test reading NamdBinVel file """
        vel = namd.NamdBinVel.read(get_fn('ala_ala_ala.vel'))
        self.assertEqual(vel.natom,33)

if __name__ == '__main__':
    unittest.main()
