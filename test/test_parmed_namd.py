"""
Tests the parmed/namd module
"""
import unittest
import os

from utils import get_fn, has_numpy
import parmed.namd as namd

class TestNamdBin(unittest.TestCase):
    """Test the NamdBinCoor class."""
    def testRead(self):
        coor = namd.NamdBinCoor.read(get_fn('ala_ala_ala.coor'))
        self.assertEqual(coor.natom,33)

class TestNamdVel(unittest.TestCase):
     """Test the NamdBinVel class."""
     def testRead(self):
        vel = namd.NamdBinVel.read(get_fn('ala_ala_ala.vel'))
        self.assertEqual(vel.natom,33)
        vel.write('test.vel')

if __name__ == '__main__':
    unittest.main()
