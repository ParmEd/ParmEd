"""
Tests the routines in the parmed.geometry module
"""
from __future__ import division
import utils

from parmed import geometry as geo
from parmed import unit as u
from parmed.utils.six.moves import zip
import math
import unittest
import numpy as np

class TestChemistryGeometry(unittest.TestCase):
    """ Tests the various routines in the geometry package """

    def assertEqualVectors(self, a, b):
        a = strip_units(a)
        b = strip_units(b)
        for x, y in zip(a, b):
            self.assertEqual(x, y)

    def testBoxLengthsVectors(self):
        """ Test converting box lengths/angles to vectors and back again """
        a, b, c = geo.box_lengths_and_angles_to_vectors(1, 1, 1, 90, 90, 90)

        self.assertEqualVectors(a, [1.0, 0.0, 0.0] * u.angstroms)
        self.assertEqualVectors(b, [0.0, 1.0, 0.0] * u.angstroms)
        self.assertEqualVectors(c, [0.0, 0.0, 1.0] * u.angstroms)

        ang = 109.475
        rad = ang * math.pi / 180
        a,b,c = geo.box_lengths_and_angles_to_vectors(50, 50, 50, ang, ang, ang)
        leng, ang = geo.box_vectors_to_lengths_and_angles(a, b, c)
        self.assertEqualVectors(leng, (50, 50, 50))
        self.assertEqualVectors(ang, (rad, rad, rad))

    def testCenterOfMass(self):
        """ Tests the center-of-mass calculator """
        almost_equal = np.testing.assert_array_almost_equal
        # Make some systems we know the answer for
        array = np.asarray([[1, 0, 0], [-1, 0, 0]])
        masses = np.asarray([1, 1])
        self.assertTrue(np.all(geo.center_of_mass(array, masses) == 0))
        # Change masses
        masses = np.asarray([2, 1])
        almost_equal(geo.center_of_mass(array, masses), np.array([1/3, 0, 0]))

def strip_units(x):
    if u.is_quantity(x):
        return x.value_in_unit_system(u.akma_unit_system)
    return x

if __name__ == '__main__':
    unittest.main()
