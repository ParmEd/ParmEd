"""
Tests the routines in the parmed.geometry module
"""
from __future__ import division
import utils

from parmed import Atom
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

    def test_distance2(self):
        """ Tests the distance2 calculation """
        a1, a2 = Atom(), Atom()
        a1.xx = a1.xy = a1.xz = 0
        a2.xx = 3
        a2.xy = 4
        a2.xz = 0
        self.assertEqual(geo.distance2(a1, a2), 25)
        # Make sure it also works for tuples and a mixture of both
        self.assertEqual(geo.distance2((0, 0, 0), (3, 4, 0)), 25)
        self.assertEqual(geo.distance2(a1, (3, 4, 0)), 25)
        self.assertEqual(geo.distance2((0, 0, 0), a2), 25)

    def test_angle(self):
        """ Tests the angle calculation """
        a1, a2, a3 = Atom(), Atom(), Atom()
        a1.xx = a1.xy = a1.xz = 0
        a2.xx, a2.xy, a2.xz = 1, 0, 0
        a3.xx, a3.xy, a3.xz = 1, 1, 0
        self.assertAlmostEqual(geo.angle(a1, a2, a3), 90)
        # Check pathological cases
        a3.xx, a3.xy, a3.xz = 2, 0, 0
        self.assertAlmostEqual(geo.angle(a1, a2, a3), 180)
        a3.xx = a3.xy = a3.xz = 0
        self.assertAlmostEqual(geo.angle(a1, a2, a3), 0)

    def test_dihedral(self):
        """ Tests calculating a torsion angle between 4 points """
        a1, a2, a3, a4 = Atom(), Atom(), Atom(), Atom()
        a1.xx, a1.xy, a1.xz = 1, 0, 0
        a2.xx, a2.xy, a2.xz = 0, 0, 0
        a3.xx, a3.xy, a3.xz = 0, 0, 1
        a4.xx, a4.xy, a4.xz = 0.1, 0.6, 1
        self.assertAlmostEqual(geo.dihedral(a1, a2, a3, a4), 80.537677791974374)
        self.assertAlmostEqual(
                geo.dihedral([1, 0, 0], [0, 0, 0], [0, 0, 1], [0.1, 0.6, 1]),
                80.537677791974374
        )

    def test_box_lengths_vectors(self):
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

    def test_center_of_mass(self):
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
