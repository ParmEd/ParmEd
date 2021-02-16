""" Tests the Vec3 object """
from unittest import TestCase, skipIf
from parmed import Vec3
from utils import has_old_vec3

@skipIf(has_old_vec3(), 'Vec3 is taken from OpenMM with an old implementation. Tests will fail')
class TestVec3(TestCase):
    """ Tests the Vec3 type """

    def test_attributes(self):
        """ Tests that the Vec3 object has x, y, and z accessors """
        vec1 = Vec3(1, 2, 3)
        self.assertEqual(vec1.x, 1)
        self.assertEqual(vec1.y, 2)
        self.assertEqual(vec1.z, 3)

    def test_negation(self):
        """ Tests that the negation unary operator works on Vec3 """
        vec1 = Vec3(1, 2, 3)
        vec1_neg = Vec3(-1, -2, -3)
        self.assertEqual(-vec1, vec1_neg)

    def test_equality(self):
        """ Tests the == operator for Vec3 """
        vec1 = Vec3(1, 2, 3)
        vec2 = Vec3(1, 2, 3)
        self.assertEqual(vec1, vec2)

    def test_addition(self):
        vec1 = Vec3(1, 2, 3)
        vec2 = Vec3(4, 5, 6)
        vec2_tup = (4, 5, 6)
        result = Vec3(5, 7, 9)
        self.assertEqual(vec1 + vec2, result)
        self.assertEqual(vec1 + vec2_tup, result)

    def test_subtraction(self):
        vec1 = Vec3(1, 2, 3)
        vec2 = Vec3(3, 2, 1)
        vec2_tup = (3, 2, 1)
        result = Vec3(-2, 0, 2)
        self.assertEqual(vec1 - vec2, result)
        self.assertEqual(vec1 - vec2_tup, result)
        self.assertEqual(vec2_tup - vec1, -result)

    def test_multiplication(self):
        vec1 = Vec3(1, 2, 3)
        factor = 2
        result = Vec3(2, 4, 6)
        self.assertEqual(vec1 * factor, result)
        self.assertEqual(factor * vec1, result)

    def test_division(self):
        vec1 = Vec3(4, 5, 6)
        factor = 2
        result = Vec3(2, 2.5, 3)
        self.assertEqual(vec1 / factor, result)
        with self.assertRaises(TypeError):
            2 / vec1
