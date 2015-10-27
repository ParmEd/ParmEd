""" A set of miscellaneous tests for ParmEd based on coverage report """

import parmed as pmd
import unittest

class InitTestCase(unittest.TestCase):
    """ Test various package __init__ code """

    def testParmedVersion(self):
        """ Tests parmed.version class """
        version = type(pmd.version)

        v2_0_0 = version(2, 0, 0)
        v2_0_0.beta = None
        v2_0_0b1 = version(2, 0, 0)
        v2_0_0b1.beta = 1
        v2_0_0b2 = version(2, 0, 0)
        v2_0_0b2.beta = 2

        self.assertGreater(v2_0_0, (1, 9, 9))
        self.assertGreater(v2_0_0, (1, 9, 9, 1))
        self.assertGreater(v2_0_0, (2, 0, 0, 1))
        self.assertLess(v2_0_0, (2, 0, 1))
        self.assertLess(v2_0_0, (2, 1, 0))
        self.assertLess(v2_0_0, (2, 0, 1, 1))
        self.assertEqual(v2_0_0, (2, 0, 0))

        self.assertGreaterEqual(v2_0_0, (1, 9, 9))
        self.assertGreaterEqual(v2_0_0, (1, 9, 9, 1))
        self.assertGreaterEqual(v2_0_0, (2, 0, 0, 1))
        self.assertGreaterEqual(v2_0_0, (2, 0, 0))
        self.assertLessEqual(v2_0_0, (2, 0, 1))
        self.assertLessEqual(v2_0_0, (2, 1, 0))
        self.assertLessEqual(v2_0_0, (2, 0, 1, 1))
        self.assertLessEqual(v2_0_0, (2, 0, 0))

        self.assertGreater(v2_0_0b1, (1, 9, 9))
        self.assertGreater(v2_0_0b1, (1, 9, 9, 1))
        self.assertGreater(v2_0_0b2, (2, 0, 0, 1))
        self.assertLess(v2_0_0b1, (2, 0, 0))
        self.assertLess(v2_0_0b1, (2, 0, 0, 2))
        self.assertLess(v2_0_0b1, (2, 0, 1))
        self.assertLess(v2_0_0b1, (2, 0, 1, 1))
        self.assertLess(v2_0_0b2, (2, 0, 0, 3))
        self.assertEqual(v2_0_0b1, (2, 0, 0, 1))
        self.assertEqual(v2_0_0b2, (2, 0, 0, 2))

        self.assertGreaterEqual(v2_0_0b1, (1, 9, 9))
        self.assertGreaterEqual(v2_0_0b1, (1, 9, 9, 1))
        self.assertGreaterEqual(v2_0_0b2, (2, 0, 0, 1))
        self.assertGreaterEqual(v2_0_0b1, (2, 0, 0, 1))
        self.assertGreaterEqual(v2_0_0b2, (2, 0, 0, 2))
        self.assertLessEqual(v2_0_0b1, (2, 0, 0))
        self.assertLessEqual(v2_0_0b1, (2, 0, 0, 2))
        self.assertLessEqual(v2_0_0b1, (2, 0, 1))
        self.assertLessEqual(v2_0_0b1, (2, 0, 1))
        self.assertLessEqual(v2_0_0b2, (2, 0, 0, 3))
        self.assertLessEqual(v2_0_0b1, (2, 0, 0, 1))
        self.assertLessEqual(v2_0_0b2, (2, 0, 0, 2))

        self.assertNotEqual(v2_0_0, (1, 9, 9))
        self.assertNotEqual(v2_0_0, (1, 9, 9, 1))
        self.assertNotEqual(v2_0_0, (2, 0, 0, 1))
        self.assertNotEqual(v2_0_0, (2, 0, 1))
        self.assertNotEqual(v2_0_0, (2, 1, 0))
        self.assertNotEqual(v2_0_0, (2, 0, 1, 1))

        self.assertNotEqual(v2_0_0b1, (1, 9, 9))
        self.assertNotEqual(v2_0_0b1, (1, 9, 9, 1))
        self.assertNotEqual(v2_0_0b2, (2, 0, 0, 1))
        self.assertNotEqual(v2_0_0b1, (2, 0, 0))
        self.assertNotEqual(v2_0_0b1, (2, 0, 0, 2))
        self.assertNotEqual(v2_0_0b1, (2, 0, 1))
        self.assertNotEqual(v2_0_0b1, (2, 0, 1))
        self.assertNotEqual(v2_0_0b2, (2, 0, 0, 3))

        self.assertEqual(repr(v2_0_0), 'version(major=2, minor=0, patchlevel=0)')
        self.assertEqual(repr(v2_0_0b1), 'version(major=2, minor=0, patchlevel=0, beta=1)')

        self.assertIs(version.__lt__(v2_0_0, 1), NotImplemented)
        self.assertIs(version.__gt__(v2_0_0, 1), NotImplemented)
        self.assertIs(version.__eq__(v2_0_0, 1), NotImplemented)
