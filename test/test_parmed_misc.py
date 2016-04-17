""" A set of miscellaneous tests for ParmEd based on coverage report """

import parmed as pmd
import unittest

class InitTestCase(unittest.TestCase):
    """ Test various package __init__ code """

    def test_parmed_version(self):
        """ Tests parmed.version class """
        version = type(pmd.version)

        v2_0_0 = version(2, 0, 0, 0)
        v2_0_0p1 = version(2, 0, 0, 1)
        v2_0_0p2 = version(2, 0, 0, 2)
        v2_0_0.git_hash = v2_0_0p1.git_hash = v2_0_0p2.git_hash = 'abc123def'
        v2_0_0.dirty = v2_0_0p1.dirty = v2_0_0p2.dirty = False

        self.assertGreater(v2_0_0, (1, 9, 9))
        self.assertGreater(v2_0_0, (1, 9, 9, 1))
        self.assertLess(v2_0_0, (2, 0, 0, 1))
        self.assertLess(v2_0_0, (2, 0, 1))
        self.assertLess(v2_0_0, (2, 1, 0))
        self.assertLess(v2_0_0, (2, 0, 1, 1))
        self.assertEqual(v2_0_0, (2, 0, 0))

        self.assertGreaterEqual(v2_0_0, (1, 9, 9))
        self.assertGreaterEqual(v2_0_0, (1, 9, 9, 1))
        self.assertLessEqual(v2_0_0, (2, 0, 0, 1))
        self.assertGreaterEqual(v2_0_0, (2, 0, 0))
        self.assertLessEqual(v2_0_0, (2, 0, 1))
        self.assertLessEqual(v2_0_0, (2, 1, 0))
        self.assertLessEqual(v2_0_0, (2, 0, 1, 1))
        self.assertLessEqual(v2_0_0, (2, 0, 0))

        self.assertGreater(v2_0_0p1, (1, 9, 9))
        self.assertGreater(v2_0_0p1, (1, 9, 9, 1))
        self.assertGreater(v2_0_0p2, (2, 0, 0, 1))
        self.assertGreater(v2_0_0p1, (2, 0, 0))
        self.assertLess(v2_0_0p1, (2, 0, 0, 2))
        self.assertLess(v2_0_0p1, (2, 0, 1))
        self.assertLess(v2_0_0p1, (2, 0, 1, 1))
        self.assertLess(v2_0_0p2, (2, 0, 0, 3))
        self.assertEqual(v2_0_0p1, (2, 0, 0, 1))
        self.assertEqual(v2_0_0p2, (2, 0, 0, 2))

        self.assertGreaterEqual(v2_0_0p1, (1, 9, 9))
        self.assertGreaterEqual(v2_0_0p1, (1, 9, 9, 1))
        self.assertGreaterEqual(v2_0_0p2, (2, 0, 0, 1))
        self.assertGreaterEqual(v2_0_0p1, (2, 0, 0, 1))
        self.assertGreaterEqual(v2_0_0p2, (2, 0, 0, 2))
        self.assertGreaterEqual(v2_0_0p1, (2, 0, 0))
        self.assertLessEqual(v2_0_0p1, (2, 0, 0, 2))
        self.assertLessEqual(v2_0_0p1, (2, 0, 1))
        self.assertLessEqual(v2_0_0p1, (2, 0, 1))
        self.assertLessEqual(v2_0_0p2, (2, 0, 0, 3))
        self.assertLessEqual(v2_0_0p1, (2, 0, 0, 1))
        self.assertLessEqual(v2_0_0p2, (2, 0, 0, 2))

        self.assertNotEqual(v2_0_0, (1, 9, 9))
        self.assertNotEqual(v2_0_0, (1, 9, 9, 1))
        self.assertNotEqual(v2_0_0, (2, 0, 0, 1))
        self.assertNotEqual(v2_0_0, (2, 0, 1))
        self.assertNotEqual(v2_0_0, (2, 1, 0))
        self.assertNotEqual(v2_0_0, (2, 0, 1, 1))

        self.assertNotEqual(v2_0_0p1, (1, 9, 9))
        self.assertNotEqual(v2_0_0p1, (1, 9, 9, 1))
        self.assertNotEqual(v2_0_0p2, (2, 0, 0, 1))
        self.assertNotEqual(v2_0_0p1, (2, 0, 0))
        self.assertNotEqual(v2_0_0p1, (2, 0, 0, 2))
        self.assertNotEqual(v2_0_0p1, (2, 0, 1))
        self.assertNotEqual(v2_0_0p1, (2, 0, 1))
        self.assertNotEqual(v2_0_0p2, (2, 0, 0, 3))

        self.assertEqual(repr(v2_0_0), 'version(major=2, minor=0, patchlevel=0, commits_ahead=0)')
        self.assertEqual(repr(v2_0_0p1), 'version(major=2, minor=0, patchlevel=0, commits_ahead=1)')

        self.assertIs(version.__lt__(v2_0_0, 1), NotImplemented)
        self.assertIs(version.__gt__(v2_0_0, 1), NotImplemented)
        self.assertIs(version.__eq__(v2_0_0, 1), NotImplemented)
