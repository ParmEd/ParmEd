from chemistry import read_PDB, load_rosetta
from utils import get_fn
import unittest
from unittest import skipIf, skipTest

try:
    from simtk.openmm.app import PDBFile
except:
    PDBFile = None
    pass
try:
    from rosetta import init, pose_from_sequence
except ImportError:
    init = pose_from_sequence = None
    pass


def _unpackLen(obj):
    return len(list(obj))


@skipTest(not init, "Cannot test load_rosetta module without PyRosetta.")
class TestRosetta(unittest.TestCase):
    """ Tests loading of a Rosetta pose object """

    def testLoadStruct(self):
        """ Test load_rosetta against read_PDB"""

        init()
        pose = pose_from_sequence(3*'A')

        struct = load_rosetta(pose)
        pdb = read_PDB(get_fn('ala_ala_ala.pdb'))

        self.assertEqual(len(struct.atoms), len(pdb.atoms))
        self.assertEqual(len(struct.residues), len(pdb.residues))

    @skipIf(not PDBFile, "Cannot compare topologies without OpenMM.")
    def testLoadedTopology(self):
        """ Test load_rosetta against OpenMM topology"""

        init()
        pose = pose_from_sequence(3*'A')

        struct = load_rosetta(pose)
        pdb = PDBFile(get_fn('ala_ala_ala.pdb'))

        self.assertEqual(_unpackLen(struct.topology.atoms()),
                         _unpackLen(pdb.topology.atoms()))

        self.assertEqual(_unpackLen(struct.topology.bonds()),
                         _unpackLen(pdb.topology.bonds()))

        self.assertEqual(_unpackLen(struct.topology.residues()),
                         _unpackLen(pdb.topology.residues()))


if __name__ == '__main__':
    unittest.main()
