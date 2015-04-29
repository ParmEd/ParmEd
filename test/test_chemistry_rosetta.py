from chemistry import read_PDB, load_rosetta
from itertools import chain
from utils import get_fn, skipIf
import unittest

try:
    from simtk.openmm.app import PDBFile
except:
    PDBFile = None
try:
    from rosetta import init, pose_from_sequence
except ImportError:
    init = pose_from_sequence = None


def _unpackLen(obj):
    return len(list(obj))


@skipIf(not init, "Cannot test load_rosetta module without PyRosetta.")
class TestRosetta(unittest.TestCase):
    """ Tests loading of a Rosetta pose object """

    def testLoadedPositions(self):
        """ Test that positions were properly loaded"""

        init()
        seq = 3*'A'
        pose = pose_from_sequence(seq)

        struct = load_rosetta(pose)

        posexyz = list(
            chain(*[[tuple(atom.xyz()) for atom in res.atoms()]
                    for res in [pose.residue(idx)
                                for idx in range(1, len(seq)+1)]]))

        structxyz = [(atom.xx, atom.xy, atom.xz) for atom in struct.atoms]

        self.assertEqual(posexyz, structxyz)

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
