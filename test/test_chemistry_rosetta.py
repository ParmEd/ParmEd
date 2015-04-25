from chemistry import read_PDB, load_rosetta
from utils import get_fn
import unittest

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


class TestRosetta(unittest.TestCase):
    """ Tests loading of a Rosetta pose object """

    def testLoadPose(self):
        if not init or not pose_from_sequence:
            raise ImportError('Could not load the PyRosetta module.')

        init()
        pose = pose_from_sequence(3*'A')
        struct = load_rosetta(pose)
        if PDBFile:
            pdb = PDBFile(get_fn('ala_ala_ala.pdb'))
            self.assertEqual(_unpackLen(struct.topology.bonds()),
                             _unpackLen(pdb.topology.bonds()))
        else:
            pdb = read_PDB(get_fn('ala_ala_ala.pdb'))

        self.assertEqual(_unpackLen(struct.topology.atoms()),
                         _unpackLen(pdb.topology.atoms()))

        self.assertEqual(_unpackLen(struct.topology.residues()),
                         _unpackLen(pdb.topology.residues()))


if __name__ == '__main__':
    if init:
        unittest.main()
