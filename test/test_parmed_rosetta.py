from parmed import read_PDB, load_rosetta
from parmed.utils.six.moves import range
from itertools import chain
from utils import get_fn
import unittest
try:
    from simtk.openmm.app import PDBFile
except:
    PDBFile = None
try:
    from rosetta import init, pose_from_sequence
except ImportError:
    init = pose_from_sequence = None

@unittest.skipIf(init is None, "Cannot test load_rosetta module without PyRosetta.")
class TestRosetta(unittest.TestCase):
    """ Tests loading of a Rosetta pose object """

    def test_loaded_positions(self):
        """ Test that positions were properly loaded"""

        init()
        seq = 3*'A'
        pose = pose_from_sequence(seq)

        struct = load_rosetta(pose)

        posexyz = list(
            chain(*[[tuple(atom.xyz()) for atom in res.atoms()]
                    for res in [pose.residue(idx) for idx in range(1, len(seq)+1)]]))

        structxyz = [(atom.xx, atom.xy, atom.xz) for atom in struct.atoms]

        self.assertEqual(posexyz, structxyz)

    def test_load_struct(self):
        """ Test load_rosetta against read_PDB"""

        init()
        pose = pose_from_sequence(3*'A')

        struct = load_rosetta(pose)
        pdb = read_PDB(get_fn('ala_ala_ala.pdb'))

        self.assertEqual(len(struct.atoms), len(pdb.atoms))
        self.assertEqual(len(struct.bonds), len(pdb.bonds))
        self.assertEqual(len(struct.residues), len(pdb.residues))

    @unittest.skipIf(PDBFile is None, "Cannot compare topologies without OpenMM.")
    def test_loaded_topology(self):
        """ Test load_rosetta against OpenMM topology"""

        init()
        pose = pose_from_sequence(3*'A')

        struct = load_rosetta(pose)
        pdb = PDBFile(get_fn('ala_ala_ala.pdb'))

        self.assertEqual(len(list(struct.topology.atoms())),
                         len(list(pdb.topology.atoms())))

        self.assertEqual(len(list(struct.topology.bonds())),
                         len(list(pdb.topology.bonds())))

        self.assertEqual(len(list(struct.topology.residues())),
                         len(list(pdb.topology.residues())))
