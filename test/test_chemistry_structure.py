"""
Tests the chemistry/structure module
"""
import chemistry.structure as structure
import unittest
from utils import get_fn

class TestChemistryStructure(unittest.TestCase):
    
    def setUp(self):
        self.pdb = get_fn('4lzt.pdb')
        self.pdbgz = get_fn('4lzt.pdb.gz')
        self.pdbbz2 = get_fn('4lzt.pdb.bz2')
        self.models = get_fn('2koc.pdb')
        self.overflow = get_fn('4lyt_vmd.pdb')

    def testAscii(self):
        """ Test PDB file parsing """
        self._check4lyt(structure.read_PDB(self.pdb))
        # The PDB file with multiple models
        pdbfile = structure.read_PDB(open(self.models))
        self.assertEqual(len(pdbfile.pdbxyz), 20)
        self.assertEqual(pdbfile.pdbxyz[0][:3], [-8.886, -5.163, 9.647])
        self.assertEqual(pdbfile.pdbxyz[19][-3:], [-12.051, 5.205, -2.146])

    def testGzip(self):
        """ Test Gzipped-PDB file parsing """
        self._check4lyt(structure.read_PDB(self.pdbgz))

    def testBzip(self):
        """ Test Bzipped-PDB file parsing """
        self._check4lyt(structure.read_PDB(self.pdbbz2))

    def testVmdOverflow(self):
        """ Test PDB file where atom and residue numbers overflow """
        pdbfile = structure.read_PDB(self.overflow)
        self.assertEqual(len(pdbfile.atoms), 110237)
        self.assertEqual(len(pdbfile.residues), 35697)
        self.assertEqual(pdbfile.box, [0, 0, 0, 90, 90, 90])

    # Private helper test functions
    def _check4lyt(self, obj):
        self.assertEqual(len(obj.pdbxyz), 1)
        self.assertEqual(obj.box,
                         [27.24, 31.87, 34.23, 88.52, 108.53, 111.89])
        self.assertEqual(len(obj.atoms), 1164)
        self.assertEqual(len(obj.residues[0]), 9)
        # Check that alternate conformations are taken into account
        total_natoms = 0
        for i, atom in enumerate(obj.atoms):
            total_natoms += 1
            for key in atom.other_locations:
                total_natoms += 1
                atom2 = atom.other_locations[key]
                self.assertEqual(atom.altloc, 'A')
                self.assertEqual(atom2.altloc, 'B')
                if i in [388, 389]:
                    # Something wrong with atoms 388 and 389 -- sum of their
                    # occupancies are 1.02 :)
                    self.assertEqual(atom2.occupancy + atom.occupancy, 1.02)
                else:
                    self.assertEqual(atom2.occupancy + atom.occupancy, 1)
        self.assertEqual(total_natoms, 1183)
        # Check the metadata
        self.assertEqual(obj.experimental, 'X-RAY DIFFRACTION')
        self.assertEqual(len(obj.residues), 274)
        self.assertEqual(obj.pmid, '9761848')
        self.assertEqual(obj.journal_authors, 'M.A.WALSH,T.R.SCHNEIDER,'
                         'L.C.SIEKER,Z.DAUTER,V.S.LAMZIN,K.S.WILSON')
        self.assertEqual(obj.journal, 'ACTA CRYSTALLOGR.,SECT.D')
        self.assertEqual(obj.year, 1998)
        self.assertEqual(obj.keywords, ['HYDROLASE', 'O-GLYCOSYL',
                         'GLYCOSIDASE'])
        self.assertEqual(obj.title, 'REFINEMENT OF TRICLINIC HEN EGG-WHITE '
                         'LYSOZYME AT ATOMIC RESOLUTION.')
        self.assertEqual(obj.doi, '10.1107/S0907444997013656')
        self.assertEqual(obj.volume, '54')
        self.assertEqual(obj.page, '522')
        # Check the TER card is picked up
        for i, residue in enumerate(obj.residues):
            if i == 128:
                self.assertTrue(residue.ter)
            else:
                self.assertFalse(residue.ter)
