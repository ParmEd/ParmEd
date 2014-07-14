"""
Tests the chemistry/system module
"""
import chemistry.system as system
import unittest
from utils import get_fn

class TestChemistrySystem(unittest.TestCase):
    
    def setUp(self):
        self.pdb = get_fn('4lzt.pdb')
        self.pdbgz = get_fn('4lzt.pdb.gz')
        self.pdbbz2 = get_fn('4lzt.pdb.bz2')
        self.models = get_fn('2koc.pdb')

    def testAscii(self):
        """ Test PDB file parsing """
        self._check4lyt(system.ChemicalSystem.load_from_pdb(self.pdb))
        # The PDB file with multiple models
        pdbfile = system.ChemicalSystem.load_from_open_pdb(open(self.models))
        self.assertEqual(pdbfile.models, 20)
        self.assertEqual(pdbfile.positions(1)[:3], [-8.886, -5.163, 9.647])
        self.assertEqual(pdbfile.positions(20)[-3:], [-12.051, 5.205, -2.146])

    def testGzip(self):
        """ Test Gzipped-PDB file parsing """
        self._check4lyt(system.ChemicalSystem.load_from_pdb(self.pdbgz))

    def testBzip(self):
        """ Test Bzipped-PDB file parsing """
        self._check4lyt(system.ChemicalSystem.load_from_pdb(self.pdbbz2))

    # Private helper test functions
    def _check4lyt(self, obj):
        self.assertEqual(obj.models, 1)
        self.assertEqual(obj.box,
                         [27.24, 31.87, 34.23, 88.52, 108.53, 111.89])
        self.assertEqual(obj.experimental, 'X-RAY DIFFRACTION')
        self.assertEqual(len(obj), 274)
        self.assertEqual(obj.pmid, '9761848')
