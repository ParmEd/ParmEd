from __future__ import print_function
import sys
import parmed as pmd
import unittest
import numpy as np
from utils import get_fn

try:
    import rdkit
    has_rdkit = True
except ImportError:
    has_rdkit = False

is_linux = sys.platform.startswith('linux')

@unittest.skipUnless(has_rdkit and is_linux, "Only test load_rdkit module on Linux")
class TestRDKit(unittest.TestCase):
    """ Tests loading of an rdkit Mol object """

    def test_load_rdkit_mol(self):
        """ test load rdkit from Mol """
        from rdkit import Chem
        m1 = Chem.MolFromSmiles('C1=CC=CN=C1')
        parm = pmd.load_rdkit(m1)
        self.assertEqual([atom.name for atom in parm.atoms], 
                         ['C1', 'C2', 'C3', 'C4', 'N1', 'C5']) 
        self.assertEqual(parm.residues[0].name, 'UNL')

    def test_load_smiles(self):
        """ test load rdkit from smiles string """
        smiles = 'C1=CC=CN=C1'

        # coordinates = False
        parm = pmd.rdkit.from_smiles(smiles, coordinates=False)
        self.assertEqual([atom.name for atom in parm.atoms], 
                         ['C1', 'C2', 'C3', 'C4', 'N1', 'C5']) 
        self.assertEqual(parm.residues[0].name, 'UNL')
        self.assertIs(parm.coordinates, None)
        self.assertIs(parm.get_coordinates(), None)

        # coordinates = True (default)
        parm = pmd.rdkit.from_smiles(smiles)
        np.testing.assert_allclose(parm.coordinates[0], [-1.076,  0.83 ,  0.011])

    def test_load_sdf(self):
        """ test load rdkit from SDF format """
        sdffile = get_fn('test.sdf')
        parmlist = pmd.rdkit.from_sdf(sdffile)
        self.assertEqual(len(parmlist[0].atoms), 34)
        self.assertEqual(len(parmlist[1].atoms), 43)
        print(parmlist[0].coordinates[0])
        np.testing.assert_almost_equal(parmlist[0].coordinates[0], [2.0000, 2.7672, 0.0000], decimal=3)
        np.testing.assert_almost_equal(parmlist[0].coordinates[-1], [9.9858, -2.8473, 0.0000], decimal=3)
        np.testing.assert_almost_equal(parmlist[1].coordinates[0], [7.0468, -1.7307, 0.0000], decimal=3)
        np.testing.assert_almost_equal(parmlist[1].coordinates[-1], [1.5269, 2.1331, 0.0000], decimal=3)
