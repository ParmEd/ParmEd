from __future__ import print_function
import parmed as pmd
import unittest

try:
    from rdkit import Chem
    has_rdkit = True
except ImportError:
    has_rdkit = False


@unittest.skipIf(not has_rdkit, "Cannot test load_rdkit module without rdkit.")
class TestRDKit(unittest.TestCase):
    """ Tests loading of an rdkit Mol object """

    def test_load_rdkit_mol(self):
        m1 = Chem.MolFromSmiles('C1=CC=CN=C1')
        parm = pmd.load_rdkit(m1)
        assert [atom.name for atom in parm.atoms] == ['C1', 'C2', 'C3', 'C4', 'N1', 'C5'] 
        assert parm.residues[0].name == 'UNL', 'resname must be UNL'

    def test_load_smile(self):
        smile = 'C1=CC=CN=C1'
        parm = pmd.rdkit.from_smile(smile)
        assert [atom.name for atom in parm.atoms] == ['C1', 'C2', 'C3', 'C4', 'N1', 'C5'] 
        assert parm.residues[0].name == 'UNL', 'resname must be UNL'

if __name__ == '__main__':
    unittest.main()
