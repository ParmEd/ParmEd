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

@unittest.skipUnless(has_rdkit, "Only test load_rdkit module on Linux")
class TestRDKit(unittest.TestCase):
    """ Tests loading of an rdkit Mol object """

    def test_load_rdkit_mol(self):
        """ test load rdkit from Mol """
        from rdkit import Chem
        m1 = Chem.MolFromSmiles('C1=CC=CN=C1')
        parm = pmd.load_rdkit(m1)
        self.assertEqual([atom.name for atom in parm.atoms], ['C1', 'C2', 'C3', 'C4', 'N1', 'C5']) 
        self.assertEqual(parm.residues[0].name, 'UNL')

    def test_load_smiles(self):
        """ test load rdkit from smiles string """
        smiles = 'C1=CC=CN=C1'

        # coordinates = False
        parm = pmd.rdkit.from_smiles(smiles, coordinates=False, hydrogens=False)
        self.assertEqual([atom.name for atom in parm.atoms], ['C1', 'C2', 'C3', 'C4', 'N1', 'C5']) 
        self.assertEqual(parm.residues[0].name, 'UNL')
        self.assertIs(parm.coordinates, None)
        self.assertIs(parm.get_coordinates(), None)

        # coordinates = True (default)
        parm = pmd.rdkit.from_smiles(smiles, coordinates=True, hydrogens=False)
        np.testing.assert_allclose(parm.coordinates[0], [-1.072,  0.829 ,  0.108])

    def test_load_smiles_explicit_hydrogen(self):
        """ test adding explict hydrogens from smiles string"""
        smiles = "CC"

        num_atoms = 8 # CH3-CH3
        parm = pmd.rdkit.from_smiles(smiles=smiles, coordinates=False, hydrogens=True)
        self.assertEqual(len(parm.atoms), num_atoms)
        self.assertEqual([atom.name for atom in parm.atoms], ['C1', 'C2', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6'])

    def test_to_mol(self):
        """ Test converting a Structure to an RDKit Mol object """
        from rdkit.Chem import Mol

        structure = pmd.load_file(get_fn("4lzt.pdb"))
        mol = structure.rdkit_mol

        self.assertIsInstance(mol, Mol)
        self.assertEqual(mol.GetNumAtoms(), len(structure.atoms))
        self.assertEqual(mol.GetNumBonds(), len(structure.bonds))
