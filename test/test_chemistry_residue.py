"""
Tests the functionality in chemistry.residue
"""

import utils
import unittest
from chemistry import residue

class TestChemistryResidue(unittest.TestCase):

    def testResidueMembers(self):
        """ Test that the main amino acid residues are members of residue """
        self.assertTrue(hasattr(residue, 'ASP'))
        self.assertTrue(hasattr(residue, 'GLU'))
        self.assertTrue(hasattr(residue, 'HIS'))
        self.assertTrue(hasattr(residue, 'GLN'))
        self.assertTrue(hasattr(residue, 'ASN'))
        self.assertTrue(hasattr(residue, 'LYS'))
        self.assertTrue(hasattr(residue, 'TYR'))
        self.assertTrue(hasattr(residue, 'TRP'))
        self.assertTrue(hasattr(residue, 'ILE'))
        self.assertTrue(hasattr(residue, 'LEU'))
        self.assertTrue(hasattr(residue, 'ALA'))
        self.assertTrue(hasattr(residue, 'GLY'))
        self.assertTrue(hasattr(residue, 'PRO'))
        self.assertTrue(hasattr(residue, 'VAL'))
        self.assertTrue(hasattr(residue, 'THR'))
        self.assertTrue(hasattr(residue, 'SER'))
        self.assertTrue(hasattr(residue, 'CYS'))
        self.assertTrue(hasattr(residue, 'MET'))
        self.assertTrue(hasattr(residue, 'ARG'))
        self.assertTrue(hasattr(residue, 'PHE'))

    def testNameLookup(self):
        """ Test that looking up residues by name works """
        self.assertIs(residue.AminoAcidResidue.get('Alanine'), residue.ALA)
        self.assertIs(residue.AminoAcidResidue.get('Glutamine'), residue.GLN)
        self.assertIs(residue.AminoAcidResidue.get('glutamate'), residue.GLU)
        self.assertIs(residue.AminoAcidResidue.get('HiStIDiNe'), residue.HIS)
        self.assertIs(residue.AminoAcidResidue.get('ISOLEUCINE'), residue.ILE)

    def testAbbrLookup(self):
        """ Test that looking up residues by abbreviation works """
        self.assertIs(residue.AminoAcidResidue.get('ALA'), residue.ALA)
        self.assertIs(residue.AminoAcidResidue.get('GLY'), residue.GLY)
        self.assertIs(residue.AminoAcidResidue.get('glu'), residue.GLU)
        self.assertIs(residue.AminoAcidResidue.get('Asp'), residue.ASP)

    def testSymbolLookup(self):
        """ Test that looking up residues by symbol works """
        self.assertIs(residue.AminoAcidResidue.get('E'), residue.GLU)
        self.assertIs(residue.AminoAcidResidue.get('D'), residue.ASP)
        self.assertIs(residue.AminoAcidResidue.get('A'), residue.ALA)

    def testReprOutput(self):
        """ Test the %r representation of the Amino Acids """
        for res in residue.AminoAcidResidue.all_residues:
            self.assertEqual('<Amino Acid Residue %s: %s [%s]>' % (res.name,
                res.abbr, res.symbol), repr(res))

    def testBadLookup(self):
        """ Test that lookups of non-existent residues fails """
        self.assertRaises(KeyError,
                lambda: residue.AminoAcidResidue.get('NoResidue'))

if __name__ == '__main__':
    unittest.main()
