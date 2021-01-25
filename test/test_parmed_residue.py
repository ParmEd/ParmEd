"""
Tests the functionality in parmed.residue
"""

import utils
import parmed as pmd
from parmed import residue
from parmed.residue import RNAResidue, DNAResidue
import unittest

class TestChemistryResidue(unittest.TestCase):

    def test_residue_members(self):
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

    def test_name_lookup(self):
        """ Test that looking up residues by name works """
        self.assertIs(residue.AminoAcidResidue.get('Alanine'), residue.ALA)
        self.assertIs(residue.AminoAcidResidue.get('Glutamine'), residue.GLN)
        self.assertIs(residue.AminoAcidResidue.get('glutamate'), residue.GLU)
        self.assertIs(residue.AminoAcidResidue.get('HiStIDiNe'), residue.HIS)
        self.assertIs(residue.AminoAcidResidue.get('ISOLEUCINE'), residue.ILE)

    def test_abbr_lookup(self):
        """ Test that looking up residues by abbreviation works """
        self.assertIs(residue.AminoAcidResidue.get('ALA'), residue.ALA)
        self.assertIs(residue.AminoAcidResidue.get('GLY'), residue.GLY)
        self.assertIs(residue.AminoAcidResidue.get('glu'), residue.GLU)
        self.assertIs(residue.AminoAcidResidue.get('Asp'), residue.ASP)

    def test_symbol_lookup(self):
        """ Test that looking up residues by symbol works """
        self.assertIs(residue.AminoAcidResidue.get('E'), residue.GLU)
        self.assertIs(residue.AminoAcidResidue.get('D'), residue.ASP)
        self.assertIs(residue.AminoAcidResidue.get('A'), residue.ALA)

    def test_repr_output(self):
        """ Test the %r representation of the Amino Acids """
        for res in residue.AminoAcidResidue.all_residues:
            self.assertEqual('<Amino Acid Residue %s: %s [%s]>' % (res.name,
                res.abbr, res.symbol), repr(res))

    def test_bad_lookup(self):
        """ Test that lookups of non-existent residues fails """
        self.assertRaises(KeyError,
                lambda: residue.AminoAcidResidue.get('NoResidue'))

    def test_has(self):
        """ Tests the `has` method of BiomolecularResidue """
        self.assertTrue(residue.AminoAcidResidue.has('E'))
        self.assertTrue(residue.AminoAcidResidue.has('GLU'))
        self.assertTrue(residue.AminoAcidResidue.has(residue.GLU))
        # Now test some that should be False
        self.assertFalse(residue.AminoAcidResidue.has(residue.A))
        self.assertFalse(residue.AminoAcidResidue.has('GUA'))
        self.assertFalse(residue.AminoAcidResidue.has('DG'))

    def test_aliases(self):
        """ Tests that the common aliases used by Amber works """
        self.assertTrue(residue.AminoAcidResidue.has('GL4'))
        self.assertTrue(residue.AminoAcidResidue.has('GLH'))
        self.assertTrue(residue.AminoAcidResidue.has('ASH'))
        self.assertIs(residue.AminoAcidResidue.get('ASH'), residue.ASP)
        self.assertIs(residue.AminoAcidResidue.get('GL4'), residue.GLU)

    def test_termini(self):
        """ Tests that decorated termini names are properly recognized """
        self.assertTrue(residue.AminoAcidResidue.has('CGLY'))
        self.assertIs(residue.AminoAcidResidue.get('CHIS'), residue.HIS)
        self.assertIs(residue.AminoAcidResidue.get('CASH'), residue.ASP)
        self.assertIs(residue.AminoAcidResidue.get('NHIS'), residue.HIS)
        self.assertIs(residue.AminoAcidResidue.get('NASH'), residue.ASP)
        with self.assertRaises(KeyError):
            residue.AminoAcidResidue.get('XASH')

class TestNucleicAcidResidues(unittest.TestCase):

    def test_residue_members(self):
        """ Tests that all of the nucleic acid residues are defined """
        self.assertTrue(hasattr(residue, 'A'))
        self.assertTrue(hasattr(residue, 'U'))
        self.assertTrue(hasattr(residue, 'G'))
        self.assertTrue(hasattr(residue, 'C'))
        self.assertTrue(hasattr(residue, 'DA'))
        self.assertTrue(hasattr(residue, 'DT'))
        self.assertTrue(hasattr(residue, 'DG'))
        self.assertTrue(hasattr(residue, 'DC'))

    def test_name_lookup(self):
        """ Test that looking up DNA/RNA residues by name works """
        self.assertIs(residue.DNAResidue.get('Guanine'), residue.DG)
        self.assertIs(residue.DNAResidue.get('ADE'), residue.DA)
        self.assertIs(residue.DNAResidue.get('thymine'), residue.DT)
        self.assertIs(residue.RNAResidue.get('Uracil'), residue.U)
        self.assertIs(residue.RNAResidue.get('G'), residue.G)
        self.assertIs(residue.RNAResidue.get('ADE'), residue.A)

    def test_bad_lookup(self):
        """ Test that lookups of non-existent nucleic acid residues fails """
        self.assertRaises(KeyError, lambda: residue.RNAResidue.get('NoResidue'))
        self.assertRaises(KeyError, lambda: residue.DNAResidue.get('NoResidue'))

    def test_termini(self):
        """ Tests that decorated DNA/RNA termini are properly recognized """
        self.assertTrue(residue.RNAResidue.has('G5'))
        self.assertTrue(residue.RNAResidue.has('G3'))
        self.assertTrue(residue.RNAResidue.has('RG5'))
        self.assertTrue(residue.RNAResidue.has('RG3'))
        self.assertTrue(residue.DNAResidue.has('DG5'))
        self.assertTrue(residue.DNAResidue.has('DG3'))
        self.assertFalse(residue.RNAResidue.has('G4'))
