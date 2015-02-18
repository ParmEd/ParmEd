"""
Tests the functionality in chemistry.residue
"""

import utils
import unittest
from chemistry import residue, Atom
from chemistry.amber import AmberParm
from chemistry.exceptions import BondError
get_fn = utils.get_fn

class TestChemistryResidueTemplate(unittest.TestCase):
    """ Tests the ResidueTemplate class """

    def setUp(self):
        self.templ = templ = residue.ResidueTemplate('ACE')
        templ.add_atom(Atom(name='HH31', type='HC'))
        templ.add_atom(Atom(name='CH3', type='CT'))
        templ.add_atom(Atom(name='HH32', type='HC'))
        templ.add_atom(Atom(name='HH33', type='HC'))
        templ.add_atom(Atom(name='C', type='C'))
        templ.add_atom(Atom(name='O', type='O'))
        templ.tail = templ.atoms[4]

    def testAddAtom(self):
        """ Tests the ResidueTemplate.add_atom function """
        templ = self.templ
        self.assertEqual(len(templ), 6)
        for i in range(6):
            self.assertIs(templ[i], templ.atoms[i])
        for x, y in zip(templ, templ.atoms):
            self.assertIs(x, y)
        for i, atom in enumerate(templ):
            self.assertIs(atom, templ.atoms[i])
        self.assertIs(templ.head, None)
        self.assertIs(templ.tail, templ[-2])
        self.assertRaises(ValueError, lambda: templ.add_atom(Atom(name='C')))

    def testAddBondsAtoms(self):
        """ Tests the ResidueTemplate.add_bond function w/ indices """
        templ = self.templ
        a1, a2, a3, a4, a5, a6 = templ.atoms
        a7 = Atom(name='Unimportant', type='black hole')
        templ.add_bond(a1, a2)
        templ.add_bond(a2, a3)
        templ.add_bond(a3, a4)
        templ.add_bond(a2, a5)
        templ.add_bond(a5, a6)
        self.assertRaises(RuntimeError, lambda: templ.add_bond(a1, a7))
        self.assertRaises(BondError, lambda: templ.add_bond(a1, a1))
        self.assertIn(a1, a2.bond_partners)
        self.assertIn(a2, a1.bond_partners)
        self.assertIn(a3, a2.bond_partners)
        self.assertIn(a2, a3.bond_partners)
        self.assertIn(a5, a2.bond_partners)
        self.assertIn(a2, a5.bond_partners)
        self.assertIn(a5, a6.bond_partners)
        self.assertIn(a6, a5.bond_partners)
        self.assertEqual(len(templ.bonds), 5)

    def testAddBondsIdx(self):
        """ Tests the ResidueTemplate.add_bond function w/ atoms """
        templ = self.templ
        a1, a2, a3, a4, a5, a6 = range(6)
        a7 = Atom(name='Unimportant', type='black hole')
        templ.add_bond(a1, a2)
        templ.add_bond(a2, a3)
        templ.add_bond(a3, a4)
        templ.add_bond(a2, a5)
        templ.add_bond(a5, a6)
        self.assertRaises(RuntimeError, lambda: templ.add_bond(a1, a7))
        self.assertRaises(BondError, lambda: templ.add_bond(a1, a1))
        a1, a2, a3, a4, a5, a6 = templ.atoms
        self.assertIn(a1, a2.bond_partners)
        self.assertIn(a2, a1.bond_partners)
        self.assertIn(a3, a2.bond_partners)
        self.assertIn(a2, a3.bond_partners)
        self.assertIn(a5, a2.bond_partners)
        self.assertIn(a2, a5.bond_partners)
        self.assertIn(a5, a6.bond_partners)
        self.assertIn(a6, a5.bond_partners)
        self.assertEqual(len(templ.bonds), 5)

    def testFromResidue(self):
        """ Tests the ResidueTemplate.from_residue function """
        # Grab this residue from an amber prmtop file
        struct = AmberParm(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        for res in struct.residues:
            self._check_arbitrary_res(struct, res)

    def _check_arbitrary_res(self, struct, res):
        orig_indices = [a.idx for a in res]
        templ = residue.ResidueTemplate.from_residue(res)
        # Make sure we didn't clobber any of the atoms in res
        for i, atom in zip(orig_indices, res.atoms):
            self.assertIs(atom.list, struct.atoms)
            self.assertEqual(atom.idx, i)
        # Make sure that we have the same number of atoms in the residue as the
        # source
        self.assertEqual(len(res), len(templ))
        for a1, a2 in zip(res, templ):
            self.assertIsInstance(a1, Atom)
            self.assertIsInstance(a2, Atom)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.xx, a2.xx)
            self.assertEqual(a1.xy, a2.xy)
            self.assertEqual(a1.xz, a2.xz)
        # Make sure we have the correct number of bonds in the residue
        bondset = set()
        for atom in res:
            for bond in atom.bonds:
                if bond.atom1 in res and bond.atom2 in res:
                    bondset.add(bond)
        self.assertGreater(len(bondset), 0)
        self.assertEqual(len(bondset), len(templ.bonds))
        # Make sure that each atom has the correct number of bonds
        for i, atom in enumerate(res):
            for bond in atom.bonds:
                try:
                    id1 = res.atoms.index(bond.atom1)
                    id2 = res.atoms.index(bond.atom2)
                except ValueError:
                    if bond.atom1 in res:
                        oatom = bond.atom2
                        idx = res.atoms.index(bond.atom1)
                    else:
                        oatom = bond.atom1
                        idx = res.atoms.index(bond.atom2)
                    if oatom.residue.idx == res.idx - 1:
                        self.assertIs(templ.head, templ[idx])
                    elif oatom.residue.idx == res.idx + 1:
                        self.assertIs(templ.tail, templ[idx])
                    elif oatom.residue.idx == res.idx:
                        self.assertTrue(False) # Should never hit
                    else:
                        # Should only happen with CYX for amber prmtop...
                        self.assertEqual(res.name, 'CYX')
                        self.assertEqual(atom.name, 'SG')
                        if bond.atom1 in res:
                            self.assertIn(templ[idx], templ.connections)
                else:
                    self.assertIn(templ[id1], templ[id2].bond_partners)
                    self.assertIn(templ[id2], templ[id1].bond_partners)

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
