"""
Tests the functionality in chemistry.modeller
"""
from chemistry import Atom
from chemistry.modeller import ResidueTemplate
from chemistry.amber import AmberParm, AmberOFFLibrary
from chemistry.exceptions import BondError
import unittest
import utils
get_fn = utils.get_fn

class TestResidueTemplate(unittest.TestCase):
    """ Tests the ResidueTemplate class """

    def setUp(self):
        self.templ = templ = ResidueTemplate('ACE')
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
        templ = ResidueTemplate.from_residue(res)
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
        # Make sure that our coordinates come as a numpy array
        if utils.has_numpy():
            self.assertIsInstance(templ.coordinates, utils.numpy.ndarray)
            self.assertEqual(templ.coordinates.shape, (len(templ)*3,))

class TestAmberOFFLibrary(unittest.TestCase):
    """ Tests the AmberOFFLibrary class """

    def testReadInternal(self):
        """ Tests reading Amber amino12 OFF library (internal residues) """
        offlib = AmberOFFLibrary.parse(get_fn('amino12.lib'))
        self.assertEqual(len(offlib), 28)
        for name, res in offlib.items():
            self.assertIsInstance(res, ResidueTemplate)
            self.assertEqual(name, res.name)
            self.assertEqual(res.head.name, 'N')
            self.assertEqual(res.tail.name, 'C')
        # Check two residues in particular: ALA and CYX
        ala = offlib['ALA']
        self.assertEqual(len(ala), 10)
        self.assertEqual(len(ala.bonds), 9)
        self.assertIn(ala[0], ala[1].bond_partners)
        self.assertIn(ala[0], ala[2].bond_partners)
        self.assertIn(ala[2], ala[3].bond_partners)
        self.assertIn(ala[2], ala[4].bond_partners)
        self.assertIn(ala[2], ala[8].bond_partners)
        self.assertIn(ala[4], ala[5].bond_partners)
        self.assertIn(ala[4], ala[6].bond_partners)
        self.assertIn(ala[4], ala[7].bond_partners)
        self.assertIn(ala[8], ala[9].bond_partners)
        self.assertAlmostEqual(ala[0].xx, 3.325770)
        self.assertAlmostEqual(ala[0].xy, 1.547909)
        self.assertAlmostEqual(ala[0].xz, -1.607204E-06)
        self.assertAlmostEqual(ala[1].xx, 3.909407)
        self.assertAlmostEqual(ala[1].xy, 0.723611)
        self.assertAlmostEqual(ala[1].xz, -2.739882E-06)
        self.assertAlmostEqual(ala[2].xx, 3.970048)
        self.assertAlmostEqual(ala[2].xy, 2.845795)
        self.assertAlmostEqual(ala[2].xz, -1.311163E-07)
        self.assertAlmostEqual(ala[3].xx, 3.671663)
        self.assertAlmostEqual(ala[3].xy, 3.400129)
        self.assertAlmostEqual(ala[3].xz, -0.889820)
        self.assertAlmostEqual(ala[4].xx, 3.576965)
        self.assertAlmostEqual(ala[4].xy, 3.653838)
        self.assertAlmostEqual(ala[4].xz, 1.232143)
        self.assertAlmostEqual(ala[5].xx, 3.877484)
        self.assertAlmostEqual(ala[5].xy, 3.115795)
        self.assertAlmostEqual(ala[5].xz, 2.131197)
        self.assertAlmostEqual(ala[6].xx, 4.075059)
        self.assertAlmostEqual(ala[6].xy, 4.623017)
        self.assertAlmostEqual(ala[6].xz, 1.205786)
        self.assertAlmostEqual(ala[7].xx, 2.496995)
        self.assertAlmostEqual(ala[7].xy, 3.801075)
        self.assertAlmostEqual(ala[7].xz, 1.241379)
        self.assertAlmostEqual(ala[8].xx, 5.485541)
        self.assertAlmostEqual(ala[8].xy, 2.705207)
        self.assertAlmostEqual(ala[8].xz, -4.398755E-06)
        self.assertAlmostEqual(ala[9].xx, 6.008824)
        self.assertAlmostEqual(ala[9].xy, 1.593175)
        self.assertAlmostEqual(ala[9].xz, -8.449768E-06)
        # now cyx
        cyx = offlib['CYX']
        self.assertEqual(len(cyx), 10)
        self.assertEqual(len(cyx.bonds), 9)
        self.assertIn(cyx[0], cyx[1].bond_partners)
        self.assertIn(cyx[0], cyx[2].bond_partners)
        self.assertIn(cyx[2], cyx[3].bond_partners)
        self.assertIn(cyx[2], cyx[4].bond_partners)
        self.assertIn(cyx[2], cyx[8].bond_partners)
        self.assertIn(cyx[4], cyx[5].bond_partners)
        self.assertIn(cyx[4], cyx[6].bond_partners)
        self.assertIn(cyx[4], cyx[7].bond_partners)
        self.assertIn(cyx[8], cyx[9].bond_partners)
        self.assertAlmostEqual(cyx[0].xx, 3.325770)
        self.assertAlmostEqual(cyx[0].xy, 1.547909)
        self.assertAlmostEqual(cyx[0].xz, -1.607204E-06)
        self.assertAlmostEqual(cyx[1].xx, 3.909407)
        self.assertAlmostEqual(cyx[1].xy, 0.723611)
        self.assertAlmostEqual(cyx[1].xz, -2.739882E-06)
        self.assertAlmostEqual(cyx[2].xx, 3.970048)
        self.assertAlmostEqual(cyx[2].xy, 2.845795)
        self.assertAlmostEqual(cyx[2].xz, -1.311163E-07)
        self.assertAlmostEqual(cyx[3].xx, 3.671663)
        self.assertAlmostEqual(cyx[3].xy, 3.400129)
        self.assertAlmostEqual(cyx[3].xz, -0.889820)
        self.assertAlmostEqual(cyx[4].xx, 3.576965)
        self.assertAlmostEqual(cyx[4].xy, 3.653838)
        self.assertAlmostEqual(cyx[4].xz, 1.232143)
        self.assertAlmostEqual(cyx[5].xx, 2.496995)
        self.assertAlmostEqual(cyx[5].xy, 3.801075)
        self.assertAlmostEqual(cyx[5].xz, 1.241379)
        self.assertAlmostEqual(cyx[6].xx, 3.877484)
        self.assertAlmostEqual(cyx[6].xy, 3.115795)
        self.assertAlmostEqual(cyx[6].xz, 2.131197)
        self.assertAlmostEqual(cyx[7].xx, 4.309573)
        self.assertAlmostEqual(cyx[7].xy, 5.303523)
        self.assertAlmostEqual(cyx[7].xz, 1.366036)
        self.assertAlmostEqual(cyx[8].xx, 5.485541)
        self.assertAlmostEqual(cyx[8].xy, 2.705207)
        self.assertAlmostEqual(cyx[8].xz, -4.398755E-06)
        self.assertAlmostEqual(cyx[9].xx, 6.008824)
        self.assertAlmostEqual(cyx[9].xy, 1.593175)
        self.assertAlmostEqual(cyx[9].xz, -8.449768E-06)
        # Check connections
        self.assertEqual(len(cyx.connections), 1)
        self.assertEqual(cyx.connections[0].name, 'SG')

    def testReadNTerm(self):
        """ Test reading N-terminal amino acid Amber OFF library """
        offlib = AmberOFFLibrary.parse(get_fn('aminont12.lib'))
        self.assertEqual(len(offlib), 24)
        for name, res in offlib.items():
            self.assertIsInstance(res, ResidueTemplate)
            self.assertEqual(name, res.name)
            self.assertIs(res.head, None)
            self.assertEqual(res.tail.name, 'C')

    def testReadCTerm(self):
        """ Test reading C-terminal amino acid Amber OFF library """
        offlib = AmberOFFLibrary.parse(get_fn('aminoct12.lib'))
        self.assertEqual(len(offlib), 26)
        for name, res in offlib.items():
            self.assertIsInstance(res, ResidueTemplate)
            self.assertEqual(name, res.name)
            self.assertIs(res.head.name, 'N')

    def testReadSolvents(self):
        """ Test reading solvent Amber OFF lib (multi-res units) """
        offlib = AmberOFFLibrary.parse(get_fn('solvents.lib'))
        self.assertEqual(len(offlib), 24)
