"""
Tests the functionality in chemistry.modeller
"""
from chemistry import Atom
from chemistry.exceptions import AmberOFFWarning
from chemistry.modeller import (ResidueTemplate, ResidueTemplateContainer,
                                PROTEIN, SOLVENT)
from chemistry.amber import AmberParm, AmberOFFLibrary
from chemistry.exceptions import BondError
import unittest
import utils
import warnings
get_fn = utils.get_fn

try:
    import cStringIO as StringIO
except ImportError:
    # Must be Python 3
    import io as StringIO
try:
    from itertools import izip as zip
except ImportError:
    pass # Must by py3

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
            self.assertIs(res.type, PROTEIN)
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
            self.assertIs(res.type, PROTEIN)

    def testReadCTerm(self):
        """ Test reading C-terminal amino acid Amber OFF library """
        offlib = AmberOFFLibrary.parse(get_fn('aminoct12.lib'))
        self.assertEqual(len(offlib), 26)
        for name, res in offlib.items():
            self.assertIsInstance(res, ResidueTemplate)
            self.assertEqual(name, res.name)
            self.assertIs(res.head.name, 'N')
            self.assertIs(res.type, PROTEIN)

    def testReadSolvents(self):
        """ Test reading solvent Amber OFF lib (multi-res units) """
        # Turn off warnings... the solvents.lib file is SO broken.
        warnings.filterwarnings('ignore', module='.', category=AmberOFFWarning)
        offlib = AmberOFFLibrary.parse(get_fn('solvents.lib'))
        self.assertEqual(len(offlib), 24)
        for name, res in offlib.items():
            if 'BOX' in name:
                self.assertIsInstance(res, ResidueTemplateContainer)
                # Make sure all residues have the same features as the first
                for r in res:
                    self.assertIs(r.type, SOLVENT)
                    for a1, a2 in zip(r, res[0]):
                        self.assertEqual(a1.name, a2.name)
                        self.assertEqual(a1.type, a2.type)
                        self.assertEqual(a1.charge, a2.charge)
                        self.assertEqual(a1.atomic_number, a2.atomic_number)
                        self.assertEqual(len(a1.bond_partners),
                                         len(a2.bond_partners))
                        set1 = set([x.name for x in a1.bond_partners])
                        set2 = set([x.name for x in a2.bond_partners])
                        self.assertEqual(set1, set2)
                        if a1 is not a2:
                            self.assertNotEqual(a1.xx, a2.xx)
                            self.assertNotEqual(a1.xy, a2.xy)
                            self.assertNotEqual(a1.xz, a2.xz)
            else:
                self.assertIs(res.type, SOLVENT)
                self.assertIsInstance(res, ResidueTemplate)
        # Check a few solvent boxes in particular
        chcl3 = offlib['CHCL3BOX']
        self.assertEqual(len(chcl3), 1375)
        self.assertEqual(chcl3.box[0], 56.496)
        self.assertEqual(chcl3.box[1], 56.496)
        self.assertEqual(chcl3.box[2], 56.496)
        self.assertAlmostEqual(chcl3.box[3], 90, places=4)
        self.assertAlmostEqual(chcl3.box[4], 90, places=4)
        self.assertAlmostEqual(chcl3.box[5], 90, places=4)
        # Check some positions (but obviously not all)
        self.assertAlmostEqual(chcl3[0][0].xx, -22.675111)
        self.assertAlmostEqual(chcl3[0][0].xy, -13.977137)
        self.assertAlmostEqual(chcl3[0][0].xz, -21.470579)
        self.assertAlmostEqual(chcl3[1][0].xx, -9.668111)
        self.assertAlmostEqual(chcl3[1][0].xy, -15.097137)
        self.assertAlmostEqual(chcl3[1][0].xz, -18.569579)

    def testReadWriteInternal(self):
        """ Tests reading/writing of Amber OFF internal AA libs """
        offlib = AmberOFFLibrary.parse(get_fn('amino12.lib'))
        outfile = StringIO.StringIO()
        AmberOFFLibrary.write(offlib, outfile)
        outfile.seek(0)
        offlib2 = AmberOFFLibrary.parse(outfile)
        self._check_read_written_libs(offlib, offlib2)

    def testReadWriteCTerm(self):
        """ Tests reading/writing of Amber OFF C-terminal AA libs """
        offlib = AmberOFFLibrary.parse(get_fn('aminoct12.lib'))
        outfile = StringIO.StringIO()
        AmberOFFLibrary.write(offlib, outfile)
        outfile.seek(0)
        offlib2 = AmberOFFLibrary.parse(outfile)
        self._check_read_written_libs(offlib, offlib2)

    def testReadWriteNTerm(self):
        """ Tests reading/writing of Amber OFF N-terminal AA libs """
        offlib = AmberOFFLibrary.parse(get_fn('aminont12.lib'))
        outfile = StringIO.StringIO()
        AmberOFFLibrary.write(offlib, outfile)
        outfile.seek(0)
        offlib2 = AmberOFFLibrary.parse(outfile)
        self._check_read_written_libs(offlib, offlib2)

    def _check_read_written_libs(self, offlib, offlib2):
        # Check that offlib and offlib2 are equivalent
        self.assertEqual(len(offlib), len(offlib2))
        self.assertEqual(offlib.keys(), offlib2.keys())
        for key in offlib.keys():
            r1 = offlib[key]
            r2 = offlib2[key]
            # Check residues
            self.assertEqual(len(r1), len(r2))
            self.assertIs(r1.type, r2.type)
            # Check head and tail
            if r1.head is None or r2.head is None:
                self.assertIs(r1.head, None)
                self.assertIs(r2.head, None)
            else:
                self.assertEqual(r1.head.name, r2.head.name)
                self.assertEqual(r1.head.type, r2.head.type)
                self.assertEqual(r1.head.idx, r2.head.idx)
            if r1.tail is None or r2.tail is None:
                self.assertIs(r1.tail, None)
                self.assertIs(r2.tail, None)
            else:
                self.assertEqual(r1.tail.name, r2.tail.name)
                self.assertEqual(r1.tail.type, r2.tail.type)
                self.assertEqual(r1.tail.idx, r2.tail.idx)
            # Check atom properties
            for a1, a2 in zip(r1, r2):
                self.assertEqual(a1.name, a2.name)
                self.assertEqual(a1.type, a2.type)
                self.assertAlmostEqual(a1.charge, a2.charge)
                self.assertAlmostEqual(a1.xx, a2.xx, places=4)
                self.assertAlmostEqual(a1.xy, a2.xy, places=4)
                self.assertAlmostEqual(a1.xz, a2.xz, places=4)
                self.assertEqual(a1.vx, a2.vx)
                self.assertEqual(a1.vy, a2.vy)
                self.assertEqual(a1.vz, a2.vz)
            # Check bonds
            self.assertEqual(len(r1.bonds), len(r2.bonds))
            for b1, b2 in zip(r1.bonds, r2.bonds):
                self.assertEqual(b1.atom1.name, b2.atom1.name)
                self.assertEqual(b1.atom2.name, b2.atom2.name)
