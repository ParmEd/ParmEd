"""
Tests the functionality in parmed.modeller
"""
from __future__ import division

from copy import copy
import numpy as np
try:
    import pandas as pd
except ImportError:
    pd = None
try:
    import networkx as nx
except ImportError:
    nx = None
import os
import parmed as pmd
from parmed import Atom, read_PDB, Structure
from parmed.amber import AmberParm, AmberOFFLibrary
from parmed.exceptions import AmberWarning, Mol2Error
from parmed.modeller import (ResidueTemplate, ResidueTemplateContainer,
                             PatchTemplate,
                             PROTEIN, SOLVENT, StandardBiomolecularResidues)
from parmed.formats import Mol2File, PDBFile
from parmed.geometry import distance2
from parmed.exceptions import MoleculeError
from parmed.utils import find_atom_pairs
from parmed.utils.six import iteritems
from parmed.utils.six.moves import zip, range, StringIO
from parmed.tools import changeRadii
import random
import sys
import unittest
import utils
import warnings
get_fn = utils.get_fn

class TestResidueTemplate(unittest.TestCase):
    """ Tests the ResidueTemplate class """

    def setUp(self):
        self.templ = templ = ResidueTemplate('ACE')
        templ.add_atom(Atom(name='HH31', type='HC', atomic_number=1))
        templ.add_atom(Atom(name='CH3', type='CT', atomic_number=6))
        templ.add_atom(Atom(name='HH32', type='HC', atomic_number=1))
        templ.add_atom(Atom(name='HH33', type='HC', atomic_number=1))
        templ.add_atom(Atom(name='C', type='C', atomic_number=6))
        templ.add_atom(Atom(name='O', type='O', atomic_number=8))
        templ.type = PROTEIN
        templ.tail = templ.atoms[4]
        repr(PROTEIN) # make sure it works, don't care what its value is
        assert str(PROTEIN) == 'PROTEIN'

    def test_chemical_formula(self):
        """ Test computation of empirical chemical formula. """
        a1, a2, a3, a4, a5, a6 = self.templ.atoms
        self.templ.add_bond(a1, a2)
        self.templ.add_bond(a2, a3)
        self.templ.add_bond(a2, a4)
        self.templ.add_bond(a2, a5)
        self.templ.add_bond(a5, a6)
        self.assertEqual(self.templ.empirical_chemical_formula, 'C2H3O')

    @unittest.skipIf(pd is None, "Cannot test without pandas")
    def test_data_frame(self):
        """ Test converting ResidueTemplate to a DataFrame """
        df = self.templ.to_dataframe()
        self.assertEqual(df.shape, (6, 20))
        self.assertAlmostEqual(df.charge.sum(), 0)
        self.assertEqual(df.atomic_number.sum(), 23)

    @unittest.skipIf(nx is None, "Cannot test without networkx")
    def test_to_networkx(self):
        """ Test converting ResidueTemplate to NetworkX graph """
        a1, a2, a3, a4, a5, a6 = self.templ.atoms
        self.templ.add_bond(a1, a2)
        self.templ.add_bond(a2, a3)
        self.templ.add_bond(a3, a4)
        self.templ.add_bond(a5, a6)
        G1 = self.templ.to_networkx()
        assert not nx.is_connected(G1)
        self.templ.add_bond(a2, a5)
        G2 = self.templ.to_networkx()
        assert nx.is_connected(G2)

    def test_add_atom(self):
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
        self.assertRaises(IndexError, lambda: templ['NOAT'])
        # Make sure we can print an atom's repr when it is in a ResidueTemplate
        repr(templ.atoms[0])

    def test_to_structure(self):
        """ Tests the ResidueTemplate.to_structure function """
        a1, a2, a3, a4, a5, a6 = self.templ.atoms
        self.templ.add_bond(a1, a2)
        self.templ.add_bond(a2, a3)
        self.templ.add_bond(a3, a4)
        self.templ.add_bond(a2, a5)
        self.templ.add_bond(a5, a6)
        struct = self.templ.to_structure()

        self.assertIsInstance(struct, Structure)
        self.assertEqual(len(struct.atoms), 6)
        self.assertEqual(len(struct.residues), 1)
        self.assertEqual(len(struct.bonds), 5)

    def test_delete_bond(self):
        """ Tests the ResidueTemplate.delete_bond function """
        templ = self.templ

        a1, a2, a3, a4, a5, a6 = self.templ.atoms
        self.templ.add_bond(a1, a2)
        self.templ.add_bond(a2, a3)
        self.templ.add_bond(a3, a4)
        self.templ.add_bond(a2, a5)
        self.templ.add_bond(a5, a6)

        bond = templ.bonds[0]
        self.assertEqual(len(templ.bonds), 5)
        templ.delete_bond(bond)
        self.assertEqual(len(templ.bonds), 4)
        # Make sure we can't delete the bond again
        self.assertRaises(ValueError, lambda: templ.delete_bond(bond))

    def test_copy(self):
        """ Tests ResidueTemplate __copy__ functionality """
        for atom in self.templ:
            atom.xx = random.random() * 20 - 10
            atom.xy = random.random() * 20 - 10
            atom.xz = random.random() * 20 - 10
        templcopy = copy(self.templ)
        self.assertIsNot(templcopy, self.templ)
        self.assertEqual(len(templcopy.atoms), len(self.templ.atoms))
        self.assertEqual(len(templcopy.bonds), len(self.templ.bonds))
        for a1, a2 in zip(templcopy.atoms, self.templ):
            self.assertIsNot(a1, a2)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.xx, a2.xx)
            self.assertEqual(a1.xy, a2.xy)
            self.assertEqual(a1.xz, a2.xz)
        for b1, b2 in zip(templcopy.bonds, self.templ.bonds):
            self.assertIsNot(b1, b2)
            self.assertIsNot(b1.atom1, b2.atom1)
            self.assertIsNot(b1.atom2, b2.atom2)
            self.assertEqual(b1.atom1.name, b2.atom1.name)
            self.assertEqual(b1.atom2.name, b2.atom2.name)
        self.assertIs(self.templ.head, None)
        self.assertIs(templcopy.head, None)
        self.assertIs(self.templ.tail, self.templ[4])
        self.assertIs(templcopy.tail, templcopy[4])
        self.assertIs(templcopy.type, self.templ.type)
        self.assertEqual(templcopy.name, self.templ.name)

        # Give ResidueTemplate a head atom and a connection
        self.templ.head = self.templ[0]
        self.templ.connections.append(self.templ[1])
        templcopy = copy(self.templ)
        self.assertIs(self.templ.head, self.templ[0])
        self.assertIs(templcopy.head, templcopy[0])
        self.assertEqual(len(self.templ.connections), 1)
        self.assertEqual(len(templcopy.connections), 1)
        for a1, a2 in zip(self.templ.connections, templcopy.connections):
            self.assertIsNot(a1, a2)
            self.assertEqual(a1.name, a2.name) # Enough to verify correctness
        self.assertIs(self.templ.first_patch, None)
        self.assertIs(self.templ.last_patch, None)
        self.assertIs(templcopy.first_patch, None)
        self.assertIs(templcopy.last_patch, None)

        # Now give ResidueTemplate patches (can be blank, doesn't matter)
        self.templ.first_patch = ResidueTemplate()
        self.templ.last_patch = ResidueTemplate()
        templcopy = copy(self.templ)
        self.assertIsNot(self.templ.first_patch, None)
        self.assertIsNot(self.templ.last_patch, None)
        self.assertIs(self.templ.first_patch, templcopy.first_patch)
        self.assertIs(self.templ.last_patch, templcopy.last_patch)

    def test_fix_charge(self):
        """ Tests charge fixing for ResidueTemplate """
        self.assertRaises(ValueError, lambda:
                ResidueTemplate().fix_charges())
        charges = [random.random()*2 - 2 for a in self.templ]
        for a, charge in zip(self.templ, charges):
            a.charge = charge
        net_charge = sum(charges)
        # The odds of 6 random numbers adding to an exact integer is miniscule
        self.assertNotEqual(net_charge, round(net_charge))
        self.assertEqual(round(net_charge), int(round(net_charge)))
        # Find what the new charges *should* be
        diff = (net_charge - round(net_charge)) / 6
        fixed_charges = [x - diff for x in charges]
        self.assertAlmostEqual(sum(fixed_charges), round(net_charge), places=10)
        # Fix the charges
        self.assertEqual(sum(a.charge for a in self.templ), net_charge)
        precision = random.randint(4, 10)
        return_value = self.templ.fix_charges(precision=precision)
        self.assertAlmostEqual(sum(a.charge for a in self.templ),
                               sum(fixed_charges), places=10)
        for a, chg in zip(self.templ, fixed_charges):
            self.assertAlmostEqual(a.charge, chg, delta=2*10**-precision)
        # Check that the return value is the residue itself
        self.assertIs(return_value, self.templ)

    def test_fix_charge_2(self):
        """ Tests charge fixing to a specific value for ResidueTemplate """
        desired_charge = random.choice(range(-10, 11))
        charges = [random.random()*2 - 2 for a in self.templ]
        for a, charge in zip(self.templ, charges):
            a.charge = charge
        net_charge = sum(charges)
        # The odds of 6 random numbers adding to my exact integer are negligible
        self.assertNotEqual(net_charge, desired_charge)
        # Find what the new charges *should* be
        diff = (net_charge - desired_charge) / 6
        fixed_charges = [x - diff for x in charges]
        self.assertAlmostEqual(sum(fixed_charges), desired_charge, places=10)
        # Fix the charges
        self.assertAlmostEqual(sum(a.charge for a in self.templ), net_charge,
                               places=10)
        precision = random.randint(4, 10)
        return_value = self.templ.fix_charges(desired_charge, precision)
        self.assertAlmostEqual(sum(a.charge for a in self.templ),
                               sum(fixed_charges), places=10)
        for a, chg in zip(self.templ, fixed_charges):
            self.assertAlmostEqual(a.charge, chg, delta=2*10**precision)
        # Check that the return value is the residue itself
        self.assertIs(return_value, self.templ)

    def test_fix_charge_container(self):
        """ Tests charge fixing for ResidueTemplateContainer """
        rescont = ResidueTemplateContainer()
        self.assertRaises(ValueError, rescont.fix_charges)
        for i in range(10):
            templcopy = copy(self.templ)
            templcopy.name = '%s%d' % (templcopy.name, i)
            sumchg = 0
            for a in templcopy.atoms:
                a.charge = random.random()*2 - 2
                sumchg += a.charge
            self.assertNotEqual(sumchg, round(sumchg))
            rescont.append(templcopy)
        self.assertEqual(len(rescont), 10)
        orig_charges = [sum(a.charge for a in r) for r in rescont]
        new_charges = [round(x) for x in orig_charges]
        for res, oc, nc in zip(rescont, orig_charges, new_charges):
            self.assertNotEqual(oc, nc)
            self.assertEqual(round(oc), nc)
            self.assertEqual(sum(a.charge for a in res), oc)
        precision = random.randint(4, 10)
        retval = rescont.fix_charges(precision=precision)
        self.assertIs(retval, rescont)
        for res, oc, nc in zip(rescont, orig_charges, new_charges):
            self.assertNotEqual(oc, nc)
            self.assertEqual(round(oc), nc)
            strchgs = ['%%.%df' % precision % atom.charge for atom in res.atoms]
            # If the string charges are equal to the rounded charge to *greater
            # than the requested precision*, then it is clearly exactly equal.
            # We go 3 orders of magnitude tighter than the printed precision to
            # make sure
            self.assertAlmostEqual(round(oc), sum(float(x) for x in strchgs),
                                   places=precision+3)
            self.assertAlmostEqual(sum(a.charge for a in res), nc, places=10)
        # Make all charges alternate between 1.0 and -1.0 for the first residue,
        # then show that fixing charges does not change them
        charges = []
        for i, atom in enumerate(rescont[0].atoms):
            atom.charge = 1**(-i)
            charges.append(1**(-i))
        rescont[0].fix_charges()
        for a, c in zip(rescont[0].atoms, charges):
            self.assertEqual(a.charge, c)

    def test_delete_atom(self):
        """ Tests the ResidueTemplate.delete_atom function """
        templ = copy(self.templ)
        a1, a2, a3, a4, a5, a6 = templ.atoms
        a7 = Atom(name='Unimportant', type='black hole')
        templ.add_bond(a1, a2)
        templ.add_bond(a2, a3)
        templ.add_bond(a3, a4)
        templ.add_bond(a2, a5)
        templ.add_bond(a5, a6)
        templ.delete_atom(a6)
        #self.assertRaises(RuntimeError, lambda: templ.delete_atom(a7))
        self.assertIn(a1, a2.bond_partners)
        self.assertIn(a2, a1.bond_partners)
        self.assertIn(a3, a2.bond_partners)
        self.assertIn(a2, a3.bond_partners)
        self.assertIn(a5, a2.bond_partners)
        self.assertIn(a2, a5.bond_partners)
        self.assertNotIn(a5, a6.bond_partners)
        self.assertNotIn(a6, a5.bond_partners)
        self.assertEqual(len(templ.bonds), 4)

    @unittest.skipIf(nx is None, "Cannot test without networkx")
    def test_patch_residue(self):
        """ Tests ResidueTemplate.patch_residue function """
        templ = self.templ
        a1, a2, a3, a4, a5, a6 = templ.atoms
        templ.add_bond(a1, a2)
        templ.add_bond(a2, a3)
        templ.add_bond(a2, a4)
        templ.add_bond(a2, a5)
        templ.add_bond(a5, a6)
        self.assertIs(templ.tail, templ[4])
        patch = PatchTemplate()
        patch.delete_atoms.append(a5.name)
        patch.delete_atoms.append(a6.name)
        residue = self.templ.apply_patch(patch)
        self.assertEqual(len(residue.atoms), 4)
        self.assertEqual(len(residue.bonds), 3)
        a1, a2, a3, a4 = residue.atoms
        self.assertIn(a1, a2.bond_partners)
        self.assertIn(a2, a1.bond_partners)
        self.assertIn(a3, a2.bond_partners)
        self.assertIn(a2, a3.bond_partners)
        self.assertIn(a4, a2.bond_partners)
        self.assertIn(a2, a4.bond_partners)
        self.assertIs(residue.tail, None)

    def test_add_bonds_atoms(self):
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
        self.assertRaises(MoleculeError, lambda: templ.add_bond(a1, a1))
        self.assertIn(a1, a2.bond_partners)
        self.assertIn(a2, a1.bond_partners)
        self.assertIn(a3, a2.bond_partners)
        self.assertIn(a2, a3.bond_partners)
        self.assertIn(a5, a2.bond_partners)
        self.assertIn(a2, a5.bond_partners)
        self.assertIn(a5, a6.bond_partners)
        self.assertIn(a6, a5.bond_partners)
        self.assertEqual(len(templ.bonds), 5)

    def test_add_bonds_idx(self):
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
        self.assertRaises(MoleculeError, lambda: templ.add_bond(a1, a1))
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

    def test_from_residue(self):
        """ Tests the ResidueTemplate.from_residue function """
        # Grab this residue from an amber prmtop file
        struct = AmberParm(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        for res in struct.residues:
            self._check_arbitrary_res(struct, res)

    def test_from_residue_noorder(self):
        """ Tests the from_residue function when residue order unknown """
        struct = AmberParm(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        residues = struct.residues[:]
        del struct.residues[:]
        for res in residues:
            self.assertEqual(res.idx, -1)
        self.assertGreater(len(residues), 0)
        # Check what should be a warning
        with self.assertWarns(UserWarning):
            ResidueTemplate.from_residue(residues[1])

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
            self.assertIn(a2, templ)
            self.assertNotIn(a1, templ)
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
        self.assertIsInstance(templ.coordinates, np.ndarray)
        self.assertEqual(templ.coordinates.shape, (len(templ), 3))
        # Check that repr doesn't fail
        repr(templ)

class TestResidueTemplateContainer(unittest.TestCase):
    """ Tests the ResidueTemplateContainer class """

    def test_from_structure(self):
        """ Tests building ResidueTemplateContainer from a Structure """
        struct = AmberParm(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        cont = ResidueTemplateContainer.from_structure(struct)
        for res, sres in zip(cont, struct.residues):
            self.assertIsInstance(res, ResidueTemplate)
            self.assertEqual(len(res), len(sres))
            for a1, a2 in zip(res, sres):
                self.assertEqual(a1.name, a2.name)
                self.assertEqual(a1.type, a2.type)
                self.assertEqual(a1.charge, a2.charge)
                self.assertEqual(a1.xx, a2.xx)
                self.assertEqual(a1.xy, a2.xy)
                self.assertEqual(a1.xz, a2.xz)
        # Check accessor
        self.assertIs(cont[0], cont[cont[0].name])

    def test_to_library(self):
        """ Tests converting a ResidueTemplateContainer to a library/dict """
        lib = ResidueTemplateContainer.from_structure(
                AmberParm(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        ).to_library()
        self.assertIsInstance(lib, dict)
        self.assertEqual(len(lib.keys()), 23)
        refset = set(["NSER", "ASP", "LYS", "ILE", "HID", "LEU", "THR", "SER",
                      "PHE", "VAL", "ALA", "GLY", "ASH", "TRP", "GLU", "CYX",
                      "PRO", "MET", "TYR", "GLN", "ASN", "ARG", "CALA"])
        self.assertEqual(set(lib.keys()), refset)

class TestResidueTemplateSaver(utils.FileIOTestCase):
    " Tests the .save method on ResidueTemplate and ResidueTemplateContainer "

    def setUp(self):
        # Create a ResidueTemplate
        utils.FileIOTestCase.setUp(self)
        self.ace = templ = ResidueTemplate('ACE')
        templ.add_atom(Atom(name='HH31', type='HC'))
        templ.add_atom(Atom(name='CH3', type='CT'))
        templ.add_atom(Atom(name='HH32', type='HC'))
        templ.add_atom(Atom(name='HH33', type='HC'))
        templ.add_atom(Atom(name='C', type='C'))
        templ.add_atom(Atom(name='O', type='O'))
        templ.add_bond(0, 1)
        templ.add_bond(1, 2)
        templ.add_bond(1, 3)
        templ.add_bond(1, 4)
        templ.add_bond(4, 5)
        templ.type = PROTEIN
        templ.tail = templ.atoms[4]
        for a in templ.atoms:
            a.xx = random.random() * 10 - 5
            a.xy = random.random() * 10 - 5
            a.xz = random.random() * 10 - 5

        self.nme = templ = ResidueTemplate('NME')
        templ.add_atom(Atom(name='N', type='N'))
        templ.add_atom(Atom(name='H', type='H'))
        templ.add_atom(Atom(name='CH3', type='CT'))
        templ.add_atom(Atom(name='H31', type='H1'))
        templ.add_atom(Atom(name='H32', type='H1'))
        templ.add_atom(Atom(name='H33', type='H1'))
        templ.add_bond(0, 1)
        templ.add_bond(0, 2)
        templ.add_bond(2, 3)
        templ.add_bond(2, 4)
        templ.add_bond(2, 5)
        templ.head = templ.atoms[0]
        for a in templ.atoms:
            a.xx = random.random() * 10 - 5
            a.xy = random.random() * 10 - 5
            a.xz = random.random() * 10 - 5

        self.container = ResidueTemplateContainer()
        self.container.append(self.ace)
        self.container.append(self.nme)

    def test_residue_template_mol2_save(self):
        """ Tests ResidueTemplate.save() method for Mol2 file """
        # Check saving mol2 files by keyword
        self.ace.save(self.get_fn('test', written=True), format='mol2')
        self.assertTrue(Mol2File.id_format(self.get_fn('test', written=True)))
        x = Mol2File.parse(self.get_fn('test', written=True))
        self._check_templates(x, self.ace, preserve_headtail=False)
        self.nme.save(self.get_fn('test', written=True), format='mol2')
        self.assertTrue(Mol2File.id_format(self.get_fn('test', written=True)))
        x = Mol2File.parse(self.get_fn('test', written=True))
        self._check_templates(x, self.nme, preserve_headtail=False)
        # Check saving mol2 files by filename extension
        self.ace.save(self.get_fn('test.mol2', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol2', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol2', written=True))
        self._check_templates(x, self.ace, preserve_headtail=False)
        self.nme.save(self.get_fn('test.mol2', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol2', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol2', written=True))
        self._check_templates(x, self.nme, preserve_headtail=False)
        # Now try zipped
        self.ace.save(self.get_fn('test.mol2.gz', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol2.gz', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol2.gz', written=True))
        self._check_templates(x, self.ace, preserve_headtail=False)
        self.nme.save(self.get_fn('test.mol2.gz', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol2.gz', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol2.gz', written=True))
        self._check_templates(x, self.nme, preserve_headtail=False)
        self.ace.save(self.get_fn('test.mol2.bz2', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol2.bz2', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol2.bz2', written=True))
        self._check_templates(x, self.ace, preserve_headtail=False)
        self.nme.save(self.get_fn('test.mol2.bz2', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol2.bz2', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol2.bz2', written=True))
        self._check_templates(x, self.nme, preserve_headtail=False)

    def test_residue_template_pdb_save(self):
        """ Tests ResidueTemplate.save() method for PDB file """
        # Check saving pdb files by keyword
        self.ace.save(self.get_fn('test', written=True), format='pdb')
        self.assertTrue(PDBFile.id_format(self.get_fn('test', written=True)))
        x = PDBFile.parse(self.get_fn('test', written=True))
        self.assertEqual(len(x.atoms), len(self.ace.atoms))
        for a1, a2 in zip(x.atoms, self.ace.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.residue.name, a2.residue.name)

        self.nme.save(self.get_fn('test.pdb', written=True))
        self.assertTrue(PDBFile.id_format(self.get_fn('test', written=True)))
        x = PDBFile.parse(self.get_fn('test.pdb', written=True))
        self.assertEqual(len(x.atoms), len(self.nme.atoms))
        for a1, a2 in zip(x.atoms, self.nme.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.residue.name, a2.residue.name)

    def test_residue_template_mol3_save(self):
        """ Tests ResidueTemplate.save() method for Mol3 file """
        # Check saving mol3 files by keyword
        self.ace.save(self.get_fn('test', written=True), format='mol3')
        self.assertTrue(Mol2File.id_format(self.get_fn('test', written=True)))
        x = Mol2File.parse(self.get_fn('test', written=True))
        self._check_templates(x, self.ace, preserve_headtail=True)
        self.nme.save(self.get_fn('test', written=True), format='mol3')
        self.assertTrue(Mol2File.id_format(self.get_fn('test', written=True)))
        x = Mol2File.parse(self.get_fn('test', written=True))
        self._check_templates(x, self.nme, preserve_headtail=True)
        # Check saving mol3 files by filename extension
        self.ace.save(self.get_fn('test.mol3', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol3', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol3', written=True))
        self._check_templates(x, self.ace, preserve_headtail=True)
        self.nme.save(self.get_fn('test.mol3', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol3', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol3', written=True))
        self._check_templates(x, self.nme, preserve_headtail=True)
        # Now try zipped
        self.ace.save(self.get_fn('test.mol3.gz', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol3.gz', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol3.gz', written=True))
        self._check_templates(x, self.ace, preserve_headtail=True)
        self.nme.save(self.get_fn('test.mol3.gz', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol3.gz', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol3.gz', written=True))
        self._check_templates(x, self.nme, preserve_headtail=True)
        self.ace.save(self.get_fn('test.mol3.bz2', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol3.bz2', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol3.bz2', written=True))
        self._check_templates(x, self.ace, preserve_headtail=True)
        self.nme.save(self.get_fn('test.mol3.bz2', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol3.bz2', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol3.bz2', written=True))
        self._check_templates(x, self.nme, preserve_headtail=True)

    def test_residue_template_off_save(self):
        """ Tests ResidueTemplate.save() method for OFF lib file """
        # Check saving OFF files by keyword
        self.ace.save(self.get_fn('test', written=True), format='offlib')
        self.assertTrue(AmberOFFLibrary.id_format(self.get_fn('test', written=True)))
        x = AmberOFFLibrary.parse(self.get_fn('test', written=True))['ACE']
        self._check_templates(x, self.ace, preserve_headtail=True)
        self.nme.save(self.get_fn('test', written=True), format='offlib')
        self.assertTrue(AmberOFFLibrary.id_format(self.get_fn('test', written=True)))
        x = AmberOFFLibrary.parse(self.get_fn('test', written=True))['NME']
        self._check_templates(x, self.nme, preserve_headtail=True)
        # Check saving OFF files by filename extension
        self.ace.save(self.get_fn('test.lib', written=True))
        self.assertTrue(AmberOFFLibrary.id_format(self.get_fn('test.lib', written=True)))
        x = AmberOFFLibrary.parse(self.get_fn('test.lib', written=True))['ACE']
        self._check_templates(x, self.ace, preserve_headtail=True)
        self.nme.save(self.get_fn('test.off', written=True))
        self.assertTrue(AmberOFFLibrary.id_format(self.get_fn('test.off', written=True)))
        x = AmberOFFLibrary.parse(self.get_fn('test.off', written=True))['NME']
        self._check_templates(x, self.nme, preserve_headtail=True)
        # Now try zipped
        self.ace.save(self.get_fn('test.lib.gz', written=True))
        self.assertTrue(AmberOFFLibrary.id_format(self.get_fn('test.lib.gz', written=True)))
        x = AmberOFFLibrary.parse(self.get_fn('test.lib.gz', written=True))['ACE']
        self._check_templates(x, self.ace, preserve_headtail=True)
        self.nme.save(self.get_fn('test.off.gz', written=True))
        self.assertTrue(AmberOFFLibrary.id_format(self.get_fn('test.off.gz', written=True)))
        x = AmberOFFLibrary.parse(self.get_fn('test.off.gz', written=True))['NME']
        self._check_templates(x, self.nme, preserve_headtail=True)
        self.ace.save(self.get_fn('test.lib.bz2', written=True))
        self.assertTrue(AmberOFFLibrary.id_format(self.get_fn('test.lib.bz2', written=True)))
        x = AmberOFFLibrary.parse(self.get_fn('test.lib.bz2', written=True))['ACE']
        self._check_templates(x, self.ace, preserve_headtail=True)
        self.nme.save(self.get_fn('test.off.bz2', written=True))
        self.assertTrue(AmberOFFLibrary.id_format(self.get_fn('test.off.bz2', written=True)))
        x = AmberOFFLibrary.parse(self.get_fn('test.off.bz2', written=True))['NME']
        self._check_templates(x, self.nme, preserve_headtail=True)

    def test_residue_template_save_bad_format(self):
        """ Tests proper exceptions for bad format types to ResidueTemplate """
        for bad_fname in ('noextension', 'bad_extension.txt', 'bad_extension.bz2'):
            with self.assertRaises(ValueError):
                self.nme.save(self.get_fn(bad_fname, written=True))
        with self.assertRaises(ValueError):
            self.nme.save(self.get_fn('test.mol2', written=True), format='SOMETHING')

    def test_residue_template_container_mol2_save(self):
        """ Tests ResidueTemplateContainer.save() method for mol2 files """
        # Check saving mol2 files by keyword
        self.container.save(self.get_fn('test', written=True), format='mol2')
        self.assertTrue(Mol2File.id_format(self.get_fn('test', written=True)))
        x = Mol2File.parse(self.get_fn('test', written=True))
        self.assertIsInstance(x, ResidueTemplateContainer)
        self._check_templates(x[0], self.ace, preserve_headtail=False)
        self._check_templates(x[1], self.nme, preserve_headtail=False)
        # Make sure it is a multi-@<MOLECULE> mol2 file (so it can't be loaded
        # as a Structure)
        self.assertRaises(Mol2Error, lambda:
                Mol2File.parse(self.get_fn('test', written=True), structure=True))
        # Check saving mol2 files by filename extension
        self.container.save(self.get_fn('test.mol2', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol2', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol2', written=True))
        self.assertIsInstance(x, ResidueTemplateContainer)
        self._check_templates(x[0], self.ace, preserve_headtail=False)
        self._check_templates(x[1], self.nme, preserve_headtail=False)
        # Now try zipped
        self.container.save(self.get_fn('test.mol2.gz', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol2.gz', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol2.gz', written=True))
        self.assertIsInstance(x, ResidueTemplateContainer)
        self._check_templates(x[0], self.ace, preserve_headtail=False)
        self._check_templates(x[1], self.nme, preserve_headtail=False)
        self.container.save(self.get_fn('test.mol2.bz2', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol2.bz2', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol2.bz2', written=True))
        self.assertIsInstance(x, ResidueTemplateContainer)
        self._check_templates(x[0], self.ace, preserve_headtail=False)
        self._check_templates(x[1], self.nme, preserve_headtail=False)

    def test_residue_template_container_mol3_save(self):
        """ Tests ResidueTemplateContainer.save() method for mol3 files """
        # Check saving mol3 files by keyword
        self.container.save(self.get_fn('test', written=True), format='mol3')
        self.assertTrue(Mol2File.id_format(self.get_fn('test', written=True)))
        x = Mol2File.parse(self.get_fn('test', written=True))
        self.assertIsInstance(x, ResidueTemplateContainer)
        self._check_templates(x[0], self.ace, preserve_headtail=True)
        self._check_templates(x[1], self.nme, preserve_headtail=True)
        # Make sure it is a multi-@<MOLECULE> mol3 file (so it can't be loaded
        # as a Structure)
        with self.assertRaises(Mol2Error):
            Mol2File.parse(self.get_fn('test', written=True), structure=True)
        # Check saving mol3 files by filename extension
        self.container.save(self.get_fn('test.mol3', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol3', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol3', written=True))
        self.assertIsInstance(x, ResidueTemplateContainer)
        self._check_templates(x[0], self.ace, preserve_headtail=True)
        self._check_templates(x[1], self.nme, preserve_headtail=True)
        # Now try zipped
        self.container.save(self.get_fn('test.mol3.gz', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol3.gz', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol3.gz', written=True))
        self.assertIsInstance(x, ResidueTemplateContainer)
        self._check_templates(x[0], self.ace, preserve_headtail=True)
        self._check_templates(x[1], self.nme, preserve_headtail=True)
        self.container.save(self.get_fn('test.mol3.bz2', written=True))
        self.assertTrue(Mol2File.id_format(self.get_fn('test.mol3.bz2', written=True)))
        x = Mol2File.parse(self.get_fn('test.mol3.bz2', written=True))
        self.assertIsInstance(x, ResidueTemplateContainer)
        self._check_templates(x[0], self.ace, preserve_headtail=True)
        self._check_templates(x[1], self.nme, preserve_headtail=True)

    def test_residue_template_container_off_save(self):
        """ Tests ResidueTemplateContainer.save() method for Amber OFF files """
        # Check saving Amber OFF files by keyword
        self.container.save(self.get_fn('test', written=True), format='offlib')
        self.assertTrue(AmberOFFLibrary.id_format(self.get_fn('test', written=True)))
        x = AmberOFFLibrary.parse(self.get_fn('test', written=True))
        self._check_templates(x['ACE'], self.ace, preserve_headtail=True)
        self._check_templates(x['NME'], self.nme, preserve_headtail=True)
        # Check saving mol2 files by filename extension
        self.container.save(self.get_fn('test.off', written=True))
        self.assertTrue(AmberOFFLibrary.id_format(self.get_fn('test.off', written=True)))
        x = AmberOFFLibrary.parse(self.get_fn('test.off', written=True))
        self._check_templates(x['ACE'], self.ace, preserve_headtail=True)
        self._check_templates(x['NME'], self.nme, preserve_headtail=True)
        # Now try zipped
        self.container.save(self.get_fn('test.lib.gz', written=True))
        self.assertTrue(AmberOFFLibrary.id_format(self.get_fn('test.lib.gz', written=True)))
        x = AmberOFFLibrary.parse(self.get_fn('test.lib.gz', written=True))
        self._check_templates(x['ACE'], self.ace, preserve_headtail=True)
        self._check_templates(x['NME'], self.nme, preserve_headtail=True)
        self.container.save(self.get_fn('test.lib.bz2', written=True))
        self.assertTrue(AmberOFFLibrary.id_format(self.get_fn('test.lib.bz2', written=True)))
        x = AmberOFFLibrary.parse(self.get_fn('test.lib.bz2', written=True))
        self._check_templates(x['ACE'], self.ace, preserve_headtail=True)
        self._check_templates(x['NME'], self.nme, preserve_headtail=True)

    def test_bad_save(self):
        """ Test error handling in ResidueTemplateContainer.save """
        self.assertRaises(ValueError, lambda:
                self.container.save('unrecognized.file.bz2'))
        self.assertRaises(ValueError, lambda:
                self.container.save('unrecognized.file'))
        self.assertRaises(ValueError, lambda:
                self.container.save('something', format='NOPE'))

    def _check_templates(self, templ1, templ2, preserve_headtail=True):
        self.assertEqual(len(templ1.atoms), len(templ2.atoms))
        self.assertEqual(len(templ1.bonds), len(templ2.bonds))
        if preserve_headtail:
            if templ1.head is None or templ2.head is None:
                self.assertIs(templ1.head, None)
                self.assertIs(templ2.head, None)
            else:
                self.assertIsNot(templ1.head, None)
                self.assertIsNot(templ2.head, None)
                self.assertEqual(templ1.head.idx, templ2.head.idx)
            if templ1.tail is None or templ2.tail is None:
                self.assertIs(templ1.tail, None)
                self.assertIs(templ2.tail, None)
            else:
                self.assertIsNot(templ1.tail, None)
                self.assertIsNot(templ2.tail, None)
                self.assertEqual(templ1.tail.idx, templ2.tail.idx)
        elif preserve_headtail is not None:
            self.assertIs(templ1.head, None)
            self.assertIs(templ1.tail, None)
        bs1 = set()
        bs2 = set()
        for b1, b2 in zip(templ1.bonds, templ2.bonds):
            bs1.add(tuple(sorted([b1.atom1.idx, b1.atom2.idx])))
            bs2.add(tuple(sorted([b2.atom1.idx, b2.atom2.idx])))
        self.assertEqual(bs1, bs2)

class TestAmberOFFLibrary(utils.FileIOTestCase):
    """ Tests the AmberOFFLibrary class """

    def test_from_library(self):
        """ Tests ResidueTemplateContainer.from_library functionality """
        offlib = AmberOFFLibrary.parse(get_fn('amino12.lib'))
        lib = ResidueTemplateContainer.from_library(offlib)
        self.assertIsInstance(lib, ResidueTemplateContainer)
        self.assertEqual(len(lib), len(offlib))
        for r1, r2 in zip(offlib, lib):
            r1 = offlib[r1]
            self.assertIs(r1, r2)
        lib2 = ResidueTemplateContainer.from_library(offlib, copy=True)
        for r1, r2 in zip(offlib, lib2):
            r1 = offlib[r1]
            self.assertIsNot(r1, r2)

    def test_read_internal(self):
        """ Tests reading Amber amino12 OFF library (internal residues) """
        offlib = AmberOFFLibrary.parse(get_fn('amino12.lib'))
        self.assertEqual(len(offlib), 28)
        for name, res in iteritems(offlib):
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

    def test_read_n_term(self):
        """ Test reading N-terminal amino acid Amber OFF library """
        offlib = AmberOFFLibrary.parse(get_fn('aminont12.lib'))
        self.assertEqual(len(offlib), 24)
        for name, res in iteritems(offlib):
            self.assertIsInstance(res, ResidueTemplate)
            self.assertEqual(name, res.name)
            self.assertIs(res.head, None)
            self.assertEqual(res.tail.name, 'C')
            self.assertIs(res.type, PROTEIN)

    def test_read_c_term(self):
        """ Test reading C-terminal amino acid Amber OFF library """
        offlib = AmberOFFLibrary.parse(get_fn('aminoct12.lib'))
        self.assertEqual(len(offlib), 26)
        for name, res in iteritems(offlib):
            self.assertIsInstance(res, ResidueTemplate)
            self.assertEqual(name, res.name)
            self.assertEqual(res.head.name, 'N')
            self.assertIs(res.type, PROTEIN)

    def test_read_solvents(self):
        """ Test reading solvent Amber OFF lib (multi-res units) """
        # Turn off warnings... the solvents.lib file is SO broken.
        offlib = AmberOFFLibrary.parse(get_fn('solvents.lib'))
        self.assertEqual(len(offlib), 24)
        for name, res in iteritems(offlib):
            self.assertEqual(res.name, name)
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
        self.assertEqual(chcl3.name, 'CHCL3BOX')
        self.assertEqual(len(chcl3), 1375)
        self.assertEqual(chcl3.box[0], 56.496)
        self.assertEqual(chcl3.box[1], 56.496)
        self.assertEqual(chcl3.box[2], 56.496)
        for res in chcl3:
            self.assertEqual(res.name, 'CL3')
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
        # We can't convert the solvents library to a ResidueTemplateContainer,
        # since it contains ResidueTemplateContainer instances
        self.assertRaises(ValueError, lambda:
                ResidueTemplateContainer.from_library(offlib))

    def test_read_write_internal(self):
        """ Tests reading/writing of Amber OFF internal AA libs """
        offlib = AmberOFFLibrary.parse(get_fn('amino12.lib'))
        outfile = StringIO()
        AmberOFFLibrary.write(offlib, outfile)
        outfile.seek(0)
        offlib2 = AmberOFFLibrary.parse(outfile)
        self._check_read_written_libs(offlib, offlib2)

    def test_read_write_c_term(self):
        """ Tests reading/writing of Amber OFF C-terminal AA libs """
        offlib = AmberOFFLibrary.parse(get_fn('aminoct12.lib'))
        outfile = StringIO()
        AmberOFFLibrary.write(offlib, outfile)
        outfile.seek(0)
        offlib2 = AmberOFFLibrary.parse(outfile)
        self._check_read_written_libs(offlib, offlib2)

    def test_read_write_n_term(self):
        """ Tests reading/writing of Amber OFF N-terminal AA libs """
        offlib = AmberOFFLibrary.parse(get_fn('aminont12.lib'))
        outfile = StringIO()
        AmberOFFLibrary.write(offlib, outfile)
        outfile.seek(0)
        offlib2 = AmberOFFLibrary.parse(outfile)
        self._check_read_written_libs(offlib, offlib2)

    def test_read_write_solvent_lib(self):
        """ Tests reading/writing of Amber OFF solvent libs """
        offlib = AmberOFFLibrary.parse(get_fn('solvents.lib'))
        outfile = StringIO()
        AmberOFFLibrary.write(offlib, outfile)
        outfile.seek(0)
        offlib2 = AmberOFFLibrary.parse(outfile)

    def test_bad_off_files(self):
        """ Tests error checking in OFF library files """
        self.assertRaises(ValueError, lambda:
                AmberOFFLibrary.parse(get_fn('trx.prmtop')))
        with open(self.get_fn('test.off', written=True), 'w') as f:
            with open(get_fn('amino12.lib'), 'r') as ff:
                for i in range(10):
                    f.write(ff.readline())
        with self.assertRaises(RuntimeError):
            AmberOFFLibrary.parse(self.get_fn('test.off', written=True))

    @unittest.skipIf(pd is None, "Cannot test without pandas")
    def test_data_frame(self):
        """ Test converting ResidueTemplate to a DataFrame """
        offlib = AmberOFFLibrary.parse(get_fn('amino12.lib'))
        df = offlib['ALA'].to_dataframe()
        self.assertEqual(df.shape, (10, 26))
        self.assertAlmostEqual(df.charge.sum(), 0)
        self.assertEqual(df.atomic_number.sum(), 38)

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

class TestAmberOFFLeapCompatibility(utils.FileIOTestCase):
    """ Tests the AmberOFFLibrary classes written in LEaP """

    def setUp(self):
        super().setUp()
        self.tleap = utils.which('tleap')
        self.cwd = os.getcwd()
        os.chdir(self._temporary_directory.name)

    def tearDown(self):
        os.chdir(self.cwd)
        super().tearDown()

    @unittest.skipIf(utils.which('tleap') is None, "Cannot test without tleap")
    def test_amber_amino_internal(self):
        """ Test that the internal AA OFF library writes work with LEaP """
        # First create the parm to test against... we are in "writes" right now
        offlib = AmberOFFLibrary.parse(get_fn('amino12.lib'))
        AmberOFFLibrary.write(offlib, 'testinternal.lib')
        f = open('tleap_orig.in', 'w')
        f.write("""\
source "%s"
l = sequence {ALA ARG ASH ASN ASP CYM CYS CYX GLH GLN GLU GLY HID HIE HIP \
              HYP ILE LEU LYN LYS MET PHE PRO SER THR TRP TYR VAL}
set default PBRadii mbondi2
savePDB l alphabet.pdb
saveAmberParm l alphabet.parm7 alphabet.rst7
quit
""" % get_fn('leaprc.ff12SB'))
        f.close()
        # Now create the leaprc for our new files
        f = open('tleap_new.in', 'w')
        f.write("""\
loadAmberParams parm10.dat
loadAmberParams frcmod.ff12SB
loadOFF testinternal.lib
l = sequence {ALA ARG ASH ASN ASP CYM CYS CYX GLH GLN GLU GLY HID HIE HIP \
              HYP ILE LEU LYN LYS MET PHE PRO SER THR TRP TYR VAL}
savePDB l alphabet2.pdb
saveAmberParm l alphabet2.parm7 alphabet2.rst7
quit
""")
        f.close()
        os.system('tleap -f tleap_orig.in > tleap_orig.out 2>&1')
        os.system('tleap -f tleap_new.in > tleap_new.out 2>&1')
        # Compare the resulting files
        pdb1 = read_PDB('alphabet.pdb')
        pdb2 = read_PDB('alphabet2.pdb')
        parm1 = AmberParm('alphabet.parm7', 'alphabet.rst7')
        parm2 = AmberParm('alphabet2.parm7', 'alphabet2.rst7')
        # Since there are some specific parts of the leaprc that affect default
        # radii, change it here intentionally
        changeRadii(parm1, 'mbondi2').execute()
        changeRadii(parm2, 'mbondi2').execute()
        self._check_corresponding_files(pdb1, pdb2, parm1, parm2)

    @unittest.skipIf(utils.which('tleap') is None, "Cannot test without tleap")
    def test_amber_amino_termini(self):
        """ Test that the terminal AA OFF library writes work with LEaP """
        offlib_nter = AmberOFFLibrary.parse(get_fn('aminont12.lib'))
        offlib_cter = AmberOFFLibrary.parse(get_fn('aminoct12.lib'))
        AmberOFFLibrary.write(offlib_nter, 'testnt.lib')
        AmberOFFLibrary.write(offlib_cter, 'testct.lib')
        # Test all pairs a random set of 10 pairs
        keys1 = [random.choice(list(offlib_nter.keys())) for i in range(10)]
        keys2 = [random.choice(list(offlib_cter.keys())) for i in range(10)]
        for key1, key2 in zip(keys1, keys2):
            f = open('tleap_orig.in', 'w')
            f.write("""\
source "%s"
l = sequence {%s %s}
savePDB l alphabet.pdb
saveAmberParm l alphabet.parm7 alphabet.rst7
quit
""" % (get_fn('leaprc.ff12SB'), key1, key2))
            f.close()
            f = open('tleap_new.in', 'w')
            f.write("""\
loadAmberParams parm10.dat
loadAmberParams frcmod.ff12SB
loadOFF testct.lib
loadOFF testnt.lib
l = sequence {%s %s}
savePDB l alphabet2.pdb
saveAmberParm l alphabet2.parm7 alphabet2.rst7
quit
""" % (key1, key2))
            f.close()
            os.system('tleap -f tleap_orig.in > tleap_orig.out 2>&1')
            os.system('tleap -f tleap_new.in > tleap_new.out 2>&1')
            # Compare the resulting files
            pdb1 = read_PDB('alphabet.pdb')
            pdb2 = read_PDB('alphabet2.pdb')
            parm1 = AmberParm('alphabet.parm7', 'alphabet.rst7')
            parm2 = AmberParm('alphabet2.parm7', 'alphabet2.rst7')
            # Since there are some specific parts of the leaprc that affect
            # default radii, change it here intentionally
            changeRadii(parm1, 'mbondi2').execute()
            changeRadii(parm2, 'mbondi2').execute()
            self._check_corresponding_files(pdb1, pdb2, parm1, parm2, False)

    def _check_corresponding_files(self, pdb1, pdb2, parm1, parm2, tree=True):
        self.assertEqual(len(pdb1.atoms), len(pdb2.atoms))
        self.assertEqual(len(parm1.atoms), len(parm2.atoms))
        self.assertEqual(len(parm1.bonds), len(parm2.bonds))
        for a1, a2 in zip(pdb1.atoms, pdb2.atoms):
            self.assertEqual(a1.name, a2.name, 'pdb1 atoms:\n{}\npdb2 atoms:\n{}\n'.format(pdb1.atoms, pdb2.atoms))
            self.assertEqual(a1.atomic_number, a2.atomic_number)
        for a1, a2 in zip(parm1.atoms, parm2.atoms):
            # Check EVERYTHING
            self.assertIsNot(a1, a2)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.nb_idx, a2.nb_idx)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.atom_type.rmin, a2.atom_type.rmin)
            self.assertEqual(a1.atom_type.epsilon, a2.atom_type.epsilon)
            self.assertEqual(a1.solvent_radius, a2.solvent_radius)
            self.assertEqual(a1.screen, a2.screen)
            # Ugh. OFF libs are inconsistent
            if tree:
                self.assertEqual(a1.tree, a2.tree)
            self.assertEqual(len(a1.bonds), len(a2.bonds))
            self.assertEqual(len(a1.angles), len(a2.angles))
            self.assertEqual(len(a1.dihedrals), len(a2.dihedrals))
            set1 = set([a.name for a in a1.bond_partners])
            set2 = set([a.name for a in a2.bond_partners])
            self.assertEqual(set1, set2)
            set1 = set([a.name for a in a1.angle_partners])
            set2 = set([a.name for a in a2.angle_partners])
            self.assertEqual(set1, set2)
            set1 = set([a.name for a in a1.dihedral_partners])
            set2 = set([a.name for a in a2.dihedral_partners])
            self.assertEqual(set1, set2)
            # Check residue properties
            self.assertEqual(a1.residue.name, a2.residue.name)

class TestSlice(unittest.TestCase):
    '''Test slicing ResidueTemplate'''

    def test_slice_from_array_like(self):
        """ Test slicing by a tuple/list"""
        residue = pmd.load_file(get_fn('aminont12.lib'))['NALA']

        names = ['CA', 'CB', 'C', 'N']
        for op in [list, tuple]:
            atomlist = residue[op(names)]
            self.assertEqual(len(names), len(atomlist))
            for atom in atomlist:
                self.assertEqual(atom.name, residue[atom.name].name)

        indices = [0, 4, 7, 3]
        for op in [list, tuple]:
            atomlist = residue[op(indices)]
            self.assertEqual(len(indices), len(atomlist))
            for atom in atomlist:
                self.assertEqual(atom.name, residue[atom.name].name)

class TestBondDetermination(utils.FileIOTestCase):
    """ Tests assigning bonds to structures """

    def test_standard_residue_database(self):
        """ Tests the availability of standard residue templates """
        for resname in ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY',
                        'HIS', 'HYP', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                        'SER', 'THR', 'TRP', 'TYR', 'VAL', 'DA', 'DT', 'DG',
                        'DC', 'A', 'U', 'G', 'C', 'ACE', 'NME'):
            self.assertIn(resname, StandardBiomolecularResidues)
            self.assertIsInstance(StandardBiomolecularResidues[resname],
                                  ResidueTemplate)

    def test_simple_bond_assignment(self):
        """ Tests the assignment of bonds to simple Structure instances """
        for name, res in iteritems(StandardBiomolecularResidues):
            s = Structure()
            for a in res.atoms:
                s.add_atom(copy(a), name, 1)
            s.assign_bonds()
            # Now make sure we have the same bond
            self.assertEqual(len(res.atoms), len(s.atoms))
            self.assertEqual(len(res.bonds), len(s.bonds))
            for sa, ra in zip(res.atoms, s.atoms):
                self.assertEqual(sa.name, ra.name)
                self.assertEqual(len(sa.bond_partners), len(ra.bond_partners))
                self.assertEqual({a.name for a in sa.bond_partners},
                                 {a.name for a in ra.bond_partners})

    def test_headtail_bond_assignment(self):
        """ Tests assignment of bonds to polymeric Structures """
        s = read_PDB(get_fn('ava.pdb'))
        # Make sure PDB files have their bond automatically determined
        self.assertGreater(len(s.bonds), 0)
        self.assertIn(s.view[0, 'C'].atoms[0], s.view[1, 'N'].atoms[0].bond_partners)
        self.assertIn(s.view[1, 'C'].atoms[0], s.view[2, 'N'].atoms[0].bond_partners)
        self.assertIn(s.view[2, 'C'].atoms[0], s.view[3, 'N'].atoms[0].bond_partners)
        self.assertIn(s.view[3, 'C'].atoms[0], s.view[4, 'N'].atoms[0].bond_partners)
        # TER card -- nobody bonded to *next* residue
        for a in s.residues[4].atoms:
            for ba in a.bond_partners:
                self.assertNotEqual(ba.residue.idx, 5)
        self.assertIn(s.view[5, 'C'].atoms[0], s.view[6, 'N'].atoms[0].bond_partners)
        self.assertIn(s.view[6, 'C'].atoms[0], s.view[7, 'N'].atoms[0].bond_partners)
        self.assertIn(s.view[7, 'C'].atoms[0], s.view[8, 'N'].atoms[0].bond_partners)
        self.assertIn(s.view[8, 'C'].atoms[0], s.view[9, 'N'].atoms[0].bond_partners)
        # Make sure we have exactly 86 bonds defined
        self.assertEqual(len(s.bonds), 86)

    def test_connect_parsing(self):
        """ Tests processing of PDB CONECT records and see that it adds bonds """
        s = read_PDB(get_fn('4lzt.pdb'))
        # Check that the CONECT record bond specs are actual bonds
        self.assertIn(s.view[5, 'SG'].atoms[0], s.view[126, 'SG'].atoms[0].bond_partners)

    def test_bond_distance_assignment(self):
        """ Tests assignment of disulfides if no CONECT records available """
        fn = self.get_fn('test.pdb', written=True)
        with open(get_fn('4lzt.pdb'), 'r') as f, open(fn, 'w') as fw:
            for line in f:
                if line.startswith('CONECT'): continue
                fw.write(line)
        s = read_PDB(fn)
        # Check that the disulfide is present even without CONECT records
        self.assertIn(s.view[5, 'SG'].atoms[0], s.view[126, 'SG'].atoms[0].bond_partners)

    def test_pairlist(self):
        """ Tests pairlist builder """
        s = Structure()
        for i in range(2000):
            s.add_atom(Atom('XYZ'), 'RES', i)
        s.coordinates = np.random.rand(2000, 3) * 20 - 10
        pairs = find_atom_pairs(s, 5.0)
        self.assertTrue(any(len(x) > 0 for x in pairs))
        for a1 in s.atoms:
            for a2 in s.atoms:
                if a1 is a2: continue
                if distance2(a1, a2) < 25:
                    self.assertIn(a1, pairs[a2.idx])
                    self.assertIn(a2, pairs[a1.idx])

    def test_element_override(self):
        """ Tests that templates improve element information """
        f = StringIO("""\
ATOM      1  CA   CA A   1      -0.641  26.272   5.586  1.00 24.68
""")
        s = read_PDB(f)
        self.assertEqual(s[0].atomic_number, pmd.periodic_table.AtomicNum['Ca'])
