"""
Tests the fancy indexing and slicing capabilities of Structure
"""
from collections import defaultdict
import parmed as pmd
from parmed.utils.six import iteritems
from parmed.utils.six.moves import range, zip
import random
import unittest
import utils

parm = pmd.load_file(utils.get_fn('trx.prmtop'))
pdb1 = pmd.load_file(utils.get_fn('4lzt.pdb'))
pdb2 = pmd.load_file(utils.get_fn('1kip.cif'))

class TestStructureSlicing(unittest.TestCase):
    """ Tests the fancy slicing/indexing of Structure """

    def testIntIndex(self):
        """ Tests simple Structure indexing (integer) """
        for i, atom in enumerate(parm.atoms):
            self.assertIs(atom, parm[i])
        for i, atom in enumerate(pdb1.atoms):
            self.assertIs(atom, pdb1[i])
        for i, atom in enumerate(pdb2.atoms):
            self.assertIs(atom, pdb2[i])

    def testTwoIntIndex(self):
        """ Tests simple Structure indexing w/ residue and atom (int, int) """
        for i, res in enumerate(parm.residues):
            for j, atom in enumerate(res):
                self.assertIs(res[j], parm[i,j])
        for i, res in enumerate(pdb1.residues):
            for j, atom in enumerate(res):
                self.assertIs(res[j], pdb1[i,j])
        for i, res in enumerate(pdb2.residues):
            for j, atom in enumerate(res):
                self.assertIs(res[j], pdb2[i,j])

    def testThreeSimpleIndex(self):
        """ Tests simple indexing with chain, residue, and atom ID """
        chains = defaultdict(pmd.TrackedList)
        for res in pdb2.residues:
            chains[res.chain].append(res)
        for chain_name, chain in iteritems(chains):
            for i, res in enumerate(chain):
                for j, atom in enumerate(res):
                    self.assertIs(atom, pdb2[chain_name, i, j])

    def testSingleAtomSlice(self):
        """ Test that non-trivial, single-atom selections give a Structure """
        sel1 = parm['@1']
        sel2 = parm[:1]
        self.assertIsInstance(sel1, pmd.Structure)
        self.assertIsInstance(sel2, pmd.Structure)
        self.assertEqual(len(sel1.atoms), 1)
        self.assertEqual(len(sel2.atoms), 1)

    def testEmptyStructure(self):
        """ Tests empty selection of a Structure (view) """
        sel1 = parm[10000000:]
        sel2 = pdb1['@NOTHING']
        sel3 = pdb2['G',:,:]
        sel4 = parm['NOTHING',:]
        self.assertFalse(sel1)
        self.assertFalse(sel2)
        self.assertFalse(sel3)
        self.assertFalse(sel4)
        self.assertIsInstance(sel1, pmd.amber.AmberParm)
        self.assertIsInstance(sel2, pmd.Structure)
        self.assertIsInstance(sel3, pmd.Structure)
        self.assertIsInstance(sel4, pmd.amber.AmberParm)

    def testSimpleSlice(self):
        """ Tests simple atom slicing """
        sl11 = parm[:10]
        sl21 = pdb1[:10]
        sl31 = pdb2[:10]
        sl12 = parm[10:20]
        sl22 = pdb1[10:20]
        sl32 = pdb2[10:20]
        sl13 = parm[5:100:5] # 19 selections -- 100 not inclusive
        sl23 = pdb1[5:100:5]
        sl33 = pdb2[5:100:5]
        # Check that the slices are all the correct # of atoms
        self.assertEqual(len(sl11.atoms), 10)
        self.assertEqual(len(sl21.atoms), 10)
        self.assertEqual(len(sl31.atoms), 10)
        self.assertEqual(len(sl12.atoms), 10)
        self.assertEqual(len(sl22.atoms), 10)
        self.assertEqual(len(sl32.atoms), 10)
        self.assertEqual(len(sl13.atoms), 19)
        self.assertEqual(len(sl23.atoms), 19)
        self.assertEqual(len(sl33.atoms), 19)

        # Check that the resulting types of the slices are correct
        self.assertIsInstance(sl11, pmd.amber.AmberParm)
        self.assertIsInstance(sl21, pmd.Structure)
        self.assertIsInstance(sl31, pmd.Structure)
        self.assertIsInstance(sl12, pmd.amber.AmberParm)
        self.assertIsInstance(sl22, pmd.Structure)
        self.assertIsInstance(sl32, pmd.Structure)
        self.assertIsInstance(sl13, pmd.amber.AmberParm)
        self.assertIsInstance(sl23, pmd.Structure)
        self.assertIsInstance(sl33, pmd.Structure)

        # Check that the atoms sliced out are correct
        for atom in sl11.atoms:
            patom = parm.atoms[atom.idx]
            self.assertEqual(patom.name, atom.name)
            self.assertEqual(patom.type, atom.type)
            self.assertEqual(patom.charge, atom.charge)
            self.assertEqual(patom.nb_idx, atom.nb_idx)
            self.assertEqual(patom.rmin, atom.rmin)
            self.assertEqual(patom.epsilon, atom.epsilon)
            self.assertEqual(patom.sigma, atom.sigma)

        for atom in sl22.atoms:
            patom = pdb1.atoms[atom.idx+10]
            self.assertEqual(patom.name, atom.name)
            self.assertEqual(patom.xx, atom.xx)
            self.assertEqual(patom.xy, atom.xy)
            self.assertEqual(patom.xz, atom.xz)

        for atom in sl33.atoms:
            patom = pdb2.atoms[atom.idx*5+5]
            self.assertEqual(patom.name, atom.name)
            self.assertEqual(patom.xx, atom.xx)
            self.assertEqual(patom.xy, atom.xy)
            self.assertEqual(patom.xz, atom.xz)

    def testMaskArray(self):
        """ Tests Structure selection using a mask array """
        mask = [a.name in ('CA', 'CB') for a in parm.atoms]
        sel = parm[mask]
        self.assertEqual(len(sel.atoms), 207)
        for atom in sel.atoms:
            self.assertIn(atom.name, ('CA', 'CB'))

    def testSelectionArray(self):
        """ Tests Structure selection using indices array """
        indices = [random.randint(0, len(pdb1.atoms)-1) for i in range(100)]
        sel = pdb1[indices]
        self.assertEqual(len(sel.atoms), len(set(indices)))
        self.assertGreater(len(sel.atoms), 0)
        for i, val in enumerate(sorted(set(indices))):
            self.assertEqual(sel.atoms[i].name, pdb1.atoms[val].name)
            self.assertEqual(sel.atoms[i].xx, pdb1.atoms[val].xx)
            self.assertEqual(sel.atoms[i].xy, pdb1.atoms[val].xy)
            self.assertEqual(sel.atoms[i].xz, pdb1.atoms[val].xz)

    def testResidueAtomSelection(self):
        """ Tests combined residue,atom slicing/selections """
        sel = pdb1[10:20,:5] # First five atoms of residues 10-19
        self.assertEqual(len(sel.atoms), 49) # 1 residue has only 4 atoms
        c = 0
        for res in pdb1.residues[10:20]:
            for i in range(min(len(res), 5)):
                self.assertEqual(sel.atoms[c].name, res[i].name)
                self.assertEqual(sel.atoms[c].xx, res[i].xx)
                self.assertEqual(sel.atoms[c].xy, res[i].xy)
                self.assertEqual(sel.atoms[c].xz, res[i].xz)
                c += 1

        sel = pdb1[[0,2,3,5], :2] # First 2 atoms of residues 0, 2, 3, and 5
        self.assertEqual(len(sel.atoms), 8)
        c = 0
        for i, res in enumerate([0, 2, 3, 5]):
            res = pdb1.residues[res]
            self.assertEqual(sel.atoms[c].name, res[0].name)
            self.assertEqual(sel.atoms[c].xx, res[0].xx)
            self.assertEqual(sel.atoms[c].xy, res[0].xy)
            self.assertEqual(sel.atoms[c].xz, res[0].xz)
            c += 1
            self.assertEqual(sel.atoms[c].name, res[1].name)
            self.assertEqual(sel.atoms[c].xx, res[1].xx)
            self.assertEqual(sel.atoms[c].xy, res[1].xy)
            self.assertEqual(sel.atoms[c].xz, res[1].xz)
            c += 1

        sel = pdb1[8,:]
        self.assertEqual(len(sel.atoms), len(pdb1.residues[8]))
        for a1, a2 in zip(sel.atoms, pdb1.residues[8]):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.xx, a2.xx)
            self.assertEqual(a1.xy, a2.xy)
            self.assertEqual(a1.xz, a2.xz)

        self.assertIs(pdb1[8,4], pdb1.residues[8][4])

        sel = pdb1[['ALA', 'GLY'], :]
        for atom in sel.atoms:
            self.assertIn(atom.residue.name, ('ALA', 'GLY'))

        nalagly = sum(r.name in ('ALA', 'GLY') for r in pdb1.residues)
        natom = sum(len(r) for r in pdb1.residues if r.name in ('ALA', 'GLY'))
        self.assertEqual(nalagly, len(sel.residues))
        self.assertEqual(natom, len(sel.atoms))

    def testChainResidueAtomSelection(self):
        """ Tests combined chain,residue,atom slicing/selections """
        sel = pdb2['A',:,:] # All of chain A
        for a in sel.atoms:
            self.assertEqual(a.residue.chain, 'A')
        # Get all chains
        chainA = pdb2['A',:,:]
        chainB = pdb2['B',:,:]
        chainC = pdb2['C',:,:]
        # Now try some different kinds of indexing
        sel = pdb2['B',0,0]
        self.assertIs(sel, pdb2.atoms[826])
        sel = pdb2[['A','B'], :10:2, 0]
        self.assertEqual(len(sel.atoms), 2*5*1)

class TestStructureViewSlicing(unittest.TestCase):
    """ Tests the fancy slicing/indexing of Structure """

    def testIntIndex(self):
        """ Tests simple Structure indexing (integer) (view) """
        for i, atom in enumerate(parm.atoms):
            self.assertIs(atom, parm.view[i])
        for i, atom in enumerate(pdb1.atoms):
            self.assertIs(atom, pdb1.view[i])
        for i, atom in enumerate(pdb2.atoms):
            self.assertIs(atom, pdb2.view[i])

    def testTwoIntIndex(self):
        """ Tests simple Structure view w/ residue and atom (int, int) """
        for i, res in enumerate(parm.residues):
            for j, atom in enumerate(res):
                self.assertIs(res[j], parm.view[i,j])
        for i, res in enumerate(pdb1.residues):
            for j, atom in enumerate(res):
                self.assertIs(res[j], pdb1.view[i,j])
        for i, res in enumerate(pdb2.residues):
            for j, atom in enumerate(res):
                self.assertIs(res[j], pdb2.view[i,j])

    def testThreeSimpleIndex(self):
        """ Tests simple indexing with chain, residue, and atom ID (view) """
        chains = defaultdict(pmd.TrackedList)
        for res in pdb2.residues:
            chains[res.chain].append(res)
        for chain_name, chain in iteritems(chains):
            for i, res in enumerate(chain):
                for j, atom in enumerate(res):
                    self.assertIs(atom, pdb2.view[chain_name, i, j])

    def testSingleAtomSlice(self):
        """ Test that non-trivial, single-atom selections give StructureView """
        sel1 = parm.view['@1']
        sel2 = parm.view[:1]
        self.assertIsInstance(sel1, pmd.structure.StructureView)
        self.assertIsInstance(sel2, pmd.structure.StructureView)
        self.assertEqual(len(sel1.atoms), 1)
        self.assertEqual(len(sel2.atoms), 1)
        self.assertIs(sel1.atoms[0], parm.atoms[0])
        self.assertIs(sel2.atoms[0], parm.atoms[0])

    def testEmptyStructure(self):
        """ Tests empty selection of a Structure (view) """
        sel1 = parm.view[10000000:]
        sel2 = pdb1.view['@NOTHING']
        sel3 = pdb2.view['G',:,:]
        sel4 = parm.view['NOTHING',:]
        self.assertFalse(sel1)
        self.assertFalse(sel2)
        self.assertFalse(sel3)
        self.assertFalse(sel4)
        self.assertIsInstance(sel1, pmd.structure.StructureView)
        self.assertIsInstance(sel2, pmd.structure.StructureView)
        self.assertIsInstance(sel3, pmd.structure.StructureView)
        self.assertIsInstance(sel4, pmd.structure.StructureView)

    def testSimpleSlice(self):
        """ Tests simple atom slicing (view) """
        sl11 = parm.view[:10]
        sl21 = pdb1.view[:10]
        sl31 = pdb2.view[:10]
        sl12 = parm.view[10:20]
        sl22 = pdb1.view[10:20]
        sl32 = pdb2.view[10:20]
        sl13 = parm.view[5:100:5] # 19 selections -- 100 not inclusive
        sl23 = pdb1.view[5:100:5]
        sl33 = pdb2.view[5:100:5]
        # Check that the slices are all the correct # of atoms
        self.assertEqual(len(sl11.atoms), 10)
        self.assertEqual(len(sl21.atoms), 10)
        self.assertEqual(len(sl31.atoms), 10)
        self.assertEqual(len(sl12.atoms), 10)
        self.assertEqual(len(sl22.atoms), 10)
        self.assertEqual(len(sl32.atoms), 10)
        self.assertEqual(len(sl13.atoms), 19)
        self.assertEqual(len(sl23.atoms), 19)
        self.assertEqual(len(sl33.atoms), 19)

        # Check that the resulting types of the slices are correct
        self.assertIsInstance(sl11, pmd.structure.StructureView)
        self.assertIsInstance(sl21, pmd.structure.StructureView)
        self.assertIsInstance(sl31, pmd.structure.StructureView)
        self.assertIsInstance(sl12, pmd.structure.StructureView)
        self.assertIsInstance(sl22, pmd.structure.StructureView)
        self.assertIsInstance(sl32, pmd.structure.StructureView)
        self.assertIsInstance(sl13, pmd.structure.StructureView)
        self.assertIsInstance(sl23, pmd.structure.StructureView)
        self.assertIsInstance(sl33, pmd.structure.StructureView)

        # Check that the atoms sliced out are correct
        for atom in sl11.atoms:
            self.assertIs(atom, parm.atoms[atom.idx])
        for atom in sl12.atoms:
            self.assertIs(atom, parm.atoms[atom.idx])
        for atom in sl13.atoms:
            self.assertIs(atom, parm.atoms[atom.idx])
        for atom in sl21.atoms:
            self.assertIs(atom, pdb1.atoms[atom.idx])
        for atom in sl22.atoms:
            self.assertIs(atom, pdb1.atoms[atom.idx])
        for atom in sl23.atoms:
            self.assertIs(atom, pdb1.atoms[atom.idx])
        for atom in sl31.atoms:
            self.assertIs(atom, pdb2.atoms[atom.idx])
        for atom in sl32.atoms:
            self.assertIs(atom, pdb2.atoms[atom.idx])
        for atom in sl33.atoms:
            self.assertIs(atom, pdb2.atoms[atom.idx])

        # Check that iterations over views go over the atoms
        for atom in sl11:
            self.assertIs(atom, parm.atoms[atom.idx])

    def testMaskArray(self):
        """ Tests Structure selection using a mask array (view) """
        mask = [a.name in ('CA', 'CB') for a in parm.atoms]
        sel = parm.view[mask]
        self.assertIsInstance(sel, pmd.structure.StructureView)
        self.assertEqual(len(sel.atoms), 207)
        for atom in sel.atoms:
            self.assertIn(atom.name, ('CA', 'CB'))

    def testSelectionArray(self):
        """ Tests Structure selection using indices array (view) """
        indices = [random.randint(0, len(pdb1.atoms)-1) for i in range(100)]
        sel = pdb1.view[indices]
        self.assertEqual(len(sel.atoms), len(set(indices)))
        self.assertGreater(len(sel.atoms), 0)
        for i, val in enumerate(sorted(set(indices))):
            self.assertIs(sel.atoms[i], pdb1.atoms[sel.atoms[i].idx])

    def testResidueAtomSelection(self):
        """ Tests combined residue,atom slicing/selections (view) """
        sel = pdb1.view[10:20,:5] # First five atoms of residues 10-19
        self.assertIsInstance(sel, pmd.structure.StructureView)
        self.assertEqual(len(sel.atoms), 49) # 1 residue has only 4 atoms
        c = 0
        for res in pdb1.residues[10:20]:
            for i in range(min(len(res), 5)):
                self.assertIs(sel.atoms[c], res[i])
                c += 1

        sel = pdb1.view[[0,2,3,5], :2] # First 2 atoms of residues 0, 2, 3, and 5
        self.assertIsInstance(sel, pmd.structure.StructureView)
        self.assertEqual(len(sel.atoms), 8)
        c = 0
        for i, res in enumerate([0, 2, 3, 5]):
            res = pdb1.residues[res]
            self.assertIs(sel.atoms[c], res[0])
            c += 1
            self.assertIs(sel.atoms[c], res[1])
            c += 1

        sel = pdb1.view[8,:]
        self.assertIsInstance(sel, pmd.structure.StructureView)
        self.assertEqual(len(sel.atoms), len(pdb1.residues[8]))
        for a1, a2 in zip(sel.atoms, pdb1.residues[8]):
            self.assertIs(a1, a2)

        self.assertIs(pdb1[8,4], pdb1.residues[8][4])

        sel = pdb1.view[['ALA', 'GLY'], :]
        self.assertIsInstance(sel, pmd.structure.StructureView)
        for atom in sel.atoms:
            self.assertIn(atom.residue.name, ('ALA', 'GLY'))
            self.assertIs(atom.list, pdb1.atoms)

        nalagly = sum(r.name in ('ALA', 'GLY') for r in pdb1.residues)
        natom = sum(len(r) for r in pdb1.residues if r.name in ('ALA', 'GLY'))
        self.assertEqual(nalagly, len(sel.residues))
        self.assertEqual(natom, len(sel.atoms))

    def testChainResidueAtomSelection(self):
        """ Tests combined chain,residue,atom slicing/selections (view) """
        sel = pdb2.view['A',:,:] # All of chain A
        self.assertIsInstance(sel, pmd.structure.StructureView)
        for a in sel.atoms:
            self.assertEqual(a.residue.chain, 'A')
            self.assertIs(a.list, pdb2.atoms)
        # Get all chains
        chainA = pdb2.view['A',:,:]
        chainB = pdb2.view['B',:,:]
        chainC = pdb2.view['C',:,:]
        # Now try some different kinds of indexing
        sel = pdb2.view['B',0,0]
        self.assertIs(sel, pdb2.atoms[826])
        sel = pdb2.view[['A','B'], :10:2, 0]
        self.assertIsInstance(sel, pmd.structure.StructureView)
        self.assertEqual(len(sel.atoms), 2*5*1)
