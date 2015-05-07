"""
Tests the chemistry/structure module
"""
from __future__ import division

import chemistry.structure as structure
from chemistry.topologyobjects import *
from chemistry.utils.six.moves import range, zip
from copy import copy
import random
import string
import unittest
import os

class TestStructureAPI(unittest.TestCase):
    """ Tests the underlying Structure API """

    def setUp(self):
        s = self.s = structure.Structure()
        s.add_atom(Atom(), 'ALA', 1, 'A')
        s.add_atom(Atom(), 'ALA', 1, 'A')
        s.add_atom(Atom(), 'ALA', 1, 'A')
        s.add_atom(Atom(), 'ALA', 1, 'A')
        s.add_atom(Atom(), 'GLY', 2, 'A')
        s.add_atom(Atom(), 'GLY', 3, 'A')
        s.add_atom(Atom(), 'GLY', 3, 'B')
        s.add_atom(Atom(), 'GLY', 3, 'B')
        s.add_atom(Atom(), 'GLY', 3, 'B')

    def testAddAtom(self):
        """ Tests the Structure.add_atom method """
        # Check that we have the expected number of residues and atoms
        s = self.s
        self.assertEqual(len(s.atoms), 9)
        self.assertEqual(len(s.residues), 4)
        self.assertEqual(len(s.residues[0]), 4)
        self.assertEqual(len(s.residues[1]), 1)
        self.assertEqual(len(s.residues[2]), 1)
        self.assertEqual(len(s.residues[3]), 3)

        for residue in s.residues:
            self.assertEqual(residue.atoms[-1].idx - residue.atoms[0].idx + 1,
                             len(residue))

    def testAddAtomToResidue(self):
        """ Tests the Structure.add_atom_to_residue method """
        s = self.s
        res = s.residues[1]
        s.add_atom_to_residue(Atom(name='TOK'), res)
        s = self.s
        self.assertEqual(len(s.atoms), 10)
        self.assertEqual(len(s.residues), 4)
        self.assertEqual(len(s.residues[0]), 4)
        self.assertEqual(len(s.residues[1]), 2)
        self.assertEqual(len(s.residues[2]), 1)
        self.assertEqual(len(s.residues[3]), 3)

        for residue in s.residues:
            self.assertEqual(residue.atoms[-1].idx - residue.atoms[0].idx + 1,
                             len(residue))
        self.assertEqual(s.atoms[5].name, 'TOK')

        # Now try to add on to the end
        res = s.residues[-1]
        s.add_atom_to_residue(Atom(name='TOK2'), res)
        self.assertEqual(len(s.atoms), 11)
        self.assertEqual(len(s.residues), 4)
        self.assertEqual(len(s.residues[0]), 4)
        self.assertEqual(len(s.residues[1]), 2)
        self.assertEqual(len(s.residues[2]), 1)
        self.assertEqual(len(s.residues[3]), 4)

        for residue in s.residues:
            self.assertEqual(residue.atoms[-1].idx - residue.atoms[0].idx + 1,
                             len(residue))
        self.assertEqual(s.atoms[5].name, 'TOK')
        self.assertEqual(s.atoms[-1].name, 'TOK2')

class TestStructureAdd(unittest.TestCase):
    """ Tests the addition property of a System """

    def _createStructure(self, parametrized, novalence=False):
        """ Create a random Structure with random attributes

        Parameters
        ----------
        parametrized : bool
            If True, add at least two of all kinds of parameters to the
            generated random structure. If False, just fill in the atoms and
            residues and some random valence terms, but no "types"
        novalence : bool, optional
            If True, no valence terms will be added. Default is False. This is
            set to False if parametrized is True
        """
        if parametrized: novalence = False
        # Generate random atom and parameter types
        atom_types = [AtomType(''.join(random.sample(string.uppercase, 3)),
                               i, random.random()*16+1, random.randint(1, 8))
                      for i in range(random.randint(8, 20))]
        bond_types = [BondType(random.random()*2, random.random()*100)
                      for i in range(random.randint(10, 20))]
        angle_types = [AngleType(random.random()*50, random.random()*120)
                       for i in range(random.randint(10, 20))]
        dihed_types = [DihedralType(random.random()*10, random.randint(1, 6),
                                    random.choice([0, 180]))
                       for i in range(random.randint(10, 20))]
        rb_types = [RBTorsionType(*[random.random()*10 for i in range(6)])]
        imp_types = [ImproperType(random.random()*100, random.choice([0, 180]))
                     for i in range(random.randint(10, 20))]
        cmap_types = [CmapType(24, [random.random()*5 for i in range(24*24)])
                      for i in range(random.randint(5, 10))]
        oop_types = [OutOfPlaneBendType(random.random()*100)
                     for i in range(random.randint(10, 20))]
        strbnd_types = [StretchBendType(random.random()*10, random.random()*10,
                                        random.random()*2, random.random()*2,
                                        random.random()*120)
                        for i in range(random.randint(10, 20))]
        ang1, ang2 = list(range(-180,180,36)), list(range(-180,180,18))
        tortor_types = [TorsionTorsionType((10, 20), ang1[:], ang2[:],
                                [random.random()*10 for j in range(200)])
                        for i in range(random.randint(5, 10))]
        for typ in atom_types:
            typ.set_lj_params(random.random()*2, random.random()*2)

        struct = structure.Structure()
        # Add atoms in residues
        for res in range(random.randint(20, 30)):
            resname = ''.join(random.sample(string.uppercase, 3))
            resid = res + 1
            for i in range(random.randint(10, 25)):
                name = ''.join(random.sample(string.uppercase, 4))
                if parametrized:
                    typ = random.choice(atom_types)
                    type = str(typ)
                    mass = typ.mass
                    atomic_number = typ.atomic_number
                else:
                    type = ''.join(random.sample(string.uppercase, 3))
                    mass = random.random() * 16 + 1
                    atomic_number = random.randint(1, 8)
                charge = random.random() * 2 - 1
                radii = random.random() * 2
                screen = random.random() * 2
                atom = Atom(atomic_number=atomic_number, type=type,
                            charge=charge, mass=mass, radii=radii,
                            screen=screen)
                if parametrized:
                    atom.atom_type = typ
                struct.add_atom(atom, resname, resid)
        if novalence:
            return struct
        # Possibly add parameter type lists
        if parametrized:
            struct.bond_types.extend([copy(x) for x in bond_types])
            struct.bond_types.claim()
            struct.angle_types.extend([copy(x) for x in angle_types])
            struct.angle_types.claim()
            struct.dihedral_types.extend([copy(x) for x in dihed_types])
            struct.dihedral_types.claim()
            struct.rb_torsion_types.extend([copy(x) for x in rb_types])
            struct.rb_torsion_types.claim()
            struct.urey_bradley_types.extend([copy(x) for x in bond_types])
            struct.urey_bradley_types.claim()
            struct.improper_types.extend([copy(x) for x in imp_types])
            struct.improper_types.claim()
            struct.cmap_types.extend([copy(x) for x in cmap_types])
            struct.cmap_types.claim()
            struct.trigonal_angle_types.extend([copy(x) for x in angle_types])
            struct.trigonal_angle_types.claim()
            struct.out_of_plane_bend_types.extend([copy(x) for x in oop_types])
            struct.out_of_plane_bend_types.claim()
            struct.pi_torsion_types.extend([copy(x) for x in dihed_types])
            struct.pi_torsion_types.claim()
            struct.stretch_bend_types.extend([copy(x) for x in strbnd_types])
            struct.stretch_bend_types.claim()
            struct.torsion_torsion_types.extend([copy(x) for x in tortor_types])
            struct.torsion_torsion_types.claim()
            struct.adjust_types.extend([NonbondedExceptionType(0.5, 0.5, 0.6, 0.6, 0.7)
                                        for i in range(random.randint(10, 20))])
            struct.adjust_types.claim()
        # Add valence terms with optional 
        for i in range(random.randint(40, 50)):
            struct.bonds.append(Bond(*random.sample(struct.atoms, 2)))
            if parametrized:
                struct.bonds[-1].type = random.choice(struct.bond_types)
        for i in range(random.randint(35, 45)):
            struct.angles.append(Angle(*random.sample(struct.atoms, 3)))
            if parametrized:
                struct.angles[-1].type = random.choice(struct.angle_types)
        for i in range(random.randint(35, 45)):
            struct.urey_bradleys.append(UreyBradley(*random.sample(struct.atoms, 2)))
            if parametrized:
                struct.urey_bradleys[-1].type = random.choice(struct.urey_bradley_types)
        for i in range(random.randint(30, 40)):
            struct.dihedrals.append(Dihedral(*random.sample(struct.atoms, 4),
                                             improper=random.choice([True, False])))
            if parametrized:
                struct.dihedrals[-1].type = random.choice(struct.dihedral_types)
        for i in range(random.randint(30, 40)):
            struct.rb_torsions.append(Dihedral(*random.sample(struct.atoms, 4)))
            if parametrized:
                struct.rb_torsions[-1].type = random.choice(struct.rb_torsion_types)
        for i in range(random.randint(10, 20)):
            struct.impropers.append(Improper(*random.sample(struct.atoms, 4)))
            if parametrized:
                struct.impropers[-1].type = random.choice(struct.improper_types)
        for i in range(random.randint(25, 35)):
            struct.cmaps.append(Cmap(*random.sample(struct.atoms, 5)))
            if parametrized:
                struct.cmaps[-1].type = random.choice(struct.cmaps)
        for i in range(random.randint(30, 40)):
            struct.trigonal_angles.append(TrigonalAngle(*random.sample(struct.atoms, 4)))
            if parametrized:
                struct.trigonal_angles[-1].type = random.choice(struct.trigonal_angle_types)
        for i in range(random.randint(30, 40)):
            struct.out_of_plane_bends.append(TrigonalAngle(*random.sample(struct.atoms, 4)))
            if parametrized:
                struct.out_of_plane_bends[-1].type = random.choice(struct.out_of_plane_bend_types)
        for i in range(random.randint(30, 40)):
            struct.stretch_bends.append(StretchBend(*random.sample(struct.atoms, 3)))
            if parametrized:
                struct.stretch_bends[-1].type = random.choice(struct.stretch_bend_types)
        for i in range(random.randint(20, 30)):
            struct.pi_torsions.append(PiTorsion(*random.sample(struct.atoms, 6)))
            if parametrized:
                struct.pi_torsions[-1].type = random.choice(struct.pi_torsion_types)
        for i in range(random.randint(10, 20)):
            struct.torsion_torsions.append(TorsionTorsion(*random.sample(struct.atoms, 5)))
            if parametrized:
                struct.torsion_torsions[-1].type = random.choice(struct.torsion_torsion_types)
        # Now use some lesser-used features
        for i in range(random.randint(5, 10)):
            struct.acceptors.append(AcceptorDonor(*random.sample(struct.atoms, 2)))
            struct.donors.append(AcceptorDonor(*random.sample(struct.atoms, 2)))
            struct.groups.append(Group(*random.sample(range(1, 11), 3)))
            struct.chiral_frames.append(ChiralFrame(*random.sample(struct.atoms, 2),
                                                    chirality=random.choice([-1, 1])))
            struct.multipole_frames.append(MultipoleFrame(random.choice(struct.atoms),
                                                          0, 1, 2, 3))
        for i in range(random.randint(20, 30)):
            struct.adjusts.append(NonbondedException(*random.sample(struct.atoms, 2)))
            if parametrized:
                struct.adjusts[-1].type = random.choice(struct.adjust_types)
        struct.prune_empty_terms()
        struct.unchange()
        return struct

    def testAddParametrized(self):
        """ Tests addition of two parametrized Structure instances """
        s1 = self._createStructure(parametrized=True)
        s2 = self._createStructure(parametrized=True)
        s = s1 + s2

    def testAddNotParametrized(self):
        """ Tests addition of two non-parametrized Structure instances """
        s1 = self._createStructure(parametrized=False)
        s2 = self._createStructure(parametrized=False)

    def testAddNoValence(self):
        """ Tests addition of two minimal Structure instances """
        s1 = self._createStructure(parametrized=False, novalence=True)
        s2 = self._createStructure(parametrized=False, novalence=True)

    def testMultiplyParametrized(self):
        """ Tests replicating a parametrized Structure instance """
        struct = self._createStructure(parametrized=True)

    def testMultiplyNotParametrized(self):
        """ Tests replicating a non-parametrized Structure instance """
        struct = self._createStructure(parametrized=False)

    def testAddNoValence(self):
        """ Tests addition of two minimal Structure instances """
        struct = self._createStructure(parametrized=False, novalence=True)

if __name__ == '__main__':
    unittest.main()
