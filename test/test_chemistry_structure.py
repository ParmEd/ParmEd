"""
Tests the chemistry/structure module
"""
from __future__ import division

import chemistry.structure as structure
from chemistry.topologyobjects import *
from chemistry.utils.six import integer_types
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
                struct.cmaps[-1].type = random.choice(struct.cmap_types)
        for i in range(random.randint(30, 40)):
            struct.trigonal_angles.append(TrigonalAngle(*random.sample(struct.atoms, 4)))
            if parametrized:
                struct.trigonal_angles[-1].type = random.choice(struct.trigonal_angle_types)
        for i in range(random.randint(30, 40)):
            struct.out_of_plane_bends.append(OutOfPlaneBend(*random.sample(struct.atoms, 4)))
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
        struct.update_dihedral_exclusions()
        return struct

    def _check_mult(self, s, s1, multfac):
        self.assertIsNot(s, s1)
        # Make sure that s is really the sum of s1 and s
        self.assertEqual(len(s.atoms), len(s1.atoms) * multfac)
        self.assertEqual(len(s.residues), len(s1.residues) * multfac)

        def cmp_atoms(a1, a2):
            self.assertIsNot(a1, a2)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.atom_type, a2.atom_type)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.radii, a2.radii)
            self.assertEqual(a1.screen, a2.screen)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.insertion_code, a2.residue.insertion_code)
            self.assertEqual(len(a1.bond_partners), len(a2.bond_partners))
            self.assertEqual(len(a1.angle_partners), len(a2.angle_partners))
            self.assertEqual(len(a1.dihedral_partners), len(a2.dihedral_partners))
            self.assertEqual(len(a1.bonds), len(a2.bonds))
            self.assertEqual(len(a1.angles), len(a2.angles))
            self.assertEqual(len(a1.dihedrals), len(a2.dihedrals))
            self.assertEqual(len(a1.impropers), len(a2.impropers))

        for a1, a2 in zip(s.atoms, s1.atoms * multfac):
            cmp_atoms(a1, a2)
        for r1, r2 in zip(s.residues, s1.residues * multfac):
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.chain, r2.chain)
            self.assertEqual(r1.insertion_code, r2.insertion_code)
        self.assertEqual(len(s.bonds), len(s1.bonds)*multfac)
        self.assertEqual(len(s.angles), len(s1.angles)*multfac)
        self.assertEqual(len(s.dihedrals), len(s1.dihedrals)*multfac)
        self.assertEqual(len(s.urey_bradleys), len(s1.urey_bradleys)*multfac)
        self.assertEqual(len(s.impropers), len(s1.impropers)*multfac)
        self.assertEqual(len(s.rb_torsions), len(s1.rb_torsions)*multfac)
        self.assertEqual(len(s.cmaps), len(s1.cmaps)*multfac)
        self.assertEqual(len(s.stretch_bends), len(s1.stretch_bends)*multfac)
        self.assertEqual(len(s.trigonal_angles), len(s1.trigonal_angles)*multfac)
        self.assertEqual(len(s.out_of_plane_bends), len(s1.out_of_plane_bends)*multfac)
        self.assertEqual(len(s.torsion_torsions), len(s1.torsion_torsions)*multfac)
        self.assertEqual(len(s.acceptors), len(s1.acceptors)*multfac)
        self.assertEqual(len(s.donors), len(s1.donors)*multfac)
        self.assertEqual(len(s.pi_torsions), len(s1.pi_torsions)*multfac)
        self.assertEqual(len(s.chiral_frames), len(s1.chiral_frames)*multfac)
        self.assertEqual(len(s.multipole_frames), len(s1.multipole_frames)*multfac)
        self.assertEqual(len(s.groups), len(s1.groups)*multfac)
        self.assertEqual(len(s.adjusts), len(s1.adjusts)*multfac)
        # Check types
        self.assertEqual(len(s.bond_types), len(s1.bond_types))
        self.assertEqual(len(s.angle_types), len(s1.angle_types))
        self.assertEqual(len(s.dihedral_types), len(s1.dihedral_types))
        self.assertEqual(len(s.urey_bradley_types), len(s1.urey_bradley_types))
        self.assertEqual(len(s.improper_types), len(s1.improper_types))
        self.assertEqual(len(s.rb_torsion_types), len(s1.rb_torsion_types))
        self.assertEqual(len(s.cmap_types), len(s1.cmap_types))
        self.assertEqual(len(s.stretch_bend_types), len(s1.stretch_bend_types))
        self.assertEqual(len(s.trigonal_angle_types), len(s1.trigonal_angle_types))
        self.assertEqual(len(s.out_of_plane_bend_types), len(s1.out_of_plane_bend_types))
        self.assertEqual(len(s.torsion_torsion_types), len(s1.torsion_torsion_types))
        self.assertEqual(len(s.pi_torsion_types), len(s1.pi_torsion_types))
        self.assertEqual(len(s.adjust_types), len(s1.adjust_types))
        # Check all valence terms
        def chk_valence(val1, val2):
            self.assertIs(type(val1[0]), type(val2[0]))
            self.assertEqual(len(val1), len(val2))
            attrs = [attr for attr in dir(val1[0]) if attr.startswith('atom')]
            for v1, v2 in zip(val1, val2):
                at1 = [getattr(v1, attr) for attr in attrs]
                at2 = [getattr(v2, attr) for attr in attrs]
                self.assertIsNot(v1, v2)
                for a1, a2 in zip(at1, at2):
                    cmp_atoms(a1, a2)
                if hasattr(v1, 'type'):
                    self.assertEqual(v1.type, v2.type)
                    if not isinstance(v1.type, integer_types):
                        if v1.type is None:
                            self.assertIs(v2.type, None)
                        else:
                            self.assertIsNot(v1.type, v2.type)
        chk_valence(s.bonds, s1.bonds*multfac)
        chk_valence(s.angles, s1.angles*multfac)
        chk_valence(s.dihedrals, s1.dihedrals*multfac)
        chk_valence(s.rb_torsions, s1.rb_torsions*multfac)
        chk_valence(s.urey_bradleys, s1.urey_bradleys*multfac)
        chk_valence(s.impropers, s1.impropers*multfac)
        chk_valence(s.cmaps, s1.cmaps*multfac)
        chk_valence(s.trigonal_angles, s1.trigonal_angles*multfac)
        chk_valence(s.out_of_plane_bends, s1.out_of_plane_bends*multfac)
        chk_valence(s.pi_torsions, s1.pi_torsions*multfac)
        chk_valence(s.torsion_torsions, s1.torsion_torsions*multfac)
        chk_valence(s.stretch_bends, s1.stretch_bends*multfac)
        chk_valence(s.chiral_frames, s1.chiral_frames*multfac)
        chk_valence(s.multipole_frames, s1.multipole_frames*multfac)
        chk_valence(s.donors, s1.donors*multfac)
        chk_valence(s.acceptors, s1.acceptors*multfac)
        chk_valence(s.groups, s1.groups*multfac)

    def _check_sum(self, s, s1, s2):
        self.assertIsNot(s, s1)
        self.assertIsNot(s, s2)
        # Make sure that s is really the sum of s1 and s2
        self.assertEqual(len(s.atoms), len(s1.atoms) + len(s2.atoms))
        self.assertEqual(len(s.residues), len(s1.residues) + len(s2.residues))

        def cmp_atoms(a1, a2):
            self.assertIsNot(a1, a2)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.atom_type, a2.atom_type)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.radii, a2.radii)
            self.assertEqual(a1.screen, a2.screen)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.insertion_code, a2.residue.insertion_code)
            self.assertEqual(len(a1.bond_partners), len(a2.bond_partners))
            self.assertEqual(len(a1.angle_partners), len(a2.angle_partners))
            self.assertEqual(len(a1.dihedral_partners), len(a2.dihedral_partners))
            self.assertEqual(len(a1.bonds), len(a2.bonds))
            self.assertEqual(len(a1.angles), len(a2.angles))
            self.assertEqual(len(a1.dihedrals), len(a2.dihedrals))
            self.assertEqual(len(a1.impropers), len(a2.impropers))

        for a1, a2 in zip(s.atoms, s1.atoms + s2.atoms):
            cmp_atoms(a1, a2)
        for r1, r2 in zip(s.residues, s1.residues + s2.residues):
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.chain, r2.chain)
            self.assertEqual(r1.insertion_code, r2.insertion_code)
        self.assertEqual(len(s.bonds), len(s1.bonds)+len(s2.bonds))
        self.assertEqual(len(s.angles), len(s1.angles)+len(s2.angles))
        self.assertEqual(len(s.dihedrals), len(s1.dihedrals)+len(s2.dihedrals))
        self.assertEqual(len(s.urey_bradleys), len(s1.urey_bradleys)+len(s2.urey_bradleys))
        self.assertEqual(len(s.impropers), len(s1.impropers)+len(s2.impropers))
        self.assertEqual(len(s.rb_torsions), len(s1.rb_torsions)+len(s2.rb_torsions))
        self.assertEqual(len(s.cmaps), len(s1.cmaps)+len(s2.cmaps))
        self.assertEqual(len(s.stretch_bends), len(s1.stretch_bends)+len(s2.stretch_bends))
        self.assertEqual(len(s.trigonal_angles), len(s1.trigonal_angles)+len(s2.trigonal_angles))
        self.assertEqual(len(s.out_of_plane_bends), len(s1.out_of_plane_bends)+len(s2.out_of_plane_bends))
        self.assertEqual(len(s.torsion_torsions), len(s1.torsion_torsions)+len(s2.torsion_torsions))
        self.assertEqual(len(s.acceptors), len(s1.acceptors)+len(s2.acceptors))
        self.assertEqual(len(s.donors), len(s1.donors)+len(s2.donors))
        self.assertEqual(len(s.pi_torsions), len(s1.pi_torsions)+len(s2.pi_torsions))
        self.assertEqual(len(s.chiral_frames), len(s1.chiral_frames)+len(s2.chiral_frames))
        self.assertEqual(len(s.multipole_frames), len(s1.multipole_frames)+len(s2.multipole_frames))
        self.assertEqual(len(s.groups), len(s1.groups)+len(s2.groups))
        self.assertEqual(len(s.adjusts), len(s1.adjusts)+len(s2.adjusts))
        # Check types
        self.assertEqual(len(s.bond_types), len(s1.bond_types) +
                len(s2.bond_types))
        self.assertEqual(len(s.angle_types),
                         len(s1.angle_types)+len(s2.angle_types))
        self.assertEqual(len(s.dihedral_types),
                         len(s1.dihedral_types)+len(s2.dihedral_types))
        self.assertEqual(len(s.urey_bradley_types),
                         len(s1.urey_bradley_types)+len(s2.urey_bradley_types))
        self.assertEqual(len(s.improper_types),
                         len(s1.improper_types)+len(s2.improper_types))
        self.assertEqual(len(s.rb_torsion_types),
                         len(s1.rb_torsion_types)+len(s2.rb_torsion_types))
        self.assertEqual(len(s.cmap_types),
                         len(s1.cmap_types)+len(s2.cmap_types))
        self.assertEqual(len(s.stretch_bend_types),
                         len(s1.stretch_bend_types)+len(s2.stretch_bend_types))
        self.assertEqual(len(s.trigonal_angle_types),
                         len(s1.trigonal_angle_types)+len(s2.trigonal_angle_types))
        self.assertEqual(len(s.out_of_plane_bend_types),
                         len(s1.out_of_plane_bend_types)+len(s2.out_of_plane_bend_types))
        self.assertEqual(len(s.torsion_torsion_types),
                         len(s1.torsion_torsion_types)+len(s2.torsion_torsion_types))
        self.assertEqual(len(s.pi_torsion_types),
                         len(s1.pi_torsion_types)+len(s2.pi_torsion_types))
        self.assertEqual(len(s.adjust_types), len(s1.adjust_types)+len(s2.adjust_types))
        # Check all valence terms
        def chk_valence(val1, val2):
            self.assertIs(type(val1[0]), type(val2[0]))
            self.assertEqual(len(val1), len(val2))
            attrs = [attr for attr in dir(val1[0]) if attr.startswith('atom')]
            for v1, v2 in zip(val1, val2):
                at1 = [getattr(v1, attr) for attr in attrs]
                at2 = [getattr(v2, attr) for attr in attrs]
                self.assertIsNot(v1, v2)
                for a1, a2 in zip(at1, at2):
                    cmp_atoms(a1, a2)
                if hasattr(v1, 'type'):
                    self.assertEqual(v1.type, v2.type)
                    if not isinstance(v1.type, integer_types):
                        self.assertIsNot(v1.type, v2.type)
        chk_valence(s.bonds, s1.bonds+s2.bonds)
        chk_valence(s.angles, s1.angles+s2.angles)
        chk_valence(s.dihedrals, s1.dihedrals+s2.dihedrals)
        chk_valence(s.rb_torsions, s1.rb_torsions+s2.rb_torsions)
        chk_valence(s.urey_bradleys, s1.urey_bradleys+s2.urey_bradleys)
        chk_valence(s.impropers, s1.impropers+s2.impropers)
        chk_valence(s.cmaps, s1.cmaps+s2.cmaps)
        chk_valence(s.trigonal_angles, s1.trigonal_angles+s2.trigonal_angles)
        chk_valence(s.out_of_plane_bends, s1.out_of_plane_bends+s2.out_of_plane_bends)
        chk_valence(s.pi_torsions, s1.pi_torsions+s2.pi_torsions)
        chk_valence(s.torsion_torsions, s1.torsion_torsions+s2.torsion_torsions)
        chk_valence(s.stretch_bends, s1.stretch_bends+s2.stretch_bends)
        chk_valence(s.chiral_frames, s1.chiral_frames+s2.chiral_frames)
        chk_valence(s.multipole_frames, s1.multipole_frames+s2.multipole_frames)
        chk_valence(s.donors, s1.donors+s2.donors)
        chk_valence(s.acceptors, s1.acceptors+s2.acceptors)
        chk_valence(s.groups, s1.groups+s2.groups)

    def testAddParametrized(self):
        """ Tests addition of two parametrized Structure instances """
        s1 = self._createStructure(parametrized=True)
        s2 = self._createStructure(parametrized=True)
        self.assertTrue(bool(s1.bond_types))
        self.assertTrue(bool(s2.bond_types))
        s = s1 + s2
        self._check_sum(s, s1, s2)

    def testIAdd(self):
        """ Tests in-place addition of two Structure instances """
        s1 = self._createStructure(parametrized=True)
        s2 = self._createStructure(parametrized=True)
        s1cp = copy(s1)
        s1 += s2
        self._check_sum(s1, s1cp, s2)

    def testAddNotParametrized(self):
        """ Tests addition of two non-parametrized Structure instances """
        s1 = self._createStructure(parametrized=False)
        s2 = self._createStructure(parametrized=False)
        self.assertFalse(bool(s1.bond_types))
        self.assertFalse(bool(s2.bond_types))
        s = s1 + s2
        self.assertIsNot(s, s1)
        self.assertIsNot(s, s2)
        # Make sure that s is really the sum of s1 and s2
        self.assertEqual(len(s.atoms), len(s1.atoms) + len(s2.atoms))
        self.assertEqual(len(s.residues), len(s1.residues) + len(s2.residues))

        def cmp_atoms(a1, a2):
            self.assertIsNot(a1, a2)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.atom_type, a2.atom_type)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.radii, a2.radii)
            self.assertEqual(a1.screen, a2.screen)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.insertion_code, a2.residue.insertion_code)
            self.assertEqual(len(a1.bond_partners), len(a2.bond_partners))
            self.assertEqual(len(a1.angle_partners), len(a2.angle_partners))
            self.assertEqual(len(a1.dihedral_partners), len(a2.dihedral_partners))
            self.assertEqual(len(a1.bonds), len(a2.bonds))
            self.assertEqual(len(a1.angles), len(a2.angles))
            self.assertEqual(len(a1.dihedrals), len(a2.dihedrals))
            self.assertEqual(len(a1.impropers), len(a2.impropers))

        for a1, a2 in zip(s.atoms, s1.atoms + s2.atoms):
            cmp_atoms(a1, a2)
        for r1, r2 in zip(s.residues, s1.residues + s2.residues):
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.chain, r2.chain)
            self.assertEqual(r1.insertion_code, r2.insertion_code)
        self.assertEqual(len(s.bonds), len(s1.bonds)+len(s2.bonds))
        self.assertEqual(len(s.angles), len(s1.angles)+len(s2.angles))
        self.assertEqual(len(s.dihedrals), len(s1.dihedrals)+len(s2.dihedrals))
        self.assertEqual(len(s.urey_bradleys), len(s1.urey_bradleys)+len(s2.urey_bradleys))
        self.assertEqual(len(s.impropers), len(s1.impropers)+len(s2.impropers))
        self.assertEqual(len(s.rb_torsions), len(s1.rb_torsions)+len(s2.rb_torsions))
        self.assertEqual(len(s.cmaps), len(s1.cmaps)+len(s2.cmaps))
        self.assertEqual(len(s.stretch_bends), len(s1.stretch_bends)+len(s2.stretch_bends))
        self.assertEqual(len(s.trigonal_angles), len(s1.trigonal_angles)+len(s2.trigonal_angles))
        self.assertEqual(len(s.out_of_plane_bends), len(s1.out_of_plane_bends)+len(s2.out_of_plane_bends))
        self.assertEqual(len(s.torsion_torsions), len(s1.torsion_torsions)+len(s2.torsion_torsions))
        self.assertEqual(len(s.acceptors), len(s1.acceptors)+len(s2.acceptors))
        self.assertEqual(len(s.donors), len(s1.donors)+len(s2.donors))
        self.assertEqual(len(s.pi_torsions), len(s1.pi_torsions)+len(s2.pi_torsions))
        self.assertEqual(len(s.chiral_frames), len(s1.chiral_frames)+len(s2.chiral_frames))
        self.assertEqual(len(s.multipole_frames), len(s1.multipole_frames)+len(s2.multipole_frames))
        self.assertEqual(len(s.groups), len(s1.groups)+len(s2.groups))
        self.assertEqual(len(s.adjusts), len(s1.adjusts)+len(s2.adjusts))
        self.assertEqual(len(s.adjust_types), len(s1.adjust_types)+len(s2.adjust_types))
        # Check all valence terms
        def chk_valence(val1, val2):
            self.assertIs(type(val1[0]), type(val2[0]))
            self.assertEqual(len(val1), len(val2))
            attrs = [attr for attr in dir(val1[0]) if attr.startswith('atom')]
            for v1, v2 in zip(val1, val2):
                at1 = [getattr(v1, attr) for attr in attrs]
                at2 = [getattr(v2, attr) for attr in attrs]
                self.assertIsNot(v1, v2)
                for a1, a2 in zip(at1, at2):
                    cmp_atoms(a1, a2)
        chk_valence(s.bonds, s1.bonds+s2.bonds)
        chk_valence(s.angles, s1.angles+s2.angles)
        chk_valence(s.dihedrals, s1.dihedrals+s2.dihedrals)
        chk_valence(s.rb_torsions, s1.rb_torsions+s2.rb_torsions)
        chk_valence(s.urey_bradleys, s1.urey_bradleys+s2.urey_bradleys)
        chk_valence(s.impropers, s1.impropers+s2.impropers)
        chk_valence(s.cmaps, s1.cmaps+s2.cmaps)
        chk_valence(s.trigonal_angles, s1.trigonal_angles+s2.trigonal_angles)
        chk_valence(s.out_of_plane_bends, s1.out_of_plane_bends+s2.out_of_plane_bends)
        chk_valence(s.pi_torsions, s1.pi_torsions+s2.pi_torsions)
        chk_valence(s.torsion_torsions, s1.torsion_torsions+s2.torsion_torsions)
        chk_valence(s.stretch_bends, s1.stretch_bends+s2.stretch_bends)
        chk_valence(s.chiral_frames, s1.chiral_frames+s2.chiral_frames)
        chk_valence(s.multipole_frames, s1.multipole_frames+s2.multipole_frames)
        chk_valence(s.donors, s1.donors+s2.donors)
        chk_valence(s.acceptors, s1.acceptors+s2.acceptors)
        chk_valence(s.groups, s1.groups+s2.groups)

    def testAddNoValence(self):
        """ Tests addition of two minimal Structure instances """
        s1 = self._createStructure(parametrized=False, novalence=True)
        s2 = self._createStructure(parametrized=False, novalence=True)
        s = s1 + s2
        self.assertIsNot(s, s1)
        self.assertIsNot(s, s2)
        # Make sure that s is really the sum of s1 and s2
        self.assertEqual(len(s.atoms), len(s1.atoms) + len(s2.atoms))
        self.assertEqual(len(s.residues), len(s1.residues) + len(s2.residues))

        def cmp_atoms(a1, a2):
            self.assertIsNot(a1, a2)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.atom_type, a2.atom_type)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.radii, a2.radii)
            self.assertEqual(a1.screen, a2.screen)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.insertion_code, a2.residue.insertion_code)
            self.assertEqual(len(a1.bond_partners), len(a2.bond_partners))
            self.assertEqual(len(a1.angle_partners), len(a2.angle_partners))
            self.assertEqual(len(a1.dihedral_partners), len(a2.dihedral_partners))
            self.assertEqual(len(a1.bonds), len(a2.bonds))
            self.assertEqual(len(a1.angles), len(a2.angles))
            self.assertEqual(len(a1.dihedrals), len(a2.dihedrals))
            self.assertEqual(len(a1.impropers), len(a2.impropers))

        for a1, a2 in zip(s.atoms, s1.atoms + s2.atoms):
            cmp_atoms(a1, a2)
        for r1, r2 in zip(s.residues, s1.residues + s2.residues):
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.chain, r2.chain)
            self.assertEqual(r1.insertion_code, r2.insertion_code)

    def testMultiplyParametrized(self):
        """ Tests replicating a parametrized Structure instance """
        s1 = self._createStructure(parametrized=True)
        multfac = random.randint(2, 5)
        s2 = s1 * multfac
        self.assertIsNot(s1, s2)
        self.assertEqual(len(s2.atoms), len(s1.atoms) * multfac)
        self._check_mult(s2, s1, multfac)

    def testMultiplyNotParametrized(self):
        """ Tests replicating a non-parametrized Structure instance """
        s1 = self._createStructure(parametrized=False)
        multfac = random.randint(2, 5)
        s2 = s1 * multfac
        self.assertIsNot(s1, s2)
        self.assertEqual(len(s2.atoms), len(s1.atoms) * multfac)
        self._check_mult(s2, s1, multfac)

    def testMultNoValence(self):
        """ Tests addition of two minimal Structure instances """
        s1 = self._createStructure(parametrized=False, novalence=True)
        multfac = random.randint(2, 5)
        s = s1 * multfac
        self.assertIsNot(s1, s)
        self.assertEqual(len(s.atoms), len(s1.atoms) * multfac)
        self.assertIsNot(s, s1)
        # Make sure that s is really the sum of s1 and s
        self.assertEqual(len(s.atoms), len(s1.atoms) * multfac)
        self.assertEqual(len(s.residues), len(s1.residues) * multfac)

        def cmp_atoms(a1, a2):
            self.assertIsNot(a1, a2)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.atom_type, a2.atom_type)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.radii, a2.radii)
            self.assertEqual(a1.screen, a2.screen)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.insertion_code, a2.residue.insertion_code)
            self.assertEqual(len(a1.bond_partners), len(a2.bond_partners))
            self.assertEqual(len(a1.angle_partners), len(a2.angle_partners))
            self.assertEqual(len(a1.dihedral_partners), len(a2.dihedral_partners))
            self.assertEqual(len(a1.bonds), len(a2.bonds))
            self.assertEqual(len(a1.angles), len(a2.angles))
            self.assertEqual(len(a1.dihedrals), len(a2.dihedrals))
            self.assertEqual(len(a1.impropers), len(a2.impropers))

        for a1, a2 in zip(s.atoms, s1.atoms * multfac):
            cmp_atoms(a1, a2)
        for r1, r2 in zip(s.residues, s1.residues * multfac):
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.chain, r2.chain)
            self.assertEqual(r1.insertion_code, r2.insertion_code)

if __name__ == '__main__':
    unittest.main()
