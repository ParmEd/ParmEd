"""
Unittests for the classes and interfaces defined in parmed.topologyobjects

By Jason Swails
"""
from __future__ import division

import math
import numpy as np
from parmed.exceptions import MoleculeError, ParameterError
from parmed.amber.readparm import AmberFormat
from parmed.constants import DEG_TO_RAD
import parmed.topologyobjects as topologyobjects
from parmed.topologyobjects import _ListItem, _FourAtomTerm, _strip_units
from parmed.topologyobjects import *
from parmed.structure import Structure
import parmed.unit as u
from parmed.utils.six import PY3
from parmed.utils.six.moves import range, zip
from copy import copy
import unittest
from utils import get_fn, has_openmm
import random
import warnings

class TestTopologyObjects(unittest.TestCase):
    """
    This test set is responsible for testing the classes in
    parmed.topologyobjects
    """

    #=============================================

    def test_list_item1(self):
        """ Test the _ListItem base class in standard containers """
        my_container = list()
        objects = [_ListItem() for i in range(100)]
        # Assign it a list container
        for obj in objects:
            self.assertEqual(obj.idx, -1)
            obj.list = my_container
            my_container.append(obj)
        for i, obj in enumerate(my_container):
            self.assertEqual(obj.idx, i)
        # Delete 10 random, but unique, objects from my_container
        choices = set()
        while len(choices) < 10:
            choices.add(random.choice(range(100)))
        choices = sorted(list(choices), reverse=True)
        for idx in choices:
            my_container.pop(idx)
        # Now check that the index of those 10 items are -1
        for i in choices:
            self.assertEqual(objects[i].idx, -1)
        # Now check that the index of the remaining objects are correctly
        # updated
        for i, obj in enumerate(my_container):
            self.assertEqual(obj.idx, i)
        # Test some corner cases
        class NewList(list): pass
        objects = NewList([_ListItem() for i in range(100)])
        objects.needs_indexing = True
        for i, obj in enumerate(objects):
            obj.list = objects
            self.assertEqual(i, obj.idx)
        objects[-1].list = NewList()
        objects[-1].list.needs_indexing = True
        self.assertEqual(objects[-1].idx, -1)

    #=============================================

    def test_list_item2(self):
        """ Test the _ListItem base class in the TrackedList container """
        my_container = TrackedList()
        objects = [_ListItem() for i in range(100)]
        # Assign it a list container
        for obj in objects:
            obj.list = my_container
            my_container.append(obj)
        for i, obj in enumerate(my_container):
            self.assertEqual(obj.idx, i)
        # Delete 10 random, but unique, objects from my_container
        choices = set()
        while len(choices) < 10:
            choices.add(random.choice(range(100)))
        choices = sorted(list(choices), reverse=True)
        for idx in choices:
            my_container.pop(idx)
        # Now check that the index of those 10 items are -1
        for i in choices:
            self.assertEqual(objects[i].idx, -1)
        # Now check that the index of the remaining objects are correctly
        # updated
        for i, obj in enumerate(my_container):
            self.assertEqual(obj.idx, i)
        # Now add those 10 objects back to the end of the list, and make sure
        # that the indexing is still correct, then delete those same items using
        # "del" instead of "pop" and make sure it still works
        for i in choices:
            my_container.append(objects[i])
        objects = [obj for obj in my_container] # Re-sync ordering
        for i in choices:
            del my_container[i]
        for i in choices:
            self.assertEqual(objects[i].idx, -1)
        for i, obj in enumerate(my_container):
            self.assertEqual(obj.idx, i)
        # Now try deleting a slice
        objects = [obj for obj in my_container]
        del my_container[::2]
        self.assertEqual(len(my_container), 45)
        for i in range(90):
            if i % 2 == 0:
                self.assertEqual(objects[i].idx, -1)
            else:
                self.assertEqual(objects[i].idx, i//2)
        for i, obj in enumerate(my_container):
            self.assertEqual(obj.idx, i)

    #=============================================

    def test_four_atom_term(self):
        """ Tests the _FourAtomTerm base class """
        a1, a2, a3, a4, a5 = Atom(), Atom(), Atom(), Atom(), Atom()
        fat = _FourAtomTerm(a1, a2, a3, a4)
        # Verify ordering
        self.assertIs(fat.atom1, a1)
        self.assertIs(fat.atom2, a2)
        self.assertIs(fat.atom3, a3)
        self.assertIs(fat.atom4, a4)
        # Verify container
        self.assertIn(a1, fat)
        self.assertIn(a2, fat)
        self.assertIn(a3, fat)
        self.assertIn(a4, fat)
        self.assertNotIn(a5, fat)
        fat.delete()
        self.assertIs(fat.atom1, None)
        self.assertIs(fat.atom2, None)
        self.assertIs(fat.atom3, None)
        self.assertIs(fat.atom4, None)

    #=============================================

    def test_atom(self):
        """ Tests the Atom object """
        a1 = Atom(atomic_number=6, name='C1', type='CT', charge=-0.1,
                  mass=12.01, nb_idx=1, solvent_radius=1.8, tree='M')
        a2 = Atom(atomic_number=6, name='C2', type='CT', charge=0.1,
                  mass=12.01, nb_idx=1, solvent_radius=1.8, tree='M')
        a3 = Atom(atomic_number=6, name='C3', type='CT', charge=0.0,
                  mass=12.01, nb_idx=1, solvent_radius=1.8, tree='M')
        a4 = Atom(atomic_number=6, name='C4', type='CT', charge=-0.1,
                  mass=12.01, nb_idx=1, solvent_radius=1.8, tree='M')
        a5 = Atom(atomic_number=6, name='C2', type='CT', charge=0.1,
                  mass=12.01, nb_idx=1, solvent_radius=1.8, tree='M')
        self.assertRaises(IndexError, a1.nonbonded_exclusions)
        # Make sure the atom attributes are transferred correctly (only check
        # the first atom)
        self.assertEqual(a1.atomic_number, 6)
        self.assertEqual(a1.element, 6)
        self.assertEqual(a1.element_name, 'C')
        self.assertEqual(a1.name, 'C1')
        self.assertEqual(a1.type, 'CT')
        self.assertEqual(a1.charge, -0.1)
        self.assertEqual(a1.mass, 12.01)
        self.assertEqual(a1.tree, 'M')
        self.assertEqual(a1.nb_idx, 1)
        self.assertEqual(a1.solvent_radius, 1.8)
        self.assertEqual(a1.idx, -1)
        self.assertEqual(a1.marked, 0)
        self.assertEqual(a1.screen, 0)

        bonds = list()
        angles = list()
        torsions = list()
        tortors = list()
        bonds.append(Bond(a1, a2))
        bonds.append(Bond(a2, a3))
        bonds.append(Bond(a3, a4))
        bonds.append(Bond(a4, a5))
        angles.append(Angle(a1, a2, a3))
        angles.append(Angle(a2, a3, a4))
        angles.append(Angle(a3, a4, a5))
        torsions.append(Dihedral(a1, a2, a3, a4))
        torsions.append(Dihedral(a2, a3, a4, a5))
        tortors.append(TorsionTorsion(a1, a2, a3, a4, a5))

        self.assertIn(a1, bonds[0])
        self.assertNotIn(a1, bonds[1])
        self.assertIn(a1, angles[0])
        self.assertNotIn(a1, angles[1])
        self.assertIn(a1, torsions[0])
        self.assertNotIn(a1, torsions[1])
        self.assertIn(a1, tortors[0])

        self.assertIn(bonds[0], a1.bonds)
        self.assertIn(bonds[0], a2.bonds)
        self.assertIn(bonds[1], a2.bonds)
        self.assertIn(bonds[1], a3.bonds)

        self.assertIn(angles[0], a1.angles)
        self.assertIn(angles[0], a2.angles)
        self.assertIn(angles[1], a2.angles)
        self.assertIn(angles[0], a3.angles)
        self.assertIn(angles[1], a3.angles)
        self.assertIn(angles[2], a3.angles)

        self.assertIn(torsions[0], a1.dihedrals)
        self.assertIn(torsions[0], a2.dihedrals)
        self.assertIn(torsions[0], a3.dihedrals)
        self.assertIn(torsions[0], a4.dihedrals)
        self.assertNotIn(torsions[0], a5.dihedrals)
        self.assertNotIn(torsions[1], a1.dihedrals)
        self.assertIn(torsions[1], a2.dihedrals)
        self.assertIn(torsions[1], a3.dihedrals)
        self.assertIn(torsions[1], a4.dihedrals)
        self.assertIn(torsions[1], a5.dihedrals)

        self.assertIn(tortors[0], a1.tortors)
        self.assertIn(tortors[0], a2.tortors)
        self.assertIn(tortors[0], a3.tortors)
        self.assertIn(tortors[0], a4.tortors)
        self.assertIn(tortors[0], a5.tortors)

        self.assertIn(a1, a2.bond_partners)
        self.assertIn(a2, a1.bond_partners)
        # Those in bond partners should not be in other partners
        self.assertNotIn(a1, a2.angle_partners)
        self.assertNotIn(a1, a2.dihedral_partners)

        self.assertIn(a1, a3.angle_partners)
        self.assertIn(a3, a1.angle_partners)
        self.assertNotIn(a1, a3.dihedral_partners)
        self.assertNotIn(a1, a3.bond_partners)
        self.assertIn(a1, a4.dihedral_partners)

        # Test it in a list
        atoms = [a1, a2, a3, a4, a5]
        for i, atom in enumerate(atoms):
            atom.list = atoms
            self.assertEqual(atom.idx, i)

        # Test exception handling
        atoms = AtomList()
        atoms.extend([a1, a2, a3, a4, a5])
        self.assertEqual(a1.nonbonded_exclusions(), [1, 2, 3, 4])
        self.assertEqual(a2.nonbonded_exclusions(), [2, 3, 4])
        self.assertEqual(a3.nonbonded_exclusions(), [3, 4])
        self.assertEqual(a4.nonbonded_exclusions(), [4])
        self.assertEqual(a5.nonbonded_exclusions(), [])
        self.assertEqual(a5.nonbonded_exclusions(only_greater=False), [0, 1, 2, 3])

        # Test L-J handling
        lja1 = Atom(atomic_number=6, name='C1', type='CT', charge=-0.1,
                    mass=12.01, nb_idx=1, solvent_radius=1.8, tree='M')
        self.assertEqual(lja1.sigma, 0.0)
        lja1.atom_type = AtomType('CT', 1, 12.01, 6)
        lja1.atom_type.set_lj_params(1.0, 2.0, 1.1, 2.1)
        self.assertEqual(lja1.epsilon, 1.0)
        self.assertEqual(lja1.uepsilon, 1.0*u.kilocalories_per_mole)
        self.assertEqual(lja1.rmin, 2.0)
        self.assertEqual(lja1.urmin, 2.0 * u.angstroms)
        self.assertEqual(lja1.sigma, 2.0*2**(-1/6)*2)
        self.assertEqual(lja1.usigma, 2.0*2**(-1/6)*2 * u.angstroms)
        lja2 = Atom(atomic_number=6, name='C2', type='CT', charge=0.1,
                    mass=12.01, nb_idx=1, solvent_radius=1.8, tree='M')
        lja2.sigma = 0.1 * u.nanometers
        self.assertEqual(lja2.sigma, 1.0)
        self.assertEqual(lja2.sigma_14, 1.0)
        lja2.sigma_14 = 0.15 * u.nanometers
        self.assertEqual(lja2.sigma_14, 1.5)
        self.assertEqual(lja2.rmin_14, 1.5*2**(1/6)/2)
        self.assertAlmostEqual(lja2.rmin, 2**(1/6)/2)
        lja3 = Atom(atomic_number=6, name='C3', type='CT', charge=0.0,
                    mass=12.01, nb_idx=1, solvent_radius=1.8, tree='M')
        lja3.atom_type = AtomType('CX', 2, 12.01, 6)
        lja3.atom_type.set_lj_params(2.0, 3.0)
        self.assertEqual(lja3.sigma_14, 3*2**(-1/6)*2)
#       lja4 = Atom(atomic_number=6, name='C4', type='CT', charge=-0.1,
#                   mass=12.01, nb_idx=1, solvent_radius=1.8, tree='M')
#       lja5 = Atom(atomic_number=6, name='C2', type='CT', charge=0.1,
#                   mass=12.01, nb_idx=1, solvent_radius=1.8, tree='M')

        a1.element = 7
        self.assertEqual(a1.atomic_number, 7)
        self.assertEqual(a1.element, 7)

        # Test L-J setters/getters
        a = Atom()
        self.assertEqual(a.rmin, 0)
        self.assertEqual(a.sigma, 0)
        self.assertEqual(a.epsilon, 0)
        self.assertEqual(a.rmin_14, 0)
        self.assertEqual(a.sigma_14, 0)
        self.assertEqual(a.epsilon_14, 0)
        a.rmin = 1.0
        self.assertEqual(a.rmin, 1.0)
        self.assertEqual(a.sigma, 2**(5/6))
        self.assertEqual(a.rmin_14, 1)
        self.assertEqual(a.sigma_14, 2**(5/6))
        self.assertEqual(a.epsilon, 0)
        self.assertEqual(a.epsilon_14, 0)
        a.rmin_14 = 2.0
        self.assertEqual(a.rmin, 1.0)
        self.assertEqual(a.sigma, 2**(5/6))
        self.assertEqual(a.rmin_14, 2)
        self.assertEqual(a.sigma_14, 2*2**(5/6))
        self.assertEqual(a.epsilon, 0)
        self.assertEqual(a.epsilon_14, 0)
        a.epsilon = 0.5
        self.assertEqual(a.rmin, 1.0)
        self.assertEqual(a.sigma, 2**(5/6))
        self.assertEqual(a.rmin_14, 2)
        self.assertEqual(a.sigma_14, 2*2**(5/6))
        self.assertEqual(a.epsilon, 0.5)
        self.assertEqual(a.epsilon_14, 0.5)
        a.epsilon_14 = 0.25
        self.assertEqual(a.rmin, 1.0)
        self.assertEqual(a.sigma, 2**(5/6))
        self.assertEqual(a.rmin_14, 2)
        self.assertEqual(a.sigma_14, 2*2**(5/6))
        self.assertEqual(a.epsilon, 0.5)
        self.assertEqual(a.epsilon_14, 0.25)

        # Test behavior of charge so that it falls back to the charge on the
        # atom type if it's not set on the atom
        a = Atom(name='CA', type='CX', mass=12.01, atomic_number=6)
        self.assertEqual(a.charge, 0)
        at = AtomType('CX', 1, 12.01, 6, charge=0.5)
        a.atom_type = at
        self.assertEqual(a.charge, 0.5)
        a.charge = 0.3
        self.assertEqual(a.charge, 0.3)
        self.assertEqual(at.charge, 0.5)

        # Test that units get stripped on assignment and that unit-ed attributes have units
        a = Atom(name='CA')
        a.charge = (0.1 * u.elementary_charge).in_units_of(u.coulomb)
        self.assertAlmostEqual(a.charge, 0.1)
        self.assertAlmostEqual(a.ucharge.value_in_unit(u.elementary_charge), 0.1)
        a.rmin = 1
        self.assertEqual(a.rmin, 1)
        self.assertEqual(a.urmin, 1*u.angstroms)
        a.sigma = 1
        self.assertEqual(a.sigma, 1)
        self.assertEqual(a.usigma, 1*u.angstroms)
        a.epsilon = 1
        self.assertEqual(a.epsilon, 1)
        self.assertEqual(a.uepsilon, 1*u.kilocalories_per_mole)
        a.rmin_14 = 1.1
        self.assertAlmostEqual(a.rmin_14, 1.1)
        self.assertAlmostEqual(a.urmin_14.value_in_unit(u.angstroms), 1.1)
        a.sigma_14 = 1.1
        self.assertAlmostEqual(a.sigma_14, 1.1)
        self.assertAlmostEqual(a.usigma_14.value_in_unit(u.angstroms), 1.1)
        a.epsilon = 1.1
        self.assertAlmostEqual(a.epsilon, 1.1)
        self.assertAlmostEqual(a.uepsilon.value_in_unit(u.kilocalories_per_mole), 1.1)

    #=============================================

    def test_atom_type(self):
        """ Test the AtomType API """
        at1 = AtomType(1, None, 12.01, 6)
        at2 = AtomType(None, 1, 12.01, 6)
        self.assertEqual(at1, at2)
        self.assertEqual(at1, 1)
        self.assertEqual(at2, 1)
        at1.set_lj_params(0.5, 1.0)
        self.assertNotEqual(at1, at2)
        # Now check out strings
        at3 = AtomType("CA", 1, 12.01, 6)
        at4 = AtomType("CA", 1, 12.01, 6)
        self.assertEqual(at3, 'CA')
        self.assertEqual(at4, (1, 'CA'))
        # bond_type -- this is like OpenMM "class"
        self.assertEqual(at3.bond_type, 'CA')
        at3.bond_type = 'CN'
        self.assertEqual(at3.name, 'CA')
        self.assertEqual(at3.bond_type, 'CN')
        # Try the sigma_14 setter
        at3.sigma_14 = 1.0
        self.assertEqual(at3.sigma_14, 1.0)
        self.assertEqual(at3.rmin_14, 2**(-5/6))
        # Check charge
        at5 = AtomType('CX', 2, 12.01, 6, charge=0.5)
        self.assertEqual(at5.charge, 0.5)
        at4.charge = 0.3
        self.assertEqual(at4.charge, 0.3)

    #=============================================

    def test_extra_point(self):
        """ Test ExtraPoint API and exclusion handling """
        # Analyze the exception parameters for bonding pattern
        #
        # E1 -- A1 -- A2 -- A3 -- A4 -- A5 -- E5
        #             |     |     |
        #             E2    E3    E4

        ep1 = ExtraPoint(name='E1', type='EP', atomic_number=0, weights=[1, 2])
        ep2 = ExtraPoint(name='E2', type='EP', atomic_number=0)
        ep3 = ExtraPoint(name='E3', type='EP', atomic_number=0)
        ep4 = ExtraPoint(name='E4', type='EP', atomic_number=0)
        ep5 = ExtraPoint(name='E5', type='EP', atomic_number=0)
        self.assertIs(ep1.parent, None)
        self.assertEqual(ep1.bond_partners, [])
        self.assertEqual(ep1.angle_partners, [])
        self.assertEqual(ep1.dihedral_partners, [])
        self.assertEqual(ep1.tortor_partners, [])
        self.assertEqual(ep1.exclusion_partners, [])
        a1 = Atom(name='A1', type='AX', atomic_number=6)
        a2 = Atom(name='A2', type='AY', atomic_number=6)
        a3 = Atom(name='A3', type='AZ', atomic_number=7)
        a4 = Atom(name='A4', type='AX', atomic_number=6)
        a5 = Atom(name='A5', type='AY', atomic_number=6)
        bond_type = BondType(10.0, 1.0)
        bond_type2 = BondType(10.0, 2.0)
        bond_type3 = BondType(10.0, 0.5)
        bond_type4 = BondType(10.0, math.sqrt(2))
        angle_type = AngleType(10.0, 90)
        bonds = [Bond(a1, ep1, type=bond_type),
                 Bond(ep2, a2, type=bond_type),
                 Bond(a3, ep3, type=bond_type3),
                 Bond(a4, ep4, type=bond_type)]
        self.assertRaises(ValueError, lambda:
                ThreeParticleExtraPointFrame.from_weights(a2, a1, a3, 1.0, 1.0)
        )
        bonds.extend([Bond(a1, a2, type=bond_type),
                      Bond(a4, a3, type=bond_type4),
                      Bond(a3, a2, type=bond_type4),
                      Bond(a4, a5, type=bond_type2),
                      Bond(a5, ep5, type=bond_type)])
        self.assertRaises(ValueError, lambda:
                ThreeParticleExtraPointFrame.from_weights(a2, a1, a3, 1.0, 1.0)
        )
        angles = [Angle(a1, a2, a3, type=angle_type),
                  Angle(a2, a3, a4, type=angle_type),
                  Angle(a3, a4, a5)]
        dihedrals = [Dihedral(a1, a2, a3, a4), Dihedral(a2, a3, a4, a5)]
        self.assertEqual(ep1.weights, [1, 2])
        self.assertEqual(ep2.weights, None)
        self.assertIs(ep1.parent, a1)
        self.assertIs(ep2.parent, a2)
        self.assertIs(ep3.parent, a3)
        self.assertIs(ep4.parent, a4)
        self.assertIs(ep5.parent, a5)
        # Make sure exclusions are properly handled with EPs
        # EP1
        self.assertIn(a1, ep1.bond_partners)
        self.assertIn(a2, ep1.bond_partners)
        self.assertIn(ep2, ep1.bond_partners)
        self.assertIn(a3, ep1.angle_partners)
        self.assertIn(ep3, ep1.angle_partners)
        self.assertIn(a4, ep1.dihedral_partners)
        self.assertIn(ep4, ep1.dihedral_partners)
        # A1
        self.assertIn(ep1, a1.bond_partners)
        self.assertIn(a2, a1.bond_partners)
        self.assertIn(ep2, a1.bond_partners)
        self.assertIn(a3, a1.angle_partners)
        self.assertIn(ep3, a1.angle_partners)
        self.assertIn(a4, a1.dihedral_partners)
        self.assertIn(ep4, a1.dihedral_partners)
        self.assertNotIn(ep3, a1.bond_partners)
        self.assertNotIn(a5, a1.bond_partners+a1.angle_partners+a1.dihedral_partners)
        self.assertNotIn(a5, a1.exclusion_partners)
        # EP2
        self.assertIn(ep1, ep2.bond_partners)
        self.assertIn(a1, ep2.bond_partners)
        self.assertIn(a2, ep2.bond_partners)
        self.assertIn(a3, ep2.bond_partners)
        self.assertIn(ep3, ep2.bond_partners)
        self.assertIn(a4, ep2.angle_partners)
        self.assertIn(ep4, ep2.angle_partners)
        self.assertIn(a5, ep2.dihedral_partners)
        self.assertIn(ep5, ep2.dihedral_partners)
        # A2
        self.assertIn(ep1, a2.bond_partners)
        self.assertIn(ep2, a2.bond_partners)
        self.assertIn(ep3, a2.bond_partners)
        self.assertIn(a1, a2.bond_partners)
        self.assertIn(a3, a2.bond_partners)
        self.assertIn(a4, a2.angle_partners)
        self.assertIn(ep4, a2.angle_partners)
        self.assertIn(ep5, a2.dihedral_partners)
        self.assertIn(a5, a2.dihedral_partners)
        # EP3
        self.assertIn(ep2, ep3.bond_partners)
        self.assertIn(ep4, ep3.bond_partners)
        self.assertIn(a2, ep3.bond_partners)
        self.assertIn(a4, ep3.bond_partners)
        self.assertIn(a3, ep3.bond_partners)
        self.assertIn(ep1, ep3.angle_partners)
        self.assertIn(a1, ep3.angle_partners)
        self.assertIn(a5, ep3.angle_partners)
        self.assertIn(ep5, ep3.angle_partners)
        # A3
        self.assertIn(ep2, a3.bond_partners)
        self.assertIn(ep4, a3.bond_partners)
        self.assertIn(a2, a3.bond_partners)
        self.assertIn(a4, a3.bond_partners)
        self.assertIn(ep3, a3.bond_partners)
        self.assertIn(ep1, a3.angle_partners)
        self.assertIn(a1, a3.angle_partners)
        self.assertIn(a5, a3.angle_partners)
        self.assertIn(ep5, a3.angle_partners)
        # EP4
        self.assertIn(ep3, ep4.bond_partners)
        self.assertIn(a3, ep4.bond_partners)
        self.assertIn(a4, ep4.bond_partners)
        self.assertIn(a5, ep4.bond_partners)
        self.assertIn(ep5, ep4.bond_partners)
        self.assertIn(ep2, ep4.angle_partners)
        self.assertIn(a2, ep4.angle_partners)
        self.assertIn(ep1, ep4.dihedral_partners)
        self.assertIn(a1, ep4.dihedral_partners)
        # A4
        self.assertIn(ep3, a4.bond_partners)
        self.assertIn(a3, a4.bond_partners)
        self.assertIn(ep4, a4.bond_partners)
        self.assertIn(a5, a4.bond_partners)
        self.assertIn(ep5, a4.bond_partners)
        self.assertIn(ep2, a4.angle_partners)
        self.assertIn(a2, a4.angle_partners)
        self.assertIn(ep1, a4.dihedral_partners)
        self.assertIn(a1, a4.dihedral_partners)
        # EP5
        self.assertIn(ep4, ep5.bond_partners)
        self.assertIn(a4, ep5.bond_partners)
        self.assertIn(a5, ep5.bond_partners)
        self.assertIn(a3, ep5.angle_partners)
        self.assertIn(ep3, ep5.angle_partners)
        self.assertIn(a2, ep5.dihedral_partners)
        self.assertIn(ep2, ep5.dihedral_partners)
        # A5
        self.assertIn(ep4, a5.bond_partners)
        self.assertIn(a4, a5.bond_partners)
        self.assertIn(ep5, a5.bond_partners)
        self.assertIn(a3, a5.angle_partners)
        self.assertIn(ep3, a5.angle_partners)
        self.assertIn(a2, a5.dihedral_partners)
        self.assertIn(ep2, a5.dihedral_partners)
        # Test exclusions now
        a1.exclude(a5)
        self.assertIn(a1, a5.exclusion_partners)
        self.assertIn(a5, a1.exclusion_partners)
        self.assertIn(ep5, a1.exclusion_partners)
        self.assertIn(ep1, a5.exclusion_partners)
        self.assertIn(ep1, ep5.exclusion_partners)
        self.assertIn(ep5, ep1.exclusion_partners)
        self.assertIn(a1, ep5.exclusion_partners)
        self.assertIn(a5, ep1.exclusion_partners)
        # Test proper behavior of exclusion pruning
        a1.exclude(a2)
        self.assertIn(a1, a2._exclusion_partners)
        self.assertIn(a2, a1._exclusion_partners)
        a1.prune_exclusions()
        a2.prune_exclusions()
        self.assertNotIn(a2, a1._exclusion_partners)
        self.assertNotIn(a1, a2._exclusion_partners)
        # Now try TorsionTorsion
        tortor = TorsionTorsion(a1, a2, a3, a4, a5)
        self.assertIn(a1, a5.tortor_partners)
        self.assertIn(a5, a1.tortor_partners)
        self.assertIn(ep1, a5.tortor_partners)
        self.assertIn(ep5, a1.tortor_partners)
        self.assertIn(ep1, ep5.tortor_partners)
        self.assertIn(a1, ep5.tortor_partners)
        self.assertIn(a5, ep1.tortor_partners)
        self.assertIn(ep5, ep1.tortor_partners)
        # Exception handling
        self.assertRaises(MoleculeError, lambda: a1.exclude(a1))

        # Check the frame type. EP1 is TwoParticleExtraPointFrame
        ft1 = ep1.frame_type
        ft2 = ep2.frame_type
        ft3 = ep3.frame_type
        ft4 = ep4.frame_type
        ft5 = ep5.frame_type
        self.assertIsInstance(ft1, TwoParticleExtraPointFrame)
        self.assertIsInstance(ft2, ThreeParticleExtraPointFrame)
        self.assertIsInstance(ft3, ThreeParticleExtraPointFrame)
        self.assertIsInstance(ft4, ThreeParticleExtraPointFrame)
        self.assertIsInstance(ft5, TwoParticleExtraPointFrame)
        self.assertTrue(ft1.inside)
        self.assertFalse(ft2.inside)
        self.assertFalse(ft3.inside)
        self.assertFalse(ft4.inside)
        self.assertTrue(ft5.inside)
        # Now allow coordinates to dictate inside or outside. Make EP1 "outside"
        # and EP5 "inside"
        ep1.xx, ep1.xy, ep1.xz = -1.0, 0.0, 0.0
        a1.xx, a1.xy, a1.xz = 0.0, 0.0, 0.0
        a2.xx, a2.xy, a2.xz = 1.0, 0.0, 0.0
        ep5.xx, ep5.xy, ep5.xz = 9.0, 0.0, 0.0
        a5.xx, a5.xy, a5.xz = 10.0, 0.0, 0.0
        a4.xx, a4.xy, a4.xz = 8.0, 0.0, 0.0
        ep1._frame_type = ep5._frame_type = None
        ft1 = ep1.frame_type
        ft5 = ep5.frame_type
        self.assertFalse(ft1.inside)
        self.assertTrue(ft5.inside)
        self.assertIsInstance(ft5, TwoParticleExtraPointFrame)
        # Check the frame type API
        self.assertEqual(ft1.get_atoms(), (a1, a2))
        self.assertEqual(ft5.get_atoms(), (a5, a4))
        self.assertEqual(ft1.get_weights(), (2, -1))
        self.assertEqual(ft5.get_weights(), (0.5, 0.5))

        # Make EP3 inside
        a2.xx, a2.xy, a2.xz = 0.0, 0.0, 0.0
        a3.xx, a3.xy, a3.xz = 1.0, 1.0, 0.0
        a4.xx, a4.xy, a4.xz = 2.0, 0.0, 0.0
        ep3.xx, ep3.xy, ep3.xz = 1.0, 0.5, 0.0
        ep3._frame_type = None
        self.assertTrue(ep3.frame_type.inside)
        self.assertEqual(ep3.frame_type.get_atoms(), (a3, a4, a2))
        w1, w2, w3 = ep3.frame_type.get_weights()
        self.assertAlmostEqual(ep3.xx, a3.xx*w1+a2.xx*w2+a4.xx*w3)
        self.assertAlmostEqual(ep3.xy, a3.xy*w1+a2.xy*w2+a4.xy*w3)
        self.assertAlmostEqual(ep3.xz, a3.xz*w1+a2.xz*w2+a4.xz*w3)

        # Make EP3 outside
        ep3._frame_type = None
        ep3.xx, ep3.xy, ep3.xz = 1.0, 1.5, 0.0
        self.assertFalse(ep3.frame_type.inside)
        self.assertEqual(ep3.frame_type.get_atoms(), (a3, a4, a2))
        w1, w2, w3 = ep3.frame_type.get_weights()
        self.assertAlmostEqual(ep3.xx, a3.xx*w1+a2.xx*w2+a4.xx*w3)
        self.assertAlmostEqual(ep3.xy, a3.xy*w1+a2.xy*w2+a4.xy*w3)
        self.assertAlmostEqual(ep3.xz, a3.xz*w1+a2.xz*w2+a4.xz*w3)

        self.assertAlmostEqual(
                ThreeParticleExtraPointFrame.from_weights(a3, a2, a4, w2, w3), 0.5
        )
        self.assertRaises(ValueError, lambda:
                ThreeParticleExtraPointFrame.from_weights(a3, a2, a4, w2, w3-0.1)
        )

        # Check error checking for EP frame bonding patterns
        bonds.append(Bond(a1, a3))
        # ft1 is outdated, and should fail now
        self.assertRaises(RuntimeError, ft1.get_atoms)
        self.assertRaises(ValueError, ft1.get_weights)
        self.assertRaises(RuntimeError, ft3.get_atoms)
        self.assertRaises(ValueError, ft3.get_weights)
        # Create a duplicate bond and check ft5
        duplicate_bond = Bond(a4, a5, type=bond_type2)
        ep5._frame_type = None
        del a5.xx
        self.assertRaises(RuntimeError, lambda: ep5.frame_type)
        # Delete the duplicate bonds
        duplicate_bond.delete(); bonds.pop().delete()
        # Check that from_weights works again
        self.assertAlmostEqual(
                ThreeParticleExtraPointFrame.from_weights(a3, a2, a4, w2, w3), 0.5
        )
        types = []
        for bond in bonds:
            types.append(bond.type)
            bond.type = None
        self.assertRaises(ParameterError, lambda:
                ThreeParticleExtraPointFrame.from_weights(a3, a2, a4, w2, w3)
        )
        for bond, type in zip(bonds, types):
            if bond.atom1 is a3 and bond.atom2 is a2:
                continue
            bond.type = type
        # Still missing 1 bond type
        self.assertRaises(ParameterError, lambda:
                ThreeParticleExtraPointFrame.from_weights(a3, a2, a4, w2, w3)
        )
        for bond, type in zip(bonds, types):
            bond.type = type
        angles[1].type = None

    #=============================================

    def test_tip4p_ep_pattern(self):
        """ Tests handling of typical TIP4P-style extra points """
        o = Atom(name='O', atomic_number=8)
        h1 = Atom(name='H1', atomic_number=1)
        h2 = Atom(name='H2', atomic_number=1)
        ep = ExtraPoint(name='EP', atomic_number=0)
        bt = BondType(10.0, 1.0)
        epbt = BondType(10.0, 0.5)
        bonds = [Bond(o, h1, bt), Bond(o, h2, bt), Bond(o, ep, epbt), Bond(h1, h2)]
        self.assertRaises(ParameterError, lambda:
                ThreeParticleExtraPointFrame.from_weights(o, h1, h2, 1.0, 1.0)
        )
        self.assertAlmostEqual(
                ThreeParticleExtraPointFrame.from_weights(o, h1, h2, 0.5, 0.5,
                                                          dp1=1.0, dp2=1.0, d12=1.0),
                abs(math.cos(math.acos(0.5)*0.5))
        )
        self.assertAlmostEqual(
                ThreeParticleExtraPointFrame.from_weights(o, h1, h2, 0.5, 0.5,
                                                          dp1=1.0, dp2=1.0, theteq=45),
                abs(math.cos(math.pi/8))
        )
        self.assertRaises(ValueError, lambda:
                ThreeParticleExtraPointFrame.from_weights(o, h1, h2, 0.5, 0.5,
                                                          dp1=1.0, dp2=1.1, theteq=45)
        )

    #=============================================

    def test_copy_atom(self):
        """ Tests Atom's __copy__ ability """
        atom = Atom(name='CA', type='CX', atomic_number=6, mass=12.01)
        atom2 = Atom(name='CB', type='CY', atomic_number=6, mass=12.01)
        atype = AtomType('CX', 1, 12.01, 6)
        atype.set_lj_params(1.0, 1.2)
        atom.atom_type = atype
        bond = Bond(atom, atom2)
        # Copy the atom -- it should *not* inherit any valence terms
        copy1 = copy(atom)
        # Black hole conditions
        atom.xx = atom.xy = atom.xz = 0.0
        atom2.xx = atom2.xy = atom2.xz = 0.0
        atom.vx = atom.vy = atom.vz = 0.0
        atom2.vx = atom2.vy = atom2.vz = 0.0
        copy2 = copy(atom)
        # Now check that all attributes were transferred, they are different
        # objects, and that they did not inherit any bonds or bond_partners
        self.assertEqual(copy1.name, atom.name)
        self.assertEqual(copy1.mass, atom.mass)
        self.assertEqual(copy1.atomic_number, atom.atomic_number)
        self.assertEqual(copy1.type, atom.type)
        self.assertEqual(copy1.rmin, atom.rmin)
        self.assertEqual(copy1.epsilon, atom.epsilon)
        self.assertEqual(copy1.rmin_14, atom.rmin_14)
        self.assertEqual(copy1.epsilon_14, atom.epsilon_14)
        self.assertRaises(AttributeError, lambda: copy1.xx)
        self.assertRaises(AttributeError, lambda: copy1.vx)
        self.assertIs(copy1.atom_type, atom.atom_type) # shallow copy
        self.assertIsNot(copy1, atom)
        # Check copy 2
        self.assertEqual(copy2.name, atom.name)
        self.assertEqual(copy2.mass, atom.mass)
        self.assertEqual(copy2.atomic_number, atom.atomic_number)
        self.assertEqual(copy2.type, atom.type)
        self.assertEqual(copy2.xx, atom.xx)
        self.assertEqual(copy2.xy, atom.xy)
        self.assertEqual(copy2.xz, atom.xz)
        self.assertEqual(copy2.vx, atom.vx)
        self.assertEqual(copy2.vy, atom.vy)
        self.assertEqual(copy2.vz, atom.vz)
        self.assertEqual(copy2.rmin, atom.rmin)
        self.assertEqual(copy2.epsilon, atom.epsilon)
        self.assertEqual(copy2.rmin_14, atom.rmin_14)
        self.assertEqual(copy2.epsilon_14, atom.epsilon_14)
        self.assertIs(copy2.atom_type, atom.atom_type)
        self.assertIsNot(copy2, atom)
        self.assertIsNot(copy1, copy2)

    #=============================================

    def test_ljparams(self):
        """ Tests handling of Lennard-Jones parameters in Atom and AtomType """
        atom = Atom(name='CA', type='CX', atomic_number=6, mass=12.01)
        atype = AtomType('CX', 1, 12.01, 6)
        atype.set_lj_params(1.0, 1.2)
        self.assertEqual(atype.epsilon, 1.0)
        self.assertEqual(atype.rmin, 1.2)
        self.assertEqual(atype.epsilon_14, 1.0)
        self.assertEqual(atype.rmin_14, 1.2)
        self.assertEqual(atype.sigma, 1.2*2**(-1/6)*2)
        self.assertEqual(atype.sigma_14, atype.sigma)
        # Now try setting sigma and make sure it also changes rmin
        atype.sigma = 1.2
        self.assertEqual(atype.rmin, 1.2*2**(1/6)/2)
        atom.atom_type = atype
        self.assertEqual(atom.epsilon, atype.epsilon)
        self.assertEqual(atom.sigma, atype.sigma)
        self.assertEqual(atom.rmin, atype.rmin)

    #=============================================

    def test_bond(self):
        """ Tests the Bond and BondType classes """
        a1, a2, a3, a4 = Atom(), Atom(), Atom(), Atom()
        bond = Bond(a1, a2)
        # Test Bond as a container
        self.assertIn(a1, bond)
        self.assertIn(a2, bond)
        self.assertNotIn(a3, bond)
        self.assertNotIn(a4, bond)
        self.assertIn(bond, a1.bonds)
        self.assertIn(bond, a2.bonds)
        self.assertEqual(repr(bond), '<Bond %r--%r; type=None>' % (a1, a2))
        # Delete the bond and make sure it is deleted from the atoms, too
        bond.delete()
        self.assertEqual(len(a1.bonds), 0)
        self.assertNotIn(bond, a1.bonds)
        self.assertNotIn(a2, a1.bonds)
        self.assertNotIn(a1, a2.bonds)
        # Add multiple bonds
        bond1 = Bond(a1, a2)
        bond2 = Bond(a2, a3)
        bond3 = Bond(a3, a4)
        bond4 = Bond(a1, a2) # Duplicate
        bond1.delete()
        self.assertIn(a1, a2.bond_partners) # Duplicated bond
        bond4.delete()
        self.assertNotIn(a1, a2.bond_partners)
        # Now test measurements
        self.assertIs(bond1.measure(), None)
        self.assertIs(bond1.umeasure(), None)
        self.assertIs(bond1.energy(), None)
        self.assertIs(bond1.uenergy(), None)
        self.assertIs(bond2.measure(), None)
        self.assertIs(bond2.energy(), None)
        a2.xx = a2.xy = a2.xz = 0.0
        a3.xx, a3.xy, a3.xz = 1.0, 2.0, 3.0
        self.assertEqual(bond2.measure(), math.sqrt(14))
        self.assertEqual(bond2.umeasure(), math.sqrt(14) * u.angstroms)
        # Now test the bond types
        bond_types = TrackedList()
        bond_types.append(BondType(10.0, 1.0, bond_types))
        bond_types.append(BondType(12.0, 1.1, bond_types))
        bond_types.append(BondType(10.0, 1.0, bond_types))
        self.assertEqual(repr(bond_types[0]), '<BondType; k=10.000, req=1.000>')
        bond1 = Bond(a1, a2, bond_types[0])
        bond4 = Bond(a1, a2, bond_types[0])
        bond2.type = bond_types[1]
        bond3.type = bond_types[2]
        self.assertIs(bond1.type, bond4.type)
        self.assertIsNot(bond1.type, bond3.type)
        self.assertEqual(bond1.type, bond3.type)
        self.assertEqual(bond1.type.idx, 0)
        self.assertEqual(bond4.type.idx, 0)
        self.assertEqual(bond2.type.idx, 1)
        self.assertEqual(bond3.type.idx, 2)
        self.assertIs(bond1.type.list, bond_types)
        # Now test energy calculations
        ene = bond2.type.k*(bond2.measure() - bond2.type.req)**2
        self.assertAlmostEqual(bond2.energy(), ene)
        self.assertAlmostEqual(bond2.uenergy(), ene * u.kilocalories_per_mole)
        # Test the BondTypes.__copy__ method
        cp = copy(bond_types[0])
        self.assertIsNot(cp, bond_types[0])
        self.assertIs(cp.list, None)
        self.assertEqual(cp.idx, -1)
        self.assertEqual(cp.k, bond_types[0].k)
        self.assertEqual(cp.req, bond_types[0].req)
        # Make sure BondType is hashable
        bt1 = BondType(10.0, 1.0)
        bt2 = BondType(10.01, 1.0)
        bt3 = BondType(10.0, 1.01)
        bt4 = BondType(11.0, 1.1)
        self.assertEqual(hash(bt1), hash(BondType(10.0, 1.0)))
        self.assertEqual(hash(bt2), hash(BondType(10.01, 1.0)))
        self.assertEqual(hash(bt3), hash(BondType(10.0, 1.01)))
        self.assertEqual(hash(bt4), hash(BondType(11.0, 1.1)))
        # Check storage in a dict
        d = {bt1: 1, bt2: 2, bt3: 3, bt4: 4}
        self.assertEqual(1, d[BondType(10.0, 1.0)])
        self.assertEqual(2, d[BondType(10.01, 1.0)])
        self.assertEqual(3, d[BondType(10.0, 1.01)])
        self.assertEqual(4, d[BondType(11.0, 1.1)])
        # Error handling
        self.assertRaises(MoleculeError, lambda: Bond(a1, a1))
        self.assertRaises(MoleculeError, lambda: a1.bond_to(a1))
        self.assertRaises(KeyError, lambda: d[BondType(10.001, 1.0)])
        # Units
        bt = BondType(10.0 * u.kilojoules_per_mole / u.angstroms**2, 0.1 * u.nanometers)
        self.assertEqual(bt.k, 10.0 / 4.184)
        self.assertEqual(bt.req, 1)
        bt.k = 11
        self.assertEqual(bt.k, 11)
        bt.req = 1.1
        self.assertEqual(bt.req, 1.1)
        self.assertEqual(bt.uk, 11 * u.kilocalories_per_mole / u.angstroms**2)
        self.assertEqual(bt.ureq, 1.1 * u.angstroms)

    #=============================================

    def test_angle(self):
        """ Tests the Angle and AngleType classes """
        a1, a2, a3, a4 = Atom(), Atom(), Atom(), Atom()
        b1 = Bond(a1, a2)
        b2 = Bond(a2, a3)
        b3 = Bond(a3, a4)
        b4 = Bond(a4, a2)
        ang1 = Angle(a1, a2, a3)
        ang2 = Angle(a2, a3, a4)
        self.assertEqual(repr(ang1), '<Angle %r--%r--%r; type=None>' % (a1, a2, a3))
        # Now check measurements
        self.assertIs(ang1.measure(), None)
        self.assertIs(ang1.energy(), None)
        self.assertIs(ang1.umeasure(), None)
        self.assertIs(ang1.uenergy(), None)
        a1.xx = a1.xy = a1.xz = 0
        a2.xx, a2.xy, a2.xz = 0, 1, 0
        a3.xx, a3.xy, a3.xz = 1, 1, 0
        self.assertAlmostEqual(ang1.measure(), 90)
        self.assertAlmostEqual(ang1.umeasure().value_in_unit(u.degrees), 90)
        angle_types = TrackedList()
        # Make a set of angle types
        angle_types.append(AngleType(50.0, 109.5, angle_types))
        angle_types.append(AngleType(30.0, 120.0, angle_types))
        angle_types.append(AngleType(50.0, 109.5, angle_types))
        self.assertEqual(repr(angle_types[0]), '<AngleType; k=50.000, theteq=109.500>')
        # Assign the angle types to the angles (assigning to one new one)
        ang3 = Angle(a3, a4, a2, angle_types[2])
        ang1.type = angle_types[0]
        ang2.type = angle_types[1]
        # Now test energy
        self.assertAlmostEqual(ang1.energy(), ang1.type.k *
                               (ang1.type.theteq*DEG_TO_RAD - ang1.measure()*DEG_TO_RAD) ** 2)
        # Test angle as a container
        self.assertIn(a1, ang1)
        self.assertIn(a2, ang1)
        self.assertIn(a3, ang1)
        self.assertNotIn(a4, ang1)
        self.assertIn(b1, ang1)
        self.assertIn(b2, ang1)
        self.assertNotIn(b3, ang1)
        self.assertIn(b2, ang2)
        self.assertIn(b3, ang2)
        self.assertNotIn(b4, ang2)
        # Check that this angle is properly registered in each atom, and that
        # two atoms are not both angle and bond partners
        self.assertIn(a2, a1.bond_partners)
        self.assertIn(a3, a1.angle_partners)
        self.assertIn(a3, a2.bond_partners)
        self.assertNotIn(a2, a1.angle_partners)
        self.assertIn(ang1, a1.angles)
        self.assertIn(ang1, a2.angles)
        self.assertIn(ang1, a3.angles)
        self.assertIn(ang2, a2.angles)
        # Now check the types
        self.assertIs(ang1.type, angle_types[0])
        self.assertIsNot(ang3.type, ang1.type)
        self.assertEqual(ang3.type, ang1.type)
        self.assertNotEqual(ang1.type, ang2.type)
        self.assertEqual(angle_types[0].idx, 0)
        self.assertEqual(angle_types[1].idx, 1)
        self.assertEqual(angle_types[2].idx, 2)
        # Test the AngleType.__copy__ method
        cp = copy(angle_types[0])
        self.assertIsNot(cp, angle_types[0])
        self.assertIs(cp.list, None)
        self.assertEqual(cp.idx, -1)
        self.assertEqual(cp.k, angle_types[0].k)
        self.assertEqual(cp.theteq, angle_types[0].theteq)
        # Make sure BondType is hashable
        at1 = AngleType(10.0, 1.0)
        at2 = AngleType(10.01, 1.0)
        at3 = AngleType(10.0, 1.01)
        at4 = AngleType(11.0, 1.1)
        self.assertEqual(hash(at1), hash(AngleType(10.0, 1.0)))
        self.assertEqual(hash(at2), hash(AngleType(10.01, 1.0)))
        self.assertEqual(hash(at3), hash(AngleType(10.0, 1.01)))
        self.assertEqual(hash(at4), hash(AngleType(11.0, 1.1)))
        # Check storage in a dict
        d = {at1: 1, at2: 2, at3: 3, at4: 4}
        self.assertEqual(1, d[AngleType(10.0, 1.0)])
        self.assertEqual(2, d[AngleType(10.01, 1.0)])
        self.assertEqual(3, d[AngleType(10.0, 1.01)])
        self.assertEqual(4, d[AngleType(11.0, 1.1)])
        # Error handling
        self.assertRaises(KeyError, lambda: d[AngleType(10.1, 1.0)])
        self.assertRaises(KeyError, lambda: d[BondType(10.0, 1.0)])
        self.assertRaises(MoleculeError, lambda: Angle(a1, a2, a1))
        self.assertRaises(MoleculeError, lambda: a1.angle_to(a1))

    #=============================================

    def test_dihedral(self):
        """ Tests the Dihedral and DihedralType classes """
        atoms = AtomList()
        atoms.extend([Atom(list=atoms) for i in range(10)])
        bonds = []
        # Sequential bonds
        for i in range(len(atoms)-2):
            bonds.append(Bond(atoms[i], atoms[i+1]))
        bonds.append(Bond(atoms[-1], atoms[1])) # to make an improper
        bonds.append(Bond(atoms[0], atoms[-1])) # stupid bond
        dihed_types = TrackedList()
        dihed_types.append(DihedralType(5.0, 2, 0.0, 1.2, 2.0, dihed_types))
        dihed_types.append(DihedralType(1.0, 3, 180.0, 1.2, 2.0, dihed_types))
        dihed_types.append(DihedralType(2.0, 4, 180.0, 1.2, 2.0, dihed_types))
        dihed_types.append(DihedralType(10.0, 1, 180.0, 0., 0., dihed_types))
        d1 = Dihedral(atoms[0], atoms[-1], atoms[1], atoms[2], improper=True, type=dihed_types[3])
        d2 = Dihedral(atoms[0], atoms[1], atoms[2], atoms[3], type=dihed_types[0])
        d3 = Dihedral(atoms[0], atoms[1], atoms[2], atoms[3], ignore_end=True, type=dihed_types[1])
        d4 = Dihedral(atoms[1], atoms[2], atoms[3], atoms[4], type=dihed_types[2])
        self.assertTrue(d1.same_atoms([0, len(atoms)-1, 1, 2]))
        self.assertFalse(d1.same_atoms([0, 1, 2, 3]))
        self.assertRaises(TypeError, lambda: d1.same_atoms([0, 1, 2]))
        self.assertEqual(repr(d1), '<Dihedral [imp]; %r--%r--%r--%r; type=%r>' %
                (atoms[0], atoms[-1], atoms[1], atoms[2], dihed_types[3]))
        self.assertEqual(repr(d2), '<Dihedral; %r--%r--%r--%r; type=%r>' %
                (atoms[0], atoms[1], atoms[2], atoms[3], dihed_types[0]))
        self.assertEqual(repr(d3), '<Dihedral [ign]; %r--%r--%r--%r; type=%r>' %
                (atoms[0], atoms[1], atoms[2], atoms[3], dihed_types[1]))
        # Test container properties
        self.assertIn(bonds[-2], d1)
        self.assertIn(bonds[0], d1)
        self.assertIn(bonds[1], d1)
        self.assertNotIn(bonds[-1], d1)
        self.assertIn(bonds[0], d2)
        self.assertIn(bonds[1], d2)
        self.assertIn(bonds[2], d2)
        self.assertNotIn(bonds[3], d2)
        for i, a in enumerate(atoms):
            if i < 4:
                self.assertIn(a, d2)
                self.assertIn(a, d3)
                if i > 0:
                    self.assertIn(a, d4)
                else:
                    self.assertNotIn(a, d4)
            else:
                if i == 4:
                    self.assertIn(a, d4)
                else:
                    self.assertNotIn(a, d4)
                self.assertNotIn(a, d2)
                self.assertNotIn(a, d3)
        # Check the various attributes of the different dihedrals
        self.assertTrue(d1.improper and d1.ignore_end)
        self.assertTrue(d3.ignore_end and not d3.improper)
        self.assertFalse(d2.ignore_end or d2.improper)
        self.assertFalse(d4.ignore_end or d4.improper)
        self.assertEqual(d1.type.idx, 3)
        self.assertEqual(d2.type.idx, 0)
        self.assertEqual(d3.type.idx, 1)
        self.assertEqual(d4.type.idx, 2)
        self.assertEqual(d1.type.phi_k, 10.0)
        self.assertEqual(d1.type.per, 1)
        self.assertEqual(d1.type.phase, 180)
        self.assertEqual(d1.type.scee, 0.)
        self.assertEqual(d1.type.scnb, 0.)

        self.assertTrue(d1.same_atoms((0, 9, 1, 2)))
        self.assertFalse(d1.same_atoms((0, 1, 9, 2)))
        self.assertTrue(d1.same_atoms(d1))
        self.assertFalse(d1.same_atoms(d2))
        self.assertTrue(d2.same_atoms((0, 1, 2, 3)))
        self.assertTrue(d2.same_atoms((3, 2, 1, 0)))
        self.assertFalse(d2.same_atoms((2, 3, 1, 0)))
        self.assertNotIn(atoms[0], atoms[1].dihedral_partners)
        d2.delete()
        self.assertIn(atoms[0], atoms[3].dihedral_partners)
        self.assertIn(atoms[3], atoms[0].dihedral_partners)
        # Test the measurement capabilities
        self.assertIs(d3.measure(), None)
        self.assertIs(d3.energy(), None)
        self.assertIs(d3.umeasure(), None)
        self.assertIs(d3.uenergy(), None)
        a1, a2, a3, a4 = atoms[:4]
        a1.xx, a1.xy, a1.xz = 1, 0, 0
        a2.xx, a2.xy, a2.xz = 0, 0, 0
        a3.xx, a3.xy, a3.xz = 0, 0, 1
        a4.xx, a4.xy, a4.xz = 0.1, 0.6, 1
        self.assertAlmostEqual(d3.measure(), 80.537677791974374)
        self.assertAlmostEqual(d3.umeasure(), 80.537677791974374 * u.degrees)
        ene = d3.type.phi_k*(1+math.cos(d3.type.per*d3.measure()*DEG_TO_RAD-d3.type.phase*DEG_TO_RAD))
        self.assertAlmostEqual(d3.energy(), ene)
        self.assertEqual(d3.uenergy(), ene*u.kilocalories_per_mole)

        # Now delete dihedral 3
        d3.delete()
        self.assertNotIn(atoms[0], atoms[3].dihedral_partners)
        self.assertNotIn(atoms[3], atoms[0].dihedral_partners)
        # Test the DihedralType.__copy__ method
        cp = copy(dihed_types[0])
        self.assertIsNot(cp, dihed_types[0])
        self.assertIsNot(cp, dihed_types[0])
        self.assertIs(cp.list, None)
        self.assertEqual(cp.idx, -1)
        self.assertEqual(cp.phi_k, dihed_types[0].phi_k)
        self.assertEqual(cp.per, dihed_types[0].per)
        self.assertEqual(cp.phase, dihed_types[0].phase)
        self.assertEqual(cp.scee, dihed_types[0].scee)
        self.assertEqual(cp.scnb, dihed_types[0].scnb)
        # Test hashing
        self.assertEqual(hash(dihed_types[0]), hash(DihedralType(5, 2, 0, 1.2, 2.0)))
        self.assertEqual(hash(dihed_types[1]), hash(DihedralType(1.0, 3, 180.0, 1.2, 2.0, dihed_types)))
        self.assertEqual(hash(dihed_types[2]), hash(DihedralType(2.0, 4, 180.0, 1.2, 2.0, dihed_types)))
        self.assertEqual(hash(dihed_types[3]), hash(DihedralType(10.0, 1, 180.0, 0., 0., dihed_types)))
        d = dict(zip(dihed_types, range(len(dihed_types))))
        for i, dt in enumerate(dihed_types):
            self.assertEqual(d[dt], i)
        # Unit handling
        dt = DihedralType(5 * u.kilojoules_per_mole, 5, math.pi / 4 * u.radians)
        self.assertEqual(dt.phi_k, 5 / 4.184)
        self.assertEqual(dt.per, 5)
        self.assertEqual(dt.phase, 45)
        self.assertEqual(dt.uphi_k, (5 / 4.184) * u.kilocalories_per_mole)
        self.assertEqual(dt.uphase, 45 * u.degrees)
        # Exception handling
        self.assertRaises(KeyError, lambda: d[BondType(10.0, 1.0)])
        self.assertRaises(MoleculeError, lambda: Dihedral(atoms[0], atoms[1], atoms[2], atoms[0]))
        self.assertRaises(MoleculeError, lambda: atoms[0].dihedral_to(atoms[0]))

    #=============================================

    def test_dihedral_type_list(self):
        """ Tests the DihedralTypeList class """
        dihed_types = TrackedList()
        dihed_types.append(DihedralTypeList(list=dihed_types))
        dihed_types[0].append(DihedralType(5.0, 2, 0.0, 1.2, 2.0))
        dihed_types[0].append(DihedralType(1.0, 3, 180.0, 1.2, 2.0))
        dihed_types[0].append(DihedralType(2.0, 4, 180.0, 1.2, 2.0))
        dihed_types[0].append(DihedralType(10.0, 1, 180.0, 0., 0.))
        self.assertRaises(TypeError, lambda: dihed_types[0].append(Atom()))
        self.assertRaises(ParameterError, lambda: dihed_types[0].append(DihedralType(1.0, 2, 0.0, 1.2, 2.0)))
        self.assertIs(dihed_types[0].list, dihed_types)
        self.assertEqual(len(dihed_types), 1)
        self.assertEqual(dihed_types[0].idx, 0)
        self.assertEqual(len(dihed_types[0]), 4)
        self.assertEqual(repr(dihed_types[0]), '<DihedralTypes %s>' % list.__repr__(dihed_types[0]))
        # Now try DihedralTypeList.from_rbtorsion
        dtl = DihedralTypeList.from_rbtorsion(RBTorsionType(1, 2, 3, -4, -1, 0))

        # Now test DihedralTypeList.__copy__
        cp = copy(dihed_types[0])
        self.assertIsNot(cp, dihed_types[0])
        self.assertIs(cp.list, None)
        self.assertEqual(cp.idx, -1)
        self.assertEqual(len(cp), 4)
        i = 0
        for t1, t2 in zip(cp, dihed_types[0]):
            self.assertIsNot(t1, t2)
            self.assertEqual(t1.phi_k, t2.phi_k)
            self.assertEqual(t1.per, t2.per)
            self.assertEqual(t1.phase, t2.phase)
            self.assertEqual(t1.scee, t2.scee)
            self.assertEqual(t1.scnb, t2.scnb)
            i += 1
        self.assertEqual(i, 4)
        # Make sure DihedralTypeList is hashable
        hash(dihed_types[0])
        dtcp = copy(dihed_types[0])
        self.assertIsNot(dtcp, dihed_types)
        self.assertEqual(hash(dtcp), hash(dihed_types[0]))

        # check dihedral angle and energy
        atoms = AtomList()
        atoms.extend([Atom(list=atoms) for i in range(4)])

        f = [  0.,   -19.73, -13.12,   1.53,  -6.56]
        #
        # Coordinates HOOH
        #
        a1, a2, a3, a4 = atoms[:4]
        a1.xx, a1.xy, a1.xz = 0.000000,    0.000000,   -0.000000
        a2.xx, a2.xy, a2.xz = 0.000000,   -0.000000,    0.966700
        a3.xx, a3.xy, a3.xz = 1.425871,    0.000000,    1.237188
        a4.xx, a4.xy, a4.xz = 1.532691,    0.871896,    1.640791

        dihedral_list = DihedralTypeList()
        dihedral_list.append(DihedralType(f[0], 0, 180.0))
        dihedral_list.append(DihedralType(f[1], 1, 180.0))
        dihedral_list.append(DihedralType(f[2], 2, 180.0))
        dihedral_list.append(DihedralType(f[3], 3, 180.0))
        dihedral_list.append(DihedralType(f[4], 4, 180.0))
        d = Dihedral(atoms[0], atoms[1], atoms[2], atoms[3], ignore_end=True, type=dihedral_list)
        self.assertAlmostEqual(d.measure(), 113.362320189)
        self.assertAlmostEqual(d.energy(), -28.2654260649)

    #=============================================

    @unittest.skipUnless(has_openmm, "Cannot test without OpenMM")
    def test_rb_torsion_type_conversion_openmm(self):
        """ Test energetics/forces of converted R-B torsion """
        s1 = Structure()
        s2 = Structure()
        def add_atoms(struct):
            struct.add_atom(Atom(name='A1', atomic_number=6), 'ALA', 1)
            struct.add_atom(Atom(name='A2', atomic_number=6), 'ALA', 1)
            struct.add_atom(Atom(name='A3', atomic_number=6), 'ALA', 1)
            struct.add_atom(Atom(name='A4', atomic_number=6), 'ALA', 1)
        add_atoms(s1)
        add_atoms(s2)
        crd = np.random.rand(4, 3)
        s1.coordinates = crd.copy()
        s1.rb_torsions.append(Dihedral(*s1.atoms))
        dt = RBTorsionType(1, 2, 3, -4, -1, 0)
        s1.rb_torsions[0].type = dt
        s1.rb_torsion_types.append(dt)
        s1.rb_torsion_types.claim()

        s2.dihedrals.append(Dihedral(*s2.atoms))
        dt = DihedralTypeList.from_rbtorsion(dt)
        s2.dihedral_types.append(dt)
        s2.dihedrals[0].type = dt
        s2.coordinates = crd.copy()
        s2.dihedral_types.claim()

    #=============================================

    def test_rb_torsion_type(self):
        """ Tests the RBTorsionType class """
        rb_types = TrackedList()
        rb_types.append(RBTorsionType(10, 20, 30, 40, 50, 60))
        rb_types.append(RBTorsionType(11, 21, 31, 41, 51, 61))
        rb_types.append(RBTorsionType(12, 22, 32, 42, 52, 62, list=rb_types))
        self.assertEqual(repr(rb_types[0]),
                '<RBTorsionType; c0=10.000; c1=20.000; c2=30.000; c3=40.000; c4=50.000;'
                ' c5=60.000; scee=1.0; scnb=1.0>'
        )
        rb_types.claim()
        for i, rb_typ in enumerate(rb_types):
            self.assertEqual(i, rb_typ.idx)
            self.assertEqual(rb_typ.c0, 10+i)
            self.assertEqual(rb_typ.c1, 20+i)
            self.assertEqual(rb_typ.c2, 30+i)
            self.assertEqual(rb_typ.c3, 40+i)
            self.assertEqual(rb_typ.c4, 50+i)
            self.assertEqual(rb_typ.c5, 60+i)
        # Test RBTorsion.__copy__
        cp = copy(rb_types[0])
        self.assertIsNot(cp, rb_types[0])
        self.assertIs(cp.list, None)
        self.assertEqual(cp.idx, -1)
        self.assertEqual(cp.c0, 10)
        self.assertEqual(cp.c1, 20)
        self.assertEqual(cp.c2, 30)
        self.assertEqual(cp.c3, 40)
        self.assertEqual(cp.c4, 50)
        self.assertEqual(cp.c5, 60)
        # Check that RBTorsionType is hashable
        for rb_type in rb_types:
            cp = copy(rb_type)
            self.assertIsNot(rb_type, cp)
            self.assertEqual(hash(rb_type), hash(cp))
            self.assertIn(rb_type, {cp})
        # Check unit handling
        rb = RBTorsionType(1*u.kilocalories_per_mole, 2*u.kilocalories_per_mole,
                           3*u.kilocalories_per_mole, 4*u.kilocalorie_per_mole,
                           5*u.kilocalories_per_mole, 6*u.kilocalorie_per_mole)
        self.assertEqual(rb.c0, 1)
        self.assertEqual(rb.c1, 2)
        self.assertEqual(rb.c2, 3)
        self.assertEqual(rb.c3, 4)
        self.assertEqual(rb.c4, 5)
        self.assertEqual(rb.c5, 6)
        rb.c0 = rb.c1 = rb.c2 = rb.c3 = rb.c4 = rb.c5 = 1
        self.assertEqual(rb.c0, 1)
        self.assertEqual(rb.c1, 1)
        self.assertEqual(rb.c2, 1)
        self.assertEqual(rb.c3, 1)
        self.assertEqual(rb.c4, 1)
        self.assertEqual(rb.c5, 1)
        self.assertEqual(rb.uc0, 1 * u.kilocalories_per_mole)
        self.assertEqual(rb.uc1, 1 * u.kilocalories_per_mole)
        self.assertEqual(rb.uc2, 1 * u.kilocalories_per_mole)
        self.assertEqual(rb.uc3, 1 * u.kilocalories_per_mole)
        self.assertEqual(rb.uc4, 1 * u.kilocalories_per_mole)
        self.assertEqual(rb.uc5, 1 * u.kilocalories_per_mole)

    #=============================================

    def test_urey_bradley(self):
        """ Tests the Urey-Bradley term """
        atoms = TrackedList()
        atoms.extend([Atom(), Atom(), Atom()])
        b1 = Bond(atoms[0], atoms[1])
        b2 = Bond(atoms[1], atoms[2])
        b3 = Bond(atoms[0], atoms[2])
        a1 = Angle(atoms[0], atoms[1], atoms[2])
        u1 = UreyBradley(atoms[0], atoms[2], BondType(10.0, 1.0))
        self.assertEqual(repr(u1), '<UreyBradley %r--%r; type=%r>' % (atoms[0], atoms[2], u1.type))
        # Test containers
        self.assertIn(atoms[0], u1)
        self.assertIn(atoms[2], u1)
        self.assertNotIn(atoms[1], u1)
        self.assertIn(b1, u1)
        self.assertIn(b2, u1)
        self.assertNotIn(b3, u1)
        self.assertEqual(u1.type.k, 10)
        self.assertEqual(u1.type.req, 1)
        # Test urey-bradley measurement and energy
        self.assertIs(u1.measure(), None)
        self.assertIs(u1.energy(), None)
        atoms[0].xx = atoms[0].xy = atoms[0].xz = 0
        atoms[2].xx = atoms[2].xy = atoms[2].xz = 1
        self.assertAlmostEqual(u1.measure(), math.sqrt(3))
        self.assertAlmostEqual(u1.energy(), u1.type.k*(u1.measure()-u1.type.req)**2)
        # Test error checking
        self.assertRaises(MoleculeError, lambda: UreyBradley(atoms[0], atoms[0]))

    #=============================================

    def test_improper(self):
        """ Tests the CHARMM improper torsion term and type """
        atoms = AtomList()
        atoms.extend([Atom(), Atom(), Atom(), Atom()])
        b1 = Bond(atoms[0], atoms[1])
        b2 = Bond(atoms[0], atoms[2])
        b3 = Bond(atoms[0], atoms[3])
        imp_types = TrackedList()
        imp_types.append(ImproperType(10.0, 180.0, list=imp_types))
        imp_types.append(ImproperType(10.0, 180.0, list=imp_types))
        imp = Improper(atoms[0], atoms[1], atoms[2], atoms[3], imp_types[0])
        imp2 = Improper(atoms[0], atoms[3], atoms[1], atoms[2], imp_types[1])
        self.assertEqual(repr(imp), '<Improper; %r--(%r,%r,%r); type=%r>' %
                         (atoms[0], atoms[1], atoms[2], atoms[3], imp_types[0]))
        for a in atoms:
            self.assertIn(a, imp)
        self.assertEqual(imp.type.idx, 0)
        self.assertEqual(imp.type.psi_k, 10.0)
        self.assertEqual(imp.type.psi_eq, 180.0)
        imp.type.psi_eq = 3 * 180 / math.pi
        imp.type.psi_k = 1
        self.assertEqual(imp.type.upsi_k, 1 * u.kilocalories_per_mole / u.radians**2)
        self.assertEqual(imp.type.upsi_eq, 3 * 180 / math.pi * u.degrees)
        self.assertTrue(imp.same_atoms((0, 1, 2, 3)))
        self.assertTrue(imp.same_atoms((0, 2, 3, 1)))
        self.assertFalse(imp.same_atoms((3, 2, 1, 0)))
        self.assertTrue(imp.same_atoms(imp2))
        self.assertTrue(imp.same_atoms(imp))
        # Test ImproperType.__copy__
        cp = copy(imp_types[0])
        self.assertIsNot(cp, imp_types[0])
        self.assertIs(cp.list, None)
        self.assertEqual(cp.idx, -1)
        self.assertEqual(cp.psi_k, imp_types[0].psi_k)
        self.assertEqual(cp.psi_eq, imp_types[0].psi_eq)
        # Now test that same_atoms detects *not-same* improper atoms
        imp3 = Improper(atoms[1], atoms[0], atoms[2], atoms[3], imp_types[0])
        self.assertFalse(imp.same_atoms(imp3))
        self.assertFalse(imp3.same_atoms(imp))
        self.assertRaises(MoleculeError, lambda: imp3.same_atoms([0, 1, 2]))
        for imp_type in imp_types:
            cp = copy(imp_type)
            self.assertIsNot(cp, imp_type)
            self.assertEqual(hash(cp), hash(imp_type))
            self.assertIn(cp, {imp_type})

    #=============================================

    def test_cmap(self):
        """ Tests the coupled-torsion CMAP terms used by CHARMM """
        atoms = TrackedList()
        atoms.extend([Atom(list=atoms) for i in range(10)])
        bonds = [Bond(atoms[i], atoms[i+1]) for i in range(9)]
        cg1 = list(range(36))
        cg2 = list(reversed(range(36)))
        cmap_types = TrackedList()
        cmap_types.append(CmapType(6, cg1, list=cmap_types))
        cmap_types.append(CmapType(6, cg2, list=cmap_types))
        cmap_types.append(CmapType(6, cg1[:], list=cmap_types))
        cmap1 = Cmap.extended(atoms[0], atoms[1], atoms[2], atoms[3],
                              atoms[1], atoms[2], atoms[3], atoms[4],
                              cmap_types[0])
        cmap2 = Cmap(atoms[3], atoms[4], atoms[5], atoms[6], atoms[7],
                     cmap_types[1])
        cmap3 = Cmap(atoms[5], atoms[6], atoms[7], atoms[8], atoms[9],
                     cmap_types[2])
        self.assertEqual(repr(cmap1), '<Cmap; %r--%r--%r--%r--%r; type=%r>' %
                         tuple(atoms[:5] + [cmap_types[0]]))
        # Check illegal CmapType assignment
        self.assertRaises(TypeError, lambda: CmapType(5, cg1))
        self.assertRaises(NotImplementedError, lambda:
                Cmap.extended(atoms[0], atoms[1], atoms[2], atoms[3],
                              atoms[4], atoms[5], atoms[6], atoms[7])
        )
        self.assertRaises(MoleculeError, lambda:
                cmap1.same_atoms((1, 2, 3, 4))
        )
        # Check container functionality
        for i, atom in enumerate(atoms):
            if 0 <= i < 5:
                self.assertIn(atom, cmap1)
            if 3 <= i < 7:
                self.assertIn(atom, cmap2)
            if 5 <= i < 9:
                self.assertIn(atom, cmap3)
        self.assertEqual(cmap_types[0], cmap_types[2])
        self.assertNotEqual(cmap_types[0], cmap_types[1])
        for i, bond in enumerate(bonds):
            if 0 <= i < 4:
                self.assertIn(bond, cmap1)
            if 3 <= i < 7:
                self.assertIn(bond, cmap2)
            if 5 <= i < 9:
                self.assertIn(bond, cmap3)
        # Check the same_atoms method
        self.assertTrue(cmap1.same_atoms(cmap1))
        self.assertFalse(cmap2.same_atoms(cmap1))
        self.assertTrue(cmap1.same_atoms((0, 1, 2, 3, 4)))
        self.assertTrue(cmap1.same_atoms((4, 3, 2, 1, 0)))
        self.assertFalse(cmap1.same_atoms((0, 2, 1, 3, 4)))
        # Check Cmap grid indexing
        for i, val in enumerate(cg1):
            self.assertEqual(cmap_types[0].grid[i], val)
            i1 = i // 6
            i2 = i % 6
            self.assertEqual(cmap_types[0].grid[i1,i2], val)
        transpose = cmap_types[0].grid.T
        for i, val in enumerate(cg1):
            i1 = i // 6
            i2 = i % 6
            self.assertEqual(transpose[i2, i1], val)
        # Transpose of the transpose should be the original
        self.assertEqual(transpose.T, cmap_types[0].grid)
        # Check caching
        self.assertIs(transpose, cmap_types[0].grid.T)
        # Check switch_range functionality
        sg = cmap_types[0].grid.switch_range()
        for i in range(6):
            for j in range(6):
                self.assertEqual(sg[i,j], cmap_types[0].grid[(i+3)%6,(j+3)%6])
        # Check the CmapType.__copy__ functionality
        cp = copy(cmap_types[0])
        self.assertIsNot(cp, cmap_types[0])
        self.assertIs(cp.list, None)
        self.assertEqual(cp.idx, -1)
        self.assertEqual(cp.resolution, cmap_types[0].resolution)
        self.assertEqual(cp.grid, cmap_types[0].grid)
        # Check hashability of CmapType
        for cmap_type in cmap_types:
            cp = copy(cmap_type)
            self.assertIsNot(cp, cmap_type)
            self.assertIsNot(cp.grid, cmap_type.grid)
            self.assertEqual(hash(cp), hash(cmap_type))
            self.assertIn(cp, {cmap_type})
            self.assertIn(cp.grid, {cmap_type.grid})
        # Check error handling
        self.assertRaises(IndexError, lambda: cmap_types[0].grid[10,10])
        def err():
            cmap_types[0].grid[10, 10] = 10.0
        self.assertRaises(IndexError, err)

    #=============================================

    def test_cmap_grid(self):
        """ Tests the CmapGrid API """
        raw = list(range(36))
        cg1 = CmapType(6, list(range(36))).grid
        for i in range(36):
            self.assertEqual(cg1[i], i)
        for i in range(6):
            for j in range(6):
                self.assertEqual(cg1[i,j], i*6+j)
        # Now check *setting* elements
        for i, j in enumerate(reversed(range(36))):
            cg1[i] = j
        for i in range(36):
            self.assertEqual(cg1[i], 35-i)
        for i in range(6):
            for j in range(6):
                cg1[i,j] = i*10+j*6
        for i in range(6):
            for j in range(6):
                self.assertEqual(cg1[i,j], i*10+j*6)
        # Now try setting slices
        cg1[:10] = 0
        for i in range(10):
            self.assertEqual(cg1[i], 0)
        for i in range(10, 36):
            self.assertNotEqual(cg1[i], 0)
        # Now set slice to an array (akin to numpy broadcasting)
        cg1[:10] = [(x+4)*10 for x in range(10)]
        for i in range(10):
            self.assertEqual(cg1[i], i*10+40)
        # Test error handling for improper broadcasting
        def err():
            cg1[:10] = [1, 2, 3]
        self.assertRaises(ValueError, err)
        # Test string methods
        self.assertEqual(repr(cg1), '<_CmapGrid: 6x6>')

        # Test comparison methods
        cg2 = CmapType(7, list(range(49))).grid
        self.assertNotEqual(cg1, cg2)
        cg1[:] = list(range(36))
        cg3 = copy(cg1)
        self.assertEqual(cg1, cg3)
        cg3[-1] += 0.01
        self.assertNotEqual(cg1, cg3)
        self.assertNotEqual(cg1, 1)

    #=============================================

    def test_trigonal_angle(self):
        """ Tests the trigonal angle term used in the AMOEBA force field """
        atoms = TrackedList()
        atoms.extend([Atom(list=atoms) for x in range(8)])
        bonds = []
        bonds.append(Bond(atoms[0], atoms[1]))
        bonds.append(Bond(atoms[2], atoms[1]))
        bonds.append(Bond(atoms[3], atoms[1]))
        bonds.append(Bond(atoms[0], atoms[2]))
        t1 = TrigonalAngle(atoms[0], atoms[1], atoms[2], atoms[3])
        t2 = TrigonalAngle(atoms[4], atoms[5], atoms[6], atoms[7],
                           AngleType(50.0, 90.0))
        t1.type = AngleType(50.0, 90.0)
        self.assertIsNot(t1.type, t2.type)
        self.assertEqual(t1.type, t2.type)
        self.assertIn(atoms[0], t1)
        self.assertIn(atoms[1], t1)
        self.assertIn(atoms[2], t1)
        self.assertIn(atoms[3], t1)
        self.assertIn(atoms[4], t2)
        self.assertIn(atoms[5], t2)
        self.assertIn(atoms[6], t2)
        self.assertIn(atoms[7], t2)
        self.assertIn(atoms[0], t1)
        self.assertIn(bonds[0], t1)
        self.assertIn(bonds[1], t1)
        self.assertIn(bonds[2], t1)
        self.assertNotIn(bonds[3], t1)
        self.assertRaises(MoleculeError, lambda:
                          TrigonalAngle(atoms[0], atoms[1], atoms[2], atoms[1]))
        self.assertEqual(t1.type.k, 50)
        self.assertEqual(t1.type.theteq, 90.0)
        self.assertEqual(repr(t1), '<TrigonalAngle; %r--(%r,%r,%r); type=%r>' %
                (atoms[1], atoms[0], atoms[2], atoms[3], t1.type))
        self.assertEqual(repr(t2), '<TrigonalAngle; %r--(%r,%r,%r); type=%r>' %
                (atoms[5], atoms[4], atoms[6], atoms[7], t2.type))

    #=============================================

    def test_oopbend_angle(self):
        """ Tests the out-of-plane bend term used in the AMOEBA force field """
        atoms = TrackedList()
        atoms.extend([Atom(list=atoms) for x in range(8)])
        bonds = []
        bonds.append(Bond(atoms[0], atoms[1]))
        bonds.append(Bond(atoms[2], atoms[1]))
        bonds.append(Bond(atoms[3], atoms[1]))
        bonds.append(Bond(atoms[0], atoms[2]))
        t1 = OutOfPlaneBend(atoms[0], atoms[1], atoms[2], atoms[3])
        t2 = OutOfPlaneBend(atoms[4], atoms[5], atoms[6], atoms[7], OutOfPlaneBendType(50.0))
        t1.type = OutOfPlaneBendType(50.0)
        self.assertIsNot(t1.type, t2.type)
        self.assertEqual(t1.type, t2.type)
        self.assertIn(atoms[0], t1)
        self.assertIn(atoms[1], t1)
        self.assertIn(atoms[2], t1)
        self.assertIn(atoms[3], t1)
        self.assertIn(atoms[4], t2)
        self.assertIn(atoms[5], t2)
        self.assertIn(atoms[6], t2)
        self.assertIn(atoms[7], t2)
        self.assertIn(bonds[0], t1)
        self.assertIn(bonds[1], t1)
        self.assertIn(bonds[2], t1)
        self.assertNotIn(bonds[3], t1)
        self.assertRaises(MoleculeError, lambda: OutOfPlaneBend(atoms[0], atoms[1], atoms[2], atoms[1]))
        self.assertEqual(t1.type.k, 50)
        self.assertEqual(repr(t1), '<OutOfPlaneBend; %r--(%r,%r,%r); type=%r>' %
                                   (atoms[1], atoms[0], atoms[2], atoms[3], t1.type))
        cp = OutOfPlaneBendType(50.0)
        self.assertIsNot(cp, t1.type)
        self.assertEqual(hash(cp), hash(t1.type))
        self.assertIn(cp, {t1.type})

    #=============================================

    def test_pitorsion(self):
        """ Tests the pi-torsion term used in the AMOEBA force field """
        atoms = TrackedList()
        atoms.extend([Atom(list=atoms) for i in range(6)])
        bonds = [Bond(atoms[0], atoms[2]), Bond(atoms[1], atoms[2]),
                 Bond(atoms[2], atoms[3]), Bond(atoms[3], atoms[4]),
                 Bond(atoms[3], atoms[5])]
        pit = PiTorsion(atoms[0], atoms[1], atoms[2], atoms[3], atoms[4],
                        atoms[5], type=DihedralType(10.0, 2, 180.0))
        for atom in atoms:
            self.assertIn(atom, pit)
        for bond in bonds:
            self.assertIn(bond, pit)
        self.assertNotIn(Bond(atoms[0], atoms[1]), pit)
        self.assertEqual(pit.type.phi_k, 10)
        self.assertEqual(pit.type.per, 2)
        self.assertEqual(pit.type.phase, 180)
        self.assertEqual(repr(pit), '<PiTorsion; (%r,%r)--%r--%r--(%r,%r); type=%r>' %
                tuple(atoms + [pit.type]))

    #=============================================

    def test_stretchbend(self):
        """ Tests the stretch-bend term and type """
        atoms = TrackedList()
        atoms.extend([Atom(list=atoms) for i in range(3)])
        bonds = [Bond(atoms[0], atoms[1]), Bond(atoms[1], atoms[2])]
        strbnd = StretchBend(atoms[0], atoms[1], atoms[2], StretchBendType(10.0, 11.0, 1.1, 1.2, 109.0))
        strbnds = TrackedList()
        strbnds.append(strbnd.type)
        strbnd.type.list = strbnds
        for obj in atoms + bonds:
            self.assertIn(obj, strbnd)
        self.assertEqual(strbnd.type.idx, 0)
        self.assertIs(strbnd.type.list, strbnds)
        self.assertEqual(strbnd.type.k1, 10)
        self.assertEqual(strbnd.type.uk1, 10*u.kilocalories_per_mole / (u.radians*u.angstroms))
        self.assertEqual(strbnd.type.k2, 11)
        self.assertEqual(strbnd.type.uk2, 11*u.kilocalories_per_mole / (u.radians*u.angstroms))
        self.assertEqual(strbnd.type.req1, 1.1)
        self.assertEqual(strbnd.type.ureq1, 1.1 * u.angstroms)
        self.assertEqual(strbnd.type.req2, 1.2)
        self.assertEqual(strbnd.type.ureq2, 1.2 * u.angstroms)
        self.assertEqual(strbnd.type.theteq, 109.0)
        self.assertEqual(strbnd.type.utheteq, 109.0 * u.degrees)
        # Test the StretchBendType.__copy__ method
        cp = copy(strbnd.type)
        self.assertIsNot(cp, strbnd.type)
        self.assertIs(cp.list, None)
        self.assertEqual(cp.idx, -1)
        self.assertEqual(cp.k1, 10)
        self.assertEqual(cp.k2, 11)
        self.assertEqual(cp.req1, 1.1)
        self.assertEqual(cp.req2, 1.2)
        self.assertEqual(cp.theteq, 109.0)
        self.assertEqual(repr(strbnd), '<StretchBend; %r--%r--%r; type=%r>' %
                tuple(atoms + [strbnd.type]))
        cp = copy(strbnd.type)
        self.assertIsNot(cp, strbnd.type)
        self.assertEqual(hash(cp), hash(strbnd.type))
        self.assertIn(cp, {strbnd.type})

    #=============================================

    def test_torsion_torsion(self):
        """ Tests the torsion-torsion term used in the AMOEBA force field """
        data = AmberFormat(get_fn('amoeba.parm7')).parm_data
        pre = 'AMOEBA_TORSION_TORSION_TORTOR_TABLE_'
        atoms = TrackedList()
        atoms.extend([Atom(list=atoms) for i in range(15)])
        bonds = [Bond(atoms[i], atoms[i+1]) for i in range(14)]
        tortor_types = TrackedList()
        tortor_types.extend([
                TorsionTorsionType((25, 25), data[pre+'01_ANGLE1'],
                    data[pre+'01_ANGLE2'],
                    data[pre+'01_FUNC'],
                    list=tortor_types),
                TorsionTorsionType((25, 25),
                    data[pre+'01_ANGLE1'], data[pre+'01_ANGLE2'],
                    data[pre+'01_FUNC'], data[pre+'01_DFUNC_DANGLE1'],
                    data[pre+'01_DFUNC_DANGLE2'],
                    data[pre+'01_D2FUNC_DANGLE1_DANGLE2'], list=tortor_types),
                TorsionTorsionType((25, 25),
                    data[pre+'03_ANGLE1'], data[pre+'03_ANGLE2'],
                    data[pre+'03_FUNC'], data[pre+'03_DFUNC_DANGLE1'],
                    data[pre+'03_DFUNC_DANGLE2'],
                    data[pre+'03_D2FUNC_DANGLE1_DANGLE2'], list=tortor_types),
        ])
        self.assertEqual(repr(tortor_types[0]), '<TorsionTorsionType; 25x25>')
        tortor1 = TorsionTorsion(atoms[0], atoms[1], atoms[2], atoms[3],
                                 atoms[4], type=tortor_types[0])
        tortor2 = TorsionTorsion(atoms[5], atoms[6], atoms[7], atoms[8],
                                 atoms[9], type=tortor_types[1])
        tortor3 = TorsionTorsion(atoms[10], atoms[11], atoms[12], atoms[13],
                                 atoms[14], type=tortor_types[2])
        # Check some error handling
        self.assertRaises(ValueError, lambda:
                TorsionTorsionType((25,),
                    data[pre+'03_ANGLE1'], data[pre+'03_ANGLE2'],
                    data[pre+'03_FUNC'], data[pre+'03_DFUNC_DANGLE1'],
                    data[pre+'03_DFUNC_DANGLE2'],
                    data[pre+'03_D2FUNC_DANGLE1_DANGLE2']
                )
        )
        self.assertRaises(ValueError, lambda:
                TorsionTorsionType((26, 25),
                    data[pre+'03_ANGLE1'], data[pre+'03_ANGLE2'],
                    data[pre+'03_FUNC'], data[pre+'03_DFUNC_DANGLE1'],
                    data[pre+'03_DFUNC_DANGLE2'],
                    data[pre+'03_D2FUNC_DANGLE1_DANGLE2'])
        )

        # Check the container properties of the torsion-torsion
        for i, atom in enumerate(atoms):
            if i < 5:
                self.assertIn(atom, tortor1)
                self.assertNotIn(atom, tortor2)
                self.assertNotIn(atom, tortor3)
            elif i < 10:
                self.assertNotIn(atom, tortor1)
                self.assertIn(atom, tortor2)
                self.assertNotIn(atom, tortor3)
            else:
                self.assertNotIn(atom, tortor1)
                self.assertNotIn(atom, tortor2)
                self.assertIn(atom, tortor3)

        omitted_bonds = 0
        for i, bond in enumerate(bonds):
            if bond.atom1.idx < 5 and bond.atom2.idx < 5:
                self.assertIn(bond, tortor1)
                self.assertNotIn(bond, tortor2)
                self.assertNotIn(bond, tortor3)
            elif 5 <= bond.atom1.idx < 10 and 5 <= bond.atom2.idx < 10:
                self.assertNotIn(bond, tortor1)
                self.assertIn(bond, tortor2)
                self.assertNotIn(bond, tortor3)
            elif 10 <= bond.atom1.idx < 15 and 10 <= bond.atom2.idx < 15:
                self.assertNotIn(bond, tortor1)
                self.assertNotIn(bond, tortor2)
                self.assertIn(bond, tortor3)
            else:
                omitted_bonds += 1
                self.assertNotIn(bond, tortor1)
                self.assertNotIn(bond, tortor2)
                self.assertNotIn(bond, tortor3)
        self.assertTrue(omitted_bonds > 0)
        self.assertEqual(tortor1.type.idx, 0)
        self.assertEqual(tortor2.type.idx, 1)
        self.assertEqual(tortor3.type.idx, 2)
        # Now check the _TorTorTable API
        self.assertRaises(TypeError, lambda:
                topologyobjects._TorTorTable([0, 90, 180], [0, 90, 180], [1, 2])
        )
        self.assertEqual(tortor1.type.f, tortor2.type.f)
        self.assertNotEqual(tortor1.type.f, tortor3.type.f)
        # Check the element accessors
        c = 0
        for a1 in data[pre+'01_ANGLE1']:
            for a2 in data[pre+'01_ANGLE2']:
                self.assertEqual(tortor1.type.f[(a1, a2)], data[pre+'01_FUNC'][c])
                self.assertEqual(tortor1.type.f[a1, a2], data[pre+'01_FUNC'][c])
                c += 1
        # Check the TorsionTorsionType.__copy__ method
        cp = copy(tortor_types[-1])
        self.assertIsNot(cp, tortor_types[-1])
        self.assertIsNot(cp.dfda1, tortor_types[-1].dfda1)
        self.assertIsNot(cp.dfda2, tortor_types[-1].dfda2)
        self.assertIsNot(cp.d2fda1da2, tortor_types[-1].d2fda1da2)
        self.assertIsNot(cp.dfda1.data, tortor_types[-1].dfda1.data)
        self.assertIsNot(cp.dfda2.data, tortor_types[-1].dfda2.data)
        self.assertIsNot(cp.d2fda1da2.data, tortor_types[-1].d2fda1da2.data)
        self.assertEqual(cp, tortor_types[-1])
        # Exception handling
        self.assertRaises(MoleculeError, lambda:
                TorsionTorsion(atoms[0], atoms[1], atoms[2], atoms[3], atoms[0])
        )
        self.assertRaises(MoleculeError, lambda: atoms[0].tortor_to(atoms[0]))
        # Make sure that if a torsion-torsion table has different angles, but
        # same energies, it is considered non-equal
        ttt = topologyobjects._TorTorTable(
                [x+0.001 for x in data[pre+'01_ANGLE1']],
                data[pre+'01_ANGLE2'], data[pre+'01_FUNC']
        )
        self.assertNotEqual(tortor1.type.f, ttt)
        self.assertNotEqual(ttt, tortor1.type.f)

    #=============================================

    def test_chiral_frame(self):
        """ Tests the chiral frame object used in the AMOEBA force field """
        atom1 = Atom()
        atom2 = Atom()
        cf1 = ChiralFrame(atom1, atom2, 1)
        cf2 = ChiralFrame(atom1, atom2, -1)
        self.assertRaises(ValueError, lambda: ChiralFrame(atom1, atom2, 2))
        self.assertIn(atom1, cf1)
        self.assertIn(atom2, cf1)
        self.assertIn(atom1, cf2)
        self.assertIn(atom2, cf2)
        self.assertEqual(repr(cf1), '<ChiralFrame; %r--%r, direction=1>' % (atom1, atom2))

    #=============================================

    def test_docstrings(self):
        """ Running the topologyobjects docstring examples/tests """
        import doctest
        results = doctest.testmod(topologyobjects)
        self.assertEqual(results.failed, 0)

    #=============================================

    def test_residue(self):
        """ Tests the Residue object """
        atoms = AtomList()
        atoms.extend([Atom(name='CA') for i in range(10)])
        res = Residue('ALA', segid='SYS')
        for atom in atoms:
            res.add_atom(atom)
        self.assertEqual(len(res), len(atoms))
        for atom in atoms:
            self.assertIs(atom.residue, res)
        for x, y in zip(atoms, res):
            self.assertIs(x, y)
        for atom in atoms:
            self.assertIn(atom, res)

        # Check segid assignment
        self.assertEqual(res.segid, 'SYS')

        # __repr__ testing
        self.assertEqual(repr(atoms[0]), '<Atom CA [0]; In ALA -1>')
        self.assertEqual(repr(res), '<Residue ALA[-1]; segid=SYS>')
        self.assertEqual(repr(res), '<Residue ALA[-1]; segid=SYS>')
        res.chain = 'B'
        res.insertion_code = 'B'
        self.assertEqual(repr(res),
                '<Residue ALA[-1]; chain=B; insertion_code=B; segid=SYS>')

        # Try deleting all of our atoms
        while len(res):
            lenres = len(res)
            atom = random.choice(atoms)
            is_inside = atom in res
            res.delete_atom(atom)
            self.assertEqual(lenres-len(res), is_inside)
            self.assertIs(atom.residue, None)
        self.assertEqual(len(res), 0)
        for atom in atoms:
            self.assertIs(atom.residue, None)
        self.assertTrue(res.is_empty())

    #=============================================

    def test_atom_list(self):
        """ Tests the AtomList class """
        atoms = AtomList()
        res = Residue('ALA')
        atoms.extend([Atom(name='CA') for i in range(15)])
        # Test __repr__ before adding residue
        self.assertEqual(repr(atoms[0]), '<Atom CA [0]>')
        for atom in atoms:
            res.add_atom(atom)
            self.assertIs(atom.list, atoms)
        # Test the indexing
        for i, atom in enumerate(atoms):
            self.assertEqual(atom.idx, i)
        atoms.changed = False
        i = 0
        # Test both pop and remove
        while atoms:
            if i % 2 == 0:
                atom = atoms.pop()
            else:
                atom = atoms[-1]
                atoms.remove(atom)
            self.assertEqual(atom.idx, -1)
            self.assertIs(atom.residue, None)
            self.assertIs(atom.list, None)
            self.assertTrue(atoms.changed)
            atoms.changed = False
            i += 1
        atoms.extend([Atom() for i in range(15)])
        self.assertTrue(res.is_empty())
        for i, atom in enumerate(atoms):
            self.assertEqual(atom.idx, i)
            self.assertIs(atom.residue, None)
        atoms.changed = False
        # Tests the insertion
        atoms.insert(2, Atom(name='TOK'))
        self.assertTrue(atoms.changed)
        self.assertEqual(len(atoms), 16)
        for i, atom in enumerate(atoms):
            self.assertEqual(atom.idx, i)
            self.assertIs(atom.list, atoms)
        self.assertEqual(atoms[2].name, 'TOK')
        atoms.changed = False
        # Tests appending
        atoms.append(Atom(name='TOK2'))
        self.assertTrue(atoms.changed)
        self.assertEqual(len(atoms), 17)
        for i, atom in enumerate(atoms):
            self.assertEqual(atom.idx, i)
        self.assertEqual(atoms[2].name, 'TOK')
        self.assertEqual(atoms[-1].name, 'TOK2')
        # Test atom comparisons
        self.assertGreater(atoms[1], atoms[0])
        self.assertLess(atoms[0], atoms[1])
        self.assertGreaterEqual(atoms[0], atoms[0])
        self.assertGreaterEqual(atoms[1], atoms[0])
        self.assertLessEqual(atoms[0], atoms[0])
        self.assertLessEqual(atoms[0], atoms[1])
        # Check some error handling
        self.assertRaises(ValueError, lambda: atoms.remove(atoms.pop()))
        self.assertRaises(RuntimeError, atoms.assign_nbidx_from_types)
        self.assertIs(atoms.__iadd__([Atom(), Atom()]), NotImplemented)
        self.assertRaises(IndexError, lambda: atoms.find_original_index(1000))

    #=============================================

    def test_residue_list(self):
        """ Tests the ResidueList class """
        atoms = AtomList()
        reslist = ResidueList()
        for i in range(40):
            new_at = Atom()
            if i < 10:
                reslist.add_atom(new_at, 'A', 10)
            elif i < 21:
                reslist.add_atom(new_at, 'A', 11)
            elif i < 33:
                reslist.add_atom(new_at, 'B', 20)
            else:
                reslist.add_atom(new_at, 'C', 21)
            atoms.append(new_at)
            self.assertEqual(atoms[-1].idx, i)

        self.assertEqual(len(reslist), 4)
        self.assertEqual(len(reslist[0]), 10)
        self.assertEqual(len(reslist[1]), 11)
        self.assertEqual(len(reslist[2]), 12)
        self.assertEqual(len(reslist[3]), 7)
        del atoms[:10]
        self.assertEqual(len(reslist), 4)
        self.assertTrue(reslist[0].is_empty())
        reslist.prune()
        self.assertEqual(len(reslist), 3) # Got rid of empty residue.
        self.assertEqual(len(atoms), sum([len(r) for r in reslist]))

        # Check ordering of residues
        self.assertGreater(reslist[1], reslist[0])
        self.assertFalse(reslist[1] > reslist[1])
        self.assertLess(reslist[0], reslist[1])
        self.assertFalse(reslist[0] < reslist[0])

        self.assertGreaterEqual(reslist[0], reslist[0])
        self.assertGreaterEqual(reslist[1], reslist[0])
        self.assertFalse(reslist[0] >= reslist[1])
        self.assertLessEqual(reslist[0], reslist[0])
        self.assertLessEqual(reslist[0], reslist[1])
        self.assertFalse(reslist[1] <= reslist[0])
        # Delete a couple residues
        del reslist[:2]
        self.assertEqual(len(reslist), 1)
        del reslist[0]
        self.assertEqual(len(reslist), 0)
        del reslist[10:20]
        self.assertEqual(len(reslist), 0)

    #=============================================

    def test_tracked_list(self):
        """ Test the TrackedList class """
        items = TrackedList()
        items.extend([Atom() for i in range(100)])
        for atom in items:
            self.assertEqual(atom.idx, -1)
            self.assertIs(atom.list, None)
        items.claim()
        for i, atom in enumerate(items):
            self.assertEqual(atom.idx, i)
        self.assertTrue(items.changed)
        self.assertFalse(items.needs_indexing)
        items.changed = False
        self.assertFalse(items.changed)
        all_atoms = items[:]
        for atom in all_atoms:
            self.assertIsNot(atom.list, all_atoms)
            self.assertIn(atom, items)
        # Test the repr
        self.assertEqual(repr(items), 'TrackedList([\n\t' +
                '\n\t'.join(repr(a) for a in items[:24]) + '\n\t...\n\t' +
                '\n\t'.join(repr(a) for a in items[-5:]) + '\n])')
        # Delete the whole list in random order
        while items:
            atom = items.pop(random.choice(range(len(items))))
            for item in items:
                self.assertIsNot(item, atom)
            self.assertEqual(atom.idx, -1)
            self.assertTrue(items.changed)
            items.changed = False
        # Now test the remove method
        self.assertFalse(items.changed)
        items.append(Atom())
        self.assertTrue(items.changed)
        items.changed = False
        items.remove(items[0])
        self.assertTrue(items.changed)
        # Delete slice
        while all_atoms:
            del all_atoms[:10]
        self.assertEqual(len(all_atoms), 0)
        # Try a tracked list of an immutable built-in type
        ints = TrackedList(range(9))
        ints.append(9)
        self.assertTrue(ints.needs_indexing)
        ints.claim()
        self.assertFalse(ints.needs_indexing)
        before_len = len(ints)
        ints.prune_unused() # Should be a no-op
        self.assertEqual(len(ints), before_len)
        self.assertEqual(repr(ints), """TrackedList([
\t0
\t1
\t2
\t3
\t4
\t5
\t6
\t7
\t8
\t9
])""")

    #=============================================

    def test_nonbonded_exception(self):
        """ Tests the NonbondedException and NonbondedExceptionType objects """
        a1 = Atom(name='CA', type='CX', atomic_number=6)
        a2 = Atom(name='CB', type='CT', atomic_number=6)
        a3 = Atom(name='DU', type='DU', atomic_number=0)
        nbe = NonbondedException(a1, a2)
        self.assertIn(a1, nbe)
        self.assertIn(a2, nbe)
        self.assertNotIn(a3, nbe)
        self.assertEqual(repr(nbe), '<NonbondedException; %r and %r>' % (a1, a2))
        # Now add the type
        nbet = NonbondedExceptionType(1.0, 2.0, chgscale=3.0)
        nbet2 = NonbondedExceptionType(1.1, 2.0, chgscale=3.0)
        self.assertEqual(nbet.sigma, 2**(-1/6))
        self.assertEqual(nbet.usigma, 2**(-1/6) * u.angstroms)
        self.assertEqual(nbet.urmin, 1.0 * u.angstroms)
        self.assertEqual(nbet.uepsilon, 2.0 * u.kilocalories_per_mole)
        self.assertEqual(nbet, nbet)
        self.assertNotEqual(nbet, nbet2)
        nbet2.sigma = 1.0
        self.assertEqual(nbet2.rmin*2**(-1/6), nbet2.sigma)
        nbe.type = nbet
        self.assertEqual(repr(nbe), '<NonbondedException; %r and %r, type=%r>' % (a1, a2, nbet))
        cp = copy(nbet)
        self.assertIsNot(cp, nbet)
        self.assertEqual(hash(cp), hash(nbet))
        self.assertIn(cp, {nbet})

    #=============================================

    def test_acceptor_donor(self):
        """ Test the AcceptorDonor API """
        a1 = Atom(name='O', atomic_number=8)
        a2 = Atom(name='H', atomic_number=1)
        a3 = Atom(name='H2', atomic_number=1)
        ad = AcceptorDonor(a1, a2)
        self.assertIn(a1, ad)
        self.assertIn(a2, ad)
        self.assertNotIn(a3, ad)
        self.assertEqual(repr(ad), '<AcceptorDonor; %r %r>' % (a1, a2))

    #=============================================

    def test_helper_functions(self):
        """ Test the helper private functions in topologyobjects.py """
        self.assertEqual(_strip_units(10), 10)
        self.assertAlmostEqual(_strip_units(math.pi*u.radians), 180)
        self.assertAlmostEqual(_strip_units(1*u.kilojoule_per_mole), 1/4.184)
        self.assertEqual(_strip_units(10*u.kilojoule_per_mole, u.kilocalorie_per_mole),
                         10/4.184)
        self.assertAlmostEqual(_strip_units(180*u.degree, u.radians), math.pi)
