""" Tests the functionality in parmed/parameters.py """
from __future__ import division, print_function

import parmed as pmd
import unittest
import utils

class TestParameterSet(unittest.TestCase):
    """ Tests ParameterSet """

    def test_duplicate_bond_type(self):
        """ Tests handling of duplicate bond type in ParameterSet """
        struct = pmd.Structure()
        struct.add_atom(pmd.Atom('CA', type='CX'), 'ALA', 1)
        struct.add_atom(pmd.Atom('CB', type='CY'), 'ALA', 1)
        struct.add_atom(pmd.Atom('CD', type='CX'), 'GLY', 2)
        struct.add_atom(pmd.Atom('CE', type='CY'), 'GLY', 2)
        struct.bond_types.append(pmd.BondType(10.0, 1.0))
        struct.bond_types.append(pmd.BondType(11.0, 1.1))
        struct.bond_types.claim()
        struct.bonds.append(
                pmd.Bond(struct[0], struct[1], type=struct.bond_types[0])
        )
        struct.bonds.append(
                pmd.Bond(struct[2], struct[3], type=struct.bond_types[1])
        )
        self.assertRaises(pmd.exceptions.ParameterError, lambda:
                pmd.ParameterSet.from_structure(struct,
                            allow_unequal_duplicates=False)
        )

    def test_duplicate_angle_type(self):
        """ Tests handling of duplicate angle type in ParameterSet """
        struct = pmd.Structure()
        struct.add_atom(pmd.Atom('CA', type='CX'), 'ALA', 1)
        struct.add_atom(pmd.Atom('CB', type='CY'), 'ALA', 1)
        struct.add_atom(pmd.Atom('CC', type='CZ'), 'ALA', 1)
        struct.add_atom(pmd.Atom('CD', type='CX'), 'GLY', 2)
        struct.add_atom(pmd.Atom('CE', type='CY'), 'GLY', 2)
        struct.add_atom(pmd.Atom('CF', type='CZ'), 'GLY', 2)
        struct.angle_types.append(pmd.AngleType(10.0, 120.0))
        struct.angle_types.append(pmd.AngleType(11.0, 109.0))
        struct.angle_types.claim()
        struct.angles.append(
                pmd.Angle(struct[0], struct[1], struct[2],
                          type=struct.angle_types[0])
        )
        struct.angles.append(
                pmd.Angle(struct[3], struct[4], struct[5],
                          type=struct.angle_types[1])
        )
        self.assertRaises(pmd.exceptions.ParameterError, lambda:
                pmd.ParameterSet.from_structure(struct,
                            allow_unequal_duplicates=False)
        )
