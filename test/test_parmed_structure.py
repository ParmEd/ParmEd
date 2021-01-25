"""
Tests the parmed/structure module
"""
from __future__ import division

import bz2
import gzip
from copy import copy
import numpy as np
import os
import parmed as pmd
from parmed.exceptions import CharmmWarning, ParameterWarning
import parmed.structure as structure
from parmed.topologyobjects import *
import parmed.unit as u
from parmed.utils.six.moves import StringIO
from parmed.utils.six import integer_types
from parmed.utils.six.moves import range, zip
from parmed.utils import PYPY
import random
import string
import unittest
from utils import create_random_structure, get_fn, FileIOTestCase, has_openmm, app, HAS_GROMACS
import warnings

class TestStructureAPI(unittest.TestCase):
    """ Tests the underlying Structure API """

    def setUp(self):
        s = self.s = structure.Structure()
        s.add_atom(Atom(atomic_number=6), 'ALA', 1, 'A')
        s.add_atom(Atom(atomic_number=1), 'ALA', 1, 'A')
        s.add_atom(Atom(atomic_number=7), 'ALA', 1, 'A')
        s.add_atom(Atom(atomic_number=8), 'ALA', 1, 'A')
        s.add_atom(Atom(atomic_number=16), 'GLY', 2, 'A')
        s.add_atom(Atom(atomic_number=99), 'GLY', 3, 'A')
        s.add_atom(Atom(atomic_number=6), 'GLY', 3, 'B')
        s.add_atom(Atom(atomic_number=6), 'GLY', 3, 'B')
        s.add_atom(Atom(atomic_number=6), 'GLY', 3, 'B')

    @unittest.skipUnless(has_openmm, 'Cannot test without OpenMM')
    def test_gb_assignment(self):
        """ Tests GB parameter assignment """
        # GBneck1
        params = self.s._get_gb_parameters(app.GBn)
        # Radii are assigned elsewhere -- check screen
        self.assertEqual(params[0][1], 0.48435382330)
        self.assertEqual(params[1][1], 1.09085413633)
        self.assertEqual(params[2][1], 0.700147318409)
        self.assertEqual(params[3][1], 1.06557401132)
        self.assertEqual(params[4][1], 0.602256336067)
        self.assertEqual(params[5][1], 0.5)
        # GBneck2
        self.s.bonds.append(Bond(self.s[0], self.s[1]))
        params = self.s._get_gb_parameters(app.GBn2)
        self.assertEqual(params[0][1], 1.058554)
        self.assertEqual(params[0][2], 0.733756)
        self.assertEqual(params[0][3], 0.506378)
        self.assertEqual(params[0][4], 0.205844)
        self.assertEqual(params[1][1], 1.425952)
        self.assertEqual(params[1][2], 0.788440)
        self.assertEqual(params[1][3], 0.798699)
        self.assertEqual(params[1][4], 0.437334)
        self.assertEqual(params[2][1], 0.733599)
        self.assertEqual(params[2][2], 0.503364)
        self.assertEqual(params[2][3], 0.316828)
        self.assertEqual(params[2][4], 0.192915)
        self.assertEqual(params[3][1], 1.061039)
        self.assertEqual(params[3][2], 0.867814)
        self.assertEqual(params[3][3], 0.876635)
        self.assertEqual(params[3][4], 0.387882)
        self.assertEqual(params[4][1], -0.703469)
        self.assertEqual(params[4][2], 0.867814)
        self.assertEqual(params[4][3], 0.876635)
        self.assertEqual(params[4][4], 0.387882)
        self.assertEqual(params[5][1], 0.5)
        self.assertEqual(params[5][2], 1.0)
        self.assertEqual(params[5][3], 0.8)
        self.assertEqual(params[5][4], 4.85)

    def test_combining_rule(self):
        """ Tests the Structure.combining_rule attribute """
        self.assertEqual(self.s.combining_rule, 'lorentz')
        self.s.combining_rule = 'geometric'
        self.assertEqual(self.s.combining_rule, 'geometric')
        def fail():
            self.s.combining_rule = 'badchoice'
        self.assertRaises(ValueError, fail)

    def test_add_atom(self):
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

    def test_split_duplicate_detection(self):
        """ Tests detection of duplicate molecules in split """
        parm = pmd.load_file(get_fn('ash.parm7')) * 5
        self.assertEqual(len(parm.split()), 1)
        self.assertEqual(len(parm.split()[0][1]), 5)
        parm.residues[0].name = 'ACE1'
        x = parm.split()
        self.assertEqual(len(x), 2)
        self.assertEqual(len(x[0][1]), 1)
        self.assertEqual(len(x[1][1]), 4)
        parm.residues[0].name = 'ACE'
        self.assertEqual(len(parm.split()), 1)
        parm[0].name += '1'
        x = parm.split()
        self.assertEqual(len(x), 2)
        self.assertEqual(len(x[0][1]), 1)
        self.assertEqual(len(x[1][1]), 4)

    def test_split_single_atom(self):
        """ Tests splitting of a structure with single atom residues """
        s = structure.Structure()
        na = pmd.periodic_table.AtomicNum['Na']
        cl = pmd.periodic_table.AtomicNum['Cl']
        li = pmd.periodic_table.AtomicNum['Li']
        f = pmd.periodic_table.AtomicNum['F']
        s.add_atom(Atom(name='NA', type='Na+', atomic_number=na), 'NA', 1, 'A')
        s.add_atom(Atom(name='NA', type='Na+', atomic_number=na), 'NA', 2, 'A')
        s.add_atom(Atom(name='LI', type='Li+', atomic_number=li), 'LI', 3, 'A')
        s.add_atom(Atom(name='LI', type='Li+', atomic_number=li), 'LI', 4, 'A')
        s.add_atom(Atom(name='CL', type='Cl-', atomic_number=cl), 'CL', 5, 'A')
        s.add_atom(Atom(name='CL', type='Cl-', atomic_number=cl), 'CL', 6, 'A')
        s.add_atom(Atom(name='F', type='F-', atomic_number=f), 'F', 7, 'B')
        s.add_atom(Atom(name='F', type='F-', atomic_number=f), 'F', 8, 'B')
        x = s.split()
        self.assertEqual(len(x), 4)
        self.assertEqual(set(len(_[1]) for _ in x), {2})

    def test_bool(self):
        """ Tests bool-ness of Structure """
        self.assertTrue(bool(self.s))
        self.assertFalse(bool(structure.Structure()))

    def test_add_atom_to_residue(self):
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

        # Error handling
        self.assertRaises(ValueError, lambda:
                s.add_atom_to_residue(Atom(name='AXE'), Residue('ALA')))

    def test_box_handling(self):
        """ Tests that Structure.box is always the expected type """
        s = create_random_structure(parametrized=False)
        self.assertIs(s.box, None)
        self.assertIs(s.get_box(), None)
        self.assertRaises(IndexError, lambda: s.get_box(0))
        s.box = [10, 10, 10, 90, 90, 90]
        self.assertIsInstance(s.box, np.ndarray)
        self.assertEqual(s.box[0], 10)
        self.assertEqual(s.box[1], 10)
        self.assertEqual(s.box[2], 10)
        self.assertEqual(s.box[3], 90)
        self.assertEqual(s.box[4], 90)
        self.assertEqual(s.box[5], 90)
        s.box = np.array([10, 10, 10, 90, 90, 90], dtype=np.float64)
        self.assertIsInstance(s.box, np.ndarray)
        self.assertEqual(s.box[0], 10)
        self.assertEqual(s.box[1], 10)
        self.assertEqual(s.box[2], 10)
        self.assertEqual(s.box[3], 90)
        self.assertEqual(s.box[4], 90)
        self.assertEqual(s.box[5], 90)
        s.box = [10*u.angstroms, 10*u.angstroms, 10*u.angstroms,
                 90*u.degrees, 90*u.degrees, 90*u.degrees]
        self.assertIsInstance(s.box, np.ndarray)
        self.assertEqual(s.box[0], 10)
        self.assertEqual(s.box[1], 10)
        self.assertEqual(s.box[2], 10)
        self.assertEqual(s.box[3], 90)
        self.assertEqual(s.box[4], 90)
        self.assertEqual(s.box[5], 90)
        s.box = [[10*u.angstroms, 1*u.nanometers, 10*u.angstroms,
                  90*u.degrees, 90*u.degrees, 90*u.degrees],
                 [11*u.angstroms, 1*u.nanometers, 11*u.angstroms,
                  91*u.degrees, 91*u.degrees, 91*u.degrees]]
        self.assertIsInstance(s.box, np.ndarray)
        self.assertEqual(s.box[0], 10)
        self.assertEqual(s.box[1], 10)
        self.assertEqual(s.box[2], 10)
        self.assertEqual(s.box[3], 90)
        self.assertEqual(s.box[4], 90)
        self.assertEqual(s.box[5], 90)
        box = s.get_box(0)
        self.assertEqual(box[0], 10)
        self.assertEqual(box[1], 10)
        self.assertEqual(box[2], 10)
        self.assertEqual(box[3], 90)
        self.assertEqual(box[4], 90)
        self.assertEqual(box[5], 90)
        box = s.get_box(1)
        self.assertEqual(box[0], 11)
        self.assertEqual(box[1], 10)
        self.assertEqual(box[2], 11)
        self.assertEqual(box[3], 91)
        self.assertEqual(box[4], 91)
        self.assertEqual(box[5], 91)
        box = s.get_box('all')
        self.assertEqual(box[0][0], 10)
        self.assertEqual(box[0][1], 10)
        self.assertEqual(box[0][2], 10)
        self.assertEqual(box[0][3], 90)
        self.assertEqual(box[0][4], 90)
        self.assertEqual(box[0][5], 90)
        self.assertEqual(box[1][0], 11)
        self.assertEqual(box[1][1], 10)
        self.assertEqual(box[1][2], 11)
        self.assertEqual(box[1][3], 91)
        self.assertEqual(box[1][4], 91)
        self.assertEqual(box[1][5], 91)
        self.assertRaises(IndexError, lambda: s.get_box(3))
        # Now test box vectors
        s.box_vectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] * u.nanometers
        np.testing.assert_equal(s.box, [10, 10, 10, 90, 90, 90])

    def test_bad_box_handling(self):
        """ Tests error handling when Structure.box is improperly assigned """
        s = create_random_structure(parametrized=True)
        def wrong_number_of_args():
            s.box = [0, 1, 2, 3, 4]
        def wrong_number_of_args2():
            s.box = [0, 1, 2, 3, 4, 5, 6]
        self.assertRaises(ValueError, wrong_number_of_args)
        self.assertRaises(ValueError, wrong_number_of_args2)
        try:
            wrong_number_of_args()
        except ValueError as err:
            self.assertIn('6', str(err))
        # Try wrong units
        for i in range(6):
            box = [10, 10, 10, 90, 90, 90]
            box[i] *= u.liters
            def func():
                s.box = box
            self.assertRaises(TypeError, func)

    def test_coordinates(self):
        """ Tests coordinate handling in Structure """
        s = create_random_structure(parametrized=False)
        self.assertIs(s.coordinates, None)
        self.assertIs(s.positions, None)
        natom = len(s.atoms)
        # Make sure coordinates will be generated from Atom.xx, xy, xz if not
        # otherwise generated
        xyz = np.random.random((natom, 3))
        for a, x in zip(s.atoms, xyz):
            a.xx, a.xy, a.xz = x
        self.assertEqual(s.coordinates.shape, (natom, 3))
        np.testing.assert_equal(s.coordinates, xyz[:,:])
        self.assertIs(s._coordinates, None)
        # Now set multiple frames
        xyz = np.random.random((5, natom, 3)).tolist()
        s.coordinates = xyz
        self.assertIsInstance(s.coordinates, np.ndarray)
        self.assertEqual(s.coordinates.shape, (natom, 3))
        for a, x in zip(s.atoms, xyz[0]):
            self.assertEqual(a.xx, x[0])
            self.assertEqual(a.xy, x[1])
            self.assertEqual(a.xz, x[2])
        np.testing.assert_equal(s.get_coordinates('all'), xyz)
        for i in range(5):
            np.testing.assert_equal(s.get_coordinates(i), xyz[i])
        # Now try setting with units
        xyz = u.Quantity(np.random.random((3, natom, 3)), u.nanometers)
        s.coordinates = xyz
        self.assertIsInstance(s.coordinates, np.ndarray)
        self.assertEqual(s.coordinates.shape, (natom, 3))
        for a, x in zip(s.atoms, xyz[0]._value):
            self.assertEqual(a.xx, x[0]*10)
            self.assertEqual(a.xy, x[1]*10)
            self.assertEqual(a.xz, x[2]*10)
        # Now check setting None
        s.coordinates = None
        for a in s.atoms:
            self.assertFalse(hasattr(a, 'xx'))
            self.assertFalse(hasattr(a, 'xy'))
            self.assertFalse(hasattr(a, 'xz'))
        self.assertIs(s.coordinates, None)
        # Now check setting flattened arrays
        s.coordinates = np.random.random((natom, 3))
        self.assertEqual(s.coordinates.shape, (natom, 3))
        s.coordinates = np.random.random(natom*3)
        self.assertEqual(s.coordinates.shape, (natom, 3))
        s.coordinates = np.random.random(natom*3*10)
        self.assertEqual(s.coordinates.shape, (natom, 3))
        # Now check other iterables
        old_crds = s.coordinates
        s.coordinates = (random.random() for i in range(3*len(s.atoms)))
        self.assertEqual(s.coordinates.shape, (natom, 3))
        diff = (old_crds - s.coordinates).ravel()**2
        self.assertGreater(diff.sum(), 0.01)
        # Check setting the positions attribute
        new_crds = np.random.rand(len(s.atoms), 3) * u.nanometers
        s.positions = new_crds
        np.testing.assert_allclose(new_crds*10, s.coordinates)
        new_crds2 = np.random.rand(len(s.atoms)*3) * u.nanometers
        s.positions = new_crds2
        np.testing.assert_allclose(new_crds2*10,
                s.coordinates.reshape(new_crds2.shape))
        def bad():
            s.positions = [1, 2, 3]
        self.assertRaises(ValueError, bad)
        # Try deleting the coordinates for the first atom
        del s[0].xx, s[0].xy, s[0].xz
        self.assertIs(s.coordinates, None)
        # Make sure it doesn't affect the rest
        for atom in s.view[1:].atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))
        # Now make sure that setting coordinates to an empty iterable wipes out
        # all xx/y/z attributes
        s.coordinates = []
        self.assertIs(s.coordinates, None)
        # Restore coordinates
        s.coordinates = np.random.rand(3, len(s.atoms), 3)
        del s[0].xx
        self.assertIs(s.get_coordinates(), None)
        s.coordinates = xyz = np.random.rand(3, len(s.atoms), 3)
        del s.atoms[-1]
        np.testing.assert_equal(s.get_coordinates().flatten(), xyz[0,:-1,:].flatten())
        s.coordinates = np.random.rand(3, len(s.atoms), 3)
        s[0].xx = 10
        self.assertEqual(s.get_coordinates().shape, (1, len(s.atoms), 3))
        s.coordinates = xyz = np.random.rand(3, len(s.atoms), 3)
        s[0].xx = 10
        np.testing.assert_equal(s.get_coordinates(0)[1:,:], xyz[0,1:,:])
        s.coordinates = xyz = np.random.rand(3, len(s.atoms), 3)
        s[0].xx = 10
        self.assertRaises(IndexError, lambda: s.get_coordinates(1))

    def test_velocities(self):
        """ Tests coordinate handling in Structure """
        s = create_random_structure(parametrized=False)
        self.assertIs(s.velocities, None)
        x = np.random.rand(len(s.atoms), 3)
        s.velocities = x * u.nanometers/u.picosecond
        np.testing.assert_allclose(x*10, s.velocities)

    def test_coordinate_set_to_empty_list(self):
        """ Tests behavior of setting coordinates to an empty iterable """
        s = create_random_structure(parametrized=True)
        xyz = np.random.random((len(s.atoms), 3))
        s.coordinates = xyz
        # Make sure that the x, y, and z attributes of atoms are equal to the
        # given coordinates
        for atom, pos in zip(s.atoms, xyz):
            self.assertEqual(atom.xx, pos[0])
            self.assertEqual(atom.xy, pos[1])
            self.assertEqual(atom.xz, pos[2])
        # Now set coordinates to an empty list
        s.coordinates = []
        self.assertIs(s.coordinates, None)
        for atom in s.atoms:
            self.assertFalse(hasattr(atom, 'xx'))
            self.assertFalse(hasattr(atom, 'xy'))
            self.assertFalse(hasattr(atom, 'xz'))

    def test_strip(self):
        """ Tests the Structure.strip method """
        s = create_random_structure(parametrized=True)
        nres = len(s.residues)
        natom_per_res = [len(res) for res in s.residues]
        s.strip(':1-5')
        self.assertEqual(len(s.atoms), sum(natom_per_res) - sum(natom_per_res[:5]))
        # Bad usage of strip
        mask = pmd.amber.AmberMask(create_random_structure(parametrized=True), ':1')
        self.assertRaises(TypeError, lambda: s.strip(mask))
        self.assertRaises(TypeError, lambda: s.strip(10))

    def test_strip_array(self):
        """ Tests Structure.strip with a mask array """
        s = create_random_structure(parametrized=True)
        natom = len(s.atoms)
        maskarray = [random.choice([0, 1]) for i in range(natom)]
        s.strip(maskarray)
        self.assertEqual(len(s.atoms), natom-sum(maskarray))
        # Error checking
        self.assertRaises(ValueError, lambda: s.strip([1, 0]))

    def test_strip_with_coords(self):
        """ Tests the Structure.strip method when it has coordinates """
        s = create_random_structure(parametrized=True)
        coords = np.random.rand(10, len(s.atoms), 3)
        s.coordinates = coords
        nres = len(s.residues)
        natom_per_res = [len(res) for res in s.residues]
        s.strip(':1-5')
        self.assertEqual(len(s.atoms), sum(natom_per_res) - sum(natom_per_res[:5]))
        self.assertEqual(s.get_coordinates().shape, (10, len(s.atoms), 3))
        # No longer equal -- it's truncated
        self.assertNotEqual(coords.shape, (10, len(s.atoms), 3))
        # Check that the coordinates left are the same as the original
        # coordinates corresponding to the ones that did *not* get stripped
        n = sum(natom_per_res[:5])
        self.assertTrue((coords[:,n:,:] == s.get_coordinates()).all())

    def test_helpers(self):
        """ Test private helper functions in parmed/structure.py """
        ang, deg = u.angstroms, u.degrees
        strip_units = structure._strip_box_units
        self.assertEqual(strip_units([1, 2, 3]), [1, 2, 3])
        self.assertEqual(strip_units([1*ang, 2*ang, 3*ang]), [1, 2, 3])
        self.assertEqual(strip_units([1*ang, 90*deg]), [1, 90])
        self.assertEqual(strip_units([[1*ang, 90*deg], [2*ang, 109*deg]]),
                         [[1, 90], [2, 109]])
        self.assertRaises(TypeError, lambda: strip_units(['this', 'is', 'bad']))

    def test_radii_assignment(self):
        """ Test radius assignment in Structure """
        parm = pmd.load_file(get_fn('trx.prmtop'))
        criteria_satisfied = [0, 0, 0, 0, 0, 0]
        for atom in parm.atoms:
            if atom.atomic_number == 1:
                oa = atom.bond_partners[0]
                an, mass = oa.atomic_number, oa.mass
                oa.atomic_number = 8
                oa.mass = pmd.periodic_table.Mass['O']
                self.assertEqual(structure._mbondi(atom), 0.8)
                oa.atomic_number = 2
                oa.mass = pmd.periodic_table.Mass['He']
                self.assertEqual(structure._mbondi(atom), 1.2)
                oa.atomic_number = an
                oa.mass = mass
                criteria_satisfied[0] = 1
            if atom.name.startswith('OD') and atom.residue.name == 'ASP':
                self.assertEqual(structure._mbondi3(atom), 1.4)
                criteria_satisfied[1] = 1
            if atom.name.startswith('OXT'):
                self.assertEqual(structure._mbondi3(atom), 1.4)
                criteria_satisfied[2] = 1
            if atom.name.startswith('OE') and atom.residue.name == 'GLU':
                self.assertEqual(structure._mbondi3(atom), 1.4)
                criteria_satisfied[3] = 1
            if atom.name.startswith('HH') and atom.residue.name == 'ARG':
                self.assertEqual(structure._mbondi3(atom), 1.17)
                criteria_satisfied[4] = 1
            if atom.name.startswith('HE') and atom.residue.name == 'ARG':
                self.assertEqual(structure._mbondi3(atom), 1.17)
                criteria_satisfied[5] = 1
        # Make sure we ran all checks
        self.assertEqual(set(criteria_satisfied), {1})

    @unittest.skipUnless(has_openmm, 'Cannot test without OpenMM')
    def test_screen_assignment(self):
        """ Testing screen parameter assignment """
        self.assertEqual(
                structure._gb_rad_screen(Atom(atomic_number=9), app.GBn2)[1],
                0.88
        )
        self.assertEqual(
                structure._gb_rad_screen(Atom(atomic_number=15), app.GBn2)[1],
                0.86
        )
        self.assertEqual(
                structure._gb_rad_screen(Atom(atomic_number=16), app.GBn2)[1],
                0.96
        )
        self.assertEqual(
                structure._gb_rad_screen(Atom(atomic_number=2), app.GBn2)[1],
                0.8
        )

    @unittest.skipUnless(has_openmm, 'Cannot test without OpenMM')
    def test_omm_topology(self):
        """ Tests the creation of the OpenMM Topology instance """
        # Make sure an empty structure returns an empty topology
        self.assertEqual(len(list(pmd.Structure().topology.atoms())), 0)

    def test_join_dihedrals(self):
        """ Tests the Structure.join_dihedrals method """
        s = self.s
        dt1 = DihedralType(10, 2, 0)
        dt2 = DihedralType(12, 1, 180)
        dt3 = DihedralType(13, 3, 0)
        s.dihedrals.append(Dihedral(s[0], s[1], s[2], s[3], type=dt1))
        s.dihedrals.append(Dihedral(s[3], s[2], s[1], s[0], type=dt2))
        s.dihedrals.append(Dihedral(s[0], s[1], s[2], s[3], type=dt3))
        s.dihedral_types.extend([dt1, dt2, dt3])
        s.dihedral_types.claim()
        # Now join dihedrals
        s.join_dihedrals()
        self.assertEqual(len(s.dihedrals), 1)
        self.assertEqual(len(s.dihedral_types), 1)
        self.assertIsInstance(s.dihedral_types[0], pmd.DihedralTypeList)
        self.assertEqual(len(s[0].dihedral_partners), 3)
        self.assertEqual(len(s[0].dihedrals), 1)

    def test_prune_empty(self):
        """ Tests the prune_empty_terms function """
        s = create_random_structure(parametrized=True)
        # Bonds
        nterms = len(s.bonds)
        s.bonds[-1].delete()
        s.prune_empty_terms()
        self.assertEqual(len(s.bonds), nterms-1)
        # Angles
        nterms = len(s.angles)
        s.angles[-1].delete()
        s.prune_empty_terms()
        self.assertEqual(len(s.angles), nterms-1)
        # Dihedrals
        nterms = len(s.dihedrals)
        s.dihedrals[-1].delete()
        s.prune_empty_terms()
        self.assertEqual(len(s.dihedrals), nterms-1)
        # R-B torsions
        nterms = len(s.rb_torsions)
        s.rb_torsions[-1].delete()
        s.prune_empty_terms()
        self.assertEqual(len(s.rb_torsions), nterms-1)
        # Urey-Bradleys
        nterms = len(s.urey_bradleys)
        s.urey_bradleys[-1].delete()
        s.prune_empty_terms()
        self.assertEqual(len(s.urey_bradleys), nterms-1)
        # Impropers
        nterms = len(s.impropers)
        s.impropers[-1].delete()
        s.prune_empty_terms()
        self.assertEqual(len(s.impropers), nterms-1)
        # CMAPs
        nterms = len(s.cmaps)
        s.cmaps[-1].delete()
        s.prune_empty_terms()
        self.assertEqual(len(s.cmaps), nterms-1)
        # Trigonal angles
        nterms = len(s.trigonal_angles)
        s.trigonal_angles[-1].atom1 = s.trigonal_angles[-1].atom2 = \
                s.trigonal_angles[-1].atom3 = s.trigonal_angles[-1].atom4 = None
        s.prune_empty_terms()
        self.assertEqual(len(s.trigonal_angles), nterms-1)
        # OOP bends
        nterms = len(s.out_of_plane_bends)
        s.out_of_plane_bends[-1].atom1 = s.out_of_plane_bends[-1].atom2 = \
                s.out_of_plane_bends[-1].atom3 = s.out_of_plane_bends[-1].atom4 = None
        s.prune_empty_terms()
        self.assertEqual(len(s.out_of_plane_bends), nterms-1)
        # Pi-torsions
        nterms = len(s.pi_torsions)
        s.pi_torsions[-1].atom1 = s.pi_torsions[-1].atom2 = \
                s.pi_torsions[-1].atom3 = s.pi_torsions[-1].atom4 = \
                s.pi_torsions[-1].atom5 = s.pi_torsions[-1].atom6 = None
        s.prune_empty_terms()
        self.assertEqual(len(s.pi_torsions), nterms-1)
        # Stretch-bends
        nterms = len(s.stretch_bends)
        s.stretch_bends[-1].atom1 = s.stretch_bends[-1].atom2 = \
                s.stretch_bends[-1].atom3 = None
        s.prune_empty_terms()
        self.assertEqual(len(s.stretch_bends), nterms-1)
        # Torsion-torsions
        nterms = len(s.torsion_torsions)
        s.torsion_torsions[-1].delete()
        s.prune_empty_terms()
        self.assertEqual(len(s.torsion_torsions), nterms-1)
        # Chiral frames
        nterms = len(s.chiral_frames)
        s.chiral_frames[-1].atom1 = s.chiral_frames[-1].atom2 = None
        s.prune_empty_terms()
        self.assertEqual(len(s.chiral_frames), nterms-1)
        # Adjusts
        nterms = len(s.adjusts)
        s.adjusts[-1].atom1 = s.adjusts[-1].atom2 = None
        s.prune_empty_terms()
        self.assertEqual(len(s.adjusts), nterms-1)

class TestStructureAdd(unittest.TestCase):
    """ Tests the addition property of a System """

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
            self.assertEqual(a1.solvent_radius, a2.solvent_radius)
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
            self.assertEqual(a1.solvent_radius, a2.solvent_radius)
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

    def test_add_parametrized(self):
        """ Tests addition of two parametrized Structure instances """
        s1 = create_random_structure(parametrized=True)
        s2 = create_random_structure(parametrized=True)
        s1.coordinates = np.random.random((4, len(s1.atoms), 3))
        s2.coordinates = np.random.random((2, len(s2.atoms), 3))
        s1.box = [[10, 10, 10, 90, 90, 90], [11, 11, 11, 90, 90, 90]]
        s2.box = [[20, 20, 20, 90, 90, 90], [21, 21, 21, 90, 90, 90]]
        self.assertTrue(bool(s1.bond_types))
        self.assertTrue(bool(s2.bond_types))
        s = s1 + s2
        self._check_sum(s, s1, s2)
        self.assertEqual(s.get_coordinates('all').shape, (2, len(s.atoms), 3))

    def test_add_to_empty_structure(self):
        """ Tests addition to empty Structure """
        s1 = create_random_structure(parametrized=True)
        s2 = structure.Structure()
        s2 += s1
        self._check_sum(s2, structure.Structure(), s1)

    def test_iadd(self):
        """ Tests in-place addition of two Structure instances """
        s1 = create_random_structure(parametrized=True)
        s2 = create_random_structure(parametrized=True)
        s1cp = copy(s1)
        s1 += s2
        self._check_sum(s1, s1cp, s2)

    def test_add_not_parametrized(self):
        """ Tests addition of two non-parametrized Structure instances """
        s1 = create_random_structure(parametrized=False)
        s2 = create_random_structure(parametrized=False)
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
            self.assertEqual(a1.solvent_radius, a2.solvent_radius)
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

    def test_add_no_valence(self):
        """ Tests addition of two minimal Structure instances """
        s1 = create_random_structure(parametrized=False, novalence=True)
        s2 = create_random_structure(parametrized=False, novalence=True)
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
            self.assertEqual(a1.solvent_radius, a2.solvent_radius)
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

    def test_multiply_parametrized(self):
        """ Tests replicating a parametrized Structure instance """
        s1 = create_random_structure(parametrized=True)
        multfac = random.randint(2, 5)
        s2 = s1 * multfac
        self.assertIsNot(s1, s2)
        self.assertEqual(len(s2.atoms), len(s1.atoms) * multfac)
        self._check_mult(s2, s1, multfac)
        # Check commutativity
        s3 = multfac * s1
        self.assertIsNot(s1, s3)
        self.assertEqual(len(s3.atoms), len(s1.atoms) * multfac)
        self._check_mult(s3, s1, multfac)
        # Check when coordinates exist
        crd = np.random.rand(5, len(s1.atoms), 3)
        s1.coordinates = crd
        np.testing.assert_equal((s1*multfac).get_coordinates(),
                                np.hstack([crd for i in range(multfac)]))

    def test_multiply_not_parametrized(self):
        """ Tests replicating a non-parametrized Structure instance """
        s1 = create_random_structure(parametrized=False)
        multfac = random.randint(2, 5)
        s2 = s1 * multfac
        self.assertIsNot(s1, s2)
        self.assertEqual(len(s2.atoms), len(s1.atoms) * multfac)
        self._check_mult(s2, s1, multfac)
        # Check commutativity
        s3 = multfac * s1
        self.assertIsNot(s1, s3)
        self.assertEqual(len(s3.atoms), len(s1.atoms) * multfac)
        self._check_mult(s3, s1, multfac)

    def test_mult_no_valence(self):
        """ Tests addition of two minimal Structure instances """
        s1 = create_random_structure(parametrized=False, novalence=True)
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
            self.assertEqual(a1.solvent_radius, a2.solvent_radius)
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

class TestStructureSave(FileIOTestCase):
    """ Tests the universal "save" function in Structure """

    def setUp(self):
        self.sys1 = pmd.load_file(get_fn('ala3_solv.psf'))
        self.sys1.coordinates = pmd.load_file(get_fn('ala3_solv.crd')).coordinates
        self.sys1.box = [3.271195e1, 3.299596e1, 3.300715e1, 90, 90, 90]
        self.sys1.load_parameters(
                pmd.charmm.CharmmParameterSet(
                        get_fn('par_all36_prot.prm'),
                        get_fn('toppar_water_ions.str')
                )
        )
        self.sys2 = pmd.load_file(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        self.sys3 = pmd.load_file(get_fn(os.path.join('01.1water', 'topol.top')),
                                  xyz=get_fn(os.path.join('01.1water', 'conf.gro')))
        self.sys4 = pmd.load_file(get_fn('ala_ala_ala.psf'))
        self.sys4.coordinates = np.random.rand(len(self.sys4.atoms), 3)
        super(TestStructureSave, self).setUp()

    def test_save_pdb(self):
        """ Test saving various Structure instances as a PDB """
        self.sys1.save(self.get_fn('test.pdb', written=True))
        self.sys2.save(self.get_fn('test2.pdb', written=True))
        self.sys3.save(self.get_fn('test3.pdb', written=True))
        self.sys4.save(self.get_fn('test4.pdb', written=True))
        stringio_file = StringIO()
        self.sys4.save(stringio_file, format='pdb')
        stringio_file.seek(0)
        self.assertTrue(pmd.formats.PDBFile.id_format(stringio_file))
        stringio_file.seek(0)
        x1 = pmd.formats.PDBFile.parse(self.get_fn('test.pdb', written=True))
        x2 = pmd.formats.PDBFile.parse(self.get_fn('test2.pdb', written=True))
        x3 = pmd.formats.PDBFile.parse(self.get_fn('test3.pdb', written=True))
        x4 = pmd.formats.PDBFile.parse(self.get_fn('test4.pdb', written=True))
        x5 = pmd.formats.PDBFile.parse(stringio_file)
        self.assertEqual([a.name for a in self.sys1.atoms], [a.name for a in x1.atoms])
        self.assertEqual([a.name for a in self.sys2.atoms], [a.name for a in x2.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x3.atoms])
        self.assertEqual([a.name for a in self.sys4.atoms], [a.name for a in x4.atoms])
        self.assertEqual([a.name for a in self.sys4.atoms], [a.name for a in x5.atoms])
        # Make sure atom types as integers are preserved
        self.sys4.load_parameters(
                pmd.charmm.CharmmParameterSet(
                    get_fn('top_all22_prot.inp'),
                    get_fn('par_all22_prot.inp')
                )
        )
        for a in self.sys4.atoms:
            a.type = int(a.atom_type)
        self.sys4.save(self.get_fn('test5.pdb', written=True))
        for a in self.sys4.atoms:
            self.assertIsInstance(a.type, int)

        # raise if specifying file-like object and format is None
        stringio_file = StringIO()
        self.assertRaises(RuntimeError, lambda: self.sys4.save(stringio_file))

    def test_save_cif(self):
        """ Test saving various Structure instances as a PDBx/mmCIF """
        self.sys1.save(self.get_fn('test.cif', written=True))
        self.sys2.save(self.get_fn('test2.cif', written=True))
        self.sys3.save(self.get_fn('test3.cif', written=True))
        stringio_file = StringIO()
        self.sys3.save(stringio_file, format='cif')
        stringio_file.seek(0)
        x1 = pmd.formats.CIFFile.parse(self.get_fn('test.cif', written=True))
        x2 = pmd.formats.CIFFile.parse(self.get_fn('test2.cif', written=True))
        x3 = pmd.formats.CIFFile.parse(self.get_fn('test3.cif', written=True))
        x4 = pmd.formats.CIFFile.parse(stringio_file)
        self.assertEqual([a.name for a in self.sys1.atoms], [a.name for a in x1.atoms])
        self.assertEqual([a.name for a in self.sys2.atoms], [a.name for a in x2.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x3.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x4.atoms])
        # Try a gzip and a bzip2 file
        self.sys1.save(self.get_fn('test.cif.bz2', written=True))
        self.sys1.save(self.get_fn('test.cif.gz', written=True))
        # Make sure they're compressed
        bz2.BZ2File(self.get_fn('test.cif.bz2', written=True), 'r').read()
        gzip.open(self.get_fn('test.cif.gz', written=True), 'r').read()
        # Make sure they're the right format
        self.assertTrue(pmd.formats.CIFFile.id_format(self.get_fn('test.cif.gz', written=True)))
        self.assertTrue(pmd.formats.CIFFile.id_format(self.get_fn('test.cif.bz2', written=True)))

    def test_save_mol2(self):
        """ Test saving various Structure instances as Mol2 files """
        self.sys1.save(self.get_fn('test.mol2', written=True))
        self.sys2.save(self.get_fn('test2.mol2', written=True))
        self.sys3.save(self.get_fn('test3.mol2', written=True))
        stringio_file = StringIO()
        self.sys3.save(stringio_file, format='mol2')
        stringio_file.seek(0)
        x1 = pmd.formats.Mol2File.parse(self.get_fn('test.mol2', written=True), structure=True)
        x2 = pmd.formats.Mol2File.parse(self.get_fn('test2.mol2', written=True), structure=True)
        x3 = pmd.formats.Mol2File.parse(self.get_fn('test3.mol2', written=True), structure=True)
        x4 = pmd.formats.Mol2File.parse(stringio_file, structure=True)
        self.assertEqual([a.name for a in self.sys1.atoms], [a.name for a in x1.atoms])
        self.assertEqual([a.name for a in self.sys2.atoms], [a.name for a in x2.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x3.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x4.atoms])
        self.assertEqual(len(self.sys1.bonds), len(x1.bonds))
        self.assertEqual(len(self.sys2.bonds), len(x2.bonds))
        self.assertEqual(len(self.sys3.bonds), len(x3.bonds))

    def test_save_mol3(self):
        """ Test saving various Structure instances as Mol3 files """
        self.sys1.save(self.get_fn('test.mol3', written=True))
        stringio_file = StringIO()
        self.sys1.save(stringio_file, format='mol3')
        stringio_file.seek(0)
        x1 = pmd.formats.Mol2File.parse(self.get_fn('test.mol3', written=True), structure=True)
        x2 = pmd.formats.Mol2File.parse(stringio_file, structure=True)
        self.assertEqual([a.name for a in self.sys1.atoms], [a.name for a in x1.atoms])
        self.assertEqual([a.name for a in self.sys1.atoms], [a.name for a in x2.atoms])
        self.assertEqual(len(self.sys1.bonds), len(x1.bonds))
        with open(self.get_fn('test.mol3', written=True), 'r') as f:
            for line in f:
                if line.startswith('@<TRIPOS>HEADTAIL'):
                    break
            else:
                self.assertFalse(True)

    def test_save_amber_parm(self):
        """ Test saving various Structure instances as Amber prmtop files """
        self.sys1.save(self.get_fn('test.parm7', written=True))
        self.sys2.save(self.get_fn('test2.parm7', written=True))
        self.sys3.save(self.get_fn('test3.parm7', written=True))
        stringio_file = StringIO()
        self.sys3.save(stringio_file, format='amber')
        stringio_file.seek(0)
        x1 = pmd.amber.LoadParm(self.get_fn('test.parm7', written=True))
        x2 = pmd.amber.LoadParm(self.get_fn('test2.parm7', written=True))
        x3 = pmd.amber.LoadParm(self.get_fn('test3.parm7', written=True))
        x4 = pmd.amber.LoadParm(stringio_file)
        self.assertIsInstance(x1, pmd.amber.ChamberParm)
        self.assertIsInstance(x2, pmd.amber.AmberParm)
        self.assertIsInstance(x3, pmd.amber.AmberParm)
        self.assertIsInstance(x4, pmd.amber.AmberParm)
        # Check equivalence of topologies
        self.assertEqual([a.name for a in self.sys1.atoms], [a.name for a in x1.atoms])
        self.assertEqual([a.name for a in self.sys2.atoms], [a.name for a in x2.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x3.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x4.atoms])
        self.assertEqual(len(self.sys1.bonds), len(x1.bonds))
        self.assertEqual(len(self.sys2.bonds), len(x2.bonds))
        self.assertEqual(len(self.sys3.bonds), len(x3.bonds))
        self.assertEqual(len(self.sys1.angles), len(x1.angles))
        self.assertEqual(len(self.sys2.angles), len(x2.angles))
        self.assertEqual(len(self.sys3.angles), len(x3.angles))
        self.assertEqual(sum([len(d.type) for d in self.sys1.dihedrals]), len(x1.dihedrals))
        self.assertEqual(len(self.sys2.dihedrals), len(x2.dihedrals))
        self.assertEqual(sum([len(d.type) for d in self.sys3.dihedrals]), len(x3.dihedrals))
        self.assertEqual(len(self.sys1.impropers), len(x1.impropers))
        self.assertEqual(len(self.sys2.impropers), len(x2.impropers))
        self.assertEqual(len(self.sys3.impropers), len(x3.impropers))
        self.assertEqual(len(self.sys1.cmaps), len(x1.cmaps))
        self.assertEqual(len(self.sys2.cmaps), len(x2.cmaps))
        self.assertEqual(len(self.sys3.cmaps), len(x3.cmaps))
        for a1, a2 in zip(self.sys1.atoms, x1.atoms):
            self.assertAlmostEqual(a1.rmin, a2.rmin)
            self.assertAlmostEqual(abs(a1.epsilon), abs(a2.epsilon))
        for a1, a2 in zip(self.sys2.atoms, x2.atoms):
            self.assertAlmostEqual(a1.rmin, a2.rmin)
            self.assertAlmostEqual(abs(a1.epsilon), abs(a2.epsilon))
        for a1, a2 in zip(self.sys3.atoms, x3.atoms):
            self.assertAlmostEqual(a1.rmin, a2.rmin)
            self.assertAlmostEqual(abs(a1.epsilon), abs(a2.epsilon))
        # Now try the Amoeba topology
        parm = pmd.load_file(get_fn('nma.parm7'))
        parm.save(self.get_fn('test4.parm7', written=True))
        self.assertIsInstance(pmd.load_file(self.get_fn('test4.parm7', written=True)), pmd.amber.AmoebaParm)

    @unittest.skipUnless(HAS_GROMACS, 'Cannot test without GROMACS')
    def test_save_amber_parm2(self):
        """ Test saving AmberParm with custom exceptions """
        parm = pmd.load_file(os.path.join(get_fn('04.Ala'), 'topol.top'),
                             xyz=os.path.join(get_fn('04.Ala'), 'conf.gro'))
        fn = self.get_fn('test.parm7', written=True)
        parm.save(fn, overwrite=True)
        self.assertIs(type(pmd.load_file(fn)), pmd.amber.AmberParm)
        # Now modify one of the exceptions
        parm.adjust_types.append(copy(parm.adjust_types[0]))
        parm.adjust_types[-1].rmin = 0.1
        parm.adjust_types[-1].epsilon = 0.1
        for adjust in parm.adjusts:
            if adjust.atom1 not in adjust.atom2.dihedral_partners:
                continue
            break
        else:
            assert False, 'No pair types that are 1-4 pairs'
        adjust.type = parm.adjust_types[-1]
        parm.save(fn, overwrite=True)
        # Now it should be a ChamberParm
        self.assertIs(type(pmd.load_file(fn)), pmd.amber.ChamberParm)

    def test_save_psf(self):
        """ Test saving various Structure instances as CHARMM PSF files """
        self.sys1.save(self.get_fn('test.psf', written=True))
        self.sys2.save(self.get_fn('test2.psf', written=True))
        self.sys3.save(self.get_fn('test3.psf', written=True))
        stringio_file = StringIO()
        self.sys3.save(stringio_file, format='psf')
        stringio_file.seek(0)
        x1 = pmd.charmm.CharmmPsfFile(self.get_fn('test.psf', written=True))
        x2 = pmd.charmm.CharmmPsfFile(self.get_fn('test2.psf', written=True))
        x3 = pmd.charmm.CharmmPsfFile(self.get_fn('test3.psf', written=True))
        x4 = pmd.charmm.CharmmPsfFile(stringio_file)
        # PSF files save "improper periodic" torsions under the improper list,
        # only moving them over to the dihedral list once parameters have been
        # assigned and the fact that it's a periodic improper torsion becomes
        # known. Since no parameters are assigned in this section, we need to
        # have a way of making sure that we don't get a false positive based on
        # this difference. Add methods to determine the numbers of proper and
        # improper torsions
        def _propers(struct):
            # Only uniques
            nnormal = 0
            added = set()
            for dih in struct.dihedrals:
                a1, a2, a3, a4 = dih.atom1, dih.atom2, dih.atom3, dih.atom4
                if dih.improper: continue
                if (a1, a2, a3, a4) in added or (a4, a3, a2, a1) in added:
                    continue
                nnormal += 1
                added.add((a1, a2, a3, a4))
            return nnormal
        def _impropers(struct):
            return sum(1 for dih in struct.dihedrals if dih.improper) + len(struct.impropers)
        # Check equivalence of topologies
        self.assertEqual([a.name for a in self.sys1.atoms], [a.name for a in x1.atoms])
        self.assertEqual([a.name for a in self.sys2.atoms], [a.name for a in x2.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x3.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x4.atoms])
        self.assertEqual(len(self.sys1.bonds), len(x1.bonds))
        self.assertEqual(len(self.sys2.bonds), len(x2.bonds))
        self.assertEqual(len(self.sys3.bonds), len(x3.bonds))
        self.assertEqual(len(self.sys1.angles), len(x1.angles))
        self.assertEqual(len(self.sys2.angles), len(x2.angles))
        self.assertEqual(len(self.sys3.angles), len(x3.angles))
        self.assertEqual(_propers(self.sys1), len(x1.dihedrals))
        self.assertEqual(_propers(self.sys2), len(x2.dihedrals))
        self.assertEqual(_propers(self.sys3), len(x3.dihedrals))
        self.assertEqual(_impropers(self.sys1), len(x1.impropers))
        self.assertEqual(_impropers(self.sys2), len(x2.impropers))
        self.assertEqual(_impropers(self.sys3), len(x3.impropers))
        self.assertEqual(len(self.sys1.cmaps), len(x1.cmaps))
        self.assertEqual(len(self.sys2.cmaps), len(x2.cmaps))
        self.assertEqual(len(self.sys3.cmaps), len(x3.cmaps))

    def test_save_charmm_crd(self):
        """ Test saving various Structure instances as CHARMM coord files """
        self.sys1.save(self.get_fn('test.crd', written=True))
        self.sys2.save(self.get_fn('test2.crd', written=True))
        self.sys3.save(self.get_fn('test3.crd', written=True))
        x1 = pmd.charmm.CharmmCrdFile(self.get_fn('test.crd', written=True))
        x2 = pmd.charmm.CharmmCrdFile(self.get_fn('test2.crd', written=True))
        x3 = pmd.charmm.CharmmCrdFile(self.get_fn('test3.crd', written=True))

        np.testing.assert_allclose(self.sys1.coordinates,
                x1.coordinates.reshape(self.sys1.coordinates.shape))
        np.testing.assert_allclose(self.sys2.coordinates,
                x2.coordinates.reshape(self.sys2.coordinates.shape))
        np.testing.assert_allclose(self.sys3.coordinates,
                x3.coordinates.reshape(self.sys3.coordinates.shape))

    def test_save_gromacs(self):
        """ Test saving various Structure instances as GROMACS top files """
        self.sys1.save(self.get_fn('test.top', written=True))
        self.sys2.save(self.get_fn('test2.top', written=True))
        self.sys3.save(self.get_fn('test3.top', written=True))
        stringio_file = StringIO()
        self.sys3.save(stringio_file, format='gromacs')
        stringio_file.seek(0)
        x1 = pmd.gromacs.GromacsTopologyFile(self.get_fn('test.top', written=True))
        x2 = pmd.gromacs.GromacsTopologyFile(self.get_fn('test2.top', written=True))
        x3 = pmd.gromacs.GromacsTopologyFile(self.get_fn('test3.top', written=True))
        x4 = pmd.gromacs.GromacsTopologyFile(stringio_file)
        # Check equivalence of topologies
        self.assertEqual([a.name for a in self.sys1.atoms], [a.name for a in x1.atoms])
        self.assertEqual([a.name for a in self.sys2.atoms], [a.name for a in x2.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x3.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x4.atoms])
        self.assertEqual(len(self.sys1.bonds), len(x1.bonds))
        self.assertEqual(len(self.sys2.bonds), len(x2.bonds))
        self.assertEqual(len(self.sys3.bonds), len(x3.bonds))
        self.assertEqual(len(self.sys1.angles), len(x1.angles))
        self.assertEqual(len(self.sys2.angles), len(x2.angles))
        self.assertEqual(len(self.sys3.angles), len(x3.angles))
        self.assertEqual(len(self.sys1.dihedrals), len(x1.dihedrals))
        self.assertEqual(len(self.sys2.dihedrals), len(x2.dihedrals))
        self.assertEqual(len(self.sys3.dihedrals), len(x3.dihedrals))
        self.assertEqual(len(self.sys1.impropers), len(x1.impropers))
        self.assertEqual(len(self.sys2.impropers), len(x2.impropers))
        self.assertEqual(len(self.sys3.impropers), len(x3.impropers))
        self.assertEqual(len(self.sys1.cmaps), len(x1.cmaps))
        self.assertEqual(len(self.sys2.cmaps), len(x2.cmaps))
        self.assertEqual(len(self.sys3.cmaps), len(x3.cmaps))

    def test_save_psf2(self):
        """ Test saving PSF file for unparametrized system """
        url = 'http://ambermd.org/tutorials/advanced/tutorial1/files/polyAT.pdb'
        pmd.load_file(url).save(self.get_fn('test.psf', written=True))

    def test_save_gro(self):
        """ Test saving various Structure instances as a PDB """
        self.sys1.save(self.get_fn('test.gro', written=True))
        self.sys2.save(self.get_fn('test2.gro', written=True))
        self.sys3.save(self.get_fn('test3.gro', written=True))
        stringio_file = StringIO()
        self.sys3.save(stringio_file, format='gro')
        stringio_file.seek(0)
        x1 = pmd.gromacs.GromacsGroFile.parse(self.get_fn('test.gro', written=True))
        x2 = pmd.gromacs.GromacsGroFile.parse(self.get_fn('test2.gro', written=True))
        x3 = pmd.gromacs.GromacsGroFile.parse(self.get_fn('test3.gro', written=True))
        x4 = pmd.gromacs.GromacsGroFile.parse(stringio_file)
        self.assertEqual([a.name for a in self.sys1.atoms], [a.name for a in x1.atoms])
        self.assertEqual([a.name for a in self.sys2.atoms], [a.name for a in x2.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x3.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x4.atoms])

    def test_save_rst7(self):
        """ Test saving various Structure instances as Amber ASCII restarts """
        f1 = self.get_fn('test.rst7', written=True)
        f2 = self.get_fn('test1.restrt', written=True)
        f3 = self.get_fn('test2.inpcrd', written=True)
        f4 = self.get_fn('test3.amberrst', written=True)
        stringio_file = StringIO()
        self.sys1.save(f1)
        self.sys2.save(f2)
        self.sys3.save(f3)
        self.sys1.save(f4, format='rst7')
        self.sys1.save(stringio_file, format='rst7')
        stringio_file.seek(0)

        self.assertTrue(pmd.amber.AmberAsciiRestart.id_format(f1))
        self.assertTrue(pmd.amber.AmberAsciiRestart.id_format(f2))
        self.assertTrue(pmd.amber.AmberAsciiRestart.id_format(f3))
        self.assertTrue(pmd.amber.AmberAsciiRestart.id_format(f4))
        self.assertTrue(pmd.amber.AmberAsciiRestart.id_format(stringio_file))

        np.testing.assert_allclose(self.sys1.coordinates,
                                   pmd.load_file(f1).coordinates[0],
                                   atol=1e-6)
        np.testing.assert_allclose(self.sys2.coordinates,
                                   pmd.load_file(f2).coordinates[0],
                                   atol=1e-6)
        np.testing.assert_allclose(self.sys3.coordinates,
                                   pmd.load_file(f3).coordinates[0],
                                   atol=1e-6)
        np.testing.assert_allclose(self.sys1.coordinates,
                                   pmd.load_file(f4).coordinates[0],
                                   atol=1e-6)

    @unittest.skipIf(PYPY, 'NetCDF tests cannot run on pypy yet')
    def test_save_ncrst7(self):
        """ Test saving various Structure instances as Amber NetCDF restarts """
        f1 = self.get_fn('test.ncrst', written=True)
        f2 = self.get_fn('test1.ncrst', written=True)
        f3 = self.get_fn('test2.ncrestart', written=True)
        self.sys1.save(f1)
        self.sys2.save(f2)
        self.sys3.save(f3, format='ncrst')

        self.assertTrue(pmd.amber.NetCDFRestart.id_format(f1))
        self.assertTrue(pmd.amber.NetCDFRestart.id_format(f2))
        self.assertTrue(pmd.amber.NetCDFRestart.id_format(f3))

        np.testing.assert_allclose(self.sys1.coordinates,
                                   pmd.load_file(f1).coordinates[0])
        np.testing.assert_allclose(self.sys2.coordinates,
                                   pmd.load_file(f2).coordinates[0])
        np.testing.assert_allclose(self.sys3.coordinates,
                                   pmd.load_file(f3).coordinates[0])

    def test_save_pqr(self):
        """ Test saving various Structure instances as PQR files """
        f1 = self.get_fn('test', written=True)
        f2 = self.get_fn('test.pqr', written=True)
        f3 = self.get_fn('test2.pqr', written=True)
        self.sys1.save(f1, format='pqr')
        self.sys2.save(f2)
        self.sys3.save(f3)

        self.assertTrue(pmd.formats.PQRFile.id_format(f1))
        self.assertTrue(pmd.formats.PQRFile.id_format(f2))
        self.assertTrue(pmd.formats.PQRFile.id_format(f3))

        x1 = pmd.formats.PQRFile.parse(f1)
        x2 = pmd.formats.PQRFile.parse(f2)
        x3 = pmd.formats.PQRFile.parse(f3)

        self.assertEqual([a.name for a in self.sys1.atoms], [a.name for a in x1.atoms])
        self.assertEqual([a.name for a in self.sys2.atoms], [a.name for a in x2.atoms])
        self.assertEqual([a.name for a in self.sys3.atoms], [a.name for a in x3.atoms])

        np.testing.assert_allclose(np.array([a.charge for a in self.sys1.atoms]),
                                   np.array([a.charge for a in x1.atoms]))

        np.testing.assert_allclose(np.array([a.charge for a in self.sys2.atoms]),
                                   np.array([a.charge for a in x2.atoms]))

        np.testing.assert_allclose(np.array([a.charge for a in self.sys3.atoms]),
                                   np.array([a.charge for a in x3.atoms]))

    def test_error_handling(self):
        """ Tests handling of bad input """
        self.assertRaises(ValueError, lambda:
                self.sys1.save('somefile', format='NOFMT'))
        self.assertRaises(ValueError, lambda: self.sys1.save('somefile.nofmt'))

    def test_overwrite(self):
        """ Test overwrite option of Structure.save """
        open(self.get_fn('test.pdb', written=True), 'w').close()
        with self.assertRaises(IOError):
            self.sys1.save(self.get_fn('test.pdb', written=True))
        # Should not raise
        self.sys1.save(self.get_fn('test.pdb', written=True), overwrite=True)
        pdb = pmd.load_file(self.get_fn('test.pdb', written=True))
        self.assertEqual(len(pdb.atoms), len(self.sys1.atoms))
