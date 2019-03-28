"""
Unittests for serializing various objects in ParmEd
"""
from __future__ import division

from io import BytesIO
import numpy as np
import os
import parmed as pmd
from parmed.utils.six.moves import range, zip
try:
    import cPickle as pickle
except ImportError:
    import pickle
import random
try:
    from string import uppercase
except ImportError:
    from string import ascii_uppercase as uppercase
import unittest
import utils
from utils import HAS_GROMACS

class TestParmedSerialization(unittest.TestCase):
    """ Tests ParmEd serialization """

    _equal_atoms = utils.equal_atoms

    def test_atom_serialization(self):
        """ Tests the serialization of Atom """
        atom = pmd.Atom(atomic_number=random.randint(1, 100),
                        name=random.choice(uppercase)+random.choice(uppercase),
                        type=random.choice(uppercase)+random.choice(uppercase),
                        charge=random.random()*2-1, mass=random.random()*30+1,
                        nb_idx=random.randint(1, 20),
                        solvent_radius=random.random()*2,
                        screen=random.random()*2, tree='M',
                        join=random.random()*2, irotat=random.random(),
                        occupancy=random.random(), bfactor=random.random()*10,
                        altloc=random.choice(uppercase), rmin=random.random()*2,
                        epsilon=random.random()/2, rmin14=random.random()*2,
                        epsilon14=random.random()/2)
        atom.xx, atom.xy, atom.xz = (random.random()*100-50 for i in range(3))
        atom.number = random.randint(1, 100)
        atom.vx, atom.vy, atom.vz = (random.random()*100-50 for i in range(3))
        atom.multipoles = np.random.rand(10) * 10

        fobj = BytesIO()
        pickle.dump(atom, fobj)
        fobj.seek(0)
        unpickled = pickle.load(fobj)

        self.assertIsInstance(unpickled, pmd.Atom)
        self._equal_atoms(unpickled, atom)

    def test_bond_serialization(self):
        """ Tests the serialization of Bond """
        struct = utils.create_random_structure(True)
        bond = struct.bonds[0]
        fobj = BytesIO()
        pickle.dump(bond, fobj)
        fobj.seek(0)
        unpickled = pickle.load(fobj)

        self.assertIsInstance(bond, pmd.Bond)

    def test_bondtype_serialization(self):
        """ Tests the serialization of BondType """
        struct = utils.create_random_structure(True)
        bt = struct.bond_types[0]

        fobj = BytesIO()
        pickle.dump(bt, fobj)
        fobj.seek(0)

        unpickled = pickle.load(fobj)

        self.assertEqual(unpickled, bt)
        self.assertIsNot(unpickled, bt)

    def test_residue_serialization(self):
        """ Tests the serialization of Residue """
        struct = utils.create_random_structure(parametrized=True)
        res = struct.residues[0]

        fobj = BytesIO()
        pickle.dump(res, fobj)
        fobj.seek(0)
        unpickled = pickle.load(fobj)

        self.assertEqual(len(res.atoms), len(unpickled.atoms))
        for a1, a2 in zip(res, unpickled):
            self._equal_atoms(a1, a2)
            self.assertIs(a1.residue, res)
            self.assertIs(a2.residue, unpickled)

    def test_structure_serialization(self):
        """ Tests the serialization of Structure """
        structure = utils.create_random_structure(parametrized=True)
        # Make sure we copy over exclusions
        structure.atoms[0].exclude(structure.atoms[10])
        fobj = BytesIO()
        pickle.dump(structure, fobj)
        fobj.seek(0)
        unpickled = pickle.load(fobj)

        self._compare_structures(unpickled, structure)

    def test_fortran_format_serialization(self):
        """ Tests the serialization of FortranFormat """
        fmt = pmd.amber.FortranFormat('8I10')
        unpickled = pickle.loads(pickle.dumps(fmt))

        self.assertEqual(fmt.format, unpickled.format)
        self.assertEqual(fmt.strip_strings, unpickled.strip_strings)
        self.assertIs(fmt.type, unpickled.type)
        self.assertEqual(fmt.nitems, unpickled.nitems)
        self.assertEqual(fmt.itemlen, unpickled.itemlen)
        self.assertEqual(fmt.fmt, unpickled.fmt)

    def test_amberformat_serialization(self):
        """ Tests the serialization of AmberFormat """
        amber = pmd.load_file(utils.get_fn('cSPCE.mdl'))
        unpickled = pickle.loads(pickle.dumps(amber))

        self.assertEqual(set(amber.parm_data.keys()), set(unpickled.parm_data.keys()))
        self.assertEqual(amber.flag_list, unpickled.flag_list)
        self.assertEqual(set(amber.formats.keys()), set(unpickled.formats.keys()))
        for k1 in amber.parm_data.keys():
            self.assertEqual(amber.parm_data[k1], unpickled.parm_data[k1])
            self.assertEqual(amber.formats[k1], unpickled.formats[k1])

        self.assertEqual(amber.charge_flag, unpickled.charge_flag)
        self.assertEqual(amber.version, unpickled.version)
        self.assertEqual(amber.name, unpickled.name)

    def test_amberparm_serialization(self):
        """ Tests the serialization of AmberParm """
        structure = pmd.load_file(utils.get_fn('ash.parm7'))
        unpickled = pickle.loads(pickle.dumps(structure))
        self.assertFalse(structure.unknown_functional)
        self.assertFalse(structure.unknown_functional)

        self._compare_structures(unpickled, structure)

        self.assertEqual(set(structure.parm_data.keys()), set(unpickled.parm_data.keys()))
        self.assertEqual(structure.flag_list, unpickled.flag_list)
        self.assertEqual(set(structure.formats.keys()), set(unpickled.formats.keys()))
        for k1 in structure.parm_data.keys():
            self.assertEqual(structure.parm_data[k1], unpickled.parm_data[k1])
            self.assertEqual(structure.formats[k1], unpickled.formats[k1])

        self.assertEqual(structure.charge_flag, unpickled.charge_flag)
        self.assertEqual(structure.version, unpickled.version)
        self.assertEqual(structure.name, unpickled.name)
        self.assertIs(pmd.amber.AmberParm, type(unpickled))

        for key in 'pointers LJ_types LJ_radius LJ_depth'.split():
            self.assertEqual(hasattr(structure, key),
                             hasattr(unpickled, key))
            if hasattr(structure, key):
                self.assertEqual(getattr(structure, key), getattr(unpickled, key))

        # Now check that unknown_functional gets properly deserialized
        structure.unknown_functional = True
        self.assertTrue(pickle.loads(pickle.dumps(structure)).unknown_functional)

    def test_chamberparm_serialization(self):
        """ Tests the serialization of ChamberParm """
        structure = pmd.load_file(utils.get_fn('ala_ala_ala.parm7'),
                                  utils.get_fn('ala_ala_ala.rst7'))
        unpickled = pickle.loads(pickle.dumps(structure))

        self._compare_structures(unpickled, structure)

        self.assertEqual(set(structure.parm_data.keys()), set(unpickled.parm_data.keys()))
        self.assertEqual(structure.flag_list, unpickled.flag_list)
        self.assertEqual(set(structure.formats.keys()), set(unpickled.formats.keys()))
        for k1 in structure.parm_data.keys():
            self.assertEqual(structure.parm_data[k1], unpickled.parm_data[k1])
            self.assertEqual(structure.formats[k1], unpickled.formats[k1])

        self.assertEqual(structure.charge_flag, unpickled.charge_flag)
        self.assertEqual(structure.version, unpickled.version)
        self.assertEqual(structure.name, unpickled.name)
        self.assertIs(pmd.amber.ChamberParm, type(unpickled))

    def test_amoebaparm_serialization(self):
        """ Tests the serialization of AmoebaParm """
        structure = pmd.load_file(utils.get_fn('nma.parm7'),
                                  utils.get_fn('nma.rst7'))
        unpickled = pickle.loads(pickle.dumps(structure))

        self._compare_structures(unpickled, structure)

        self.assertEqual(set(structure.parm_data.keys()), set(unpickled.parm_data.keys()))
        self.assertEqual(structure.flag_list, unpickled.flag_list)
        self.assertEqual(set(structure.formats.keys()), set(unpickled.formats.keys()))
        for k1 in structure.parm_data.keys():
            self.assertEqual(structure.parm_data[k1], unpickled.parm_data[k1])
            self.assertEqual(structure.formats[k1], unpickled.formats[k1])

        self.assertEqual(structure.charge_flag, unpickled.charge_flag)
        self.assertEqual(structure.version, unpickled.version)
        self.assertEqual(structure.name, unpickled.name)
        self.assertIs(pmd.amber.AmoebaParm, type(unpickled))

    def test_pdb_serialization(self):
        """ Tests the serialization of a parsed PDB file """
        structure = pmd.load_file(utils.get_fn('4lzt.pdb'))
        unpickled = pickle.loads(pickle.dumps(structure))

        self._compare_structures(unpickled, structure)

        # Check metadata
        for key in ('experimental', 'journal', 'authors', 'keywords', 'doi',
                    'pmid', 'journal_authors', 'volume', 'title', 'year',
                    'resolution', 'related_entries', 'space_group'):
            self.assertEqual(getattr(structure, key), getattr(unpickled, key))

    def test_pdbtraj_serialization(self):
        """ Tests the serialization of a Structure with trajectory """
        structure = pmd.load_file(utils.get_fn('2koc.pdb'))
        unpickled = pickle.loads(pickle.dumps(structure))

        self._compare_structures(unpickled, structure)

        # Check metadata
        for key in ('experimental', 'journal', 'authors', 'keywords', 'doi',
                    'pmid', 'journal_authors', 'volume', 'title', 'year',
                    'resolution', 'related_entries'):
            self.assertEqual(getattr(structure, key), getattr(unpickled, key))

        self.assertGreater(structure.get_coordinates().shape[0], 1)

    def test_parm_velocities_serialization(self):
        """ Tests the serialization of a Structure with velocities """
        structure = pmd.load_file(utils.get_fn('tip4p.parm7'), utils.get_fn('tip4p.rst7'))
        unpickled = pickle.loads(pickle.dumps(structure))
        self._compare_structures(unpickled, structure)

    @unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without GROMACS")
    def test_gromacstop_serialization(self):
        """ Tests the serialization of a GromacsTopologyFile """
        structure = pmd.load_file(os.path.join(utils.get_fn('03.AlaGlu'), 'topol.top'),
                                  xyz=os.path.join(utils.get_fn('03.AlaGlu'), 'conf.gro'))
        unpickled = pickle.loads(pickle.dumps(structure))
        self._compare_structures(unpickled, structure)
        self.assertEqual(structure.defaults, unpickled.defaults)
        self._compare_parametersets(structure.parameterset, unpickled.parameterset)

    @unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without GROMACS")
    def test_gromacscharmm_serialization(self):
        """ Tests the serialization of a CHARMM FF Gromacs topology """
        structure = pmd.load_file(utils.get_fn('1aki.charmm27.solv.top'))
        unpickled = pickle.loads(pickle.dumps(structure))
        self._compare_structures(unpickled, structure)
        self.assertEqual(structure.defaults, unpickled.defaults)
        self._compare_parametersets(structure.parameterset, unpickled.parameterset)

    def test_charmmpsf_serialization(self):
        """ Tests the serialization of a CHARMM PSF file """
        structure = pmd.load_file(utils.get_fn('ala_ala_ala.psf'))
        unpickled = pickle.loads(pickle.dumps(structure))
        self._compare_structures(unpickled, structure)

    def _compare_structures(self, unpickled, structure):

        self.assertEqual(len(unpickled.residues), len(structure.residues))
        for r1, r2 in zip(unpickled.residues, structure.residues):
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.idx, r2.idx)
            for a1, a2 in zip(r1, r2):
                self._equal_atoms(a1, a2)

        def cmp_alists(alist1, alist2):
            self.assertEqual(len(alist1), len(alist2))
            for a1, a2 in zip(alist1, alist2):
                self._equal_atoms(a1, a2)

        for a1, a2 in zip(unpickled, structure):
            self._equal_atoms(a1, a2)
            self.assertEqual(a1.idx, a2.idx)
            self.assertEqual(len(a1.bonds), len(a2.bonds))
            self.assertEqual(len(a1.angles), len(a2.angles))
            self.assertEqual(len(a1.dihedrals), len(a2.dihedrals))
            self.assertEqual(len(a1.impropers), len(a2.impropers))
            cmp_alists(a1.bond_partners, a2.bond_partners)
            cmp_alists(a1.angle_partners, a2.angle_partners)
            cmp_alists(a1.dihedral_partners, a2.dihedral_partners)
            cmp_alists(a1.tortor_partners, a2.tortor_partners)
            cmp_alists(a1.exclusion_partners, a2.exclusion_partners)

        # Check coordinates
        if structure.get_coordinates() is None:
            self.assertIs(unpickled.get_coordinates(), None)
        else:
            np.testing.assert_equal(structure.get_coordinates(), unpickled.get_coordinates())
        # Check unit cell
        if structure.box is None:
            self.assertIs(unpickled.box, None)
            self.assertIs(unpickled.box_vectors, None)
        else:
            np.testing.assert_equal(structure.box, unpickled.box)
            self.assertEqual(structure.box_vectors, unpickled.box_vectors)

        # Check velocities
        if structure.velocities is None:
            self.assertIs(unpickled.velocities, None)
        else:
            np.testing.assert_equal(structure.velocities, unpickled.velocities)

        # Check nrexcl and combining_rule
        self.assertEqual(structure.nrexcl, unpickled.nrexcl)
        self.assertEqual(structure.combining_rule, unpickled.combining_rule)

        # Make sure all of the type arrays are equivalent
        def cmp_type_arrays(arr1, arr2):
            self.assertEqual(len(arr1), len(arr2))
            for x1, x2 in zip(arr1, arr2):
                self.assertEqual(x1, x2)

        cmp_type_arrays(structure.bond_types, unpickled.bond_types)
        cmp_type_arrays(structure.angle_types, unpickled.angle_types)
        cmp_type_arrays(structure.dihedral_types, unpickled.dihedral_types)
        cmp_type_arrays(structure.improper_types, unpickled.improper_types)
        cmp_type_arrays(structure.urey_bradley_types, unpickled.urey_bradley_types)
        cmp_type_arrays(structure.rb_torsion_types, unpickled.rb_torsion_types)
        cmp_type_arrays(structure.cmap_types, unpickled.cmap_types)
        cmp_type_arrays(structure.trigonal_angle_types, unpickled.trigonal_angle_types)
        cmp_type_arrays(structure.out_of_plane_bend_types, unpickled.out_of_plane_bend_types)
        cmp_type_arrays(structure.stretch_bend_types, unpickled.stretch_bend_types)
        cmp_type_arrays(structure.torsion_torsion_types, unpickled.torsion_torsion_types)
        cmp_type_arrays(structure.pi_torsion_types, unpickled.pi_torsion_types)
        cmp_type_arrays(structure.adjust_types, unpickled.adjust_types)

        # Make sure all of the connectivity arrays are equivalent
        def cmp_top_arrays(arr1, arr2):
            self.assertEqual(len(arr1), len(arr2))
            for t1, t2 in zip(arr1, arr2):
                self.assertIs(type(t1), type(t2))
                atoms = [attr for attr in dir(t1) if attr.startswith('atom')]
                for a in atoms:
                    self._equal_atoms(getattr(t1, a), getattr(t2, a))
                if hasattr(t1, 'type'):
                    self.assertEqual(t1.type, t2.type)

        cmp_top_arrays(structure.bonds, unpickled.bonds)
        cmp_top_arrays(structure.angles, unpickled.angles)
        cmp_top_arrays(structure.dihedrals, unpickled.dihedrals)
        cmp_top_arrays(structure.impropers, unpickled.impropers)
        cmp_top_arrays(structure.urey_bradleys, unpickled.urey_bradleys)
        cmp_top_arrays(structure.rb_torsions, unpickled.rb_torsions)
        cmp_top_arrays(structure.cmaps, unpickled.cmaps)
        cmp_top_arrays(structure.trigonal_angles, unpickled.trigonal_angles)
        cmp_top_arrays(structure.out_of_plane_bends, unpickled.out_of_plane_bends)
        cmp_top_arrays(structure.pi_torsions, unpickled.pi_torsions)
        cmp_top_arrays(structure.stretch_bends, unpickled.stretch_bends)
        cmp_top_arrays(structure.torsion_torsions, unpickled.torsion_torsions)
        cmp_top_arrays(structure.chiral_frames, unpickled.chiral_frames)
        cmp_top_arrays(structure.multipole_frames, unpickled.multipole_frames)
        cmp_top_arrays(structure.adjusts, unpickled.adjusts)
        cmp_top_arrays(structure.groups, unpickled.groups)

    def _compare_parametersets(self, set1, set2):

        def cmp_lists(l1, l2, reversible=False):
            self.assertEqual(len(l1), len(l2))
            self.assertEqual(set(l1.keys()), set(l2.keys()))
            for k in l1:
                self.assertEqual(l1[k], l2[k])
                if reversible:
                    if pmd.NoUreyBradley is l1[k] or pmd.NoUreyBradley is l2[k]:
                        self.assertIs(l1[k], pmd.NoUreyBradley)
                        self.assertIs(l2[k], pmd.NoUreyBradley)
                    self.assertIs(l1[k], l1[tuple(reversed(k))])
                    self.assertIs(l2[k], l2[tuple(reversed(k))])

        # Now compare all of the parameter lists
        cmp_lists(set1.bond_types, set2.bond_types, True)
        cmp_lists(set1.angle_types, set2.angle_types, True)
        cmp_lists(set1.dihedral_types, set2.dihedral_types, True)
        cmp_lists(set1.improper_types, set2.improper_types)
        cmp_lists(set1.improper_periodic_types, set2.improper_periodic_types)
        cmp_lists(set1.cmap_types, set2.cmap_types, True)
        cmp_lists(set1.atom_types, set2.atom_types)
        cmp_lists(set1.urey_bradley_types, set2.urey_bradley_types, True)
        cmp_lists(set1.nbfix_types, set2.nbfix_types)
        self.assertEqual(set1.combining_rule, set2.combining_rule)
