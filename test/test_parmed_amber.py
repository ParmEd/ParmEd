"""
Tests the functionality in the parmed.amber package
"""
from __future__ import print_function, division

from copy import copy
import glob
import math
import numpy as np
import os
import re
import sys
import parmed as pmd
from parmed.amber import (
    readparm, asciicrd, mask, parameters, mdin, FortranFormat, titratable_residues, AmberOFFLibrary
)
from parmed.exceptions import (
    AmberWarning, MoleculeError, AmberError, MaskError, InputError, ParameterWarning
)
from parmed.modeller import ResidueTemplateContainer
from parmed import topologyobjects, load_file, Structure
from parmed.tools import change
import parmed.unit as u
from parmed.utils import PYPY
from parmed.utils.six import string_types, iteritems
from parmed.utils.six.moves import range, zip, StringIO
import pickle
import random
import saved_outputs as saved
import shutil
import unittest
from utils import (
    get_fn, FileIOTestCase, equal_atoms, create_random_structure, HAS_GROMACS,
    diff_files, has_openmm
)
import warnings
try:
    from string import letters
except ImportError:
    from string import ascii_letters as letters

def _picklecycle(obj):
    return pickle.loads(pickle.dumps(obj))


class TestReadParm(FileIOTestCase):
    """ Tests the various Parm file classes """

    def test_fortran_format(self):
        """ Tests the FortranFormat object """
        fmt = FortranFormat('(F8.5)')
        self.assertEqual(fmt.nitems, 1)
        self.assertEqual(fmt.itemlen, 8)
        self.assertEqual(fmt.num_decimals, 5)
        self.assertEqual(fmt.fmt, '%8.5F')
        fmt = FortranFormat('(E8.5)')
        self.assertEqual(fmt.nitems, 1)
        self.assertEqual(fmt.itemlen, 8)
        self.assertEqual(fmt.num_decimals, 5)
        self.assertEqual(fmt.fmt, '%8.5E')
        self.assertEqual(repr(fmt), '<FortranFormat: (E8.5)>')
        file = StringIO()
        fmt.write(10, file)
        file.seek(0)
        self.assertEqual(file.read(), '1.00000E+01\n')
        file = StringIO()
        fmt = FortranFormat('a80')
        fmt.write('this', file)
        file.seek(0)
        self.assertEqual(file.read(), 'this' + ' '*76 + '\n')
        fmt = FortranFormat('6a12', strip_strings=False)
        fmt.read = fmt._read_nostrip
        stuff = fmt.read(' '*12 + 'abcde ' + ' '*18)
        self.assertEqual(stuff, [' '*12, 'abcde' + ' '*7, ' '*12])
        # Test hashability of FortranFormat
        obj1 = object()
        obj2 = object()
        d = {FortranFormat('6a12') : obj1, FortranFormat('6a12', False) : obj2}
        self.assertIs(d[FortranFormat('6a12')], obj1)
        self.assertIs(d[FortranFormat('6a12', False)], obj2)

    def test_amber_format(self):
        """ Test some general functionality of the AmberFormat class """
        parm = readparm.AmberFormat(self.get_fn('ash.parm7'))
        after = parm.flag_list[-1]
        parm.add_flag('NEW_FLAG', '10i6', num_items=20, after=after, comments='This is a comment')
        self.assertEqual(parm.flag_list[-1], 'NEW_FLAG')
        self.assertEqual(parm.parm_comments['NEW_FLAG'], ['This is a comment'])
        with self.assertRaises(AmberError):
            parm.add_flag('NEW_FLAG2', '10i6')

    def test_optimized_reader(self):
        """ Check that the optimized reader imports correctly """
        from parmed.amber import _rdparm

    def test_nbfix_from_structure(self):
        """ Tests AmberParm.from_structure with NBFIXes """
        s = Structure()
        at1 = topologyobjects.AtomType('A1', 1, 12.01, atomic_number=6)
        at2 = topologyobjects.AtomType('A2', 2, 12.01, atomic_number=6)
        at3 = topologyobjects.AtomType('A3', 3, 12.01, atomic_number=6)
        at4 = topologyobjects.AtomType('A4', 4, 12.01, atomic_number=6)
        at1.set_lj_params(0.5, 1.0)
        at2.set_lj_params(0.6, 1.1)
        at3.set_lj_params(0.7, 1.2)
        at4.set_lj_params(0.8, 1.3)

        # Add some NBFIXes
        at1.add_nbfix('A2', 0.8, 4.0)
        at2.add_nbfix('A1', 0.8, 4.0)
        at3.add_nbfix('A4', 0.9, 4.1)
        at4.add_nbfix('A3', 0.9, 4.1)

        # Add a handful of atoms
        s.add_atom(topologyobjects.Atom(name='AA', type='A1', atomic_number=6), 'LJR', 1)
        s.add_atom(topologyobjects.Atom(name='AB', type='A2', atomic_number=6), 'LJR', 2)
        s.add_atom(topologyobjects.Atom(name='AC', type='A3', atomic_number=6), 'LJR', 3)
        s.add_atom(topologyobjects.Atom(name='AD', type='A4', atomic_number=6), 'LJR', 4)

        # Assign the types
        s[0].atom_type = at1
        s[1].atom_type = at2
        s[2].atom_type = at3
        s[3].atom_type = at4

        self.assertTrue(s.has_NBFIX())

        # Convert to Amber topology file
        parm = readparm.AmberParm.from_structure(s)

        self.assertTrue(parm.has_NBFIX())
        np.testing.assert_allclose(parm.parm_data['LENNARD_JONES_ACOEF'],
                np.array([2048.0, 0.27487790694400016, 7713.001578629537,
                    7605.122117724271, 14202.299844691719, 25564.24320523959,
                    13860.025454471592, 25302.038907727154, 1.1579610995721004,
                    76343.16532934578])
        )
        np.testing.assert_allclose(parm.parm_data['LENNARD_JONES_BCOEF'],
                np.array([64.0, 2.097152000000001, 136.05588480000006,
                    134.1529115728351, 191.8764421334576, 267.54416639999994,
                    187.2522338751463, 264.8000511276929, 4.3578162,
                    494.2652416000001])
        )

    def test_load_parm(self):
        """ Test the arbitrary parm loader """
        parm = readparm.LoadParm(self.get_fn('trx.prmtop'))
        parm2 = readparm.AmberParm(self.get_fn('trx.prmtop'))
        self.assertIs(parm.view_as(readparm.AmberParm), parm)
        for key in parm.parm_data:
            self.assertEqual(parm.parm_data[key], parm2.parm_data[key])
        parm.box = [2*u.nanometer, 2*u.nanometer, 2*u.nanometer,
                    math.pi/2*u.radian, math.pi/2*u.radian, math.pi/2*u.radian]
        np.testing.assert_allclose(parm.box, [20, 20, 20, 90, 90, 90])
        # Now check that the box info is set properly
        crd3 = load_file(self.get_fn('solv2.rst7'))
        parm3 = readparm.LoadParm(self.get_fn('solv2.parm7'), xyz=crd3.coordinates)
        np.testing.assert_equal(parm3.box[:3], parm3.parm_data['BOX_DIMENSIONS'][1:])
        self.assertEqual(parm3.box[3], parm3.parm_data['BOX_DIMENSIONS'][0])

    def test_gzipped_parm(self):
        """ Check that gzipped prmtop files can be parsed correctly """
        parm = readparm.LoadParm(self.get_fn('small.parm7.gz'))
        self.assertEqual(parm.ptr('natom'), 864)

    def test_bzipped_parm(self):
        """ Check that bzip2ed prmtop files can be parsed correctly """
        parm = readparm.LoadParm(self.get_fn('small.parm7.bz2'))
        self.assertEqual(parm.ptr('natom'), 864)

    def test_atomic_number_setting(self):
        """ Make sure that if ATOMIC_NUMBER is -1, the mass sets the element """
        fn = self.get_fn('test.parm7', written=True)
        parm = readparm.LoadParm(get_fn('ash.parm7'))
        # Turn off atomic numbers
        change(parm, 'ATOMIC_NUMBER', ':*', -1).execute()
        parm.write_parm(fn)
        parm2 = readparm.LoadParm(fn)
        parm = readparm.LoadParm(get_fn('ash.parm7'))
        self.assertEqual(len(parm.atoms), len(parm2.atoms))
        for a1, a2 in zip(parm.atoms, parm2.atoms):
            self.assertEqual(a1.atomic_number, a2.atomic_number)

    @unittest.skipUnless(has_openmm, 'Cannot test without OpenMM')
    def test_change_detection(self):
        """ Test the is_changed function on AmberParm """
        parm = readparm.AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertFalse(parm.is_changed())
        # Find the OpenMM Topology
        top = parm.topology
        # Delete the last bond
        del parm.bonds[-1]
        # Make sure our parm is changed
        self.assertTrue(parm.is_changed())
        # Make sure our OMM topology changes correspondingly
        self.assertIsNot(top, parm.topology)

    def test_deprecations(self):
        """ Test proper deprecation of old/renamed AmberParm features """
        rst7_name = get_fn('ash.rst7')
        parm7_name = get_fn('ash.parm7')
        with self.assertWarns(DeprecationWarning):
            readparm.AmberParm(parm7_name, rst7_name=rst7_name)
        with self.assertWarns(DeprecationWarning):
            readparm.AmberParm(parm7_name, rst7_name, rst7_name=rst7_name)
        parm = readparm.AmberParm(get_fn('ash.parm7'), rst7_name=get_fn('ash.rst7'))
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))

    def test_molecule_error_detection(self):
        """ Tests noncontiguous molecule detection """
        parm = readparm.AmberParm(get_fn('things.parm7'))
        self.assertRaises(MoleculeError, lambda: parm.rediscover_molecules(fix_broken=False))

    def test_recalculate_lj(self):
        """ Test the AmberParm.recalculate_LJ() method """
        parm = readparm.AmberParm(get_fn('things.parm7'))
        self._recalculate_lj_test(parm)

    def test_recalculate_lj_after_pickling(self):
        """ Test the AmberParm.recalculate_LJ() method with a pickled object """
        parm = readparm.AmberParm(get_fn('things.parm7'))
        parm_p = _picklecycle(parm)
        self._recalculate_lj_test(parm_p)

    def _recalculate_lj_test(self, parm):
        """ run the tests for AmberParm.recalculate_LJ() """
        orig_LJ_A = np.array(parm.parm_data['LENNARD_JONES_ACOEF'])
        orig_LJ_B = np.array(parm.parm_data['LENNARD_JONES_BCOEF'])
        parm.recalculate_LJ()
        np.testing.assert_allclose(orig_LJ_A, np.array(parm.parm_data['LENNARD_JONES_ACOEF']))
        np.testing.assert_allclose(orig_LJ_B, np.array(parm.parm_data['LENNARD_JONES_BCOEF']))

    def test_detect_nbfix(self):
        """ Tests NBFIX detection for AmberParm """
        parm = readparm.AmberParm(get_fn('ash.parm7'))
        self._detect_nbfix_test(parm)

    def test_detect_nbfix_after_pickling(self):
        """ Tests NBFIX detection for AmberParm with a pickled object """
        parm = readparm.AmberParm(get_fn('ash.parm7'))
        parm_p = _picklecycle(parm)
        self._detect_nbfix_test(parm_p)

    def _detect_nbfix_test(self, parm):
        """ run the tests for NBFIX detection for AmberParm """
        self.assertFalse(parm.has_NBFIX())
        parm.parm_data['LENNARD_JONES_BCOEF'][0] = 0.0
        self.assertTrue(parm.has_NBFIX())

    def test_dihedral_reorder(self):
        """ Tests dihedral reordering if first atom in 3rd or 4th spot """
        parm = readparm.AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        parm.strip('@1-8')
        parm.remake_parm()
        # This will result in torsion terms in which the first atom in the
        # system is the third atom in the torsion. So make sure that AmberParm
        # recognizes this and reorders the indices appropriately to avoid a 0 in
        # the 3rd or 4th locations
        it = iter(parm.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'])
        for i, j, k, l, m in zip(it, it, it, it, it):
            self.assertNotEqual(k, 0)
            self.assertNotEqual(l, 0)

    def test_amber_gas_parm(self):
        """ Test the AmberParm class with a non-periodic (gas-phase) prmtop """
        parm = readparm.AmberParm(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        gasparm = readparm.AmberParm(get_fn('trx.prmtop'))
        gasparm.load_rst7(get_fn('trx.inpcrd'))
        self.assertFalse(gasparm.chamber)
        self.assertFalse(gasparm.has_cmap)
        self.assertEqual(gasparm.combining_rule, 'lorentz')

        self.assertEqual([a.xx for a in gasparm.atoms], [a.xx for a in parm.atoms])
        self.assertEqual([a.xy for a in gasparm.atoms], [a.xy for a in parm.atoms])
        self.assertEqual([a.xz for a in gasparm.atoms], [a.xz for a in parm.atoms])

        # Now run the tests for the prmtop
        self._standard_parm_tests(parm)
        self.assertFalse(parm.chamber)
        self.assertFalse(parm.amoeba)
        self.assertRaises(KeyError, lambda: parm.parm_data['BOX_DIMENSIONS'])
        self.assertEqual(parm.ptr('ifbox'), 0)

        # Now check that IFBOX is set to 3 if we set the box to something non-orthogonal and
        # non-octahedral
        parm.box = [10, 10, 10, 90, 60, 90]
        self.assertEqual(parm.ptr('ifbox'), 3)

        # Now check the restart file
        rst = readparm.Rst7.open(get_fn('trx.inpcrd'))
        coords = rst.coordinates
        vels = rst.vels
        for i, atom in enumerate(gasparm.atoms):
            np.testing.assert_allclose(coords[0,i], [atom.xx, atom.xy, atom.xz])
            np.testing.assert_allclose(vels[0,i], [atom.vx, atom.vy, atom.vz])
        # Check join_dihedrals when one DT does not exist
        parm.dihedrals[-1].type = None
        parm.join_dihedrals()

    def test_amber_mdin(self):
        """ Tests the Amber Mdin class """
        inp = mdin.Mdin()
        inp.change('cntrl', 'ntb', 2)

    def test_remake_parm(self):
        """ Tests the rebuilding of the AmberParm raw data structures """
        parm = readparm.AmberParm(get_fn('trx.prmtop'))
        parm2 = readparm.AmberParm(get_fn('trx.prmtop'))
        parm.remake_parm()
        self.assertEqual(parm.flag_list, parm2.flag_list)
        for flag in parm.flag_list:
            for x1, x2 in zip(parm.parm_data[flag], parm2.parm_data[flag]):
                if isinstance(x1, string_types) or isinstance(x2, string_types):
                    self.assertEqual(x1, x2)
                else:
                    self.assertAlmostEqual(x1, x2)
        self.assertEqual(parm.combining_rule, 'lorentz')
        self.assertEqual(parm2.combining_rule, 'lorentz')

    def test_remake_chamber_parm(self):
        """ Tests the rebuilding of the ChamberParm raw data structures """
        parm = readparm.ChamberParm(get_fn('ala_ala_ala.parm7'))
        parm2 = readparm.ChamberParm(get_fn('ala_ala_ala.parm7'))
        parm.remake_parm()
        self.assertEqual(set(parm.flag_list), set(parm2.flag_list))
        for flag in parm.flag_list:
            for x1, x2 in zip(parm.parm_data[flag], parm2.parm_data[flag]):
                if isinstance(x1, string_types) or isinstance(x2, string_types):
                    self.assertEqual(x1, x2)
                else:
                    self.assertAlmostEqual(x1, x2)
        self.assertEqual(parm.combining_rule, 'lorentz')
        self.assertEqual(parm2.combining_rule, 'lorentz')

    def test_amber_solv_parm(self):
        """ Test the AmberParm class with a periodic prmtop """
        parm = readparm.AmberParm(get_fn('solv.prmtop'), get_fn('solv.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        self._standard_parm_tests(parm)
        self._solv_pointer_tests(parm)
        self.assertFalse(parm.chamber)
        self.assertFalse(parm.amoeba)
        self.assertEqual(parm.ptr('ifbox'), 1)
        orig_box = parm.box
        orig_boxdims = parm.parm_data['BOX_DIMENSIONS'][:]
        # Strip the solvent and make sure that the box info sticks around
        parm.strip(':WAT,Cl-')
        self.assertEqual(parm.ptr('ifbox'), 1)
        np.testing.assert_equal(orig_box, parm.box)
        self.assertEqual(parm.parm_data['BOX_DIMENSIONS'], orig_boxdims)
        self.assertEqual(parm.parm_data['SOLVENT_POINTERS'], [163, 2, 3])
        self.assertEqual(parm.parm_data['ATOMS_PER_MOLECULE'], [2603, 12])
        # Make sure that setting box to None gets rid of all boxy stuff
        parm.box = None
        self.assertEqual(parm.ptr('ifbox'), 0)
        self.assertIs(parm.box, None)
        self.assertNotIn('BOX_DIMENSIONS', parm.parm_data)
        self.assertNotIn('SOLVENT_POINTERS', parm.parm_data)
        self.assertNotIn('ATOMS_PER_MOLECULE', parm.parm_data)

    def test_amber_cmap(self):
        """ Test the AmberParm class with a CMAP-containing topology """
        def check_parm_for_cmap(parm):
            self.assertEqual(len(parm.cmaps), 18)
            self.assertEqual(len(parm.cmap_types), 9)
            self.assertTrue(parm.has_cmap)
        parm = readparm.AmberParm(get_fn('amber-parm-with-cmap.parm7'))
        check_parm_for_cmap(parm)
        # Make sure serializing/deserializing works
        parm2 = pickle.loads(pickle.dumps(parm))
        check_parm_for_cmap(parm2)
        f = StringIO()
        parm.write_parm(f)
        f.seek(0)
        check_parm_for_cmap(readparm.AmberParm(f))
        f = StringIO()
        parm2.write_parm(f)
        f.seek(0)
        check_parm_for_cmap(readparm.AmberParm(f))

    def test_chamber_gas_parm(self):
        """Test the ChamberParm class with a non-periodic (gas phase) prmtop"""
        parm = readparm.ChamberParm(get_fn('ala_ala_ala.parm7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertTrue(parm.chamber)
        self.assertTrue(parm.has_cmap)
        self.assertFalse(parm.amoeba)
        self.assertEqual(parm.ptr('ifbox'), 0)
        # Make sure that a corrupted data section is properly caught
        fmt = readparm.AmberFormat(get_fn('ala_ala_ala.parm7'))
        del fmt.parm_data['CHARMM_UREY_BRADLEY'][-1]
        self.assertRaises(AmberError, lambda: readparm.ChamberParm.from_rawdata(fmt))

    def test_chamber_solv_parm(self):
        """ Test the ChamberParm class with a periodic prmtop """
        parm = readparm.ChamberParm(get_fn('dhfr_cmap_pbc.parm7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        self._standard_parm_tests(parm)
        self._solv_pointer_tests(parm)
        self.assertTrue(parm.chamber)
        self.assertTrue(parm.has_cmap)
        self.assertEqual(parm.ptr('ifbox'), 1)

    def test_chamber_eliminate_cmap(self):
        """ Tests that CMAP flags are properly disposed of when they are deleted """
        parm = readparm.ChamberParm(get_fn('ala_ala_ala.parm7'))
        for cmap in parm.cmaps:
            cmap.delete()
        del parm.cmaps[:]
        del parm.cmap_types[:]
        parm.remake_parm()
        for flag in parm.parm_data:
            self.assertFalse(flag.startswith('CHARMM_CMAP'))
        self.assertFalse(parm.has_cmap)

    def test_amoeba_big(self):
        """ Test the AmoebaParm class with a large system """
        parm = readparm.AmoebaParm(get_fn('amoeba.parm7'))
        self.assertEqual(parm.ptr('natom'), len(parm.atoms))
        self.assertEqual([a.name for a in parm.atoms], parm.parm_data['ATOM_NAME'])
        self.assertEqual([a.type for a in parm.atoms], parm.parm_data['AMBER_ATOM_TYPE'])
        self.assertTrue(parm.amoeba)
        for attr in ['bonds', 'angles', 'urey_bradleys', 'angles',
                     'trigonal_angles', 'out_of_plane_bends', 'dihedrals',
                     'pi_torsions', 'stretch_bends', 'torsion_torsions',
                     'chiral_frames', 'multipole_frames', 'adjusts',
                     'adjust_types', 'bond_types', 'angle_types',
                     'trigonal_angle_types', 'out_of_plane_bend_types',
                     'dihedral_types', 'pi_torsion_types', 'stretch_bend_types',
                     'torsion_torsion_types']:
            self.assertTrue(hasattr(parm, attr))
        # Check that TORSION_TORSION data is restored properly
        used_tortors = sorted(set(tt.type.idx+1 for tt in parm.torsion_torsions))
        tortordata = dict()
        for flag, data in iteritems(parm.parm_data):
            if 'TORSION_TORSION' not in flag: continue
            tortordata[flag] = np.array(data)

        parm.remake_parm()

        myre = re.compile('AMOEBA_TORSION_TORSION_TORTOR_TABLE_(\d\d)')
        for flag, data in iteritems(parm.parm_data):
            if 'TORSION_TORSION' not in flag: continue
            rematch = myre.match(flag)
            if rematch:
                idx = int(rematch.groups()[0]) - 1
                tortornum = used_tortors[idx]
                oldflag = flag.replace(rematch.groups()[0], '%02d' % tortornum)
                np.testing.assert_equal(np.array(data), tortordata[oldflag])
            elif flag == 'AMOEBA_TORSION_TORSION_LIST':
                # We need to adjust for the unused types that we deleted
                for i, (new, old) in enumerate(zip(data, tortordata[flag])):
                    if i % 6 == 5:
                        self.assertEqual(used_tortors[new-1], old)
                    else:
                        self.assertEqual(new, old)
            elif flag == 'AMOEBA_TORSION_TORSION_NUM_PARAMS':
                self.assertEqual(tortordata[flag], [8])
                self.assertEqual(data, [3])
            else:
                np.testing.assert_equal(tortordata[flag], np.array(data))

        # Now check getting rid of all torsion-torsions
        del parm.torsion_torsions[:], parm.torsion_torsion_types[:]
        parm.remake_parm()
        for flag in parm.parm_data:
            self.assertFalse(flag.startswith('AMOEBA_TORSION_TORSION'))

    def test_amoeba_small(self):
        """ Test the AmoebaParm class w/ small system (not all terms) """
        parm = readparm.AmoebaParm(get_fn('nma.parm7'), get_fn('nma.rst'))
        rst7 = readparm.BeemanRestart(get_fn('nma.rst'))
        self.assertEqual(3*rst7.natom, len(rst7.coordinates.flatten()))
        self.assertEqual(rst7.natom, parm.ptr('natom'))
        self.assertFalse(parm.torsion_torsions)
        self.assertTrue(parm.amoeba)
        coords = np.random.rand(len(parm.atoms), 3)
        with self.assertRaises(TypeError):
            parm.initialize_topology(xyz=get_fn('ala3_solv.parm7'))
        parm.initialize_topology(xyz=coords, box=[1, 1, 1, 90, 90, 90])
        np.testing.assert_allclose(parm.coordinates, coords)
        np.testing.assert_equal(parm.box, [1, 1, 1, 90, 90, 90])
        self.assertEqual(saved.AMOEBA_SMALL_MDIN, parm.mdin_skeleton())
        self.assertEqual(len(parm.urey_bradleys), 818)
        self.assertEqual(len(parm.urey_bradley_types), 1)
        # Now make sure that the future layout of the stretch bend force
        # constants will work (i.e., when it is split into 2 fields)
        parm.add_flag(
            'AMOEBA_STRETCH_BEND_FORCE_CONSTANT_1',
            str(parm.formats['AMOEBA_STRETCH_BEND_FORCE_CONSTANT']),
            data=parm.parm_data['AMOEBA_STRETCH_BEND_FORCE_CONSTANT'][:],
            after='AMOEBA_STRETCH_BEND_FORCE_CONSTANT',
        )
        parm.add_flag(
            'AMOEBA_STRETCH_BEND_FORCE_CONSTANT_2',
            str(parm.formats['AMOEBA_STRETCH_BEND_FORCE_CONSTANT']),
            data=parm.parm_data['AMOEBA_STRETCH_BEND_FORCE_CONSTANT'][:],
            after='AMOEBA_STRETCH_BEND_FORCE_CONSTANT_1',
        )
        parm.delete_flag('AMOEBA_STRETCH_BEND_FORCE_CONSTANT')
        parm2 = readparm.AmoebaParm.from_rawdata(parm)
        self.assertEqual(len(parm.stretch_bends), len(parm2.stretch_bends))
        for sb1, sb2 in zip(parm.stretch_bends, parm2.stretch_bends):
            self.assertEqual(sb1.atom1.idx, sb2.atom1.idx)
            self.assertEqual(sb1.atom2.idx, sb2.atom2.idx)
            self.assertEqual(sb1.atom3.idx, sb2.atom3.idx)
            self.assertEqual(sb1.type, sb2.type)
        # Now walk through, deleting each of the valence terms and make sure
        # that gets rid of the corresponding parm sections
        del parm.bonds[:], parm.bond_types[:]
        parm.remake_parm()
        self.assertNotIn('AMOEBA_REGULAR_BOND_NUM_PARAMS', parm.parm_data)
        self.assertNotIn('AMOEBA_REGULAR_BOND_FORCE_CONSTANT', parm.parm_data)
        self.assertNotIn('AMOEBA_REGULAR_BOND_EQUIL_VALUE', parm.parm_data)
        self.assertNotIn('AMOEBA_REGULAR_BOND_FTAB_DEGREE', parm.parm_data)
        self.assertNotIn('AMOEBA_REGULAR_BOND_FTAB_COEFFS', parm.parm_data)
        self.assertNotIn('AMOEBA_REGULAR_BOND_NUM_LIST', parm.parm_data)
        self.assertNotIn('AMOEBA_REGULAR_BOND_LIST', parm.parm_data)
        del parm.angles[:], parm.angle_types[:]
        parm.remake_parm()
        self.assertNotIn('AMOEBA_REGULAR_ANGLE_NUM_PARAMS', parm.parm_data)
        self.assertNotIn('AMOEBA_REGULAR_ANGLE_FORCE_CONSTANT', parm.parm_data)
        self.assertNotIn('AMOEBA_REGULAR_ANGLE_EQUIL_VALUE', parm.parm_data)
        self.assertNotIn('AMOEBA_REGULAR_ANGLE_FTAB_DEGREE', parm.parm_data)
        self.assertNotIn('AMOEBA_REGULAR_ANGLE_FTAB_COEFFS', parm.parm_data)
        self.assertNotIn('AMOEBA_REGULAR_ANGLE_NUM_LIST', parm.parm_data)
        self.assertNotIn('AMOEBA_REGULAR_ANGLE_LIST', parm.parm_data)
        del parm.stretch_bends[:], parm.stretch_bend_types[:]
        parm.remake_parm()
        self.assertNotIn('AMOEBA_STRETCH_BEND_FORCE_CONSTANT', parm.parm_data)
        self.assertNotIn('AMOEBA_STRETCH_BEND_FORCE_CONSTANT_1', parm.parm_data)
        self.assertNotIn('AMOEBA_STRETCH_BEND_FORCE_CONSTANT_2', parm.parm_data)
        self.assertNotIn('AMOEBA_STRETCH_BEND_BOND1_EQUIL_VALUE', parm.parm_data)
        self.assertNotIn('AMOEBA_STRETCH_BEND_BOND2_EQUIL_VALUE', parm.parm_data)
        self.assertNotIn('AMOEBA_STRETCH_BEND_ANGLE_EQUIL_VALUE', parm.parm_data)
        self.assertNotIn('AMOEBA_STRETCH_BEND_NUM_LIST', parm.parm_data)
        self.assertNotIn('AMOEBA_STRETCH_BEND_LIST', parm.parm_data)
        del parm.multipole_frames[:], parm.chiral_frames[:]
        parm.remake_parm()
        self.assertNotIn('AMOEBA_CHIRAL_FRAME_NUM_LIST', parm.parm_data)
        self.assertNotIn('AMOEBA_CHIRAL_FRAME_LIST', parm.parm_data)
        self.assertNotIn('AMOEBA_FRAME_DEF_NUM_LIST', parm.parm_data)
        self.assertNotIn('AMOEBA_FRAME_DEF_LIST', parm.parm_data)
        del parm.urey_bradleys[:], parm.urey_bradley_types[:]
        parm.remake_parm()
        self.assertNotIn('AMOEBA_UREY_BRADLEY_BOND_NUM_PARAMS', parm.parm_data)
        self.assertNotIn('AMOEBA_UREY_BRADLEY_BOND_FORCE_CONSTANT', parm.parm_data)
        self.assertNotIn('AMOEBA_UREY_BRADLEY_BOND_EQUIL_VALUE', parm.parm_data)
        self.assertNotIn('AMOEBA_UREY_BRADLEY_BOND_FTAB_DEGREE', parm.parm_data)
        self.assertNotIn('AMOEBA_UREY_BRADLEY_BOND_FTAB_COEFFS', parm.parm_data)
        self.assertNotIn('AMOEBA_UREY_BRADLEY_BOND_NUM_LIST', parm.parm_data)
        self.assertNotIn('AMOEBA_UREY_BRADLEY_BOND_LIST', parm.parm_data)

    def test_beeman_restart(self):
        """ Tests the BeemanRestart class """
        rst = load_file(get_fn('formbox_amoeba.rst'))
        rst2 = load_file(get_fn('nma.rst'))

        self.assertEqual(rst.natom, rst.parm_data['ATOMIC_COORDS_NUM_LIST'][0])
        self.assertEqual(rst2.natom, rst2.parm_data['ATOMIC_COORDS_NUM_LIST'][0])

        np.testing.assert_equal(rst.coordinates.flatten(), rst.parm_data['ATOMIC_COORDS_LIST'])
        np.testing.assert_equal(rst2.coordinates.flatten(), rst2.parm_data['ATOMIC_COORDS_LIST'])

        # Make the number of atoms in each case a little bit bigger
        oldcrds = rst.coordinates, rst2.coordinates
        rst.natom = rst.natom + 20
        rst2.natom = rst2.natom + 20

        np.testing.assert_equal(oldcrds[0].flatten().tolist() + [0 for i in range(60)],
                                rst.coordinates.flatten())
        np.testing.assert_equal(oldcrds[1].flatten().tolist() + [0 for i in range(60)],
                                rst2.coordinates.flatten())

        # Set coordinates
        rst.coordinates = np.arange(rst.natom*3, dtype=np.float64)
        np.testing.assert_equal(rst.coordinates, np.arange(rst.natom*3).reshape((1, rst.natom, 3)))
        np.testing.assert_equal(np.arange(rst.natom*3), rst.parm_data['ATOMIC_COORDS_LIST'])

        # Make sure accelerations, and old_accelerations are not present in the
        # restart that does not contain them
        self.assertRaises(AttributeError, lambda: rst2.velocities)
        self.assertRaises(AttributeError, lambda: rst2.accelerations)
        self.assertRaises(AttributeError, lambda: rst2.old_accelerations)

        # Make sure they do exist in the restart that *does* have them
        np.testing.assert_equal(rst.velocities.flatten(),
                                rst.parm_data['ATOMIC_VELOCITIES_LIST'])
        np.testing.assert_equal(rst.accelerations.flatten(),
                                rst.parm_data['ATOMIC_ACCELERATIONS_LIST'])
        np.testing.assert_equal(rst.old_accelerations.flatten(),
                                rst.parm_data['OLD_ATOMIC_ACCELERATIONS_LIST'])

        # Set velocities, accelerations, and old_accelerations
        rst2.velocities = np.random.rand(1, rst2.natom, 3)
        rst2.accelerations = np.random.rand(1, rst2.natom, 3)
        rst2.old_accelerations = np.random.rand(1, rst2.natom, 3)
        np.testing.assert_equal(rst2.velocities.flatten(),
                                rst2.parm_data['ATOMIC_VELOCITIES_LIST'])
        np.testing.assert_equal(rst2.accelerations.flatten(),
                                rst2.parm_data['ATOMIC_ACCELERATIONS_LIST'])
        np.testing.assert_equal(rst2.old_accelerations.flatten(),
                                rst2.parm_data['OLD_ATOMIC_ACCELERATIONS_LIST'])

        # Reduce number of atoms
        rst = load_file(get_fn('formbox_amoeba.rst'))
        rst_bak = load_file(get_fn('formbox_amoeba.rst'))
        rst2 = load_file(get_fn('nma.rst'))
        rst2_bak = load_file(get_fn('nma.rst'))

        rst.natom = rst.natom - 20
        np.testing.assert_equal(rst_bak.coordinates[:,:-20,:], rst.coordinates)
        np.testing.assert_equal(rst_bak.velocities[:,:-20,:], rst.velocities)
        np.testing.assert_equal(rst_bak.accelerations[:,:-20,:], rst.accelerations)
        np.testing.assert_equal(rst_bak.old_accelerations[:,:-20,:], rst.old_accelerations)

        rst2.natom = rst2.natom - 20
        np.testing.assert_equal(rst2_bak.coordinates[:,:-20,:], rst2.coordinates)
        self.assertRaises(AttributeError, lambda: rst2.velocities)
        self.assertRaises(AttributeError, lambda: rst2.accelerations)
        self.assertRaises(AttributeError, lambda: rst2.old_accelerations)

        # Check the error checking
        def bad_setting(trial):
            if trial == 1:
                rst.coordinates = [1]
            elif trial == 2:
                rst.velocities = [1]
            elif trial == 3:
                rst.accelerations = [1]
            elif trial == 4:
                rst.old_accelerations = [1]
            elif trial == 5:
                rst.box = [1]

        self.assertRaises(ValueError, lambda: bad_setting(1))
        self.assertRaises(ValueError, lambda: bad_setting(2))
        self.assertRaises(ValueError, lambda: bad_setting(3))
        self.assertRaises(ValueError, lambda: bad_setting(4))
        self.assertRaises(ValueError, lambda: bad_setting(5))

        # Overwrite existing vels
        rst.velocities = np.random.rand(1, rst.natom, 3)
        rst.accelerations = np.random.rand(1, rst.natom, 3)
        rst.old_accelerations = np.random.rand(1, rst.natom, 3)
        rst.box = [1, 2, 3, 90, 90, 90]
        np.testing.assert_equal(rst.velocities.flatten(),
                                rst.parm_data['ATOMIC_VELOCITIES_LIST'])
        np.testing.assert_equal(rst.accelerations.flatten(),
                                rst.parm_data['ATOMIC_ACCELERATIONS_LIST'])
        np.testing.assert_equal(rst.old_accelerations.flatten(),
                                rst.parm_data['OLD_ATOMIC_ACCELERATIONS_LIST'])
        np.testing.assert_equal(rst.box, rst.parm_data['UNIT_CELL_PARAMETERS'])

    def test_1012(self):
        """ Test that 10-12 prmtop files are recognized properly """
        parm = readparm.AmberParm(get_fn('ff91.parm7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        self._standard_parm_tests(parm, has1012=True)

    def test_getitem_nbtables(self):
        """ Tests that nonbond tables are set correctly following slicing """
        parm = load_file(get_fn('ash.parm7'), get_fn('ash.rst7'))
        new = parm[list(range(len(parm.atoms)))]
        self.assertEqual(len(parm.parm_data['NONBONDED_PARM_INDEX']),
                         len(new.parm_data['NONBONDED_PARM_INDEX']))
        self.assertEqual(len(parm.parm_data['LENNARD_JONES_ACOEF']),
                         len(new.parm_data['LENNARD_JONES_ACOEF']))
        self.assertEqual(len(parm.parm_data['LENNARD_JONES_BCOEF']),
                         len(new.parm_data['LENNARD_JONES_BCOEF']))

    def test_bad_arguments(self):
        """ Test proper error handling for bad AmberParm arguments """
        with self.assertRaises(TypeError):
            readparm.AmberParm(get_fn('trx.prmtop'), xyz=get_fn('trx.prmtop'))
        struct = create_random_structure(True)
        struct.unknown_functional = True
        self.assertRaises(TypeError, lambda: readparm.AmberParm.from_structure(struct))
        try:
            readparm.AmberParm.from_structure(struct)
        except TypeError as e:
            self.assertEqual('Cannot instantiate an AmberParm from unknown functional', str(e))
        struct.unknown_functional = False
        self.assertRaises(TypeError, lambda: readparm.AmberParm.from_structure(struct))
        try:
            readparm.AmberParm.from_structure(struct)
        except TypeError as e:
            self.assertTrue(str(e).startswith('AmberParm does not support all'))
        try:
            readparm.ChamberParm.from_structure(struct)
        except TypeError as e:
            self.assertTrue(str(e).startswith('ChamberParm does not support all'))
        with self.assertRaises(TypeError):
            readparm.AmberParm.from_structure(load_file(get_fn('ala_ala_ala.parm7')))
        try:
            readparm.AmberParm.from_structure(load_file(get_fn('ala_ala_ala.parm7')))
        except TypeError as e:
            self.assertTrue(str(e).endswith('Try ChamberParm'))

    def test_corrupt_parms(self):
        """ Test proper error detection on bad topology files """
        # First test proper error handling of a truncated section
        parm = readparm.AmberFormat(get_fn('ash.parm7'))
        atom_names = parm.parm_data['ATOM_NAME'][:]
        del parm.parm_data['ATOM_NAME'][0]
        self.assertRaises(AmberError, lambda: readparm.AmberParm.from_rawdata(parm))
        self.assertRaises(AmberError, lambda: readparm.AmoebaParm(get_fn('ash.parm7')))
        parm.parm_data['ATOM_NAME'] = atom_names
        parm.add_flag('AMOEBA_FORCEFIELD', '1I10', data=[0])
        self.assertRaises(AmberError, lambda: readparm.AmoebaParm.from_rawdata(parm))
        self.assertRaises(AmberError, lambda: parm.add_flag('RESIDUE_LABEL', '20a4', num_items=10))
        with self.assertRaises(IndexError):
            parm.add_flag('NEW_FLAG', '20a4', num_items=10, after='NO_FLAG')

    def test_amber_parm_box_xyz_args(self):
        """ Test passing coord and box arrays to AmberParm """
        crds = load_file(get_fn('ff14ipq.rst7'))
        parm = load_file(get_fn('ff14ipq.parm7'), xyz=crds.coordinates, box=crds.box)
        np.testing.assert_allclose(crds.coordinates[0], parm.coordinates)
        np.testing.assert_allclose(crds.box, parm.box)
        parm = readparm.AmberParm(get_fn('ff14ipq.parm7'), xyz=crds.coordinates, box=crds.box)
        np.testing.assert_allclose(crds.coordinates[0], parm.coordinates)
        np.testing.assert_allclose(crds.box, parm.box)
        # Also check .from_rawdata
        parm.velocities = np.random.rand(3, len(parm.atoms))
        parm2 = readparm.AmberParm.from_rawdata(parm)
        for a1, a2 in zip(parm.atoms, parm2.atoms):
            equal_atoms(self, a1, a2)
        np.testing.assert_allclose(parm2.coordinates, parm.coordinates)
        np.testing.assert_allclose(parm2.velocities, parm.velocities)
        np.testing.assert_allclose(parm2.box, parm.box)

    @unittest.skipUnless(HAS_GROMACS, 'Cannot test without Gromacs')
    def test_amber_parm_from_structure(self):
        """ Tests AmberParm.from_structure """
        aparm = load_file(get_fn('ash.parm7'), get_fn('ash.rst7'))
        cparm = load_file(get_fn('ala_ala_ala.parm7'), get_fn('ala_ala_ala.rst7'))
        acopy = readparm.AmberParm.from_structure(aparm, copy=True)
        ccopy = readparm.ChamberParm.from_structure(cparm, copy=True)
        anocopy = readparm.AmberParm.from_structure(aparm, copy=False)
        cnocopy = readparm.ChamberParm.from_structure(cparm, copy=False)

        self.assertEqual(len(aparm.atoms), len(acopy.atoms))
        self.assertEqual(len(aparm.bonds), len(acopy.bonds))
        self.assertEqual(len(cparm.atoms), len(ccopy.atoms))
        self.assertEqual(len(cparm.bonds), len(ccopy.bonds))
        self.assertEqual(len(aparm.atoms), len(anocopy.atoms))
        self.assertEqual(len(aparm.bonds), len(anocopy.bonds))
        self.assertEqual(len(cparm.atoms), len(cnocopy.atoms))
        self.assertEqual(len(cparm.bonds), len(cnocopy.bonds))

        self.assertIsNot(aparm.atoms, acopy.atoms)
        self.assertIsNot(aparm.bonds, acopy.bonds)
        self.assertIsNot(aparm.angles, acopy.angles)
        self.assertIsNot(aparm.dihedrals, acopy.dihedrals)
        self.assertIsNot(aparm.coordinates, acopy.coordinates)

        self.assertIs(aparm.atoms, anocopy.atoms)
        self.assertIs(aparm.bonds, anocopy.bonds)
        self.assertIs(aparm.angles, anocopy.angles)
        self.assertIs(aparm.dihedrals, anocopy.dihedrals)

        self.assertIsNot(cparm.atoms, ccopy.atoms)
        self.assertIsNot(cparm.bonds, ccopy.bonds)
        self.assertIsNot(cparm.angles, ccopy.angles)
        self.assertIsNot(cparm.dihedrals, ccopy.dihedrals)
        self.assertIsNot(cparm.coordinates, ccopy.coordinates)

        self.assertIs(cparm.atoms, cnocopy.atoms)
        self.assertIs(cparm.bonds, cnocopy.bonds)
        self.assertIs(cparm.angles, cnocopy.angles)
        self.assertIs(cparm.dihedrals, cnocopy.dihedrals)

        # Check that non-orthogonal unit cell conversions are properly handled
        tmp = load_file(get_fn(os.path.join('02.6water', 'topol.top')),
                        xyz=get_fn(os.path.join('02.6water', 'conf.gro')))
        tmp.box = [3, 3, 3, 109, 109, 90]
        parm = readparm.AmberParm.from_structure(tmp)
        np.testing.assert_equal(parm.box, tmp.box)
        self.assertEqual(parm.ptr('ifbox'), 3)
        self.assertEqual(parm.parm_data['BOX_DIMENSIONS'], [109, 3, 3, 3])

        # Check that a loaded structure without periodicities is properly warned
        # against
        tmp = load_file(get_fn(os.path.join('04.Ala', 'topol.top')),
                        xyz=get_fn(os.path.join('04.Ala', 'conf.gro')))
        # Test that a periodicity of zero is properly handled by ParmEd and
        # converted to a dummy term with a force constant of 0 (and ensure it
        # warns)
        tmp.dihedral_types[0][0].per = 0
        with self.assertWarns(AmberWarning):
            readparm.ChamberParm.from_structure(tmp)
        parm = readparm.AmberParm.from_structure(tmp)
        self.assertEqual(parm.dihedral_types[0].per, 0)
        self.assertAlmostEqual(parm.dihedral_types[0].phi_k, 0.27)
        self.assertEqual(parm.dihedral_types[0].phase, tmp.dihedral_types[0][0].phase)
        parm = readparm.ChamberParm.from_structure(tmp)
        self.assertEqual(parm.dihedral_types[0].per, 1)
        self.assertEqual(parm.dihedral_types[0].phi_k, 0)
        self.assertEqual(parm.dihedral_types[0].phase, tmp.dihedral_types[0][0].phase)
        def assign_box_badly():
            readparm.AmberParm(get_fn('ash.parm7')).box = [1, 2, 3, 4]
        self.assertRaises(ValueError, assign_box_badly)

    def test_old_parm_format(self):
        """ Test reading old Amber prmtop file format """
        self.assertTrue(readparm.AmberParm.id_format(get_fn('old.prmtop')))
        parm = load_file(get_fn('old.prmtop'), get_fn('old.inpcrd'))
        self.assertIsInstance(parm, readparm.AmberParm)
        self._standard_parm_tests(parm)

        # Now check the "slow" reader
        parm = readparm.AmberFormat()
        parm.rdparm(self.get_fn('old.prmtop'), slow=True)
        parm = parm.view_as(readparm.AmberParm)
        self._standard_parm_tests(parm)
        self.assertRaises(TypeError, lambda: readparm.AmberFormat().rdparm_slow(object()))

    # Tests for individual prmtops
    def _standard_parm_tests(self, parm, has1012=False):
        self.assertEqual(parm.ptr('natom'), len(parm.atoms))
        self.assertEqual(parm.ptr('nres'), len(parm.residues))
        self.assertEqual(parm.ptr('nbonh'), len(list(parm.bonds_inc_h)))
        self.assertEqual(parm.ptr('nbona'), len(list(parm.bonds_without_h)))
        self.assertEqual(parm.ptr('ntheth'), len(list(parm.angles_inc_h)))
        self.assertEqual(parm.ptr('ntheta'), len(list(parm.angles_without_h)))
        self.assertEqual(parm.ptr('nphih'), len(list(parm.dihedrals_inc_h)))
        self.assertEqual(parm.ptr('nphia'), len(list(parm.dihedrals_without_h)))
        self.assertEqual([a.name for a in parm.atoms],
                         parm.parm_data['ATOM_NAME'])
        self.assertEqual([a.type for a in parm.atoms],
                         parm.parm_data['AMBER_ATOM_TYPE'])
        if has1012:
            self.assertTrue(parm.has_1012())
        else:
            self.assertFalse(parm.has_1012())

    def _solv_pointer_tests(self, parm):
        self.assertEqual(parm.ptr('nspm'),
                         parm.parm_data['SOLVENT_POINTERS'][1])
        self.assertEqual(parm.ptr('nspm'),
                         len(parm.parm_data['ATOMS_PER_MOLECULE']))
        self.assertEqual(parm.ptr('natom'),
                         sum(parm.parm_data['ATOMS_PER_MOLECULE']))

    def _extensive_checks(self, parm):
        # Check the __contains__ methods of the various topologyobjects
        atoms = parm.atoms
        for bond in parm.bonds:
            self.assertEqual(sum([a in bond for a in atoms]), 2)
        for angle in parm.angles:
            self.assertEqual(sum([a in angle for a in atoms]), 3)
            self.assertEqual(sum([b in angle for b in parm.bonds]), 2)
        for dihedral in parm.dihedrals:
            self.assertEqual(sum([a in dihedral for a in atoms]), 4)
            self.assertEqual(sum([b in dihedral for b in parm.bonds]), 3)
        for residue in parm.residues:
            self.assertTrue(all([a in residue for a in residue.atoms]))
            self.assertEqual(sum([a in residue for a in atoms]),
                             len(residue))
        if not parm.chamber: return
        # Chamber tests now
        for ub in parm.urey_bradleys:
            self.assertEqual(sum([a in ub for a in atoms]), 2)
            self.assertEqual(sum([b in ub for b in parm.bonds]), 2)
        for imp in parm.impropers:
            self.assertEqual(sum([a in imp for a in atoms]), 4)
            self.assertEqual(sum([b in imp for b in parm.bonds]), 3)
        if parm.has_cmap:
            for cmap in parm.cmaps:
                self.assertEqual(sum([a in cmap for a in atoms]), 5)
                self.assertEqual(sum([b in cmap for b in parm.bonds]), 4)

def _num_unique_types(dct):
    return len(set(id(item) for _, item in iteritems(dct)))

def _num_unique_dtypes(dct):
    used_types = set()
    num = 0
    for _, x in iteritems(dct):
        if id(x) in used_types: continue
        used_types.add(id(x))
        num += len(x)
    return num

class TestParameterFiles(FileIOTestCase):
    """ Tests Amber parameter and frcmod files """

    @unittest.skipIf(os.getenv('AMBERHOME') is None, 'Cannot test w/out Amber')
    def test_find_amber_files(self):
        """ Tests the Amber file finder helper function """
        finder = parameters._find_amber_file
        self.assertEqual(finder(__file__, False), __file__)
        self.assertRaises(FileNotFoundError, lambda: finder('nofile', False))
        # Check looking in oldff
        self.assertRaises(FileNotFoundError, lambda: finder('rna.amberua.lib', False))
        self.assertEqual(finder('rna.amberua.lib', True),
            os.path.join(os.getenv('AMBERHOME'), 'dat', 'leap', 'lib', 'oldff', 'rna.amberua.lib')
        )

    def test_file_detection_frcmod(self):
        """ Tests the detection of Amber frcmod files """
        for fname in glob.glob(os.path.join(get_fn('parm'), 'frcmod.*')):
            self.assertTrue(parameters.AmberParameterSet.id_format(fname))
        # Now try creating a bunch of non-frcmod files to test the file ID
        # discrimination
        fn = self.get_fn('test.frcmod', written=True)
        with open(fn, 'w') as f:
            f.write('\n\n\n\n\n\n')
        self.assertFalse(parameters.AmberParameterSet.id_format(fn))

    def test_file_detection_parm(self):
        """ Tests the detection of Amber parm.dat files """
        for fname in glob.glob(os.path.join(get_fn('parm'), 'parm*.dat')):
            self.assertTrue(parameters.AmberParameterSet.id_format(fname))
        # Now try a bunch of slightly-off files to test discrimination
        for fname in glob.glob(os.path.join(get_fn('noparm'), '*')):
            self.assertFalse(parameters.AmberParameterSet.id_format(fname))

    def test_frcmod_parsing(self):
        """ Tests parsing an Amber frcmod file """
        self._check_ff99sb(
                parameters.AmberParameterSet(
                    os.path.join(get_fn('parm'), 'frcmod.ff99SB')
                )
        )
        self._check_ff99sb(
                parameters.AmberParameterSet(
                    os.path.join(get_fn('parm'), 'frcmod.1')
                )
        )

    def test_parm_dat_bad_equivalencing(self):
        """ Test handling of erroneous atom equivalencing in parm.dat files """
        with self.assertWarns(AmberWarning):
            parameters.AmberParameterSet(
                    os.path.join(get_fn('parm'), 'parmAM1.dat')
            )
        params = parameters.AmberParameterSet(os.path.join(get_fn('parm'), 'parmAM1.dat'))
        # Make sure CA and C have different types, even though they are
        # explicitly equivalenced
        self.assertEqual(params.atom_types['C'].rmin, 1.9127)
        self.assertEqual(params.atom_types['C'].epsilon, 0.086)
        self.assertEqual(params.atom_types['CA'].rmin, 1.9061)
        self.assertEqual(params.atom_types['CA'].epsilon, 0.086)

    def test_frcmod_with_tabstops(self):
        """ Test parsing an Amber frcmod file with tabs instead of spaces """
        params = parameters.AmberParameterSet(
                os.path.join(get_fn('parm'), 'all_modrna08.frcmod')
        )
        self.assertEqual(len(params.atom_types), 38) # Ugh! Duplicates??  Really??
        self.assertEqual(params.bond_types[('C', 'CM')],
                         topologyobjects.BondType(449.9, 1.466)) # OVERWRITING IN THE SAME FILE??

    def test_nonconsecutive_torsions(self):
        """ Test proper warning of non-consecutive multi-term dihedrals """
        with self.assertWarns(ParameterWarning):
            parameters.AmberParameterSet(os.path.join(get_fn('parm'), 'parm14ipq.dat'))

    def _check_ff99sb(self, params):
        self.assertEqual(_num_unique_types(params.atom_types), 0)
        self.assertEqual(_num_unique_types(params.bond_types), 0)
        self.assertEqual(_num_unique_types(params.angle_types), 0)
        self.assertEqual(_num_unique_types(params.dihedral_types), 4)
        self.assertEqual(_num_unique_dtypes(params.dihedral_types), 16)
        self.assertEqual(_num_unique_types(params.improper_periodic_types), 0)
        # Small enough to check all of the parameters
        self.assertEqual(params.dihedral_types[('C','N','CT','C')][0],
                         topologyobjects.DihedralType(0, 4, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('C','N','CT','C')][1],
                         topologyobjects.DihedralType(0.42, 3, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('C','N','CT','C')][2],
                         topologyobjects.DihedralType(0.27, 2, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('C','N','CT','C')][3],
                         topologyobjects.DihedralType(0, 1, 0, 1.2, 2.0))

        self.assertEqual(params.dihedral_types[('N','CT','C','N')][0],
                         topologyobjects.DihedralType(0, 4, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('N','CT','C','N')][1],
                         topologyobjects.DihedralType(0.55, 3, 180, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('N','CT','C','N')][2],
                         topologyobjects.DihedralType(1.58, 2, 180, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('N','CT','C','N')][3],
                         topologyobjects.DihedralType(0.45, 1, 180, 1.2, 2.0))

        self.assertEqual(params.dihedral_types[('CT','CT','N','C')][0],
                         topologyobjects.DihedralType(0, 4, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('CT','CT','N','C')][1],
                         topologyobjects.DihedralType(0.4, 3, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('CT','CT','N','C')][2],
                         topologyobjects.DihedralType(2.0, 2, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('CT','CT','N','C')][3],
                         topologyobjects.DihedralType(2.0, 1, 0, 1.2, 2.0))

        self.assertEqual(params.dihedral_types[('CT','CT','C','N')][0],
                         topologyobjects.DihedralType(0, 4, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('CT','CT','C','N')][1],
                         topologyobjects.DihedralType(0.4, 3, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('CT','CT','C','N')][2],
                         topologyobjects.DihedralType(0.2, 2, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('CT','CT','C','N')][3],
                         topologyobjects.DihedralType(0.2, 1, 0, 1.2, 2.0))

    def test_parm_parsing(self):
        """ Tests parsing an Amber parm.dat file """
        params = parameters.AmberParameterSet(
                os.path.join(get_fn('parm'), 'parm10.dat')
        )
        self.assertEqual(_num_unique_types(params.atom_types), 63)
        self.assertEqual(_num_unique_types(params.bond_types), 151)
        self.assertEqual(_num_unique_types(params.angle_types), 400)
        self.assertEqual(_num_unique_dtypes(params.dihedral_types), 275)
        self.assertEqual(_num_unique_types(params.improper_periodic_types), 59)
        # Check a couple random atom types
        self.assertEqual(params.atom_types['C'].mass, 12.01)
        self.assertEqual(params.atom_types['C'].atomic_number, 6)
        self.assertEqual(params.atom_types['H3'].mass, 1.008)
        self.assertEqual(params.atom_types['H3'].atomic_number, 1)
        self.assertEqual(params.atom_types['EP'].atomic_number, 0)
        self.assertEqual(params.atom_types['EP'].mass, 0)
        self.assertEqual(params.atom_types['N*'].mass, 14.01)
        self.assertEqual(params.atom_types['N*'].atomic_number, 7)
        # Check a couple random bond types
        self.assertEqual(params.bond_types[('OW', 'HW')].req, 0.9572)
        self.assertEqual(params.bond_types[('OW', 'HW')].k, 553)
        self.assertEqual(params.bond_types[('C', 'C')].req, 1.525)
        self.assertEqual(params.bond_types[('C', 'C')].k, 310)
        self.assertEqual(params.bond_types[('OH', 'P')].req, 1.61)
        self.assertEqual(params.bond_types[('OH', 'P')].k, 230)
        self.assertEqual(params.bond_types[('C4', 'N*')].req, 1.365)
        self.assertEqual(params.bond_types[('C4', 'N*')].k, 448)
        # Check a couple random angle types
        self.assertEqual(params.angle_types[('HW', 'OW', 'HW')].theteq, 104.52)
        self.assertEqual(params.angle_types[('HW', 'OW', 'HW')].k, 100)
        self.assertEqual(params.angle_types[('CC', 'NA', 'P')].theteq, 125.10)
        self.assertEqual(params.angle_types[('CC', 'NA', 'P')].k, 76.7)
        self.assertEqual(params.angle_types[('HA', 'CM', 'CT')].theteq, 120.00)
        self.assertEqual(params.angle_types[('HA', 'CM', 'CT')].k, 50.0)
        self.assertEqual(params.angle_types[('C4', 'N*', 'CT')].theteq, 121.2)
        self.assertEqual(params.angle_types[('C4', 'N*', 'CT')].k, 70)
        self.assertEqual(params.angle_types[('EP', 'S', 'S')].theteq, 96.70)
        self.assertEqual(params.angle_types[('EP', 'S', 'S')].k, 150)
        # Check a couple random dihedral types
        d = params.dihedral_types[('X', 'C', 'C', 'X')]
        self.assertEqual(len(d), 1)
        self.assertEqual(d[0].phi_k, 14.5/4)
        self.assertEqual(d[0].per, 2)
        self.assertEqual(d[0].phase, 180)
        d = params.dihedral_types[('CT', 'OS', 'CT', 'CI')]
        self.assertEqual(len(d), 2)
        self.assertEqual(d[0].phi_k, 0.383)
        self.assertEqual(d[0].per, 3)
        self.assertEqual(d[0].phase, 0)
        self.assertEqual(d[1].phi_k, 0.1)
        self.assertEqual(d[1].per, 2)
        self.assertEqual(d[1].phase, 180)
        d = params.dihedral_types[('N', 'CT', 'CT', 'OH')]
        self.assertEqual(len(d), 4)
        self.assertEqual(d[0].phi_k, 0)
        self.assertEqual(d[0].per, 1)
        self.assertEqual(d[0].phase, 0)
        self.assertEqual(d[1].phi_k, 1.49)
        self.assertEqual(d[1].per, 2)
        self.assertEqual(d[1].phase, 0)
        self.assertEqual(d[2].phi_k, 0.156)
        self.assertEqual(d[2].per, 3)
        self.assertEqual(d[2].phase, 0)
        self.assertEqual(d[3].phi_k, 0)
        self.assertEqual(d[3].per, 4)
        self.assertEqual(d[3].phase, 0)
        d = params.dihedral_types[('EP', 'S', 'S', 'EP')]
        self.assertEqual(len(d), 1)
        self.assertEqual(d[0].phi_k, 0)
        self.assertEqual(d[0].per, 3)
        self.assertEqual(d[0].phase, 0)
        # Check a couple random improper types
        d = params.improper_periodic_types[('O', 'X', 'C', 'X')]
        self.assertEqual(d.phi_k, 10.5)
        self.assertEqual(d.per, 2)
        self.assertEqual(d.phase, 180)
        d = params.improper_periodic_types[('CB', 'CK', 'N*', 'CT')]
        self.assertEqual(d.phi_k, 1.0)
        self.assertEqual(d.per, 2)
        self.assertEqual(d.phase, 180)
        d = params.improper_periodic_types[('CC', 'CR', 'NA', 'P')]
        self.assertEqual(d.phi_k, 1.1)
        self.assertEqual(d.per, 2)
        self.assertEqual(d.phase, 180)
        # Check some of the L-J parameters
        self.assertEqual(params.atom_types['H'].rmin, 0.6)
        self.assertEqual(params.atom_types['H'].epsilon, 0.0157)
        self.assertEqual(params.atom_types['N'].rmin, 1.824)
        self.assertEqual(params.atom_types['N'].epsilon, 0.17)
        self.assertEqual(params.atom_types['I'].rmin, 2.35)
        self.assertEqual(params.atom_types['I'].epsilon, 0.4)
        self.assertEqual(params.atom_types['C*'].rmin, 1.908)
        self.assertEqual(params.atom_types['C*'].epsilon, 0.086)
        # Now check some of the equivalenced atom types
        self.assertEqual(params.atom_types['NA'].rmin, 1.824)
        self.assertEqual(params.atom_types['NA'].epsilon, 0.17)
        self.assertEqual(params.atom_types['NY'].rmin, 1.824)
        self.assertEqual(params.atom_types['NY'].epsilon, 0.17)
        self.assertEqual(params.atom_types['CA'].rmin, 1.908)
        self.assertEqual(params.atom_types['CA'].epsilon, 0.086)
        self.assertEqual(params.atom_types['CP'].rmin, 1.908)
        self.assertEqual(params.atom_types['CP'].epsilon, 0.086)

    def test_parm_parsing_ljedit(self):
        """ Tests parsing an Amber parm.dat file with an LJEDIT section """
        params = parameters.AmberParameterSet(os.path.join(get_fn('parm'), 'parm14ipq.dat'))
        self.assertEqual(_num_unique_types(params.atom_types), 74)
        self.assertEqual(_num_unique_types(params.bond_types), 217)
        self.assertEqual(_num_unique_types(params.angle_types), 724)
        self.assertEqual(_num_unique_dtypes(params.dihedral_types), 1848)
        self.assertEqual(_num_unique_types(params.improper_periodic_types), 102)
        self.assertEqual(_num_unique_types(params.nbfix_types), 6)
        # Check a couple dihedral types, since this file has them disordered
        d = params.dihedral_types[('TN', 'TG', 'C', 'N')]
        self.assertEqual(len(d), 4)
        self.assertEqual(d[0].phi_k, 0.04031)
        self.assertEqual(d[0].per, 4)
        self.assertEqual(d[0].phase, 0)
        self.assertEqual(d[1].phi_k, 0.06853)
        self.assertEqual(d[1].per, 3)
        self.assertEqual(d[1].phase, 180)
        self.assertEqual(d[2].phi_k, 0.19829)
        self.assertEqual(d[2].per, 2)
        self.assertEqual(d[2].phase, 180)
        self.assertEqual(d[3].phi_k, 1.46258)
        self.assertEqual(d[3].per, 1)
        self.assertEqual(d[3].phase, 180)
        # Check the nbfix types
        self.assertEqual(params.nbfix_types[('O3', 'OW')][0], math.sqrt(0.162750*0.21))
        self.assertEqual(params.nbfix_types[('O3', 'OW')][1], 1.775931+1.8605)
        self.assertEqual(params.nbfix_types[('OA', 'OW')][0], math.sqrt(0.162750*0.2104))
        self.assertEqual(params.nbfix_types[('OA', 'OW')][1], 1.775931+1.66)
        # Check inside an frcmod file
        params = parameters.AmberParameterSet(os.path.join(get_fn('parm'), 'frcmod.2'))
        self.assertEqual(_num_unique_types(params.atom_types), 4)
        self.assertEqual(_num_unique_types(params.bond_types), 0)
        self.assertEqual(_num_unique_types(params.angle_types), 0)
        self.assertEqual(_num_unique_types(params.dihedral_types), 0)
        self.assertEqual(params.nbfix_types[('HC', 'OH')][0], math.sqrt(0.0150*0.2))
        self.assertEqual(params.nbfix_types[('HC', 'OH')][1], 1.377+1.721)

    @unittest.skipIf(os.getenv('AMBERHOME') is None, 'Cannot test w/out Amber')
    def test_load_leaprc(self):
        """ Tests loading a leaprc file to define a force field """
        fn = os.path.join(os.getenv('AMBERHOME'), 'dat', 'leap', 'cmd', 'leaprc.protein.ff14SB')
        frcmod_fn = os.path.join(os.getenv('AMBERHOME'), 'dat', 'leap', 'parm', 'frcmod.tip4p')
        params = parameters.AmberParameterSet.from_leaprc(fn)
        params.load_parameters(frcmod_fn)
        self.assertEqual(params.atom_types['H'].atomic_number, 1)
        self.assertEqual(params.atom_types['3C'].atomic_number, 6)
        self.assertEqual(params.atom_types['EP'].atomic_number, 0)
        self.assertTrue(params.residues)
        with open(fn) as f, open(frcmod_fn) as frcmod:
            params = parameters.AmberParameterSet.from_leaprc(f)
            params.load_parameters(frcmod)
        self.assertEqual(params.atom_types['H'].atomic_number, 1)
        self.assertEqual(params.atom_types['3C'].atomic_number, 6)
        self.assertEqual(params.atom_types['EP'].atomic_number, 0)
        self.assertTrue(params.residues)
        # Now make sure it accepts search_oldff=True
        params = parameters.AmberParameterSet.from_leaprc(fn, search_oldff=True)
        params.load_parameters(frcmod_fn)
        self.assertEqual(params.atom_types['H'].atomic_number, 1)
        self.assertEqual(params.atom_types['3C'].atomic_number, 6)
        self.assertEqual(params.atom_types['EP'].atomic_number, 0)
        self.assertTrue(params.residues)

    def test_load_leaprc_filenames_with_spaces(self):
        """ Tests loading a leaprc file with filenames containing spaces """
        fn1 = self.get_fn('leaprc', written=True)
        fn2 = self.get_fn('amino 12.lib', written=True)
        fn3 = os.path.join(get_fn('parm'), 'parm10.dat')
        fn4 = self.get_fn('parm 10.dat', written=True)
        with open(fn1, 'w') as f:
            f.write('loadOFF "%s"\n' % fn2)
            f.write('loadAmberParams "%s"\n' % fn4)
        shutil.copy(get_fn('amino12.lib'), fn2)
        shutil.copy(fn3, fn4)
        params = parameters.AmberParameterSet.from_leaprc(fn1)
        with open(fn1, 'w') as f:
            f.write('loadOFF %s\n' % fn2.replace(' ', r'\ '))
            f.write('loadAmberParams %s\n' % fn4.replace(' ', r'\ '))
        params = parameters.AmberParameterSet.from_leaprc(fn1)

    def test_load_leaprc_with_mol2(self):
        """ Tests loading a leaprc file with loadMol2 files """
        fn1 = self.get_fn('leaprc', written=True)
        with open(fn1, 'w') as f:
            f.write('DAN = loadMol2 "%s"\n' % get_fn('tripos1.mol2'))
            f.write('GPN = loadMol3 "%s"\n' % get_fn('tripos9.mol2'))
        params = parameters.AmberParameterSet.from_leaprc(fn1)
        self.assertEqual(len(params.residues), 2)
        self.assertIn('DAN', params.residues)
        self.assertIn('GPN', params.residues)
        # Now make sure we warn about mult-residue mol2 files
        with open(fn1, 'w') as f:
            f.write('SOME = loadMol2 "%s"\n' % get_fn('multimol.mol2'))
        with self.assertWarns(AmberWarning):
            parameters.AmberParameterSet.from_leaprc(fn1)
        params = parameters.AmberParameterSet.from_leaprc(fn1)
        self.assertEqual(len(params.residues), 200)
        self.assertIn('ZINC00000016_1', params.residues)

    def test_parm_set_parsing(self):
        """ Tests parsing a set of Amber parameter files """
        params = parameters.AmberParameterSet(
            os.path.join(get_fn('parm'), 'parm99.dat'),
            [
                os.path.join(get_fn('parm'), 'frcmod.ff99SB'),
                os.path.join(get_fn('parm'), 'frcmod.parmbsc0'),
            ],
        )
        self._check_paramset(params)
        params = parameters.AmberParameterSet(
            open(os.path.join(get_fn('parm'), 'parm99.dat')),
            open(os.path.join(get_fn('parm'), 'frcmod.ff99SB')),
            open(os.path.join(get_fn('parm'), 'frcmod.parmbsc0')),
        )
        self._check_paramset(params)

    def _check_paramset(self, params):
        self.assertGreater(_num_unique_types(params.atom_types), 0)
        self.assertGreater(_num_unique_types(params.bond_types), 0)
        self.assertGreater(_num_unique_types(params.angle_types), 0)
        self.assertGreater(_num_unique_types(params.dihedral_types), 0)
        self.assertGreater(_num_unique_types(params.improper_periodic_types), 0)
        # Check that parameters were properly overridden. parm99.dat defines
        # C-N-CT-C torsion as follows:
        #
        # C -N -CT-C    1    0.850       180.000          -2.
        # C -N -CT-C    1    0.800         0.000           1.
        #
        # whereas ff99SB defines that torsion as:
        #
        # C -N -CT-C    1    0.00          0.0            -4.
        # C -N -CT-C    1    0.42          0.0            -3.
        # C -N -CT-C    1    0.27          0.0            -2.
        # C -N -CT-C    1    0.00          0.0             1.
        #
        # Since ff99SB is loaded last, this should be the one that is stored
        self.assertEqual(len(params.dihedral_types[('C','N','CT','C')]), 4)
        self.assertEqual(params.dihedral_types[('C','N','CT','C')][0],
                         topologyobjects.DihedralType(0, 4, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('C','N','CT','C')][1],
                         topologyobjects.DihedralType(0.42, 3, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('C','N','CT','C')][2],
                         topologyobjects.DihedralType(0.27, 2, 0, 1.2, 2.0))
        self.assertEqual(params.dihedral_types[('C','N','CT','C')][3],
                         topologyobjects.DihedralType(0, 1, 0, 1.2, 2.0))

    def test_private_functions(self):
        """ Test some of the private utility functions in AmberParameterSet """
        # Check on AmberParm._truncate_array
        parm = readparm.AmberParm(get_fn('ash.parm7'), get_fn('ash.rst7'))
        parm._truncate_array('ATOM_NAME', 2)
        self.assertEqual(len(parm.parm_data['ATOM_NAME']), 2)

    def test_load_lib(self):
        """ Test parsing Amber .lib files within a set of parameter files """
        params = parameters.AmberParameterSet(
                os.path.join(get_fn('parm'), 'parm10.dat'),
                os.path.join(get_fn('parm'), 'frcmod.ff14SB'),
                get_fn('amino12.lib'),
                get_fn('aminoct12.lib')
        )
        self.assertTrue(params.atom_types)
        self.assertTrue(params.bond_types)
        self.assertTrue(params.angle_types)
        self.assertTrue(params.dihedral_types)
        self.assertTrue(params.improper_periodic_types)
        self.assertTrue(params.residues)
        params = parameters.AmberParameterSet(
            [
                os.path.join(get_fn('parm'), 'parm10.dat'),
                os.path.join(get_fn('parm'), 'frcmod.ff14SB'),
            ],
            [
                get_fn('amino12.lib'),
                get_fn('aminoct12.lib'),
            ]
        )
        self.assertTrue(params.atom_types)
        self.assertTrue(params.bond_types)
        self.assertTrue(params.angle_types)
        self.assertTrue(params.dihedral_types)
        self.assertTrue(params.improper_periodic_types)
        self.assertTrue(params.residues)

    @unittest.skipIf(os.getenv('AMBERHOME') is None, 'Cannot test w/out Amber')
    def test_load_lib_with_blank_lines(self):
        """ Tests parsing of .lib files with blank lines """
        fn = os.path.join(os.getenv('AMBERHOME'), 'dat', 'leap', 'lib', 'all_aminoAM1.lib')
        self.assertTrue(AmberOFFLibrary.id_format(fn))
        lib = AmberOFFLibrary.parse(fn)
        self.assertEqual(len(lib), 27)
        self.assertEqual(lib['ALA'].atoms[1].name, 'H')
        self.assertEqual(lib['ALA'].atoms[1].charge, 0.423221)
        self.assertEqual(lib['VAL'].atoms[8].name, 'HG12')
        self.assertEqual(lib['VAL'].atoms[8].charge, 0.062124)

    def test_lib_with_box(self):
        """ Tests handling of OFF files with multiple residues and a box """
        solvents = AmberOFFLibrary.parse(get_fn('solvents.lib'))
        for res in solvents['TIP3PBOX']:
            self.assertEqual(len(res.bonds), 3)
            for atom in res.atoms:
                self.assertEqual(len(atom.bonds), 2)
        self.assertIsNot(solvents['TIP3PBOX'].box, None)
        # Now create one
        struct = readparm.LoadParm(get_fn('cyclohexane.parm7'), get_fn('cyclohexane.md.rst7'))
        self.assertIsNot(struct.box, None)
        reslib = ResidueTemplateContainer.from_structure(struct)
        reslib.name = 'CYCHBOX'
        self.assertIsNot(reslib.box, None)
        np.testing.assert_equal(struct.box, reslib.box)
        # Now write an OFF file
        fn = self.get_fn('test.lib', written=True)
        AmberOFFLibrary.write(dict(CYCHBOX=reslib), fn)
        # Now read it and make sure I have the appropriate bonds
        lib2 = AmberOFFLibrary.parse(fn)
        # All residues should have exactly the same number of bonds
        self.assertEqual(len({len(x.bonds) for x in lib2['CYCHBOX']}), 1)
        self.assertGreater(len(lib2['CYCHBOX'][0].bonds), 0)
        np.testing.assert_allclose(lib2['CYCHBOX'].box, struct.box, atol=0.001)

    @unittest.skipIf(os.getenv('AMBERHOME') is None, 'Cannot test w/out Amber')
    def test_lib_without_residueconnect(self):
        """ Test parsing OFF library files without RESIDUECONNECT """
        fn = os.path.join(os.getenv('AMBERHOME'), 'dat', 'leap', 'lib',
                          'lipid14.lib')
        self.assertTrue(AmberOFFLibrary.id_format(fn))
        lib = AmberOFFLibrary.parse(fn)
        self.assertIs(lib['LA'].head, lib['LA'].tail) # weird...
        # Nucleic acid caps
        fn = os.path.join(os.getenv('AMBERHOME'), 'dat', 'leap', 'lib',
                          'cph_nucleic_caps.lib')
        self.assertTrue(AmberOFFLibrary.id_format(fn))
        lib = AmberOFFLibrary.parse(fn)
        self.assertEqual(len(lib), 2)

    def test_glycam_parsing(self):
        """ Tests reading GLYCAM parameter files (weird dihedrals) """
        fn = os.path.join(get_fn('parm'), 'GLYCAM_06j.dat')
        params = parameters.AmberParameterSet(fn)
        self.assertEqual(len(params.dihedral_types[('Oy', 'Cy', 'Os', 'CT')]), 3)

class TestCoordinateFiles(FileIOTestCase):
    """ Tests the various coordinate file classes """

    def test_mdcrd(self):
        """ Test the ASCII trajectory file parsing """
        mdcrd = asciicrd.AmberMdcrd(get_fn('tz2.truncoct.crd'),
                                    natom=5827, hasbox=True, mode='r')
        self.assertEqual(mdcrd.frame, 10)
        self.assertIsInstance(mdcrd.coordinates, np.ndarray)

        runsum = 0
        for i in range(10):
            arr1 = mdcrd.coordinates[i]
            runsum += arr1.sum()
        self.assertAlmostEqual(runsum, 7049.817, places=3)

        pos = mdcrd.positions
        self.assertEqual(len(pos), 5827)
        self.assertTrue(u.is_quantity(pos))
        for (px, py, pz), (x, y, z) in zip(pos, mdcrd.coordinates[0]):
            self.assertEqual(px, x*u.angstroms)
            self.assertEqual(py, y*u.angstroms)
            self.assertEqual(pz, z*u.angstroms)

    def test_error_handling(self):
        """ Tests error handling of Amber ASCII coordinate files """
        self.assertRaises(ValueError, lambda:
                asciicrd.AmberMdcrd(get_fn('tz2.truncoct.crd'), mode='x',
                                    natom=5827, hasbox=True)
        )
        self.assertRaises(NotImplementedError, lambda:
                asciicrd._AmberAsciiCoordinateFile(get_fn('tz2.truncoct.crd'),
                            mode='r', natom=5827, hasbox=True)
        )
        class NewType(asciicrd._AmberAsciiCoordinateFile):
            DEFAULT_TITLE = 'default'
            CRDS_PER_LINE = 6
        self.assertRaises(NotImplementedError, lambda:
                NewType(get_fn('tz2.truncoct.crd'), mode='r', natom=5827,
                        hasbox=True)
        )

    def test_restart(self):
        """ Test the ASCII restart file parsing """
        restart = asciicrd.AmberAsciiRestart(get_fn('tz2.ortho.rst7'), 'r')
        self.assertEqual(restart.natom, 5293)
        self.assertTrue(restart.hasbox)
        self.assertFalse(restart.hasvels)
        self.assertIsInstance(restart.coordinates, np.ndarray)
        crdsum = restart.coordinates.sum()
        self.assertAlmostEqual(crdsum, 301623.26028240257, places=4)

        pos = restart.positions
        self.assertEqual(len(pos), 5293)
        self.assertTrue(u.is_quantity(pos))
        for (px, py, pz), (x, y, z) in zip(pos, restart.coordinates[0]):
            self.assertEqual(px, x*u.angstroms)
            self.assertEqual(py, y*u.angstroms)
            self.assertEqual(pz, z*u.angstroms)

    def test_write_restart(self):
        """ Test writing Amber restart files """
        self._check_restarts_with_atoms(10)
        self._check_restarts_with_atoms(11)

    def test_restart_error_handling(self):
        """ Test Amber ASCII restart file error handling """
        fn = self.get_fn('test_file', written=True)
        with open(fn, 'w') as f:
            f.write('Some arbitrary title\n')
            f.write('%5d%15e7\n' % (10,10))
            for j in range(5):
                for i in range(6):
                    f.write('%12.7f' % (random.random()*10-5))
                f.write('\n')
            f.write('\n\n\n\n\n\n\ntoo many lines!\n')
        self.assertRaises(RuntimeError, lambda: load_file(fn))

        # Try illegal stuff with coordinate setting
        rst = asciicrd.AmberAsciiRestart(fn, 'w', natom=10)
        crd = np.random.rand(10, 3)
        self.assertRaises(RuntimeError, lambda: rst.coordinates)
        def set_bad_num_crd():
            rst.coordinates = np.random.rand(20, 3)
        self.assertRaises(ValueError, set_bad_num_crd)
        rst.coordinates = crd
        np.testing.assert_equal(rst.coordinates.squeeze(), crd)
        def write_crd_twice():
            rst.coordinates = np.random.rand(10, 3)
        self.assertRaises(RuntimeError, write_crd_twice)
        rst.close()
        def set_crd_on_old_file():
            load_file(fn).coordinates = np.random.rand(10, 3)
        self.assertRaises(RuntimeError, set_crd_on_old_file)

        # Try illegal stuff with velocities
        rst = asciicrd.AmberAsciiRestart(fn, 'w', natom=10)
        crd = np.random.rand(10, 3)
        vel = np.random.rand(10, 3)
        self.assertRaises(RuntimeError, lambda: rst.velocities)
        def set_vel_early():
            rst.velocities = vel
        self.assertRaises(RuntimeError, set_vel_early)
        rst.coordinates = crd
        np.testing.assert_equal(rst.coordinates.squeeze(), crd)
        def set_bad_num_vel():
            rst.velocities = np.random.rand(20, 3)
        self.assertRaises(ValueError, set_bad_num_vel)
        rst.velocities = vel
        np.testing.assert_equal(rst.velocities.squeeze(), vel)
        def write_vel_twice():
            rst.velocities = np.random.rand(10, 3)
        self.assertRaises(RuntimeError, write_vel_twice)
        rst.close()
        def set_vel_on_old_file():
            load_file(fn).velocities = np.random.rand(10, 3)
        self.assertRaises(RuntimeError, set_vel_on_old_file)

    def _check_restarts_with_atoms(self, natom):
        # Write a file with coordinates. Then read it back in and make sure
        # everything is OK
        fn = self.get_fn('test.rst7', written=True)
        restart = asciicrd.AmberAsciiRestart(fn, 'w', natom=natom, title='nose', time=100)
        crd = np.random.rand(natom, 3) * 5 - 10
        restart.coordinates = crd
        restart.close()
        check = load_file(fn)
        np.testing.assert_allclose(crd, check.coordinates.squeeze())
        self.assertIs(check.velocities, None)
        self.assertIs(check.box, None)

        # Write a file with coordinates and velocities. Then read it back in and
        # make sure everything is OK
        restart = asciicrd.AmberAsciiRestart(fn, 'w', natom=natom, title='nose', time=100)
        crd = np.random.rand(natom, 3) * 5 - 10
        vel = np.random.rand(natom, 3) * 5 - 10
        restart.coordinates = crd
        restart.velocities = vel
        restart.close()
        check = load_file(fn)
        np.testing.assert_allclose(crd, check.coordinates.squeeze(), atol=1e-4)
        np.testing.assert_allclose(vel, check.velocities.squeeze(), atol=1e-4)
        self.assertIs(check.box, None)

        # Write a file with coordinates and box. Then read it back in and
        # make sure everything is OK
        restart = asciicrd.AmberAsciiRestart(fn, 'w', natom=natom, title='nose', time=100)
        crd = np.random.rand(natom, 3) * 5 - 10
        box = [20, 20, 20, 90, 90, 90]
        restart.coordinates = crd
        restart.box = box
        restart.close()
        check = load_file(fn)
        np.testing.assert_allclose(crd, check.coordinates.squeeze(), atol=1e-4)
        np.testing.assert_equal(check.box, [20, 20, 20, 90, 90, 90])
        self.assertIs(check.velocities, None)

        # Write a file with coordinates and velocities. Then read it back in and
        # make sure everything is OK
        restart = asciicrd.AmberAsciiRestart(fn, 'w', natom=natom, title='nose', time=100)
        crd = np.random.rand(natom, 3) * 5 - 10
        vel = np.random.rand(natom, 3) * 5 - 10
        box = [20, 20, 20, 90, 90, 90]
        restart.coordinates = crd
        restart.velocities = vel
        restart.box = box
        restart.close()
        check = load_file(fn)
        np.testing.assert_allclose(crd, check.coordinates.squeeze(), atol=1e-4)
        np.testing.assert_allclose(vel, check.velocities.squeeze(), atol=1e-4)
        np.testing.assert_equal(box, check.box)

    def test_auto_detection(self):
        """ Tests ASCII coordinate file autodetections """
        fn = self.get_fn('test_file', written=True)
        with open(fn, 'w') as f:
            f.write('Some arbitrary title\n')
            f.write('%5d\n' % -1)
            for i in range(6):
                f.write('%12.7f' % random.random())
            f.write('\n')
        self.assertFalse(asciicrd.AmberAsciiRestart.id_format(fn))

        with open(fn, 'w') as f:
            f.write('Some arbitrary title\n')
            f.write('%5d%15e7\n' % (10,10))
            for j in range(5):
                for i in range(6):
                    f.write('%12.7f' % (random.random()*10-5))
                f.write('\n')
        self.assertTrue(asciicrd.AmberAsciiRestart.id_format(fn))

        with open(fn, 'w') as f:
            f.write('Some arbitrary title\n')
            f.write('%5d' % 10)
            for j in range(5):
                for i in range(6):
                    f.write('%12.7f' % (random.random()*10-5))
                f.write('\n')
        self.assertFalse(asciicrd.AmberAsciiRestart.id_format(fn))

        with open(fn, 'w') as f:
            f.write('Some arbitrary title\n')
            f.write('%5d\n' % 10)
            for j in range(5):
                for i in range(6):
                    f.write('%12.7f' % (random.random()*10-5))
                f.write('\n')
            f.write('\n\n\n\n')
        self.assertTrue(asciicrd.AmberAsciiRestart.id_format(fn))
        rst = load_file(fn)
        self.assertEqual(rst.natom, 10)
        self.assertEqual(rst.time, 0)
        self.assertTrue(
                np.all(np.logical_and(-5<=rst.coordinates, rst.coordinates<=5))
        )

        with open(fn, 'w') as f:
            f.write('Some arbitrary title\n')
            f.write('%5d\n' % 10)
            for j in range(5):
                for i in range(5):
                    f.write('%12.7f' % (random.random()*10-5))
                f.write('   123456790')
                f.write('\n')
        self.assertFalse(asciicrd.AmberAsciiRestart.id_format(fn))

        with open(fn, 'w') as f:
            f.write('Some arbitrary title\n')
            f.write('%5d\n' % 10)
            for j in range(5):
                f.write('   1.345679 ')
                for i in range(5):
                    f.write('%12.7f' % (random.random()*10-5))
                f.write('\n')
        self.assertFalse(asciicrd.AmberAsciiRestart.id_format(fn))

        with open(fn, 'w') as f:
            f.write('Some arbitrary title\n')
            f.write('%5d\n' % 10)
            for j in range(5):
                for i in range(4):
                    f.write('%12.7f' % (random.random()*10-5))
                f.write('   1.345679 ')
                f.write('%12.7f' % (random.random()*10-5))
                f.write('\n')
        self.assertFalse(asciicrd.AmberAsciiRestart.id_format(fn))

        with open(fn, 'w') as f:
            f.write('Some arbitrary title\n')
            f.write('%5d\n' % 10)
            for j in range(5):
                for i in range(5):
                    f.write('%12.7f' % (random.random()*10-5))
                f.write('   1.345a79 ')
                f.write('\n')
        self.assertFalse(asciicrd.AmberAsciiRestart.id_format(fn))

        with open(fn, 'w') as f:
            f.write('Only one line\n')
        self.assertFalse(asciicrd.AmberAsciiRestart.id_format(fn))

    @unittest.skipIf(PYPY, 'Test does not yet run under pypy')
    def test_parsing_netcdf_file_without_unitcells(self):
        """Parsing Amber netcdf file that does not have box info"""
        fn, tn = get_fn('tz2.nc'), get_fn('tz2.parm7')
        parm = pmd.load_file(tn, fn)
        self.assertEqual(parm.get_coordinates().shape, (101, 223, 3))
        self.assertEqual(len(parm.atoms), 223)
        self.assertIs(parm.box, None)


class TestAmberMask(unittest.TestCase):
    """ Test the Amber mask parser """

    def test_mask(self):
        """ Test the Amber mask parser """
        parm = readparm.AmberParm(get_fn('trx.prmtop'))
        mask_res1 = mask.AmberMask(parm, ':1')
        mask_res2 = mask.AmberMask(parm, ':1:2')
        mask_res3 = mask.AmberMask(parm, ':1-3')
        mask_resala = mask.AmberMask(parm, ':ALA')
        mask_carbonat = mask.AmberMask(parm, '@/C')
        mask_atnum = mask.AmberMask(parm, '@1-10')
        mask_atnum2 = mask.AmberMask(parm, '@1-10,21-30')
        mask_atname = mask.AmberMask(parm, '@CA')
        mask_atname2 = mask.AmberMask(parm, '@CA@CB')
        mask_atname3 = mask.AmberMask(parm, '@CA,1,CB,2-3')
        mask_resat = mask.AmberMask(parm, ':ALA@CA')
        mask_attyp = mask.AmberMask(parm, '@%CT')
        mask_wldcrd = mask.AmberMask(parm, '@H=')
        mask_wldcrd2 = mask.AmberMask(parm, '*&(@H*)')
        mask_wldcrd3 = mask.AmberMask(parm, ':AL=')
        mask_wldcrd4 = mask.AmberMask(parm, ':A=A')
        mask_wldcrd5 = mask.AmberMask(parm, ':1,*')
        mask_wldcrd6 = mask.AmberMask(parm, '@*')
        mask_wldcrd7 = mask.AmberMask(parm, '@CA|@*')
        mult_op = mask.AmberMask(parm, '(!:1&@CA)|!:2')
        mask_res12 = mask.AmberMask(parm, ':1 | :2 ')

        # Check all of the masks
        self.assertEqual(sum(mask_res1.Selection(prnlev=9)), 13)
        for idx in mask_res1.Selected():
            self.assertEqual(parm.atoms[idx].residue.idx, 0)
        self.assertEqual(list(range(13)), list(mask_res1.Selected()))
        sel = mask_res1.Selection()
        for atom in parm.atoms:
            if atom.residue.idx == 0:
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)
        self.assertEqual(sum(mask_res1.Selection(invert=True))+13,
                         len(parm.atoms))

        self.assertEqual(sum(mask_res2.Selection()), len(parm.residues[0])+
                                                     len(parm.residues[1]))
        _ = (0, 1)
        for idx in mask_res1.Selected():
            self.assertIn(parm.atoms[idx].residue.idx, _)

        for atom, sel in zip(parm.atoms, mask_res3.Selection()):
            if atom.residue.idx < 3:
                self.assertEqual(sel, 1)
            else:
                self.assertEqual(sel, 0)

        self.assertEqual(sum(mask_resala.Selection()), 121)
        for idx in mask_resala.Selected():
            self.assertEqual(parm.atoms[idx].residue.name, 'ALA')
        sel = mask_resala.Selection()
        for atom in parm.atoms:
            if atom.residue.name == 'ALA':
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)

        for idx in mask_carbonat.Selected():
            self.assertEqual(parm.atoms[idx].atomic_number, 6)
        sel = mask_carbonat.Selection()
        for atom in parm.atoms:
            if atom.atomic_number == 6:
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)

        self.assertEqual(sum(mask_atnum.Selection()), 10)
        self.assertEqual(sum(mask_atnum2.Selection()), 20)

        self.assertEqual(sum(mask_atname.Selection()), 108)
        for idx in mask_atname.Selected():
            self.assertEqual(parm.atoms[idx].name, 'CA')
        sel = mask_atname.Selection()
        for atom in parm.atoms:
            if atom.name == 'CA':
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)

        sel = mask_atname2.Selection()
        for atom in parm.atoms:
            if atom.name == 'CA' or atom.name == 'CB':
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)

        for atom, val in zip(parm.atoms, mask_atname3.Selection()):
            if atom.name in ('CA', 'CB'):
                self.assertEqual(val, 1)
            elif atom.idx < 3:
                self.assertEqual(val, 1)
            else:
                self.assertEqual(val, 0)

        self.assertEqual(sum(mask_resat.Selection()), 12)
        for idx in mask_resat.Selected():
            self.assertEqual(parm.atoms[idx].name, 'CA')
            self.assertEqual(parm.atoms[idx].residue.name, 'ALA')
        sel = mask_resat.Selection()
        for atom in parm.atoms:
            if atom.residue.name == 'ALA' and atom.name == 'CA':
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)

        self.assertEqual(sum(mask_attyp.Selection()), 341)
        for idx in mask_attyp.Selected():
            self.assertEqual(parm.atoms[idx].type, 'CT')
        sel = mask_attyp.Selection()
        for atom in parm.atoms:
            if atom.type == 'CT':
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)

        for atom, sel in zip(parm.atoms, mask_wldcrd.Selection()):
            if sel:
                self.assertTrue(atom.name.startswith('H'))
            else:
                self.assertFalse(atom.name.startswith('H'))

        for atom, sel in zip(parm.atoms, mask_wldcrd2.Selection()):
            if sel:
                self.assertTrue(atom.name.startswith('H'))
            else:
                self.assertFalse(atom.name.startswith('H'))

        for atom, sel in zip(parm.atoms, mask_wldcrd3.Selection()):
            if sel:
                self.assertTrue(atom.residue.name.startswith('AL'))
            else:
                self.assertFalse(atom.residue.name.startswith('AL'))

        for atom, sel in zip(parm.atoms, mask_wldcrd4.Selection()):
            if sel:
                self.assertEqual(atom.residue.name[0], 'A')
                self.assertEqual(atom.residue.name[2], 'A')
            else:
                if len(atom.residue.name) >= 3:
                    self.assertFalse(atom.residue.name[0] == 'A' and
                                     atom.residue.name[2] == 'A')

        for sel in mask_wldcrd5.Selection():
            self.assertEqual(sel, 1)

        for sel in mask_wldcrd6.Selection():
            self.assertEqual(sel, 1)

        for sel in mask_wldcrd7.Selection():
            self.assertEqual(sel, 1)

        for atom, val in zip(parm.atoms, mult_op.Selection()):
            if atom.residue.idx == 1:
                if atom.name == 'CA':
                    self.assertEqual(val, 1)
                else:
                    self.assertEqual(val, 0)
            else:
                self.assertEqual(val, 1)

        self.assertEqual({parm.atoms[i].residue.idx for i in mask_res12.Selected()}, {0, 1})

    def test_illegal_masks(self):
        """ Test bad mask strings """
        parm = load_file(get_fn('ash.parm7'), get_fn('ash.rst7'))
        parm_nocor = load_file(get_fn('ash.parm7'))

        mask1 = mask.AmberMask(parm, '(@1)=')
        mask2 = mask.AmberMask(parm_nocor, ':2<@5')
        mask3 = mask.AmberMask(parm, '@1<10')
        mask4 = mask.AmberMask(parm, 'C')
        mask5 = mask.AmberMask(parm, 'C&@1')
        mask6 = mask.AmberMask(parm, ':(T=HIS)')
        mask7 = mask.AmberMask(parm, ':1,%^')
        mask8 = mask.AmberMask(parm, '(\xc3)')
        mask9 = mask.AmberMask(parm, '@&:1')
        mask10 = mask.AmberMask(parm, '(((:10)&@CA)|:1')
        mask11 = mask.AmberMask(parm, '!!:1')
        mask12 = mask.AmberMask(parm, ':1<:10.a')
        mask13 = mask.AmberMask(parm, '(:A=A&@CA))')
        mask14 = mask.AmberMask(parm, ':ALA(&)')
        mask15 = mask.AmberMask(parm, '(:ALA<)|:3')
        mask16 = mask.AmberMask(parm, ':ALA<3')
        mask17 = mask.AmberMask(parm, ':ALA<')
        mask18 = mask.AmberMask(parm, '@1-20,40-xyz')
        mask19 = mask.AmberMask(parm, ':1-xyz')
        mask20 = mask.AmberMask(parm, ':1-4,40-xyz')
        mask21 = mask.AmberMask(parm, '@C%')
        mask22 = mask.AmberMask(parm, ':C%')
        mask23 = mask.AmberMask(parm, '@/Fk') # Looking for all Fakeiums

        self.assertRaises(MaskError, mask1.Selection)
        self.assertRaises(MaskError, mask2.Selection)
        self.assertRaises(MaskError, mask3.Selection)
        self.assertRaises(MaskError, mask4.Selection)
        self.assertRaises(MaskError, mask5.Selection)
        self.assertRaises(MaskError, mask6.Selection)
        self.assertRaises(MaskError, mask7.Selection)
        self.assertRaises(MaskError, mask8.Selection)
        self.assertRaises(MaskError, mask9.Selection)
        self.assertRaises(MaskError, mask10.Selection)
        self.assertRaises(MaskError, mask11.Selection)
        self.assertRaises(MaskError, mask12.Selection)
        self.assertRaises(MaskError, mask13.Selection)
        self.assertRaises(MaskError, mask14.Selection)
        self.assertRaises(MaskError, mask15.Selection)
        self.assertRaises(MaskError, mask16.Selection)
        self.assertRaises(MaskError, mask17.Selection)
        self.assertRaises(MaskError, mask18.Selection)
        self.assertRaises(MaskError, mask19.Selection)
        self.assertRaises(MaskError, mask20.Selection)
        self.assertRaises(MaskError, mask21.Selection)
        self.assertRaises(MaskError, mask22.Selection)
        self.assertRaises(MaskError, mask23.Selection)
        # These should never be triggered, but they are good safeguards
        self.assertRaises(MaskError, lambda:
                mask23._binop('/', mask._mask(len(parm.atoms)), mask._mask(len(parm.atoms))))
        self.assertRaises(MaskError, lambda: mask23._priority('/'))

        # Test the internal mask interface
        pmask = mask._mask(20)
        pmask2 = mask._mask(21)
        self.assertEqual(len(pmask), 20)
        self.assertEqual(pmask.pop(), 0)
        self.assertEqual(len(pmask), 20) # Pop does not change length
        self.assertRaises(MaskError, lambda: pmask.append(10)) # append disabled
        self.assertRaises(MaskError, lambda: pmask.extend([0,0]))
        self.assertRaises(MaskError, lambda: pmask.remove(0)) # remove disabled
        self.assertRaises(MaskError, lambda: pmask.And(pmask2))
        self.assertRaises(MaskError, lambda: pmask.Or(pmask2))

    def test_compound_mask(self):
        """ Tests compound/complex Amber selection masks """
        parm = readparm.AmberParm(get_fn('trx.prmtop'))
        mask1 = mask.AmberMask(parm, ':1-6@CA,C,O,N')
        mask2 = mask.AmberMask(parm, ':1-6,ALA,10@%CT | @O')
        mask3 = mask.AmberMask(parm, '(:1-9,ALA,GLY@CA,H*)|@%CT')

        sel = mask1.Selection()
        for atom in parm.atoms:
            if (atom.residue.idx < 6 and atom.name in ('CA', 'C', 'O', 'N')):
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)

        sel = mask2.Selection()
        for atom in parm.atoms:
            if (atom.residue.idx < 6 or atom.residue.name == 'ALA' or
                    atom.residue.idx == 9) and atom.type == 'CT':
                self.assertEqual(sel[atom.idx], 1)
            elif atom.name == 'O':
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)

        sel = mask3.Selection()
        for atom in parm.atoms:
            res = atom.residue
            if ((res.idx < 9 or res.name in ('ALA', 'GLY')) and
                    (atom.name == 'CA' or atom.name.startswith('H'))):
                self.assertEqual(sel[atom.idx], 1)
            elif atom.type == 'CT':
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)

    def test_distance_based_mask_pdb(self):
        """ Test distance-based mask selections on a PDB file """
        parm = load_file(get_fn('4lzt.pdb'))
        # All atoms within 5 A of residue 8
        mask1 = mask.AmberMask(parm, ':8<@5')
        sel = mask1.Selection()
        self.assertGreater(sum(sel), 0)
        for i, atom in enumerate(parm.atoms):
            for j, a2 in enumerate(parm.residues[7]):
                dx = atom.xx - a2.xx
                dy = atom.xy - a2.xy
                dz = atom.xz - a2.xz
                if dx*dx + dy*dy + dz*dz < 25:
                    self.assertTrue(sel[i])
                    break
            else:
                self.assertFalse(sel[i])

    def test_distance_based_mask(self):
        """ Test distance-based mask selections """
        parm = readparm.AmberParm(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        # All atoms within 5 A of residue 8
        mask1 = mask.AmberMask(parm, ':8<@5')
        # All atoms in any residue with at least one atom more than 10 A away from residue 1
        mask2 = mask.AmberMask(parm, ':1>:10')

        sel = mask1.Selection()
        for i, atom in enumerate(parm.atoms):
            within = 0
            for atom2 in parm.residues[7].atoms:
                dx = atom2.xx - atom.xx
                dy = atom2.xy - atom.xy
                dz = atom2.xz - atom.xz
                if (dx*dx + dy*dy + dz*dz) < 25:
                    within = 1
                    break
            self.assertEqual(sel[i], within)

        sel = mask2.Selection()
        self.assertEqual(sum(sel), 1588)
        for i, res in enumerate(parm.residues):
            within = 0
            for atom in res.atoms:
                for atom2 in parm.residues[0].atoms:
                    dx = atom2.xx - atom.xx
                    dy = atom2.xy - atom.xy
                    dz = atom2.xz - atom.xz
                    if (dx*dx + dy*dy + dz*dz) > 100:
                        within = 1
                        break
                if within: break
            for atom in res.atoms:
                self.assertEqual(sel[atom.idx], within)

    def test_mask_underscore(self):
        """ Test mask selection with atom name having an underscore """
        parm = readparm.AmberParm(get_fn('ash.parm7'))
        name = 'AT_A'
        change(parm, 'ATOM_NAME', '@1', name).execute()
        # Make sure a selection will grab this atom
        mask1 = mask.AmberMask(parm, '@%s' % name)
        self.assertEqual(list(mask1.Selected()), [0])

class TestWriteFiles(FileIOTestCase):

    def test_write_amber_parm(self):
        """ Test writing an AmberParm file """
        parm = readparm.AmberParm(self.get_fn('trx.prmtop'))
        parm.write_parm(self.get_fn('trx.prmtop', written=True))
        f1 = open(self.get_fn('trx.prmtop'), 'r')
        f2 = open(self.get_fn('trx.prmtop', written=True), 'r')
        try:
            for line1, line2 in zip(f1, f2):
                if line1.startswith('%VERSION'):
                    self.assertTrue(line2.startswith('%VERSION'))
                    continue
                self.assertEqual(line1.strip(), line2.strip())
        finally:
            f1.close()
            f2.close()

    def test_write_chamber_parm(self):
        """ Checks for correct units in improper phase in chamber prmtop """
        parm = readparm.ChamberParm(self.get_fn('test_fad.prmtop'))
        parm.write_parm(self.get_fn('test_fad.prmtop', written=True))
        self.assertTrue(
            diff_files(self.get_fn('test_fad.prmtop', written=True),
                       self.get_fn('test_fad.prmtop.save', saved=True),
                       absolute_error=1e-4)
        )

    def test_save_amber_parm(self):
        """ Test writing AmberParm file with AmberParm.save """
        parm = readparm.AmberParm(self.get_fn('trx.prmtop'))
        parm.add_flag('NEW_FLAG', '10I6', num_items=parm.ptr('nres'))
        self.assertIn('NEW_FLAG', parm.parm_data)
        self.assertIn('NEW_FLAG', parm.flag_list)
        parm.save(self.get_fn('trx.prmtop', written=True))
        parm2 = readparm.AmberParm(self.get_fn('trx.prmtop', written=True))
        self.assertIn('NEW_FLAG', parm2.parm_data)

    def test_write_pdb_with_LES_parm(self):
        """ Tests writing a PDB file with a parm created with LES in mind """
        output = StringIO()
        saved_pdb = self.get_fn('4lzt.les.pdb')
        rst7_name = self.get_fn('4lzt.les.rst7')
        parm_name = self.get_fn('4lzt.les.parm7')
        parm = pmd.load_file(parm_name, rst7_name)
        parm.write_pdb(output)
        output.seek(0)
        new_parm = pmd.read_PDB(output)
        saved_parm = pmd.load_file(saved_pdb)

        for new_atom, saved_atom  in zip(new_parm.atoms, saved_parm.atoms):
            if saved_atom.other_locations:
                self.assertTrue(new_atom.other_locations)

        # make sure the labels are added, only pick two atoms since we already
        # tested above
        output.seek(0)
        buffer = output.read()
        self.assertIn('ATOM     18  HE2ALYS     1      -0.780   9.159  10.504', buffer)
        self.assertIn('ATOM     23  NZ BLYS     1      -0.618   8.282   8.901', buffer)

    def test_amber_restart(self):
        """ Test writing an ASCII Amber restart file """
        Restart = asciicrd.AmberAsciiRestart
        box = [10, 10, 10, 90, 90, 90]
        rst = Restart(self.get_fn('testc.rst7', written=True), 'w', natom=9)
        rst.coordinates = list(range(27))
        rst.close()
        rst = Restart(self.get_fn('testcv.rst7', written=True), 'w', natom=20)
        rst.coordinates = list(range(60))
        rst.velocities = list(reversed(range(60)))
        rst.close()
        rst = Restart(self.get_fn('testcb.rst7', written=True), 'w', natom=7)
        rst.coordinates = list(range(21))
        rst.box = box[:]
        rst.close()
        rst = Restart(self.get_fn('testcvb.rst7', written=True), 'w', natom=15)
        rst.coordinates = list(range(45))
        rst.velocities = list(reversed(range(45)))
        self.assertRaises(RuntimeError, lambda: rst.box)
        rst.box = box[:]
        rst.close()
        self._check_written_restarts(box)

    def test_amber_restart_numpy(self):
        """ Test writing Amber restart file passing numpy arrays """
        Restart = asciicrd.AmberAsciiRestart
        box = np.asarray([10, 10, 10, 90, 90, 90])
        rst = Restart(self.get_fn('testc.rst7', written=True), 'w', natom=9)
        rst.coordinates = np.arange(27).reshape((9,3))
        rst.close()
        rst = Restart(self.get_fn('testcv.rst7', written=True), 'w', natom=20)
        rst.coordinates = np.arange(60).reshape((20,3))
        rst.velocities = np.asarray(list(reversed(range(60)))).reshape((20,3))
        rst.close()
        rst = Restart(self.get_fn('testcb.rst7', written=True), 'w', natom=7)
        rst.coordinates = np.arange(21).reshape((7,3))
        rst.box = box
        rst.close()
        rst = Restart(self.get_fn('testcvb.rst7', written=True), 'w', natom=15)
        rst.coordinates = np.arange(45).reshape((15,3))
        rst.velocities = np.asarray(list(reversed(range(45)))).reshape((15,3))
        rst.box = box
        rst.close()
        self._check_written_restarts(box)

    def test_amber_mdcrd(self):
        """ Test writing ASCII trajectory file """
        box = [15, 15, 15]
        Mdcrd = asciicrd.AmberMdcrd
        crd = Mdcrd(self.get_fn('testc.mdcrd', written=True), natom=15, hasbox=False,
                    mode='w', title='Test file')
        crd.add_coordinates(list(range(45)))
        crd.add_coordinates([x+1 for x in range(45)])
        crd.add_coordinates([x+2 for x in range(45)])
        crd.add_coordinates([x+3 for x in range(45)])
        crd.add_coordinates([x+4 for x in range(45)])
        crd.close()
        crd = Mdcrd(self.get_fn('testcb.mdcrd', written=True), natom=18, hasbox=True,
                    mode='w', title='Test file')
        crd.add_coordinates(list(range(54)))
        crd.add_box(box)
        crd.add_coordinates([x+1 for x in range(54)])
        crd.add_box(box)
        crd.add_coordinates([x+2 for x in range(54)])
        crd.add_box(box)
        crd.add_coordinates([x+3 for x in range(54)])
        crd.add_box(box)
        crd.add_coordinates([x+4 for x in range(54)])
        crd.add_box(box)
        crd.close()
        self._check_written_mdcrds(box)

    def test_amber_mdcrd_numpy(self):
        """ Test writing ASCII trajectory file passing numpy arrays """
        box = np.asarray([15, 15, 15])
        Mdcrd = asciicrd.AmberMdcrd
        crd = Mdcrd(self.get_fn('testc.mdcrd', written=True), natom=15, hasbox=False,
                    mode='w', title='Test file')
        coorddata = np.arange(45).reshape((15,3))
        crd.add_coordinates(coorddata)
        crd.add_coordinates(coorddata+1)
        crd.add_coordinates(coorddata+2)
        crd.add_coordinates(coorddata+3)
        crd.add_coordinates(coorddata+4)
        crd.close()
        crd = Mdcrd(self.get_fn('testcb.mdcrd', written=True), natom=18, hasbox=True,
                    mode='w', title='Test file')
        coorddata = np.arange(54).reshape((18,3))
        crd.add_coordinates(coorddata)
        crd.add_box(box)
        crd.add_coordinates(coorddata+1)
        crd.add_box(box)
        crd.add_coordinates(coorddata+2)
        crd.add_box(box)
        crd.add_coordinates(coorddata+3)
        crd.add_box(box)
        crd.add_coordinates(coorddata+4)
        crd.add_box(box)
        crd.close()
        self._check_written_mdcrds(box)

    def test_bad_file_usage(self):
        """ Check that illegal file usage results in desired exceptions """
        Restart = asciicrd.AmberAsciiRestart
        Mdcrd = asciicrd.AmberMdcrd
        box = [10, 10, 10, 90, 90, 90]
        rst = Restart(self.get_fn('testc.rst7', written=True), 'w', natom=9, hasbox=True)
        def assign(obj, stmnt):
            rst = crd = obj
            exec(stmnt)
        try:
            self.assertRaises(ValueError,
                              lambda: assign(rst, 'rst.coordinates=range(20)'))
            self.assertRaises(RuntimeError,
                              lambda: assign(rst, 'rst.box=[10]*3+[90]*3'))
            self.assertRaises(RuntimeError,
                              lambda: assign(rst, 'rst.cell_angles=[90]*3'))
            rst.coordinates = list(range(27))
            self.assertRaises(ValueError,
                              lambda: assign(rst, 'rst.cell_lengths=[0, 1]'))
            self.assertRaises(RuntimeError, lambda: rst.cell_lengths)
            self.assertRaises(RuntimeError, lambda: rst.cell_angles)
            rst.box = box
            self.assertRaises(RuntimeError, lambda:
                              assign(rst, 'rst.velocities=list(range(27))'))
            self.assertRaises(RuntimeError, lambda:
                              assign(rst, 'rst.cell_lengths=[1, 2, 3]'))
            self.assertRaises(RuntimeError, lambda:
                              assign(rst, 'rst.cell_angles=[1, 2, 3]'))
        finally:
            rst.close()
        try:
            rst = Restart(self.get_fn('testc.rst7', written=True), 'r')
            with self.assertRaises(RuntimeError):
                assign(rst, 'rst.cell_lengths=[1, 2, 3]')
            with self.assertRaises(RuntimeError):
                assign(rst, 'rst.cell_angles=[1, 2, 3]')
        finally:
            rst.close()
        crd = Mdcrd(self.get_fn('testc.mdcrd', written=True), natom=15, hasbox=True, mode='w', title='Test file')
        s = 'list(range(45))'
        s2 = 'list(range(42))'
        try:
            crd.add_coordinates(eval(s))
            self.assertRaises(RuntimeError,
                              lambda: assign(crd, 'crd.add_coordinates(%s)'%s))
            crd.add_box([10, 10, 10])
            self.assertRaises(ValueError,
                              lambda: assign(crd, 'crd.add_coordinates(%s)'%s2))
        finally:
            crd.close()

    def _check_written_restarts(self, box):
        # Now try to read them and verify the information (keep in mind that the
        # restart velocities are scaled down then back up, so you'll need to use
        # assertAlmostEqual in this case).
        rst = readparm.Rst7.open(self.get_fn('testc.rst7', written=True))
        self.assertFalse(rst.hasbox)
        self.assertIs(rst.box, None)
        self.assertFalse(rst.hasvels)
        np.testing.assert_equal(rst.coordinates.flatten(), list(range(27)))
        rst = asciicrd.AmberAsciiRestart(self.get_fn('testc.rst7', written=True))
        self.assertIs(rst.cell_lengths, None)
        self.assertIs(rst.cell_angles, None)
        rst = readparm.Rst7.open(self.get_fn('testcb.rst7', written=True))
        self.assertTrue(rst.hasbox)
        self.assertFalse(rst.hasvels)
        np.testing.assert_equal(rst.coordinates.flatten(), list(range(21)))
        np.testing.assert_equal(rst.box.flatten(), box)
        rst = readparm.Rst7.open(self.get_fn('testcv.rst7', written=True))
        self.assertTrue(rst.hasvels)
        self.assertFalse(rst.hasbox)
        np.testing.assert_equal(rst.coordinates, np.arange(60).reshape(rst.coordinates.shape))
        np.testing.assert_allclose(rst.velocities, np.array(list(reversed(range(60)))).reshape(rst.velocities.shape))
        rst = readparm.Rst7.open(self.get_fn('testcvb.rst7', written=True))
        self.assertTrue(rst.hasvels)
        self.assertTrue(rst.hasbox)
        np.testing.assert_equal(rst.coordinates,
                np.arange(45).reshape(rst.coordinates.shape))
        np.testing.assert_allclose(rst.velocities,
                np.array(list(reversed(range(45)))).reshape(rst.velocities.shape))
        np.testing.assert_equal(rst.box, box)

    def _check_written_mdcrds(self, box):
        # Now try to read them and verify the information
        crd = asciicrd.AmberMdcrd(self.get_fn('testc.mdcrd', written=True), 15, False, 'r')
        self.assertEqual(crd.title, 'Test file')
        self.assertFalse(crd.hasbox)
        for i in range(crd.frame):
            shape = crd.coordinates[i].shape
            refcrd = np.arange(45) + i
            np.testing.assert_equal(crd.coordinates[i], refcrd.reshape(shape))
        for i, array in enumerate(crd.coordinates):
            refcrd = (np.arange(45) + i).reshape(array.shape)
            np.testing.assert_equal(array, refcrd)
        crd.close()

        crd = asciicrd.AmberMdcrd(self.get_fn('testcb.mdcrd', written=True), 18, True, 'r')
        self.assertEqual(crd.title, 'Test file')
        self.assertTrue(crd.hasbox)
        for i in range(crd.frame):
            refcrd = (np.arange(54) + i).reshape(crd.coordinates[i].shape)
            np.testing.assert_equal(crd.coordinates[i], refcrd)
            np.testing.assert_equal(crd.box[i], box)
        for i, (coords, mybox) in enumerate(zip(crd.coordinates, crd.box)):
            np.testing.assert_equal(coords, (np.arange(54)+i).reshape(coords.shape))
            np.testing.assert_equal(mybox, box)

class TestObjectAPIs(unittest.TestCase):
    """ Tests various object APIs """

    def test_tracked_list(self):
        """ Tests the TrackedList object """
        mylist = topologyobjects.TrackedList(range(20))
        mylist2 = topologyobjects.TrackedList(reversed(range(20)))
        self.assertFalse(mylist.changed)
        self.assertIsInstance(mylist[0], int)
        self.assertIsInstance(mylist[0:5], list)
        self.assertIsInstance(mylist + mylist2, topologyobjects.TrackedList)
        self.assertIsInstance(mylist * 2, topologyobjects.TrackedList)
        mylist += mylist2
        self.assertTrue(mylist.changed)
        self.assertFalse(mylist2.changed)
        for i in range(20):
            self.assertEqual(i, mylist[i])
            self.assertEqual(i, mylist[39-i])
        mylist.changed = False
        mylist.append(10)
        self.assertTrue(mylist.changed)
        mylist.changed = False
        mylist.extend(mylist2)
        self.assertTrue(mylist.changed)
        self.assertIsInstance(mylist, topologyobjects.TrackedList)
        mylist.changed = False
        mylist.pop()
        self.assertTrue(mylist.changed)
        mylist.changed = False
        mylist[5] = 8
        self.assertTrue(mylist.changed)
        mylist.changed = False
        del mylist[20:]
        self.assertTrue(mylist.changed)
        mylist.changed = False
        del mylist[0]
        self.assertTrue(mylist.changed)
        mylist.changed = False
        mylist *= 2
        self.assertTrue(mylist.changed)

class TestAmberParmSlice(unittest.TestCase):
    """ Tests fancy slicing """

    def test_split(self):
        """ Tests the molecule splitting functionality """
        parm = readparm.AmberParm(get_fn('solv.prmtop'))
        parts = parm.split()
        # Make sure the sum of the parts is equal to the whole
        natom = sum(len(part[0].atoms)*len(part[1]) for part in parts)
        self.assertEqual(len(parm.atoms), natom)
        self.assertEqual(len(parts), 4) # 4 types of molecules
        self.assertEqual(len(parts[0][1]), 1)
        self.assertEqual(len(parts[1][1]), 1)
        self.assertEqual(len(parts[2][1]), 8)
        self.assertEqual(len(parts[3][1]), 9086)

    def test_split_2(self):
        """ Tests splitting distinct single-residue molecules with same name """
        parm = readparm.AmberParm(get_fn('phenol.prmtop'))
        self.assertEqual(len(parm.residues), 1)
        self.assertEqual(parm.residues[0].name, 'MOL')
        parm2 = readparm.AmberParm(get_fn('biphenyl.prmtop'))
        self.assertEqual(len(parm2.residues), 1)
        self.assertEqual(parm2.residues[0].name, 'MOL')

        comb = parm * 20 + parm2 * 30

        self.assertEqual(len(comb.residues), 50)

        parts = comb.split()
        self.assertEqual(len(parts), 2)
        self.assertEqual(len(parts[0][0].atoms), len(parm.atoms))
        self.assertEqual(len(parts[1][0].atoms), len(parm2.atoms))
        self.assertEqual(len(parts[0][1]), 20)
        self.assertEqual(len(parts[1][1]), 30)

    def test_add(self):
        """ Tests combining AmberParm instances """
        parm1 = readparm.AmberParm(get_fn('phenol.prmtop'))
        parm2 = readparm.AmberParm(get_fn('biphenyl.prmtop'))
        comb = parm1 + parm2
        self.assertEqual(len(comb.atoms), len(parm1.atoms) + len(parm2.atoms))
        for a1, a2 in zip(comb.atoms, parm1.atoms + parm2.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.solvent_radius, a2.solvent_radius)
        self.assertEqual(len(comb.residues), len(parm1.residues) + len(parm2.residues))
        for r1, r2 in zip(comb.residues, parm1.residues + parm2.residues):
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.chain, r2.chain)
        # In-place now
        parm1 += parm2
        self.assertEqual(len(parm1.atoms), len(comb.atoms))
        for a1, a2 in zip(comb.atoms, parm1.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.solvent_radius, a2.solvent_radius)
        self.assertEqual(len(parm1.residues), len(comb.residues))
        for r1, r2 in zip(comb.residues, parm1.residues):
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.chain, r2.chain)

    def test_mult(self):
        """ Tests replicating AmberParm instances """
        parm = readparm.AmberParm(get_fn('phenol.prmtop'))
        mult = parm * 5
        self.assertEqual(len(mult.atoms), 5*len(parm.atoms))
        self.assertEqual(len(mult.residues), 5*len(parm.residues))
        for i, a1 in enumerate(mult.atoms):
            a2 = parm[i%len(parm.atoms)]
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.solvent_radius, a2.solvent_radius)
        for i, r1 in enumerate(mult.residues):
            r2 = parm.residues[i%len(parm.residues)]
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.chain, r2.chain)
        # In-place now
        parm *= 5
        self.assertEqual(len(parm.atoms), len(mult.atoms))
        for a1, a2 in zip(mult.atoms, parm.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.solvent_radius, a2.solvent_radius)
        self.assertEqual(len(parm.residues), len(mult.residues))
        for r1, r2 in zip(mult.residues, parm.residues):
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.chain, r2.chain)

    def test_simple_slice(self):
        """ Tests simple slicing of AmberParm """
        parm1 = readparm.AmberParm(get_fn('trx.prmtop'))
        parm2 = readparm.AmberParm(get_fn('trx.prmtop'))
        parm2.strip('!@CA,C,O,N,HA,H')
        selection = parm1['@CA,C,O,N,HA,H']
        self.assertIs(type(parm1), readparm.AmberParm)
        self.assertIs(type(parm2), readparm.AmberParm)
        self.assertIs(type(selection), readparm.AmberParm)
        self.assertEqual(len(parm2.atoms), len(selection.atoms))
        self.assertEqual(len(parm2.residues), len(selection.residues))
        self.assertLess(len(parm2.atoms), len(parm1.atoms))
        def cmp_atoms(a1, a2):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.tree, a2.tree)
            self.assertEqual(a1.solvent_radius, a2.solvent_radius)
            self.assertEqual(a1.screen, a2.screen)
            self.assertEqual(a1.join, a2.join)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.idx, a2.residue.idx)
            self.assertEqual(a1.nb_idx, a2.nb_idx)
        for a1, a2 in zip(parm2.atoms, selection.atoms):
            cmp_atoms(a1, a2)
        # Now check valence terms
        self.assertEqual(len(parm2.bonds), len(selection.bonds))
        self.assertEqual(len(parm2.angles), len(selection.angles))
        self.assertEqual(len(parm2.dihedrals), len(selection.dihedrals))
        self.assertGreater(len(parm2.bonds), 0)
        self.assertGreater(len(parm2.angles), 0)
        self.assertGreater(len(parm2.dihedrals), 0)
        for b1, b2 in zip(parm2.bonds, selection.bonds):
            cmp_atoms(b1.atom1, b2.atom1)
            cmp_atoms(b1.atom2, b2.atom2)
            self.assertEqual(b1.type, b2.type)
        for a1, a2 in zip(parm2.angles, selection.angles):
            cmp_atoms(a1.atom1, a2.atom1)
            cmp_atoms(a1.atom2, a2.atom2)
            cmp_atoms(a1.atom3, a2.atom3)
            self.assertEqual(a1.type, a2.type)
        for d1, d2 in zip(parm2.dihedrals, selection.dihedrals):
            cmp_atoms(d1.atom1, d2.atom1)
            cmp_atoms(d1.atom2, d2.atom2)
            cmp_atoms(d1.atom3, d2.atom3)
            cmp_atoms(d1.atom4, d2.atom4)
            self.assertEqual(d1.ignore_end, d2.ignore_end)
            self.assertEqual(d1.improper, d2.improper)
            self.assertEqual(d1.type, d2.type)

class TestAmberMdin(FileIOTestCase):
    """ Tests the Mdin class.... not a good class """

    def test_mdin_API(self):
        """ Tests the Mdin object basic features """
        fn = self.get_fn('test.mdin', written=True)
        mdin1 = mdin.Mdin('sander')
        self.assertEqual(set(mdin1.valid_namelists), {'cntrl', 'ewald', 'qmmm', 'pb'})
        self.assertEqual(mdin1.title, 'mdin prepared by mdin.py')
        self.assertEqual(mdin1.verbosity, 0)
        # What the heck was this for?
        self.assertTrue(mdin1.check())
        mdin1.time()
        self.assertEqual(mdin1.cntrl_nml['dt'], 0.001)
        self.assertEqual(mdin1.cntrl_nml['nstlim'], 1000000)
        self.assertEqual(mdin1.cntrl_nml['imin'], 0)
        mdin1.SHAKE()
        self.assertEqual(mdin1.cntrl_nml['ntf'], 2)
        self.assertEqual(mdin1.cntrl_nml['ntc'], 2)
        self.assertEqual(mdin1.cntrl_nml['dt'], 0.002)
        mdin1.time()
        self.assertEqual(mdin1.cntrl_nml['dt'], 0.002)
        self.assertEqual(mdin1.cntrl_nml['nstlim'], 500000)
        mdin1.constPressure(press=100.0, taup=10.0)
        self.assertEqual(mdin1.cntrl_nml['ntb'], 2)
        self.assertEqual(mdin1.cntrl_nml['ntp'], 1)
        self.assertEqual(mdin1.cntrl_nml['pres0'], 100.0)
        self.assertEqual(mdin1.cntrl_nml['taup'], 10.0)
        mdin1.constVolume()
        self.assertEqual(mdin1.cntrl_nml['ntb'], 1)
        self.assertEqual(mdin1.cntrl_nml['ntp'], 0)
        mdin1.constTemp(ntt=3, temp=100)
        self.assertEqual(mdin1.cntrl_nml['ntt'], 3)
        self.assertEqual(mdin1.cntrl_nml['gamma_ln'], 2.0)
        self.assertEqual(mdin1.cntrl_nml['ig'], -1)
        self.assertEqual(mdin1.cntrl_nml['tautp'], 1.0)
        mdin1.constTemp(ntt=2, temp=100)
        self.assertEqual(mdin1.cntrl_nml['ntt'], 2)
        self.assertEqual(mdin1.cntrl_nml['gamma_ln'], 0)
        self.assertEqual(mdin1.cntrl_nml['ig'], -1)
        self.assertEqual(mdin1.cntrl_nml['tautp'], 1.0)
        mdin1.constpH(solvph=1.0)
        self.assertEqual(mdin1.cntrl_nml['solvph'], 1.0)
        self.assertEqual(mdin1.cntrl_nml['icnstph'], 1)
        self.assertEqual(mdin1.cntrl_nml['ntcnstph'], 10)
        self.assertEqual(mdin1.cntrl_nml['igb'], 2)
        self.assertEqual(mdin1.cntrl_nml['ntb'], 0)
        self.assertEqual(mdin1.cntrl_nml['saltcon'], 0.1)
        mdin1.restrainHeavyAtoms(10.0)
        self.assertEqual(mdin1.cntrl_nml['restraint_wt'], 10)
        self.assertEqual(mdin1.cntrl_nml['restraintmask'], '!@H=')
        mdin1.restrainBackbone(5.0)
        self.assertEqual(mdin1.cntrl_nml['restraint_wt'], 5)
        self.assertEqual(mdin1.cntrl_nml['restraintmask'], '@N,CA,C')
        mdin1.genBorn(igb=8, rgbmax=15.0)
        self.assertEqual(mdin1.cntrl_nml['igb'], 8)
        self.assertEqual(mdin1.cntrl_nml['rgbmax'], 15)
        self.assertEqual(mdin1.cntrl_nml['ntp'], 0)
        self.assertEqual(mdin1.cntrl_nml['ntb'], 0)
        mdin1.time(dt=0.004, time=2000)
        self.assertEqual(mdin1.cntrl_nml['dt'], 0.004)
        self.assertEqual(mdin1.cntrl_nml['nstlim'], 500000)
        self.assertEqual(mdin1.cntrl_nml['imin'], 0)
        mdin1.heat()
        self.assertEqual(mdin1.cntrl_nml['tempi'], 0)
        self.assertEqual(mdin1.cntrl_nml['temp0'], 300)
        self.assertEqual(mdin1.cntrl_nml['ntt'], 3)
        self.assertEqual(mdin1.cntrl_nml['ig'], -1)
        self.assertEqual(mdin1.cntrl_nml['gamma_ln'], 5.0)
        mdin1.restart()
        self.assertEqual(mdin1.cntrl_nml['irest'], 1)
        self.assertEqual(mdin1.cntrl_nml['ntx'], 5)
        mdin1.TI(clambda=0.5)
        self.assertEqual(mdin1.cntrl_nml['clambda'], 0.5)
        self.assertEqual(mdin1.cntrl_nml['icfe'], 1)
        mdin1.softcore_TI(scmask='@1-10', crgmask='@1-10', logdvdl=1000)
        self.assertEqual(mdin1.cntrl_nml['icfe'], 1)
        self.assertEqual(mdin1.cntrl_nml['ifsc'], 1)
        self.assertEqual(mdin1.cntrl_nml['scalpha'], 0.5)
        self.assertEqual(mdin1.cntrl_nml['scmask'], '@1-10')
        self.assertEqual(mdin1.cntrl_nml['crgmask'], '@1-10')
        self.assertEqual(mdin1.cntrl_nml['logdvdl'], 1000)
        mdin1.minimization(imin=5)
        self.assertEqual(mdin1.cntrl_nml['imin'], 5)
        self.assertEqual(mdin1.cntrl_nml['maxcyc'], 1)
        self.assertEqual(mdin1.cntrl_nml['ncyc'], 10)
        self.assertEqual(mdin1.cntrl_nml['ntmin'], 1)
        mdin1.change('cntrl', 'ifqnt', 1)
        mdin1.change('pb', 'istrng', 1.0)
        mdin1.change('qmmm', 'qmcharge', -1)
        mdin1.change('qmmm', 'printdipole', 1)
        mdin1.change('ewald', 'vdwmeth', 0)
        mdin1.change('ewald', 'nfft1', 50)
        mdin1.change('ewald', 'nfft2', 64)
        mdin1.change('ewald', 'nfft3', 96)
        mdin1.AddCard(title='Restraints!', cardString='10.5\nRES 1 10')
        mdin1.write(fn)
        mdin2 = mdin.Mdin(program='sander')
        mdin2.read(fn)
        for var in mdin1.cntrl_nml.keys():
            self.assertEqual(mdin1.cntrl_nml[var], mdin2.cntrl_nml[var])
        for var in mdin1.qmmm_nml.keys():
            self.assertEqual(mdin1.qmmm_nml[var], mdin2.qmmm_nml[var])
        mdin3 = mdin.Mdin(program='pmemd')
        self.assertRaises(InputError, lambda: mdin3.change('cntrl', 'ievb', 1))
        self.assertRaises(InputError, lambda: mdin.Mdin(program='charmm'))
        mdin4 = mdin.Mdin(program='sander.APBS')
        mdin4.change('cntrl', 'igb', 10)
        mdin4.change('pb', 'bcfl', 1)

class TestRst7Class(FileIOTestCase):
    """ Test the Rst7 class """

    def test_ascii(self):
        """ Test the Rst7 class reading ASCII coordinates """
        rst = readparm.Rst7.open(self.get_fn('ash.rst7'))
        np.testing.assert_equal(rst.coordinates, load_file(self.get_fn('ash.rst7')).coordinates)
        np.testing.assert_equal(readparm.Rst7(self.get_fn('ash.rst7')).coordinates, load_file(self.get_fn('ash.rst7')).coordinates)
        rst2 = readparm.Rst7.copy_from(rst)
        np.testing.assert_equal(rst.coordinates, rst2.coordinates)
        rst3 = copy(rst)
        np.testing.assert_equal(rst.coordinates, rst3.coordinates)
        self.assertIsNot(rst, rst3)
        self.assertIsNot(rst.coordinates, rst3.coordinates)

    def test_netcdf(self):
        """ Test the Rst7 class reading NetCDF coordinates """
        rst = readparm.Rst7.open(self.get_fn('ncinpcrd.rst7'))
        np.testing.assert_equal(rst.coordinates, load_file(self.get_fn('ncinpcrd.rst7')).coordinates)
        with self.assertRaises(AmberError):
            readparm.Rst7.open(self.get_fn('trx.prmtop'))
        with self.assertRaises(RuntimeError):
            readparm.Rst7().write(self.get_fn('test.nc', written=True), netcdf=True)

class TestNetCDFTrajectorywithBox(FileIOTestCase):
    """ Test trajecotry with more than 1 frame and with box """

    @unittest.skipIf(PYPY, 'Test does not yet run under pypy')
    def test_netcdf_long_trajectory(self):
        """ Test netcdf trajectory with box """
        parmfile, ncfile = self.get_fn('tz2.parm7'), self.get_fn('tz2.nc')
        parm = pmd.load_file(parmfile, xyz=ncfile, box=np.random.rand(101, 6))
        boxes = parm.get_box('all')
        self.assertEqual(boxes.shape, (101, 6))

class TestAmberTitratableResidues(FileIOTestCase):
    """ Test Amber's titration module capabilities """

    def test_line_buffer(self):
        """ Tests private _LineBuffer for cpin utilities """
        fobj = StringIO()
        fobj2 = StringIO()
        fobj3 = StringIO()
        buf = titratable_residues._LineBuffer(fobj)
        buf2 = titratable_residues._LineBuffer(fobj2)
        buf3 = titratable_residues._LineBuffer(fobj3)
        words = []
        alphabet = list(letters)
        for i in range(random.randint(20, 40)):
            word = ''.join(np.random.choice(alphabet, size=random.randint(5, 20),
                                            replace=True).tolist()
            )
            words.append(word)
            buf.add_word(word)
        buf2.add_words(words)
        buf3.add_words(words, space_delimited=True)
        buf.flush()
        buf.flush() # Make sure second ones don't do anything
        buf2.flush()
        buf3.flush()

        fobj.seek(0)
        fobj2.seek(0)
        fobj3.seek(0)
        self.assertEqual(fobj.read(), fobj2.read())
        self.assertNotEqual(fobj.read(), fobj3.read())
        fobj3.seek(0)
        self.assertEqual(fobj3.read().split(), words)

    def test_old_cpin_creation(self):
        """ Test TitratableResidueList and cpin creation at the old format """
        import cpinutil
        parm = self.get_fn('trx.prmtop')
        output = self.get_fn('test.old.cpin', written=True)
        opt = cpinutil.parser.parse_args(
            ['-igb', '2', '-p', parm,'--old-format', '-states', '0,0,1,0,1,1,0,1,0,1,1,1', '-o', output]
        )
        cpinutil.main(opt)
        self.assertTrue(
            diff_files(self.get_fn('test.old.cpin', saved=True), self.get_fn('test.old.cpin', written=True),
                       absolute_error=1e-6, spacechar='=,')
        )

    def test_cpin_creation(self):
        """ Test TitratableResidueList and cpin creation """
        import cpinutil
        parm = self.get_fn('trx.prmtop')
        output = self.get_fn('test.cpin', written=True)
        opt = cpinutil.parser.parse_args(
            ['-igb', '2', '-p', parm, '-states', '0,0,1,0,1,1,0,1,0,1,1,1', '-o', output]
        )
        cpinutil.main(opt)
        self.assertTrue(
            diff_files(self.get_fn('test.cpin', saved=True), self.get_fn('test.cpin', written=True),
                       absolute_error=1e-6, spacechar='=,')
        )

    def test_cein_creation(self):
        """ Test cein creation """
        import ceinutil
        parm = self.get_fn('mp8.prmtop')
        output = self.get_fn('mp8.cein', written=True)
        opt = ceinutil.parser.parse_args(['-igb', '2', '-p', parm, '-o', output])
        ceinutil.main(opt)
        self.assertTrue(
            diff_files(self.get_fn('mp8.cein', saved=True), self.get_fn('mp8.cein', written=True),
                       absolute_error=1e-6, spacechar='=,')
        )

    def test_cpein_creation(self):
        """ Test cpein creation """
        import cpeinutil
        parm = self.get_fn('tyx.prmtop')
        output = self.get_fn('tyx.cpein', written=True)
        opt = cpeinutil.parser.parse_args(['-igb', '2', '-p', parm, '-o', output])
        cpeinutil.main(opt)
        self.assertTrue(
            diff_files(self.get_fn('tyx.cpein', saved=True), self.get_fn('tyx.cpein', written=True),
                       absolute_error=1e-6, spacechar='=,')
        )
        parm = self.get_fn('mp8.prmtop')
        output = self.get_fn('mp8.cpein', written=True)
        opt = cpeinutil.parser.parse_args(
            ['-igb', '2', '-p', parm, '-o', output]
        )
        cpeinutil.main(opt)
        self.assertTrue(
            diff_files(self.get_fn('mp8.cpein', saved=True), self.get_fn('mp8.cpein', written=True),
                       absolute_error=1e-6, spacechar='=,')
        )

    def test_titratable_residue(self):
        """ Tests the TitratableResidue object """
        as4 = titratable_residues.AS4
        self.assertEqual(str(as4), saved.AS4_TITR_OUTPUT)
        # Test error handling for TitratableResidue
        newres = titratable_residues.TitratableResidue(
                'NWR', ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'], "ph", 7.0,
        )
        self.assertEqual(newres.pKa, 7.0)
        self.assertRaises(AmberError, lambda:
                newres.add_state([10.0, 20.0], 10.0, 10.0, 3, 7.0)
        )
        self.assertRaises(AmberError, lambda:
            newres.add_states([[1, 2, 3, 4, 5, 6, 7], [2, 3, 4,5, 6, 7]],
                              [10, 20, 30], [10, 20, 30], [3, 2, 1], [7.0, 0.0, 0.0])
        )
        self.assertRaises(AmberError, lambda: newres.cpin_pointers(10))
        newres.set_first_state(0)
        newres.set_first_state(0) # Second setting should be ignored
        self.assertRaises(AmberError, lambda: newres.set_first_state(1))
        newres.set_first_charge(0)
        newres.set_first_charge(0)
        self.assertRaises(AmberError, lambda: newres.set_first_charge(1))
