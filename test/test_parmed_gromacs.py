"""
Tests the functionality in the parmed.gromacs package
"""
from contextlib import closing
import copy
import sys
import os
import unittest
import warnings

import numpy as np

from parmed import (load_file, Structure, ExtraPoint, DihedralTypeList, Atom,
                    ParameterSet, Bond, NonbondedException, DihedralType,
                    RBTorsionType, Improper, Cmap, UreyBradley, BondType,
                    UnassignedAtomType, NonbondedExceptionType, NoUreyBradley)
from parmed.charmm import CharmmParameterSet
from parmed.exceptions import GromacsWarning, GromacsError, ParameterError
from parmed.gromacs import GromacsTopologyFile, GromacsGroFile
from parmed.gromacs._gromacsfile import GromacsFile
from parmed import gromacs as gmx, periodic_table
from parmed.topologyobjects import UnassignedAtomType
from parmed.utils.six.moves import range, zip, StringIO
from utils import (get_fn, diff_files, get_saved_fn, FileIOTestCase, HAS_GROMACS,
                   create_random_structure)
import utils

@unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without Gromacs")
class TestGromacsTop(FileIOTestCase):
    """ Tests the Gromacs topology file parser """

    def _charmm27_checks(self, top):
        # Check that the number of terms are correct
        self.assertEqual(len(top.atoms), 1960)
        self.assertEqual(len(top.bonds), 1984)
        self.assertEqual(len(top.angles), 3547)
        self.assertEqual(len(top.dihedrals), 5187)
        self.assertEqual(len(top.impropers), 351)
        self.assertEqual(len(top.cmaps), 127)
        self.assertEqual(len(top.adjusts), 5106)
        self.assertFalse(top.unknown_functional)
        # Check the first and last of most of the terms to make sure that they
        # are the same as what is defined in the topology file
        self.assertEqual(top.atoms[0].type, 'NH3')
        self.assertEqual(top.atoms[0].name, 'N')
        self.assertEqual(top.atoms[0].mass, 14.007)
        self.assertEqual(top.atoms[0].charge, -0.3)
        self.assertEqual(top.atoms[0].atomic_number, 7)
        self.assertEqual(top.atoms[0].residue.name, 'LYS')
        self.assertEqual(top.atoms[1].type, 'HC')
        self.assertEqual(top.atoms[1].name, 'H1')
        self.assertEqual(top.atoms[1].mass, 1.008)
        self.assertEqual(top.atoms[1].charge, 0.33)
        self.assertEqual(top.atoms[1].atomic_number, 1)
        self.assertEqual(top.atoms[1].residue.name, 'LYS')
        self.assertEqual(top.atoms[1958].type, 'OC')
        self.assertEqual(top.atoms[1958].name, 'OT1')
        self.assertEqual(top.atoms[1958].mass, 15.9994)
        self.assertEqual(top.atoms[1958].charge, -0.67)
        self.assertEqual(top.atoms[1958].atomic_number, 8)
        self.assertEqual(top.atoms[1958].residue.name, 'LEU')
        self.assertEqual(top.atoms[1959].type, 'OC')
        self.assertEqual(top.atoms[1959].name, 'OT2')
        self.assertEqual(top.atoms[1959].mass, 15.9994)
        self.assertEqual(top.atoms[1959].charge, -0.67)
        self.assertEqual(top.atoms[1959].atomic_number, 8)
        self.assertEqual(top.atoms[1959].residue.name, 'LEU')
        # Bonds
        self.assertIs(top.bonds[0].atom1, top.atoms[0])
        self.assertIs(top.bonds[0].atom2, top.atoms[1])
        self.assertEqual(top.bonds[0].funct, 1)
        self.assertIs(top.bonds[1983].atom1, top.atoms[1957])
        self.assertIs(top.bonds[1983].atom2, top.atoms[1959])
        self.assertEqual(top.bonds[1983].funct, 1)
        # Angles
        self.assertIs(top.angles[0].atom1, top.atoms[1])
        self.assertIs(top.angles[0].atom2, top.atoms[0])
        self.assertIs(top.angles[0].atom3, top.atoms[2])
        self.assertEqual(top.angles[0].funct, 5)
        self.assertIs(top.angles[3546].atom1, top.atoms[1958])
        self.assertIs(top.angles[3546].atom2, top.atoms[1957])
        self.assertIs(top.angles[3546].atom3, top.atoms[1959])
        self.assertEqual(top.angles[0].funct, 5)
        # Dihedrals
        self.assertIs(top.dihedrals[0].atom1, top.atoms[1])
        self.assertIs(top.dihedrals[0].atom2, top.atoms[0])
        self.assertIs(top.dihedrals[0].atom3, top.atoms[4])
        self.assertIs(top.dihedrals[0].atom4, top.atoms[5])
        self.assertEqual(top.dihedrals[0].funct, 9)
        self.assertIs(top.dihedrals[5186].atom1, top.atoms[1949])
        self.assertIs(top.dihedrals[5186].atom2, top.atoms[1947])
        self.assertIs(top.dihedrals[5186].atom3, top.atoms[1953])
        self.assertIs(top.dihedrals[5186].atom4, top.atoms[1956])
        self.assertEqual(top.dihedrals[5186].funct, 9)
        # Impropers
        self.assertIs(top.impropers[0].atom1, top.atoms[22])
        self.assertIs(top.impropers[0].atom2, top.atoms[4])
        self.assertIs(top.impropers[0].atom3, top.atoms[24])
        self.assertIs(top.impropers[0].atom4, top.atoms[23])
        self.assertEqual(top.impropers[0].funct, 2)
        self.assertIs(top.impropers[350].atom1, top.atoms[1957])
        self.assertIs(top.impropers[350].atom2, top.atoms[1942])
        self.assertIs(top.impropers[350].atom3, top.atoms[1959])
        self.assertIs(top.impropers[350].atom4, top.atoms[1958])
        self.assertEqual(top.impropers[350].funct, 2)
        # Cmaps
        self.assertIs(top.cmaps[0].atom1, top.atoms[22])
        self.assertIs(top.cmaps[0].atom2, top.atoms[24])
        self.assertIs(top.cmaps[0].atom3, top.atoms[26])
        self.assertIs(top.cmaps[0].atom4, top.atoms[38])
        self.assertIs(top.cmaps[0].atom5, top.atoms[40])
        self.assertEqual(top.cmaps[0].funct, 1)
        self.assertIs(top.cmaps[126].atom1, top.atoms[1914])
        self.assertIs(top.cmaps[126].atom2, top.atoms[1916])
        self.assertIs(top.cmaps[126].atom3, top.atoms[1918])
        self.assertIs(top.cmaps[126].atom4, top.atoms[1938])
        self.assertIs(top.cmaps[126].atom5, top.atoms[1940])
        self.assertEqual(top.cmaps[126].funct, 1)
        # Adjusts
        self.assertIs(top.adjusts[0].atom1, top.atoms[0])
        self.assertIs(top.adjusts[0].atom2, top.atoms[7])
        self.assertEqual(top.adjusts[0].funct, 1)
        self.assertIs(top.adjusts[5105].atom1, top.atoms[1952])
        self.assertIs(top.adjusts[5105].atom2, top.atoms[1953])
        self.assertEqual(top.adjusts[5105].funct, 1)

    def test_charmm27_top(self):
        """ Tests parsing a Gromacs topology with CHARMM 27 FF """
        top = GromacsTopologyFile(get_fn('1aki.charmm27.top'))
        self.assertEqual(top.combining_rule, 'lorentz')
        self.assertEqual(top.itps, ['charmm27.ff/forcefield.itp',
                                    'charmm27.ff/tip3p.itp',
                                    'charmm27.ff/ions.itp'])
        self._charmm27_checks(top)

    def test_gromacs_top_detection(self):
        """ Tests automatic file detection of GROMACS topology files """
        fn = self.get_fn('test.top', written=True)
        with open(fn, 'w') as f:
            f.write('# not a gromacs topology file\n')
        self.assertFalse(GromacsTopologyFile.id_format(fn))
        with open(fn, 'w') as f:
            pass
        self.assertFalse(GromacsTopologyFile.id_format(fn))

    def test_gromacs_file(self):
        """ Test GromacsFile helper class """
        f = StringIO('This is the first line \\\n  this is still the first line'
                     '\nThis is the second line')
        gf = GromacsFile(f)
        self.assertEqual(gf.readline(), 'This is the first line    this is '
                         'still the first line\n')
        self.assertEqual(gf.readline(), 'This is the second line')
        self.assertEqual(gf.readline(), '')
        f.seek(0)
        lines = [line for line in gf]
        self.assertEqual(lines[0], 'This is the first line    this is '
                         'still the first line\n')
        self.assertEqual(lines[1], 'This is the second line')

        # Try with comments now
        f = StringIO('This is the first line \\\n  this is still the first line'
                     ' ; and this is a comment\nThis is the second line ; and '
                     'this is also a comment')
        gf = GromacsFile(f)
        self.assertEqual(gf.readline(), 'This is the first line    this is still'
                         ' the first line \n')
        self.assertEqual(gf.readline(), 'This is the second line \n')
        f.seek(0)
        lines = [line for line in gf]
        self.assertEqual(lines[0], 'This is the first line    this is still'
                         ' the first line \n')
        self.assertEqual(lines[1], 'This is the second line \n')
        f.seek(0)
        lines = gf.readlines()
        self.assertEqual(lines[0], 'This is the first line    this is still'
                         ' the first line \n')
        self.assertEqual(lines[1], 'This is the second line \n')
        f.seek(0)
        self.assertEqual(gf.read(), 'This is the first line    this is still'
                         ' the first line \nThis is the second line \n')

        # Error handling
        self.assertRaises(IOError, lambda:
                GromacsFile('some_file that does_not exist'))

    def test_write_charmm27_top(self):
        """ Tests writing a Gromacs topology file with CHARMM 27 FF """
        top = load_file(get_fn('1aki.charmm27.top'))
        self.assertEqual(top.combining_rule, 'lorentz')
        GromacsTopologyFile.write(top, self.get_fn('1aki.charmm27.top', written=True))
        top2 = load_file(self.get_fn('1aki.charmm27.top', written=True))
        self._charmm27_checks(top)

    def test_write_with_unicode_escape(self):
        """ Tests writing a .top when sys.argv includes \n in a string"""
        top = load_file(get_fn('1aki.charmm27.top'))
        self.assertEqual(top.combining_rule, 'lorentz')
        sys.argv.append('foo\nbar')
        GromacsTopologyFile.write(top, get_fn('foobar.top'))
        with open(get_fn('foobar.top')) as fh:
            for _ in range(15):
                assert fh.readline().startswith(';')
        sys.argv.remove('foo\nbar')

    def test_moleculetype_distinction(self):
        """ Tests moleculetype distinction for different parameters """
        parm = load_file(get_fn('different_molecules.parm7'))
        parm.save(self.get_fn('different_molecules.top', written=True),
                  format='gromacs', overwrite=True)
        top = load_file(self.get_fn('different_molecules.top', written=True))
        # Make sure all atoms have the same parameters
        for a1, a2 in zip(parm.atoms, top.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertAlmostEqual(a1.charge, a2.charge, places=4)

    def _check_ff99sbildn(self, top):
        self.assertEqual(len(top.atoms), 4235)
        self.assertEqual(len(top.residues), 1046)
        self.assertEqual(sum(1 for a in top.atoms if isinstance(a, ExtraPoint)),
                         1042)
        self.assertEqual(len(top.bonds), 3192)
        self.assertEqual(len(top.angles), 1162)
        self.assertEqual(len(top.dihedrals), 179)

    def _check_equal_structures(self, top1, top2):
        def cmp_atoms(a1, a2):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.atom_type, a2.atom_type)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.idx, a2.residue.idx)

        def cmp_valence(val1, val2, typeattrs=None):
            self.assertEqual(len(val1), len(val2))
            for v1, v2 in zip(val1, val2):
                self.assertIs(type(v1), type(v2))
                attrs = [attr for attr in dir(v1) if attr.startswith('atom')]
                atoms1 = [getattr(v1, attr) for attr in attrs]
                atoms2 = [getattr(v2, attr) for attr in attrs]
                for a1, a2 in zip(atoms1, atoms2):
                    cmp_atoms(a1, a2)
                # Check the type lists
                if typeattrs is not None:
                    for attr in typeattrs:
                        self.assertAlmostEqual(getattr(v1.type, attr),
                                               getattr(v2.type, attr), places=5)
                else:
                    self.assertEqual(v1.type, v2.type)

        def cmp_dihedrals(dih1, dih2):
            self.assertEqual(len(dih1), len(dih2))
            for v1, v2 in zip(dih1, dih2):
                self.assertIs(type(v1), type(v2))
                self.assertIs(type(v1.type), type(v2.type))
                atoms1 = [v1.atom1, v1.atom2, v1.atom3, v1.atom4]
                atoms2 = [v2.atom1, v2.atom2, v2.atom3, v2.atom4]
                for a1, a2 in zip(atoms1, atoms2):
                    cmp_atoms(a1, a2)
                self.assertEqual(v1.improper, v2.improper)
                self.assertEqual(v1.ignore_end, v2.ignore_end)
                if isinstance(v1, DihedralTypeList):
                    self.assertEqual(len(v1.type), len(v2.type))
                    for vt1, vt2 in zip(v1.type, v2.type):
                        self.assertAlmostEqual(v1.type.phi_k, v2.type.phi_k, places=5)
                        self.assertAlmostEqual(v1.type.per, v2.type.per, places=5)
                        self.assertAlmostEqual(v1.type.phase, v2.type.phase, places=5)

        self.assertEqual(len(top1.atoms), len(top2.atoms))
        for a1, a2 in zip(top1.atoms, top2.atoms):
            cmp_atoms(a1, a2)
        cmp_valence(top1.bonds, top2.bonds, ['k', 'req'])
        cmp_valence(top1.angles, top2.angles, ['k', 'theteq'])
        cmp_dihedrals(top1.dihedrals, top2.dihedrals)

    def test_read_amber99SBILDN(self):
        """ Tests parsing a Gromacs topology with Amber99SBILDN and water """
        top = load_file(get_fn('ildn.solv.top'))
        self.assertEqual(top.combining_rule, 'lorentz')
        self._check_ff99sbildn(top)
        dts = top.dihedral_types[:]
        top.join_dihedrals()
        for dt1, dt2 in zip(dts, top.dihedral_types):
            self.assertIs(dt1, dt2)

    def test_write_amber99SBILDN(self):
        """ Tests writing a Gromacs topology with multiple molecules """
        top = load_file(get_fn('ildn.solv.top'))
        self.assertEqual(top.combining_rule, 'lorentz')
        fn = self.get_fn('ildn.solv.top', written=True)
        top.write(fn, combine=None)
        top2 = load_file(fn)
        self._check_ff99sbildn(top2)
        self._check_equal_structures(top, top2)
        # Turn off gen_pairs and write another version
        top.defaults.gen_pairs = 'no'
        top.write(fn)
        top3 = load_file(fn)
        self._check_ff99sbildn(top3)
        self._check_equal_structures(top, top3)

    def test_duplicate_system_names(self):
        """ Tests that Gromacs topologies never have duplicate moleculetypes """
        parm = load_file(get_fn('phenol.prmtop'))
        parm = parm * 20 + load_file(get_fn('biphenyl.prmtop')) * 20
        top = GromacsTopologyFile.from_structure(parm)
        self.assertEqual(top.combining_rule, 'lorentz')
        top.write(self.get_fn('phenol_biphenyl.top', written=True))
        top2 = GromacsTopologyFile(self.get_fn('phenol_biphenyl.top', written=True))
        self.assertEqual(len(top.residues), 40)

        # Now test this when we use "combine"
        parm = load_file(os.path.join(get_fn('12.DPPC'), 'topol3.top'))
        fn = self.get_fn('samename.top', written=True)
        parm.residues[3].name = 'SOL' # Rename a DPPC to SOL
        parm.write(fn, combine=[[0, 1]])
        parm2 = load_file(fn)
        self.assertEqual(len(parm2.atoms), len(parm.atoms))
        self.assertEqual(len(parm2.residues), len(parm2.residues))
        for a1, a2 in zip(parm.atoms, parm2.atoms):
            self._equal_atoms(a1, a2)
        for r1, r2 in zip(parm.residues, parm2.residues):
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.name, r2.name)
            for a1, a2 in zip(r1, r2):
                self._equal_atoms(a1, a2)

    def test_gromacs_top_from_structure(self):
        """ Tests the GromacsTopologyFile.from_structure constructor """
        struct = create_random_structure(True)
        with self.assertRaises(TypeError):
            GromacsTopologyFile.from_structure(struct)
        parm = load_file(get_fn('ash.parm7'))
        parm.dihedrals[0].type.scee = 8.0
        with self.assertRaises(GromacsError):
            GromacsTopologyFile.from_structure(parm)
        for dt in parm.dihedral_types: dt.scee = dt.scnb = 0
        top = GromacsTopologyFile.from_structure(parm)
        self.assertEqual(top.defaults.fudgeLJ, 1.0)
        self.assertEqual(top.defaults.fudgeQQ, 1.0)

    def test_OPLS(self):
        """ Tests the geometric combining rules in Gromacs with OPLS/AA """
        parm = load_file(os.path.join(get_fn('05.OPLS'), 'topol.top'),
                         xyz=os.path.join(get_fn('05.OPLS'), 'conf.gro'))
        self.assertEqual(parm.combining_rule, 'geometric')
        self.assertEqual(parm.defaults.comb_rule, 3)
        parm.write(self.get_fn('test.topol', written=True), combine='all')
        parm2 = load_file(self.get_fn('test.topol', written=True))
        self.assertEqual(len(parm.atoms), len(parm2.atoms))
        # Check that the charge attribute is read correctly
        self.assertEqual(parm.parameterset.atom_types['opls_001'].charge, 0.5)
        # Check the xyz argument in the constructor being coordinates
        parm2 = load_file(os.path.join(get_fn('05.OPLS'), 'topol.top'),
                          xyz=parm.coordinates, box=[10, 10, 10, 90, 90, 90])
        np.testing.assert_equal(parm2.coordinates, parm.coordinates)
        np.testing.assert_equal(parm2.box, [10, 10, 10, 90, 90, 90])
        # Check the copy constructor
        p2 = GromacsTopologyFile.from_structure(parm, copy=True)
        self.assertEqual(p2.combining_rule, 'geometric')
        self.assertEqual(p2.defaults.comb_rule, 3)
        self.assertEqual(len(p2.atoms), len(parm.atoms))
        for a1, a2 in zip(p2.atoms, parm.atoms):
            self.assertIsNot(a1, a2)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.mass, a2.mass)
        np.testing.assert_equal(p2.box, parm.box)

    def test_write_settles(self):
        """ Tests that settles is only written for water """
        fn = self.get_fn('blah.top', written=True)
        parm = load_file(os.path.join(get_fn('01.1water'), 'topol.top'))
        parm[0].atomic_number = parm[0].atom_type.atomic_number = 7
        parm.write(fn)
        with closing(GromacsFile(fn)) as f:
            for line in f:
                self.assertNotIn('settles', line)

    def test_write_extra_points(self):
        """ Test writing of GROMACS files with virtual sites """
        f = StringIO('; TIP4Pew water molecule\n#include "amber99sb.ff/forcefield.itp"\n'
                     '#include "amber99sb.ff/tip4pew.itp"\n[ system ]\nWATER\n'
                     '[ molecules ]\nSOL 1\n')
        parm = GromacsTopologyFile(f)
        fn = self.get_fn('test.top', written=True)
        parm.write(fn)
        parm2 = load_file(fn)
        self.assertEqual(len(parm.atoms), len(parm2.atoms))
        self.assertEqual(len(parm.bonds), len(parm2.bonds))

    def test_without_parametrize(self):
        """ Tests loading a Gromacs topology without parametrizing """
        parm = load_file(os.path.join(get_fn('05.OPLS'), 'topol.top'),
                         xyz=os.path.join(get_fn('05.OPLS'), 'conf.gro'), parametrize=False)
        self.assertIs(parm.atoms[0].atom_type, UnassignedAtomType)
        self.assertTrue(all(x.type is None for x in parm.bonds))
        # Now try writing it out again
        fn = self.get_fn('test.top', written=True)
        parm.write(fn)

        parm = load_file(os.path.join(get_fn('04.Ala'), 'topol.top'), parametrize=False)
        # Add an improper
        parm.impropers.append(Improper(*parm.atoms[:4]))
        parm.write(fn)

    def test_bad_top_loads(self):
        """ Test error catching in GromacsTopologyFile reading """
        fn = os.path.join(get_fn('03.AlaGlu'), 'topol.top')
        self.assertRaises(TypeError, lambda: load_file(fn, xyz=fn))
        self.assertRaises(ValueError, lambda: GromacsTopologyFile(xyz=1, box=1))
        wfn = os.path.join(get_fn('gmxtops'), 'duplicated_topol.top')
        self.assertRaises(GromacsError, lambda: GromacsTopologyFile(wfn))
        f = StringIO('\n[ defaults ]\n; not enough data\n 1\n\n')
        self.assertRaises(GromacsError, lambda: GromacsTopologyFile(f))
        f = StringIO('\n[ defaults ]\n; unsupported function type\n 2 1 yes\n')
        with self.assertWarns(GromacsWarning):
            GromacsTopologyFile(f)
        f.seek(0)
        self.assertTrue(GromacsTopologyFile(f).unknown_functional)
        fn = os.path.join(get_fn('gmxtops'), 'bad_vsites3.top')
        self.assertRaises(GromacsError, lambda: load_file(fn))
        self.assertRaises(ValueError, lambda: load_file(fn, defines=dict()))
        f = StringIO('\n[ defaults ]\n1 1 yes\n\n[ system ]\nname\n[ molecules ]\nNOMOL  2\n')
        self.assertRaises(GromacsError, lambda: GromacsTopologyFile(f))
        fn = os.path.join(get_fn('gmxtops'), 'bad_nrexcl.top')
        with self.assertWarns(GromacsWarning):
            GromacsTopologyFile(fn, defines=dict(SMALL_NREXCL=1))
        with self.assertWarns(GromacsWarning):
            GromacsTopologyFile(fn)
        with self.assertWarns(GromacsWarning), self.assertRaises(GromacsError):
            GromacsTopologyFile(wfn, defines=dict(NODUP=1))
        with self.assertRaises(GromacsError):
            GromacsTopologyFile(wfn, defines=dict(NODUP=1, BADNUM=1))
        self.assertRaises(RuntimeError, GromacsTopologyFile().parametrize)

    def test_top_parsing_missing_types(self):
        """ Test GROMACS topology files with missing types """
        fn = os.path.join(get_fn('gmxtops'), 'missing_atomtype.top')
        with self.assertWarns(GromacsWarning):
            GromacsTopologyFile(fn, parametrize=False)
        top = GromacsTopologyFile(fn, parametrize=False)
        self.assertIs(top[0].atom_type, UnassignedAtomType)
        self.assertEqual(top[0].mass, -1)
        self.assertEqual(top[0].atomic_number, -1)
        self.assertEqual(top[1].atomic_number, 1)  # taken from atom_type
        self.assertEqual(top[-1].atomic_number, 1) # taken from atom_type
        self.assertEqual(top[-1].charge, 0) # removed
        self.assertEqual(top.bonds[0].funct, 2)
        self.assertTrue(top.unknown_functional)

    def test_molecule_ordering(self):
        """ Tests non-contiguous atoms in Gromacs topology file writes """
        parm = load_file(os.path.join(get_fn('12.DPPC'), 'topol3.top'))
        parm.write(self.get_fn('topol3.top', written=True))
        parm2 = load_file(self.get_fn('topol3.top', written=True))
        self.assertEqual(len(parm.atoms), len(parm2.atoms))
        self.assertEqual(len(parm.residues), len(parm2.residues))
        for r1, r2 in zip(parm.residues, parm2.residues):
            self.assertEqual(r1.name, r2.name)
            for a1, a2 in zip(r1.atoms, r2.atoms):
                self.assertEqual(a1.name, a2.name)
                self.assertEqual(a1.type, a2.type)

    def test_copying_defaults(self):
        """ Tests that copying GromacsTopologyFile copies Defaults too """
        parm = load_file(get_fn('159.top'))
        newfile = StringIO()
        copy.copy(parm).write(newfile)
        newfile.seek(0)
        newparm = GromacsTopologyFile(newfile)
        self.assertEqual(parm.defaults, newparm.defaults)

    def test_getitem_defaults(self):
        """ Tests that GromacsTopologyFile[] sets Defaults correctly """
        parm = load_file(get_fn('159.top'))
        newfile = StringIO()
        parm[0,:].write(newfile)
        newfile.seek(0)
        newparm = GromacsTopologyFile(newfile)
        self.assertEqual(parm.defaults, newparm.defaults)

    def test_molecule_combine(self):
        """ Tests selective molecule combination in Gromacs topology files """
        parm = load_file(os.path.join(get_fn('12.DPPC'), 'topol3.top'))
        fname = self.get_fn('combined.top', written=True)
        # Make sure that combining non-adjacent molecules fails
        self.assertRaises(ValueError, lambda:
                parm.write(fname, combine=[[1, 3]]))
        self.assertRaises(ValueError, lambda:
                parm.write(fname, combine='joey'))
        self.assertRaises(TypeError, lambda:
                parm.write(fname, combine=[1, 2, 3]))
        self.assertRaises(TypeError, lambda:
                parm.write(fname, combine=1))
        parm.write(fname, combine=[[3, 4], [126, 127, 128, 129, 130]])
        with open(fname, 'r') as f:
            for line in f:
                if line.startswith('[ molecules ]'):
                    break
            molecule_list = []
            for line in f:
                if line[0] == ';': continue
                words = line.split()
                molecule_list.append((words[0], int(words[1])))
        parm2 = load_file(fname)
        self.assertEqual(molecule_list, [('DPPC', 3), ('system1', 1),
                         ('SOL', 121), ('system2', 1), ('SOL', 121)])
        self.assertEqual(len(parm2.atoms), len(parm.atoms))
        self.assertEqual(len(parm2.residues), len(parm2.residues))
        for a1, a2 in zip(parm.atoms, parm2.atoms):
            self._equal_atoms(a1, a2)
        for r1, r2 in zip(parm.residues, parm2.residues):
            self.assertEqual(len(r1), len(r2))
            for a1, a2 in zip(r1, r2):
                self._equal_atoms(a1, a2)
        # Now combine multiple DPPC and multiple waters in the same molecule
        parm.write(fname, combine=[[2, 3, 4, 5], [128, 129, 130, 131]])
        with open(fname, 'r') as f:
            for line in f:
                if line.startswith('[ molecules ]'):
                    break
            molecule_list = []
            for line in f:
                if line[0] == ';': continue
                words = line.split()
                molecule_list.append((words[0], int(words[1])))
        self.assertEqual(molecule_list, [('DPPC', 2), ('system1', 1),
                         ('SOL', 120), ('DPPC', 2), ('system2', 1), ('SOL', 120)])
        parm2 = load_file(fname)
        self.assertEqual(len(parm2.atoms), len(parm.atoms))
        self.assertEqual(len(parm2.residues), len(parm2.residues))
        for a1, a2 in zip(parm.atoms, parm2.atoms):
            self._equal_atoms(a1, a2)
        for r1, r2 in zip(parm.residues, parm2.residues):
            self.assertEqual(len(r1), len(r2))
            for a1, a2 in zip(r1, r2):
                self._equal_atoms(a1, a2)

    def test_gromacs_top_write(self):
        """ Tests the GromacsTopologyFile writer """
        def total_diheds(dlist):
            n = 0
            for d in dlist:
                if isinstance(d.type, DihedralTypeList):
                    n += len(d.type)
                elif not d.improper:
                    n += 1
            return n
        parm = load_file(get_fn('ash.parm7'))
        top = GromacsTopologyFile.from_structure(parm)
        self.assertRaises(TypeError, lambda: top.write(10))
        f = StringIO()
        self.assertRaises(ValueError, lambda: top.write(f, parameters=10))
        # Write parameters and topology to same filename
        fn = self.get_fn('test.top', written=True)
        top.write(fn, parameters=fn)
        top2 = load_file(fn)
        self.assertEqual(len(top2.atoms), len(top.atoms))
        self.assertEqual(len(top2.bonds), len(top.bonds))
        self.assertEqual(len(top2.angles), len(top.angles))
        self.assertEqual(total_diheds(top2.dihedrals), total_diheds(top.dihedrals))
        for a1, a2 in zip(top2.atoms, top.atoms):
            self.assertAlmostEqual(a1.atom_type.sigma, a2.atom_type.sigma, places=3)
            self.assertAlmostEqual(a1.atom_type.epsilon, a2.atom_type.epsilon, places=3)
            self.assertEqual(a1.atom_type.name, a2.atom_type.name)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(set(a.name for a in a1.bond_partners),
                             set(a.name for a in a2.bond_partners))
        # Now try passing open files
        with open(fn, 'w') as f:
            top.write(f, parameters=f)
        top2 = load_file(fn)
        self.assertEqual(len(top2.atoms), len(top.atoms))
        self.assertEqual(len(top2.bonds), len(top.bonds))
        self.assertEqual(len(top2.angles), len(top.angles))
        self.assertEqual(total_diheds(top2.dihedrals), total_diheds(top.dihedrals))
        for a1, a2 in zip(top2.atoms, top.atoms):
            self.assertAlmostEqual(a1.atom_type.sigma, a2.atom_type.sigma, places=3)
            self.assertAlmostEqual(a1.atom_type.epsilon, a2.atom_type.epsilon, places=3)
            self.assertEqual(a1.atom_type.name, a2.atom_type.name)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(set(a.name for a in a1.bond_partners),
                             set(a.name for a in a2.bond_partners))
        # Now try separate parameter/topology file
        fn2 = self.get_fn('test.itp', written=True)
        top.write(fn, parameters=fn2)
        top2 = load_file(fn)
        self.assertEqual(len(top2.atoms), len(top.atoms))
        self.assertEqual(len(top2.bonds), len(top.bonds))
        self.assertEqual(len(top2.angles), len(top.angles))
        self.assertEqual(total_diheds(top2.dihedrals), total_diheds(top.dihedrals))
        for a1, a2 in zip(top2.atoms, top.atoms):
            self.assertAlmostEqual(a1.atom_type.sigma, a2.atom_type.sigma, places=3)
            self.assertAlmostEqual(a1.atom_type.epsilon, a2.atom_type.epsilon, places=3)
            self.assertEqual(a1.atom_type.name, a2.atom_type.name)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(set(a.name for a in a1.bond_partners),
                             set(a.name for a in a2.bond_partners))
        # Now try separate parameter/topology/molfile files
        fn3 = self.get_fn('test_mol.itp', written=True)
        top.write(fn, parameters=fn2, molfile=fn3)
        top2 = load_file(fn)
        self.assertEqual(len(top2.atoms), len(top.atoms))
        self.assertEqual(len(top2.bonds), len(top.bonds))
        self.assertEqual(len(top2.angles), len(top.angles))
        self.assertEqual(total_diheds(top2.dihedrals), total_diheds(top.dihedrals))
        for a1, a2 in zip(top2.atoms, top.atoms):
            self.assertAlmostEqual(a1.atom_type.sigma, a2.atom_type.sigma, places=3)
            self.assertAlmostEqual(a1.atom_type.epsilon, a2.atom_type.epsilon, places=3)
            self.assertEqual(a1.atom_type.name, a2.atom_type.name)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(set(a.name for a in a1.bond_partners),
                             set(a.name for a in a2.bond_partners))
        # Now force writing pair types...
        top.defaults.gen_pairs = 'no'
        top.write(fn, parameters=fn)
        top2 = load_file(fn)
        self.assertEqual(len(top2.atoms), len(top.atoms))
        self.assertEqual(len(top2.bonds), len(top.bonds))
        self.assertEqual(len(top2.angles), len(top.angles))
        self.assertEqual(total_diheds(top2.dihedrals), total_diheds(top.dihedrals))
        for a1, a2 in zip(top2.atoms, top.atoms):
            self.assertAlmostEqual(a1.atom_type.sigma, a2.atom_type.sigma, places=3)
            self.assertAlmostEqual(a1.atom_type.epsilon, a2.atom_type.epsilon, places=3)
            self.assertEqual(a1.atom_type.name, a2.atom_type.name)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(set(a.name for a in a1.bond_partners),
                             set(a.name for a in a2.bond_partners))

        # Now force writing pair types to [ pairtypes ] (instead of in-line)
        fn2 = self.get_fn('testpairtypes.top', written=True)
        top2.write(fn2, parameters=fn2)
        top3 = load_file(fn2)
        self.assertEqual(top3.defaults.gen_pairs, 'no')
        self.assertEqual(len(top2.atoms), len(top3.atoms))
        self.assertEqual(len(top2.bonds), len(top3.bonds))
        self.assertEqual(len(top2.angles), len(top3.angles))
        self.assertEqual(total_diheds(top2.dihedrals), total_diheds(top3.dihedrals))
        for a1, a2 in zip(top2.atoms, top3.atoms):
            self.assertAlmostEqual(a1.atom_type.sigma, a2.atom_type.sigma, places=3)
            self.assertAlmostEqual(a1.atom_type.epsilon, a2.atom_type.epsilon, places=3)
            self.assertEqual(a1.atom_type.name, a2.atom_type.name)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(set(a.name for a in a1.bond_partners),
                             set(a.name for a in a2.bond_partners))
        self.assertEqual(len(top2.adjusts), len(top3.adjusts))
        for adj1, adj2 in zip(top2.adjusts, top3.adjusts):
            self.assertEqual({adj1.atom1.idx, adj1.atom2.idx}, {adj2.atom1.idx, adj2.atom2.idx})
            self.assertEqual(adj1.type, adj2.type)

        # We can't combine molecules that don't exist
        self.assertRaises(IndexError, lambda: top.write(fn, combine=[[1, 2]]))
        # Now change all angle types to urey-bradleys and make sure they're
        # written with 0s for those parameters
        psf = load_file(get_fn('ala_ala_ala.psf'))
        psf.load_parameters(
            CharmmParameterSet(get_fn('top_all22_prot.inp'), get_fn('par_all22_prot.inp'))
        )
        for atom in psf.atoms:
            self.assertIsNot(atom.atom_type, UnassignedAtomType)
        ctop = GromacsTopologyFile.from_structure(psf)
        for atom in ctop.atoms:
            self.assertIsNot(atom.atom_type, UnassignedAtomType)
            self.assertIsInstance(atom.type, str)
        ctop.write(fn, parameters=fn)
        top2 = load_file(fn)
        self.assertGreater(len(top2.urey_bradleys), 0)
        self.assertEqual(len(top2.urey_bradleys), len(ctop.urey_bradleys))

    def test_gromacs_itp_write(self):
        """ Tests the GromacsTopologyFile writer with itp option """
        def total_diheds(dlist):
            n = 0
            for d in dlist:
                if isinstance(d.type, DihedralTypeList):
                    n += len(d.type)
                elif not d.improper:
                    n += 1
            return n
        parm = load_file(get_fn('ash.parm7'))
        top = GromacsTopologyFile.from_structure(parm)
        # Write parameters and topology to same filename
        fn = self.get_fn('test.itp', written=True)
        top.write(fn, parameters=fn, itp=True)
        top2 = load_file(fn, parametrize=False)
        self.assertEqual(len(top2.atoms), len(top.atoms))
        self.assertEqual(len(top2.bonds), len(top.bonds))
        self.assertEqual(len(top2.angles), len(top.angles))
        self.assertEqual(total_diheds(top2.dihedrals), total_diheds(top.dihedrals))
        for a1, a2 in zip(top2.atoms, top.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(set(a.name for a in a1.bond_partners),
                             set(a.name for a in a2.bond_partners))
        # Now try passing open files
        with open(fn, 'w') as f:
            top.write(f, parameters=f, itp=True)
        top2 = load_file(fn, parametrize=False)
        self.assertEqual(len(top2.atoms), len(top.atoms))
        self.assertEqual(len(top2.bonds), len(top.bonds))
        self.assertEqual(len(top2.angles), len(top.angles))
        self.assertEqual(total_diheds(top2.dihedrals), total_diheds(top.dihedrals))
        for a1, a2 in zip(top2.atoms, top.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(set(a.name for a in a1.bond_partners),
                             set(a.name for a in a2.bond_partners))
        # Now try separate parameter/topology file
        fn2 = self.get_fn('test.itp', written=True)
        top.write(fn, parameters=fn2, itp=True)
        top2 = load_file(fn, parametrize=False)
        self.assertEqual(len(top2.atoms), len(top.atoms))
        self.assertEqual(len(top2.bonds), len(top.bonds))
        self.assertEqual(len(top2.angles), len(top.angles))
        self.assertEqual(total_diheds(top2.dihedrals), total_diheds(top.dihedrals))
        for a1, a2 in zip(top2.atoms, top.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(set(a.name for a in a1.bond_partners),
                             set(a.name for a in a2.bond_partners))
        # Now try separate parameter/topology/molfile files
        fn3 = self.get_fn('test_mol.itp', written=True)
        top.write(fn, parameters=fn2, molfile=fn3, itp=True)
        top2 = load_file(fn, parametrize=False)
        self.assertEqual(len(top2.atoms), len(top.atoms))
        self.assertEqual(len(top2.bonds), len(top.bonds))
        self.assertEqual(len(top2.angles), len(top.angles))
        self.assertEqual(total_diheds(top2.dihedrals), total_diheds(top.dihedrals))
        for a1, a2 in zip(top2.atoms, top.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(set(a.name for a in a1.bond_partners),
                             set(a.name for a in a2.bond_partners))

    def test_nonbond_params(self):
        """ Test the reading  and writing of the `nonbond_params` directive """
        top = load_file(get_fn('nonbond_params.top'))
        self.assertTrue(top.has_NBFIX())
        top.write(get_fn('test_nonbond_params.top'))
        top_test = load_file(get_fn('test_nonbond_params.top'))
        self.assertTrue(top_test.has_NBFIX())

        assert np.allclose(
            top.parameterset.nbfix_types[('opls_136', 'opls_135')],
            top_test.parameterset.nbfix_types[('opls_136', 'opls_135')],
        )

    def test_private_functions(self):
        """ Tests private helper functions for GromacsTopologyFile """
        Defaults = gmx.gromacstop._Defaults
        self.assertRaises(ValueError, lambda: Defaults(nbfunc=3))
        self.assertRaises(ValueError, lambda: Defaults(comb_rule=4))
        self.assertRaises(ValueError, lambda: Defaults(gen_pairs='nada'))
        self.assertRaises(ValueError, lambda: Defaults(fudgeLJ=-1.0))
        self.assertRaises(ValueError, lambda: Defaults(fudgeQQ=-1.0))
        repr(Defaults()) # To make sure it doesn't crash
        defaults = Defaults(nbfunc=2, comb_rule=3, gen_pairs='no', fudgeLJ=2.0, fudgeQQ=1.5)
        self.assertEqual(defaults[0], 2)
        self.assertEqual(defaults[1], 3)
        self.assertEqual(defaults[2], 'no')
        self.assertEqual(defaults[3], 2.0)
        self.assertEqual(defaults[4], 1.5)
        self.assertEqual(defaults[-1], 1.5)
        self.assertEqual(defaults[-2], 2.0)
        self.assertEqual(defaults[-3], 'no')
        self.assertEqual(defaults[-4], 3)
        self.assertEqual(defaults[-5], 2)
        self.assertRaises(IndexError, lambda: defaults[-6])
        self.assertRaises(IndexError, lambda: defaults[5])
        # Try setting as an array
        defaults[0] = defaults[1] = 1
        defaults[2] = 'yes'
        defaults[3] = 1.5
        defaults[4] = 2.0
        self.assertEqual(defaults.nbfunc, 1)
        self.assertEqual(defaults.comb_rule, 1)
        self.assertEqual(defaults.gen_pairs, 'yes')
        self.assertEqual(defaults.fudgeLJ, 1.5)
        self.assertEqual(defaults.fudgeQQ, 2.0)
        defaults[-5] = defaults[-4] = 1
        defaults[-3] = 'yes'
        defaults[-2] = 1.5
        defaults[-1] = 2.0
        self.assertEqual(defaults.nbfunc, 1)
        self.assertEqual(defaults.comb_rule, 1)
        self.assertEqual(defaults.gen_pairs, 'yes')
        self.assertEqual(defaults.fudgeLJ, 1.5)
        self.assertEqual(defaults.fudgeQQ, 2.0)
        # Error checking
        for idx, illegal_value in ((0, 0), (1, 0), (2, 'nada'), (3, -1), (4, -1)):
            with self.assertRaises(ValueError):
                defaults[idx] = illegal_value
        with self.assertRaises(IndexError):
            defaults[5] = 0
        
    _equal_atoms = utils.equal_atoms

@unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without Gromacs")
class TestGromacsMissingParameters(FileIOTestCase):
    """ Test handling of missing parameters """

    def setUp(self):
        self.top = load_file(get_fn('ildn.solv.top'), parametrize=False)
        FileIOTestCase.setUp(self)

    def test_missing_pairtypes(self):
        """ Tests handling of missing pairtypes parameters """
        self.top.defaults.gen_pairs = 'no'
        self.assertRaises(ParameterError, self.top.parametrize)

    def test_missing_bondtypes(self):
        """ Tests handling of missing bondtypes parameters """
        b1 = self.top.bonds[0]
        del self.top.parameterset.bond_types[(b1.atom1.type, b1.atom2.type)]
        del self.top.parameterset.bond_types[(b1.atom2.type, b1.atom1.type)]
        self.assertRaises(ParameterError, self.top.parametrize)

    def test_extra_pairs(self):
        """ Tests warning if "extra" exception pair found """
        self.top.adjusts.append(NonbondedException(self.top[0], self.top[-1]))
        with self.assertWarns(GromacsWarning):
            self.top.parametrize()

    def test_missing_angletypes(self):
        """ Tests handling of missing angletypes parameters """
        a1 = self.top.angles[0]
        key = (a1.atom1.type, a1.atom2.type, a1.atom3.type)
        del self.top.parameterset.angle_types[key]
        if key != tuple(reversed(key)):
            del self.top.parameterset.angle_types[tuple(reversed(key))]
        self.assertRaises(ParameterError, self.top.parametrize)

    def test_missing_wildcard_dihedraltypes(self):
        """ Tests handling of wild-card dihedraltypes parameters """
        def get_key(d, wc=None):
            if wc is None:
                return (d.atom1.type, d.atom2.type, d.atom3.type, d.atom4.type)
            if wc == 0:
                return ('X', d.atom2.type, d.atom3.type, d.atom4.type)
            if wc == 3:
                return (d.atom1.type, d.atom2.type, d.atom3.type, 'X')
            else:
                return ('X', d.atom2.type, d.atom3.type, 'X')
        d1 = self.top.dihedrals[0]
        for d in self.top.dihedrals:
            if get_key(d) == get_key(d1): continue
            if get_key(d, 0) == get_key(d1, 0): continue
            if get_key(d, 3) == get_key(d1, 3): continue
            if get_key(d, 0) == get_key(d1, 3): continue
            if get_key(d1, 0) == get_key(d, 3): continue
            if d.improper: continue
            if d.type is not None: continue
            d2 = d
            break
        else:
            assert False, 'Bad test parm'
        # Now make sure the two dihedrals match where only one wild-card is
        # present
        params = self.top.parameterset
        if get_key(d1) in params.dihedral_types:
            del params.dihedral_types[get_key(d1)]
            del params.dihedral_types[tuple(reversed(get_key(d1)))]
        if get_key(d2) in params.dihedral_types:
            del params.dihedral_types[get_key(d2)]
            del params.dihedral_types[tuple(reversed(get_key(d2)))]
        dt1 = DihedralTypeList([DihedralType(10, 180, 1)])
        dt2 = DihedralTypeList([DihedralType(11, 0, 2)])
        params.dihedral_types[get_key(d1, 0)] = dt1
        params.dihedral_types[tuple(reversed(get_key(d1, 0)))] = dt1
        params.dihedral_types[get_key(d2, 3)] = dt2
        params.dihedral_types[tuple(reversed(get_key(d2, 3)))] = dt2

        self.top.parametrize()
        self.assertEqual(d1.type, dt1)
        self.assertEqual(d2.type, dt2)

    def test_missing_dihedraltypes(self):
        """ Tests handling of missing dihedraltypes parameters """
        def get_key(d, wc=None):
            if wc is None:
                return (d.atom1.type, d.atom2.type, d.atom3.type, d.atom4.type)
            if wc == 0:
                return ('X', d.atom2.type, d.atom3.type, d.atom4.type)
            if wc == 3:
                return (d.atom1.type, d.atom2.type, d.atom3.type, 'X')
            else:
                return ('X', d.atom2.type, d.atom3.type, 'X')
        for d in self.top.dihedrals:
            if d.type is not None: continue
            break
        params = self.top.parameterset
        if get_key(d) in params.dihedral_types:
            del params.dihedral_types[get_key(d)]
            del params.dihedral_types[tuple(reversed(get_key(d)))]
        if get_key(d, wc=100) in params.dihedral_types:
            del params.dihedral_types[get_key(d, wc=100)]
            del params.dihedral_types[tuple(reversed(get_key(d, wc=100)))]
        self.assertRaises(ParameterError, self.top.parametrize)

    def test_missing_impropertypes(self):
        """ Tests handling of missing improper type """
        for key in set(self.top.parameterset.improper_periodic_types.keys()):
            del self.top.parameterset.improper_periodic_types[key]
        self.assertRaises(ParameterError, self.top.parametrize)

    def test_wildcard_rbtorsions(self):
        """ Tests handling of missing and wild-cards with R-B torsion types """
        def get_key(d, wc=None):
            if wc is None:
                return (d.atom1.type, d.atom2.type, d.atom3.type, d.atom4.type)
            if wc == 0:
                return ('X', d.atom2.type, d.atom3.type, d.atom4.type)
            if wc == 3:
                return (d.atom1.type, d.atom2.type, d.atom3.type, 'X')
            else:
                return ('X', d.atom2.type, d.atom3.type, 'X')
        for i, d1 in enumerate(self.top.dihedrals):
            if not d1.improper and d1.type is None:
                break
        else:
            assert False, 'Bad topology file for test'
        del self.top.dihedrals[i]
        for i, d2 in enumerate(self.top.dihedrals):
            if get_key(d1) == get_key(d2): continue
            if get_key(d1, 0) == get_key(d2, 0): continue
            if get_key(d1, 3) == get_key(d2, 3): continue
            if get_key(d1, 0) == get_key(d2, 3): continue
            if get_key(d2, 0) == get_key(d1, 3): continue
            if not d2.improper and d2.type is None:
                break
        else:
            assert False, 'Bad topology file for test'
        del self.top.dihedrals[i]
        self.top.rb_torsions.extend([d1, d2])
        self.assertRaises(ParameterError, self.top.parametrize)
        # Now assign wild-cards
        params = self.top.parameterset
        rbt = RBTorsionType(1, 2, 3, 4, 5, 6)
        params.rb_torsion_types[get_key(d1, 0)] = rbt
        params.rb_torsion_types[tuple(reversed(get_key(d1, 0)))] = rbt
        rbt2 = RBTorsionType(2, 3, 4, 5, 6, 7)
        params.rb_torsion_types[get_key(d2, 1)] = rbt2
        params.rb_torsion_types[tuple(reversed(get_key(d2, 1)))] = rbt2

        self.top.parametrize()

        self.assertEqual(d1.type, rbt)
        self.assertEqual(d2.type, rbt2)

    def test_missing_impropers(self):
        """ Test handling of missing impropers """
        self.top.impropers.append(Improper(*tuple(self.top.atoms[:4])))
        self.assertRaises(ParameterError, self.top.parametrize)

    def test_missing_cmaps(self):
        """ Test handling of missing cmaptypes """
        self.top.cmaps.append(Cmap(*tuple(self.top.atoms[:5])))
        self.assertRaises(ParameterError, self.top.parametrize)

    def test_missing_ureybradleys(self):
        """ Test handling of missing Urey-Bradley types """
        self.top.angles[0].funct = 5
        self.top.urey_bradleys.append(
            UreyBradley(self.top.angles[0].atom1, self.top.angles[0].atom3)
        )
        self.assertRaises(ParameterError, self.top.parametrize)

class TestGromacsTopHelperFunctions(FileIOTestCase):
    """ Test GROMACS helper functions """

    def setUp(self):
        self.top = GromacsTopologyFile()
        self.top.add_atom(Atom(name='C1'), 'ABC', 1)
        self.top.add_atom(Atom(name='C1'), 'DEF', 2)
        self.top.add_atom(Atom(name='C1'), 'GHI', 3)
        self.top.add_atom(Atom(name='C1'), 'JKL', 4)
        self.top.add_atom(Atom(name='C1'), 'MNO', 5)
        self.top.add_atom(Atom(name='C1'), 'PQR', 5)
        FileIOTestCase.setUp(self)

    def test_parse_pairs(self):
        """ Test GromacsTopologyFile._parse_pairs """
        with self.assertWarns(GromacsWarning):
            self.top._parse_pairs('1  2  3\n', dict(), self.top.atoms)
        self.top._parse_pairs('1  2  3\n', dict(), self.top.atoms)
        self.assertTrue(self.top.unknown_functional)

    def test_parse_angles(self):
        """ Test GromacsTopologyFile._parse_angles """
        with self.assertWarns(GromacsWarning):
            self.top._parse_angles('1  2  3  2\n', dict(), dict(), self.top.atoms)
        self.top._parse_angles('1  2  3  2\n', dict(), dict(), self.top.atoms)
        self.assertTrue(self.top.unknown_functional)

    def test_parse_dihedrals(self):
        """ Test GromacsTopologyFile._parse_dihedrals """
        with self.assertWarns(GromacsWarning):
            self.top._parse_dihedrals('1 2 3 4 100\n', dict(), dict(), self.top)
        self.assertTrue(self.top.unknown_functional)
        self.assertEqual(len(self.top.dihedrals), 1)
        dih = self.top.dihedrals[0]
        self.assertIs(dih.atom1, self.top[0])
        self.assertIs(dih.atom2, self.top[1])
        self.assertIs(dih.atom3, self.top[2])
        self.assertIs(dih.atom4, self.top[3])
        self.assertIs(dih.type, None)
        # Test multi-term dihedrals
        dt = dict()
        PMD = dict()
        self.top._parse_dihedrals('1 2 3 4 9 180 50 1', dt, PMD, self.top)
        self.assertIn(tuple(self.top.atoms[:4]), PMD)
        self.top._parse_dihedrals('1 2 3 4 9 180 40 2', dt, PMD, self.top)
        self.assertEqual(len(PMD[tuple(self.top.atoms[:4])]), 2)

    def test_parse_cmaps(self):
        """ Test GromacsTopologyFile._parse_cmaps """
        with self.assertWarns(GromacsWarning):
            self.top._parse_cmaps('1 2 3 4 5 2\n', self.top.atoms)
        self.top._parse_cmaps('1 2 3 4 5 2\n', self.top.atoms)
        self.assertTrue(self.top.unknown_functional)

    def test_parse_settles(self):
        """ Test GromacsTopologyFile._parse_settles """
        with self.assertRaises(GromacsError):
            self.top._parse_settles('whatever', self.top.atoms)
        with self.assertRaises(GromacsError):
            self.top._parse_settles('whatever', self.top.atoms[:3])
        self.top[0].atomic_number = 8
        self.top[1].atomic_number = self.top[2].atomic_number = 1
        with self.assertRaises(GromacsError):
            self.top._parse_settles('1 2 nofloat nofloat\n', self.top.atoms[:3])

    def test_parse_vsite3(self):
        """ Test GromacsTopologyFile._parse_vsites3 """
        with self.assertRaises(GromacsError):
            self.top._parse_vsites3('1 2 3 4 1 1.2 1.3\n', self.top.atoms, ParameterSet())
        with self.assertRaises(GromacsError):
            self.top._parse_vsites3('1 2 3 4 2 1.2 1.3\n', self.top.atoms, ParameterSet())
        bond = Bond(self.top[0], self.top[1])
        with self.assertRaises(GromacsError):
            self.top._parse_vsites3('1 2 3 4 1 1.2 1.2', self.top.atoms, ParameterSet())

    def test_parse_atomtypes(self):
        """ Test GromacsTopologyFile._parse_atomtypes """
        name, typ = self.top._parse_atomtypes('CX 12.01 0 A 0.1 2.0')
        self.assertEqual(name, 'CX')
        self.assertEqual(typ.atomic_number, 6)
        self.assertEqual(typ.charge, 0)
        self.assertEqual(typ.epsilon, 2.0/4.184)
        self.assertEqual(typ.sigma, 1)

    def test_parse_bondtypes(self):
        """ Test GromacsTopologyFile._parse_bondtypes """
        with self.assertWarns(GromacsWarning):
            self.top._parse_bondtypes('CA CB 2 0.1 5000')
        self.top._parse_bondtypes('CA CB 2 0.1 5000')
        self.assertTrue(self.top.unknown_functional)

    def test_parse_angletypes(self):
        """ Test GromacsTopologyFile._parse_angletypes """
        with self.assertWarns(GromacsWarning):
            self.top._parse_angletypes('CA CB CC 2 120 5000')
        self.top._parse_angletypes('CA CB CC 2 120 5000')
        self.assertTrue(self.top.unknown_functional)

    def test_parse_dihedraltypes(self):
        """ Test GromacsTopologyFile._parse_dihedraltypes """
        key, dtype, ptype, replace = self.top._parse_dihedraltypes('CA CA 9 180 50.0 2')
        self.assertEqual(key, ('X', 'CA', 'CA', 'X'))
        self.assertEqual(dtype, 'normal')
        self.assertFalse(replace)
        self.assertEqual(ptype.phase, 180)
        self.assertEqual(ptype.phi_k, 50/4.184)
        self.assertEqual(ptype.per, 2)
        with self.assertWarns(GromacsWarning):
            self.top._parse_dihedraltypes('CX CA CA CX 10 180 50.0 2')
        self.top._parse_dihedraltypes('CX CA CA CX 10 180 50.0 2')
        self.assertTrue(self.top.unknown_functional)

    def test_parse_cmaptypes(self):
        """ Test GromacsTopologyFile._parse_cmaptypes """
        self.assertRaises(GromacsError, lambda:
                self.top._parse_cmaptypes('C1 C2 C3 C4 C5 1 24 24 1 2 3 4 5'))
        self.assertRaises(GromacsError, lambda:
                self.top._parse_cmaptypes('C1 C2 C3 C4 C5 1 2 3 1 2 3 4 5 6'))

class TestGromacsGro(FileIOTestCase):
    """ Tests the Gromacs GRO file parser """

    def test_gro_velocities(self):
        """ Test parsing and writing GROMACS GRO files with velocities """
        gro = load_file(get_fn('1aki.ff99sbildn.gro'))
        self.assertIsInstance(gro, Structure)
        vel = np.random.rand(len(gro.atoms), 3)
        gro.velocities = vel
        fn = self.get_fn('test.gro', written=True)
        gro.save(fn)
        self.assertTrue(GromacsGroFile.id_format(fn))
        gro.save(fn, overwrite=True, precision=8)
        self.assertTrue(GromacsGroFile.id_format(fn))
        with open(fn) as f:
            tmp = GromacsGroFile.parse(f)
        np.testing.assert_allclose(vel, tmp.velocities, atol=1e-8)

    def test_gro_detection(self):
        """ Tests automatic detection of GROMACS GRO files """
        fn = self.get_fn('candidate.gro', written=True)
        with open(fn, 'w') as f:
            f.write('Some title\n 1000\n    aISNot a valid format\n')
        self.assertFalse(GromacsGroFile.id_format(fn))
        self.assertRaises(GromacsError, lambda: GromacsGroFile.parse(fn))
        f = StringIO('Gromacs title line\n notanumber\nsome line\n')
        self.assertRaises(GromacsError, lambda: GromacsGroFile.parse(f))

    def test_gro_elements_bonds(self):
        """ Tests improved element and bond assignment in GRO files """
        gro = GromacsGroFile.parse(
                os.path.join(get_fn('07.DHFR-Liquid-NoPBC'), 'conf.gro')
        )
        self.assertGreater(len(gro.bonds), 0)
        for atom in gro.view['NA', :].atoms:
            self.assertEqual(atom.atomic_number, periodic_table.AtomicNum['Na'])

    def test_read_gro_file(self):
        """ Tests reading GRO file """
        gro = GromacsGroFile.parse(get_fn('1aki.ff99sbildn.gro'))
        self.assertIsInstance(gro, Structure)
        self.assertEqual(len(gro.atoms), 1960)
        self.assertEqual(len(gro.residues), 129)
        self.assertAlmostEqual(gro.atoms[0].xx, 44.6)
        self.assertAlmostEqual(gro.atoms[0].xy, 49.86)
        self.assertAlmostEqual(gro.atoms[0].xz, 18.10)
        self.assertAlmostEqual(gro.atoms[1959].xx, 50.97)
        self.assertAlmostEqual(gro.atoms[1959].xy, 39.80)
        self.assertAlmostEqual(gro.atoms[1959].xz, 38.64)
        self.assertAlmostEqual(gro.box[0], 74.1008)
        self.assertAlmostEqual(gro.box[1], 74.10080712)
        self.assertAlmostEqual(gro.box[2], 74.10074585)
        self.assertAlmostEqual(gro.box[3], 70.52882666)
        self.assertAlmostEqual(gro.box[4], 109.47126278)
        self.assertAlmostEqual(gro.box[5], 70.52875398)
        # Check atomic number and mass assignment
        self.assertEqual(gro.atoms[0].atomic_number, 7)
        self.assertEqual(gro.atoms[0].mass, 14.0067)
        fn = self.get_fn('test.gro', written=True)
        # Test bad GRO files
        with open(fn, 'w') as wf, open(get_fn('1aki.charmm27.solv.gro')) as f:
            for i in range(1000):
                wf.write(f.readline())
        self.assertRaises(GromacsError, lambda: GromacsGroFile.parse(fn))
        with open(get_fn('1aki.ff99sbildn.gro')) as f:
            lines = f.readlines()
        lines[-1] = 'not a legal box line\n'
        with open(fn, 'w') as f:
            f.write(''.join(lines))
        self.assertRaises(GromacsError, lambda: GromacsGroFile.parse(fn))

    def test_write_gro_file(self):
        """ Tests writing GRO file """
        gro = GromacsGroFile.parse(get_fn('1aki.ff99sbildn.gro'))
        GromacsGroFile.write(gro, self.get_fn('1aki.ff99sbildn.gro', written=True))
        gro = load_file(self.get_fn('1aki.ff99sbildn.gro', written=True))
        self.assertIsInstance(gro, Structure)
        self.assertEqual(len(gro.atoms), 1960)
        self.assertEqual(len(gro.residues), 129)
        self.assertAlmostEqual(gro.atoms[0].xx, 44.6)
        self.assertAlmostEqual(gro.atoms[0].xy, 49.86)
        self.assertAlmostEqual(gro.atoms[0].xz, 18.10)
        self.assertAlmostEqual(gro.atoms[1959].xx, 50.97)
        self.assertAlmostEqual(gro.atoms[1959].xy, 39.80)
        self.assertAlmostEqual(gro.atoms[1959].xz, 38.64)
        self.assertAlmostEqual(gro.box[0], 74.1008)
        self.assertAlmostEqual(gro.box[1], 74.10080712)
        self.assertAlmostEqual(gro.box[2], 74.10074585)
        self.assertAlmostEqual(gro.box[3], 70.52882666)
        self.assertAlmostEqual(gro.box[4], 109.47126278)
        self.assertAlmostEqual(gro.box[5], 70.52875398)
        self.assertRaises(TypeError, lambda: GromacsGroFile.write(gro, 10))

    def test_write_gro_file_nobox(self):
        """ Test GROMACS GRO file writing without a box """
        parm = load_file(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        self.assertIs(parm.box, None)
        fn = self.get_fn('test.gro', written=True)
        parm.save(fn)
        # Make sure it has a box
        gro = load_file(fn)
        self.assertIsNot(gro.box, None)
        parm.save(fn, nobox=True, overwrite=True)
        gro = load_file(fn)
        self.assertIs(gro.box, None)

    def test_read_write_high_precision_gro_file(self):
        """ Tests reading/writing high-precision GRO files """
        gro = GromacsGroFile.parse(get_fn('1aki.ff99sbildn.gro'))
        GromacsGroFile.write(gro, self.get_fn('1aki.ff99sbildn_highprec.gro', written=True), precision=6)
        self.assertTrue(diff_files(self.get_fn('1aki.ff99sbildn_highprec.gro', saved=True),
                                   self.get_fn('1aki.ff99sbildn_highprec.gro', written=True))
        )
        gro2 = GromacsGroFile.parse(self.get_fn('1aki.ff99sbildn_highprec.gro', written=True))
        gro3 = load_file(self.get_fn('1aki.ff99sbildn_highprec.gro', written=True))
        self.assertIsInstance(gro3, Structure)
        for a1, a2, a3 in zip(gro.atoms, gro2.atoms, gro3.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a3.name, a2.name)
            self.assertEqual(a1.xx, a2.xx)
            self.assertEqual(a3.xx, a2.xx)
            self.assertEqual(a1.xy, a2.xy)
            self.assertEqual(a3.xy, a2.xy)
            self.assertEqual(a1.xz, a2.xz)
            self.assertEqual(a3.xz, a2.xz)

    def test_write_gro_matching_topology(self):
        """Tests the usage of the match_topology and combine keyword arguments
        to GromacsGroFile.write that can be used to write a Gro file with atom
        ordering that always matches a GromacsTopologyFile written from the
        same structure.
        """
        def compare_top_gro_atom_order(struct, combine):
            """From a given Structure write a Gro and topology file and check
            that the atom order is the same regardless of the value of combine
            """
            # write out files into StringIO's to trigger atom reordering
            topology = GromacsTopologyFile.from_structure(struct)
            top_file = StringIO()
            topology.write(top_file, combine=combine)
            gro_file = StringIO()
            GromacsGroFile.write(struct, gro_file, combine=combine)

            # parse StringIO's back into parmed objects
            top_file.seek(0)
            gro_file.seek(0)
            top = GromacsTopologyFile()
            top.read(top_file, parametrize=False)
            gro = GromacsGroFile.parse(gro_file, skip_bonds=True)
            for top_at, gro_at in zip(top.atoms, gro.atoms):
                assert top_at.name == gro_at.name

            return top, gro

        # Build a very simple 7 atom structure
        # Bonds will be added between 'A' and 'D' atoms later
        struct = Structure()
        for i, (letter, chain) in enumerate(zip('ABCDABD', 'ABCDEFG')):
            atom = Atom(name=letter, type=letter, atomic_number=6, mass=12.01)
            # coordinate info needed to create gro file
            atom.xx, atom.xy, atom.xz = float(i), 0., 0.
            struct.add_atom(atom, resname=letter, resnum=i, chain=chain)

        # test the possible combinations of the combine argument
        # in a system that IS NOT reordered when writing the topology
        combine_values = None, 'all', [(0, 1)]
        for combine in combine_values:
            compare_top_gro_atom_order(struct, combine=combine)

        # adding bonds will trigger atom reordering when combine=None
        struct.bonds.append(Bond(struct.atoms[0], struct.atoms[6]))
        struct.bonds.append(Bond(struct.atoms[4], struct.atoms[3]))

        for combine in combine_values:
            top, gro = compare_top_gro_atom_order(struct, combine=combine)
            if combine is None:
                # check that each 'A' particle has been paired with the correct
                # 'D' particle by checking the unique xx coordinate
                self.assertEqual(gro.atoms[0].xx, 0.)
                self.assertEqual(gro.atoms[1].xx, 6.)
                self.assertEqual(gro.atoms[4].xx, 3.)
                self.assertEqual(gro.atoms[5].xx, 4.)

        # add some extra bonds that change reordering behaviour
        struct.bonds.append(Bond(struct.atoms[0], struct.atoms[1]))
        struct.bonds.append(Bond(struct.atoms[4], struct.atoms[5]))

        for combine in combine_values:
            top, gro = compare_top_gro_atom_order(struct, combine=combine)
            if combine is None:
                # check that each 'A' particle has been paired with the correct
                # 'D' particle by checking the unique xx coordinate
                self.assertEqual(gro.atoms[0].xx, 0.)
                self.assertEqual(gro.atoms[2].xx, 6.)
                self.assertEqual(gro.atoms[4].xx, 3.)
                self.assertEqual(gro.atoms[5].xx, 4.)
