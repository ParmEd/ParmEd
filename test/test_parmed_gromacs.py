"""
Tests the functionality in the parmed.gromacs package
"""
import utils
from parmed import load_file, Structure, ExtraPoint, DihedralTypeList
from parmed.exceptions import GromacsWarning
from parmed.gromacs import GromacsTopologyFile, GromacsGroFile
from parmed import gromacs as gmx
from parmed.utils.six.moves import range, zip, StringIO
import os
import unittest
from utils import get_fn, diff_files, get_saved_fn, FileIOTestCase
import warnings

@unittest.skipIf(not os.path.exists(gmx.GROMACS_TOPDIR), "Cannot run GROMACS tests without Gromacs")
class TestGromacsTop(FileIOTestCase):
    """ Tests the Gromacs topology file parser """

    def setUp(self):
        warnings.filterwarnings('error', category=GromacsWarning)
        FileIOTestCase.setUp(self)

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

    def testCharmm27Top(self):
        """ Tests parsing a Gromacs topology with CHARMM 27 FF """
        top = GromacsTopologyFile(get_fn('1aki.charmm27.top'))
        self.assertEqual(top.combining_rule, 'lorentz')
        self.assertEqual(top.itps, ['charmm27.ff/forcefield.itp',
                                    'charmm27.ff/tip3p.itp',
                                    'charmm27.ff/ions.itp'])
        self._charmm27_checks(top)

    def testWriteCharmm27Top(self):
        """ Tests writing a Gromacs topology file with CHARMM 27 FF """
        top = load_file(get_fn('1aki.charmm27.top'))
        self.assertEqual(top.combining_rule, 'lorentz')
        GromacsTopologyFile.write(top,
                get_fn('1aki.charmm27.top', written=True))
        top2 = load_file(get_fn('1aki.charmm27.top', written=True))
        self._charmm27_checks(top)

    def _check_ff99sbildn(self, top):
        self.assertEqual(len(top.atoms), 40560)
        self.assertEqual(len(top.residues), 9779)
        self.assertEqual(len([a for a in top.atoms if isinstance(a, ExtraPoint)]),
                         9650)
        self.assertEqual(len(top.bonds), 30934)
        self.assertEqual(len(top.angles), 13197)
        self.assertEqual(len(top.dihedrals), 5613)

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

    def testReadAmber99SBILDN(self):
        """ Tests parsing a Gromacs topology with Amber99SBILDN and water """
        top = load_file(get_fn('1aki.ff99sbildn.top'))
        self.assertEqual(top.combining_rule, 'lorentz')
        self._check_ff99sbildn(top)

    def testWriteAmber99SBILDN(self):
        """ Tests writing a Gromacs topology with multiple molecules """
        top = load_file(get_fn('1aki.ff99sbildn.top'))
        self.assertEqual(top.combining_rule, 'lorentz')
        GromacsTopologyFile.write(top,
                get_fn('1aki.ff99sbildn.top', written=True),
                combine=None)
        top2 = load_file(get_fn('1aki.ff99sbildn.top', written=True))
        self._check_ff99sbildn(top2)
        self._check_equal_structures(top, top2)

    def testDuplicateSystemNames(self):
        """ Tests that Gromacs topologies never have duplicate moleculetypes """
        parm = load_file(get_fn('phenol.prmtop'))
        parm = parm * 20 + load_file(get_fn('biphenyl.prmtop')) * 20
        top = GromacsTopologyFile.from_structure(parm)
        self.assertEqual(top.combining_rule, 'lorentz')
        top.write(get_fn('phenol_biphenyl.top', written=True))
        top2 = GromacsTopologyFile(get_fn('phenol_biphenyl.top', written=True))
        self.assertEqual(len(top.residues), 40)

    def testOPLS(self):
        """ Tests the geometric combining rules in Gromacs with OPLS/AA """
        parm = load_file(os.path.join(get_fn('05.OPLS'), 'topol.top'),
                         xyz=os.path.join(get_fn('05.OPLS'), 'conf.gro'))
        self.assertEqual(parm.combining_rule, 'geometric')
        self.assertEqual(parm.defaults.comb_rule, 3)
        parm.write(get_fn('test.topol', written=True), combine='all')
        parm2 = load_file(get_fn('test.topol', written=True))
        self.assertEqual(len(parm.atoms), len(parm2.atoms))

    def testMoleculeOrdering(self):
        """ Tests non-contiguous atoms in Gromacs topology file writes """
        warnings.filterwarnings('ignore', category=GromacsWarning)
        parm = load_file(os.path.join(get_fn('12.DPPC'), 'topol3.top'))
        parm.write(get_fn('topol3.top', written=True))
        parm2 = load_file(get_fn('topol3.top', written=True))
        self.assertEqual(len(parm.atoms), len(parm2.atoms))
        self.assertEqual(len(parm.residues), len(parm2.residues))
        for r1, r2 in zip(parm.residues, parm2.residues):
            self.assertEqual(r1.name, r2.name)
            for a1, a2 in zip(r1.atoms, r2.atoms):
                self.assertEqual(a1.name, a2.name)
                self.assertEqual(a1.type, a2.type)

class TestGromacsGro(FileIOTestCase):
    """ Tests the Gromacs GRO file parser """

    def testReadGroFile(self):
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

    def testWriteGroFile(self):
        """ Tests writing GRO file """
        gro = GromacsGroFile.parse(get_fn('1aki.ff99sbildn.gro'))
        GromacsGroFile.write(gro, get_fn('1aki.ff99sbildn.gro', written=True))
        gro = load_file(get_fn('1aki.ff99sbildn.gro', written=True))
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

    def testReadWriteHighPrecisionGroFile(self):
        """ Tests reading/writing high-precision GRO files """
        gro = GromacsGroFile.parse(get_fn('1aki.ff99sbildn.gro'))
        GromacsGroFile.write(gro, get_fn('1aki.ff99sbildn_highprec.gro',
                                         written=True),
                             precision=6)
        self.assertTrue(diff_files(get_saved_fn('1aki.ff99sbildn_highprec.gro'),
                                   get_fn('1aki.ff99sbildn_highprec.gro',
                                          written=True)
                                   )
        )
        gro2 = GromacsGroFile.parse(get_fn('1aki.ff99sbildn_highprec.gro',
                                           written=True))
        gro3 = load_file(get_fn('1aki.ff99sbildn_highprec.gro', written=True))
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
