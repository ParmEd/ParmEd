"""
Tests the functionality in the chemistry.gromacs package
"""
from chemistry import load_file, Structure, ExtraPoint
from chemistry.exceptions import PreProcessorError, PreProcessorWarning
from chemistry.gromacs import GromacsTopologyFile, GromacsGroFile
from chemistry.utils.six.moves import range, zip, StringIO
import os
import unittest
from utils import get_fn, diff_files, get_saved_fn
import warnings

class TestGromacsTop(unittest.TestCase):
    """ Tests the Gromacs topology file parser """

    def setUp(self):
        try:
            os.makedirs(get_fn('writes'))
        except OSError:
            pass

    def tearDown(self):
        try:
            for f in os.listdir(get_fn('writes')):
                os.unlink(get_fn(f, written=True))
            os.rmdir(get_fn('writes'))
        except OSError:
            pass

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
        top, itps = GromacsTopologyFile.parse(get_fn('1aki.charmm27.top'),
                                              return_itps=True)
        self.assertEqual(itps, ['charmm27.ff/forcefield.itp',
                                'charmm27.ff/tip3p.itp',
                                'charmm27.ff/ions.itp'])
        self._charmm27_checks(top)

    def testWriteCharmm27Top(self):
        """ Tests writing a Gromacs topology file with CHARMM 27 FF """
        top, param = load_file(get_fn('1aki.charmm27.top'), return_params=True)
        GromacsTopologyFile.write(top,
                get_fn('1aki.charmm27.top', written=True),
                include_itps='charmm27.ff/forcefield.itp')
        top2, param2, itp = load_file(get_fn('1aki.charmm27.top', written=True),
                                      return_params=True, return_itps=True)
        self.assertEqual(itp, ['charmm27.ff/forcefield.itp'])
        self._charmm27_checks(top)

    def _check_ff99sbildn(self, top):
        self.assertEqual(len(top.atoms), 40560)
        self.assertEqual(len(top.residues), 9779)
        self.assertEqual(len([a for a in top.atoms if isinstance(a, ExtraPoint)]),
                         9650)
        self.assertEqual(len(top.bonds), 40584)
        self.assertEqual(len(top.angles), 3547)
        self.assertEqual(len(top.dihedrals), 5613)

    def testReadAmber99SBILDN(self):
        """ Tests parsing a Gromacs topology with Amber99SBILDN and water """
        top = load_file(get_fn('1aki.ff99sbildn.top'))
        self._check_ff99sbildn(top)

    def testWriteAmber99SBILDN(self):
        """ Tests writing a Gromacs topology with multiple molecules """
        top = load_file(get_fn('1aki.ff99sbildn.top'))
        GromacsTopologyFile.write(top,
                get_fn('1aki.ff99sbildn.top', written=True),
                combine=None)
        top2 = load_file(get_fn('1aki.ff99sbildn.top', written=True))
        self._check_ff99sbildn(top2)
        self.assertTrue(diff_files(get_fn('1aki.ff99sbildn.top', written=True),
                                   get_saved_fn('1aki.ff99sbildn.top')))

class TestGromacsGro(unittest.TestCase):
    """ Tests the Gromacs GRO file parser """

    def setUp(self):
        try:
            os.makedirs(get_fn('writes'))
        except OSError:
            pass

    def tearDown(self):
        try:
            for f in os.listdir(get_fn('writes')):
                os.unlink(get_fn(f, written=True))
            os.rmdir(get_fn('writes'))
        except OSError:
            pass

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

