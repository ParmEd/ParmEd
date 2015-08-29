"""
Tests the functionality in the parmed.amber package
"""
from __future__ import print_function, division

import glob
import math
import numpy as np
import os
from parmed.amber import readparm, asciicrd, mask, parameters
from parmed.exceptions import AmberWarning
from parmed import topologyobjects, load_file
from parmed.utils.six import string_types, iteritems
from parmed.utils.six.moves import range, zip
import random
import unittest
from utils import get_fn, has_numpy, FileIOTestCase
import warnings

class TestReadParm(unittest.TestCase):
    """ Tests the various Parm file classes """
    
    def testOptimizedReader(self):
        """ Check that the optimized reader imports correctly """
        from parmed.amber import _rdparm

    def testLoadParm(self):
        """ Test the arbitrary parm loader """
        parm = readparm.LoadParm(get_fn('trx.prmtop'))
        parm2 = readparm.AmberParm(get_fn('trx.prmtop'))
        for key in parm.parm_data:
            self.assertEqual(parm.parm_data[key], parm2.parm_data[key])

    def testGzippedParm(self):
        """ Check that gzipped prmtop files can be parsed correctly """
        parm = readparm.LoadParm(get_fn('small.parm7.gz'))
        self.assertEqual(parm.ptr('natom'), 864)

    def testBzippedParm(self):
        """ Check that bzip2ed prmtop files can be parsed correctly """
        parm = readparm.LoadParm(get_fn('small.parm7.bz2'))
        self.assertEqual(parm.ptr('natom'), 864)

    def testAmberGasParm(self):
        """ Test the AmberParm class with a non-periodic (gas-phase) prmtop """
        parm = readparm.AmberParm(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        gasparm = readparm.AmberParm(get_fn('trx.prmtop'))
        gasparm.load_rst7(get_fn('trx.inpcrd'))
        self.assertEqual(gasparm.combining_rule, 'lorentz')

        self.assertEqual([a.xx for a in gasparm.atoms],
                         [a.xx for a in parm.atoms])
        self.assertEqual([a.xy for a in gasparm.atoms],
                         [a.xy for a in parm.atoms])
        self.assertEqual([a.xz for a in gasparm.atoms],
                         [a.xz for a in parm.atoms])
        
        # Now run the tests for the prmtop
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertFalse(parm.chamber)
        self.assertFalse(parm.amoeba)
        self.assertRaises(KeyError, lambda: parm.parm_data['BOX_DIMENSIONS'])
        self.assertEqual(parm.ptr('ifbox'), 0)

        # Now check the restart file
        rst = readparm.Rst7.open(get_fn('trx.inpcrd'))
        coords = rst.coordinates
        vels = rst.vels
        for i, atom in enumerate(gasparm.atoms):
            np.testing.assert_allclose(coords[0,i], [atom.xx, atom.xy, atom.xz])
            np.testing.assert_allclose(vels[0,i], [atom.vx, atom.vy, atom.vz])

    def testRemakeParm(self):
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

    def testRemakeChamberParm(self):
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

    def testAmberSolvParm(self):
        """ Test the AmberParm class with a periodic prmtop """
        parm = readparm.AmberParm(get_fn('solv.prmtop'),
                                  get_fn('solv.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        self._standard_parm_tests(parm)
        self._solv_pointer_tests(parm)
        self.assertFalse(parm.chamber)
        self.assertFalse(parm.amoeba)
        self.assertEqual(parm.ptr('ifbox'), 1)

    def testChamberGasParm(self):
        """Test the ChamberParm class with a non-periodic (gas phase) prmtop"""
        parm = readparm.ChamberParm(get_fn('ala_ala_ala.parm7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertTrue(parm.chamber)
        self.assertTrue(parm.has_cmap)
        self.assertEqual(parm.ptr('ifbox'), 0)

    def testChamberSolvParm(self):
        """ Test the ChamberParm class with a periodic prmtop """
        parm = readparm.ChamberParm(get_fn('dhfr_cmap_pbc.parm7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        self._standard_parm_tests(parm)
        self._solv_pointer_tests(parm)
        self.assertTrue(parm.chamber)
        self.assertTrue(parm.has_cmap)
        self.assertEqual(parm.ptr('ifbox'), 1)

    def testAmoebaBig(self):
        """ Test the AmoebaParm class with a large system """
        parm = readparm.AmoebaParm(get_fn('amoeba.parm7'))
        self.assertEqual(parm.ptr('natom'), len(parm.atoms))
        self.assertEqual([a.name for a in parm.atoms],
                         parm.parm_data['ATOM_NAME'])
        self.assertEqual([a.type for a in parm.atoms],
                         parm.parm_data['AMBER_ATOM_TYPE'])
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

    def testAmoebaSmall(self):
        """ Test the AmoebaParm class w/ small system (not all terms) """
        parm = readparm.AmoebaParm(get_fn('nma.parm7'))
        rst7 = readparm.BeemanRestart(get_fn('nma.rst'))
        self.assertEqual(3*rst7.natom, len(rst7.coordinates))
        self.assertEqual(rst7.coordinates, rst7.parm_data['ATOMIC_COORDS_LIST'])
        self.assertEqual(rst7.natom, parm.ptr('natom'))
        self.assertFalse(parm.torsion_torsions)
        self.assertTrue(parm.amoeba)

    def test1012(self):
        """ Test that 10-12 prmtop files are recognized properly """
        parm = readparm.AmberParm(get_fn('ff91.parm7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        self._standard_parm_tests(parm, has1012=True)

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

class TestParameterFiles(unittest.TestCase):
    """ Tests Amber parameter and frcmod files """

    def testFileDetectionFrcmod(self):
        """ Tests the detection of Amber frcmod files """
        for fname in glob.glob(os.path.join(get_fn('parm'), 'frcmod.*')):
            self.assertTrue(parameters.AmberParameterSet.id_format(fname))

    def testFileDetectionParm(self):
        """ Tests the detection of Amber parm.dat files """
        for fname in glob.glob(os.path.join(get_fn('parm'), 'parm*.dat')):
            self.assertTrue(parameters.AmberParameterSet.id_format(fname))

    def testFrcmodParsing(self):
        """ Tests parsing an Amber frcmod file """
        params = parameters.AmberParameterSet(
                os.path.join(get_fn('parm'), 'frcmod.ff99SB')
        )
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

    def testParmParsing(self):
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

    def testParmParsingLJEDIT(self):
        """ Tests parsing an Amber parm.dat file with an LJEDIT section """
        params = parameters.AmberParameterSet(
                os.path.join(get_fn('parm'), 'parm14ipq.dat')
        )
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
        self.assertEqual(params.nbfix_types[('O3', 'OW')][0],
                         math.sqrt(0.162750*0.21))
        self.assertEqual(params.nbfix_types[('O3', 'OW')][1],
                         1.775931+1.8605)
        self.assertEqual(params.nbfix_types[('OA', 'OW')][0],
                         math.sqrt(0.162750*0.2104))
        self.assertEqual(params.nbfix_types[('OA', 'OW')][1],
                         1.775931+1.66)

    @unittest.skipIf(os.getenv('AMBERHOME') is None, 'Cannot test w/out Amber')
    def testLoadLeaprc(self):
        """ Tests loading a leaprc file to define a force field """
        warnings.filterwarnings('ignore', category=AmberWarning)
        params = parameters.AmberParameterSet.from_leaprc(
                os.path.join(os.getenv('AMBERHOME'), 'dat', 'leap', 'cmd',
                             'leaprc.ff14SB')
        )
        self.assertEqual(params.atom_types['H'].atomic_number, 1)
        self.assertEqual(params.atom_types['3C'].atomic_number, 6)
        self.assertEqual(params.atom_types['K+'].atomic_number, 19)
        self.assertTrue(params.residues)

    def testParmSetParsing(self):
        """ Tests parsing a set of Amber parameter files """
        params = parameters.AmberParameterSet(
                os.path.join(get_fn('parm'), 'parm99.dat'),
                os.path.join(get_fn('parm'), 'frcmod.ff99SB'),
                os.path.join(get_fn('parm'), 'frcmod.parmbsc0'),
        )
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

class TestCoordinateFiles(unittest.TestCase):
    """ Tests the various coordinate file classes """
    
    def testMdcrd(self):
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

    def testRestart(self):
        """ Test the ASCII restart file parsing """
        restart = asciicrd.AmberAsciiRestart(get_fn('tz2.ortho.rst7'), 'r')
        self.assertEqual(restart.natom, 5293)
        self.assertTrue(restart.hasbox)
        self.assertFalse(restart.hasvels)
        self.assertIsInstance(restart.coordinates, np.ndarray)
        crdsum = restart.coordinates.sum()
        self.assertAlmostEqual(crdsum, 301623.26028240257, places=4)

class TestAmberMask(unittest.TestCase):
    """ Test the Amber mask parser """
    
    def testMask(self):
        """ Test the Amber mask parser """
        parm = readparm.AmberParm(get_fn('trx.prmtop'))
        mask_res1 = mask.AmberMask(parm, ':1')
        mask_resala = mask.AmberMask(parm, ':ALA')
        mask_atname = mask.AmberMask(parm, '@CA')
        mask_resat = mask.AmberMask(parm, ':ALA@CA')
        mask_attyp = mask.AmberMask(parm, '@%CT')

        # Check all of the masks
        self.assertEqual(sum(mask_res1.Selection()), 13)
        for idx in mask_res1.Selected():
            self.assertEqual(parm.atoms[idx].residue.idx, 0)
        self.assertEqual(list(range(13)), list(mask_res1.Selected()))
        sel = mask_res1.Selection()
        for atom in parm.atoms:
            if atom.residue.idx == 0:
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)

        self.assertEqual(sum(mask_resala.Selection()), 121)
        for idx in mask_resala.Selected():
            self.assertEqual(parm.atoms[idx].residue.name, 'ALA')
        sel = mask_resala.Selection()
        for atom in parm.atoms:
            if atom.residue.name == 'ALA':
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)
        
        self.assertEqual(sum(mask_atname.Selection()), 108)
        for idx in mask_atname.Selected():
            self.assertEqual(parm.atoms[idx].name, 'CA')
        sel = mask_atname.Selection()
        for atom in parm.atoms:
            if atom.name == 'CA':
                self.assertEqual(sel[atom.idx], 1)
            else:
                self.assertEqual(sel[atom.idx], 0)

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

    def testCompoundMask(self):
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

    def testDistanceBasedMaskPDB(self):
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

    def testDistanceBasedMask(self):
        """ Test distance-based mask selections """
        parm = readparm.AmberParm(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        # All atoms within 5 A of residue 8
        mask1 = mask.AmberMask(parm, ':8<@5')
        # All atoms more than 10 A away from residue 1
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

class TestWriteFiles(FileIOTestCase):
    
    def testWriteAmberParm(self):
        """ Test writing an AmberParm file """
        parm = readparm.AmberParm(get_fn('trx.prmtop'))
        parm.write_parm(get_fn('trx.prmtop', written=True))
        f1 = open(get_fn('trx.prmtop'), 'r')
        f2 = open(get_fn('trx.prmtop', written=True), 'r')
        try:
            for line1, line2 in zip(f1, f2):
                if line1.startswith('%VERSION'):
                    self.assertTrue(line2.startswith('%VERSION'))
                    continue
                self.assertEqual(line1.strip(), line2.strip())
        finally:
            f1.close()
            f2.close()

    def testSaveAmberParm(self):
        """ Test writing AmberParm file with AmberParm.save """
        parm = readparm.AmberParm(get_fn('trx.prmtop'))
        parm.add_flag('NEW_FLAG', '10I6', num_items=parm.ptr('nres'))
        self.assertIn('NEW_FLAG', parm.parm_data)
        self.assertIn('NEW_FLAG', parm.flag_list)
        parm.save(get_fn('trx.prmtop', written=True))
        parm2 = readparm.AmberParm(get_fn('trx.prmtop', written=True))
        self.assertIn('NEW_FLAG', parm2.parm_data)

    def testAmberRestart(self):
        """ Test writing an ASCII Amber restart file """
        Restart = asciicrd.AmberAsciiRestart
        box = [10, 10, 10, 90, 90, 90]
        rst = Restart(get_fn('testc.rst7', written=True), 'w', natom=9)
        rst.coordinates = list(range(27))
        rst.close()
        rst = Restart(get_fn('testcv.rst7', written=True), 'w', natom=20)
        rst.coordinates = list(range(60))
        rst.velocities = list(reversed(range(60)))
        rst.close()
        rst = Restart(get_fn('testcb.rst7', written=True), 'w', natom=7)
        rst.coordinates = list(range(21))
        rst.box = box[:]
        rst.close()
        rst = Restart(get_fn('testcvb.rst7', written=True), 'w', natom=15)
        rst.coordinates = list(range(45))
        rst.velocities = list(reversed(range(45)))
        rst.box = box[:]
        rst.close()
        self._check_written_restarts(box)

    def testAmberRestartNumpy(self):
        """ Test writing Amber restart file passing numpy arrays """
        Restart = asciicrd.AmberAsciiRestart
        box = np.asarray([10, 10, 10, 90, 90, 90])
        rst = Restart(get_fn('testc.rst7', written=True), 'w', natom=9)
        rst.coordinates = np.arange(27).reshape((9,3))
        rst.close()
        rst = Restart(get_fn('testcv.rst7', written=True), 'w', natom=20)
        rst.coordinates = np.arange(60).reshape((20,3))
        rst.velocities = np.asarray(list(reversed(range(60)))).reshape((20,3))
        rst.close()
        rst = Restart(get_fn('testcb.rst7', written=True), 'w', natom=7)
        rst.coordinates = np.arange(21).reshape((7,3))
        rst.box = box
        rst.close()
        rst = Restart(get_fn('testcvb.rst7', written=True), 'w', natom=15)
        rst.coordinates = np.arange(45).reshape((15,3))
        rst.velocities = np.asarray(list(reversed(range(45)))).reshape((15,3))
        rst.box = box
        rst.close()
        self._check_written_restarts(box)

    def testAmberMdcrd(self):
        """ Test writing ASCII trajectory file """
        box = [15, 15, 15]
        Mdcrd = asciicrd.AmberMdcrd
        crd = Mdcrd(get_fn('testc.mdcrd', written=True), natom=15, hasbox=False,
                    mode='w', title='Test file')
        crd.add_coordinates(list(range(45)))
        crd.add_coordinates([x+1 for x in range(45)])
        crd.add_coordinates([x+2 for x in range(45)])
        crd.add_coordinates([x+3 for x in range(45)])
        crd.add_coordinates([x+4 for x in range(45)])
        crd.close()
        crd = Mdcrd(get_fn('testcb.mdcrd', written=True), natom=18, hasbox=True,
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

    def testAmberMdcrdNumpy(self):
        """ Test writing ASCII trajectory file passing numpy arrays """
        box = np.asarray([15, 15, 15])
        Mdcrd = asciicrd.AmberMdcrd
        crd = Mdcrd(get_fn('testc.mdcrd', written=True), natom=15, hasbox=False,
                    mode='w', title='Test file')
        coorddata = np.arange(45).reshape((15,3))
        crd.add_coordinates(coorddata)
        crd.add_coordinates(coorddata+1)
        crd.add_coordinates(coorddata+2)
        crd.add_coordinates(coorddata+3)
        crd.add_coordinates(coorddata+4)
        crd.close()
        crd = Mdcrd(get_fn('testcb.mdcrd', written=True), natom=18, hasbox=True,
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

    def testBadFileUsage(self):
        """ Check that illegal file usage results in desired exceptions """
        Restart = asciicrd.AmberAsciiRestart
        Mdcrd = asciicrd.AmberMdcrd
        box = [10, 10, 10, 90, 90, 90]
        rst = Restart(get_fn('testc.rst7', written=True), 'w', natom=9,
                      hasbox=True)
        def assign(obj, stmnt):
            rst = crd = obj
            exec(stmnt)
        try:
            self.assertRaises(ValueError,
                              lambda: assign(rst, 'rst.coordinates=range(20)'))
            self.assertRaises(RuntimeError,
                              lambda: assign(rst, 'rst.box=[10]*3+[90]*3'))
            rst.coordinates = list(range(27))
            rst.box = box
            self.assertRaises(RuntimeError, lambda:
                              assign(rst, 'rst.velocities=list(range(27))'))
        finally:
            rst.close()
        crd = Mdcrd(get_fn('testc.mdcrd', written=True), natom=15, hasbox=True,
                    mode='w', title='Test file')
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
        rst = readparm.Rst7.open(get_fn('testc.rst7', written=True))
        self.assertFalse(rst.hasbox)
        self.assertFalse(rst.hasvels)
        np.testing.assert_equal(rst.coordinates.flatten(), list(range(27)))
        rst = readparm.Rst7.open(get_fn('testcb.rst7', written=True))
        self.assertTrue(rst.hasbox)
        self.assertFalse(rst.hasvels)
        np.testing.assert_equal(rst.coordinates.flatten(), list(range(21)))
        np.testing.assert_equal(rst.box.flatten(), box)
        rst = readparm.Rst7.open(get_fn('testcv.rst7', written=True))
        self.assertTrue(rst.hasvels)
        self.assertFalse(rst.hasbox)
        np.testing.assert_equal(rst.coordinates,
                np.arange(60).reshape(rst.coordinates.shape))
        np.testing.assert_allclose(rst.velocities,
                np.array(list(reversed(range(60)))).reshape(rst.velocities.shape))
        rst = readparm.Rst7.open(get_fn('testcvb.rst7', written=True))
        self.assertTrue(rst.hasvels)
        self.assertTrue(rst.hasbox)
        np.testing.assert_equal(rst.coordinates,
                np.arange(45).reshape(rst.coordinates.shape))
        np.testing.assert_allclose(rst.velocities,
                np.array(list(reversed(range(45)))).reshape(rst.velocities.shape))
        np.testing.assert_equal(rst.box, box)

    def _check_written_mdcrds(self, box):
        # Now try to read them and verify the information
        crd = asciicrd.AmberMdcrd(get_fn('testc.mdcrd', written=True),
                                  15, False, 'r')
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

        crd = asciicrd.AmberMdcrd(get_fn('testcb.mdcrd', written=True),
                                  18, True, 'r')
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

    def testTrackedList(self):
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

    def testSplit(self):
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

    def testSplit2(self):
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

    def testAdd(self):
        """ Tests combining AmberParm instances """
        parm1 = readparm.AmberParm(get_fn('phenol.prmtop'))
        parm2 = readparm.AmberParm(get_fn('biphenyl.prmtop'))
        comb = parm1 + parm2
        self.assertEqual(len(comb.atoms), len(parm1.atoms) + len(parm2.atoms))
        for a1, a2 in zip(comb.atoms, parm1.atoms + parm2.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.radii, a2.radii)
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
            self.assertEqual(a1.radii, a2.radii)
        self.assertEqual(len(parm1.residues), len(comb.residues))
        for r1, r2 in zip(comb.residues, parm1.residues):
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.chain, r2.chain)

    def testMult(self):
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
            self.assertEqual(a1.radii, a2.radii)
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
            self.assertEqual(a1.radii, a2.radii)
        self.assertEqual(len(parm.residues), len(mult.residues))
        for r1, r2 in zip(mult.residues, parm.residues):
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.chain, r2.chain)

    def testSimpleSlice(self):
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
            self.assertEqual(a1.radii, a2.radii)
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

if __name__ == '__main__':
    unittest.main()
