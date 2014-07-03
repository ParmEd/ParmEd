"""
Tests the functionality in the chemistry.amber package
"""

from array import array
from chemistry.amber import readparm, asciicrd, mask
import unittest
from utils import get_fn, has_numpy

class TestReadParm(unittest.TestCase):
    
    def testLoadParm(self):
        parm = readparm.LoadParm(get_fn('trx.prmtop'))
        parm2 = readparm.AmberParm(get_fn('trx.prmtop'))
        for key in parm.parm_data:
            self.assertEqual(parm.parm_data[key], parm2.parm_data[key])

    def testAmberGasParm(self):
        parm = readparm.AmberParm(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        gasparm = readparm.AmberParm(get_fn('trx.prmtop'))
        gasparm.LoadRst7(get_fn('trx.inpcrd'))

        self.assertEqual([a.xx for a in gasparm.atom_list],
                         [a.xx for a in parm.atom_list])
        self.assertEqual([a.xy for a in gasparm.atom_list],
                         [a.xy for a in parm.atom_list])
        self.assertEqual([a.xz for a in gasparm.atom_list],
                         [a.xz for a in parm.atom_list])
        
        # Now run the tests for the prmtop
        self._standard_parm_tests(parm)
        self.assertFalse(parm.chamber)
        self.assertFalse(parm.amoeba)
        self.assertRaises(KeyError, lambda: parm.parm_data['BOX_DIMENSIONS'])
        self.assertEqual(parm.ptr('ifbox'), 0)

        # Now check the restart file
        rst = readparm.Rst7(get_fn('trx.inpcrd'))
        coords = rst.coordinates
        vels = rst.velocities
        for i, atom in enumerate(gasparm.atom_list):
            i3 = i * 3
            self.assertEqual(coords[i3  ], atom.xx)
            self.assertEqual(coords[i3+1], atom.xy)
            self.assertEqual(coords[i3+2], atom.xz)
            self.assertEqual(vels[i3  ], atom.vx)
            self.assertEqual(vels[i3+1], atom.vy)
            self.assertEqual(vels[i3+2], atom.vz)

    def testAmberSolvParm(self):
        parm = readparm.AmberParm(get_fn('solv.prmtop'),
                                  get_fn('solv.rst7'))
        self._standard_parm_tests(parm)
        self._solv_pointer_tests(parm)
        self.assertFalse(parm.chamber)
        self.assertFalse(parm.amoeba)
        self.assertEqual(parm.ptr('ifbox'), 1)

    def testChamberGasParm(self):
        parm = readparm.ChamberParm(get_fn('ala_ala_ala.parm7'))
        self._standard_parm_tests(parm)
        self.assertTrue(parm.chamber)
        self.assertTrue(parm.has_cmap)
        self.assertEqual(parm.ptr('ifbox'), 0)

    def testChamberSolvParm(self):
        parm = readparm.ChamberParm(get_fn('dhfr_cmap_pbc.parm7'))
        self._standard_parm_tests(parm)
        self._solv_pointer_tests(parm)
        self.assertTrue(parm.chamber)
        self.assertTrue(parm.has_cmap)
        self.assertEqual(parm.ptr('ifbox'), 1)

    def testAmoebaBig(self):
        parm = readparm.AmoebaParm(get_fn('amoeba.parm7'))
        self.assertEqual(parm.ptr('natom'), len(parm.atom_list))
        self.assertEqual([a.atname for a in parm.atom_list],
                         parm.parm_data['ATOM_NAME'])
        self.assertEqual([a.attype for a in parm.atom_list],
                         parm.parm_data['AMBER_ATOM_TYPE'])
        self.assertTrue(parm.amoeba)
        for attr in ['bond_list', 'angle_list', 'urey_bradley_list',
                     'angle_list', 'trigonal_angle_list', 'oopbend_list',
                     'dihedral_list', 'pitorsion_list', 'stretch_bend_list',
                     'torsion_torsion_list', 'chiral_frame_list',
                     'multipole_frame_list', 'adjust_list', 'adjust_weights',
                     'bond_type_list', 'angle_type_list',
                     'trigonal_angle_type_list', 'oopbend_type_list',
                     'dihedral_type_list', 'pitorsion_type_list',
                     'stretch_bend_type_list', 'torsion_torsion_type_list']:
            self.assertTrue(hasattr(parm, attr))

    def testAmoebaSmall(self):
        parm = readparm.AmoebaParm(get_fn('nma.parm7'))
        rst7 = readparm.BeemanRestart(get_fn('nma.rst7'))
        self.assertEqual(3*rst7.natom, len(rst7.coordinates))
        self.assertEqual(rst7.coordinates, rst7.parm_data['ATOMIC_COORDS_LIST'])
        self.assertEqual(rst7.natom, parm.ptr('natom'))
        self.assertFalse(parm.torsion_torsion_list)
        self.assertTrue(parm.amoeba)

    # Tests for individual prmtops
    def _standard_parm_tests(self, parm):
        self.assertEqual(parm.ptr('natom'), len(parm.atom_list))
        self.assertEqual(parm.ptr('nres'), len(parm.residue_list))
        self.assertEqual(parm.ptr('nbonh'), len(parm.bonds_inc_h))
        self.assertEqual(parm.ptr('nbona'), len(parm.bonds_without_h))
        self.assertEqual(parm.ptr('ntheth'), len(parm.angles_inc_h))
        self.assertEqual(parm.ptr('ntheta'), len(parm.angles_without_h))
        self.assertEqual(parm.ptr('nphih'), len(parm.dihedrals_inc_h))
        self.assertEqual(parm.ptr('nphia'), len(parm.dihedrals_without_h))
        self.assertEqual([a.atname for a in parm.atom_list],
                         parm.parm_data['ATOM_NAME'])
        self.assertEqual([a.attype for a in parm.atom_list],
                         parm.parm_data['AMBER_ATOM_TYPE'])

    def _solv_pointer_tests(self, parm):
        self.assertEqual(parm.ptr('nspm'),
                         parm.parm_data['SOLVENT_POINTERS'][1])
        self.assertEqual(parm.ptr('nspm'),
                         len(parm.parm_data['ATOMS_PER_MOLECULE']))
        self.assertEqual(parm.ptr('natom'),
                         sum(parm.parm_data['ATOMS_PER_MOLECULE']))

class TestCoordinateFiles(unittest.TestCase):
    
    def testMdcrd(self):
        mdcrd = asciicrd.AmberMdcrd(get_fn('tz2.truncoct.crd'),
                                    natom=5827, hasbox=True, mode='r')
        self.assertEqual(mdcrd.frame, 10)
        if has_numpy():
            import numpy as np
            self.assertIsInstance(mdcrd.coordinates(0), np.ndarray)
            # Hack numpy off so we test that code path, too
            asciicrd.np = None
            asciicrd.array = array
            mdcrdarray = asciicrd.AmberMdcrd(get_fn('tz2.truncoct.crd'),
                                             natom=5827, hasbox=True, mode='r')
            asciicrd.np = np
        else:
            mdcrdarray = mdcrd
        self.assertIsInstance(mdcrdarray.coordinates(0), array)

        if has_numpy():
            runsum = 0
            for i in range(10):
                arr1 = mdcrd.coordinates(i)
                arr2 = mdcrdarray.coordinates(i)
                runsum += arr1.sum()
                self.assertAlmostEqual(arr1.sum(), sum(arr2), places=3)
        else:
            runsum = 0
            for i in range(10):
                runsum += sum(mdcrdarray.coordinates(i))
        self.assertAlmostEqual(runsum, 7049.817, places=3)

    def testRestart(self):
        restart = asciicrd.AmberAsciiRestart(get_fn('tz2.ortho.rst7'), 'r')
        self.assertEqual(restart.natom, 5293)
        self.assertTrue(restart.hasbox)
        self.assertFalse(restart.hasvels)
        if has_numpy():
            import numpy as np
            # Hack numpy off so we test that code path, too
            asciicrd.np = None
            asciicrd.array = array
            restartarray = asciicrd.AmberAsciiRestart(
                                        get_fn('tz2.ortho.rst7'), 'r')
            self.assertAlmostEqual(restart.coordinates.sum(),
                                   sum(restartarray.coordinates), places=4)
            crdsum = restart.coordinates.sum()
        else:
            restartarray = restart
            crdsum = sum(restart.coordinates)
        self.assertIsInstance(restart.coordinates, np.ndarray)
        self.assertAlmostEqual(crdsum, 301623.26028240257, places=4)

class TestAmberMask(unittest.TestCase):
    
    def testMask(self):
        parm = readparm.AmberParm(get_fn('trx.prmtop'))
        mask_res1 = mask.AmberMask(parm, ':1')
        mask_resala = mask.AmberMask(parm, ':ALA')
        mask_atname = mask.AmberMask(parm, '@CA')
        mask_resat = mask.AmberMask(parm, ':ALA@CA')
        mask_attyp = mask.AmberMask(parm, '@%CT')

        # Check all of the masks
        self.assertEqual(sum(mask_res1.Selection()), 13)
        for idx in mask_res1.Selected():
            self.assertEqual(parm.atom_list[idx].residue.idx, 1)
        self.assertEqual(range(13), list(mask_res1.Selected()))

        self.assertEqual(sum(mask_resala.Selection()), 121)
        for idx in mask_resala.Selected():
            self.assertEqual(parm.atom_list[idx].residue.resname, 'ALA')
        
        self.assertEqual(sum(mask_atname.Selection()), 108)
        for idx in mask_atname.Selected():
            self.assertEqual(parm.atom_list[idx].atname, 'CA')

        self.assertEqual(sum(mask_resat.Selection()), 12)
        for idx in mask_resat.Selected():
            self.assertEqual(parm.atom_list[idx].atname, 'CA')
            self.assertEqual(parm.atom_list[idx].residue.resname, 'ALA')

        self.assertEqual(sum(mask_attyp.Selection()), 341)
        for idx in mask_attyp.Selected():
            self.assertEqual(parm.atom_list[idx].attype, 'CT')
