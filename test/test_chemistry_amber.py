"""
Tests the functionality in the chemistry.amber package
"""

from array import array
from compat24 import all
from chemistry.amber import readparm, asciicrd, mask
from chemistry import topologyobjects
import os
import unittest
from utils import get_fn, has_numpy

class TestReadParm(unittest.TestCase):
    """ Tests the various Parm file classes """
    
    def testLoadParm(self):
        """ Test the arbitrary parm loader """
        parm = readparm.LoadParm(get_fn('trx.prmtop'))
        parm2 = readparm.AmberParm(get_fn('trx.prmtop'))
        for key in parm.parm_data:
            self.assertEqual(parm.parm_data[key], parm2.parm_data[key])

    def testAmberGasParm(self):
        """ Test the AmberParm class with a non-periodic (gas-phase) prmtop """
        parm = readparm.AmberParm(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        gasparm = readparm.AmberParm(get_fn('trx.prmtop'))
        gasparm.load_rst7(get_fn('trx.inpcrd'))

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
            i3 = i * 3
            self.assertEqual(coords[i3  ], atom.xx)
            self.assertEqual(coords[i3+1], atom.xy)
            self.assertEqual(coords[i3+2], atom.xz)
            self.assertEqual(vels[i3  ], atom.vx)
            self.assertEqual(vels[i3+1], atom.vy)
            self.assertEqual(vels[i3+2], atom.vz)

    def testRemakeParm(self):
        """ Tests the rebuilding of the AmberParm raw data structures """
        parm = readparm.AmberParm(get_fn('trx.prmtop'))
        parm2 = readparm.AmberParm(get_fn('trx.prmtop'))
        parm.remake_parm()
        self.assertEqual(parm.flag_list, parm2.flag_list)
        for flag in parm.flag_list:
            self.assertEqual(parm.parm_data[flag], parm2.parm_data[flag])

    def testRemakeChamberParm(self):
        """ Tests the rebuilding of the ChamberParm raw data structures """
        parm = readparm.ChamberParm(get_fn('ala_ala_ala.parm7'))
        parm2 = readparm.ChamberParm(get_fn('ala_ala_ala.parm7'))
        parm.remake_parm()
        self.assertEqual(set(parm.flag_list), set(parm2.flag_list))
        for flag in parm.flag_list:
            self.assertEqual(parm.parm_data[flag], parm2.parm_data[flag])

    def testAmberSolvParm(self):
        """ Test the AmberParm class with a periodic prmtop """
        parm = readparm.AmberParm(get_fn('solv.prmtop'),
                                  get_fn('solv.rst7'))
        self._standard_parm_tests(parm)
        self._solv_pointer_tests(parm)
        self.assertFalse(parm.chamber)
        self.assertFalse(parm.amoeba)
        self.assertEqual(parm.ptr('ifbox'), 1)

    def testChamberGasParm(self):
        """Test the ChamberParm class with a non-periodic (gas phase) prmtop"""
        parm = readparm.ChamberParm(get_fn('ala_ala_ala.parm7'))
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertTrue(parm.chamber)
        self.assertTrue(parm.has_cmap)
        self.assertEqual(parm.ptr('ifbox'), 0)

    def testChamberSolvParm(self):
        """ Test the ChamberParm class with a periodic prmtop """
        parm = readparm.ChamberParm(get_fn('dhfr_cmap_pbc.parm7'))
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

    # Tests for individual prmtops
    def _standard_parm_tests(self, parm):
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

class TestCoordinateFiles(unittest.TestCase):
    """ Tests the various coordinate file classes """
    
    def testMdcrd(self):
        """ Test the ASCII trajectory file parsing """
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
        """ Test the ASCII restart file parsing """
        restart = asciicrd.AmberAsciiRestart(get_fn('tz2.ortho.rst7'), 'r')
        self.assertEqual(restart.natom, 5293)
        self.assertTrue(restart.hasbox)
        self.assertFalse(restart.hasvels)
        if has_numpy():
            import numpy as np
            self.assertIsInstance(restart.coordinates, np.ndarray)
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
        self.assertIsInstance(restartarray.coordinates, array)
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

class TestWriteFiles(unittest.TestCase):
    
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
        import numpy as np
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
        import numpy as np
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
        for x1, x2 in zip(rst.coordinates, range(27)):
            self.assertEqual(x1, x2)
        rst = readparm.Rst7.open(get_fn('testcb.rst7', written=True))
        self.assertTrue(rst.hasbox)
        self.assertFalse(rst.hasvels)
        for x1, x2 in zip(rst.coordinates, range(21)):
            self.assertEqual(x1, x2)
        for x1, x2 in zip(rst.box, box):
            self.assertEqual(x1, x2)
        rst = readparm.Rst7.open(get_fn('testcv.rst7', written=True))
        self.assertTrue(rst.hasvels)
        self.assertFalse(rst.hasbox)
        for x1, x2 in zip(rst.coordinates, range(60)):
            self.assertEqual(x1, x2)
        for x1, x2 in zip(rst.vels, reversed(range(60))):
            self.assertAlmostEqual(x1, x2, places=5)
        rst = readparm.Rst7.open(get_fn('testcvb.rst7', written=True))
        self.assertTrue(rst.hasvels)
        self.assertTrue(rst.hasbox)
        for x1, x2 in zip(rst.coordinates, range(45)):
            self.assertEqual(x1, x2)
        for x1, x2 in zip(rst.vels, reversed(range(45))):
            self.assertAlmostEqual(x1, x2, places=5)
        for x1, x2 in zip(rst.box, box):
            self.assertEqual(x1, x2)

    def _check_written_mdcrds(self, box):
        # Now try to read them and verify the information
        crd = asciicrd.AmberMdcrd(get_fn('testc.mdcrd', written=True),
                                  15, False, 'r')
        self.assertEqual(crd.title, 'Test file')
        self.assertFalse(crd.hasbox)
        for i in range(crd.frame):
            for x1, x2 in zip(crd.coordinates(i), [x+i for x in range(45)]):
                self.assertEqual(x1, x2)
        for i, array in enumerate(crd.coordinates()):
            for x1, x2 in zip(array, [x+i for x in range(45)]):
                self.assertEqual(x1, x2)
        crd.close()

        crd = asciicrd.AmberMdcrd(get_fn('testcb.mdcrd', written=True),
                                  18, True, 'r')
        self.assertEqual(crd.title, 'Test file')
        self.assertTrue(crd.hasbox)
        for i in range(crd.frame):
            for x1, x2 in zip(crd.coordinates(i), [x+i for x in range(54)]):
                self.assertEqual(x1, x2)
            for x1, x2 in zip(crd.box(i), box):
                self.assertEqual(x1, x2)
        for i, (coords, mybox) in enumerate(zip(crd.coordinates(), crd.box())):
            for x1, x2 in zip(coords, [x+i for x in range(45)]):
                self.assertEqual(x1, x2)
            for x1, x2 in zip(mybox, box):
                self.assertEqual(x1, x2)

class TestObjectAPIs(unittest.TestCase):
    """ Tests various object APIs """

    def testTrackedList(self):
        """ Tests the TrackedList object """
        mylist = topologyobjects.TrackedList(range(20))
        mylist2 = topologyobjects.TrackedList(reversed(range(20)))
        self.assertFalse(mylist.changed)
        self.assertIsInstance(mylist[0], int)
        self.assertIsInstance(mylist[0:5], topologyobjects.TrackedList)
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

if not has_numpy():
    del TestWriteFiles.testAmberRestartNumpy, TestWriteFiles.testAmberMdcrdNumpy

if __name__ == '__main__':
    unittest.main()
