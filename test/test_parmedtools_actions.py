"""
Tests for the various actions in ParmEd
"""
from __future__ import division, print_function

import utils
from parmed import periodic_table, gromacs, load_file
from parmed.amber import AmberParm, ChamberParm, AmoebaParm
from parmed.charmm import CharmmPsfFile
from parmed.exceptions import AmberWarning, CharmmWarning
from parmed.formats import PDBFile, CIFFile
from parmed.utils.six.moves import range, zip
from copy import copy
import os
import parmed.tools as PT
from parmed.tools import exceptions as exc
from parmed.tools import parmlist
import saved_outputs as saved
import sys
import unittest
import warnings

get_fn = utils.get_fn
get_saved_fn = utils.get_saved_fn
diff_files = utils.diff_files

gasparm = AmberParm(get_fn('trx.prmtop'))
solvparm = AmberParm(get_fn('solv.prmtop'))
gascham = ChamberParm(get_fn('ala_ala_ala.parm7'))
solvchamber = ChamberParm(get_fn('dhfr_cmap_pbc.parm7'))
amoebaparm = AmoebaParm(get_fn('nma.parm7'))

class TestNonParmActions(unittest.TestCase):
    """ Tests all actions that do not require a prmtop instance """

    def setUp(self):
        self.parm = gasparm

    def testOverwrite(self):
        """ Test setting overwrite capabilities on ParmEd interpeter """
        self.assertTrue(PT.Action.overwrite)
        a = PT.setOverwrite(self.parm, False)
        self.assertTrue(PT.Action.overwrite)
        a.execute()
        self.assertFalse(PT.Action.overwrite)
        self.assertEqual(str(a), 'Files are NOT overwritable')
        a = PT.setOverwrite(self.parm, True)
        a.execute()
        self.assertTrue(PT.Action.overwrite)
        self.assertEqual(str(a), 'Files are overwritable')
    
    def testListParms(self):
        """ Test listing of the prmtop files in the ParmEd interpreter """
        a = PT.listParms(self.parm)
        a.execute() # Should do nothing
        lines = str(a).split('\n')
        self.assertEqual(lines[0], 'Loaded topology files:')
        self.assertEqual(lines[1], '[0]\t%s (active)' % get_fn('trx.prmtop'))

    def testChamber(self):
        """ Test the chamber action with a basic protein """
        # To keep stderr clean
        warnings.filterwarnings('ignore', category=CharmmWarning,
                                module='psf')
        a = PT.chamber(self.parm, '-psf %s' % get_fn('ala_ala_ala.psf'),
                       '-top %s' % get_fn('top_all22_prot.inp'),
                       '-param %s' % get_fn('par_all22_prot.inp'),
                       '-crd %s' % get_fn('ala_ala_ala.pdb'))
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertTrue(parm.chamber)
        self.assertTrue(parm.has_cmap)
        self.assertEqual(parm.ptr('ifbox'), 0)

    def testChamberModel(self):
        """ Test the chamber action with a model compound """
        # To keep stderr clean
        warnings.filterwarnings('ignore', category=CharmmWarning,
                                module='psf')
        a = PT.chamber(self.parm, '-psf %s' % get_fn('propane.psf'),
                       '-top %s' % get_fn('top_all36_prot.rtf'),
                       '-param %s' % get_fn('par_all36_prot.prm'),
                       '-str %s' % get_fn('toppar_all36_prot_model.str'),
                       '-str %s' % get_fn('toppar_water_ions.str'),
                       '-crd %s' % get_fn('propane.pdb'))
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertTrue(parm.chamber)
        self.assertEqual(parm.ptr('ifbox'), 0)

    def testChamberGlobbing(self):
        """ Test globbing in the chamber action """
        warnings.filterwarnings('ignore', category=CharmmWarning,
                                module='psf')
        a = PT.chamber(self.parm, '-psf', get_fn('ala_ala_ala.psf'),
                       '-toppar', get_fn('*_all22_prot.inp'),
                       '-crd', get_fn('ala_ala_ala.pdb'))
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertTrue(parm.chamber)
        self.assertTrue(parm.has_cmap)
        self.assertEqual(parm.ptr('ifbox'), 0)

    def testChamberNbfix(self):
        """ Test the chamber action with a complex system using NBFIX """
        warnings.filterwarnings('ignore', category=CharmmWarning,
                                module='psf')
        a = PT.chamber(self.parm, '-psf %s' % get_fn('ala3_solv.psf'),
                       '-param %s' % get_fn('par_all36_prot.prm'),
                       '-str %s' % get_fn('toppar_water_ions.str'),
                       '-crd %s' % get_fn('ala3_solv.crd'), '-box bounding')
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)

    def testChamberBug1(self):
        """ Test chamber BFNA creation (former bug) """
        a = PT.chamber(self.parm, '-top', get_fn('top_all36_cgenff.rtf'),
                '-top', get_fn('top_bfna_nonbonded_stitched.rtf'), '-param',
                get_fn('par_bfna_nonbonded_stitched.prm'), '-param',
                get_fn('par_all36_cgenff.prm'), '-box', 'bounding', '-psf',
                get_fn('bfna_nonbonded_vmd_autopsf.psf'), 'nocondense', '-crd',
                get_fn('bfna_nonbonded_vmd_autopsf.pdb'), '-nocmap'
        )
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        psf = CharmmPsfFile(get_fn('bfna_nonbonded_vmd_autopsf.psf'))
        self.assertEqual(len(psf.atoms), len(parm.atoms))
        self.assertEqual(len(psf.residues), len(parm.residues))
        for a1, a2 in zip(psf.atoms, parm.atoms):
            self.assertEqual(a1.name[:4], a2.name[:4])
            self.assertEqual(a1.type[:4], a2.type[:4])
            self.assertAlmostEqual(a1.charge, a2.charge)
            self.assertAlmostEqual(a1.mass, a2.mass)

    def testChamberBug2(self):
        """ Test that chamber sets the box angles for triclinics correctly """
        warnings.filterwarnings('ignore', category=CharmmWarning,
                                module='psf')
        a = PT.chamber(self.parm, '-psf %s' % get_fn('ala3_solv.psf'),
                       '-param %s' % get_fn('par_all36_prot.prm'),
                       '-str %s' % get_fn('toppar_water_ions.str'),
                       '-crd %s' % get_fn('ala3_solv.crd'), '-box',
                       '33,33,33,109.475,109.475,109.475')
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        for x, y in zip(parm.parm_data['BOX_DIMENSIONS'], [109.475] + [33]*3):
            self.assertAlmostEqual(x, y)
        for x, y in zip(parm.box, [33]*3 + [109.475]*3):
            self.assertAlmostEqual(x, y)

    def testGromber(self):
        """ Test the gromber action on a small system (no coords) """
        a = PT.gromber(None, os.path.join(get_fn('03.AlaGlu'), 'topol.top'))
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertIs(parm.box, None)

    def testGromber2(self):
        """ Test the gromber action with coordinates """
        a = PT.gromber(None, os.path.join(get_fn('03.AlaGlu'), 'topol.top'),
                       os.path.join(get_fn('03.AlaGlu'), 'conf.gro'))
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertIs(parm.box, None) # AmberParm deletes the box without solvent
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))

    def testGromber3(self):
        """ Test the gromber action passing various defines """
        a = PT.gromber(None, os.path.join(get_fn('03.AlaGlu'), 'topol.top'),
                       os.path.join(get_fn('03.AlaGlu'), 'conf.gro'),
                       'define SOMEDEF=this define SOMEDEF2=that')
        stra = str(a)
        self.assertIn('SOMEDEF', stra)
        self.assertIn('SOMEDEF2', stra)
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertIs(parm.box, None)
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))

    def testGromberBox(self):
        """ Test the gromber action when a box should be defined """
        a = PT.gromber(None, os.path.join(get_fn('09.DHFR-PME'), 'topol.top'),
                       os.path.join(get_fn('09.DHFR-PME'), 'conf.gro'))
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
#       self._extensive_checks(parm)
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))
        self.assertEqual(parm.box[0], 62.23)
        self.assertEqual(parm.box[1], 62.23)
        self.assertEqual(parm.box[2], 62.23)
        self.assertEqual(parm.box[3], 90.00)
        self.assertEqual(parm.box[4], 90.00)
        self.assertEqual(parm.box[5], 90.00)

    # Copied from test_parmed_amber -- tests the prmtop file generated by the
    # "chamber" action
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
        self.assertEqual([a.type[:4] for a in parm.atoms],
                         parm.parm_data['AMBER_ATOM_TYPE'])
    
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

class TestAmberParmActions(utils.FileIOTestCase, utils.TestCaseRelative):
    """ Tests actions on Amber prmtop files """
    
    def testParmoutOutparmLoadRestrt(self):
        """ Test parmout, outparm, and loadRestrt actions on AmberParm """
        self._empty_writes()
        parm = copy(gasparm)
        PT.loadRestrt(parm, get_fn('trx.inpcrd')).execute()
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))
        PT.parmout(parm, get_fn('test.parm7', written=True)).execute()
        self.assertEqual(len(os.listdir(get_fn('writes'))), 1)
        self.assertTrue(diff_files(get_fn('trx.prmtop'),
                                   get_fn('test.parm7', written=True)))
        self._empty_writes()
        PT.parmout(parm, get_fn('test.parm7', written=True),
                         get_fn('test.rst7', written=True)).execute()
        self.assertEqual(len(os.listdir(get_fn('writes'))), 2)
        self.assertTrue(diff_files(get_fn('trx.prmtop'),
                                   get_fn('test.parm7', written=True)))
        self.assertTrue(diff_files(get_fn('trx.inpcrd'),
                                   get_fn('test.rst7', written=True),
                                   absolute_error=0.0001))
        self._empty_writes()
        PT.outparm(parm, get_fn('test.parm7', written=True)).execute()
        self.assertEqual(len(os.listdir(get_fn('writes'))), 1)
        self.assertTrue(diff_files(get_fn('trx.prmtop'),
                                   get_fn('test.parm7', written=True)))
        self._empty_writes()
        PT.outparm(parm, get_fn('test.parm7', written=True),
                         get_fn('test.rst7', written=True)).execute()
        self.assertEqual(len(os.listdir(get_fn('writes'))), 2)
        self.assertTrue(diff_files(get_fn('trx.prmtop'),
                                   get_fn('test.parm7', written=True)))
        self.assertTrue(diff_files(get_fn('trx.inpcrd'),
                                   get_fn('test.rst7', written=True),
                                   absolute_error=0.0001))

    def testWriteFrcmod(self):
        """ Test writeFrcmod on AmberParm """
        parm = gasparm
        PT.writeFrcmod(parm, get_fn('test.frcmod', written=True)).execute()
        self.assertTrue(diff_files(get_saved_fn('test.frcmod'),
                                   get_fn('test.frcmod', written=True)))

    def testWriteOffLoadRestrt(self):
        """ Test writeOFF on AmberParm """
        parm = copy(gasparm)
        PT.loadRestrt(parm, get_fn('trx.inpcrd')).execute()
        PT.writeOFF(parm, get_fn('test.off', written=True)).execute()
        if utils.has_numpy():
            self.assertTrue(diff_files(get_saved_fn('test.off'),
                                       get_fn('test.off', written=True),
                                       absolute_error=0.0001))
        else:
            self.assertTrue(diff_files(get_saved_fn('test_nonpy.off'),
                                       get_fn('test.off', written=True),
                                       absolute_error=0.0001))

    def testChangeRadii(self):
        """ Test changeRadii on AmberParm """
        parm = copy(gasparm)
        PT.changeRadii(parm, 'amber6').execute()
        self.assertEqual(parm.parm_data['RADIUS_SET'][0],
                         'amber6 modified Bondi radii (amber6)')
        for i, atom in enumerate(parm.atoms):
            radii, atomic_number = atom.radii, atom.atomic_number
            self.assertEqual(parm.parm_data['RADII'][i], radii)
            if atomic_number == 6:
                self.assertEqual(radii, 1.7)
            elif atomic_number == 7:
                self.assertEqual(radii, 1.55)
            elif atomic_number == 8 or atomic_number == 9:
                self.assertEqual(radii, 1.5)
            elif atomic_number == 14:
                self.assertEqual(radii, 2.1)
            elif atomic_number == 15:
                self.assertEqual(radii, 1.85)
            elif atomic_number == 16:
                self.assertEqual(radii, 1.8)
            elif atomic_number == 1:
                if atom.bond_partners[0].atomic_number == 6:
                    self.assertEqual(radii, 1.3)
                elif atom.bond_partners[0].atomic_number in (8, 16):
                    self.assertEqual(radii, 0.8)
                else:
                    self.assertEqual(radii, 1.2)
            else:
                self.assertEqual(radii, 1.5)

        PT.changeRadii(parm, 'bondi').execute()
        self.assertEqual(parm.parm_data['RADIUS_SET'][0], 'Bondi radii (bondi)')
        for atom in parm.atoms:
            radii, atomic_number = atom.radii, atom.atomic_number
            if atomic_number == 6:
                self.assertEqual(radii, 1.7)
            elif atomic_number == 7:
                self.assertEqual(radii, 1.55)
            elif atomic_number == 8 or atomic_number == 9:
                self.assertEqual(radii, 1.5)
            elif atomic_number == 14:
                self.assertEqual(radii, 2.1)
            elif atomic_number == 15:
                self.assertEqual(radii, 1.85)
            elif atomic_number == 16:
                self.assertEqual(radii, 1.8)
            elif atomic_number == 1:
                self.assertEqual(radii, 1.2)
            else:
                self.assertEqual(radii, 1.5)

        PT.changeRadii(parm, 'mbondi').execute()
        self.assertEqual(parm.parm_data['RADIUS_SET'][0],
                         'modified Bondi radii (mbondi)')
        for atom in parm.atoms:
            radii, atomic_number = atom.radii, atom.atomic_number
            if atomic_number == 6:
                self.assertEqual(radii, 1.7)
            elif atomic_number == 7:
                self.assertEqual(radii, 1.55)
            elif atomic_number == 8 or atomic_number == 9:
                self.assertEqual(radii, 1.5)
            elif atomic_number == 14:
                self.assertEqual(radii, 2.1)
            elif atomic_number == 15:
                self.assertEqual(radii, 1.85)
            elif atomic_number == 16:
                self.assertEqual(radii, 1.8)
            elif atomic_number == 1:
                if atom.bond_partners[0].atomic_number in (6, 7):
                    self.assertEqual(radii, 1.3)
                elif atom.bond_partners[0].atomic_number in (8, 16):
                    self.assertEqual(radii, 0.8)
                else:
                    self.assertEqual(radii, 1.2)
            else:
                self.assertEqual(radii, 1.5)

        PT.changeRadii(parm, 'mbondi2').execute()
        self.assertEqual(parm.parm_data['RADIUS_SET'][0],
                         'H(N)-modified Bondi radii (mbondi2)')
        for atom in parm.atoms:
            radii, atomic_number = atom.radii, atom.atomic_number
            if atomic_number == 6:
                self.assertEqual(radii, 1.7)
            elif atomic_number == 7:
                self.assertEqual(radii, 1.55)
            elif atomic_number == 8 or atomic_number == 9:
                self.assertEqual(radii, 1.5)
            elif atomic_number == 14:
                self.assertEqual(radii, 2.1)
            elif atomic_number == 15:
                self.assertEqual(radii, 1.85)
            elif atomic_number == 16:
                self.assertEqual(radii, 1.8)
            elif atomic_number == 1:
                if atom.bond_partners[0].atomic_number == 7:
                    self.assertEqual(radii, 1.3)
                else:
                    self.assertEqual(radii, 1.2)
            else:
                self.assertEqual(radii, 1.5)

        PT.changeRadii(parm, 'mbondi3').execute()
        self.assertEqual(parm.parm_data['RADIUS_SET'][0],
                         'ArgH and AspGluO modified Bondi2 radii (mbondi3)')
        for i, atom in enumerate(parm.atoms):
            radii, atomic_number = atom.radii, atom.atomic_number
            if atomic_number == 6:
                self.assertEqual(radii, 1.7)
            elif atomic_number == 7:
                self.assertEqual(radii, 1.55)
            elif atomic_number == 8:
                if atom.residue.name in ('ASP,GLU') and (
                            atom.name.startswith('OD') or
                            atom.name.startswith('OE')):
                    self.assertEqual(radii, 1.4)
                elif atom.name == 'OXT' or (i < parm.ptr('natom') and
                            parm.atoms[i+1].name == 'OXT'):
                    self.assertEqual(radii, 1.4)
                else:
                    self.assertEqual(radii, 1.5)
            elif atomic_number == 9:
                self.assertEqual(radii, 1.5)
            elif atomic_number == 14:
                self.assertEqual(radii, 2.1)
            elif atomic_number == 15:
                self.assertEqual(radii, 1.85)
            elif atomic_number == 16:
                self.assertEqual(radii, 1.8)
            elif atomic_number == 1:
                if atom.residue.name == 'ARG' and \
                            atom.name[:2] in ('HH', 'HE'):
                    self.assertEqual(radii, 1.17)
                elif atom.bond_partners[0].atomic_number == 7:
                    self.assertEqual(radii, 1.3)
                else:
                    self.assertEqual(radii, 1.2)
            else:
                self.assertEqual(radii, 1.5)
        # Now test bad input
        self.assertRaises(exc.ChangeRadiiError, lambda:
                          PT.changeRadii(parm, 'mbondi6').execute())

    def testChangeLJPair(self):
        """ Test changeLJPair on AmberParm """
        parm = copy(gasparm)
        PT.changeLJPair(parm, '@%N', '@%H', 1.0, 1.0).execute()
        # Figure out what type numbers each atom type belongs to
        ntype = htype = 0
        for atom in parm.atoms:
            if atom.type == 'N':
                ntype = atom.nb_idx
            elif atom.type == 'H':
                htype = atom.nb_idx
        # Make sure the A and B coefficient matrices are what I expect them to
        # be
        indexes = sorted([ntype, htype])
        acoef = parm.parm_data['LENNARD_JONES_ACOEF']
        bcoef = parm.parm_data['LENNARD_JONES_BCOEF']
        refa = gasparm.parm_data['LENNARD_JONES_ACOEF']
        refb = gasparm.parm_data['LENNARD_JONES_BCOEF']
        ntypes = parm.ptr('ntypes')
        for i in range(ntypes):
            for j in range(i, ntypes):
                idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j]
                if [i+1, j+1] == indexes:
                    self.assertEqual(acoef[idx-1], 1.0)
                    self.assertEqual(bcoef[idx-1], 2.0)
                else:
                    self.assertEqual(acoef[idx-1], refa[idx-1])
                    self.assertEqual(bcoef[idx-1], refb[idx-1])

    def testChangeLJ14Pair(self):
        """ Check that changeLJ14Pair fails on AmberParm """
        parm = copy(gasparm)
        self.assertRaises(exc.ParmError, lambda:
            PT.changeLJ14Pair(parm, '@%N', '@%H', 1.0, 1.0).execute())

    def testChange(self):
        """ Test change on AmberParm with all properties """
        parm = copy(gasparm)
        PT.change(parm, 'CHARGE', ':ALA', 0, 'quiet').execute()
        for flag in parm.parm_data:
            if flag != 'CHARGE':
                self.assertEqual(parm.parm_data[flag], gasparm.parm_data[flag])
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(parm.parm_data['CHARGE'][i], atom.charge)
            if atom.residue.name == 'ALA':
                self.assertEqual(atom.charge, 0)
            else:
                self.assertEqual(atom.charge, gasparm.atoms[i].charge)
        PT.change(parm, 'MASS', ':GLY', 10.0).execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(parm.parm_data['MASS'][i], atom.mass)
            if atom.residue.name == 'GLY':
                self.assertEqual(atom.mass, 10.0)
            else:
                self.assertEqual(atom.mass, gasparm.atoms[i].mass)
            if atom.residue.name == 'ALA':
                self.assertEqual(atom.charge, 0.0)
            else:
                self.assertEqual(atom.charge, gasparm.atoms[i].charge)
        PT.change(parm, 'ATOM_NAME', ':ASP@C', 'JMS').execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(parm.parm_data['ATOM_NAME'][i], atom.name)
            if atom.residue.name == 'ASP' and gasparm.atoms[i].name == 'C':
                self.assertEqual(atom.name, 'JMS')
            else:
                self.assertEqual(atom.name, gasparm.atoms[i].name)
        PT.change(parm, 'AMBER_ATOM_TYPE', ':GLU@N', 'RJLS').execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(parm.parm_data['AMBER_ATOM_TYPE'][i], atom.type)
            if atom.residue.name == 'GLU' and gasparm.atoms[i].name == 'N':
                self.assertEqual(atom.type, 'RJLS')
            else:
                self.assertEqual(atom.type, gasparm.atoms[i].type)
        PT.change(parm, 'ATOM_TYPE_INDEX', '@1', 4).execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(atom.nb_idx, parm.parm_data['ATOM_TYPE_INDEX'][i])
            if i == 0:
                self.assertEqual(atom.nb_idx, 4)
            else:
                self.assertEqual(atom.nb_idx, gasparm.atoms[i].nb_idx)
        PT.change(parm, 'RADII', ':1-20', 2.0, 'quiet').execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(atom.radii, parm.parm_data['RADII'][i])
            if atom.residue.idx < 20:
                self.assertEqual(atom.radii, 2.0)
            else:
                self.assertEqual(atom.radii, gasparm.atoms[i].radii)
        PT.change(parm, 'SCREEN', '*', 0.0).execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(atom.screen, parm.parm_data['SCREEN'][i])
            self.assertEqual(atom.screen, 0.0)
        PT.change(parm, 'TREE_CHAIN_CLASSIFICATION', ':1-2', 'ABC').execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(atom.tree,
                             parm.parm_data['TREE_CHAIN_CLASSIFICATION'][i])
            if atom.residue.idx < 2:
                self.assertEqual(atom.tree, 'ABC')
        # Check bad input
        self.assertRaises(exc.ParmedChangeError, lambda:
                          PT.change(parm, 'RESIDUE_LABEL', ':*', 'NaN'))

    def testPrintInfo(self):
        """ Test printInfo for all flags on AmberParm """
        for flag in gasparm.parm_data:
            act = PT.printInfo(gasparm, flag)
            vals = []
            for line in str(act).split('\n'):
                vals += line.split()
            self.assertEqual(len(vals), len(gasparm.parm_data[flag]))
            try:
                datatype = type(gasparm.parm_data[flag][0])
            except IndexError:
                continue
            for i, j in zip(vals, gasparm.parm_data[flag]):
                # printInfo prints to 5 places for floats.
                if datatype is float:
                    self.assertAlmostEqual(datatype(i), j, places=4)
                else:
                    self.assertEqual(datatype(i), j)

    def testAddChangeLJType(self):
        """ Test addLJType and changeLJSingleType on AmberParm """
        parm = copy(gasparm)
        PT.addLJType(parm, '@1').execute()
        self.assertEqual(parm.ptr('ntypes'), gasparm.ptr('ntypes') + 1)
        self.assertEqual(parm.atoms[0].nb_idx, parm.ptr('ntypes'))
        ntypes = parm.ptr('ntypes')
        ntypes2 = ntypes - 1
        orig_type = gasparm.atoms[0].nb_idx - 1
        for i in range(ntypes):
            idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+ntypes-1]
            if i == ntypes - 1:
                idx2 = gasparm.parm_data['NONBONDED_PARM_INDEX'][
                                                ntypes2*orig_type+orig_type]
            else:
                ii, jj = sorted([orig_type, i])
                idx2 = gasparm.parm_data['NONBONDED_PARM_INDEX'][ntypes2*ii+jj]
            self.assertRelativeEqual(
                            parm.parm_data['LENNARD_JONES_ACOEF'][idx-1],
                            gasparm.parm_data['LENNARD_JONES_ACOEF'][idx2-1],
                            places=7)
            self.assertRelativeEqual(
                            parm.parm_data['LENNARD_JONES_BCOEF'][idx-1],
                            gasparm.parm_data['LENNARD_JONES_BCOEF'][idx2-1],
                            places=7)
        # Ensure that the rest of the values are unchanged (exactly equal)
        for i in range(ntypes2):
            for j in range(ntypes2):
                idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j]
                idx2 = gasparm.parm_data['NONBONDED_PARM_INDEX'][ntypes2*i+j]
                self.assertEqual(
                            parm.parm_data['LENNARD_JONES_ACOEF'][idx-1],
                            gasparm.parm_data['LENNARD_JONES_ACOEF'][idx2-1])
                self.assertEqual(
                            parm.parm_data['LENNARD_JONES_BCOEF'][idx-1],
                            gasparm.parm_data['LENNARD_JONES_BCOEF'][idx2-1])
        # Now supply keywords
        parm2 = copy(gasparm)
        PT.addLJType(parm2, '@1', radius=1.0, epsilon=1.0).execute()
        PT.changeLJSingleType(parm, '@1', 1.0, 1.0).execute()
        for x, y in zip(parm.parm_data['LENNARD_JONES_ACOEF'],
                        parm2.parm_data['LENNARD_JONES_ACOEF']):
            self.assertRelativeEqual(x, y)
        for x, y in zip(parm.parm_data['LENNARD_JONES_BCOEF'],
                        parm2.parm_data['LENNARD_JONES_BCOEF']):
            self.assertRelativeEqual(x, y)
        # Now use addLJType to hack a way to turn off LJ interactions
        PT.addLJType(parm, '*', radius=0.0, epsilon=0.0).execute()
        ntypes = parm.ptr('ntypes')
        for atom in parm.atoms:
            self.assertEqual(atom.nb_idx, ntypes)
        idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*(ntypes-1)+ntypes-1]
        self.assertEqual(parm.parm_data['LENNARD_JONES_ACOEF'][idx-1], 0.0)
        self.assertEqual(parm.parm_data['LENNARD_JONES_BCOEF'][idx-1], 0.0)

    def testPrintLJTypes(self):
        """ Test printLJTypes on AmberParm """
        # Simple test
        act = PT.printLJTypes(gasparm, '@1')
        for line in str(act).split('\n'):
            if not line.startswith('ATOM'):
                continue
            words = line.split()
            self.assertTrue(words[2].startswith('N'))
            self.assertTrue(words[3].startswith('N'))
            self.assertEqual(words[7], '1')

    def testSceeScnb(self):
        """ Test scee and scnb actions on AmberParm """
        parm = copy(gasparm)
        PT.scee(parm, 1.0).execute()
        PT.scnb(parm, 1.0).execute()
        for dih in parm.dihedrals:
            self.assertEqual(dih.type.scee, 1.0)
            self.assertEqual(dih.type.scnb, 1.0)
        for x, y in zip(parm.parm_data['SCEE_SCALE_FACTOR'],
                        parm.parm_data['SCNB_SCALE_FACTOR']):
            self.assertEqual(x, 1.0)
            self.assertEqual(y, 1.0)

    def testPrintDetails(self):
        """ Test printDetails on AmberParm """
        act = PT.printDetails(gasparm, '@1')
        self.assertEqual(str(act), saved.PRINT_DETAILS)

    def testPrintFlags(self):
        """ Test printFlags on AmberParm """
        act = PT.printFlags(gasparm)
        printed_flags = set()
        for line in str(act).split('\n'):
            if line.startswith('%FLAG'):
                printed_flags.add(line.split()[1])
        self.assertEqual(printed_flags, set(gasparm.parm_data.keys()))

    def testPrintPointers(self):
        """ Test printPointers on AmberParm """
        act = PT.printPointers(gasparm)
        printed_pointers = set()
        printed_pointers.add('NEXT') # Not in printed list
        for line in str(act).split('\n'):
            try:
                pointer = line.split()[0]
                value = int(line[line.rfind('=')+1:].strip())
            except (IndexError, ValueError):
                continue
            self.assertEqual(gasparm.ptr(pointer), value)
            printed_pointers.add(pointer)
        self.assertEqual(printed_pointers, set(gasparm.pointers.keys()))

    def testPrintBonds(self):
        """ Test printBonds on AmberParm """
        act = PT.printBonds(gasparm, '@1')
        self.assertEqual(str(act), saved.PRINT_BONDS)
        act = PT.printBonds(gasparm, '@1', ':*')
        self.assertEqual(str(act), saved.PRINT_BONDS)
        act = PT.printBonds(gasparm, '*', '@1')
        self.assertEqual(str(act), saved.PRINT_BONDS)
        act = PT.printBonds(gasparm, '@1', '@3')
        self.assertEqual(str(act), saved.PRINT_BONDS_2MASKS)
        act = PT.printBonds(gasparm, '@3', '@1')
        self.assertEqual(str(act), saved.PRINT_BONDS_2MASKS)

    def testPrintAngles(self):
        """ Test printAngles on AmberParm """
        act = PT.printAngles(gasparm, '@1')
        self.assertEqual(str(act), saved.PRINT_ANGLES)
        act = PT.printAngles(gasparm, '@1', '*')
        self.assertEqual(str(act), saved.PRINT_ANGLES_2MASKS_1)
        act = PT.printAngles(gasparm, '@1', '*', '*')
        self.assertEqual(str(act), saved.PRINT_ANGLES_2MASKS_1)
        act = PT.printAngles(gasparm, '*', '*', '@1')
        self.assertEqual(str(act), saved.PRINT_ANGLES_2MASKS_1)
        act = PT.printAngles(gasparm, '*', '@1')
        self.assertEqual(str(act), saved.PRINT_ANGLES_2MASKS_2)
        act = PT.printAngles(gasparm, '*', '@1', '*')
        self.assertEqual(str(act), saved.PRINT_ANGLES_2MASKS_2)
        act = PT.printAngles(gasparm, '@1', '@5')
        self.assertEqual(str(act), saved.PRINT_ANGLES_2MASKS_3)
        act = PT.printAngles(gasparm, '*', '@5', '@1')
        self.assertEqual(str(act), saved.PRINT_ANGLES_2MASKS_3)
        act = PT.printAngles(gasparm, '@1 @5 @7')
        self.assertEqual(str(act), saved.PRINT_ANGLES_3MASKS)
        act = PT.printAngles(gasparm, '@7 @5 @1')
        self.assertEqual(str(act), saved.PRINT_ANGLES_3MASKS)

    def testPrintDihedrals(self):
        """ Test printDihedrals on AmberParm """
        act = PT.printDihedrals(gasparm, '@1')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALS)
        act = PT.printDihedrals(gasparm, '@1', '@5')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALS_2MASKS)
        act = PT.printDihedrals(gasparm, '*', '*', '@5', '@1')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALS_2MASKS)
        act = PT.printDihedrals(gasparm, '@4', '@1', '@5')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALS_3MASKS)
        act = PT.printDihedrals(gasparm, '@4', '@1', '@5', '*')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALS_3MASKS)
        act = PT.printDihedrals(gasparm, '*', '@5', '@1', '@4')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALS_3MASKS)
        act = PT.printDihedrals(gasparm, '@7', ':1@CA', '@12', '@14')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALS_4MASKS)
        act = PT.printDihedrals(gasparm, '@14', '@12', ':1@CA', '@7')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALS_4MASKS)

    def testSetMolecules(self):
        """ Test setMolecules on AmberParm """
        parm = AmberParm(get_fn('things.parm7'), get_fn('things.rst7'))
        atoms = [atom for atom in parm.atoms] # shallow copy!
        self.assertTrue(all([x is y for x,y in zip(parm.atoms,atoms)]))
        self.assertEqual(parm.ptr('IPTRES'), 29)
        self.assertEqual(parm.ptr('NSPM'), 718)
        self.assertEqual(parm.ptr('NSPSOL'), 23)
        # To keep the output clean
        warnings.filterwarnings('ignore', category=AmberWarning)
        PT.setMolecules(parm).execute()
        self.assertFalse(all([x is y for x,y in zip(parm.atoms,atoms)]))
        # Now check that setMolecules can apply another time. solute_ions seems
        # to be broken, and I can't figure out why.
        PT.setMolecules(parm).execute()

    def testNetCharge(self):
        """ Test netCharge on AmberParm """
        act = PT.netCharge(gasparm)
        chg = act.execute() # check this part of the API
        self.assertEqual(str(act), 'The net charge of :* is %.4f' % chg)
        self.assertAlmostEqual(chg, -4.0, places=6)
        chg = PT.netCharge(gasparm, ':ASP').execute()
        self.assertAlmostEqual(chg, -10.0, places=6)

    def testStrip(self):
        """ Test stripping of AmberParm """
        parm = copy(gasparm)
        PT.strip(parm, ':1').execute()
        self.assertEqual(parm.ptr('natom'), 1641)
        self.assertEqual(len(parm.atoms), 1641)
        # Good enough for here. The strip action is repeatedly tested in the
        # core Amber test suite as part of the MM/PBSA tests via ante-MMPBSA.py
        # and that part also tests that the energies come out correct as well

    def testDefineSolvent(self):
        """ Test defineSolvent on AmberParm """
        import parmed.residue as residue
        PT.defineSolvent(gasparm, 'WAT,HOH,Na+,Cl-').execute()
        self.assertEqual(residue.SOLVENT_NAMES, 'WAT HOH Na+ Cl-'.split())
        PT.defineSolvent(gasparm, 'WAT,HOH').execute()
        self.assertEqual(residue.SOLVENT_NAMES, 'WAT HOH'.split())

    def testAddExclusions(self):
        """ Test addExclusions on AmberParm """
        parm = copy(gasparm)
        in_exclusions_before = []
        for atom1 in parm.residues[0].atoms:
            all_exclusions = (atom1.bond_partners + atom1.angle_partners + 
                             atom1.dihedral_partners + atom1.exclusion_partners)
            for atom2 in parm.residues[0].atoms:
                if atom1 is atom2: continue
                in_exclusions_before.append(atom2 in all_exclusions)
        self.assertFalse(all(in_exclusions_before))
        PT.addExclusions(parm, ':1', ':1').execute()
        in_exclusions_after = []
        for atom1 in parm.residues[0].atoms:
            all_exclusions = (atom1.bond_partners + atom1.angle_partners + 
                             atom1.dihedral_partners + atom1.exclusion_partners)
            for atom2 in parm.residues[0].atoms:
                if atom1 is atom2: continue
                in_exclusions_after.append(atom2 in all_exclusions)
                if not in_exclusions_after[-1]:
                    print('%s %s not excluded' % (atom1, atom2))
        self.assertTrue(all(in_exclusions_after))

    def testAddDeleteDihedral(self):
        """ Test addDihedral and deleteDihedral on AmberParm """
        parm = copy(gasparm)
        n = PT.deleteDihedral(parm, ':ALA@N :ALA@CA :ALA@CB :ALA@HB1').execute()
        parm.remake_parm()
        self.assertEqual(gasparm.ptr('nphih') + gasparm.ptr('nphia'),
                         parm.ptr('nphih') + parm.ptr('nphia') + n)
        NALA = sum([res.name == 'ALA' for res in parm.residues])
        self.assertEqual(n, NALA)
        PT.addDihedral(parm, ':ALA@N', ':ALA@CA', ':ALA@CB', ':ALA@HB1',
                       0.1556, 3, 0, 1.2, 2.0).execute()
        parm.remake_parm()
        self.assertEqual(gasparm.ptr('nphih') + gasparm.ptr('nphia'),
                         parm.ptr('nphih') + parm.ptr('nphia'))
        PT.addDihedral(parm, ':ALA@N', ':ALA@CA', ':ALA@CB', ':ALA@HB1',
                       0.1556, 1, 0, 1.2, 2.0, type='normal').execute()
        parm.remake_parm()
        self.assertEqual(gasparm.ptr('nphih') + gasparm.ptr('nphia'),
                         parm.ptr('nphih') + parm.ptr('nphia') - n)
        num_dihedrals = 0
        num_ignore_ends = 0
        for atom in parm.atoms:
            if atom.residue.name == 'ALA' and atom.name == 'N':
                for dih in atom.dihedrals:
                    if dih.atom1 is atom:
                        if (dih.atom2.name == 'CA' and dih.atom3.name == 'CB'
                            and dih.atom4.name == 'HB1'):
                            num_dihedrals += 1
                            if dih.ignore_end:
                                num_ignore_ends += 1
                    elif dih.atom4 is atom:
                        if (dih.atom2.name == 'CB' and dih.atom3.name == 'CA' and
                            dih.atom1.name == 'HB1'):
                            num_dihedrals += 1
                            if dih.ignore_end:
                                num_ignore_ends += 1
        self.assertEqual(num_dihedrals, 2*NALA)
        self.assertEqual(num_ignore_ends, 1*NALA)

    def testSetBond(self):
        """ Test setBond on AmberParm """
        parm = copy(gasparm)
        PT.setBond(parm, ':ALA@CA', ':ALA@CB', 300.0, 1.5).execute()
        act = PT.printBonds(parm, ':ALA@CA')
        self.assertEqual(str(act), saved.SET_BOND)
        nala = sum([1 for res in parm.residues if res.name == 'ALA'])
        nbon = len(parm.bonds)
        PT.setBond(parm, ':ALA@CB', ':ALA@HA', 100.0, 1.5).execute()
        for atom in parm.atoms:
            if atom.residue.name == 'ALA' and atom.name == 'HA':
                self.assertEqual(len(atom.bonds), 2)
                for atom2 in atom.bond_partners:
                    self.assertIn(atom2.name, ('CA','CB'))
        self.assertEqual(nbon + nala, len(parm.bonds))

    def testSetAngle(self):
        """ Test setAngle on AmberParm """
        parm = copy(gasparm)
        PT.setAngle(parm, ':ALA@CA', ':ALA@CB', ':ALA@HB1', 40, 100).execute()
        act = PT.printAngles(parm, ':ALA@CB')
        self.assertEqual(str(act), saved.SET_ANGLE)
        nala = sum([1 for res in parm.residues if res.name == 'ALA'])
        nang = len(parm.angles)
        PT.setAngle(parm, ':ALA@HA', ':ALA@CB', ':ALA@HB1', 50, 120).execute()
        self.assertEqual(nang + nala, len(parm.angles))

    def testAddAtomicNumber(self):
        """ Test addAtomicNumber on AmberParm """
        parm = copy(gasparm)
        self.assertFalse('ATOMIC_NUMBER' in parm.parm_data)
        atomic_numbers = [atom.atomic_number for atom in parm.atoms]
        PT.addAtomicNumber(parm).execute()
        self.assertEqual(parm.parm_data['ATOMIC_NUMBER'], atomic_numbers)

    def testPrintLJMatrix(self):
        """ Test printLJMatrix on AmberParm """
        act = PT.printLJMatrix(gasparm, '@1')
        self.assertEqual(str(act), saved.PRINT_LJMATRIX)

    def testDeleteBond(self):
        """ Test deleteBond on AmberParm """
        parm = copy(gasparm)
        # Pick the bond we plan to delete, pick out every angle and dihedral
        # that contains that bond, and then delete it. Then make sure none of
        # the valence terms that contained that bond remain afterwards. We
        # already have a test to make sure that the __contains__ method works
        # for atoms and bonds.
        for bond in parm.atoms[0].bonds:
            if parm.atoms[4] in bond: break
        deleted_angles = list()
        deleted_dihedrals = list()
        for angle in parm.angles:
            if bond in angle: deleted_angles.append(angle)
        for dihedral in parm.dihedrals:
            if bond in dihedral: deleted_dihedrals.append(dihedral)
        PT.deleteBond(parm, '@1', '@5').execute()
        self.assertTrue(bond not in parm.bonds)
        for angle in deleted_angles:
            self.assertTrue(angle not in parm.angles)
        for dihedral in deleted_dihedrals:
            self.assertTrue(dihedral not in parm.dihedrals)

    def testSummary(self):
        """ Test summary action on AmberParm """
        parm = AmberParm(get_fn('things.parm7'), get_fn('things.rst7'))
        act = PT.summary(parm)
        self.assertEqual(str(act), saved.SUMMARY)

    def testScale(self):
        """ Test scale action on AmberParm """
        parm = copy(gasparm)
        PT.scale(parm, 'DIHEDRAL_FORCE_CONSTANT', 2.0).execute()
        self.assertEqual(
                [2*x for x in gasparm.parm_data['DIHEDRAL_FORCE_CONSTANT']],
                parm.parm_data['DIHEDRAL_FORCE_CONSTANT'])
        PT.scale(parm, 'DIHEDRAL_FORCE_CONSTANT', 0.0).execute()
        for val in parm.parm_data['DIHEDRAL_FORCE_CONSTANT']:
            self.assertEqual(val, 0)

    def testLmod(self):
        """ Test lmod action on AmberParm """
        parm = copy(gasparm)
        self.assertFalse(all(parm.parm_data['LENNARD_JONES_ACOEF']))
        PT.lmod(parm).execute()
        self.assertTrue(all(parm.parm_data['LENNARD_JONES_ACOEF']))

    def testProtStateInterpolate(self):
        """ Test changeProtState and interpolate actions on AmberParm """
        self._empty_writes()
        parm = AmberParm(get_fn('ash.parm7'))
        origparm = copy(parm)
        origparm.name = origparm.name + '_copy1'
        PT.changeProtState(parm, ':ASH', 0).execute()
        self.assertAlmostEqual(sum(parm.parm_data['CHARGE']), -1)
        self.assertAlmostEqual(sum(origparm.parm_data['CHARGE']), 0)
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(atom.charge, parm.parm_data['CHARGE'][i])
        for i, atom in enumerate(origparm.atoms):
            self.assertEqual(atom.charge, origparm.parm_data['CHARGE'][i])
        # Now set up a ParmList so we can interpolate these topology files
        parms = parmlist.ParmList()
        parms.add_parm(parm)
        parms.add_parm(origparm)
        sys.stdout = open(os.devnull, 'w')
        PT.interpolate(parms, 5, 'eleconly', startnum=2,
                       prefix=get_fn('test.parm7', written=True)).execute()
        sys.stdout = sys.__stdout__
        self.assertEqual(len(os.listdir(get_fn('writes'))), 5)
        self.assertTrue(os.path.exists(get_fn('test.parm7.2', written=True)))
        self.assertTrue(os.path.exists(get_fn('test.parm7.3', written=True)))
        self.assertTrue(os.path.exists(get_fn('test.parm7.4', written=True)))
        self.assertTrue(os.path.exists(get_fn('test.parm7.5', written=True)))
        self.assertTrue(os.path.exists(get_fn('test.parm7.6', written=True)))
        # Now check them all
        ladder = [origparm]
        ladder.append(AmberParm(get_fn('test.parm7.2', written=True)))
        ladder.append(AmberParm(get_fn('test.parm7.3', written=True)))
        ladder.append(AmberParm(get_fn('test.parm7.4', written=True)))
        ladder.append(AmberParm(get_fn('test.parm7.5', written=True)))
        ladder.append(AmberParm(get_fn('test.parm7.6', written=True)))
        ladder.append(parm)
        natom = parm.ptr('natom')
        for i in range(1, 6):
            before = ladder[i-1].parm_data['CHARGE']
            after = ladder[i+1].parm_data['CHARGE']
            this = ladder[i].parm_data['CHARGE']
            for j in range(natom):
                if before[j] < after[j]:
                    self.assertTrue(before[j] <= this[j] <= after[j])
                else:
                    self.assertTrue(before[j] >= this[j] >= after[j])

    def testAddDeletePDB(self):
        """ Test addPDB and deletePDB actions on AmberParm """
        parm = copy(gasparm)
        PT.addPDB(parm, get_fn('trx.pdb'), 'elem', 'allicodes').execute()
        self.assertTrue('RESIDUE_ICODE' in parm.flag_list)
        self.assertTrue('ATOM_ELEMENT' in parm.flag_list)
        self.assertTrue('RESIDUE_NUMBER' in parm.flag_list)
        self.assertTrue('RESIDUE_CHAINID' in parm.flag_list)
        self.assertTrue('ATOM_OCCUPANCY' in parm.flag_list)
        self.assertTrue('ATOM_BFACTOR' in parm.flag_list)
        self.assertTrue(len(parm.parm_data['RESIDUE_ICODE']), parm.ptr('nres'))
        self.assertTrue(len(parm.parm_data['ATOM_ELEMENT']), parm.ptr('natom'))
        self.assertTrue(len(parm.parm_data['RESIDUE_NUMBER']), parm.ptr('nres'))
        self.assertTrue(len(parm.parm_data['RESIDUE_CHAINID']),parm.ptr('nres'))
        self.assertTrue(len(parm.parm_data['ATOM_OCCUPANCY']),parm.ptr('natom'))
        self.assertTrue(len(parm.parm_data['ATOM_BFACTOR']), parm.ptr('natom'))
        for i in range(parm.ptr('nres')):
            self.assertEqual(parm.parm_data['RESIDUE_NUMBER'][i], i + 21)
            self.assertEqual(parm.parm_data['RESIDUE_ICODE'][i], '')
            if parm.parm_data['RESIDUE_NUMBER'][i] < 41:
                self.assertEqual(parm.parm_data['RESIDUE_CHAINID'][i], 'A')
            else:
                self.assertEqual(parm.parm_data['RESIDUE_CHAINID'][i], 'B')
        for i, atom in enumerate(parm.atoms):
            atnum = atom.atomic_number
            elem = parm.parm_data['ATOM_ELEMENT'][i].strip()
            self.assertEqual(periodic_table.Element[atnum], elem)
            self.assertEqual(atnum, periodic_table.AtomicNum[elem])
        PT.deletePDB(parm).execute()
        self.assertFalse('RESIDUE_ICODE' in parm.flag_list)
        self.assertFalse('ATOM_ELEMENT' in parm.flag_list)
        self.assertFalse('RESIDUE_NUMBER' in parm.flag_list)
        self.assertFalse('RESIDUE_CHAINID' in parm.flag_list)
        self.assertFalse('ATOM_OCCUPANCY' in parm.flag_list)
        self.assertFalse('ATOM_BFACTOR' in parm.flag_list)

    def testAddPDB2(self):
        """ Test addPDB with atypical numbering and extra residues """
        parm = load_file(get_fn('4lzt.parm7'))
        PT.addPDB(parm, get_fn('4lzt_NoNO3.pdb')).execute()
        parm.write_parm(get_fn('4lzt_pdb.parm7', written=True))
        self.assertTrue(diff_files(get_saved_fn('4lzt_pdb.parm7'),
                                   get_fn('4lzt_pdb.parm7', written=True),
                                   absolute_error=1e-6)
        )

    def testHMassRepartition(self):
        """ Test HMassRepartition on AmberParm """
        parm = copy(solvparm)
        PT.HMassRepartition(parm, 2.0).execute()
        for atom in parm.atoms:
            if atom.atomic_number == 1:
                if atom.residue.name == 'WAT':
                    self.assertAlmostEqual(atom.mass, 1.008)
                else:
                    self.assertEqual(atom.mass, 2)
        self.assertEqual(parm.parm_data['MASS'],
                         [a.mass for a in parm.atoms])
        self.assertAlmostEqual(sum(solvparm.parm_data['MASS']),
                               sum(parm.parm_data['MASS']))
        PT.HMassRepartition(parm, 3.0, 'dowater').execute()
        for atom in parm.atoms:
            if atom.atomic_number == 1:
                self.assertEqual(atom.mass, 3.0)
        self.assertAlmostEqual(sum(solvparm.parm_data['MASS']),
                               sum(parm.parm_data['MASS']))

    def testOutPDB(self):
        """ Test the outPDB action on AmberParm """
        parm = copy(gasparm)
        PT.loadRestrt(parm, get_fn('trx.inpcrd')).execute()
        PT.outPDB(parm, get_fn('outPDB1.pdb', written=True)).execute()
        f = PDBFile.parse(get_fn('outPDB1.pdb', written=True))
        self.assertEqual(len(f.atoms), len(parm.atoms))
        self.assertEqual(len(f.residues), len(parm.residues))
        for a1, a2 in zip(f.atoms, parm.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertAlmostEqual(a1.xx, a2.xx, delta=2e-3)
            self.assertAlmostEqual(a1.xy, a2.xy, delta=2e-3)
            self.assertAlmostEqual(a1.xz, a2.xz, delta=2e-3)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.idx, a2.residue.idx)

    def testOutCIF(self):
        """ Test the outCIF action on AmberParm """
        parm = copy(gasparm)
        PT.loadRestrt(parm, get_fn('trx.inpcrd')).execute()
        PT.outCIF(parm, get_fn('outPDB1.cif', written=True)).execute()
        f = CIFFile.parse(get_fn('outPDB1.cif', written=True))
        self.assertEqual(len(f.atoms), len(parm.atoms))
        self.assertEqual(len(f.residues), len(parm.residues))
        for a1, a2 in zip(f.atoms, parm.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertAlmostEqual(a1.xx, a2.xx, delta=2e-3)
            self.assertAlmostEqual(a1.xy, a2.xy, delta=2e-3)
            self.assertAlmostEqual(a1.xz, a2.xz, delta=2e-3)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.idx, a2.residue.idx)

    def testTIMerge(self):
        """ Tests the tiMerge action on AmberParm """
        parm = AmberParm(get_fn('abs.prmtop'), get_fn('abs.inpcrd'))
        PT.tiMerge(parm, ':1-3', ':4-6', ':2', ':5').execute()
        parm.write_parm(get_fn('abs_merged.prmtop', written=True))
        parm.write_rst7(get_fn('abs_merged.inpcrd', written=True))
        self.assertTrue(diff_files(get_fn('abs_merged.prmtop', written=True),
                                   get_saved_fn('abs_merged.prmtop')))
        self.assertTrue(diff_files(get_fn('abs_merged.inpcrd', written=True),
                                   get_saved_fn('abs_merged.inpcrd')))

    def testAdd1264(self):
        """ Test the add12_6_4 action on AmberParm """
        parm = AmberParm(get_fn('Mg_ti1_b.parm7'))
        PT.addLJType(parm, '@14').execute()
        PT.changeLJPair(parm, '@14', ':MG', 3.26, 0.061666).execute()
        PT.add12_6_4(parm, ':MG', watermodel='TIP4PEW',
                     polfile=get_fn('lj_1264_pol.dat')).execute()
        parm.write_parm(get_fn('Mg_ti1_b_1264.parm7', written=True))
        self.assertTrue(diff_files(get_fn('Mg_ti1_b_1264.parm7', written=True),
                                   get_saved_fn('Mg_ti1_b_1264.parm7'))
        )

    def testAdd1264_2metals(self):
        """ Test the add12_6_4 action on AmberParm with 2+ metals """
        parm1 = AmberParm(get_fn('mg_na_cl.parm7'))
        parm2 = AmberParm(get_fn('na_cl_mg.parm7'))
        PT.add12_6_4(parm1, ':MG,NA,CL', watermodel='TIP3P',
                     polfile=get_fn('lj_1264_pol.dat')).execute()
        PT.add12_6_4(parm2, ':MG,NA,CL', watermodel='TIP3P',
                     polfile=get_fn('lj_1264_pol.dat')).execute()
        self.assertEqual(str(PT.printLJMatrix(parm1, ':MG')),
                         saved.PRINTLJMATRIX_MGNACL)
        self.assertEqual(str(PT.printLJMatrix(parm2, ':MG')),
                         saved.PRINTLJMATRIX_NACLMG)

class TestChamberParmActions(utils.TestCaseRelative, utils.FileIOTestCase):
    """ Tests actions on Amber prmtop files """
    
    def testParmoutOutparmLoadRestrt(self):
        """ Test parmout, outparm, and loadRestrt actions for ChamberParm """
        self._empty_writes()
        parm = copy(gascham)
        PT.loadRestrt(parm, get_fn('ala_ala_ala.rst7')).execute()
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))
        PT.parmout(parm, get_fn('test.parm7', written=True)).execute()
        self.assertEqual(len(os.listdir(get_fn('writes'))), 1)
        self.assertTrue(diff_files(get_fn('ala_ala_ala.parm7'),
                                   get_fn('test.parm7', written=True),
                                   absolute_error=1e-6))
        self._empty_writes()
        PT.parmout(parm, get_fn('test.parm7', written=True),
                         get_fn('test.rst7', written=True)).execute()
        self.assertEqual(len(os.listdir(get_fn('writes'))), 2)
        self.assertTrue(diff_files(get_fn('ala_ala_ala.parm7'),
                                   get_fn('test.parm7', written=True),
                                   absolute_error=1e-6))
        self.assertTrue(diff_files(get_fn('ala_ala_ala.rst7'),
                                   get_fn('test.rst7', written=True),
                                   absolute_error=0.0001))
        self._empty_writes()
        PT.outparm(parm, get_fn('test.parm7', written=True)).execute()
        self.assertEqual(len(os.listdir(get_fn('writes'))), 1)
        self.assertTrue(diff_files(get_fn('ala_ala_ala.parm7'),
                                   get_fn('test.parm7', written=True),
                                   absolute_error=1e-6))
        self._empty_writes()
        PT.outparm(parm, get_fn('test.parm7', written=True),
                         get_fn('test.rst7', written=True)).execute()
        self.assertEqual(len(os.listdir(get_fn('writes'))), 2)
        self.assertTrue(diff_files(get_fn('ala_ala_ala.parm7'),
                                   get_fn('test.parm7', written=True),
                                   absolute_error=1e-6))
        self.assertTrue(diff_files(get_fn('ala_ala_ala.rst7'),
                                   get_fn('test.rst7', written=True),
                                   absolute_error=0.0001))

    def testWriteFrcmod(self):
        """ Check that writeFrcmod fails for ChamberParm """
        parm = gascham
        self.assertRaises(exc.ParmError, lambda:
                PT.writeFrcmod(parm, get_fn('x', written=True)).execute())

    def testWriteOffLoadRestrt(self):
        """ Check that writeOFF fails for ChamberParm """
        parm = copy(gascham)
        PT.loadRestrt(parm, get_fn('ala_ala_ala.rst7')).execute()
        self.assertRaises(exc.ParmError, lambda:
                PT.writeOFF(parm, get_fn('test.off', written=True)).execute())

    def testTIMerge(self):
        """ Check that tiMerge fails for ChamberParm """
        parm = copy(gascham)
        PT.loadRestrt(parm, get_fn('ala_ala_ala.rst7')).execute()
        self.assertRaises(exc.ParmError, lambda:
                PT.tiMerge(parm, ':1-3', ':4-6', ':2', ':5').execute())

    def testChangeRadii(self):
        """ Test changeRadii for ChamberParm """
        parm = copy(gascham)
        PT.changeRadii(parm, 'amber6').execute()
        self.assertEqual(parm.parm_data['RADIUS_SET'][0],
                         'amber6 modified Bondi radii (amber6)')
        for i, atom in enumerate(parm.atoms):
            radii, atomic_number = atom.radii, atom.atomic_number
            self.assertEqual(parm.parm_data['RADII'][i], radii)
            if atomic_number == 6:
                self.assertEqual(radii, 1.7)
            elif atomic_number == 7:
                self.assertEqual(radii, 1.55)
            elif atomic_number == 8 or atomic_number == 9:
                self.assertEqual(radii, 1.5)
            elif atomic_number == 14:
                self.assertEqual(radii, 2.1)
            elif atomic_number == 15:
                self.assertEqual(radii, 1.85)
            elif atomic_number == 16:
                self.assertEqual(radii, 1.8)
            elif atomic_number == 1:
                if atom.bond_partners[0].atomic_number == 6:
                    self.assertEqual(radii, 1.3)
                elif atom.bond_partners[0].atomic_number in (8, 16):
                    self.assertEqual(radii, 0.8)
                else:
                    self.assertEqual(radii, 1.2)
            else:
                self.assertEqual(radii, 1.5)

        PT.changeRadii(parm, 'bondi').execute()
        self.assertEqual(parm.parm_data['RADIUS_SET'][0], 'Bondi radii (bondi)')
        for atom in parm.atoms:
            radii, atomic_number = atom.radii, atom.atomic_number
            if atomic_number == 6:
                self.assertEqual(radii, 1.7)
            elif atomic_number == 7:
                self.assertEqual(radii, 1.55)
            elif atomic_number == 8 or atomic_number == 9:
                self.assertEqual(radii, 1.5)
            elif atomic_number == 14:
                self.assertEqual(radii, 2.1)
            elif atomic_number == 15:
                self.assertEqual(radii, 1.85)
            elif atomic_number == 16:
                self.assertEqual(radii, 1.8)
            elif atomic_number == 1:
                self.assertEqual(radii, 1.2)
            else:
                self.assertEqual(radii, 1.5)

        PT.changeRadii(parm, 'mbondi').execute()
        self.assertEqual(parm.parm_data['RADIUS_SET'][0],
                         'modified Bondi radii (mbondi)')
        for atom in parm.atoms:
            radii, atomic_number = atom.radii, atom.atomic_number
            if atomic_number == 6:
                self.assertEqual(radii, 1.7)
            elif atomic_number == 7:
                self.assertEqual(radii, 1.55)
            elif atomic_number == 8 or atomic_number == 9:
                self.assertEqual(radii, 1.5)
            elif atomic_number == 14:
                self.assertEqual(radii, 2.1)
            elif atomic_number == 15:
                self.assertEqual(radii, 1.85)
            elif atomic_number == 16:
                self.assertEqual(radii, 1.8)
            elif atomic_number == 1:
                if atom.bond_partners[0].atomic_number in (6, 7):
                    self.assertEqual(radii, 1.3)
                elif atom.bond_partners[0].atomic_number in (8, 16):
                    self.assertEqual(radii, 0.8)
                else:
                    self.assertEqual(radii, 1.2)
            else:
                self.assertEqual(radii, 1.5)

        PT.changeRadii(parm, 'mbondi2').execute()
        self.assertEqual(parm.parm_data['RADIUS_SET'][0],
                         'H(N)-modified Bondi radii (mbondi2)')
        for atom in parm.atoms:
            radii, atomic_number = atom.radii, atom.atomic_number
            if atomic_number == 6:
                self.assertEqual(radii, 1.7)
            elif atomic_number == 7:
                self.assertEqual(radii, 1.55)
            elif atomic_number == 8 or atomic_number == 9:
                self.assertEqual(radii, 1.5)
            elif atomic_number == 14:
                self.assertEqual(radii, 2.1)
            elif atomic_number == 15:
                self.assertEqual(radii, 1.85)
            elif atomic_number == 16:
                self.assertEqual(radii, 1.8)
            elif atomic_number == 1:
                if atom.bond_partners[0].atomic_number == 7:
                    self.assertEqual(radii, 1.3)
                else:
                    self.assertEqual(radii, 1.2)
            else:
                self.assertEqual(radii, 1.5)

        PT.changeRadii(parm, 'mbondi3').execute()
        self.assertEqual(parm.parm_data['RADIUS_SET'][0],
                         'ArgH and AspGluO modified Bondi2 radii (mbondi3)')
        for i, atom in enumerate(parm.atoms):
            radii, atomic_number = atom.radii, atom.atomic_number
            if atomic_number == 6:
                self.assertEqual(radii, 1.7)
            elif atomic_number == 7:
                self.assertEqual(radii, 1.55)
            elif atomic_number == 8:
                if atom.residue.name in ('ASP,GLU') and (
                            atom.name.startswith('OD') or
                            atom.name.startswith('OE')):
                    self.assertEqual(radii, 1.4)
                elif atom.name == 'OXT' or (i < parm.ptr('natom')-1 and
                            parm.atoms[i+1].name == 'OXT'):
                    self.assertEqual(radii, 1.4)
                else:
                    self.assertEqual(radii, 1.5)
            elif atomic_number == 9:
                self.assertEqual(radii, 1.5)
            elif atomic_number == 14:
                self.assertEqual(radii, 2.1)
            elif atomic_number == 15:
                self.assertEqual(radii, 1.85)
            elif atomic_number == 16:
                self.assertEqual(radii, 1.8)
            elif atomic_number == 1:
                if atom.residue.name == 'ARG' and \
                            atom.name[:2] in ('HH', 'HE'):
                    self.assertEqual(radii, 1.17)
                elif atom.bond_partners[0].atomic_number == 7:
                    self.assertEqual(radii, 1.3)
                else:
                    self.assertEqual(radii, 1.2)
            else:
                self.assertEqual(radii, 1.5)
        # Now test bad input
        self.assertRaises(exc.ChangeRadiiError, lambda:
                          PT.changeRadii(parm, 'mbondi6').execute())

    def testChangeLJPair(self):
        """ Test changeLJPair for ChamberParm """
        parm = copy(gascham)
        PT.changeLJPair(parm, '@%NH3', '@%HC', 1.0, 1.0).execute()
        # Figure out what type numbers each atom type belongs to
        ntype = htype = 0
        for atom in parm.atoms:
            if atom.type == 'NH3':
                ntype = atom.nb_idx
            elif atom.type == 'HC':
                htype = atom.nb_idx
        # Make sure the A and B coefficient matrices are what I expect them to
        # be
        indexes = sorted([ntype, htype])
        acoef = parm.parm_data['LENNARD_JONES_ACOEF']
        bcoef = parm.parm_data['LENNARD_JONES_BCOEF']
        refa = gascham.parm_data['LENNARD_JONES_ACOEF']
        refb = gascham.parm_data['LENNARD_JONES_BCOEF']
        ntypes = parm.ptr('ntypes')
        for i in range(ntypes):
            for j in range(i, ntypes):
                idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j]
                if [i+1, j+1] == indexes:
                    self.assertEqual(acoef[idx-1], 1.0)
                    self.assertEqual(bcoef[idx-1], 2.0)
                else:
                    self.assertEqual(acoef[idx-1], refa[idx-1])
                    self.assertEqual(bcoef[idx-1], refb[idx-1])

    def testChangeLJ14Pair(self):
        """ Test changeLJ14Pair for ChamberParm """
        parm = copy(gascham)
        PT.changeLJ14Pair(parm, '@%NH3', '@%HC', 1.0, 1.0).execute()
        # Figure out what type numbers each atom type belongs to
        ntype = htype = 0
        for atom in parm.atoms:
            if atom.type == 'NH3':
                ntype = atom.nb_idx
            elif atom.type == 'HC':
                htype = atom.nb_idx
        # Make sure the A and B coefficient matrices are what I expect them to
        # be
        indexes = sorted([ntype, htype])
        acoef = parm.parm_data['LENNARD_JONES_14_ACOEF']
        bcoef = parm.parm_data['LENNARD_JONES_14_BCOEF']
        refa = gascham.parm_data['LENNARD_JONES_14_ACOEF']
        refb = gascham.parm_data['LENNARD_JONES_14_BCOEF']
        ntypes = parm.ptr('ntypes')
        for i in range(ntypes):
            for j in range(i, ntypes):
                idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j]
                if [i+1, j+1] == indexes:
                    self.assertEqual(acoef[idx-1], 1.0)
                    self.assertEqual(bcoef[idx-1], 2.0)
                else:
                    self.assertEqual(acoef[idx-1], refa[idx-1])
                    self.assertEqual(bcoef[idx-1], refb[idx-1])

    def testChange(self):
        """ Test change on ChamberParm with all properties """
        parm = copy(gascham)
        PT.change(parm, 'CHARGE', ':ALA', 0, 'quiet').execute()
        for flag in parm.parm_data:
            if flag != 'CHARGE':
                self.assertEqual(parm.parm_data[flag], gascham.parm_data[flag])
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(parm.parm_data['CHARGE'][i], atom.charge)
            if atom.residue.name == 'ALA':
                self.assertEqual(atom.charge, 0)
            else:
                self.assertEqual(atom.charge, gascham.parm_data['CHARGE'][i])
        PT.change(parm, 'MASS', ':1', 10.0).execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(parm.parm_data['MASS'][i], atom.mass)
            if atom.residue.idx == 0:
                self.assertEqual(atom.mass, 10.0)
            else:
                self.assertEqual(atom.mass, gascham.parm_data['MASS'][i])
            if atom.residue.name == 'ALA':
                self.assertEqual(atom.charge, 0.0)
            else:
                self.assertEqual(atom.charge, gascham.parm_data['CHARGE'][i])
        PT.change(parm, 'ATOM_NAME', ':ALA@C', 'JMS').execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(parm.parm_data['ATOM_NAME'][i], atom.name)
            if atom.residue.name == 'ALA' and \
                        gascham.parm_data['ATOM_NAME'][i] == 'C':
                self.assertEqual(atom.name, 'JMS')
            else:
                self.assertEqual(atom.name, gascham.parm_data['ATOM_NAME'][i])
        PT.change(parm, 'AMBER_ATOM_TYPE', ':ALA@%NH3', 'RJLS').execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(parm.parm_data['AMBER_ATOM_TYPE'][i], atom.type)
            if atom.residue.name == 'ALA' and \
                        gascham.parm_data['AMBER_ATOM_TYPE'][i] == 'NH3':
                self.assertEqual(atom.type, 'RJLS')
            else:
                self.assertEqual(atom.type,
                                 gascham.parm_data['AMBER_ATOM_TYPE'][i])
        PT.change(parm, 'ATOM_TYPE_INDEX', '@1', 4).execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(atom.nb_idx, parm.parm_data['ATOM_TYPE_INDEX'][i])
            if i == 0:
                self.assertEqual(atom.nb_idx, 4)
            else:
                self.assertEqual(atom.nb_idx,
                                 gascham.parm_data['ATOM_TYPE_INDEX'][i])
        PT.change(parm, 'RADII', ':1-2', 2.0, 'quiet').execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(atom.radii, parm.parm_data['RADII'][i])
            if atom.residue.idx < 2:
                self.assertEqual(atom.radii, 2.0)
            else:
                self.assertEqual(atom.radii, gascham.parm_data['RADII'][i])
        PT.change(parm, 'SCREEN', '*', 0.0).execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(atom.screen, parm.parm_data['SCREEN'][i])
            self.assertEqual(atom.screen, 0.0)
        # Check bad input
        self.assertRaises(exc.ParmedChangeError, lambda:
                          PT.change(parm, 'RESIDUE_LABEL', ':*', 'NaN'))

    def testPrintInfo(self):
        """ Test printInfo on ChamberParm for all FLAGs """
        for flag in gascham.parm_data:
            if flag == 'FORCE_FIELD_TYPE': continue
            if flag == 'RADIUS_SET': continue
            act = PT.printInfo(gascham, flag)
            vals = []
            for line in str(act).split('\n'):
                vals += line.split()
            self.assertEqual(len(vals), len(gascham.parm_data[flag]))
            try:
                datatype = type(gascham.parm_data[flag][0])
            except IndexError:
                continue
            for i, j in zip(vals, gascham.parm_data[flag]):
                # printInfo prints to 5 places for floats.
                if datatype is float:
                    self.assertAlmostEqual(datatype(i), j, places=4)
                else:
                    self.assertEqual(datatype(i), j)

    def testAddChangeLJType(self):
        """ Test addLJType and changeLJSingleType on ChamberParm """
        parm = copy(gascham)
        PT.addLJType(parm, '@1').execute()
        self.assertEqual(parm.ptr('ntypes'), gascham.ptr('ntypes') + 1)
        self.assertEqual(parm.atoms[0].nb_idx, parm.ptr('ntypes'))
        ntypes = parm.ptr('ntypes')
        ntypes2 = ntypes - 1
        orig_type = gascham.atoms[0].nb_idx - 1
        for i in range(ntypes):
            idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+ntypes-1]
            if i == ntypes - 1:
                idx2 = gascham.parm_data['NONBONDED_PARM_INDEX'][
                                                ntypes2*orig_type+orig_type]
            else:
                ii, jj = sorted([orig_type, i])
                idx2 = gascham.parm_data['NONBONDED_PARM_INDEX'][ntypes2*ii+jj]
            self.assertRelativeEqual(
                            parm.parm_data['LENNARD_JONES_ACOEF'][idx-1],
                            gascham.parm_data['LENNARD_JONES_ACOEF'][idx2-1],
                            places=7)
            self.assertRelativeEqual(
                            parm.parm_data['LENNARD_JONES_BCOEF'][idx-1],
                            gascham.parm_data['LENNARD_JONES_BCOEF'][idx2-1],
                            places=7)
        # Ensure that the rest of the values are unchanged (exactly equal)
        for i in range(ntypes2):
            for j in range(ntypes2):
                idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j]
                idx2 = gascham.parm_data['NONBONDED_PARM_INDEX'][ntypes2*i+j]
                self.assertEqual(
                            parm.parm_data['LENNARD_JONES_ACOEF'][idx-1],
                            gascham.parm_data['LENNARD_JONES_ACOEF'][idx2-1])
                self.assertEqual(
                            parm.parm_data['LENNARD_JONES_BCOEF'][idx-1],
                            gascham.parm_data['LENNARD_JONES_BCOEF'][idx2-1])
        # Now supply keywords
        parm2 = copy(gascham)
        PT.addLJType(parm2, '@1', radius=1.0, epsilon=1.0).execute()
        PT.changeLJSingleType(parm, '@1', 1.0, 1.0).execute()
        for x, y in zip(parm.parm_data['LENNARD_JONES_ACOEF'],
                        parm2.parm_data['LENNARD_JONES_ACOEF']):
            self.assertRelativeEqual(x, y)
        for x, y in zip(parm.parm_data['LENNARD_JONES_BCOEF'],
                        parm2.parm_data['LENNARD_JONES_BCOEF']):
            self.assertRelativeEqual(x, y)
        # Now use addLJType to hack a way to turn off LJ interactions
        PT.addLJType(parm, '*', radius=0.0, epsilon=0.0).execute()
        ntypes = parm.ptr('ntypes')
        for atom in parm.atoms:
            self.assertEqual(atom.nb_idx, ntypes)
        idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*(ntypes-1)+ntypes-1]
        self.assertEqual(parm.parm_data['LENNARD_JONES_ACOEF'][idx-1], 0.0)
        self.assertEqual(parm.parm_data['LENNARD_JONES_BCOEF'][idx-1], 0.0)

    def testPrintLJTypes(self):
        """ Test printLJTypes for ChamberParm """
        # Simple test
        act = PT.printLJTypes(gascham, '@1')
        for line in str(act).split('\n'):
            if not line.startswith('ATOM'):
                continue
            words = line.split()
            self.assertTrue(words[2].startswith('N'))
            self.assertTrue(words[3].startswith('N'))
            self.assertEqual(words[7], '1')

    def testSceeScnb(self):
        """ Test scee and scnb for ChamberParm """
        parm = copy(gascham)
        PT.scee(parm, 10).execute()
        PT.scnb(parm, 10).execute()
        for dih in parm.dihedrals:
            self.assertEqual(dih.type.scee, 10)
            self.assertEqual(dih.type.scnb, 10)
        for x, y in zip(parm.parm_data['SCEE_SCALE_FACTOR'],
                        parm.parm_data['SCNB_SCALE_FACTOR']):
            self.assertEqual(x, 10)
            self.assertEqual(y, 10)

    def testPrintDetails(self):
        """ Test printDetails for ChamberParm """
        act = PT.printDetails(gascham, '@1')
        self.assertEqual(str(act), saved.PRINT_DETAILSC)

    def testPrintFlags(self):
        """ Test printFlags for ChamberParm """
        act = PT.printFlags(gascham)
        printed_flags = set()
        for line in str(act).split('\n'):
            if line.startswith('%FLAG'):
                printed_flags.add(line.split()[1])
        self.assertEqual(printed_flags, set(gascham.parm_data.keys()))

    def testPrintPointers(self):
        """ Test printPointers for ChamberParm """
        act = PT.printPointers(gascham)
        printed_pointers = set(['NEXT', 'NIMPRTYPES', 'NUBTYPES', 'CMAP',
                                'CMAP_TYPES', 'NIMPHI', 'NUB'])
        for line in str(act).split('\n'):
            try:
                pointer = line.split()[0]
                value = int(line[line.rfind('=')+1:].strip())
            except (IndexError, ValueError):
                continue
            self.assertEqual(gascham.ptr(pointer), value)
            printed_pointers.add(pointer)
        self.assertEqual(printed_pointers, set(gascham.pointers.keys()))

    def testPrintBonds(self):
        """ Test printBonds for ChamberParm """
        act = PT.printBonds(gascham, '@1')
        self.assertEqual(str(act), saved.PRINT_BONDSC)

    def testPrintAngles(self):
        """ Test printAngles for ChamberParm """
        act = PT.printAngles(gascham, '@1')
        self.assertEqual(str(act), saved.PRINT_ANGLESC)

    def testPrintDihedrals(self):
        """ Test printDihedrals for ChamberParm """
        act = PT.printDihedrals(gascham, '@1')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALSC)

    def testSetMolecules(self):
        """ Test setMolecules for ChamberParm """
        parm = copy(solvchamber)
        atoms = [atom for atom in parm.atoms] # shallow copy!
        self.assertTrue(all([x is y for x,y in zip(parm.atoms,atoms)]))
        self.assertEqual(parm.ptr('IPTRES'), 160)
        self.assertEqual(parm.ptr('NSPM'), 17857)
        self.assertEqual(parm.ptr('NSPSOL'), 2)
        # To keep the output clean
        PT.setMolecules(parm).execute()
        self.assertTrue(all([x is y for x,y in zip(parm.atoms, atoms)]))
        # Now check that setMolecules can apply another time. solute_ions seems
        # to be broken, and I can't figure out why.
        PT.setMolecules(parm).execute()

    def testNetCharge(self):
        """ Test netCharge for ChamberParm """
        act = PT.netCharge(gascham)
        chg = act.execute() # check this part of the API
        self.assertEqual(str(act), 'The net charge of :* is %.4f' % chg)
        self.assertAlmostEqual(chg, 0.0, places=6)
        chg = PT.netCharge(gascham, ':1').execute()
        self.assertAlmostEqual(chg, 1.0, places=6)
        chg = PT.netCharge(gascham, ':3').execute()
        self.assertAlmostEqual(chg, -1.0, places=6)

    def testStrip(self):
        """ Test strip action for ChamberParm """
        parm = copy(gascham)
        PT.strip(parm, ':1').execute()
        self.assertEqual(parm.ptr('natom'), 21)
        self.assertEqual(len(parm.atoms), 21)
        # Good enough for here. The strip action is repeatedly tested in the
        # core Amber test suite as part of the MM/PBSA tests via ante-MMPBSA.py
        # and that part also tests that the energies come out correct as well

    def testDefineSolvent(self):
        """ Test defineSolvent for ChamberParm """
        import parmed.residue as residue
        PT.defineSolvent(gascham, 'WAT,HOH,Na+,Cl-').execute()
        self.assertEqual(residue.SOLVENT_NAMES, 'WAT HOH Na+ Cl-'.split())
        PT.defineSolvent(gascham, 'WAT,HOH').execute()
        self.assertEqual(residue.SOLVENT_NAMES, 'WAT HOH'.split())

    def testAddExclusions(self):
        """ Test addExclusions for ChamberParm """
        parm = copy(gascham)
        in_exclusions_before = []
        for atom1 in parm.residues[0].atoms:
            all_exclusions = (atom1.bond_partners + atom1.angle_partners + 
                             atom1.dihedral_partners + atom1.exclusion_partners)
            for atom2 in parm.residues[0].atoms:
                if atom1 is atom2: continue
                in_exclusions_before.append(atom2 in all_exclusions)
        self.assertFalse(all(in_exclusions_before))
        PT.addExclusions(parm, ':1', ':1').execute()
        in_exclusions_after = []
        for atom1 in parm.residues[0].atoms:
            all_exclusions = (atom1.bond_partners + atom1.angle_partners + 
                             atom1.dihedral_partners + atom1.exclusion_partners)
            for atom2 in parm.residues[0].atoms:
                if atom1 is atom2: continue
                in_exclusions_after.append(atom2 in all_exclusions)
                if not in_exclusions_after[-1]:
                    print('%s %s not excluded' % (atom1, atom2))
        self.assertTrue(all(in_exclusions_after))

    def testAddDeleteDihedral(self):
        """ Test the addDihedral and deleteDihedral actions for ChamberParm """
        parm = copy(gascham)
        n = PT.deleteDihedral(parm, ':ALA@N :ALA@CA :ALA@CB :ALA@HB1').execute()
        parm.remake_parm()
        self.assertEqual(gascham.ptr('nphih') + gascham.ptr('nphia'),
                         parm.ptr('nphih') + parm.ptr('nphia') + n)
        NALA = sum([res.name == 'ALA' for res in parm.residues])
        self.assertEqual(n, NALA)
        PT.addDihedral(parm, ':ALA@N', ':ALA@CA', ':ALA@CB', ':ALA@HB1',
                       0.1556, 3, 0, 1.2, 2.0).execute()
        parm.remake_parm()
        self.assertEqual(gascham.ptr('nphih') + gascham.ptr('nphia'),
                         parm.ptr('nphih') + parm.ptr('nphia'))
        PT.addDihedral(parm, ':ALA@N', ':ALA@CA', ':ALA@CB', ':ALA@HB1',
                       0.1556, 1, 0, 1.2, 2.0, type='normal').execute()
        parm.remake_parm()
        self.assertEqual(gascham.ptr('nphih') + gascham.ptr('nphia'),
                         parm.ptr('nphih') + parm.ptr('nphia') - n)
        num_dihedrals = 0
        num_ignore_ends = 0
        for atom in parm.atoms:
            if atom.residue.name == 'ALA' and atom.name == 'N':
                for dih in atom.dihedrals:
                    if dih.atom1 is atom:
                        if (dih.atom2.name == 'CA' and dih.atom3.name == 'CB'
                            and dih.atom4.name == 'HB1'):
                            num_dihedrals += 1
                            if dih.ignore_end:
                                num_ignore_ends += 1
                    elif dih.atom4 is atom:
                        if (dih.atom2.name == 'CB' and dih.atom3.name == 'CA'
                            and dih.atom1.name == 'HB1'):
                            num_dihedrals += 1
                            if dih.ignore_end:
                                num_ignore_ends += 1
        self.assertEqual(num_dihedrals, 2*NALA)
        self.assertEqual(num_ignore_ends, 1*NALA)

    def testSetBond(self):
        """ Test setBond for ChamberParm """
        parm = copy(gascham)
        PT.setBond(parm, ':ALA@CA', ':ALA@CB', 300.0, 1.5).execute()
        act = PT.printBonds(parm, ':ALA@CA')
        self.assertEqual(str(act), saved.SET_BONDC)

    def testSetAngle(self):
        """ Test setAngle for ChamberParm """
        parm = copy(gascham)
        PT.setAngle(parm, ':ALA@CA', ':ALA@CB', ':ALA@HB1', 40, 100).execute()
        act = PT.printAngles(parm, ':ALA@CB')
        self.assertEqual(str(act), saved.SET_ANGLEC)

    def testAddAtomicNumber(self):
        """ Test addAtomicNumber for ChamberParm """
        parm = copy(gascham)
        self.assertFalse('ATOMIC_NUMBER' in parm.parm_data)
        atomic_numbers = [atom.atomic_number for atom in parm.atoms]
        PT.addAtomicNumber(parm).execute()
        self.assertEqual(parm.parm_data['ATOMIC_NUMBER'], atomic_numbers)

    def testPrintLJMatrix(self):
        """ Test printLJMatrix for ChamberParm """
        act = PT.printLJMatrix(gascham, '@1')
        self.assertEqual(str(act), saved.PRINT_LJMATRIXC)

    def testDeleteBond(self):
        """ Test deleteBond for ChamberParm """
        parm = copy(gascham)
        # Pick the bond we plan to delete, pick out every angle and dihedral
        # that contains that bond, and then delete it. Then make sure none of
        # the valence terms that contained that bond remain afterwards. We
        # already have a test to make sure that the __contains__ method works
        # for atoms and bonds.
        for bond in parm.atoms[10].bonds:
            if parm.atoms[12] in bond: break
        deleted_angles = list()
        deleted_dihedrals = list()
        deleted_impropers = list()
        deleted_urey_bradleys = list()
        deleted_cmaps = list()
        for angle in parm.angles:
            if bond in angle: deleted_angles.append(angle)
        for dihedral in parm.dihedrals:
            if bond in dihedral: deleted_dihedrals.append(dihedral)
        for imp in parm.impropers:
            if bond in imp: deleted_impropers.append(imp)
        for ub in parm.urey_bradleys:
            if bond in ub: deleted_urey_bradleys.append(ub)
        for cmap in parm.cmaps:
            if bond in cmap: deleted_cmaps.append(cmap)
        act = PT.deleteBond(parm, '@11', '@13', 'verbose')
        act.execute()
        self.assertTrue(bond not in parm.bonds)
        for angle in deleted_angles:
            self.assertTrue(angle not in parm.angles)
        for dihedral in deleted_dihedrals:
            self.assertTrue(all([dihedral is not d for d in parm.dihedrals]))
        for improper in deleted_impropers:
            self.assertTrue(improper not in parm.impropers)
        for ureybrad in deleted_urey_bradleys:
            self.assertTrue(ureybrad not in parm.urey_bradleys)
        self.assertFalse(parm.has_cmap)

    def testSummary(self):
        """ Test summary action for ChamberParm """
        parm = copy(solvchamber)
        parm.load_rst7(get_fn('dhfr_cmap_pbc.rst7'))
        act = PT.summary(parm)
        self.assertTrue(utils.detailed_diff(str(act), saved.SUMMARYC1,
                                            relative_error=1e-6))

    def testScale(self):
        """ Test scale action for ChamberParm """
        parm = copy(gascham)
        PT.scale(parm, 'DIHEDRAL_FORCE_CONSTANT', 2.0).execute()
        self.assertEqual(
                [2*x for x in gascham.parm_data['DIHEDRAL_FORCE_CONSTANT']],
                parm.parm_data['DIHEDRAL_FORCE_CONSTANT'])
        for dt1, dt2 in zip(parm.dihedral_types, gascham.dihedral_types):
            self.assertEqual(dt1.phi_k, dt2.phi_k*2)
        PT.scale(parm, 'DIHEDRAL_FORCE_CONSTANT', 0.0).execute()
        for val in parm.parm_data['DIHEDRAL_FORCE_CONSTANT']:
            self.assertEqual(val, 0)
        for dt in parm.dihedral_types:
            self.assertEqual(dt.phi_k, 0)

    def testInterpolate(self):
        """ Test interpolate action for ChamberParm """
        self._empty_writes()
        parm = copy(gascham)
        origparm = copy(gascham)
        origparm.name = origparm.name + '_copy1'
        PT.change(parm, 'CHARGE', ':1', 0).execute()
        self.assertAlmostEqual(sum(parm.parm_data['CHARGE']), -1)
        self.assertAlmostEqual(sum(origparm.parm_data['CHARGE']), 0)
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(atom.charge, parm.parm_data['CHARGE'][i])
        for i, atom in enumerate(origparm.atoms):
            self.assertEqual(atom.charge, origparm.parm_data['CHARGE'][i])
        # Now set up a ParmList so we can interpolate these topology files
        parms = parmlist.ParmList()
        parms.add_parm(parm)
        parms.add_parm(origparm)
        sys.stdout = open(os.devnull, 'w')
        PT.interpolate(parms, 5, 'eleconly', startnum=2,
                       prefix=get_fn('test.parm7', written=True)).execute()
        sys.stdout = sys.__stdout__
        self.assertEqual(len(os.listdir(get_fn('writes'))), 5)
        self.assertTrue(os.path.exists(get_fn('test.parm7.2', written=True)))
        self.assertTrue(os.path.exists(get_fn('test.parm7.3', written=True)))
        self.assertTrue(os.path.exists(get_fn('test.parm7.4', written=True)))
        self.assertTrue(os.path.exists(get_fn('test.parm7.5', written=True)))
        self.assertTrue(os.path.exists(get_fn('test.parm7.6', written=True)))
        # Now check them all
        ladder = [origparm]
        ladder.append(AmberParm(get_fn('test.parm7.2', written=True)))
        ladder.append(AmberParm(get_fn('test.parm7.3', written=True)))
        ladder.append(AmberParm(get_fn('test.parm7.4', written=True)))
        ladder.append(AmberParm(get_fn('test.parm7.5', written=True)))
        ladder.append(AmberParm(get_fn('test.parm7.6', written=True)))
        ladder.append(parm)
        natom = parm.ptr('natom')
        for i in range(1, 6):
            before = ladder[i-1].parm_data['CHARGE']
            after = ladder[i+1].parm_data['CHARGE']
            this = ladder[i].parm_data['CHARGE']
            for j in range(natom):
                if before[j] < after[j]:
                    self.assertTrue(before[j] <= this[j] <= after[j])
                else:
                    self.assertTrue(before[j] >= this[j] >= after[j])

    def testChangeProtState(self):
        """ Check that changeProtState fails for ChamberParm """
        parm = copy(solvchamber)
        self.assertRaises(exc.ParmError, lambda:
                PT.changeProtState(parm, ':32', 0).execute())

    def testLmod(self):
        """ Test lmod action for ChamberParm """
        parm = copy(gascham)
        parm.parm_data['LENNARD_JONES_ACOEF'][3] = 0.0
        self.assertFalse(all(parm.parm_data['LENNARD_JONES_ACOEF']))
        PT.lmod(parm).execute()
        self.assertTrue(all(parm.parm_data['LENNARD_JONES_ACOEF']))

    def testAddDeletePDB(self):
        """ Test addPDB and deletePDB actions for ChamberParm """
        parm = copy(gascham)
        PT.addPDB(parm, get_fn('ala_ala_ala.pdb'), 'elem allicodes').execute()
        self.assertTrue('RESIDUE_ICODE' in parm.flag_list)
        self.assertTrue('ATOM_ELEMENT' in parm.flag_list)
        self.assertTrue('RESIDUE_NUMBER' in parm.flag_list)
        self.assertTrue('RESIDUE_CHAINID' in parm.flag_list)
        self.assertTrue(len(parm.parm_data['RESIDUE_ICODE']), parm.ptr('nres'))
        self.assertTrue(len(parm.parm_data['ATOM_ELEMENT']), parm.ptr('natom'))
        self.assertTrue(len(parm.parm_data['RESIDUE_NUMBER']), parm.ptr('nres'))
        self.assertTrue(len(parm.parm_data['RESIDUE_CHAINID']),parm.ptr('nres'))
        for i in range(parm.ptr('nres')):
            self.assertEqual(parm.parm_data['RESIDUE_NUMBER'][i], i+1)
            self.assertEqual(parm.parm_data['RESIDUE_ICODE'][i], '')
            self.assertEqual(parm.parm_data['RESIDUE_CHAINID'][i], 'A')
        for i, atom in enumerate(parm.atoms):
            atnum = atom.atomic_number
            elem = parm.parm_data['ATOM_ELEMENT'][i].strip()
            self.assertEqual(periodic_table.Element[atnum], elem)
            self.assertEqual(atnum, periodic_table.AtomicNum[elem])
        PT.deletePDB(parm).execute()
        self.assertFalse('RESIDUE_ICODE' in parm.flag_list)
        self.assertFalse('ATOM_ELEMENT' in parm.flag_list)
        self.assertFalse('RESIDUE_NUMBER' in parm.flag_list)
        self.assertFalse('RESIDUE_CHAINID' in parm.flag_list)
        PT.addPDB(parm, get_fn('ala_ala_ala.pdb')).execute()
        self.assertFalse('RESIDUE_ICODE' in parm.flag_list)
        self.assertFalse('ATOM_ELEMENT' in parm.flag_list)
        self.assertTrue('RESIDUE_NUMBER' in parm.flag_list)
        self.assertTrue('RESIDUE_CHAINID' in parm.flag_list)

    def testHMassRepartition(self):
        """ Test HMassRepartition action for ChamberParm """
        parm = copy(solvchamber)
        PT.defineSolvent(parm, 'TIP3').execute()
        PT.HMassRepartition(parm, 2.0).execute()
        for atom in parm.atoms:
            if atom.atomic_number == 1:
                if atom.residue.name == 'TIP3':
                    self.assertAlmostEqual(atom.mass, 1.008)
                else:
                    self.assertEqual(atom.mass, 2)
        self.assertEqual(parm.parm_data['MASS'],
                         [a.mass for a in parm.atoms])
        self.assertAlmostEqual(sum(solvchamber.parm_data['MASS']),
                               sum(parm.parm_data['MASS']))
        PT.HMassRepartition(parm, 3.0, 'dowater').execute()
        for atom in parm.atoms:
            if atom.atomic_number == 1:
                self.assertEqual(atom.mass, 3.0)
        self.assertAlmostEqual(sum(solvchamber.parm_data['MASS']),
                               sum(parm.parm_data['MASS']), places=6)

    def testOutPDB(self):
        """ Test the outPDB action on ChamberParm """
        parm = copy(gascham)
        PT.loadRestrt(parm, get_fn('ala_ala_ala.rst7')).execute()
        PT.outPDB(parm, get_fn('outPDB1.pdb', written=True)).execute()
        f = PDBFile.parse(get_fn('outPDB1.pdb', written=True))
        self.assertEqual(len(f.atoms), len(parm.atoms))
        self.assertEqual(len(f.residues), len(parm.residues))
        for a1, a2 in zip(f.atoms, parm.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertAlmostEqual(a1.xx, a2.xx, delta=2e-3)
            self.assertAlmostEqual(a1.xy, a2.xy, delta=2e-3)
            self.assertAlmostEqual(a1.xz, a2.xz, delta=2e-3)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.idx, a2.residue.idx)

    def testOutCIF(self):
        """ Test the outCIF action on ChamberParm """
        parm = copy(gascham)
        PT.loadRestrt(parm, get_fn('ala_ala_ala.rst7')).execute()
        PT.outCIF(parm, get_fn('outPDB1.cif', written=True)).execute()
        f = CIFFile.parse(get_fn('outPDB1.cif', written=True))
        self.assertEqual(len(f.atoms), len(parm.atoms))
        self.assertEqual(len(f.residues), len(parm.residues))
        for a1, a2 in zip(f.atoms, parm.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertAlmostEqual(a1.xx, a2.xx, delta=2e-3)
            self.assertAlmostEqual(a1.xy, a2.xy, delta=2e-3)
            self.assertAlmostEqual(a1.xz, a2.xz, delta=2e-3)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.idx, a2.residue.idx)

class TestAmoebaParmActions(utils.TestCaseRelative, utils.FileIOTestCase):
    """ Tests actions on Amber prmtop files """
    
    def testParmoutOutparmLoadRestrt(self):
        """ Test parmout, outparm, and loadRestrt actions on AmoebaParm """
        self._empty_writes()
        parm = copy(amoebaparm)
        PT.loadRestrt(parm, get_fn('nma.rst7')).execute()
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))
        PT.parmout(parm, get_fn('test.parm7', written=True)).execute()
        self.assertEqual(len(os.listdir(get_fn('writes'))), 1)
        self.assertTrue(diff_files(get_fn('nma.parm7'),
                                   get_fn('test.parm7', written=True)))
        self._empty_writes()
        PT.parmout(parm, get_fn('test.parm7', written=True),
                         get_fn('test.rst7', written=True)).execute()
        self.assertEqual(len(os.listdir(get_fn('writes'))), 2)
        self.assertTrue(diff_files(get_fn('nma.parm7'),
                                   get_fn('test.parm7', written=True)))
        self.assertTrue(diff_files(get_fn('nma.rst7'),
                                   get_fn('test.rst7', written=True),
                                   absolute_error=0.0001))
        self._empty_writes()
        PT.outparm(parm, get_fn('test.parm7', written=True)).execute()
        self.assertEqual(len(os.listdir(get_fn('writes'))), 1)
        self.assertTrue(diff_files(get_fn('nma.parm7'),
                                   get_fn('test.parm7', written=True)))
        self._empty_writes()
        PT.outparm(parm, get_fn('test.parm7', written=True),
                         get_fn('test.rst7', written=True)).execute()
        self.assertEqual(len(os.listdir(get_fn('writes'))), 2)
        self.assertTrue(diff_files(get_fn('nma.parm7'),
                                   get_fn('test.parm7', written=True)))
        self.assertTrue(diff_files(get_fn('nma.rst7'),
                                   get_fn('test.rst7', written=True),
                                   absolute_error=0.0001))

    def testWriteFrcmod(self):
        """ Check that writeFrcmod fails for AmoebaParm """
        self.assertRaises(exc.ParmError, lambda:
                PT.writeFrcmod(amoebaparm, get_fn('x', written=True)).execute())

    def testWriteOffLoadRestrt(self):
        """ Check that writeOFF fails for AmoebaParm """
        parm = copy(amoebaparm)
        PT.loadRestrt(parm, get_fn('nma.rst7')).execute()
        self.assertRaises(exc.ParmError, lambda:
                PT.writeOFF(parm, get_fn('test.off', written=True)).execute())

    def testChangeRadii(self):
        """ Check that changeRadii fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.changeRadii(parm, 'amber6').execute())

    def testChangeLJPair(self):
        """ Check that changeLJPair fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.changeLJPair(parm, '@%NH3', '@%HC', 1.0, 1.0).execute())

    def testChangeLJ14Pair(self):
        """ Check that changeLJ14Pair fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.changeLJ14Pair(parm, '@%NH3', '@%HC', 1.0, 1.0).execute())

    def testChange(self):
        """ Test the 'change' action for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmedChangeError, lambda:
                PT.change(parm, 'CHARGE', ':ACE', 0, 'quiet').execute())
        PT.change(parm, 'MASS', ':1', 10.0).execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(parm.parm_data['MASS'][i], atom.mass)
            if atom.residue.idx == 0:
                self.assertEqual(atom.mass, 10.0)
            else:
                self.assertEqual(atom.mass, amoebaparm.parm_data['MASS'][i])
        PT.change(parm, 'ATOM_NAME', ':WAT@O', 'JMS').execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(parm.parm_data['ATOM_NAME'][i], atom.name)
            if atom.residue.name == 'WAT' and \
                        amoebaparm.parm_data['ATOM_NAME'][i] == 'O':
                self.assertEqual(atom.name, 'JMS')
            else:
                self.assertEqual(atom.name,
                                 amoebaparm.parm_data['ATOM_NAME'][i])
        PT.change(parm, 'AMBER_ATOM_TYPE', ':WAT@H*', 'RJLS').execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(parm.parm_data['AMBER_ATOM_TYPE'][i], atom.type)
            if atom.residue.name == 'WAT' and \
                        amoebaparm.parm_data['AMBER_ATOM_TYPE'][i][0] == 'H':
                self.assertEqual(atom.type, 'RJLS')
            else:
                self.assertEqual(atom.type,
                                 amoebaparm.parm_data['AMBER_ATOM_TYPE'][i])
        self.assertRaises(exc.ParmedChangeError, lambda:
                PT.change(parm, 'ATOM_TYPE_INDEX', '@1', 4).execute())
        self.assertRaises(exc.ParmedChangeError, lambda:
                PT.change(parm, 'RADII', ':1-2', 2.0, 'quiet').execute())
        self.assertRaises(exc.ParmedChangeError, lambda:
                PT.change(parm, 'SCREEN', '*', 0.0).execute())
        # Check bad input
        self.assertRaises(exc.ParmedChangeError, lambda:
                          PT.change(parm, 'RESIDUE_LABEL', ':*', 'NaN'))

    def testPrintInfo(self):
        """ Test printInfo for all flags of AmoebaParm """
        for flag in amoebaparm.parm_data:
            if flag == 'TITLE': continue
            act = PT.printInfo(amoebaparm, flag)
            vals = []
            for line in str(act).split('\n'):
                vals += line.split()
            self.assertEqual(len(vals), len(amoebaparm.parm_data[flag]))
            try:
                datatype = type(amoebaparm.parm_data[flag][0])
            except IndexError:
                continue
            for i, j in zip(vals, amoebaparm.parm_data[flag]):
                # printInfo prints to 5 places for floats.
                if datatype is float:
                    self.assertAlmostEqual(datatype(i), j, places=4)
                else:
                    self.assertEqual(datatype(i), j)

    def testAddChangeLJType(self):
        """ Check that addLJType and changeLJSingleType fail for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.addLJType(parm, '@1').execute())
        self.assertRaises(exc.ParmError, lambda:
                PT.changeLJSingleType(parm, '@1', 1.0, 1.0).execute())

    def testPrintLJTypes(self):
        """ Check that printLJTypes fails for AmoebaParm """
        self.assertRaises(exc.ParmError, lambda:
                PT.printLJTypes(amoebaparm, '@1'))

    def testSceeScnb(self):
        """ Check that scee and scnb fail for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda: PT.scee(parm, 10).execute())
        self.assertRaises(exc.ParmError, lambda: PT.scnb(parm, 10).execute())

    def testPrintDetails(self):
        """ Test printDetails for AmoebaParm """
        act = PT.printDetails(amoebaparm, ':1-2')
        self.assertEqual(str(act), saved.PRINT_DETAILSA)

    def testPrintFlags(self):
        """ Test printFlags for AmoebaParm """
        act = PT.printFlags(amoebaparm)
        printed_flags = set()
        for line in str(act).split('\n'):
            if line.startswith('%FLAG'):
                printed_flags.add(line.split()[1])
        self.assertEqual(printed_flags, set(amoebaparm.parm_data.keys()))

    def testPrintPointers(self):
        """ Test printPointers for AmoebaParm """
        act = PT.printPointers(amoebaparm)
        printed_pointers = set(['NEXT'])
        for line in str(act).split('\n'):
            try:
                pointer = line.split()[0]
                value = int(line[line.rfind('=')+1:].strip())
            except (IndexError, ValueError):
                continue
            self.assertEqual(amoebaparm.ptr(pointer), value)
            printed_pointers.add(pointer)
        self.assertEqual(printed_pointers, set(amoebaparm.pointers.keys()))

    def testPrintBonds(self):
        """ Test printBonds for AmoebaParm """
        act = PT.printBonds(amoebaparm, '@1')
        self.assertEqual(str(act), saved.PRINT_BONDSA)

    def testPrintAngles(self):
        """ Test printAngles for AmoebaParm """
        act = PT.printAngles(amoebaparm, '@1')
        self.assertEqual(str(act), saved.PRINT_ANGLESA)

    def testPrintDihedrals(self):
        """ Test printDihedrals for AmoebaParm """
        act = PT.printDihedrals(amoebaparm, '@1')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALSA)

    def testSetMolecules(self):
        """ Test setMolecules for AmoebaParm """
        parm = copy(amoebaparm)
        atoms = [atom for atom in parm.atoms] # shallow copy!
        self.assertTrue(all([x is y for x,y in zip(parm.atoms,atoms)]))
        self.assertEqual(parm.ptr('IPTRES'), 2)
        self.assertEqual(parm.ptr('NSPM'), 819)
        self.assertEqual(parm.ptr('NSPSOL'), 2)
        PT.setMolecules(parm).execute()
        self.assertEqual(parm.ptr('IPTRES'), 2)
        self.assertEqual(parm.ptr('NSPM'), 819)
        self.assertEqual(parm.ptr('NSPSOL'), 2)
        self.assertTrue(all([x is y for x,y in zip(parm.atoms, atoms)]))
        # Now check that setMolecules can apply another time
        PT.setMolecules(parm).execute()

    def testNetCharge(self):
        """ Test netCharge for AmoebaParm (charge is the monopole) """
        act = PT.netCharge(amoebaparm)
        chg = act.execute() # check the netCharge.execute return value
        self.assertEqual(str(act), 'The net charge of :* is %.4f' % chg)
        self.assertAlmostEqual(chg, 0.0)
        chg = PT.netCharge(amoebaparm, ':WAT').execute()
        self.assertAlmostEqual(chg, 0)

    def testStrip(self):
        """ Test strip action for AmoebaParm """
        parm = copy(amoebaparm)
        natoms = len(parm.atoms)
        lenres = len(parm.residues[0])
        PT.strip(parm, ':1').execute()
        self.assertEqual(parm.ptr('natom'), natoms-lenres)
        self.assertEqual(len(parm.atoms), natoms-lenres)
        # Good enough for here. The strip action is repeatedly tested in the
        # core Amber test suite as part of the MM/PBSA tests via ante-MMPBSA.py
        # and that part also tests that the energies come out correct as well

    def testDefineSolvent(self):
        """ Test defineSolvent for AmoebaParm """
        import parmed.residue as residue
        PT.defineSolvent(amoebaparm, 'WAT,HOH,Na+,Cl-').execute()
        self.assertEqual(residue.SOLVENT_NAMES, 'WAT HOH Na+ Cl-'.split())
        PT.defineSolvent(amoebaparm, 'WAT,HOH').execute()
        self.assertEqual(residue.SOLVENT_NAMES, 'WAT HOH'.split())

    def testAddExclusions(self):
        """ Check that addExclusions fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.addExclusions(parm, ':*', ':*').execute())

    def testAddDeleteDihedral(self):
        """ Check that addDihedral and deleteDihedral fail for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.deleteDihedral(parm, '@1 @2 @3 @4').execute())
        self.assertRaises(exc.ParmError, lambda:
                PT.addDihedral(parm, '@1 @2 @3 @4 0.1556 3 0 1.2 2.0',
                               type='normal').execute()
        )

    def testSetBond(self):
        """ Check that setBond fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.setBond(parm, ':ALA@CA', ':ALA@CB', 300.0, 1.5).execute())

    def testSetAngle(self):
        """ Check that setAngle fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.setAngle(parm, ':ALA@CA :ALA@CB :ALA@HB1 40 100').execute())

    def testAddAtomicNumber(self):
        """ Test addAtomicNumber for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertFalse('ATOMIC_NUMBER' in parm.parm_data)
        atomic_numbers = [atom.atomic_number for atom in parm.atoms]
        PT.addAtomicNumber(parm).execute()
        self.assertEqual(parm.parm_data['ATOMIC_NUMBER'], atomic_numbers)

    def testPrintLJMatrix(self):
        """ Check that printLJMatrix fails for AmoebaParm """
        self.assertRaises(exc.ParmError, lambda:
                PT.printLJMatrix(amoebaparm, '@1'))

    def testDeleteBond(self):
        """ Test deleteBond for AmoebaParm """
        parm = copy(amoebaparm)
        for bond in parm.atoms[0].bonds:
            if parm.atoms[1] in bond: break
        TrackedList = type(parm.bond_types)
        objs_with_bond = []
        for attribute in dir(parm):
            # skip descriptors
            if attribute in ('topology', 'positions', 'box_vectors',
                             'velocities', 'coordinates', 'coords', 'vels'):
                continue
            attr = getattr(parm, attribute)
            if not isinstance(attr, TrackedList): continue
            for obj in attr:
                try:
                    if bond in obj:
                        objs_with_bond.append(attr)
                        break
                except TypeError:
                    break
        self.assertTrue(len(objs_with_bond) > 0)
        act = PT.deleteBond(parm, '@1', '@2', 'verbose')
        act.execute()
        self.assertTrue(bond not in parm.bonds)
        for attr in objs_with_bond:
            for obj in attr:
                self.assertNotIn(bond, attr)

    def testSummary(self):
        """ Test summary action for AmoebaParm """
        parm = copy(amoebaparm)
        act = PT.summary(parm)
        self.assertEqual(str(act), saved.SUMMARYA1)
        PT.loadRestrt(parm, get_fn('nma.rst7')).execute()
        act = PT.summary(parm)
        self.assertEqual(str(act), saved.SUMMARYA2)

    def testScale(self):
        """ Test scale action for AmoebaParm """
        parm = copy(amoebaparm)
        flag = 'AMOEBA_STRETCH_BEND_FORCE_CONSTANT'
        PT.scale(parm, flag, 2.0).execute()
        self.assertEqual([2*x for x in amoebaparm.parm_data[flag]],
                         parm.parm_data[flag])
        PT.scale(parm, flag, 0.0).execute()
        for val in parm.parm_data[flag]:
            self.assertEqual(val, 0)

    def testInterpolate(self):
        """ Check that interpolate action fails for AmoebaParm """
        self.assertRaises(exc.ParmError, lambda:
                PT.interpolate(amoebaparm).execute())

    def testChangeProtState(self):
        """ Check that changeProtState fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.changeProtState(parm, ':32', 0).execute())

    def testLmod(self):
        """ Check that lmod fails for AmoebaParm """
        self.assertRaises(exc.ParmError, lambda: PT.lmod(amoebaparm).execute())

    def testAddDeletePDB(self):
        """ Test addPDB and deletePDB for AmoebaParm """
        parm = copy(amoebaparm)
        PT.addPDB(parm, get_fn('nma.pdb'), 'elem allicodes').execute()
        self.assertTrue('RESIDUE_ICODE' in parm.flag_list)
        self.assertTrue('ATOM_ELEMENT' in parm.flag_list)
        self.assertTrue('RESIDUE_NUMBER' in parm.flag_list)
        self.assertTrue('RESIDUE_CHAINID' in parm.flag_list)
        self.assertTrue(len(parm.parm_data['RESIDUE_ICODE']), parm.ptr('nres'))
        self.assertTrue(len(parm.parm_data['ATOM_ELEMENT']), parm.ptr('natom'))
        self.assertTrue(len(parm.parm_data['RESIDUE_NUMBER']), parm.ptr('nres'))
        self.assertTrue(len(parm.parm_data['RESIDUE_CHAINID']),parm.ptr('nres'))
        for i in range(parm.ptr('nres')):
            self.assertEqual(parm.parm_data['RESIDUE_NUMBER'][i], i+1)
            self.assertEqual(parm.parm_data['RESIDUE_ICODE'][i], '')
            if parm.residues[i].name == 'WAT':
                self.assertEqual(parm.parm_data['RESIDUE_CHAINID'][i], 'B')
            else:
                self.assertEqual(parm.parm_data['RESIDUE_CHAINID'][i], 'A')
        for i, atom in enumerate(parm.atoms):
            atnum = atom.atomic_number
            elem = parm.parm_data['ATOM_ELEMENT'][i].strip()
            self.assertEqual(periodic_table.Element[atnum], elem)
            self.assertEqual(atnum, periodic_table.AtomicNum[elem])
        PT.deletePDB(parm).execute()
        self.assertFalse('RESIDUE_ICODE' in parm.flag_list)
        self.assertFalse('ATOM_ELEMENT' in parm.flag_list)
        self.assertFalse('RESIDUE_NUMBER' in parm.flag_list)
        self.assertFalse('RESIDUE_CHAINID' in parm.flag_list)
        PT.addPDB(parm, get_fn('nma.pdb')).execute()
        self.assertFalse('RESIDUE_ICODE' in parm.flag_list)
        self.assertFalse('ATOM_ELEMENT' in parm.flag_list)
        self.assertTrue('RESIDUE_NUMBER' in parm.flag_list)
        self.assertTrue('RESIDUE_CHAINID' in parm.flag_list)

    def testHMassRepartition(self):
        """ Test HMassRepartition action for AmoebaParm """
        parm = copy(amoebaparm)
        PT.HMassRepartition(parm, 2.0).execute()
        for atom in parm.atoms:
            if atom.atomic_number == 1:
                if atom.residue.name == 'WAT':
                    self.assertAlmostEqual(atom.mass, 1.008)
                else:
                    self.assertEqual(atom.mass, 2)
        self.assertEqual(parm.parm_data['MASS'],
                         [a.mass for a in parm.atoms])
        self.assertAlmostEqual(sum(amoebaparm.parm_data['MASS']),
                               sum(parm.parm_data['MASS']))
        PT.HMassRepartition(parm, 3.0, 'dowater').execute()
        for atom in parm.atoms:
            if atom.atomic_number == 1:
                self.assertEqual(atom.mass, 3.0)
        self.assertAlmostEqual(sum(amoebaparm.parm_data['MASS']),
                               sum(parm.parm_data['MASS']), places=6)

    def testOutPDB(self):
        """ Test the outPDB action on AmoebaParm """
        parm = copy(amoebaparm)
        PT.loadRestrt(parm, get_fn('nma.rst7')).execute()
        PT.outPDB(parm, get_fn('outPDB1.pdb', written=True)).execute()
        f = PDBFile.parse(get_fn('outPDB1.pdb', written=True))
        self.assertEqual(len(f.atoms), len(parm.atoms))
        self.assertEqual(len(f.residues), len(parm.residues))
        for a1, a2 in zip(f.atoms, parm.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertAlmostEqual(a1.xx, a2.xx, delta=2e-3)
            self.assertAlmostEqual(a1.xy, a2.xy, delta=2e-3)
            self.assertAlmostEqual(a1.xz, a2.xz, delta=2e-3)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.idx, a2.residue.idx)

    def testOutCIF(self):
        """ Test the outCIF action on AmoebaParm """
        parm = copy(amoebaparm)
        PT.loadRestrt(parm, get_fn('nma.rst7')).execute()
        PT.outCIF(parm, get_fn('outPDB1.cif', written=True)).execute()
        f = CIFFile.parse(get_fn('outPDB1.cif', written=True))
        self.assertEqual(len(f.atoms), len(parm.atoms))
        self.assertEqual(len(f.residues), len(parm.residues))
        for a1, a2 in zip(f.atoms, parm.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertAlmostEqual(a1.xx, a2.xx, delta=2e-3)
            self.assertAlmostEqual(a1.xy, a2.xy, delta=2e-3)
            self.assertAlmostEqual(a1.xz, a2.xz, delta=2e-3)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.idx, a2.residue.idx)

    def testTIMerge(self):
        """ Check that tiMerge fails for AmoebaParm """
        parm = copy(amoebaparm)
        PT.loadRestrt(parm, get_fn('nma.rst7')).execute()
        self.assertRaises(exc.ParmError, lambda:
                PT.tiMerge(parm, ':1-3', ':4-6', ':2', ':5').execute())

class TestOtherParm(unittest.TestCase):
    """ Tests the use of other parms as the main parm """

    def testSummary(self):
        """ Tests the use of a PDB file with the summary action """
        parm = load_file(get_fn('4lzt.pdb'))
        self.assertEqual(str(PT.summary(parm)), saved.PDB_SUMMARY)

if __name__ == '__main__':
    unittest.main()
