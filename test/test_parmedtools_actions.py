"""
Tests for the various actions in ParmEd
"""
from __future__ import division, print_function

from copy import copy
import numpy as np
import os
import parmed as pmd
from parmed import periodic_table, gromacs, load_file, amber
from parmed.amber import AmberParm, ChamberParm, AmoebaParm, AmberFormat, AmberMask
from parmed.charmm import CharmmPsfFile
from parmed.exceptions import AmberWarning, CharmmWarning
from parmed.formats import PDBFile, CIFFile
from parmed.utils import PYPY
from parmed.utils.six.moves import range, zip, StringIO
from parmed.utils.six import string_types, iteritems
import parmed.unit as u
import parmed.tools as PT
from parmed.tools import exceptions as exc
from parmed.tools import parmlist
from parmed.tools.simulations import sanderapi
from parmed.tools.actions import ArgumentList
import re
import saved_outputs as saved
import sys
import unittest
from utils import (HAS_GROMACS, get_fn, get_saved_fn, diff_files,
        FileIOTestCase, TestCaseRelative, detailed_diff, has_openmm,
        create_random_structure, app)
import warnings
try:
    import pandas as pd
except ImportError:
    pd = None
try:
    import sander
except ImportError:
    sander = None

gasparm = AmberParm(get_fn('trx.prmtop'))
solvparm = AmberParm(get_fn('solv2.parm7'))
gascham = ChamberParm(get_fn('ala_ala_ala.parm7'))
solvchamber = ChamberParm(get_fn('ala3_solv.parm7'))
amoebaparm = AmoebaParm(get_fn('nma.parm7'))

# Make sure default overwrite is False
PT.Action.overwrite = False

# Make some new action subclasses
class NewActionNoUsage(PT.Action):
    needs_parm = False
    def init(self, arg_list):
        raise exc.NoArgument()

class NewActionWithUsage(PT.Action):
    needs_parm = False
    usage = '<some_argument>'
    def init(self, arg_list):
        raise exc.NoArgument()

class TestActionAPI(unittest.TestCase):
    """ Tests the Action API """

    def test_argument_types(self):
        """ Test argument type handling of Action subclass """
        self.assertRaises(TypeError, lambda: PT.actions.Action(10))

    def test_error_handling(self):
        """ Test error handling in Action initializer """
        with self.assertRaises(exc.ParmError):
            PT.changeRadii(None, 'mbondi2')
        self.assertFalse(NewActionNoUsage(None).valid)
        self.assertFalse(NewActionWithUsage(None).valid)
        with self.assertWarns(exc.UnhandledArgumentWarning):
            PT.changeRadii(gasparm, 'mbondi2', 'extra', 'arguments')
        str(NewActionNoUsage(None)) # Make sure it doesn't crash

class TestActionAPI(unittest.TestCase):
    """ Tests the Action API """

    def test_argument_types(self):
        """ Test argument type handling of Action subclass """
        self.assertRaises(TypeError, lambda: PT.actions.Action(10))

    def test_error_handling(self):
        """ Test error handling in Action initializer """
        self.assertRaises(exc.ParmError, lambda:
                PT.changeRadii(None, 'mbondi2'))

class TestNonParmActions(FileIOTestCase):
    """ Tests all actions that do not require a prmtop instance """

    def setUp(self):
        self.parm = gasparm
        FileIOTestCase.setUp(self)

    def test_overwrite(self):
        """ Test setting overwrite capabilities on ParmEd interpeter """
        a = PT.setOverwrite(self.parm, False)
        self.assertFalse(PT.Action.overwrite)
        a.execute()
        self.assertFalse(PT.Action.overwrite)
        self.assertEqual(str(a), 'Files are NOT overwritable')
        a = PT.setOverwrite(self.parm, True)
        a.execute()
        self.assertTrue(PT.Action.overwrite)
        self.assertEqual(str(a), 'Files are overwritable')
        PT.setOverwrite(None, False).execute() # Turn it back off
        with self.assertWarns(exc.SeriousParmWarning):
            PT.setOverwrite(None, 'badarg')

    def test_list_parms(self):
        """ Test listing of the prmtop files in the ParmEd interpreter """
        a = PT.listParms(self.parm)
        a.execute() # Should do nothing
        lines = str(a).split('\n')
        self.assertEqual(lines[0], 'Loaded topology files:')
        self.assertEqual(lines[1], '[0]\t%s (active)' % self.get_fn('trx.prmtop'))

    @unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without GROMACS")
    def test_chamber(self):
        """ Test the chamber action with a basic protein """
        a = PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'), '-top',
                       self.get_fn('top_all22_prot.inp'), '-param', self.get_fn('par_all22_prot.inp'))
        str(a)
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertTrue(parm.chamber)
        self.assertTrue(parm.has_cmap)
        self.assertEqual(parm.ptr('ifbox'), 0)
        # Error checking
        self.assertRaises(exc.FileDoesNotExist, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'),
                           '-top', 'foo', '-param', self.get_fn('par_all22_prot.inp'),
                           '-crd', self.get_fn('ala_ala_ala.pdb')).execute()
        )
        self.assertRaises(exc.InputError, lambda:
                PT.chamber(self.parm, '-top', self.get_fn('top_all22_prot.inp'),
                           '-param', self.get_fn('par_all22_prot.inp'), '-crd',
                           self.get_fn('ala_ala_ala.pdb')).execute()
        )
        self.assertRaises(exc.FileDoesNotExist, lambda:
                PT.chamber(self.parm, '-psf', 'foo', '-top', self.get_fn('top_all22_prot.inp'),
                           '-param', self.get_fn('par_all22_prot.inp'), '-crd',
                           self.get_fn('ala_ala_ala.pdb')).execute()
        )
        self.assertRaises(exc.FileDoesNotExist, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'),
                           '-top', self.get_fn('top_all22_prot.inp'),
                           '-param', 'foo', '-crd', self.get_fn('ala_ala_ala.pdb')).execute()
        )
        self.assertRaises(exc.FileDoesNotExist, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'), '-top',
                           self.get_fn('top_all22_prot.inp'), '-param', self.get_fn('par_all22_prot.inp'),
                           '-crd', 'foo').execute()
        )
        self.assertRaises(exc.FileDoesNotExist, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'), '-top',
                           self.get_fn('top_all22_prot.inp'), '-param', self.get_fn('par_all22_prot.inp'),
                           '-crd', self.get_fn('ala_ala_ala.pdb'), '-str', 'foo').execute()
        )
        self.assertRaises(exc.FileDoesNotExist, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'), '-top',
                           self.get_fn('top_all22_prot.inp'), '-param', self.get_fn('par_all22_prot.inp'),
                           '-crd', self.get_fn('ala_ala_ala.pdb'), '-toppar', 'foo').execute()
        )
        self.assertRaises(exc.InputError, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'),
                           '-top', self.get_fn('top_all22_prot.inp'),
                           '-param', self.get_fn('par_all22_prot.inp'),
                           '-crd', self.get_fn('ala_ala_ala.pdb'),
                           '-toppar', self.get_fn('trx.prmtop')).execute()
        )
        fn = self.get_fn('test.inp', written=True)
        with open(fn, 'w') as f:
            pass
        self.assertRaises(exc.InputError, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'), '-top',
                           self.get_fn('top_all22_prot.inp'), '-param', self.get_fn('par_all22_prot.inp'),
                           '-crd', self.get_fn('ala_ala_ala.pdb'), '-toppar', fn).execute()
        )
        self.assertRaises(exc.InputError, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'),
                           '-top', self.get_fn('top_all22_prot.inp'), '-param',
                           self.get_fn('par_all22_prot.inp'), '-crd', self.get_fn('*.pdb')).execute()
        )
        self.assertRaises(exc.InputError, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('*.psf'), '-top', self.get_fn('top_all22_prot.inp'),
                           '-param', self.get_fn('par_all22_prot.inp'), '-crd',
                           self.get_fn('ala_ala_ala.pdb')).execute()
        )
        self.assertRaises(exc.InputError, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'), '-top',
                           self.get_fn('top_all22_prot.inp'), '-param', self.get_fn('par_all22_prot.inp'),
                           '-crd', self.get_fn('ala_ala_ala.pdb'), '-box', 'a,b,c,d,e,f').execute()
        )
        self.assertRaises(exc.InputError, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'), '-top',
                           self.get_fn('top_all22_prot.inp'), '-param', self.get_fn('par_all22_prot.inp'),
                           '-crd', self.get_fn('ala_ala_ala.pdb'), '-box', '1,2,3,4').execute()
        )
        self.assertRaises(exc.InputError, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'),
                           '-top', self.get_fn('top_all22_prot.inp'),
                           '-crd', self.get_fn('ala_ala_ala.pdb')).execute()
        )
        self.assertRaises(exc.InputError, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'),
                           '-top', self.get_fn('top_all22_prot.inp'),
                           '-param', self.get_fn('par_all22_prot.inp'), '-radii', 'foobar').execute()
        )
        self.assertRaises(exc.InputError, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'),
                           '-top', self.get_fn('top_all22_prot.inp'),
                           '-param', self.get_fn('par_all22_prot.inp'), '-radii', 'foobar').execute()
        )
        self.assertRaises(exc.ChamberError, lambda:
                PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'), '-top',
                           self.get_fn('top_all22_prot.inp'), '-param', self.get_fn('par_all22_prot.inp'),
                           '-crd', self.get_fn('trx.prmtop')).execute()
        )

    def test_chamber_model(self):
        """ Test the chamber action with a model compound """
        a = PT.chamber(self.parm, '-psf', self.get_fn('propane.psf'), '-top',
                       self.get_fn('top_all36_prot.rtf'), '-param', self.get_fn('par_all36_prot.prm'),
                       '-str', self.get_fn('toppar_all36_prot_model.str'), '-str',
                       self.get_fn('toppar_water_ions.str'), '-crd', self.get_fn('propane.pdb'))
        str(a)
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
#       self._extensive_checks(parm)
        self.assertTrue(parm.chamber)
        self.assertEqual(parm.ptr('ifbox'), 0)

    def test_chamber_globbing(self):
        """ Test globbing in the chamber action """
        a = PT.chamber(self.parm, '-psf', self.get_fn('ala_ala_ala.psf'),
                       '-toppar', self.get_fn('*_all22_prot.inp'),
                       '-crd', self.get_fn('ala_ala_ala.pdb'))
        str(a)
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertTrue(parm.chamber)
        self.assertTrue(parm.has_cmap)
        self.assertEqual(parm.ptr('ifbox'), 0)

    def test_chamber_nbfix(self):
        """ Test the chamber action with a complex system using NBFIX """
        a = PT.chamber(self.parm, '-psf', self.get_fn('ala3_solv.psf'), '-toppar',
                       self.get_fn('???_all36_prot.???'), '-toppar', self.get_fn('toppar_water_ions.str'),
                       '-crd', self.get_fn('ala3_solv.crd'), '-box', 'bounding')
        a.execute()
        str(a)
        parm = a.parm
        self.assertTrue(parm.has_NBFIX())
        self._standard_parm_tests(parm)
#       self._extensive_checks(parm)

    def test_chamber_bug1(self):
        """ Test chamber BFNA creation (former bug) """
        a = PT.chamber(self.parm, '-top', self.get_fn('top_all36_cgenff.rtf'),
                '-top', self.get_fn('top_bfna_nonbonded_stitched.rtf'), '-param',
                self.get_fn('par_bfna_nonbonded_stitched.prm'), '-param',
                self.get_fn('par_all36_cgenff.prm'), '-box', 'bounding', '-psf',
                self.get_fn('bfna_nonbonded_vmd_autopsf.psf'), 'nocondense', '-crd',
                self.get_fn('bfna_nonbonded_vmd_autopsf.pdb'), '-nocmap'
        )
        str(a)
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        psf = CharmmPsfFile(self.get_fn('bfna_nonbonded_vmd_autopsf.psf'))
        self.assertEqual(len(psf.atoms), len(parm.atoms))
        self.assertEqual(len(psf.residues), len(parm.residues))
        for a1, a2 in zip(psf.atoms, parm.atoms):
            self.assertEqual(a1.name[:4], a2.name[:4])
            self.assertEqual(a1.type[:4], a2.type[:4])
            self.assertAlmostEqual(a1.charge, a2.charge)
            self.assertAlmostEqual(a1.mass, a2.mass)

    def test_chamber_bug2(self):
        """ Test that chamber sets the box angles for triclinics correctly """
        a = PT.chamber(self.parm, '-psf', self.get_fn('ala3_solv.psf'),
                       '-param', self.get_fn('par_all36_prot.prm'),
                       '-str', self.get_fn('toppar_water_ions.str'),
                       '-crd', self.get_fn('ala3_solv.crd'), '-box',
                       '33,33,33,109.475,109.475,109.475')
        a.execute()
        str(a)
        parm = a.parm
        self._standard_parm_tests(parm)
#       self._extensive_checks(parm)
        for x, y in zip(parm.parm_data['BOX_DIMENSIONS'], [109.475] + [33]*3):
            self.assertAlmostEqual(x, y)
        for x, y in zip(parm.box, [33]*3 + [109.475]*3):
            self.assertAlmostEqual(x, y)
        # Now check implied orthorhombic UC
        a = PT.chamber(self.parm, '-psf', self.get_fn('ala3_solv.psf'), '-param',
                       self.get_fn('par_all36_prot.prm'), '-str', self.get_fn('toppar_water_ions.str'),
                       '-crd', self.get_fn('ala3_solv.crd'), '-box', '33,33,33')
        a.execute()
        str(a)
        parm = a.parm
        for x, y in zip(parm.parm_data['BOX_DIMENSIONS'], [90] + [33]*3):
            self.assertEqual(x, y)
        has_wat = False
        for res in parm.residues:
            self.assertNotEqual(res.name, 'TIP3')
            if res.name == 'WAT':
                has_wat = True
                names = [a.name for a in res]
                self.assertNotIn('OH2', names)
                self.assertIn('O', names)
        self.assertTrue(has_wat)

    @unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without GROMACS")
    def test_gromber(self):
        """ Test the gromber action on a small system (no coords) """
        a = PT.gromber(None, os.path.join(self.get_fn('03.AlaGlu'), 'topol.top'))
        str(a)
        a.execute()
        parm = a.parm
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertIs(parm.box, None)

    @unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without GROMACS")
    def test_gromber2(self):
        """ Test the gromber action with coordinates """
        a = PT.gromber(None, os.path.join(self.get_fn('03.AlaGlu'), 'topol.top'),
                       os.path.join(self.get_fn('03.AlaGlu'), 'conf.gro'))
        str(a)
        a.execute()
        parm = a.parm
        parm.box = None
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertIs(parm.box, None) # AmberParm deletes the box without solvent
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))

    @unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without GROMACS")
    def test_gromber3(self):
        """ Test the gromber action passing various defines """
        a = PT.gromber(None, os.path.join(self.get_fn('03.AlaGlu'), 'topol.top'),
                       os.path.join(self.get_fn('03.AlaGlu'), 'conf.gro'),
                       'define', 'SOMEDEF=this', 'define', 'SOMEDEF2=that')
        stra = str(a)
        self.assertIn('SOMEDEF', stra)
        self.assertIn('SOMEDEF2', stra)
        a.execute()
        parm = a.parm
        parm.box = None
        self._standard_parm_tests(parm)
        self._extensive_checks(parm)
        self.assertIs(parm.box, None)
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))

    @unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without GROMACS")
    def test_gromber_box(self):
        """ Test the gromber action when a box should be defined """
        a = PT.gromber(None, self.get_fn('ala3.solv.top'), self.get_fn('ala3.solv.gro'))
        a.execute()
        str(a)
        parm = a.parm
        self._standard_parm_tests(parm)
#       self._extensive_checks(parm)
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))
        np.testing.assert_allclose(parm.box,
                [31.3585000, 31.3585000, 31.3584443,
                 60.0000468, 60.0000468, 90.0000000])

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

class TestAmberParmActions(FileIOTestCase, TestCaseRelative):
    """ Tests actions on Amber prmtop files """

    @unittest.skipIf(PYPY, 'Cannot test with NetCDF on pypy')
    def test_parmout_outparm_load_restrt(self):
        """ Test parmout, outparm, and loadRestrt actions on AmberParm """
        parm = copy(gasparm)
        lract = PT.loadRestrt(parm, self.get_fn('trx.inpcrd'))
        lract.execute()
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))
        act = PT.parmout(parm, self.get_fn('test.parm7', written=True))
        act.execute()
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 1)
        self.assertTrue(diff_files(self.get_fn('trx.prmtop'), self.get_fn('test.parm7', written=True)))
        self._empty_writes()
        PT.parmout(parm, self.get_fn('test.parm7', written=True), self.get_fn('test.rst7', written=True)).execute()
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 2)
        self.assertTrue(diff_files(self.get_fn('trx.prmtop'), self.get_fn('test.parm7', written=True)))
        self.assertTrue(diff_files(self.get_fn('trx.inpcrd'), self.get_fn('test.rst7', written=True),
                                   absolute_error=0.0001))
        self._empty_writes()
        PT.outparm(parm, self.get_fn('test.parm7', written=True)).execute()
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 1)
        self.assertTrue(diff_files(self.get_fn('trx.prmtop'), self.get_fn('test.parm7', written=True)))
        self._empty_writes()
        PT.outparm(parm, self.get_fn('test.parm7', written=True), self.get_fn('test.rst7', written=True)).execute()
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 2)
        self.assertTrue(diff_files(self.get_fn('trx.prmtop'), self.get_fn('test.parm7', written=True)))
        self.assertTrue(diff_files(self.get_fn('trx.inpcrd'), self.get_fn('test.rst7', written=True),
                                   absolute_error=0.0001))
        self._empty_writes()
        # Now try NetCDF format
        act2 = PT.outparm(parm, self.get_fn('test.parm7', written=True),
                          self.get_fn('test.ncrst', written=True), 'netcdf')
        act2.execute()
        self.assertTrue(
                amber.NetCDFRestart.id_format(self.get_fn('test.ncrst', written=True)))
        # Now make sure the string methods work
        str(act2); str(act); str(lract)
        # Make sure overwriting does not work
        self.assertRaises(exc.FileExists, act2.execute)
        # Delete the prmtop and make sure it still raises
        os.unlink(self.get_fn('test.parm7', written=True))
        self.assertRaises(exc.FileExists, act2.execute)

    def test_write_frcmod(self):
        """ Test writeFrcmod on AmberParm """
        parm = gasparm
        act = PT.writeFrcmod(parm, self.get_fn('test.frcmod', written=True))
        act.execute()
        self.assertTrue(diff_files(self.get_fn('test.frcmod', saved=True), self.get_fn('test.frcmod', written=True)))
        self.assertRaises(exc.FileExists, act.execute)
        # Make sure 10-12 prmtops fail
        parm = load_file(self.get_fn('ff91.parm7'))
        with self.assertWarns(exc.SeriousParmWarning):
            PT.writeFrcmod(load_file(self.get_fn('ff91.parm7')), self.get_fn('test2.frcmod', written=True))
        # Make sure string does not error
        str(act)

    def test_write_off_load_restrt(self):
        """ Test writeOFF on AmberParm """
        parm = copy(gasparm)
        parm.coordinates = None
        self.assertRaises(exc.ParmError, lambda:
                PT.writeOFF(parm, self.get_fn('test.off', written=True)).execute())
        PT.loadRestrt(parm, self.get_fn('trx.inpcrd')).execute()
        act = PT.writeOFF(parm, self.get_fn('test.off', written=True))
        act.execute()
        self.assertTrue(diff_files(self.get_fn('test.off', saved=True),
                                   self.get_fn('test.off', written=True),
                                   absolute_error=0.0001))
        str(act)
        self.assertRaises(exc.FileExists, act.execute)

    def test_change_radii(self):
        """ Test changeRadii on AmberParm """
        parm = copy(gasparm)
        act = PT.changeRadii(parm, 'amber6')
        act.execute()
        str(act) # Make sure it doesn't fail
        self.assertEqual(parm.parm_data['RADIUS_SET'][0],
                         'amber6 modified Bondi radii (amber6)')
        for i, atom in enumerate(parm.atoms):
            radii, atomic_number = atom.solvent_radius, atom.atomic_number
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

        # Make sure it adds back all necessary flags
        parm.delete_flag('RADIUS_SET')
        parm.delete_flag('RADII')
        parm.delete_flag('SCREEN')
        PT.changeRadii(parm, 'bondi').execute()
        self.assertEqual(parm.parm_data['RADIUS_SET'][0], 'Bondi radii (bondi)')
        for atom in parm.atoms:
            radii, atomic_number = atom.solvent_radius, atom.atomic_number
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
            radii, atomic_number = atom.solvent_radius, atom.atomic_number
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
            radii, atomic_number = atom.solvent_radius, atom.atomic_number
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
            radii, atomic_number = atom.solvent_radius, atom.atomic_number
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
        # Test that it works on non-Amber topologies
        act = PT.changeRadii(pmd.load_file(self.get_fn('ala_ala_ala.psf')), 'mbondi3')
        act.execute()

    def test_change_lj_pair(self):
        """ Test changeLJPair on AmberParm """
        parm = copy(gasparm)
        act = PT.changeLJPair(parm, '@%N', '@%H', 1.0, 1.0)
        act.execute()
        str(act)
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
        acoef = parm.parm_data['LENNARD_JONES_ACOEF'][:]
        bcoef = parm.parm_data['LENNARD_JONES_BCOEF'][:]
        refa = gasparm.parm_data['LENNARD_JONES_ACOEF'][:]
        refb = gasparm.parm_data['LENNARD_JONES_BCOEF'][:]
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
        # Check handling with no atoms
        PT.changeLJPair(parm, '@NONE', '@%H', 1.0, 1.0).execute()
        # Make sure nothing changed
        np.testing.assert_equal(acoef, parm.parm_data['LENNARD_JONES_ACOEF'])
        np.testing.assert_equal(bcoef, parm.parm_data['LENNARD_JONES_BCOEF'])
        # Check error handling
        self.assertRaises(exc.ChangeLJPairError, lambda:
                PT.changeLJPair(parm, ':*', '@%H', 1.0, 1.0).execute())
        self.assertRaises(exc.ChangeLJPairError, lambda:
                PT.changeLJPair(parm, '@%H', ':*', 1.0, 1.0).execute())

    def test_change_lj_14_pair(self):
        """ Check that changeLJ14Pair fails on AmberParm """
        parm = copy(gasparm)
        self.assertRaises(exc.ParmError, lambda:
            PT.changeLJ14Pair(parm, '@%N', '@%H', 1.0, 1.0).execute())

    def test_change(self):
        """ Test change on AmberParm with all properties """
        parm = copy(gasparm)
        act = PT.change(parm, 'CHARGE', ':ALA', 0, 'quiet')
        act.execute()
        self.assertEqual(len(str(act).split('\n')), 1)
        for flag in parm.parm_data:
            if flag != 'CHARGE':
                self.assertEqual(parm.parm_data[flag], gasparm.parm_data[flag])
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(parm.parm_data['CHARGE'][i], atom.charge)
            if atom.residue.name == 'ALA':
                self.assertEqual(atom.charge, 0)
            else:
                self.assertEqual(atom.charge, gasparm.atoms[i].charge)
        act = PT.change(parm, 'MASS', ':GLY', 10.0)
        self.assertGreater(len(str(act).split('\n')), 1)
        act.execute()
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
            self.assertEqual(atom.solvent_radius, parm.parm_data['RADII'][i])
            if atom.residue.idx < 20:
                self.assertEqual(atom.solvent_radius, 2.0)
            else:
                self.assertEqual(atom.solvent_radius,
                                 gasparm.atoms[i].solvent_radius)
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
        # Check if we try to add strings that are too long
        with self.assertWarns(exc.ParmWarning):
            PT.change(parm, 'ATOM_NAME', ':*', 'LALALALA').execute()
        # Make sure the strings get truncated
        PT.change(parm, 'ATOM_NAME', ':*', 'LALALALA').execute()
        for atom in parm.atoms:
            self.assertEqual(atom.name, 'LALA')
            self.assertEqual(parm.parm_data['ATOM_NAME'][atom.idx], 'LALA')
        # Check bad input
        self.assertRaises(exc.ParmedChangeError, lambda:
                          PT.change(parm, 'RESIDUE_LABEL', ':*', 'NaN'))

    def test_print_info(self):
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
        with self.assertWarns(exc.SeriousParmWarning):
            repr(PT.printInfo(gasparm, 'NOTAFLAG'))

    def test_add_change_lj_type(self):
        """ Test addLJType and changeLJSingleType on AmberParm """
        parm = copy(gasparm)
        act = PT.addLJType(parm, '@1')
        act.execute()
        str(act) # Make sure it doesn't crash
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

    def test_print_lj_types(self):
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
        # Integer
        act = PT.printLJTypes(gasparm, 1)
        for line in str(act).split('\n'):
            if not line.startswith('ATOM'):
                continue
            words = line.split()
            self.assertTrue(words[2].startswith('N'))
            self.assertTrue(words[3].startswith('N'))
            self.assertEqual(words[7], '1')
        # Integer tuple
        act = PT.printLJTypes(gasparm, '1,2')
        for line in str(act).split('\n'):
            if not line.startswith('ATOM'):
                continue
            words = line.split()
            self.assertIn(words[2][0], 'NH')
            self.assertIn(words[3][0], 'NH')
            self.assertIn(words[7], ('1', '2'))
        # Default prints out every atom
        i = 0
        for line in str(PT.printLJTypes(gasparm)).split('\n'):
            if not line.startswith('ATOM'):
                continue
            i += 1
        self.assertEqual(i, len(gasparm.atoms))
        # Now check ranges and such
        act = PT.printLJTypes(gasparm, '1-3,5')
        it = (a for a in gasparm.atoms if a.nb_idx in (1, 2, 3, 5))
        for line in str(act).split('\n'):
            if not line.startswith('ATOM'):
                continue
            words = line.split()
            a = next(it)
            self.assertEqual(int(words[1]), a.idx+1)
            self.assertEqual(words[2], a.name)
            self.assertEqual(words[3], a.type)
            self.assertEqual(int(words[7]), a.nb_idx)
            i += 1
        # Check that it still works when nothing selected
        str(PT.printLJTypes(gasparm, '@NOATOM'))
        # Check illegal selections
        self.assertRaises(exc.ParmError, lambda:
                PT.printLJTypes(gasparm, 0))
        self.assertRaises(exc.ParmError, lambda:
                PT.printLJTypes(gasparm, '0-0'))

    def test_scee_scnb(self):
        """ Test scee and scnb actions on AmberParm """
        parm = copy(gasparm)
        act1 = PT.scee(parm, 1.0)
        act2 = PT.scnb(parm, 1.0)
        act1.execute()
        act2.execute()
        str(act1), str(act2)
        for dih in parm.dihedrals:
            self.assertEqual(dih.type.scee, 1.0)
            self.assertEqual(dih.type.scnb, 1.0)
        for x, y in zip(parm.parm_data['SCEE_SCALE_FACTOR'],
                        parm.parm_data['SCNB_SCALE_FACTOR']):
            self.assertEqual(x, 1.0)
            self.assertEqual(y, 1.0)

    def test_print_details(self):
        """ Test printDetails on AmberParm """
        act = PT.printDetails(gasparm, '@1')
        self.assertEqual(str(act), saved.PRINT_DETAILS)
        # Check out different ways of selecting parms
        act = PT.printDetails(gasparm, '@1', parm=0)
        self.assertEqual(str(act), saved.PRINT_DETAILS)
        act = PT.printDetails(gasparm, '@1', parm=str(gasparm))
        self.assertEqual(str(act), saved.PRINT_DETAILS)
        with self.assertWarns(exc.SeriousParmWarning):
            PT.printDetails(gasparm, '@1', parm='DNE')
        with self.assertWarns(exc.SeriousParmWarning):
            PT.printDetails(gasparm, '@1', parm=10)

    def test_print_flags(self):
        """ Test printFlags on AmberParm """
        act = PT.printFlags(gasparm)
        printed_flags = set()
        for line in str(act).split('\n'):
            if line.startswith('%FLAG'):
                printed_flags.add(line.split()[1])
        self.assertEqual(printed_flags, set(gasparm.parm_data.keys()))

    def test_print_pointers(self):
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

    def test_print_bonds(self):
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

    def test_print_bonds_with_measurements(self):
        """ Test printBonds on AmberParm with measurements """
        act = PT.printBonds(AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7')), '*')
        self.assertEqual(str(act), saved.PRINT_BONDS_MEASURE)

    def test_print_angles(self):
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
        act = PT.printAngles(gasparm, *'@1 @5 @7'.split())
        self.assertEqual(str(act), saved.PRINT_ANGLES_3MASKS)
        act = PT.printAngles(gasparm, *'@7 @5 @1'.split())
        self.assertEqual(str(act), saved.PRINT_ANGLES_3MASKS)

    def test_print_angles_with_measurements(self):
        """ Test printBonds on AmberParm with measurements """
        act = PT.printAngles(AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7')), '*')
        self.assertEqual(str(act), saved.PRINT_ANGLES_MEASURE)

    def test_print_dihedrals(self):
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
        # Now make sure other combos all work
        str(PT.printDihedrals(gasparm, ':*'))
        str(PT.printDihedrals(gasparm, ':1-10', ':2-11'))
        str(PT.printDihedrals(gasparm, ':*', ':1-10', ':1-10', ':*'))
        str(PT.printDihedrals(gasparm, ':*', ':*', ':1-10', ':1-10'))

    def test_print_angles_with_measurements(self):
        """ Test printBonds on AmberParm with measurements """
        act = PT.printDihedrals(AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7')), '*')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALS_MEASURE)

    def test_set_molecules(self):
        """ Test setMolecules on AmberParm """
        parm = AmberParm(self.get_fn('things.parm7'), self.get_fn('things.rst7'))
        atoms = [atom for atom in parm.atoms] # shallow copy!
        self.assertTrue(all([x is y for x,y in zip(parm.atoms,atoms)]))
        self.assertEqual(parm.ptr('IPTRES'), 29)
        self.assertEqual(parm.ptr('NSPM'), 718)
        self.assertEqual(parm.ptr('NSPSOL'), 23)
        # To keep the output clean
        PT.setMolecules(parm).execute()
        self.assertFalse(all([x is y for x,y in zip(parm.atoms,atoms)]))
        # Now check that setMolecules can apply another time.
        PT.setMolecules(parm, solute_ions=False).execute()

        # Now check that solute_ions keyword works as expected
        parm = AmberParm(self.get_fn('ff14ipq.parm7'), self.get_fn('ff14ipq.rst7'))
        self.assertEqual(parm.parm_data['SOLVENT_POINTERS'], [15, 926, 12])
        PT.setMolecules(parm, solute_ions=False).execute()
        self.assertEqual(parm.parm_data['SOLVENT_POINTERS'], [5, 926, 2])
        act = PT.setMolecules(parm, solute_ions=True)
        act.execute()
        str(act)
        self.assertEqual(parm.parm_data['SOLVENT_POINTERS'], [15, 926, 12])
        with self.assertWarns(exc.SeriousParmWarning):
            PT.setMolecules(parm, solute_ions='foo')

        # Make sure we warn if reordering was required
        with self.assertWarns(exc.ParmWarning):
            PT.setMolecules(pmd.load_file(self.get_fn('things.parm7')), solute_ions=True).execute()

    def test_net_charge(self):
        """ Test netCharge on AmberParm """
        act = PT.netCharge(gasparm)
        chg = act.execute() # check this part of the API
        self.assertEqual(str(act), 'The net charge of :* is %.4f' % chg)
        self.assertAlmostEqual(chg, -4.0, places=6)
        chg = PT.netCharge(gasparm, ':ASP').execute()
        self.assertAlmostEqual(chg, -10.0, places=6)

    def test_strip(self):
        """ Test stripping of AmberParm """
        parm = copy(gasparm)
        parm.box = [10, 10, 10, 90, 90, 90]
        act = PT.strip(parm, ':1')
        act.execute()
        self.assertEqual(parm.ptr('natom'), 1641)
        self.assertEqual(len(parm.atoms), 1641)
        self.assertEqual(str(act), "Removing mask ':1' (%d atoms) "
                         "from the topology file." % (len(gasparm.residues[0])))
        np.testing.assert_equal(parm.box, [10, 10, 10, 90, 90, 90])
        # Now test nobox
        act = PT.strip(parm, ':1', 'nobox')
        act.execute()
        str(act)
        self.assertIs(parm.box, None)

    def test_define_solvent(self):
        """ Test defineSolvent on AmberParm """
        import parmed.residue as residue
        PT.defineSolvent(gasparm, 'WAT,HOH,Na+,Cl-').execute()
        self.assertEqual(residue.SOLVENT_NAMES, 'WAT HOH Na+ Cl-'.split())
        act = PT.defineSolvent(gasparm, 'WAT,HOH,')
        act.execute()
        str(act)
        self.assertEqual(residue.SOLVENT_NAMES, 'WAT HOH'.split())

    def test_add_exclusions(self):
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
        act = PT.addExclusions(parm, ':1', ':1')
        act.execute()
        str(act)
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
        # Make sure the exclusions are correct when reading in a new version
        # with exceptions defined
        parm.remake_parm()
        parm = AmberParm.from_rawdata(parm)
        in_exclusions_before = []
        for atom1 in parm.residues[0].atoms:
            all_exclusions = (atom1.bond_partners + atom1.angle_partners +
                             atom1.dihedral_partners + atom1.exclusion_partners)
            for atom2 in parm.residues[0].atoms:
                if atom1 is atom2: continue
                in_exclusions_before.append(atom2 in all_exclusions)
        self.assertTrue(all(in_exclusions_before))

    def test_add_delete_dihedral(self):
        """ Test addDihedral and deleteDihedral on AmberParm """
        parm = copy(gasparm)
        act = PT.deleteDihedral(parm, *':ALA@N :ALA@CA :ALA@CB :ALA@HB1'.split())
        n = act.execute()
        str(act)
        str(PT.deleteDihedral(parm, '@NONE', '@NONE', '@NONE', '@NONE'))
        parm.remake_parm()
        self.assertEqual(gasparm.ptr('nphih') + gasparm.ptr('nphia'),
                         parm.ptr('nphih') + parm.ptr('nphia') + n)
        NALA = sum([res.name == 'ALA' for res in parm.residues])
        self.assertEqual(n, NALA)
        act = PT.addDihedral(parm, ':ALA@N', ':ALA@CA', ':ALA@CB', ':ALA@HB1',
                       0.1556, 3, 0, 1.2, 2.0)
        act.execute()
        str(act)
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
        # Try adding a new improper
        ntyp = len(parm.dihedral_types)
        for d in parm.dihedrals:
            if d.improper:
                dt = d.type
                break
        else:
            assert False, 'No natural impropers found'
        PT.addDihedral(parm, '@1', '@20', '@30', '@40', dt.phi_k,
            dt.per, dt.phase, dt.scee, dt.scnb, type='improper').execute()
        self.assertEqual(len(parm.dihedral_types), ntyp)
        self.assertIs(parm.dihedrals[-1].type, dt)
        # Error checking
        self.assertRaises(exc.SetParamError,
            PT.addDihedral(parm, '@1-2', '@3', '@4', '@5', 1, 1, 10).execute)
        self.assertRaises(exc.InputError, lambda:
            PT.addDihedral(parm, '@1', '@2', '@3', '@4', 1, 1, 10, type='badtype').execute()
        )
        self.assertRaises(exc.SetParamError, lambda:
            PT.addDihedral(parm, '@1,3', '@2,3', '@4-5', '@6-7', 1, 1, 10).execute())
        self.assertRaises(exc.DeleteDihedralError, lambda:
            PT.deleteDihedral(parm, '@1', '@2-3', '@4', '@5').execute())
        with self.assertWarns(exc.SeriousParmWarning):
            PT.deleteDihedral(parm, '@1,3', '@2-3', '@4-5', '@6-7').execute()
        self.assertEqual(PT.deleteDihedral(parm, '@1', '@25', '@35', '@45').execute(), 0)

    def test_set_bond(self):
        """ Test setBond on AmberParm """
        parm = copy(gasparm)
        PT.setBond(parm, ':ALA@CA', ':ALA@CB', 300.0, 1.5).execute()
        act = PT.printBonds(parm, ':ALA@CA')
        self.assertEqual(str(act), saved.SET_BOND)
        nala = sum([1 for res in parm.residues if res.name == 'ALA'])
        nbon = len(parm.bonds)
        act = PT.setBond(parm, ':ALA@CB', ':ALA@HA', 100.0, 1.5)
        act.execute()
        str(act)
        for atom in parm.atoms:
            if atom.residue.name == 'ALA' and atom.name == 'HA':
                self.assertEqual(len(atom.bonds), 2)
                for atom2 in atom.bond_partners:
                    self.assertIn(atom2.name, ('CA','CB'))
        self.assertEqual(nbon + nala, len(parm.bonds))
        # Add a new bond with an existing bond type and make sure a new type is
        # not added
        ntyp = len(parm.bond_types)
        PT.setBond(parm, '@1', '@20', parm.bond_types[0].k,
                   parm.bond_types[0].req).execute()
        self.assertEqual(ntyp, len(parm.bond_types))
        self.assertIs(parm.bonds[-1].type, parm.bond_types[0])
        # Error checking
        self.assertRaises(exc.SetParamError, lambda:
                PT.setBond(parm, '@1', '@2-3', 300.0, 1.5).execute())

    def test_set_angle(self):
        """ Test setAngle on AmberParm """
        parm = copy(gasparm)
        act = PT.setAngle(parm, ':ALA@CA', ':ALA@CB', ':ALA@HB1', 40, 100)
        act.execute()
        str(act)
        act = PT.printAngles(parm, ':ALA@CB')
        self.assertEqual(str(act), saved.SET_ANGLE)
        nala = sum([1 for res in parm.residues if res.name == 'ALA'])
        nang = len(parm.angles)
        PT.setAngle(parm, ':ALA@HA', ':ALA@CB', ':ALA@HB1', 50, 120).execute()
        self.assertEqual(nang + nala, len(parm.angles))
        # Add a new angle with an existing type and make sure it uses that
        angtyp = parm.angle_types[0]
        ntyp = len(parm.angle_types)
        PT.setAngle(parm, '@1', '@20', '@30', angtyp.k, angtyp.theteq).execute()
        self.assertEqual(len(parm.angle_types), ntyp)
        self.assertIs(parm.angles[-1].type, parm.angle_types[0])
        # Error checking
        self.assertRaises(exc.SetParamError,
                PT.setAngle(parm, '@1', '@2', '@3-4', 50, 120).execute)
        self.assertRaises(exc.SetParamError,
                PT.setAngle(parm, '@1', '@2-3', '@4', 50, 120).execute)

    def test_add_atomic_number(self):
        """ Test addAtomicNumber on AmberParm """
        parm = copy(gasparm)
        self.assertFalse('ATOMIC_NUMBER' in parm.parm_data)
        atomic_numbers = [atom.atomic_number for atom in parm.atoms]
        act = PT.addAtomicNumber(parm)
        act.execute()
        str(act)
        self.assertEqual(parm.parm_data['ATOMIC_NUMBER'], atomic_numbers)
        str(PT.addAtomicNumber(parm)) # Already present

    def test_print_lj_matrix(self):
        """ Test printLJMatrix on AmberParm """
        act = PT.printLJMatrix(gasparm, '@1')
        self.assertEqual(str(act), saved.PRINT_LJMATRIX)
        act = PT.printLJMatrix(gasparm, gasparm[0].nb_idx)
        self.assertEqual(str(act), saved.PRINT_LJMATRIX)

    def test_delete_bond(self):
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
        act = PT.deleteBond(parm, '@1', '@5', 'verbose')
        str(act)
        act.execute()
        self.assertTrue(bond not in parm.bonds)
        for angle in deleted_angles:
            self.assertTrue(angle not in parm.angles)
        for dihedral in deleted_dihedrals:
            self.assertTrue(dihedral not in parm.dihedrals)
        # Nothing to do, make sure it doesn't fail, and does nothing
        act = PT.deleteBond(parm, '@1', '@20')
        nbnd = len(parm.bonds)
        str(act)
        act.execute()
        self.assertEqual(nbnd, len(parm.bonds))

    def test_summary(self):
        """ Test summary action on AmberParm """
        parm = AmberParm(self.get_fn('things.parm7'), self.get_fn('things.rst7'))
        act = PT.summary(parm)
        self.assertEqual(str(act), saved.SUMMARY)

    def test_scale(self):
        """ Test scale action on AmberParm """
        parm = copy(gasparm)
        PT.scale(parm, 'DIHEDRAL_FORCE_CONSTANT', 2.0).execute()
        self.assertEqual(
                [2*x for x in gasparm.parm_data['DIHEDRAL_FORCE_CONSTANT']],
                parm.parm_data['DIHEDRAL_FORCE_CONSTANT'])
        PT.scale(parm, 'DIHEDRAL_FORCE_CONSTANT', 0.0).execute()
        for val in parm.parm_data['DIHEDRAL_FORCE_CONSTANT']:
            self.assertEqual(val, 0)

    def test_lmod(self):
        """ Test lmod action on AmberParm """
        parm = copy(gasparm)
        self.assertFalse(all(parm.parm_data['LENNARD_JONES_ACOEF']))
        act = PT.lmod(parm)
        act.execute()
        str(act)
        self.assertTrue(all(parm.parm_data['LENNARD_JONES_ACOEF']))

    def test_prot_state_interpolate(self):
        """ Test changeProtState and interpolate actions on AmberParm """
        self._empty_writes()
        parm = AmberParm(self.get_fn('ash.parm7'))
        origparm = copy(parm)
        origparm.name = origparm.name + '_copy1'
        act = PT.changeProtState(parm, ':ASH', 0)
        act.execute()
        str(act)
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
        # First do some error checking
        self.assertRaises(exc.ArgumentError, lambda:
                PT.interpolate(parms, -1, 'eleconly',
                    prefix=self.get_fn('test.parm7', written=True)).execute())
        sys.stdout = open(os.devnull, 'w')
        act = PT.interpolate(parms, 5, 'eleconly', startnum=2,
                       prefix=self.get_fn('test.parm7', written=True))
        str(act)
        act.execute()
        sys.stdout = sys.__stdout__
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 5)
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.2', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.3', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.4', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.5', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.6', written=True)))
        # Now check them all
        ladder = [origparm]
        ladder.append(AmberParm(self.get_fn('test.parm7.2', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.3', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.4', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.5', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.6', written=True)))
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
        # Now check that if no atoms are selected, it is handled properly
        parm = pmd.load_file(self.get_fn('ash.parm7'))
        act = PT.changeProtState(parm, ':GLH', 1)
        str(act)
        act.execute()
        # Now check error handling if the selected residue is not protonable
        self.assertRaises(exc.ChangeStateError,
                PT.changeProtState(parm, ':ACE', 1).execute)
        self.assertRaises(exc.ChangeStateError,
                PT.changeProtState(parm, ':ASH', 10).execute)
        self.assertRaises(exc.ChangeStateError,
                PT.changeProtState(parm, ':*', 1).execute)
        self.assertRaises(exc.ChangeStateError,
                PT.changeProtState(parm, ':ASH@CA', 1).execute)
        self.assertRaises(exc.ChangeStateError,
                PT.changeProtState(parm, '@13-25', 1).execute)
        self.assertRaises(exc.NonexistentParm, lambda:
                PT.interpolate(load_file(self.get_fn('ash.parm7')), 10, 'eleconly').execute())
        # Make a list of 3 parms and make sure ambiguity causes failure
        act = PT.parm(parm, copy=0)
        act.execute()
        PT.changeProtState(act.parm_list, ':ASH', 1).execute()
        PT.parm(act.parm_list, copy=0).execute()
        self.assertEqual(len(act.parm_list), 3)
        self.assertRaises(exc.AmbiguousParmError, lambda:
                PT.interpolate(act.parm_list, 10, 'eleconly').execute())
        act.parm_list.add_parm(self.get_fn('trx.prmtop'))
        self.assertRaises(exc.IncompatibleParmsError, lambda:
                PT.interpolate(act.parm_list, 10, 'eleconly', parm2=0).execute())
        act.parm_list[0].atoms[0].name = 'FOO'
        with self.assertWarns(exc.SeriousParmWarning):
            PT.interpolate(act.parm_list, 10, 'eleconly', parm=0, parm2=1).execute()

    def test_prot_state_interpolate_2(self):
        """ Test interpolate actions on AmberParm with more parm selection """
        self._empty_writes()
        parm = AmberParm(self.get_fn('ash.parm7'))
        origparm = copy(parm)
        origparm.name = origparm.name + '_copy1'
        act = PT.changeProtState(parm, ':ASH', 0)
        act.execute()
        str(act)
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
        parms.add_parm(load_file(self.get_fn('trx.prmtop')))
        sys.stdout = open(os.devnull, 'w')
        PT.interpolate(parms, 5, 'eleconly', parm=1, parm2=0, startnum=2,
                       prefix=self.get_fn('test.parm7', written=True)).execute()
        sys.stdout = sys.__stdout__
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 5)
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.2', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.3', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.4', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.5', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.6', written=True)))
        # Now check them all
        ladder = [origparm]
        ladder.append(AmberParm(self.get_fn('test.parm7.2', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.3', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.4', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.5', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.6', written=True)))
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

    def test_prot_state_interpolate_3(self):
        """ Test interpolate actions on AmberParm with parm name selection """
        self._empty_writes()
        parm = AmberParm(self.get_fn('ash.parm7'))
        origparm = copy(parm)
        origparm.name = origparm.name + '_copy1'
        act = PT.changeProtState(parm, ':ASH', 0)
        act.execute()
        str(act)
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
        parms.add_parm(load_file(self.get_fn('trx.prmtop')))
        sys.stdout = open(os.devnull, 'w')
        PT.interpolate(parms, 5, 'eleconly', parm=1, parm2=self.get_fn('ash.parm7'),
                       startnum=2, prefix=self.get_fn('test.parm7', written=True)).execute()
        sys.stdout = sys.__stdout__
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 5)
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.2', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.3', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.4', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.5', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.6', written=True)))
        # Now check them all
        ladder = [origparm]
        ladder.append(AmberParm(self.get_fn('test.parm7.2', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.3', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.4', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.5', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.6', written=True)))
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

    def test_add_delete_pdb(self):
        """ Test addPDB and deletePDB actions on AmberParm """
        parm = copy(gasparm)
        PT.addPDB(parm, self.get_fn('trx.pdb'), 'elem', 'allicodes').execute()
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
        # Test reading Amber topology file with insertion codes and PDB
        # information stored in it
        parm.write_parm(self.get_fn('addpdb.parm7', written=True))
        nparm = AmberParm(self.get_fn('addpdb.parm7', written=True))
        self.assertEqual(len(nparm.parm_data['RESIDUE_ICODE']),
                         len(nparm.residues))
        nparm = AmberFormat()
        nparm.rdparm(self.get_fn('addpdb.parm7', written=True), slow=True)
        self.assertEqual(len(nparm.parm_data['RESIDUE_ICODE']),
                         len(parm.residues))
        # Test deletePDB
        act = PT.deletePDB(parm)
        act.execute()
        str(act)
        self.assertFalse('RESIDUE_ICODE' in parm.flag_list)
        self.assertFalse('ATOM_ELEMENT' in parm.flag_list)
        self.assertFalse('RESIDUE_NUMBER' in parm.flag_list)
        self.assertFalse('RESIDUE_CHAINID' in parm.flag_list)
        self.assertFalse('ATOM_OCCUPANCY' in parm.flag_list)
        self.assertFalse('ATOM_BFACTOR' in parm.flag_list)
        act = PT.deletePDB(parm)
        act.execute()
        str(act)
        # Add PDB information, but set a couple new names that should be OK.
        # Make sure we don't warn
        for res in parm.residues:
            if res.name == 'ASP':
                res.name = 'AS4'
            elif res.name == 'GLU':
                res.name = 'GL4'
            elif res.name == 'HIS':
                res.name = 'HID'
            elif res.name == 'LYS':
                res.name = 'LYN'
            elif res.name == 'CYS':
                res.name = 'CYM'
        parm.remake_parm()
        PT.addPDB(parm, self.get_fn('trx.pdb'), 'elem', 'allicodes').execute()
        # Make sure addPDB works with nucleic acids, too
        parm2 = AmberParm(self.get_fn('gaucu.parm7'))
        PT.addPDB(parm2, self.get_fn('gaucu.pdb'), 'elem', 'allicodes').execute()
        self.assertTrue('RESIDUE_ICODE' in parm2.flag_list)
        self.assertTrue('ATOM_ELEMENT' in parm2.flag_list)
        self.assertTrue('RESIDUE_NUMBER' in parm2.flag_list)
        self.assertTrue('RESIDUE_CHAINID' in parm2.flag_list)
        self.assertTrue('ATOM_OCCUPANCY' in parm2.flag_list)
        self.assertTrue('ATOM_BFACTOR' in parm2.flag_list)
        self.assertTrue(len(parm2.parm_data['RESIDUE_ICODE']), parm2.ptr('nres'))
        self.assertTrue(len(parm2.parm_data['ATOM_ELEMENT']), parm2.ptr('natom'))
        self.assertTrue(len(parm2.parm_data['RESIDUE_NUMBER']), parm2.ptr('nres'))
        self.assertTrue(len(parm2.parm_data['RESIDUE_CHAINID']),parm2.ptr('nres'))
        self.assertTrue(len(parm2.parm_data['ATOM_OCCUPANCY']),parm2.ptr('natom'))
        self.assertTrue(len(parm2.parm_data['ATOM_BFACTOR']), parm2.ptr('natom'))
        PT.deletePDB(parm2).execute()
        # Now tweak some residue name
        parm2.residues[0].name = 'FOO'
        parm2.remake_parm()
        self.assertRaises(exc.AddPDBWarning, lambda:
                PT.addPDB(parm2, self.get_fn('gaucu.pdb'), 'elem', 'strict').execute())
        # Now mix-and-match PDB files that do not have the same number of atoms
        parm2.residues[0].name = 'G5'
        parm2.strip(':5')
        self.assertRaises(exc.AddPDBError, lambda:
                PT.addPDB(parm2, self.get_fn('gaucu.pdb'), 'elem').execute())

    def test_add_pdb2(self):
        """ Test addPDB with atypical numbering and extra residues """
        parm = load_file(self.get_fn('4lzt.parm7'))
        act = PT.addPDB(parm, self.get_fn('4lzt_NoNO3.pdb'), 'strict')
        act.execute()
        str(act)
        parm.write_parm(self.get_fn('4lzt_pdb.parm7', written=True))
        self.assertTrue(diff_files(get_saved_fn('4lzt_pdb.parm7'),
                                   self.get_fn('4lzt_pdb.parm7', written=True),
                                   absolute_error=1e-6)
        )
        act = PT.addPDB(parm, self.get_fn('4lzt_NoNO3.pdb'))
        str(act)
        act.execute()

    def test_h_mass_repartition(self):
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
        self.assertRaises(exc.HMassRepartitionError, lambda:
                PT.HMassRepartition(parm, 100.0).execute())

    @unittest.skipUnless(has_openmm, 'Cannot test without OpenMM')
    def test_openmm_action(self):
        """ Tests the OpenMM action for AmberParm """
        parm = AmberParm(self.get_fn('ash.parm7'))
        mdin = self.get_fn('mdin', written=True)
        with open(mdin, 'w') as f:
            f.write('''\
Basic MD simulation
 &cntrl
    imin=0, nstlim=10, dt=0.001, ntb=0, igb=5,
    ntwr=2, ntwx=2, ntpr=2, ntwv=-1, ioutfm=1,
    cut=500.0, tempi=100,
 /
''')
        traj = self.get_fn('ash.nc', written=True)
        mdout = self.get_fn('ash.mdout', written=True)
        mdinfo = self.get_fn('ash.mdinfo', written=True)
        restart = self.get_fn('ash.restrt', written=True)
        script = self.get_fn('ash.py', written=True)
        self.assertRaises(exc.SimulationError, lambda:
                PT.OpenMM(parm, '-O', '-i', mdin, '-o', mdout, '-inf', mdinfo,
                          '-r', restart, '-x', traj, '-p', self.get_fn('ash.parm7'),
                          script=script).execute()
        )
        act = PT.OpenMM(parm, '-O', '-i', mdin, '-c', self.get_fn('ash.rst7'), '-o',
                        mdout, '-inf', mdinfo, '-r', restart, '-x', traj,
                        '-p', self.get_fn('ash.parm7'), script=script)
        str(act)
        act.execute()
        self.assertTrue(os.path.exists(script))
        self.assertTrue(pmd.amber.AmberAsciiRestart.id_format(restart))
        self.assertTrue(pmd.amber.NetCDFTraj.id_format(traj))
        # Now check that the attributes are as expected
        if pd is not None:
            df = pd.read_csv(mdout)
            self.assertEqual(df.shape, (5, 6))
            np.testing.assert_allclose(df['Time (ps)'],
                    [0.002, 0.004, 0.006, 0.008, 0.010])
        # Now check that the coordinates have moved
        t = pmd.load_file(traj)
        self.assertEqual(t.coordinates.shape, (5, len(parm.atoms), 3))
        diff = pmd.load_file(self.get_fn('ash.rst7')).coordinates - t.coordinates[4]
        self.assertTrue((np.abs(diff) > 1e-3).any())

        # Now run the python script and make sure it does the same thing
        os.system('cd "%s" && "%s" "%s"' % (self._temporary_directory.name, sys.executable, script))
        if pd is not None:
            df = pd.read_csv(mdout)
            self.assertEqual(df.shape, (5, 6))
            np.testing.assert_allclose(df['Time (ps)'],
                    [0.002, 0.004, 0.006, 0.008, 0.010])
        # Now check that the coordinates have moved
        t = pmd.load_file(traj)
        self.assertEqual(t.coordinates.shape, (5, len(parm.atoms), 3))
        diff = pmd.load_file(self.get_fn('ash.rst7')).coordinates - t.coordinates[4]
        self.assertTrue((np.abs(diff) > 1e-3).any())

    @unittest.skipUnless(has_openmm, 'Cannot test energy function without OMM')
    def test_energy_openmm(self):
        """ Tests the energy action with OpenMM """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        f = StringIO()
        PT.energy.output = f
        act = PT.energy(parm, 'omm', 'decompose', igb=5, saltcon=0.1)
        str(act)
        act.execute()
        f.seek(0)
        info = f.read()
        ene = float(re.findall(r'TOTAL\s+=\s+([-\d\.]+)', info)[0])
        self.assertLess(abs(ene + 23.01), 0.05)

    @unittest.skipIf(sander is None, 'Cannot test energy function without pysander')
    def test_energy_sander_gb(self):
        """ Tests gb energy action with sander """
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        f = StringIO()
        PT.energy.output = f
        PT.energy(parm, igb=5, saltcon=0.1).execute()
        f.seek(0)
        info = f.read()
        ene = float(re.findall(r'TOTAL\s+=\s+([-\d\.]+)', info)[0])
        self.assertLess(abs(ene + 23.01), 0.05)
        self.assertRaises(exc.SimulationError, lambda: PT.energy(parm, igb=100).execute())
        self.assertRaises(exc.SimulationError, lambda: PT.energy(parm, cutoff=-2.0).execute())
        self.assertRaises(exc.SimulationError, lambda: PT.energy(parm, saltcon=-0.2).execute())
        PT.energy(parm, igb=0).execute()
        parm.coordinates = None
        self.assertRaises(exc.SimulationError, lambda: PT.energy(parm).execute())

    @unittest.skipIf(sander is None, 'Cannot test energy function without pysander')
    def test_energy_sander_explicit_water(self):
        """ Tests energy action with sander """
        parm = AmberParm(self.get_fn('solv2.parm7'), self.get_fn('solv2.rst7'))
        f = StringIO()
        PT.energy.output = f
        PT.energy(parm).execute()
        f.seek(0)
        info = f.read()
        ene = float(re.findall(r'TOTAL\s+=\s+([-\d\.]+)', info)[0])
        self.assertLess(abs(ene + 12785.68), 0.05)
        self.assertRaises(exc.SimulationError, lambda: PT.energy(parm, cutoff=-2.0).execute())

    @unittest.skipIf(sander is None, 'Cannot test amber minimization without pysander')
    def test_minimize_from_action_tools(self):
        """ Tests the minimize action with pysander and scipy """
        # just want to make sure those minimizations runnable
        parm7 = self.get_fn('ala_ala_ala.parm7')
        rst7 = self.get_fn('ala_ala_ala.rst7')
        parm = pmd.load_file(parm7, rst7)
        original_coordinates = parm.coordinates

        for igb in (0, 1, 2, 5, 6, 7, 8):
            arg_list = dict(igb=igb, maxcyc=10)
            parm.coordinates = original_coordinates
            pmd.tools.minimize(parm, **arg_list).execute()

        def test_coordinates_is_None():
            parm.coordinates = None
            pmd.tools.minimize(parm, igb=8, maxcyc=10).execute()
        self.assertRaises(exc.SimulationError, test_coordinates_is_None)

    @unittest.skipIf(True, "OpenMM minimization test is currently failing") # Need to figure out why this is segfaulting in parallel
    @unittest.skipUnless(has_openmm, 'Cannot test minimize function without OpenMM')
    def test_minimize_openmm(self):
        """ Tests the minimize action with OpenMM """
        self._check_emin_omm(AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7')), 0)
        self._check_emin_omm(AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7')), 1)
        self._check_emin_omm(AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7')), 2)
        self._check_emin_omm(AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7')), 5)
        parm = AmberParm(self.get_fn('ash.parm7'), self.get_fn('ash.rst7'))
        PT.changeRadii(parm, 'mbondi3').execute()
        self._check_emin_omm(parm, 7)
        self._check_emin_omm(parm, 8)

    def _check_emin_omm(self, parm, igb):
        if igb in (0, 6):
            IS = None
        elif igb == 1:
            IS = app.HCT
        elif igb == 2:
            IS = app.OBC1
        elif igb == 5:
            IS = app.OBC2
        elif igb == 7:
            IS = app.GBn
        elif igb == 8:
            IS = app.GBn2
        else:
            assert False, 'illegal input'
        system = parm.createSystem(implicitSolvent=IS, implicitSolventSaltConc=0.1*u.molar)
        starting_e = 0
        for _, e in pmd.openmm.energy_decomposition_system(parm, system):
            starting_e += e
        script = self.get_fn('minimize.py', written=True)
        act = PT.minimize(parm, 'omm', script=script, maxcyc=10, igb=igb, saltcon=0.1)
        str(act)
        act.execute()
        end_e = 0
        for _, e in pmd.openmm.energy_decomposition_system(parm, system):
            end_e += e
        self.assertLess(end_e, starting_e)

    def test_out_pdb(self):
        """ Test the outPDB action on AmberParm """
        parm = copy(gasparm)
        PT.loadRestrt(parm, self.get_fn('trx.inpcrd')).execute()
        PT.outPDB(parm, self.get_fn('outPDB1.pdb', written=True)).execute()
        f = PDBFile.parse(self.get_fn('outPDB1.pdb', written=True))
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

    def test_out_cif(self):
        """ Test the outCIF action on AmberParm """
        parm = copy(gasparm)
        PT.loadRestrt(parm, self.get_fn('trx.inpcrd')).execute()
        PT.outCIF(parm, self.get_fn('outPDB1.cif', written=True)).execute()
        f = CIFFile.parse(self.get_fn('outPDB1.cif', written=True))
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

    def test_ti_merge1(self):
        """ Tests the tiMerge action on AmberParm with gas-phase parm """
        timask1re = re.compile(r'''timask1 *= *["']([\d,@:]+)["']''', re.I)
        timask2re = re.compile(r'''timask2 *= *["']([\d,@:]+)["']''', re.I)
        scmask1re = re.compile(r'''scmask1 *= *["']([\d,@:]+)["']''', re.I)
        scmask2re = re.compile(r'''scmask2 *= *["']([\d,@:]+)["']''', re.I)
        def get_masks(info):
            return dict(timask1=timask1re.findall(info)[0],
                        timask2=timask2re.findall(info)[0],
                        scmask1=scmask1re.findall(info)[0],
                        scmask2=scmask2re.findall(info)[0])

        output = StringIO()
        PT.tiMerge.output = output
        parm = AmberParm(self.get_fn('ava_aaa.parm7'), self.get_fn('ava_aaa.rst7'))
        act = PT.tiMerge(parm, ':1-5', ':6-10', ':3', ':8')
        str(act)
        # Make sure topology objects are the same
        a = parm[0]
        act.execute()
        self.assertIs(a, parm[0])
        output.seek(0)
        info = output.read()
        masks = get_masks(info)
        # Make sure the timasks are what we expect them to be
        timask1 = AmberMask(parm, masks['timask1'])
        timask2 = AmberMask(parm, masks['timask2'])
        scmask1 = AmberMask(parm, masks['scmask1'])
        scmask2 = AmberMask(parm, masks['scmask2'])
        self.assertEqual(timask1.Selection(),
                         [int(a.residue.idx == 2) for a in parm.atoms])
        self.assertEqual(timask2.Selection(),
                         [int(a.residue.idx == 5) for a in parm.atoms])
        self.assertEqual(scmask1.Selection(),
                         [int(a.residue.idx == 2) for a in parm.atoms])
        self.assertEqual(scmask2.Selection(),
                         [int(a.residue.idx == 5) for a in parm.atoms])
        self.assertEqual(len(parm.residues), 6) # Kept first 5 and mutant ALA

        output.truncate() # Reset the buffer for the next test

        # Check some error processing
        self.assertRaises(exc.TiMergeError, lambda:
                PT.tiMerge(load_file(self.get_fn('ava_aaa.parm7'), self.get_fn('ava_aaa.rst7')),
                           ':1-5', ':6-10', ':8', ':8').execute())
        self.assertRaises(exc.TiMergeError, lambda:
                PT.tiMerge(load_file(self.get_fn('ava_aaa.parm7'), self.get_fn('ava_aaa.rst7')),
                           ':1-5', ':6-10', ':3', ':3').execute())
        self.assertRaises(exc.TiMergeError, lambda:
                PT.tiMerge(load_file(self.get_fn('ava_aaa.parm7'), self.get_fn('ava_aaa.rst7')),
                           ':1-5', ':5-10', ':3', ':8').execute())
        self.assertRaises(exc.TiMergeError,
                PT.tiMerge(load_file(self.get_fn('ava_aaa.parm7')), ':1-5', ':6-10',
                    ':3', ':8').execute)
        parm = load_file(self.get_fn('ava_aaa.parm7'), self.get_fn('ava_aaa.rst7'))
        parm.residues[5][0].xx += 1 # Move away so it's not recognized as same
        self.assertRaises(exc.TiMergeError, lambda:
                PT.tiMerge(parm, ':1-5', ':6-10', ':3', ':8').execute())

    def test_ti_merge2(self):
        """ Tests the tiMerge action on AmberParm with solvated parm """
        parm = AmberParm(self.get_fn('ava_aaa.solv.parm7'), self.get_fn('ava_aaa.solv.rst7'))
        timask1re = re.compile(r'''timask1 *= *["']([\d,@:]+)["']''', re.I)
        timask2re = re.compile(r'''timask2 *= *["']([\d,@:]+)["']''', re.I)
        scmask1re = re.compile(r'''scmask1 *= *["']([\d,@:]+)["']''', re.I)
        scmask2re = re.compile(r'''scmask2 *= *["']([\d,@:]+)["']''', re.I)
        def get_masks(info):
            return dict(timask1=timask1re.findall(info)[0],
                        timask2=timask2re.findall(info)[0],
                        scmask1=scmask1re.findall(info)[0],
                        scmask2=scmask2re.findall(info)[0])

        output = StringIO()
        PT.tiMerge.output = output
        act = PT.tiMerge(parm, ':1-5', ':6-10', ':3', ':8', ':11', ':12')
        str(act)
        # Make sure topology objects are the same
        a = parm[0]
        act.execute()
        self.assertIs(a, parm[0])
        output.seek(0)
        info = output.read()
        masks = get_masks(info)
        # Make sure the timasks are what we expect them to be
        timask1 = AmberMask(parm, masks['timask1'])
        timask2 = AmberMask(parm, masks['timask2'])
        scmask1 = AmberMask(parm, masks['scmask1'])
        scmask2 = AmberMask(parm, masks['scmask2'])
        self.assertEqual(timask1.Selection(),
                         [int(a.residue.idx in (2, 6)) for a in parm.atoms])
        self.assertEqual(timask2.Selection(),
                         [int(a.residue.idx in (5, 7)) for a in parm.atoms])
        self.assertEqual(scmask1.Selection(),
                         [int(a.residue.idx in (2, 6)) for a in parm.atoms])
        self.assertEqual(scmask2.Selection(),
                         [int(a.residue.idx in (5, 7)) for a in parm.atoms])

        output.truncate() # Reset the buffer for the next test
        # Error handling checking
        parm = AmberParm(self.get_fn('ava_aaa.solv.parm7'), self.get_fn('ava_aaa.solv.rst7'))
        self.assertRaises(exc.TiMergeError, lambda:
                PT.tiMerge(parm, ':1-5', ':11', ':3', ':11').execute())

    def test_add12_6_4(self):
        """ Test the add12_6_4 action on AmberParm """
        parm = AmberParm(self.get_fn('Mg_ti1_b.parm7'))
        PT.addLJType(parm, '@14').execute()
        PT.changeLJPair(parm, '@14', ':MG', 3.26, 0.061666).execute()
        act = PT.add12_6_4(parm, ':MG', watermodel='TIP4PEW',
                     polfile=self.get_fn('lj_1264_pol.dat'))
        act.execute()
        str(act)
        parm.write_parm(self.get_fn('Mg_ti1_b_1264.parm7', written=True))
        self.assertTrue(diff_files(self.get_fn('Mg_ti1_b_1264.parm7', written=True),
                                   get_saved_fn('Mg_ti1_b_1264.parm7'))
        )
        # Error handling
        self.assertRaises(exc.LJ12_6_4Error, lambda:
                PT.add12_6_4(parm, ':MG', watermodel='FOO',
                     polfile=self.get_fn('lj_1264_pol.dat')).execute()
        )
        self.assertRaises(exc.LJ12_6_4Error, lambda:
                PT.add12_6_4(parm, ':MG', watermodel='FOO',
                    c4file='BAR').execute()
        )

    def test_add12_6_4_c4file(self):
        """ Test add12_6_4 action on AmberParm specifying c4file """
        from parmed.tools.add1264 import DEFAULT_C4_PARAMS
        fn = self.get_fn('c4file', written=True)
        with open(fn, 'w') as f:
            for items in iteritems(DEFAULT_C4_PARAMS['TIP4PEW']):
                f.write('%s %s\n' % items)
        parm = AmberParm(self.get_fn('Mg_ti1_b.parm7'))
        PT.addLJType(parm, '@14').execute()
        PT.changeLJPair(parm, '@14', ':MG', 3.26, 0.061666).execute()
        act = PT.add12_6_4(parm, ':MG', c4file=fn,
                     polfile=self.get_fn('lj_1264_pol.dat'))
        act.execute()
        str(act)
        parm.write_parm(self.get_fn('Mg_ti1_b_1264.parm7', written=True))
        self.assertTrue(diff_files(self.get_fn('Mg_ti1_b_1264.parm7', written=True),
                                   get_saved_fn('Mg_ti1_b_1264.parm7'))
        )

    def test_add_12_6_4_2metals(self):
        """ Test the add12_6_4 action on AmberParm with 2+ metals """
        parm1 = AmberParm(self.get_fn('mg_na_cl.parm7'))
        parm2 = AmberParm(self.get_fn('na_cl_mg.parm7'))
        PT.add12_6_4(parm1, ':MG,NA,CL',
                     polfile=self.get_fn('lj_1264_pol.dat')).execute()
        PT.add12_6_4(parm2, ':MG,NA,CL', watermodel='TIP3P',
                     polfile=self.get_fn('lj_1264_pol.dat')).execute()
        self.assertEqual(str(PT.printLJMatrix(parm1, ':MG')),
                         saved.PRINTLJMATRIX_MGNACL)
        self.assertEqual(str(PT.printLJMatrix(parm2, ':MG')),
                         saved.PRINTLJMATRIX_NACLMG)

    @unittest.skipIf(PYPY, 'NetCDF support does not work on PYPY yet')
    def test_write_coordinates(self):
        """ Test writeCoordinates method """
        parm = copy(gasparm)
        PT.loadCoordinates(parm, self.get_fn('trx.inpcrd')).execute()
        basefn = self.get_fn('test', written=True)
        # NetCDF trajectory
        act = PT.writeCoordinates(parm, basefn + '.nc')
        act.execute()
        str(act) # Make sure it doesn't fail
        self.assertTrue(pmd.amber.NetCDFTraj.id_format(basefn + '.nc'))
        PT.writeCoordinates(parm, basefn + '_nc', 'netcdftraj').execute()
        self.assertTrue(pmd.amber.NetCDFTraj.id_format(basefn + '_nc'))
        # NetCDF restart
        PT.writeCoordinates(parm, basefn + '.ncrst').execute()
        self.assertTrue(pmd.amber.NetCDFRestart.id_format(basefn + '.ncrst'))
        PT.writeCoordinates(parm, basefn + '_ncrst', 'netcdf').execute()
        self.assertTrue(pmd.amber.NetCDFRestart.id_format(basefn + '_ncrst'))
        # PDB
        PT.writeCoordinates(parm, basefn + '.pdb').execute()
        self.assertTrue(pmd.formats.PDBFile.id_format(basefn + '.pdb'))
        PT.writeCoordinates(parm, basefn + '_pdb', 'pdb').execute()
        self.assertTrue(pmd.formats.PDBFile.id_format(basefn + '_pdb'))
        # CIF
        PT.writeCoordinates(parm, basefn + '.cif').execute()
        self.assertTrue(pmd.formats.CIFFile.id_format(basefn + '.cif'))
        PT.writeCoordinates(parm, basefn + '_cif', 'CIF').execute()
        self.assertTrue(pmd.formats.CIFFile.id_format(basefn + '_cif'))
        # ASCII restart
        PT.writeCoordinates(parm, basefn + '.rst7').execute()
        self.assertTrue(pmd.amber.AmberAsciiRestart.id_format(basefn + '.rst7'))
        PT.writeCoordinates(parm, basefn + '_rst7', 'restart').execute()
        self.assertTrue(pmd.amber.AmberAsciiRestart.id_format(basefn + '_rst7'))
        # ASCII trajectory
        PT.writeCoordinates(parm, basefn + '.mdcrd').execute()
        self.assertTrue(pmd.amber.AmberMdcrd.id_format(basefn + '.mdcrd'))
        PT.writeCoordinates(parm, basefn + '_mdcrd', 'mdcrd').execute()
        self.assertTrue(pmd.amber.AmberMdcrd.id_format(basefn + '_mdcrd'))
        # Mol2 file
        PT.writeCoordinates(parm, basefn + '.mol2').execute()
        self.assertTrue(pmd.formats.Mol2File.id_format(basefn + '.mol2'))
        PT.writeCoordinates(parm, basefn + '_mol2', 'mol2').execute()
        self.assertTrue(pmd.formats.Mol2File.id_format(basefn + '_mol2'))
        # ASCII restart by default
        PT.writeCoordinates(parm, basefn + '.noext').execute()
        self.assertTrue(pmd.amber.AmberAsciiRestart.id_format(basefn + '.noext'))
        self.assertRaises(exc.InputError, lambda:
                PT.writeCoordinates(parm, basefn, 'noformat').execute())
        # Make sure overwriting is correctly handled
        self.assertRaises(exc.FileExists, act.execute)

    def test_check_validity(self):
        """ Tests the checkValidity test """
        act = PT.checkValidity(gasparm)
        str(act)
        act.execute()

class TestChamberParmActions(FileIOTestCase, TestCaseRelative):
    """ Tests actions on Amber prmtop files """

    def test_parmout_outparm_load_restrt(self):
        """ Test parmout, outparm, and loadCoordinates actions for ChamberParm """
        self._empty_writes()
        parm = copy(gascham)
        act = PT.loadCoordinates(parm, self.get_fn('ala_ala_ala.rst7'))
        act.execute()
        str(act) # Make sure it doesn't fail
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))
        PT.parmout(parm, self.get_fn('test.parm7', written=True)).execute()
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 1)
        self.assertTrue(diff_files(get_saved_fn('ala_ala_ala.parm7'),
                                   self.get_fn('test.parm7', written=True),
                                   absolute_error=1e-6))
        self._empty_writes()
        PT.parmout(parm, self.get_fn('test.parm7', written=True),
                         self.get_fn('test.rst7', written=True)).execute()
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 2)
        self.assertTrue(diff_files(get_saved_fn('ala_ala_ala.parm7'),
                                   self.get_fn('test.parm7', written=True),
                                   absolute_error=1e-6))
        self.assertTrue(diff_files(self.get_fn('ala_ala_ala.rst7'),
                                   self.get_fn('test.rst7', written=True),
                                   absolute_error=0.0001))
        self._empty_writes()
        PT.outparm(parm, self.get_fn('test.parm7', written=True)).execute()
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 1)
        self.assertTrue(diff_files(get_saved_fn('ala_ala_ala.parm7'),
                                   self.get_fn('test.parm7', written=True),
                                   absolute_error=1e-6))
        self._empty_writes()
        PT.outparm(parm, self.get_fn('test.parm7', written=True),
                         self.get_fn('test.rst7', written=True)).execute()
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 2)
        self.assertTrue(diff_files(get_saved_fn('ala_ala_ala.parm7'),
                                   self.get_fn('test.parm7', written=True),
                                   absolute_error=1e-6))
        self.assertTrue(diff_files(self.get_fn('ala_ala_ala.rst7'),
                                   self.get_fn('test.rst7', written=True),
                                   absolute_error=0.0001))
        # Check loadCoordinates error handling
        self.assertRaises(exc.ParmError, lambda:
                PT.loadCoordinates(parm, self.get_fn('trx.prmtop')).execute())

    def test_write_frcmod(self):
        """ Check that writeFrcmod fails for ChamberParm """
        parm = gascham
        self.assertRaises(exc.ParmError, lambda:
                PT.writeFrcmod(parm, self.get_fn('x', written=True)).execute())

    def test_write_off_load_restrt(self):
        """ Check that writeOFF fails for ChamberParm """
        parm = copy(gascham)
        PT.loadRestrt(parm, self.get_fn('ala_ala_ala.rst7')).execute()
        self.assertRaises(exc.ParmError, lambda:
                PT.writeOFF(parm, self.get_fn('test.off', written=True)).execute())

    def test_ti_merge(self):
        """ Check that tiMerge joins CHARMM-specific terms in ChamberParm """
        timask1re = re.compile(r'''timask1 *= *["']([\d,@:]+)["']''', re.I)
        timask2re = re.compile(r'''timask2 *= *["']([\d,@:]+)["']''', re.I)
        scmask1re = re.compile(r'''scmask1 *= *["']([\d,@:]+)["']''', re.I)
        scmask2re = re.compile(r'''scmask2 *= *["']([\d,@:]+)["']''', re.I)
        def get_masks(info):
            return dict(timask1=timask1re.findall(info)[0],
                        timask2=timask2re.findall(info)[0],
                        scmask1=scmask1re.findall(info)[0],
                        scmask2=scmask2re.findall(info)[0])
        # First convert from CHARMM-GUI PSF file
        a = PT.chamber(parmlist.ParmList(), '-psf', self.get_fn('ava_aaa.psf'),
                '-top', self.get_fn('top_all36_prot.rtf'), '-param',
                self.get_fn('par_all36_prot.prm'), '-crd', self.get_fn('ava_full.pdb'))
        a.execute()
        PT.tiMerge(a.parm, ':1-3', ':4-6', ':2', ':5').execute()
        self.assertEqual(len(a.parm.residues), 4) # Removed 2 redundant res.
        # Error checking
        parm = copy(gascham)
        PT.loadRestrt(parm, self.get_fn('ala_ala_ala.rst7')).execute()
        self.assertRaises(exc.ParmError, lambda:
                PT.tiMerge(parm, ':1-3', ':4-6', ':2', ':5').execute())

    def test_change_radii(self):
        """ Test changeRadii for ChamberParm """
        parm = copy(gascham)
        PT.changeRadii(parm, 'amber6').execute()
        self.assertEqual(parm.parm_data['RADIUS_SET'][0],
                         'amber6 modified Bondi radii (amber6)')
        for i, atom in enumerate(parm.atoms):
            radii, atomic_number = atom.solvent_radius, atom.atomic_number
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
            radii, atomic_number = atom.solvent_radius, atom.atomic_number
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
            radii, atomic_number = atom.solvent_radius, atom.atomic_number
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
            radii, atomic_number = atom.solvent_radius, atom.atomic_number
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
            radii, atomic_number = atom.solvent_radius, atom.atomic_number
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

    def test_change_lj_pair(self):
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

    def test_change_lj_14_pair(self):
        """ Test changeLJ14Pair for ChamberParm """
        parm = copy(gascham)
        act = PT.changeLJ14Pair(parm, '@%NH3', '@%HC', 1.0, 1.0)
        act.execute()
        str(act)
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
        # Check handling with no atoms
        PT.changeLJ14Pair(parm, '@NONE', '@%H', 1.0, 1.0).execute()
        # Make sure nothing changed
        np.testing.assert_equal(acoef, parm.parm_data['LENNARD_JONES_14_ACOEF'])
        np.testing.assert_equal(bcoef, parm.parm_data['LENNARD_JONES_14_BCOEF'])
        # Check error handling
        self.assertRaises(exc.ChangeLJPairError, lambda:
                PT.changeLJ14Pair(parm, ':*', '@%H', 1.0, 1.0).execute())
        self.assertRaises(exc.ChangeLJPairError, lambda:
                PT.changeLJ14Pair(parm, '@%H', ':*', 1.0, 1.0).execute())

    def test_change(self):
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
            self.assertEqual(atom.solvent_radius, parm.parm_data['RADII'][i])
            if atom.residue.idx < 2:
                self.assertEqual(atom.solvent_radius, 2.0)
            else:
                self.assertEqual(atom.solvent_radius, gascham.parm_data['RADII'][i])
        PT.change(parm, 'SCREEN', '*', 0.0).execute()
        for i, atom in enumerate(parm.atoms):
            self.assertEqual(atom.screen, parm.parm_data['SCREEN'][i])
            self.assertEqual(atom.screen, 0.0)
        # Check bad input
        self.assertRaises(exc.ParmedChangeError, lambda:
                          PT.change(parm, 'RESIDUE_LABEL', ':*', 'NaN'))

    def test_print_info(self):
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
        # Check this on non-AmberParm Structure classes
        self.assertRaises(exc.ParmError, lambda:
                PT.printInfo(load_file(self.get_fn('2koc.pdb'))))

    def test_add_change_lj_type(self):
        """ Test addLJType and changeLJSingleType on ChamberParm """
        parm = copy(gascham)
        act = PT.addLJType(parm, '@1')
        act.execute()
        str(act)
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
        act = PT.addLJType(parm2, '@1', radius=1.0, epsilon=1.0)
        act.execute()
        str(act)
        act = PT.changeLJSingleType(parm, '@1', 1.0, 1.0)
        act.execute()
        str(act)
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
        # Make sure an empty selection still works
        act = PT.changeLJSingleType(parm, '@NOATOM', 1, 2)
        act.execute()
        str(act)
        # Turn the first atom into a new type, then make sure changeLJSingleType
        # fails if multiple atom types are selected
        PT.addLJType(parm, '@1', radius=1, epsilon=0.5).execute()
        self.assertRaises(exc.LJ_TypeError, PT.changeLJSingleType(parm, '@1-2',
                                        1, 2).execute)

    def test_print_lj_types(self):
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

    def test_scee_scnb(self):
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

    def test_print_details(self):
        """ Test printDetails for ChamberParm """
        act = PT.printDetails(gascham, '@1')
        self.assertEqual(str(act), saved.PRINT_DETAILSC)

    def test_print_flags(self):
        """ Test printFlags for ChamberParm """
        act = PT.printFlags(gascham)
        printed_flags = set()
        for line in str(act).split('\n'):
            if line.startswith('%FLAG'):
                printed_flags.add(line.split()[1])
        self.assertEqual(printed_flags, set(gascham.parm_data.keys()))

    def test_print_pointers(self):
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

    def test_print_bonds(self):
        """ Test printBonds for ChamberParm """
        act = PT.printBonds(gascham, '@1')
        self.assertEqual(str(act), saved.PRINT_BONDSC)

    def test_print_angles(self):
        """ Test printAngles for ChamberParm """
        act = PT.printAngles(gascham, '@1')
        self.assertEqual(str(act), saved.PRINT_ANGLESC)

    def test_print_dihedrals(self):
        """ Test printDihedrals for ChamberParm """
        act = PT.printDihedrals(gascham, '@1')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALSC)

    def test_set_molecules(self):
        """ Test setMolecules for ChamberParm """
        parm = copy(solvchamber)
        atoms = [atom for atom in parm.atoms] # shallow copy!
        self.assertTrue(all([x is y for x,y in zip(parm.atoms,atoms)]))
        self.assertEqual(parm.ptr('IPTRES'), 3)
        self.assertEqual(parm.ptr('NSPM'), 942)
        self.assertEqual(parm.ptr('NSPSOL'), 2)
        # To keep the output clean
        PT.setMolecules(parm).execute()
        self.assertTrue(all([x is y for x,y in zip(parm.atoms, atoms)]))
        # Now check that setMolecules can apply another time.
        PT.setMolecules(parm, solute_ions=False).execute()

    def test_net_charge(self):
        """ Test netCharge for ChamberParm """
        act = PT.netCharge(gascham)
        chg = act.execute() # check this part of the API
        self.assertEqual(str(act), 'The net charge of :* is %.4f' % chg)
        self.assertAlmostEqual(chg, 0.0, places=6)
        chg = PT.netCharge(gascham, ':1').execute()
        self.assertAlmostEqual(chg, 1.0, places=6)
        chg = PT.netCharge(gascham, ':3').execute()
        self.assertAlmostEqual(chg, -1.0, places=6)

    def test_strip(self):
        """ Test strip action for ChamberParm """
        parm = copy(gascham)
        act = PT.strip(parm, ':1')
        str(act)
        act.execute()
        self.assertEqual(parm.ptr('natom'), 21)
        self.assertEqual(len(parm.atoms), 21)
        # Good enough for here. The strip action is repeatedly tested in the
        # core Amber test suite as part of the MM/PBSA tests via ante-MMPBSA.py
        # and that part also tests that the energies come out correct as well

    def test_define_solvent(self):
        """ Test defineSolvent for ChamberParm """
        import parmed.residue as residue
        PT.defineSolvent(gascham, 'WAT,HOH,Na+,Cl-').execute()
        self.assertEqual(residue.SOLVENT_NAMES, 'WAT HOH Na+ Cl-'.split())
        PT.defineSolvent(gascham, 'WAT,HOH').execute()
        self.assertEqual(residue.SOLVENT_NAMES, 'WAT HOH'.split())

    def test_add_exclusions(self):
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

    def test_add_delete_dihedral(self):
        """ Test the addDihedral and deleteDihedral actions for ChamberParm """
        parm = copy(gascham)
        n = PT.deleteDihedral(parm, *':ALA@N :ALA@CA :ALA@CB :ALA@HB1'.split()).execute()
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

    def test_set_bond(self):
        """ Test setBond for ChamberParm """
        parm = copy(gascham)
        PT.setBond(parm, ':ALA@CA', ':ALA@CB', 300.0, 1.5).execute()
        act = PT.printBonds(parm, ':ALA@CA')
        self.assertEqual(str(act), saved.SET_BONDC)

    def test_set_angle(self):
        """ Test setAngle for ChamberParm """
        parm = copy(gascham)
        PT.setAngle(parm, ':ALA@CA', ':ALA@CB', ':ALA@HB1', 40, 100).execute()
        act = PT.printAngles(parm, ':ALA@CB')
        self.assertEqual(str(act), saved.SET_ANGLEC)

    def test_add_atomic_number(self):
        """ Test addAtomicNumber for ChamberParm """
        parm = copy(gascham)
        self.assertFalse('ATOMIC_NUMBER' in parm.parm_data)
        atomic_numbers = [atom.atomic_number for atom in parm.atoms]
        PT.addAtomicNumber(parm).execute()
        self.assertEqual(parm.parm_data['ATOMIC_NUMBER'], atomic_numbers)

    def test_print_lj_matrix(self):
        """ Test printLJMatrix for ChamberParm """
        act = PT.printLJMatrix(gascham, '@1')
        self.assertEqual(str(act), saved.PRINT_LJMATRIXC)

    def test_delete_bond(self):
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
        str(act)
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

    def test_delete_bond2(self):
        """ Test deleteBond with different parameter types """
        parm = copy(gascham)
        # Pick the bond we plan to delete, pick out every angle and dihedral
        # that contains that bond, and then delete it. Then make sure none of
        # the valence terms that contained that bond remain afterwards. We
        # already have a test to make sure that the __contains__ method works
        # for atoms and bonds.
        bond = parm.atoms[1].bonds[0]
        deleted_angles = list()
        deleted_dihedrals = list()
        deleted_impropers = list()
        deleted_urey_bradleys = list()
        deleted_cmaps = list()
        for angle in parm.angles:
            if bond in angle: deleted_angles.append(angle)
        for dihedral in parm.dihedrals:
            if bond in dihedral: deleted_dihedrals.append(dihedral)
            # Move these to the r-b torsion list
            parm.rb_torsions.append(dihedral)
        # Remove them from the dihedral list
        for i in reversed(range(len(parm.dihedrals))):
            if parm.dihedrals[i] in parm.rb_torsions:
                del parm.dihedrals[i]
        for imp in parm.impropers:
            if bond in imp: deleted_impropers.append(imp)
        for ub in parm.urey_bradleys:
            if bond in ub: deleted_urey_bradleys.append(ub)
        for cmap in parm.cmaps:
            if bond in cmap: deleted_cmaps.append(cmap)
        act = PT.deleteBond(parm, '@1', '@2', 'verbose')
        str(act)
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

    def test_summary(self):
        """ Test summary action for ChamberParm """
        parm = load_file(self.get_fn('dhfr_cmap_pbc.parm7'))
        parm.load_rst7(self.get_fn('dhfr_cmap_pbc.rst7'))
        act = PT.summary(parm)
        self.assertTrue(detailed_diff(str(act), saved.SUMMARYC1, relative_error=1e-6))
        act = PT.summary(load_file(self.get_fn('2koc.pdb')))
        repr(act)

    def test_scale(self):
        """ Test scale action for ChamberParm """
        parm = copy(gascham)
        act = PT.scale(parm, 'DIHEDRAL_FORCE_CONSTANT', 2.0)
        act.execute()
        str(act)
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
        # Error handling
        self.assertRaises(exc.ArgumentError, lambda:
                PT.scale(parm, 'NOTAFLAG', 10.0).execute())
        self.assertRaises(exc.ArgumentError, lambda:
                PT.scale(parm, 'ATOM_NAME', 10.0).execute())

    def test_interpolate(self):
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
                       prefix=self.get_fn('test.parm7', written=True)).execute()
        sys.stdout = sys.__stdout__
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 5)
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.2', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.3', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.4', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.5', written=True)))
        self.assertTrue(os.path.exists(self.get_fn('test.parm7.6', written=True)))
        # Now check them all
        ladder = [origparm]
        ladder.append(AmberParm(self.get_fn('test.parm7.2', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.3', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.4', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.5', written=True)))
        ladder.append(AmberParm(self.get_fn('test.parm7.6', written=True)))
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

    def test_change_prot_state(self):
        """ Check that changeProtState fails for ChamberParm """
        parm = copy(solvchamber)
        self.assertRaises(exc.ParmError, lambda:
                PT.changeProtState(parm, ':32', 0).execute())

    def test_lmod(self):
        """ Test lmod action for ChamberParm """
        parm = copy(gascham)
        parm.parm_data['LENNARD_JONES_ACOEF'][3] = 0.0
        self.assertFalse(all(parm.parm_data['LENNARD_JONES_ACOEF']))
        PT.lmod(parm).execute()
        self.assertTrue(all(parm.parm_data['LENNARD_JONES_ACOEF']))

    def test_add_delete_pdb(self):
        """ Test addPDB and deletePDB actions for ChamberParm """
        parm = copy(gascham)
        PT.addPDB(parm, self.get_fn('ala_ala_ala.pdb'), 'elem', 'allicodes').execute()
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
        PT.addPDB(parm, self.get_fn('ala_ala_ala.pdb')).execute()
        self.assertFalse('RESIDUE_ICODE' in parm.flag_list)
        self.assertFalse('ATOM_ELEMENT' in parm.flag_list)
        self.assertTrue('RESIDUE_NUMBER' in parm.flag_list)
        self.assertTrue('RESIDUE_CHAINID' in parm.flag_list)

    def test_h_mass_repartition(self):
        """ Test HMassRepartition action for ChamberParm """
        parm = copy(solvchamber)
        PT.defineSolvent(parm, 'TIP3').execute()
        act = PT.HMassRepartition(parm, 2.0)
        act.execute()
        str(act)
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
        act = PT.HMassRepartition(parm, 3.0, 'dowater')
        act.execute()
        str(act)
        for atom in parm.atoms:
            if atom.atomic_number == 1:
                self.assertEqual(atom.mass, 3.0)
        self.assertAlmostEqual(sum(solvchamber.parm_data['MASS']),
                               sum(parm.parm_data['MASS']), places=6)

    def test_out_pdb(self):
        """ Test the outPDB action on ChamberParm """
        parm = copy(gascham)
        PT.loadRestrt(parm, self.get_fn('ala_ala_ala.rst7')).execute()
        PT.outPDB(parm, self.get_fn('outPDB1.pdb', written=True)).execute()
        f = PDBFile.parse(self.get_fn('outPDB1.pdb', written=True))
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

    def test_out_cif(self):
        """ Test the outCIF action on ChamberParm """
        parm = copy(gascham)
        PT.loadRestrt(parm, self.get_fn('ala_ala_ala.rst7')).execute()
        PT.outCIF(parm, self.get_fn('outPDB1.cif', written=True)).execute()
        f = CIFFile.parse(self.get_fn('outPDB1.cif', written=True))
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

class TestAmoebaParmActions(FileIOTestCase, TestCaseRelative):
    """ Tests actions on Amber prmtop files """

    def test_parmout_outparm_load_restrt(self):
        """ Test parmout, outparm, and loadRestrt actions on AmoebaParm """
        self._empty_writes()
        parm = copy(amoebaparm)
        PT.loadRestrt(parm, self.get_fn('nma.rst7')).execute()
        for atom in parm.atoms:
            self.assertTrue(hasattr(atom, 'xx'))
            self.assertTrue(hasattr(atom, 'xy'))
            self.assertTrue(hasattr(atom, 'xz'))
        PT.parmout(parm, self.get_fn('test.parm7', written=True)).execute()
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 1)
        self.assertTrue(diff_files(self.get_fn('nma.parm7'),
                                   self.get_fn('test.parm7', written=True)))
        self._empty_writes()
        PT.parmout(parm, self.get_fn('test.parm7', written=True),
                         self.get_fn('test.rst7', written=True)).execute()
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 2)
        self.assertTrue(diff_files(self.get_fn('nma.parm7'),
                                   self.get_fn('test.parm7', written=True)))
        self.assertTrue(diff_files(self.get_fn('nma.rst7'),
                                   self.get_fn('test.rst7', written=True),
                                   absolute_error=0.0001))
        self._empty_writes()
        PT.outparm(parm, self.get_fn('test.parm7', written=True)).execute()
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 1)
        self.assertTrue(diff_files(self.get_fn('nma.parm7'),
                                   self.get_fn('test.parm7', written=True)))
        self._empty_writes()
        PT.outparm(parm, self.get_fn('test.parm7', written=True),
                         self.get_fn('test.rst7', written=True)).execute()
        self.assertEqual(len(os.listdir(self._temporary_directory.name)), 2)
        self.assertTrue(diff_files(self.get_fn('nma.parm7'),
                                   self.get_fn('test.parm7', written=True)))
        self.assertTrue(diff_files(self.get_fn('nma.rst7'),
                                   self.get_fn('test.rst7', written=True),
                                   absolute_error=0.0001))

    def test_write_frcmod(self):
        """ Check that writeFrcmod fails for AmoebaParm """
        with self.assertRaises(exc.ParmError):
            PT.writeFrcmod(amoebaparm, self.get_fn('x', written=True)).execute()

    def test_write_off_load_restrt(self):
        """ Check that writeOFF fails for AmoebaParm """
        parm = copy(amoebaparm)
        PT.loadRestrt(parm, self.get_fn('nma.rst7')).execute()
        self.assertRaises(exc.ParmError, lambda:
                PT.writeOFF(parm, self.get_fn('test.off', written=True)).execute())

    def test_change_radii(self):
        """ Check that changeRadii fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.changeRadii(parm, 'amber6').execute())

    def test_change_lj_pair(self):
        """ Check that changeLJPair fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.changeLJPair(parm, '@%NH3', '@%HC', 1.0, 1.0).execute())

    def test_change_lj_14_pair(self):
        """ Check that changeLJ14Pair fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.changeLJ14Pair(parm, '@%NH3', '@%HC', 1.0, 1.0).execute())

    def test_change(self):
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
        # Change atomic number
        PT.change(parm, 'ATOMIC_NUMBER', '@1', 10).execute()
        self.assertEqual(parm.atoms[0].atomic_number, 10)
        self.assertEqual(parm.parm_data['AMOEBA_ATOMIC_NUMBER'][0], 10)
        # Now make sure it adds the AMOEBA_ATOMIC_NUMBER section
        parm.delete_flag('AMOEBA_ATOMIC_NUMBER')
        PT.change(parm, 'ATOMIC_NUMBER', '@1', 1).execute()
        self.assertIn('AMOEBA_ATOMIC_NUMBER', parm.parm_data)
        self.assertEqual(parm[0].atomic_number, 1)
        self.assertEqual(parm.parm_data['AMOEBA_ATOMIC_NUMBER'][0], 1)
        # Check some error-handling
        self.assertRaises(exc.ParmedChangeError, lambda:
                PT.change(parm, 'ATOM_TYPE_INDEX', '@1', 4).execute())
        self.assertRaises(exc.ParmedChangeError, lambda:
                PT.change(parm, 'RADII', ':1-2', 2.0, 'quiet').execute())
        self.assertRaises(exc.ParmedChangeError, lambda:
                PT.change(parm, 'SCREEN', '*', 0.0).execute())
        # Check bad input
        self.assertRaises(exc.ParmedChangeError, lambda:
                          PT.change(parm, 'RESIDUE_LABEL', ':*', 'NaN'))
        # Make sure string casting always works
        act = PT.change(parm, 'MASS', '@NOTHING', 0, 'quiet')
        str(act)
        previous_masses = parm.parm_data['MASS'][:]
        act.execute()
        np.testing.assert_equal(previous_masses, parm.parm_data['MASS'])

    def test_print_info(self):
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

    def test_add_change_lj_type(self):
        """ Check that addLJType and changeLJSingleType fail for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.addLJType(parm, '@1').execute())
        self.assertRaises(exc.ParmError, lambda:
                PT.changeLJSingleType(parm, '@1', 1.0, 1.0).execute())

    def test_print_lj_types(self):
        """ Check that printLJTypes fails for AmoebaParm """
        self.assertRaises(exc.ParmError, lambda:
                PT.printLJTypes(amoebaparm, '@1'))

    def test_scee_scnb(self):
        """ Check that scee and scnb fail for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda: PT.scee(parm, 10).execute())
        self.assertRaises(exc.ParmError, lambda: PT.scnb(parm, 10).execute())

    def test_print_details(self):
        """ Test printDetails for AmoebaParm """
        act = PT.printDetails(amoebaparm, ':1-2')
        self.assertEqual(str(act), saved.PRINT_DETAILSA)

    def test_print_flags(self):
        """ Test printFlags for AmoebaParm """
        act = PT.printFlags(amoebaparm)
        printed_flags = set()
        for line in str(act).split('\n'):
            if line.startswith('%FLAG'):
                printed_flags.add(line.split()[1])
        self.assertEqual(printed_flags, set(amoebaparm.parm_data.keys()))

    def test_print_pointers(self):
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

    def test_print_bonds(self):
        """ Test printBonds for AmoebaParm """
        act = PT.printBonds(amoebaparm, '@1')
        self.assertEqual(str(act), saved.PRINT_BONDSA)

    def test_print_angles(self):
        """ Test printAngles for AmoebaParm """
        act = PT.printAngles(amoebaparm, '@1')
        self.assertEqual(str(act), saved.PRINT_ANGLESA)

    def test_print_dihedrals(self):
        """ Test printDihedrals for AmoebaParm """
        act = PT.printDihedrals(amoebaparm, '@1')
        self.assertEqual(str(act), saved.PRINT_DIHEDRALSA)
        str(PT.printDihedrals(amoebaparm, ':*', '@10-12', ':*', ':*'))
        str(PT.printDihedrals(amoebaparm, ':*', ':1-20', ':*', ':*'))

    def test_set_molecules(self):
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

    def test_net_charge(self):
        """ Test netCharge for AmoebaParm (charge is the monopole) """
        act = PT.netCharge(amoebaparm)
        chg = act.execute() # check the netCharge.execute return value
        self.assertEqual(str(act), 'The net charge of :* is %.4f' % chg)
        self.assertAlmostEqual(chg, 0.0)
        chg = PT.netCharge(amoebaparm, ':WAT').execute()
        self.assertAlmostEqual(chg, 0)

    def test_strip(self):
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

    def test_define_solvent(self):
        """ Test defineSolvent for AmoebaParm """
        import parmed.residue as residue
        PT.defineSolvent(amoebaparm, 'WAT,HOH,Na+,Cl-').execute()
        self.assertEqual(residue.SOLVENT_NAMES, 'WAT HOH Na+ Cl-'.split())
        PT.defineSolvent(amoebaparm, 'WAT,HOH').execute()
        self.assertEqual(residue.SOLVENT_NAMES, 'WAT HOH'.split())

    def test_add_exclusions(self):
        """ Check that addExclusions fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.addExclusions(parm, ':*', ':*').execute())

    def test_add_delete_dihedral(self):
        """ Check that addDihedral and deleteDihedral fail for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
            PT.deleteDihedral(parm, *'@1 @2 @3 @4'.split()).execute())
        self.assertRaises(exc.ParmError, lambda:
            PT.addDihedral(parm, '@1 @2 @3 @4 0.1556 3 0 1.2 2.0', type='normal').execute()
        )

    def test_set_bond(self):
        """ Check that setBond fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
            PT.setBond(parm, ':ALA@CA', ':ALA@CB', 300.0, 1.5).execute())

    def test_set_angle(self):
        """ Check that setAngle fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
            PT.setAngle(parm, ':ALA@CA :ALA@CB :ALA@HB1 40 100').execute())

    def test_add_atomic_number(self):
        """ Test addAtomicNumber for AmoebaParm """
        parm = copy(amoebaparm)
        atomic_numbers = [atom.atomic_number for atom in parm.atoms]
        PT.addAtomicNumber(parm).execute()
        self.assertEqual(parm.parm_data['AMOEBA_ATOMIC_NUMBER'], atomic_numbers)

    def test_print_lj_matrix(self):
        """ Check that printLJMatrix fails for AmoebaParm """
        self.assertRaises(exc.ParmError, lambda:
                PT.printLJMatrix(amoebaparm, '@1'))

    def test_delete_bond(self):
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
        str(act)
        act.execute()
        self.assertTrue(bond not in parm.bonds)
        for attr in objs_with_bond:
            for obj in attr:
                self.assertNotIn(bond, attr)

    def test_summary(self):
        """ Test summary action for AmoebaParm """
        parm = copy(amoebaparm)
        act = PT.summary(parm)
        self.assertEqual(str(act), saved.SUMMARYA1)
        PT.loadRestrt(parm, self.get_fn('nma.rst7')).execute()
        act = PT.summary(parm)
        self.assertEqual(str(act), saved.SUMMARYA2)

    def test_scale(self):
        """ Test scale action for AmoebaParm """
        parm = copy(amoebaparm)
        flag = 'AMOEBA_STRETCH_BEND_FORCE_CONSTANT'
        PT.scale(parm, flag, 2.0).execute()
        self.assertEqual([2*x for x in amoebaparm.parm_data[flag]],
                         parm.parm_data[flag])
        PT.scale(parm, flag, 0.0).execute()
        for val in parm.parm_data[flag]:
            self.assertEqual(val, 0)

    def test_interpolate(self):
        """ Check that interpolate action fails for AmoebaParm """
        self.assertRaises(exc.ParmError, lambda:
                PT.interpolate(amoebaparm).execute())

    def test_change_prot_state(self):
        """ Check that changeProtState fails for AmoebaParm """
        parm = copy(amoebaparm)
        self.assertRaises(exc.ParmError, lambda:
                PT.changeProtState(parm, ':32', 0).execute())

    def test_lmod(self):
        """ Check that lmod fails for AmoebaParm """
        self.assertRaises(exc.ParmError, lambda: PT.lmod(amoebaparm).execute())

    def test_add_delete_pdb(self):
        """ Test addPDB and deletePDB for AmoebaParm """
        parm = copy(amoebaparm)
        PT.addPDB(parm, self.get_fn('nma.pdb'), 'elem', 'allicodes').execute()
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
        parm.write_parm(self.get_fn('amoeba_pdb.parm7', written=True))
        parm = AmoebaParm(self.get_fn('amoeba_pdb.parm7', written=True))
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
        PT.addPDB(parm, self.get_fn('nma.pdb')).execute()
        self.assertFalse('RESIDUE_ICODE' in parm.flag_list)
        self.assertFalse('ATOM_ELEMENT' in parm.flag_list)
        self.assertTrue('RESIDUE_NUMBER' in parm.flag_list)
        self.assertTrue('RESIDUE_CHAINID' in parm.flag_list)

    def test_h_mass_repartition(self):
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

    def test_out_pdb(self):
        """ Test the outPDB action on AmoebaParm """
        parm = copy(amoebaparm)
        PT.loadRestrt(parm, self.get_fn('nma.rst7')).execute()
        PT.outPDB(parm, self.get_fn('outPDB1.pdb', written=True)).execute()
        f = PDBFile.parse(self.get_fn('outPDB1.pdb', written=True))
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

    def test_out_cif(self):
        """ Test the outCIF action on AmoebaParm """
        parm = copy(amoebaparm)
        PT.loadRestrt(parm, self.get_fn('nma.rst7')).execute()
        PT.outCIF(parm, self.get_fn('outPDB1.cif', written=True)).execute()
        f = CIFFile.parse(self.get_fn('outPDB1.cif', written=True))
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

    def test_ti_merge(self):
        """ Check that tiMerge fails for AmoebaParm """
        parm = copy(amoebaparm)
        PT.loadRestrt(parm, self.get_fn('nma.rst7')).execute()
        self.assertRaises(exc.ParmError, lambda:
                PT.tiMerge(parm, ':1-3', ':4-6', ':2', ':5').execute())

class TestOtherParm(FileIOTestCase):
    """ Tests the use of other parms as the main parm """

    def test_summary(self):
        """ Tests the use of a PDB file with the summary action """
        parm = load_file(self.get_fn('4lzt.pdb'))
        self.assertEqual(str(PT.summary(parm)), saved.PDB_SUMMARY)

    def test_print_bonds(self):
        """ Tests printBonds on a PSF file with no parameters """
        parm = load_file(self.get_fn('ala_ala_ala.psf'))
        repr(PT.printBonds(parm, ':1-2'))

    def test_print_angles(self):
        """ Tests printAngles on a PSF file with no parameters """
        parm = load_file(self.get_fn('ala_ala_ala.psf'))
        repr(PT.printAngles(parm, ':1'))
        repr(PT.printAngles(parm, ':1', ':1'))
        repr(PT.printAngles(parm, ':1', ':1', ':1'))

    def test_print_dihedrals(self):
        """ Tests printDihedrals on a PSF file with no parameters """
        parm = load_file(self.get_fn('ala_ala_ala.psf'))
        repr(PT.printDihedrals(parm, ':1'))
        repr(PT.printDihedrals(parm, ':*', ':1'))
        repr(PT.printDihedrals(parm, ':*', ':*', ':1'))
        repr(PT.printDihedrals(parm, ':*', ':*', ':*', ':*'))
        parm.dihedrals[0].improper = True
        parm.dihedrals[1].ignore_end = True
        repr(PT.printDihedrals(parm, ':1'))
        repr(PT.printDihedrals(parm, ':*', ':1'))
        repr(PT.printDihedrals(parm, ':*', ':*', ':1'))
        repr(PT.printDihedrals(parm, ':*', ':*', ':*', ':*'))

    @unittest.skipUnless(HAS_GROMACS, 'Cannot test without GROMACS')
    def test_parm(self):
        """ Tests the parm action on a series of topology types and listParms """
        parms = parmlist.ParmList()
        # Make sure listParms works on an empty list
        repr(PT.listParms(parms))
        # First add an AmberParm
        act = PT.parm(parms, self.get_fn('ash.parm7'))
        act.execute()
        str(act)
        self.assertEqual(len(parms), 1)
        self.assertIs(parms.parm, parms[-1])
        # Next add a PDB and make sure it is the active parm
        PT.parm(parms, self.get_fn('2koc.pdb')).execute()
        self.assertEqual(len(parms), 2)
        self.assertIs(parms.parm, parms[-1])
        # Next add a GROMACS topology file
        PT.parm(parms, self.get_fn('ildn.solv.top')).execute()
        self.assertEqual(len(parms), 3)
        self.assertIs(parms.parm, parms[-1])
        # Next copy the amber parm
        act = PT.parm(parms, copy=0)
        act.execute()
        str(act)
        self.assertEqual(len(parms), 4)
        self.assertIs(parms[-1], parms.parm)
        self.assertEqual(len(parms[0].atoms), len(parms[3].atoms))
        # Next copy the GROMACS topology file parm
        act = PT.parm(parms, copy=self.get_fn('ildn.solv.top'))
        act.execute()
        str(act)
        self.assertEqual(len(parms), 5)
        self.assertIs(parms[-1], parms.parm)
        self.assertEqual(len(parms[2].atoms), len(parms[4].atoms))

        # Check setting new active parm by name and index
        act = PT.parm(parms, select=0)
        str(act)
        act.execute()
        self.assertEqual(len(parms), 5)
        self.assertIs(parms.parm, parms[0])
        act = PT.parm(parms, select=self.get_fn('2koc.pdb'))
        str(act)
        act.execute()
        self.assertEqual(len(parms), 5)
        self.assertIs(parms.parm, parms[1])

        # Error handling
        self.assertRaises(exc.ParmError, lambda: PT.parm(parms).execute())
        self.assertRaises(exc.ParmError, lambda:
                PT.parm(parms, copy=0, select=1).execute())
        self.assertRaises(exc.ParmError, lambda:
                PT.parm(parms, 'notafile').execute())
        act = PT.parm(parms, select='notafile')
        str(act)
        self.assertWarns(exc.SeriousParmWarning, act.execute)
        act = PT.parm(parms, select=100)
        str(act)
        self.assertWarns(exc.SeriousParmWarning, act.execute)
        act = PT.parm(parms, copy='notafile')
        str(act)
        self.assertWarns(exc.SeriousParmWarning, act.execute)
        act = PT.parm(parms, copy=100)
        str(act)
        self.assertWarns(exc.SeriousParmWarning, act.execute)

        # Catch permissions error, but only on Linux
        import platform
        if platform.system() == 'Windows': return
        fn = self.get_fn('test.pdb', written=True)
        with open(fn, 'w') as fw, open(self.get_fn('ash.parm7')) as fr:
            fw.write(fr.read())
        # Change the permissions to remove read perms
        os.chmod(fn, int('333', 8))
        act = PT.parm(parms, fn)
        self.assertWarns(exc.SeriousParmWarning, act.execute)

        # Now check that listParms works
        info = repr(PT.listParms(parms))
        self.assertIn('(active)', info)

    def test_h_mass_repartition(self):
        """ Tests HMassRepartition on arbitrary Structure instances """
        from parmed import periodic_table
        S = periodic_table.AtomicNum['S']
        Sm = periodic_table.Mass['S']
        struct = pmd.Structure()
        struct.add_atom(pmd.Atom(name='H1', atomic_number=1, mass=1.001), 'H2', 1)
        struct.add_atom(pmd.Atom(name='H2', atomic_number=1, mass=1.001), 'H2', 1)
        struct.add_atom(pmd.Atom(name='H1', atomic_number=1, mass=1.001), 'HSH', 2)
        struct.add_atom(pmd.Atom(name='H2', atomic_number=1, mass=1.001), 'HSH', 2)
        struct.add_atom(pmd.Atom(name='S', atomic_number=S, mass=Sm), 'HSH', 2)
        struct.bonds.append(pmd.Bond(struct[0], struct[1]))
        struct.bonds.append(pmd.Bond(struct[2], struct[3]))
        struct.bonds.append(pmd.Bond(struct[2], struct[4]))
        struct.bonds.append(pmd.Bond(struct[3], struct[4]))
        total_mass = sum(a.mass for a in struct.atoms)

        with self.assertWarns(exc.ParmWarning):
            PT.HMassRepartition(struct).execute()
        PT.HMassRepartition(struct).execute()
        self.assertEqual(struct[0].mass, 1.001)
        self.assertEqual(struct[1].mass, 1.001)
        self.assertEqual(struct[2].mass, 3.024)
        self.assertEqual(struct[3].mass, 3.024)
        self.assertAlmostEqual(sum(a.mass for a in struct.atoms), total_mass)
        # Now what happens if we make our Hs *too* heavy?
        with self.assertRaises(exc.HMassRepartitionError):
            PT.HMassRepartition(struct, 100).execute()

    def test_delete_bond(self):
        """ Tests deleteBond on arbitrary Structure instances """
        from parmed import periodic_table
        struct = create_random_structure(parametrized=True)
        act = PT.deleteBond(struct, '@%d' % (struct.bonds[0].atom1.idx+1),
                '@%d' % (struct.bonds[0].atom2.idx+1))
        str(act)
        act.execute()
