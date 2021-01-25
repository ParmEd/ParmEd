"""
Tests parmed.formats package
"""
from __future__ import division
import utils

from copy import copy
import numpy as np
import parmed as pmd
from parmed import (amber, charmm, exceptions, formats, gromacs, residue, Structure, read_PDB, Atom,
                    read_CIF, download_PDB, download_CIF, topologyobjects, write_PDB, write_CIF)
from parmed.symmetry import Symmetry
from parmed.modeller import ResidueTemplate, ResidueTemplateContainer
from parmed.utils import PYPY
from parmed.utils.six import iteritems, add_metaclass
from parmed.utils.six.moves import zip, StringIO, range
import random
import os
import sys
import unittest
from utils import get_fn, diff_files, run_all_tests, is_jenkins, HAS_GROMACS, FileIOTestCase
import warnings

def reset_stringio(io):
    """ Resets a StringIO instance to "empty-file" state """
    io.seek(0)
    io.truncate()
    return io

try:
    import rdkit
    has_rdkit = True
except ImportError:
    has_rdkit = False

is_linux = sys.platform.startswith('linux')

class TestFileLoader(FileIOTestCase):
    """ Tests the automatic file loader """

    def test_load_blank_file(self):
        """ Makes sure that a blank file does not match any id_format """
        from parmed.formats.registry import PARSER_REGISTRY
        fn = self.get_fn('test', written=True)
        with open(fn, 'w'):
            pass
        for name, cls in iteritems(PARSER_REGISTRY):
            self.assertFalse(cls.id_format(fn))

    def test_load_off(self):
        """ Tests automatic loading of OFF files """
        off = formats.load_file(get_fn('amino12.lib'))
        self.assertIsInstance(off, dict)
        for key, item in iteritems(off):
            self.assertIsInstance(item, ResidueTemplate)

    def test_load_amber_prmtop(self):
        """ Tests automatic loading of AmberParm object """
        parm = formats.load_file(get_fn('trx.prmtop'))
        self.assertIsInstance(parm, amber.AmberParm)

    def test_load_amoeba_prmtop(self):
        """ Tests automatic loading of AmoebaParm object """
        parm = formats.load_file(get_fn('amoeba.parm7'))
        self.assertIsInstance(parm, amber.AmoebaParm)

    def test_load_chamber_prmtop(self):
        """ Tests automatic loading of ChamberParm object """
        parm = formats.load_file(get_fn('ala_ala_ala.parm7'))
        self.assertIsInstance(parm, amber.ChamberParm)

    def test_load_raw_amber_format(self):
        """ Tests automatic loading of RISM mdl (AmberFormat) object """
        parm = formats.load_file(get_fn('cSPCE.mdl'))
        self.assertIsInstance(parm, amber.AmberFormat)
        self.assertNotIsInstance(parm, amber.AmberParm)

    def test_load_amber_restart_ascii(self):
        """ Tests automatic loading of Amber ASCII restart file """
        parm = formats.load_file(get_fn('trx.inpcrd'))
        self.assertIsInstance(parm, amber.AmberAsciiRestart)

    def test_load_amber_restart_ascii_as_structure(self):
        """ Tests automatic loading of Amber ASCII restart file to Structure """
        parm = pmd.load_file(get_fn('ala3_solv.rst7'), structure=True)
        inpcrd = pmd.load_file(get_fn('ala3_solv.rst7'))
        self.assertIsInstance(parm, Structure)
        np.testing.assert_almost_equal(parm.box, inpcrd.box)
        np.testing.assert_almost_equal(parm.coordinates, inpcrd.coordinates[0])
        # dummy testing to assign box
        # issue #778
        parm.box = [0.]*6

    def test_load_amber_traj_ascii(self):
        """ Tests automatic loading of Amber mdcrd file """
        crd = formats.load_file(get_fn('tz2.truncoct.crd'), natom=5827,
                                hasbox=True)
        self.assertIsInstance(crd, amber.AmberMdcrd)
        self.assertRaises(TypeError, lambda:
                formats.load_file(get_fn('tz2.truncoct.crd')))

    def test_load_charmm_psf(self):
        """ Tests automatic loading of CHARMM PSF file """
        parm = formats.load_file(get_fn('ala_ala_ala.psf'))
        self.assertIsInstance(parm, charmm.CharmmPsfFile)

    def test_load_charmm_crd(self):
        """ Tests automatic loading of CHARMM crd file """
        crd = formats.load_file(get_fn('dhfr_min_charmm.crd'))
        self.assertIsInstance(crd, charmm.CharmmCrdFile)

    def test_load_charmm_rst(self):
        """ Tests automatic loading of CHARMM restart file """
        crd = formats.load_file(get_fn('sample-charmm.rst'))
        self.assertIsInstance(crd, charmm.CharmmRstFile)

    @unittest.skipIf(PYPY, 'Test does not yet run under pypy')
    def test_load_amber_restart_netcdf(self):
        """ Tests automatic loading of Amber NetCDF restart file """
        crd = formats.load_file(get_fn('ncinpcrd.rst7'))
        self.assertIsInstance(crd, amber.NetCDFRestart)
        self.assertFalse(amber.NetCDFTraj.id_format(get_fn('ncinpcrd.rst7')))
        self.assertFalse(amber.NetCDFTraj.id_format(get_fn('WMI_Lear.nc')))

    @unittest.skipIf(PYPY, 'Test does not yet run under pypy')
    def test_load_traj_netcdf(self):
        """ Tests automatic loading of Amber NetCDF trajectory file """
        crd = formats.load_file(get_fn('tz2.truncoct.nc'))
        self.assertIsInstance(crd, amber.NetCDFTraj)
        self.assertFalse(amber.NetCDFRestart.id_format(get_fn('tz2.truncoct.nc')))
        self.assertFalse(amber.NetCDFRestart.id_format(get_fn('WMI_Lear.nc')))

    def test_load_pdb(self):
        """ Tests automatic loading of PDB files """
        pdb = formats.load_file(get_fn('4lzt.pdb'))
        self.assertIsInstance(pdb, Structure)
        self.assertEqual(len(pdb.atoms), 1164)

    def test_load_reduced_pdb(self):
        """ Tests automatic loading of PDB files generated by reduce """
        self.assertTrue(formats.PDBFile.id_format(get_fn('reduce.pdb')))
        pdb = formats.load_file(get_fn('reduce.pdb'))
        self.assertIsInstance(pdb, Structure)
        self.assertEqual(len(pdb.atoms), 49)
        self.assertEqual(len(pdb.residues), 1)
        self.assertEqual(pdb.residues[0].name, 'SAM')

    def test_load_pdb_with_negative_resnum(self):
        """ Tests negative residue numbers in PDB writing """
        # Make a random structure
        struct = read_PDB(get_fn('4lzt.pdb'))
        for i, residue in enumerate(struct.residues):
            residue.number = i - 2
        for i, atom in enumerate(struct.atoms):
            atom.number = i - 2
        mypdb = self.get_fn('negative_indexes.pdb', written=True)
        struct.save(mypdb, renumber=False)
        struct2 = read_PDB(mypdb)
        self.assertEqual(len(struct.atoms), len(struct2.atoms))
        self.assertEqual(len(struct.residues), len(struct2.residues))
        # Now make sure the numbers are still negative
        for i, atom in enumerate(struct2.atoms):
            self.assertEqual(atom.number, i-2)
        for i, residue in enumerate(struct2.residues):
            self.assertEqual(residue.number, i-2)

    def test_load_cif(self):
        """ Tests automatic loading of PDBx/mmCIF files """
        cif = formats.load_file(get_fn('4LZT.cif'))
        self.assertIsInstance(cif, Structure)
        self.assertEqual(len(cif.atoms), 1164)

    def test_load_mol2(self):
        """ Tests automatic loading of mol2 and mol3 files """
        mol2 = formats.load_file(get_fn('test_multi.mol2'))
        self.assertIsInstance(mol2, ResidueTemplateContainer)
        mol3 = formats.load_file(get_fn('tripos9.mol2'))
        self.assertIsInstance(mol3, ResidueTemplate)
        # Check mol2 file where last 4 or 5 columns do not exist in ATOM
        mol2 = formats.load_file(get_fn('tripos2.mol2'))
        self.assertIsInstance(mol2, ResidueTemplate)
        for atom in mol2.atoms:
            self.assertEqual(atom.charge, 0)
            self.assertEqual(atom.residue.name, 'UNK')
        # Check mol2 where last several columns do not exist in SUBSTRUCTURE
        mol2 = formats.load_file(get_fn('tripos4.mol2'), structure=True)
        self.assertIsInstance(mol2, Structure)
        # Check bad file detection
        fn = self.get_fn('junk_file', written=True)
        with open(fn, 'w') as f:
            f.write('\n')
        self.assertFalse(formats.mol2.Mol2File.id_format(fn))
        with open(fn, 'w') as f:
            f.write('junkity junk junk\n')
            f.write('not a mol2 file\n')
        self.assertRaises(exceptions.Mol2Error, lambda:
                formats.mol2.Mol2File.parse(fn)
        )
        self.assertRaises(exceptions.Mol2Error, lambda:
                formats.mol2.Mol2File.parse(get_fn('error.mol2'))
        )

    @unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without GROMACS")
    def test_load_gromacs_topology(self):
        """ Tests automatic loading of Gromacs topology file """
        top = formats.load_file(get_fn('1aki.charmm27.top'))
        self.assertIsInstance(top, gromacs.GromacsTopologyFile)

    def test_load_gro(self):
        """ Tests automatic loading of Gromacs GRO file """
        gro = formats.load_file(get_fn('1aki.ff99sbildn.gro'))
        self.assertIsInstance(gro, Structure)

    def test_load_pqr(self):
        """ Tests automatic loading of PQR files """
        pqr = formats.load_file(get_fn('adk_open.pqr'))
        self.assertIsInstance(pqr, Structure)
        self.assertEqual(len(pqr.atoms), 3341)
        self.assertEqual(len(pqr.residues), 214)
        self.assertAlmostEqual(sum(a.charge for a in pqr.atoms), -4, places=4)
        self.assertEqual(pqr.atoms[0].charge, -0.30)
        self.assertEqual(pqr.atoms[0].solvent_radius, 1.85)
        self.assertEqual(pqr.atoms[0].atomic_number, 7)
        self.assertEqual(pqr.atoms[35].charge, -0.8)
        self.assertEqual(pqr.atoms[-1].charge, -0.67)
        self.assertEqual(pqr.atoms[-1].solvent_radius, 1.7)
        self.assertEqual(pqr.atoms[-1].atomic_number, 8)
        self.assertIsInstance(random.choice(pqr.residues).number, int)
        self.assertIsInstance(random.choice(pqr.atoms).number, int)

    def test_misdetect_pqr(self):
        """ Check that PQR autodetection does not identify a PDB file """
        pdb = read_PDB(get_fn('3p4a.pdb'))
        fname = self.get_fn('3p4a_chainA.pdb', written=True)
        pdb['A',:,:].save(fname)
        self.assertFalse(formats.PQRFile.id_format(fname))

    def test_bad_loads(self):
        """ Test exception handling when non-recognized files are loaded """
        self.assertRaises(exceptions.FormatNotFound, lambda:
                formats.load_file(get_fn('../test_parmed_formats.py')))
        self.assertRaises(IOError, lambda: formats.load_file('no_file'))

    def test_structure_keyword(self):
        """ Tests that the structure argument is special-cased in load_file """
        mol2 = formats.load_file(get_fn('tripos9.mol2'), structure=True)
        self.assertIsInstance(mol2, Structure)
        pdb = formats.load_file(get_fn('4lzt.pdb'), structure=True)

    def test_negative_residue_number(self):
        """ Tests automatic detection of PDB file with negative residue #s """
        pdb = formats.load_file(get_fn('1kx5.pdb'))
        self.assertTrue(any(res.number < 0 for res in pdb.residues))

    def test_dbref_keyword(self):
        """ Tests automatic detection of PDB file with DBREF record(s) """
        self.assertTrue(formats.PDBFile.id_format(get_fn('3p49.pdb')))

    def test_natom_hasbox_keywords(self):
        """ Tests that the hasbox/natom arguments are special-cased in load_file """
        crd = formats.load_file(get_fn('tz2.truncoct.crd'), natom=5827,
                                hasbox=True)
        self.assertIsInstance(crd, amber.AmberMdcrd)
        # Does not currently run under pypy
        if not PYPY:
            crd = formats.load_file(get_fn('tz2.truncoct.nc'), natom=5827,
                                    hasbox=True)
            self.assertIsInstance(crd, amber.NetCDFTraj)
        crd = formats.load_file(get_fn('trx.prmtop'), natom=5827,
                                hasbox=True)
        self.assertIsInstance(crd, amber.AmberParm)

    @unittest.skipUnless(has_rdkit and is_linux, "Only test load_rdkit module on Linux")
    def test_load_sdf(self):
        """ test load sdf format via rdkit """
        sdffile = get_fn('test.sdf')
        # structure = False
        parmlist = pmd.load_file(sdffile)
        self.assertIsInstance(parmlist, list)
        self.assertEqual(len(parmlist[0].atoms), 34)
        self.assertEqual(len(parmlist[1].atoms), 43)
        np.testing.assert_almost_equal(parmlist[0].coordinates[0], [2.0000, 2.7672, 0.0000], decimal=3)
        np.testing.assert_almost_equal(parmlist[0].coordinates[-1], [9.9858, -2.8473, 0.0000], decimal=3)
        np.testing.assert_almost_equal(parmlist[1].coordinates[0], [7.0468, -1.7307, 0.0000], decimal=3)
        np.testing.assert_almost_equal(parmlist[1].coordinates[-1], [1.5269, 2.1331, 0.0000], decimal=3)
        # structure = True
        parm = pmd.load_file(sdffile, structure=True)
        self.assertIsInstance(parm, Structure)
        self.assertEqual(len(parm.atoms), 34)
        np.testing.assert_almost_equal(parm.coordinates[0], [2.0000, 2.7672, 0.0000], decimal=3)
        np.testing.assert_almost_equal(parm.coordinates[-1], [9.9858, -2.8473, 0.0000], decimal=3)

class TestPDBStructure(FileIOTestCase):

    def setUp(self):
        self.pdb = get_fn('4lzt.pdb')
        self.pdbgz = get_fn('4lzt.pdb.gz')
        self.pdbbz2 = get_fn('4lzt.pdb.bz2')
        self.models = get_fn('2koc.pdb')
        self.overflow = get_fn('4lyt_vmd.pdb')
        self.simple = get_fn('ala_ala_ala.pdb')
        self.format_test = get_fn('SCM_A.pdb')
        self.overflow2 = get_fn('overflow.pdb')
        self.ATOMLINE = "ATOM  %5s %4s%1s%3s %1s%4s%-2s  %8s%8s%8s%6s%6s          %-2s%2s\n"
        self.ANISOULINE = "ANISOU%5s %-4s%1s%-4s%1s%4s%-2s%7s%7s%7s%7s%7s%7s      %2s%-2s\n"
        super().setUp()

    def test_pdb_anisou_inscode(self):
        """ Tests that PDB files with ANISOU records on inscodes work """
        formats.PDBFile.parse(get_fn('1gdu.pdb'))

    def test_pdb_format_detection(self):
        """ Tests PDB file detection from contents """
        fn = self.get_fn('test.pdb', written=True)
        pdbtext1 = "%-5s%d    %10.6f%10.6f%10.6f     %10.5f\n" + self.ATOMLINE
        with open(fn, 'w') as f:
            f.write(pdbtext1 % ('ORIGX', 1, 10, 10, 10, 10, 1, 'CA', '', 'ALA',
                'A', 1, '', '   1.000', '   1.000', '   1.000', '  1.00',
                '  1.00', '', ''))
        self.assertTrue(formats.PDBFile.id_format(fn))
        # ORIGX4 is not a valid keyword
        with open(fn, 'w') as f:
            f.write(pdbtext1 % ('ORIGX', 4, 10, 10, 10, 10, 1, 'CA', '', 'ALA',
                'A', 1, '', '   1.000', '   1.000', '   1.000', '  1.00',
                '  1.00', '', ''))
        self.assertFalse(formats.PDBFile.id_format(fn))
        # Safely catch if coordinates are not floats
        with open(fn, 'w') as f:
            f.write(pdbtext1 % ('ORIGX', 1, 10, 10, 10, 10, 1, 'CA', '', 'ALA',
                'A', 1, '', '   a.000', '   1.000', '   1.000', '  1.00',
                '  1.00', '', ''))
        self.assertFalse(formats.PDBFile.id_format(fn))
        # Safely catch if occupancy and b-factor are not floats
        with open(fn, 'w') as f:
            f.write(pdbtext1 % ('ORIGX', 1, 10, 10, 10, 10, 1, 'CA', '', 'ALA',
                'A', 1, '', '   1.000', '   1.000', '   1.000', '  a.00',
                '  1.00', '', ''))
        self.assertFalse(formats.PDBFile.id_format(fn))
        with open(fn, 'w') as f:
            f.write(pdbtext1 % ('ORIGX', 1, 10, 10, 10, 10, 1, 'CA', '', 'ALA',
                'A', 1, '', '   1.000', '   1.000', '   1.000', '  1.00',
                '  a.00', '', ''))
        self.assertFalse(formats.PDBFile.id_format(fn))
        # Safely catch if element is wrong
        with open(fn, 'w') as f:
            f.write(pdbtext1 % ('ORIGX', 1, 10, 10, 10, 10, 1, 'CA', '', 'ALA',
                'A', 1, '', '   1.000', '   1.000', '   1.000', '  1.00',
                '  1.00', 'C1', ''))
        self.assertFalse(formats.PDBFile.id_format(fn))
        # Safely catch blank file
        with open(fn, 'w') as f:
            pass
        self.assertFalse(formats.PDBFile.id_format(fn))

    def test_ascii(self):
        """ Test PDB file parsing """
        self._check4lzt(read_PDB(self.pdb))
        # The PDB file with multiple models
        pdbfile = read_PDB(open(self.models))
        all_crds = pdbfile.get_coordinates('all')
        self.assertEqual(all_crds.shape[0], 20)
        np.testing.assert_allclose(all_crds[0][0], [-8.886, -5.163, 9.647])
        np.testing.assert_allclose(all_crds[19][-1], [-12.051, 5.205, -2.146])

    @unittest.skipUnless(is_jenkins(), 'PDB blocks Travis from downloading files')
    def test_download(self):
        """ Tests downloading PDB files """
        self._check4lzt(download_PDB('4lzt'))
        # Check proper argument handling
        self.assertRaises(ValueError, lambda: download_PDB('not a PDB ID'))
        self.assertRaises(IOError, lambda: download_PDB('@#63'))

    @unittest.skipUnless(is_jenkins(), 'PDB blocks Travis from downloading files')
    def test_download_save(self):
        """ Tests downloading PDB files and saving a copy """
        fname = self.get_fn('downloaded.pdb', written=True)
        self._check4lzt(download_PDB('4lzt', saveto=fname))
        self._check4lzt(read_PDB(fname))

    def test_positions(self):
        """ Tests that positions are Vec3's with units """
        from parmed import unit as u
        from parmed import Vec3
        pdbfile = read_PDB(open(self.models))
        self.assertIsInstance(pdbfile.positions[0], u.Quantity)
        self.assertIsInstance(pdbfile.positions[0].value_in_unit(u.angstroms), Vec3)

    def test_gzip(self):
        """ Test Gzipped-PDB file parsing """
        self._check4lzt(read_PDB(self.pdbgz))

    def test_bzip(self):
        """ Test Bzipped-PDB file parsing """
        self._check4lzt(read_PDB(self.pdbbz2))

    @unittest.skipUnless(run_all_tests, 'Skipping large tests')
    def test_vmd_overflow(self):
        """ Test PDB file where atom and residue numbers overflow """
        pdbfile = read_PDB(self.overflow)
        self.assertEqual(len(pdbfile.atoms), 110237)
        self.assertEqual(len(pdbfile.residues), 35697)
        np.testing.assert_allclose(pdbfile.box, [0, 0, 0, 90, 90, 90])

    @unittest.skipUnless(run_all_tests, 'Skipping large tests')
    def test_regular_overflow(self):
        """ Test PDB file where atom number goes to ***** after 99999 """
        pdbfile = read_PDB(self.overflow2)
        self.assertEqual(len(pdbfile.atoms), 114277)
        self.assertEqual(len(pdbfile.residues), 25044)
        for i, atom in enumerate(pdbfile.atoms):
            self.assertEqual(atom.number, i+1)
            self.assertEqual(atom.idx, i)

    def test_residue_overflow(self):
        """ Tests PDB file where residue number overflows """
        fn = self.get_fn('test.pdb', written=True)
        pdbtext = self.ATOMLINE * 7
        with open(fn, 'w') as f:
            f.write('CRYST1%9.3f%9.3f%9.3f\n' % (10, 10, 10))
            f.write(pdbtext %
                (1, 'CA', ' ', 'RE1', 'A', 9999, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CA', ' ', 'RE2', 'A', hex(10000)[2:], '', 1, 1, 1, 1, 1, '', '',
                 3, 'CA', ' ', 'RE3', 'A', hex(10001)[2:], '', 1, 1, 1, 1, 1, '', '',
                 4, 'CA', ' ', 'RE4', 'A', hex(10002)[2:], '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE5', 'A', hex(10003)[2:], '', 1, 1, 1, 1, 1, '', '',
                 6, 'CA', ' ', 'RE6', 'A', 'ffff', '', 1, 1, 1, 1, 1, '', '',
                 7, 'CA', ' ', 'RE7', 'A', '****', '', 1, 1, 1, 1, 1, '', '')
            )
        # Check the parsing
        pdb = formats.PDBFile.parse(fn)
        np.testing.assert_equal(pdb.box, [10, 10, 10, 90, 90, 90])
        self.assertEqual(len(pdb.residues), 7)
        self.assertEqual(len(pdb.atoms), 7)
        self.assertEqual(pdb.residues[0].number, 9999)
        # Parser no longer tries to retain gaps in overflow PDB files. Make sure atoms are
        # numbered seqeuentially
        self.assertEqual(pdb.residues[1].number, 10000)
        self.assertEqual(pdb.residues[2].number, 10001)
        self.assertEqual(pdb.residues[3].number, 10002)
        self.assertEqual(pdb.residues[4].number, 10003)
        self.assertEqual(pdb.residues[5].number, 10004)
        self.assertEqual(pdb.residues[6].number, 10005)
        # Non-numerical residue numbers are only checked for uniqueness to compare with previous
        # residue to see if a new one should be created
        with open(fn, 'w') as f:
            f.write(pdbtext %
                (1, 'CA', ' ', 'RE1', 'A', 9999, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CA', ' ', 'RE2', 'A', hex(10000)[2:], '', 1, 1, 1, 1, 1, '', '',
                 3, 'CA', ' ', 'RE3', 'A', hex(10001)[2:], '', 1, 1, 1, 1, 1, '', '',
                 4, 'CA', ' ', 'RE4', 'A', hex(10002)[2:], '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE5', 'A', hex(10003)[2:], '', 1, 1, 1, 1, 1, '', '',
                 6, 'CA', ' ', 'RE6', 'A', 'ffff', '', 1, 1, 1, 1, 1, '', '',
                 7, 'CA', ' ', 'RE7', 'A', '>:-O', '', 1, 1, 1, 1, 1, '', '')
            )
        pdb = formats.PDBFile.parse(fn)
        self.assertEqual(len(pdb.residues), 7)
        self.assertEqual(len(pdb.atoms), 7)
        # Now check if the residue sequence field is simply expanded
        with open(fn, 'w') as f:
            f.write(pdbtext %
                (1, 'CA', ' ', 'RE1', 'A', 9999, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CA', ' ', 'RE2', 'A', 1000, '0', 1, 1, 1, 1, 1, '', '',
                 3, 'CA', ' ', 'RE3', 'A', 1000, '1', 1, 1, 1, 1, 1, '', '',
                 4, 'CA', ' ', 'RE4', 'A', 1000, '2', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE5', 'A', 9999, '9', 1, 1, 1, 1, 1, '', '',
                 6, 'CA', ' ', 'RE6', 'A', 1000, '00', 1, 1, 1, 1, 1, '', '',
                 7, 'CA', ' ', 'RE7', 'A', 1000, '01', 1, 1, 1, 1, 1, '', '')
            )
        # Check the parsing
        pdb = formats.PDBFile.parse(fn)
        self.assertEqual(len(pdb.residues), 7)
        self.assertEqual(len(pdb.atoms), 7)
        self.assertEqual(pdb.residues[0].number, 9999)
        self.assertEqual(pdb.residues[1].number, 10000)
        self.assertEqual(pdb.residues[2].number, 10001)
        self.assertEqual(pdb.residues[3].number, 10002)
        self.assertEqual(pdb.residues[4].number, 99999)
        self.assertEqual(pdb.residues[5].number, 10000)
        self.assertEqual(pdb.residues[6].number, 10000)

        # Check proper residue distinguishing for overflowed case. NR marks
        # atoms that *should* be turned into a new residue because there's a
        # name repeat or a new residue name
        with open(fn, 'w') as f:
            f.write(pdbtext %
                (1, 'CA', ' ', 'RE1', 'A', 9999, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CA', ' ', 'RE2', 'A', hex(10000)[2:], '', 1, 1, 1, 1, 1, '', '',
                 3, 'CA', ' ', 'RE3', 'A', 'ffff', '', 1, 1, 1, 1, 1, '', '',
                 4, 'CB', ' ', 'RE3', 'A', '****', '', 1, 1, 1, 1, 1, '', '',
                 5, 'CB', ' ', 'RE3', 'A', '****', '', 1, 1, 1, 1, 1, '', '', # NR
                 6, 'XX', ' ', 'RE4', 'A', '****', '', 1, 1, 1, 1, 1, '', '', # NR
                 7, 'EP', ' ', 'RE4', 'A', '****', '', 1, 1, 1, 1, 1, '', '')
            )
        # In the above, the PDB parser should recognize at the very least that there are at least
        # 5 distinct residues. Atoms 4 and 5 have the same residue name *and* atom name (and the
        # same residue number) so it's reasonable to treat them as the same or different residues
        pdb = formats.PDBFile.parse(fn)
        self.assertGreaterEqual(len(pdb.residues), 5)
        self.assertEqual([len(r) for r in pdb.residues[:3]], [1, 1, 1])
        self.assertIsInstance(pdb.atoms[6], topologyobjects.ExtraPoint)

    def test_atom_number_overflow(self):
        """ Tests PDB file where residue number overflows """
        fn = self.get_fn('test.pdb', written=True)
        pdbtext = self.ATOMLINE * 7
        with open(fn, 'w') as f:
            f.write(pdbtext %
                (99999, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 hex(100000)[2:], 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 hex(100001)[2:], 'CA', ' ', 'RE3', 'A', 3, '', 1, 1, 1, 1, 1, '', '',
                 hex(100002)[2:], 'CA', ' ', 'RE4', 'A', 4, '', 1, 1, 1, 1, 1, '', '',
                 hex(100003)[2:], 'CA', ' ', 'RE5', 'A', 5, '', 1, 1, 1, 1, 1, '', '',
                 'fffff', 'CA', ' ', 'RE6', 'A', 6, '', 1, 1, 1, 1, 1, '', '',
                 '*****', 'CA', ' ', 'RE7', 'A', 7, '', 1, 1, 1, 1, 1, '', '')
            )
        # Check the parsing
        pdb = formats.PDBFile.parse(fn)
        self.assertEqual(len(pdb.residues), 7)
        self.assertEqual(len(pdb.atoms), 7)
        self.assertEqual(pdb[0].number, 99999)
        # Parser no longer tries to retain gaps in overflow PDB files. Make sure atoms are
        # numbered seqeuentially
        self.assertEqual(pdb[1].number, 100000)
        self.assertEqual(pdb[2].number, 100001)
        self.assertEqual(pdb[3].number, 100002)
        self.assertEqual(pdb[4].number, 100003)
        self.assertEqual(pdb[5].number, 100004)
        self.assertEqual(pdb[6].number, 100005)
        with open(fn, 'w') as f:
            f.write(pdbtext %
                (99999, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 hex(100000)[2:], 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 hex(100001)[2:], 'CA', ' ', 'RE3', 'A', 3, '', 1, 1, 1, 1, 1, '', '',
                 hex(100002)[2:], 'CA', ' ', 'RE4', 'A', 4, '', 1, 1, 1, 1, 1, '', '',
                 hex(100003)[2:], 'CA', ' ', 'RE5', 'A', 5, '', 1, 1, 1, 1, 1, '', '',
                 'fffff', 'CA', ' ', 'RE6', 'A', 6, '', 1, 1, 1, 1, 1, '', '',
                 '*a***', 'CA', ' ', 'RE7', 'A', 7, '', 1, 1, 1, 1, 1, '', '')
            )
        # Check the parsing. PDB file no longer cares if atom number is not numeric
        pdb = formats.PDBFile.parse(fn)
        self.assertEqual(len(pdb.residues), 7)
        self.assertEqual(len(pdb.atoms), 7)
        self.assertEqual(pdb.atoms[-1].name, 'CA')
        self.assertEqual(pdb.atoms[-1].residue.name, 'RE7')

    def test_pdb_with_models(self):
        """ Test parsing of PDB files with multiple models """
        fn = self.get_fn('test.pdb', written=True)
        # Test working version
        with open(fn, 'w') as f:
            f.write("MODEL        1\n")
            f.write(self.ATOMLINE*7 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'CC', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write('ENDMDL\n')
            f.write("MODEL        2\n")
            f.write(self.ATOMLINE*7 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'CC', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write('ENDMDL\n')
        pdb = formats.PDBFile.parse(fn)
        self.assertEqual(pdb.get_coordinates().shape[0], 2)

        # Make sure it still works WITHOUT ENDMDL (to be permissive)
        with open(fn, 'w') as f:
            f.write("MODEL        1\n")
            f.write(self.ATOMLINE*7 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'CC', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write("MODEL        2\n")
            f.write(self.ATOMLINE*7 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'CC', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write('ENDMDL\n')
        with self.assertWarns(exceptions.PDBWarning):
            formats.PDBFile.parse(fn)
        pdb = formats.PDBFile.parse(fn)
        self.assertEqual(pdb.get_coordinates().shape[0], 2)

        with open(fn, 'w') as f:
            f.write("MODEL        1\n")
            f.write(self.ATOMLINE*7 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'CC', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write('ENDMDL\n')
            f.write("MODEL        2\n")
            f.write(self.ATOMLINE*8 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'CC', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 8, 'CD', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write('ENDMDL\n')
        self.assertRaises(exceptions.PDBError, lambda: formats.PDBFile.parse(fn))

        with open(fn, 'w') as f:
            f.write("MODEL        1\n")
            f.write(self.ATOMLINE*8 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'CC', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 8, 'CD', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write('ENDMDL\n')
            f.write("MODEL        2\n")
            f.write(self.ATOMLINE*7 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'CC', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write('ENDMDL\n')
        self.assertRaises(exceptions.PDBError, lambda: formats.PDBFile.parse(fn))

        with open(fn, 'w') as f:
            f.write("MODEL        1\n")
            f.write(self.ATOMLINE*7 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'CC', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write('ENDMDL\n')
            f.write("MODEL        2\n")
            f.write(self.ATOMLINE*7 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'XX', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write('ENDMDL\n')
        self.assertRaises(exceptions.PDBError, lambda: formats.PDBFile.parse(fn))

        # Make sure it still error checking works without ENDMDL
        with open(fn, 'w') as f:
            f.write("MODEL        1\n")
            f.write(self.ATOMLINE*7 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'CC', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write("MODEL        2\n")
            f.write(self.ATOMLINE*6 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write("MODEL        3\n")
        self.assertRaises(exceptions.PDBError, lambda: formats.PDBFile.parse(fn))

        # Make sure it still error checking works without ENDMDL on last MODEL
        with open(fn, 'w') as f:
            f.write("MODEL        1\n")
            f.write(self.ATOMLINE*7 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'CC', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write("MODEL        2\n")
            f.write(self.ATOMLINE*7 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 7, 'CC', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
            f.write("MODEL        3\n")
            f.write(self.ATOMLINE*6 %
                (1, 'CA', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 2, 'CB', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 3, 'CC', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 4, 'CD', ' ', 'RE1', 'A', 1, '', 1, 1, 1, 1, 1, '', '',
                 5, 'CA', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '',
                 6, 'CB', ' ', 'RE2', 'A', 2, '', 1, 1, 1, 1, 1, '', '')
            )
        self.assertRaises(exceptions.PDBError, lambda: formats.PDBFile.parse(fn))

        with open(fn, 'w') as f:
            f.write('ENDMDL\n')
        self.assertRaises(exceptions.PDBError, lambda: formats.PDBFile.parse(fn))

    def test_anisou_error_handling(self):
        """ Tests error detection/handling for bad ANISOU records in PDBs """
        fn = self.get_fn('test.pdb', written=True)
        with open(fn, 'w') as f:
            f.write(self.ATOMLINE % (1, 'CA', '', 'ALA', 'A', 1, '', 1, 1, 1, 1, 1, '', ''))
            f.write(self.ANISOULINE % (1, 'CA', '', 'ALA', 'A',  1, '', 10000, 20000, 30000, 40000, 50000, 60000, '', ''))
        pdb = formats.PDBFile.parse(fn)
        self.assertEqual(len(pdb.atoms), 1)
        np.testing.assert_equal(pdb.atoms[0].anisou, [1, 2, 3, 4, 5, 6])
        # Now test error handling
        # Bad atom number
        with open(fn, 'w') as f:
            f.write(self.ATOMLINE % (1, 'CA', '', 'ALA', 'A', 1, '',
                                     1, 1, 1, 1, 1, '', ''))
            f.write(self.ANISOULINE % ('a', 'CA', '', 'ALA', 'A',  1, '',
                                       10000, 20000, 30000, 40000, 50000, 60000,
                                       '', ''))
        with self.assertWarns(exceptions.PDBWarning):
            formats.PDBFile.parse(fn)
        self.assertIs(formats.PDBFile.parse(fn).atoms[0].anisou, None)
        # Bad residue number
        with open(fn, 'w') as f:
            f.write(self.ATOMLINE % (1, 'CA', '', 'ALA', 'A', 1, '', 1, 1, 1, 1, 1, '', ''))
            f.write(
                self.ANISOULINE % (
                    1, 'CA', '', 'ALA', 'A',  'a', '', 10000, 20000, 30000, 40000, 50000, 60000, '', ''
                )
            )
        with self.assertWarns(exceptions.PDBWarning):
            formats.PDBFile.parse(fn)
        self.assertIs(formats.PDBFile.parse(fn).atoms[0].anisou, None)
        # Bad U11
        with open(fn, 'w') as f:
            f.write(self.ATOMLINE % (1, 'CA', '', 'ALA', 'A', 1, '', 1, 1, 1, 1, 1, '', ''))
            f.write(
                self.ANISOULINE % (1, 'CA', '', 'ALA', 'A',  1, '', 'a', 20000, 30000, 40000, 50000, 60000, '', '')
            )
        with self.assertWarns(exceptions.PDBWarning):
            formats.PDBFile.parse(fn)
        self.assertIs(formats.PDBFile.parse(fn).atoms[0].anisou, None)
        # Orphaned ANISOU
        with open(fn, 'w') as f:
            f.write(self.ANISOULINE % (1, 'CA', '', 'ALA', 'A',  1, '',
                                       10000, 20000, 30000, 40000, 50000, 60000,
                                       '', ''))
        # Non-matching ANISOU
        with open(fn, 'w') as f:
            f.write(self.ATOMLINE % (1, 'CA', '', 'ALA', 'A', 1, '', 1, 1, 1, 1, 1, '', ''))
            f.write(self.ANISOULINE % (1, 'CB', '', 'ALA', 'A',  1, '',
                                       10000, 20000, 30000, 40000, 50000, 60000, '', ''))
        with self.assertWarns(exceptions.PDBWarning):
            formats.PDBFile.parse(fn)
        self.assertIs(formats.PDBFile.parse(fn).atoms[0].anisou, None)

    def test_pdb_write_simple(self):
        """ Test PDB file writing on a very simple input structure """
        pdbfile = read_PDB(self.simple)
        self.assertEqual(len(pdbfile.atoms), 33)
        self.assertEqual(len(pdbfile.residues), 3)
        output = StringIO()
        pdbfile.write_pdb(output)
        output.seek(0)
        pdbfile2 = read_PDB(output)
        self.assertEqual(len(pdbfile2.atoms), 33)
        self.assertEqual(len(pdbfile2.residues), 3)
        self._compareInputOutputPDBs(pdbfile, pdbfile2)
        # Test that passing the coordinates attribute works
        output = StringIO()
        pdbfile.write_pdb(output, coordinates=pdbfile.get_coordinates('all'))
        output.seek(0)
        pdbfile2 = read_PDB(output)
        self.assertEqual(len(pdbfile2.atoms), 33)
        self.assertEqual(len(pdbfile2.residues), 3)
        self._compareInputOutputPDBs(pdbfile, pdbfile2)
        # Check some input parameter checking
        output = StringIO()
        self.assertRaises(ValueError, lambda:
                pdbfile.write_pdb(output, altlocs='illegal')
        )
        output = StringIO()
        self.assertRaises(TypeError, lambda:
                pdbfile.write_pdb(output, coordinates=[0, 1, 2])
        )

    def test_write_long_names(self):
        """ Tests writing long atom and residue names in PDB """
        struct = Structure()
        atom = Atom(name='CBDEF', atomic_number=7, altloc='A')
        oatom = Atom(name='CBDEF', atomic_number=7, altloc='B')
        atom.xx, atom.xy, atom.xz = 1, 1, 1
        oatom.xx, oatom.xy, oatom.xz = 2, 2, 2
        struct.add_atom(atom, 'RESIDUE', 1, 'A')
        struct.atoms[-1].other_locations['B'] = oatom
        fn = self.get_fn('test.pdb', written=True)
        struct.write_pdb(fn)
        pdb = formats.PDBFile.parse(fn)
        self.assertEqual(pdb.atoms[0].name, 'CBDE')
        self.assertEqual(pdb.residues[0].name, 'RES')
        self.assertEqual(len(pdb.atoms[0].other_locations), 1)
        self.assertEqual(pdb.atoms[0].other_locations['B'].name, 'CBDE')

    def test_pdb_write_models(self):
        """ Test PDB file writing from NMR structure with models """
        pdbfile = read_PDB(self.models)
        self.assertEqual(pdbfile.get_coordinates('all').shape, (20, 451, 3))
        self.assertEqual(len(pdbfile.atoms), 451)
        output = StringIO()
        pdbfile.write_pdb(output)
        output.seek(0)
        pdbfile2 = read_PDB(output)
        self.assertEqual(len(pdbfile2.atoms), 451)
        self.assertEqual(pdbfile2.get_coordinates('all').shape, (20, 451, 3))
        np.testing.assert_allclose(pdbfile2.get_coordinates('all'),
                                   pdbfile.get_coordinates('all'))
        self._compareInputOutputPDBs(pdbfile, pdbfile2)

    def test_ter_cards(self):
        """ Tests that the addition of TER cards is correct in PDB writing """
        pdbfile = read_PDB(get_fn('ala_ala_ala.pdb'))
        pdbfile *= 5 # This should make us need 5 TER cards
        fn = self.get_fn('test.pdb', written=True)
        pdbfile.write_pdb(fn)
        with open(fn, 'r') as f:
            self.assertEqual(sum([l.startswith('TER') for l in f]), 5)
        # Make sure TER cards *don't* get added for water
        pdbfile = read_PDB(get_fn('4lzt.pdb'))
        pdbfile.write_pdb(fn)
        with open(fn, 'r') as f:
            for line in f:
                if line.startswith('TER'):
                    break
            else:
                assert False, 'No TER card found!'
            # Make sure the rest of the atoms after this are HETATM
            has_hetatms = False
            for line in f:
                self.assertFalse(line.startswith('ATOM'))
                has_hetatms = has_hetatms or line.startswith('HETATM')
            self.assertTrue(has_hetatms)

    def test_ter_copy(self):
        """ Test that copying a Structure preserves TER card presence """
        pdbfile = read_PDB(get_fn('ala_ala_ala.pdb')) * 5
        fn = self.get_fn('test.pdb', written=True)
        pdbfile.write_pdb(fn)
        parsed = read_PDB(fn)
        self.assertEqual(sum([r.ter for r in parsed.residues]), 5)
        self.assertEqual(sum([r.ter for r in copy(parsed).residues]), 5)

    def test_ter_not_increase_tercount(self):
        s = """
ATOM      1  N   CYX L   1      57.464  29.769  15.871  1.00 25.39           N
ATOM      2  SG  CYX L   1      56.982  27.807  18.150  1.00 14.20           S
TER       3      CYX L   1
ATOM      4  N   CYX H   2      36.233  17.035  12.739  1.00 10.49           N
ATOM      5  SG  CYX H   2      36.833  15.443  15.640  1.00 15.60           S
"""
        parm = pmd.read_PDB(StringIO(s))
        buf = StringIO()
        parm.write_pdb(buf)
        buf.seek(0)
        content = buf.read()
        assert "TER       3      CYX L   1" in content

        buf = StringIO()
        parm.write_pdb(buf, increase_tercount=False)
        buf.seek(0)
        content = buf.read()
        assert "TER       2      CYX L   1" in content

    def test_pdb_big_coordinates(self):
        """ Test proper PDB coordinate parsing for large coordinates """
        pdbfile = read_PDB(get_fn('bigz.pdb'))
        self.assertAlmostEqual(pdbfile.coordinates[0,0], -100.024)
        self.assertAlmostEqual(pdbfile.coordinates[0,1], -100.103)
        self.assertAlmostEqual(pdbfile.coordinates[0,2], -100.101)

    def test_pdb_write_xtal(self):
        """ Test PDB file writing from a Xtal structure """
        pdbfile = read_PDB(self.pdb)
        self._check4lzt(pdbfile)
        output = StringIO()
        pdbfile.write_pdb(output, renumber=False)
        output.seek(0)
        pdbfile2 = read_PDB(output)
        self._check4lzt(pdbfile2, check_meta=False)
        self._compareInputOutputPDBs(pdbfile, pdbfile2)
        output = reset_stringio(output)
        pdbfile.write_pdb(output)
        output.seek(0)
        pdbfile3 = read_PDB(output)
        self._check4lzt(pdbfile3, check_meta=False)
        self._compareInputOutputPDBs(pdbfile, pdbfile3, True)
        # Now check that renumbering is done correctly. 4lzt skips residues 130
        # through 200
        for res1, res2 in zip(pdbfile.residues, pdbfile3.residues):
            if res1.idx < 129:
                self.assertEqual(res1.number, res2.number)
            elif res1.idx < 135:
                self.assertEqual(res1.number, res2.number + 71)
            else:
                # Some residue numbers are skipped in the water numbering
                self.assertGreaterEqual(res1.number, res2.number + 71 + 794)

    def test_pdb_write_altloc_options(self):
        """ Test PDB file writing with different altloc options """
        pdbfile = read_PDB(self.pdb)
        self._check4lzt(pdbfile)
        output = StringIO()
        pdbfile.write_pdb(output, renumber=False, altlocs='all')
        output.seek(0)
        pdbfile2 = read_PDB(output)
        self._check4lzt(pdbfile2, check_meta=False)
        self._compareInputOutputPDBs(pdbfile, pdbfile2)
        # Check that 'first' option works
        output = reset_stringio(output)
        pdbfile.write_pdb(output, renumber=False, altlocs='first')
        output.seek(0)
        pdbfile3 = read_PDB(output)
        self._check4lzt(pdbfile3, check_meta=False, has_altloc=False)
        self._compareInputOutputPDBs(pdbfile, pdbfile3, altloc_option='first')
        # Check that the 'occupancy' option works
        output = reset_stringio(output)
        pdbfile.write_pdb(output, renumber=False, altlocs='occupancy')
        output.seek(0)
        pdbfile4 = read_PDB(output)
        self._check4lzt(pdbfile4, check_meta=False, has_altloc=False)
        self._compareInputOutputPDBs(pdbfile, pdbfile4, altloc_option='occupancy')
        # Double-check 'first' vs. 'occupancy'. Residue 85 (SER) has a conformer
        # A that has an occupancy of 0.37 and conformer B with occupancy 0.63
        self.assertEqual(pdbfile3.residues[84][4].xx, -4.162)
        self.assertEqual(pdbfile4.residues[84][4].xx, -4.157)

    def test_pdb_write_standard_names(self):
        """ Test PDB file writing converting to standard names """
        parm = formats.load_file(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        output = StringIO()
        parm.write_pdb(output, standard_resnames=True)
        output.seek(0)
        pdb = read_PDB(output)
        for res in pdb.residues:
            self.assertEqual(
                    residue.AminoAcidResidue.get(res.name).abbr, res.name
            )

    def test_pdb_write_standard_names_water(self):
        """ Test water residue name translation in PDB writing """
        parm = formats.load_file(get_fn('nma.pdb'))
        resname_set = set(res.name for res in parm.residues)
        self.assertIn('WAT', resname_set)
        self.assertNotIn('HOH', resname_set)
        assert 'HOH' not in resname_set
        output = StringIO()
        parm.write_pdb(output, standard_resnames=True)
        output.seek(0)
        pdb = read_PDB(output)
        resname_set = set(res.name for res in pdb.residues)
        self.assertNotIn('WAT', resname_set)
        self.assertIn('HOH', resname_set)

    def test_pdb_write_hetatoms(self):
        """Tests HETATM/ATOM tag writing"""
        structure = Structure()
        a = Atom(name='CA', atomic_number=6)
        structure.add_atom(copy(a), 'ASH', 2, 'A')
        structure.add_atom(copy(a), 'DG', 2, 'A')
        structure.add_atom(copy(a), 'T', 2, 'A')
        structure.add_atom(copy(a), 'MOL', 2, 'A')

        coordinates = np.zeros((len(structure.atoms), 3))

        output = StringIO()

        tests = [{"use_hetatoms": True, "tags": ['ATOM', 'ATOM', 'ATOM', 'HETATM']},
                 {"use_hetatoms": False, "tags": ['ATOM', 'ATOM', 'ATOM', 'ATOM']}]

        for test in tests:
            output.seek(0)
            structure.write_pdb(output, use_hetatoms=test["use_hetatoms"], coordinates=coordinates)
            output.seek(0)

            for tag in test["tags"]:
                assert output.readline().startswith(tag)

    def test_anisou_read(self):
        """ Tests that read_PDB properly reads ANISOU records """
        pdbfile = read_PDB(self.pdb)
        aniso1 = [2066, 1204, 1269, 44, 126, 191] # first atom's ANISOU record
        aniso2 = [2090, 1182, 921, 46, 64, 60]    # second atom's ANISOU record
        aniso3 = [3057, 3932, 5304, 126, -937, -661] # last atom's ANISOU
        self.assertEqual(len(aniso1), len(pdbfile.atoms[0].anisou))
        for x, y in zip(aniso1, pdbfile.atoms[0].anisou):
            self.assertEqual(x/10000, y)
        self.assertEqual(len(aniso2), len(pdbfile.atoms[1].anisou))
        for x, y in zip(aniso2, pdbfile.atoms[1].anisou):
            self.assertEqual(x/10000, y)
        self.assertEqual(len(aniso3), len(pdbfile.atoms[-1].anisou))
        for x, y in zip(aniso3, pdbfile.atoms[-1].anisou):
            self.assertEqual(x/10000, y)

    def test_anisou_write(self):
        """ Tests that write_pdb properly writes ANISOU records """
        def check_aniso(pdbfile):
            aniso1 = [2066, 1204, 1269, 44, 126, 191]
            aniso2 = [2090, 1182, 921, 46, 64, 60]
            aniso3 = [3057, 3932, 5304, 126, -937, -661]
            self.assertEqual(len(aniso1), len(pdbfile.atoms[0].anisou))
            for x, y in zip(aniso1, pdbfile.atoms[0].anisou):
                self.assertEqual(x/10000, y)
            self.assertEqual(len(aniso2), len(pdbfile.atoms[1].anisou))
            for x, y in zip(aniso2, pdbfile.atoms[1].anisou):
                self.assertEqual(x/10000, y)
            self.assertEqual(len(aniso3), len(pdbfile.atoms[-1].anisou))
            for x, y in zip(aniso3, pdbfile.atoms[-1].anisou):
                self.assertEqual(x/10000, y)
        pdbfile = read_PDB(self.pdb)
        check_aniso(pdbfile)
        output = StringIO()
        pdbfile.write_pdb(output)
        output.seek(0)
        pdbfile2 = read_PDB(output)
        # Should have no anisou records, since by default they are not written
        for atom in pdbfile2.atoms:
            self.assertIs(atom.anisou, None)
        output = reset_stringio(output)
        pdbfile.write_pdb(output, renumber=False, write_anisou=True)
        output.seek(0)
        # This one should have anisou records
        pdbfile3 = read_PDB(output)
        self._compareInputOutputPDBs(pdbfile, pdbfile3)
        for a1, a2 in zip(pdbfile.atoms, pdbfile3.atoms):
            self.assertEqual(a1.anisou.shape, a2.anisou.shape)
            for x, y in zip(a1.anisou, a2.anisou):
                self.assertAlmostEqual(x, y, delta=1e-4)
            self.assertEqual(len(a1.other_locations), len(a2.other_locations))
            for key in sorted(a1.other_locations.keys()):
                oa1 = a1.other_locations[key]
                oa2 = a2.other_locations[key]
                self.assertEqual(oa1.anisou.shape, oa2.anisou.shape)
                for x, y in zip(oa1.anisou, oa2.anisou):
                    self.assertAlmostEqual(x, y, delta=1e-4)

    def test_pdb_write_format(self):
        """ Test PDB atom names are properly justified per PDB standard """
        pdbfile = read_PDB(self.format_test)
        f = self.get_fn('pdb_format_test.pdb', written=True)
        pdbfile.write_pdb(f, write_anisou=True)
        self.assertTrue(diff_files(self.get_fn('SCM_A_formatted.pdb', saved=True), f))

    def test_pdb_multimodel_parsing_bug_820(self):
        """ Test model failing in parsing due to bug #820 in GitHub """
        # Just make sure it does not raise an exception
        self.assertEqual(len(read_PDB(get_fn('1aaf.pdb')).atoms), 893)

    def test_pdb_write_symmetry_data(self):
        """ Tests writing PDB file with symmetry data """
        def assert_remark_290(parm, remark_290_lines):
            output = StringIO()
            parm.write_pdb(output)
            output.seek(0)
            buffer = output.read()
            for line in remark_290_lines.split():
                self.assertTrue(line in buffer)

        # 4lzt
        pdbfile = get_fn('4lzt.pdb')
        parm = pmd.load_file(pdbfile)
        remark_290_lines = """
REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000
REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000
REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000
"""
        assert_remark_290(parm, remark_290_lines)

        # 2idg
        parm = read_PDB(get_fn('2igd.pdb'))
        remark_290_lines = """
REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000
REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000
REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000
REMARK 290   SMTRY1   2 -1.000000  0.000000  0.000000       17.52500
REMARK 290   SMTRY2   2  0.000000 -1.000000  0.000000        0.00000
REMARK 290   SMTRY3   2  0.000000  0.000000  1.000000       21.18500
REMARK 290   SMTRY1   3 -1.000000  0.000000  0.000000        0.00000
REMARK 290   SMTRY2   3  0.000000  1.000000  0.000000       20.25000
REMARK 290   SMTRY3   3  0.000000  0.000000 -1.000000       21.18500
REMARK 290   SMTRY1   4  1.000000  0.000000  0.000000       17.52500
REMARK 290   SMTRY2   4  0.000000 -1.000000  0.000000       20.25000
REMARK 290   SMTRY3   4  0.000000  0.000000 -1.000000        0.00000
"""
        assert_remark_290(parm, remark_290_lines)

        self.assertRaises(ValueError, lambda: Symmetry(np.arange(100).reshape(10, 10)))

    def test_segid_handling(self):
        """ Test handling of CHARMM-specific SEGID identifier (r/w) """
        pdbfile = read_PDB(self.overflow2, skip_bonds=True) # Big file... skip bond check
        allsegids = set(['PROA', 'PROB', 'CARA', 'CARE', 'CARC', 'CARD', 'CARB',
                         'MEMB', 'TIP3', 'POT', 'CLA'])
        foundsegids = set()
        for residue in pdbfile.residues:
            foundsegids.add(residue.segid)
        self.assertEqual(foundsegids, allsegids)
        self.assertEqual(pdbfile.atoms[0].residue.segid, 'PROA')
        self.assertEqual(pdbfile.atoms[5161].residue.segid, 'PROA')
        self.assertEqual(pdbfile.atoms[5162].residue.segid, 'PROB')
        self.assertEqual(pdbfile.atoms[-1].residue.segid, 'CLA')
        f = self.get_fn('pdb_segid_test.pdb', written=True)
        pdbfile.write_pdb(f, charmm=True)
        pdbfile2 = read_PDB(f)
        for residue in pdbfile2.residues:
            self.assertTrue(residue.segid)
            self.assertEqual(residue.segid, pdbfile.residues[residue.idx].segid)

    def test_private_functions(self):
        """ Tests the private helper functions in parmed/formats/pdb.py """
        self.assertEqual(formats.pdb._standardize_resname('ASH'), ('ASP', False))
        self.assertEqual(formats.pdb._standardize_resname('CYX'), ('CYS', False))
        self.assertEqual(formats.pdb._standardize_resname('RA'), ('A', False))
        self.assertEqual(formats.pdb._standardize_resname('DG'), ('DG', False))
        self.assertEqual(formats.pdb._standardize_resname('BLA'), ('BLA', True))
        self.assertEqual(formats.pdb._standardize_resname('WAT'), ('HOH', True))
        self.assertEqual(formats.pdb._standardize_resname('TIP3'), ('HOH', True))
        # Make sure standard residues return themselves
        for res in residue.AminoAcidResidue.all_residues:
            self.assertEqual(formats.pdb._standardize_resname(res.abbr), (res.abbr, False))
        for res in residue.DNAResidue.all_residues:
            self.assertEqual(formats.pdb._standardize_resname(res.abbr), (res.abbr, False))
        for res in residue.RNAResidue.all_residues:
            self.assertEqual(formats.pdb._standardize_resname(res.abbr), (res.abbr, False))

    def test_deprecations(self):
        """ Test functions that raise deprecation warnings """
        fn = self.get_fn('blah', written=True)
        parm = formats.load_file(get_fn('ash.parm7'), get_fn('ash.rst7'))
        with self.assertWarns(DeprecationWarning):
            write_PDB(parm, fn)
        with self.assertWarns(DeprecationWarning):
            write_CIF(parm, fn)

    def test_link(self):
        """ Tests proper handling and processing of LINK records in PDB files """
        parm = pmd.load_file(get_fn('5qk8.pdb'))
        # This PDB file has 47 LINK records, all have symmetry operations 1555
        self.assertEqual(len(parm.links), 47) # 47 LINK records in this PDB file
        # Here are the first few lines of those LINK records:
        # LINK         O   ALA A  96                MG    MG A 301     1555   1555  2.09  
        # LINK         OE2 GLU A 112                MG    MG A 302     1555   1555  2.06  
        # Check these by hand
        link1, link2 = parm.links[:2]
        self.assertEqual(link1.atom1.name, 'O')
        self.assertEqual(link1.atom2.name, 'MG')
        self.assertEqual(link1.atom1.residue.name, 'ALA')
        self.assertEqual(link1.atom1.residue.number, 96)
        self.assertEqual(link1.atom2.residue.name, 'MG')
        self.assertEqual(link1.atom2.residue.number, 301)
        self.assertEqual(link1.symmetry_op1, '1555')
        self.assertEqual(link1.symmetry_op2, '1555')

        self.assertEqual(link2.atom1.name, 'OE2')
        self.assertEqual(link2.atom2.name, 'MG')
        self.assertEqual(link2.atom1.residue.name, 'GLU')
        self.assertEqual(link2.atom1.residue.number, 112)
        self.assertEqual(link2.atom2.residue.name, 'MG')
        self.assertEqual(link2.atom2.residue.number, 302)
        self.assertEqual(link2.symmetry_op1, '1555')
        self.assertEqual(link2.symmetry_op2, '1555')

        # Now test writing
        written_file = self.get_fn('link.pdb', written=True)
        parm.write_pdb(written_file, write_links=True, renumber=False)
        parm2 = pmd.load_file(written_file)
        self.assertEqual(len(parm2.links), 47)

    # Private helper test functions
    def _compareInputOutputPDBs(self, pdbfile, pdbfile2, reordered=False,
                                altloc_option='all'):
        # Now go through all atoms and compare their attributes
        for a1, a2 in zip(pdbfile.atoms, pdbfile2.atoms):
            if altloc_option in ('first', 'all'):
                self.assertEqual(a1.occupancy, a2.occupancy)
                a1idx = a1.idx
            elif altloc_option == 'occupancy':
                a, occ = a1, a1.occupancy
                for key, oa in iteritems(a1.other_locations):
                    if oa.occupancy > occ:
                        occ = oa.occupancy
                        a = oa
                a1idx = a1.idx
                a1 = a # This is the atom we want to compare with
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.bfactor, a2.bfactor)
            self.assertEqual(a1.altloc, a2.altloc)
            self.assertEqual(a1idx, a2.idx)
            if altloc_option == 'all':
                self.assertEqual(set(a1.other_locations.keys()),
                                 set(a2.other_locations.keys()))
            self.assertEqual(a1.xx, a2.xx)
            self.assertEqual(a1.xy, a2.xy)
            self.assertEqual(a1.xz, a2.xz)
            if altloc_option != 'all':
                # There should be no alternate locations unless we keep them all
                self.assertEqual(len(a2.other_locations), 0)
            if not reordered:
                self.assertEqual(a1.number, a2.number)
            # Search all alternate locations as well
            for k1, k2 in zip(sorted(a1.other_locations.keys()),
                              sorted(a2.other_locations.keys())):
                self.assertEqual(k1, k2)
                oa1 = a1.other_locations[k1]
                oa2 = a2.other_locations[k2]
                self.assertEqual(oa1.atomic_number, oa2.atomic_number)
                self.assertEqual(oa1.name, oa2.name)
                self.assertEqual(oa1.type, oa2.type)
                self.assertEqual(oa1.mass, oa2.mass)
                self.assertEqual(oa1.charge, oa2.charge)
                self.assertEqual(oa1.occupancy, oa2.occupancy)
                self.assertEqual(oa1.bfactor, oa2.bfactor)
                self.assertEqual(oa1.altloc, oa2.altloc)
                self.assertEqual(oa1.idx, oa2.idx)
                if not reordered:
                    self.assertEqual(oa1.number, oa2.number)
        # Now compare all residues
        for r1, r2 in zip(pdbfile.residues, pdbfile2.residues):
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.idx, r2.idx)
            self.assertEqual(r1.ter, r2.ter)
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.insertion_code, r2.insertion_code)
            if not reordered:
                self.assertEqual(r1.number, r2.number)

    def _check4lzt(self, obj, check_meta=True, has_altloc=True):
        self.assertEqual(obj.get_coordinates('all').shape[0], 1)
        np.testing.assert_allclose(obj.box,
                         [27.24, 31.87, 34.23, 88.52, 108.53, 111.89])
        self.assertEqual(obj.space_group, 'P 1')
        self.assertEqual(len(obj.atoms), 1164)
        self.assertEqual(len(obj.residues[0]), 9)
        # Check that alternate conformations are taken into account
        total_natoms = 0
        for i, atom in enumerate(obj.atoms):
            total_natoms += 1
            for key in atom.other_locations:
                total_natoms += 1
                atom2 = atom.other_locations[key]
                self.assertEqual(atom.altloc, 'A')
                self.assertEqual(atom2.altloc, 'B')
                if i in [388, 389]:
                    # Sum of atom 388/389 occupancies is 1.02
                    self.assertEqual(atom2.occupancy + atom.occupancy, 1.02)
                else:
                    # Other atoms occupancy sums are 1 exactly
                    self.assertEqual(atom2.occupancy + atom.occupancy, 1)
        if has_altloc:
            self.assertEqual(total_natoms, 1183)
            self.assertEqual(len(obj.atoms), 1164)
        else:
            self.assertEqual(total_natoms, 1164) # 19 atoms have altlocs
        # Check the metadata
        if check_meta:
            self.assertEqual(obj.experimental, 'X-RAY DIFFRACTION')
            self.assertEqual(len(obj.residues), 274)
            self.assertEqual(obj.pmid, '9761848')
            self.assertEqual(obj.journal_authors, 'M.A.WALSH,T.R.SCHNEIDER,'
                             'L.C.SIEKER,Z.DAUTER,V.S.LAMZIN,K.S.WILSON')
            self.assertEqual(obj.journal, 'ACTA CRYSTALLOGR.,SECT.D')
            self.assertEqual(obj.year, 1998)
            self.assertEqual(obj.keywords, ['HYDROLASE', 'O-GLYCOSYL',
                             'GLYCOSIDASE'])
            self.assertEqual(obj.title, 'REFINEMENT OF TRICLINIC HEN EGG-WHITE '
                             'LYSOZYME AT ATOMIC RESOLUTION.')
            self.assertEqual(obj.doi, '10.1107/S0907444997013656')
            self.assertEqual(obj.volume, '54')
            self.assertEqual(obj.page, '522')
            self.assertEqual(obj.resolution, 0.95)
        # Check the TER card is picked up
        for i, residue in enumerate(obj.residues):
            if i == 128:
                self.assertTrue(residue.ter)
            else:
                self.assertFalse(residue.ter)

class TestParmedPQRStructure(FileIOTestCase):
    """ Tests the PQR parser and writer """

    def setUp(self):
        self.ATOMLINE = "ATOM  %5s %4s %3s %1s %4s %7s %7s %7s %5s %5s\n"
        self.ATOMLINE2 = "ATOM  %5s %4s %3s %1s %4s %7s %7s %7s %5s %5s   %s %s\n"
        self.ATOMLINE3 = "ATOM  %5s %4s %3s %4s %7s %7s %7s %5s %5s\n"
        super().setUp()

    def test_pqr_parsing(self):
        """ Tests parsing a PQR file """
        fn = self.get_fn('test.pqr', written=True)
        self._check_adk_pqr(formats.PQRFile.parse(get_fn('adk_open.pqr')))
        with open(get_fn('adk_open.pqr'), 'r') as f:
            self._check_adk_pqr(formats.PQRFile.parse(f))
        # Test some variants of PQR files
        with open(fn, 'w') as f:
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000',
                                      '  1.000', '', '')
            )
        self.assertRaises(ValueError, lambda: formats.PQRFile.parse(fn))
        with open(fn, 'w') as f:
            f.write(self.ATOMLINE3 % (1, 'H1', 'HOH', 1, '  0.000', '  0.000',
                                      '  0.000', ' 0.500', ' 0.950'))
            f.write(self.ATOMLINE3 % (2, 'H2', 'HOH', 1, '  1.000', '  0.000',
                                      '  0.000', ' 0.500', ' 0.950'))
            f.write(self.ATOMLINE3 % (3, 'O', 'HOH', 1, '  0.500', '  1.000',
                                      '  0.000', ' 0.000', ' 1.500'))
            f.write(self.ATOMLINE3 % (4, 'EP', 'HOH', 1, '  0.500', '  0.500',
                                      '  0.000', '-1.000', ' 0.000'))
        pqr = formats.PQRFile.parse(fn)
        self.assertIsInstance(pqr.atoms[3], topologyobjects.ExtraPoint)

    def test_pqr_with_cryst1(self):
        """ Tests parsing PQR files with CRYST1 record """
        fn = self.get_fn('test.pqr', written=True)
        with open(fn, 'w') as f:
            f.write('CRYST1   10.0   10.0   10.0   109.47  109.47   109.47\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000',
                                      '  1.000', ' -0.500', ' 1.200'))
        pqr = formats.PQRFile.parse(fn)
        np.testing.assert_equal(pqr.box, [10, 10, 10, 109.47, 109.47, 109.47])
        with open(fn, 'w') as f:
            f.write('CRYST1   10.0   10.0   10.0\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000',
                                      '  1.000', ' -0.500', ' 1.200'))
        pqr = formats.PQRFile.parse(fn)
        np.testing.assert_equal(pqr.box, [10, 10, 10, 90, 90, 90])

    def test_pqr_parsing_with_models(self):
        """ Tests parsing PQR files with multiple models """
        fn = self.get_fn('test.pqr', written=True)
        with open(fn, 'w') as f:
            f.write('MODEL   1\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000', '  1.000', ' -0.500', ' 1.200'))
            f.write('ENDMDL\n')
            f.write('MODEL   2\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  2.000', '  2.000', '  2.000', ' -0.500', ' 1.200'))
            f.write('ENDMDL\n')
        pqr = formats.PQRFile.parse(fn)
        self.assertEqual(pqr.get_coordinates().shape, (2, 1, 3))
        with open(fn, 'w') as f:
            f.write('MODEL   1\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000', '  1.000', ' -0.500', ' 1.200'))
            f.write('MODEL   2\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  2.000', '  2.000', '  2.000', ' -0.500', ' 1.200'))
        with self.assertWarns(exceptions.PDBWarning):
            formats.PQRFile.parse(fn)
        pqr = formats.PQRFile.parse(fn)
        self.assertEqual(pqr.get_coordinates().shape, (2, 1, 3))
        # Check some error catching
        with open(fn, 'w') as f:
            f.write('MODEL   1\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000', '  1.000', ' -0.500', ' 1.200'))
            f.write('ENDMDL\n')
            f.write('MODEL   2\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  2.000', '  2.000', '  2.000', ' -0.500', ' 1.200'))
            f.write(self.ATOMLINE3 % (2, 'CB', 'ALA', 1, '  2.100', '  2.100', '  2.100', ' -0.500', ' 1.200'))
            f.write('ENDMDL\n')
        self.assertRaises(exceptions.PDBError, lambda: formats.PQRFile.parse(fn))
        with open(fn, 'w') as f:
            f.write('MODEL   1\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000', '  1.000', ' -0.500', ' 1.200'))
            f.write(self.ATOMLINE3 % (2, 'CB', 'ALA', 1, '  2.100', '  2.100', '  2.100', ' -0.500', ' 1.200'))
            f.write('ENDMDL\n')
            f.write('MODEL   2\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  2.000', '  2.000', '  2.000', ' -0.500', ' 1.200'))
            f.write('ENDMDL\n')
        self.assertRaises(exceptions.PDBError, lambda: formats.PQRFile.parse(fn))
        with open(fn, 'w') as f:
            f.write('MODEL   1\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000', '  1.000', ' -0.500', ' 1.200'))
            f.write('ENDMDL\n')
            f.write('MODEL   2\n')
            f.write(self.ATOMLINE3 % (1, 'CB', 'ALA', 1, '  2.000', '  2.000', '  2.000', ' -0.500', ' 1.200'))
            f.write('ENDMDL\n')
        self.assertRaises(exceptions.PDBError, lambda: formats.PQRFile.parse(fn))
        with open(fn, 'w') as f:
            f.write('ENDMDL\n')
        self.assertRaises(exceptions.PDBError, lambda: formats.PQRFile.parse(fn))
        with open(fn, 'w') as f:
            f.write('MODEL   1\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000', '  1.000', ' -0.500', ' 1.200'))
            f.write('MODEL   2\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000', '  1.000', ' -0.500', ' 1.200'))
            f.write(self.ATOMLINE3 % (2, 'CB', 'ALA', 1, '  2.100', '  2.100', '  2.100', ' -0.500', ' 1.200'))
            f.write('MODEL   3\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  2.000', '  2.000', '  2.000', ' -0.500', ' 1.200'))
        self.assertRaises(exceptions.PDBError, lambda: formats.PQRFile.parse(fn))
        with open(fn, 'w') as f:
            f.write('MODEL   1\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000', '  1.000', ' -0.500', ' 1.200'))
            f.write(self.ATOMLINE3 % (2, 'CB', 'ALA', 1, '  2.100', '  2.100', '  2.100', ' -0.500', ' 1.200'))
            f.write('MODEL   2\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000', '  1.000', ' -0.500', ' 1.200'))
            f.write('MODEL   3\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  2.000', '  2.000', '  2.000', ' -0.500', ' 1.200'))
        self.assertRaises(exceptions.PDBError, lambda: formats.PQRFile.parse(fn))
        with open(fn, 'w') as f:
            f.write('MODEL   1\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000', '  1.000', ' -0.500', ' 1.200'))
            f.write(self.ATOMLINE3 % (2, 'CB', 'ALA', 1, '  2.100', '  2.100', '  2.100', ' -0.500', ' 1.200'))
            f.write('MODEL   2\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  1.000', '  1.000', '  1.000', ' -0.500', ' 1.200'))
            f.write(self.ATOMLINE3 % (2, 'CB', 'ALA', 1, '  2.100', '  2.100', '  2.100', ' -0.500', ' 1.200'))
            f.write('MODEL   3\n')
            f.write(self.ATOMLINE3 % (1, 'CA', 'ALA', 1, '  2.000', '  2.000', '  2.000', ' -0.500', ' 1.200'))
        self.assertRaises(exceptions.PDBError, lambda: formats.PQRFile.parse(fn))

    def _check_adk_pqr(self, pqr):
        self.assertIsInstance(pqr, Structure)
        self.assertEqual(len(pqr.atoms), 3341)
        self.assertEqual(len(pqr.residues), 214)
        self.assertAlmostEqual(sum(a.charge for a in pqr.atoms), -4, places=4)
        self.assertEqual(pqr.atoms[0].charge, -0.30)
        self.assertEqual(pqr.atoms[0].solvent_radius, 1.85)
        self.assertEqual(pqr.atoms[0].atomic_number, 7)
        self.assertEqual(pqr.atoms[35].charge, -0.8)
        self.assertEqual(pqr.atoms[-1].charge, -0.67)
        self.assertEqual(pqr.atoms[-1].solvent_radius, 1.7)
        self.assertEqual(pqr.atoms[-1].atomic_number, 8)

    def test_pqr_writer(self):
        """ Tests writing a PQR file with charges and radii """
        parm = formats.load_file(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        fn = self.get_fn('test.pqr', written=True)
        # Create multiple models
        coords = []
        coords.append(parm.coordinates)
        coords.append(parm.coordinates + 1)
        coords.append(parm.coordinates + 2)
        parm.coordinates = np.vstack(coords)
        formats.PQRFile.write(parm, fn, renumber=True)
        pqr = formats.PQRFile.parse(fn)
        self.assertEqual(len(parm.atoms), len(pqr.atoms))
        for a1, a2 in zip(parm.atoms, pqr.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.idx, a2.residue.idx)
            self.assertAlmostEqual(a1.charge, a2.charge)
            self.assertAlmostEqual(a1.solvent_radius, a2.solvent_radius)
        self.assertEqual(pqr.get_coordinates().shape[0], 3)
        np.testing.assert_allclose(pqr.get_coordinates(0),
                                   parm.get_coordinates(0), atol=2e-3)
        # Pass coordinates explicitly
        formats.PQRFile.write(parm, fn, coordinates=np.vstack(coords[:2]))
        self.assertEqual(formats.PQRFile.parse(fn).get_coordinates().shape, (2, len(parm.atoms), 3))
        with self.assertRaises(TypeError):
            formats.PQRFile.write(parm, fn, coordinates=[1, 2, 3])

    def test_pqr_standard_resnames(self):
        """ Test standard residue name replacement in PQR writing """
        struct = Structure()
        a = Atom(name='CA', atomic_number=6, charge=0.5, solvent_radius=1.2)
        struct.add_atom(a, 'ASH', 2, 'A')
        fobj = StringIO()
        formats.PQRFile.write(struct, fobj, standard_resnames=True,
                              coordinates=[1, 1, 1], renumber=False)
        fobj.seek(0)
        self.assertEqual(formats.PQRFile.parse(fobj).residues[0].name, 'ASP')

    def test_pqr_with_element(self):
        """ Tests reading a PQR file that has an element column """
        self.assertTrue(formats.PQRFile.id_format(get_fn('elem.pqr')))
        pqr = formats.PQRFile.parse(get_fn('elem.pqr'))
        self.assertEqual(len(pqr.atoms), 458)
        self.assertEqual(len(pqr.residues), 14)
        self.assertEqual(pqr.atoms[0].charge, -0.9526)
        self.assertEqual(pqr.atoms[-1].solvent_radius, 0.8)

    def test_pqr_format_detection(self):
        """ Tests PDB file detection from contents """
        fn = self.get_fn('test.pdb', written=True)
        pdbtext1 = "%s%d    %9.6f %9.6f %9.6f     %10.5f\n" + self.ATOMLINE
        with open(fn, 'w') as f:
            f.write(pdbtext1 % ('ORIGX', 1, 10, 10, 10, 10, 1, 'CA', 'ALA',
                'A', 1, '   1.000', '  1.000', '  1.000', '%.4f' % -0.5,
                '%.4f' % 1.5))
        self.assertTrue(formats.PQRFile.id_format(fn))
        # Check failures
        with open(fn, 'w') as f:
            f.write(pdbtext1 % ('ORIGX', 8, 10, 10, 10, 10, 1, 'CA', 'ALA',
                'A', 1, '   1.000', '  1.000', '  1.000', '%.4f' % -0.5,
                '%.4f' % 1.5))
        self.assertFalse(formats.PQRFile.id_format(fn))
        with open(fn, 'w') as f:
            f.write(pdbtext1 % ('ORIGX', 1, 10, 10, 10, 10, 1, 'CA', 'ALA',
                'A', 1, '   1.000', '  1.000', '  1.000', '', ''))
        self.assertFalse(formats.PQRFile.id_format(fn))
        with open(fn, 'w') as f:
            f.write(pdbtext1 % ('ORIGX', 1, 10, 10, 10, 10, 1, 'CA', 'ALA',
                'A', 1, '   1.000', '  1.000', '  a.000', '%.4f' % -0.5,
                '%.4f' % 1.5))
        self.assertFalse(formats.PQRFile.id_format(fn))
        with open(fn, 'w') as f:
            pass
        self.assertFalse(formats.PQRFile.id_format(fn))

class TestCIFStructure(FileIOTestCase):

    def setUp(self):
        self.lztpdb = get_fn('4lzt.pdb')
        self.lzt = get_fn('4LZT.cif')
        self.largecif = get_fn('1ffk.cif')
        super().setUp()

    def test_write_cif(self):
        """ Test CIF writing capabilities """
        cif = read_CIF(self.lzt)
        written = self.get_fn('test.cif', written=True)
        with self.assertRaises(TypeError):
            cif.write_cif(written, coordinates=[1, 2, 3])
        cif.write_cif(written, renumber=False, write_anisou=True)
        cif2 = read_CIF(written)
        # cif and cif2 should have equivalent atom properties (basically,
        # everything should be the same except the metadata)
        self.assertEqual(len(cif.atoms), len(cif2.atoms))
        self.assertEqual(len(cif.residues), len(cif2.residues))

        # Check residue properties
        for res1, res2 in zip(cif.residues, cif2.residues):
            self.assertEqual(len(res1), len(res2))
            self.assertEqual(res1.name, res2.name)
            self.assertEqual(res1.chain, res2.chain)
            self.assertEqual(res1.insertion_code, res2.insertion_code)
            self.assertEqual(res1.number, res2.number)

        # Check atom properties
        for a1, a2 in zip(cif.atoms, cif2.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.number, a2.number)
            self.assertEqual(a1.element, a2.element)
            self.assertEqual(a1.altloc, a2.altloc)
            self.assertEqual(a1.xx, a2.xx)
            self.assertEqual(a1.xy, a2.xy)
            self.assertEqual(a1.xz, a2.xz)
            self.assertEqual(len(a1.anisou), 6)
            self.assertEqual(len(a2.anisou), 6)
            for x, y in zip(a1.anisou, a2.anisou):
                self.assertEqual(x, y)

        # Check box
        self.assertEqual(len(cif.box), len(cif2.box))
        for x, y in zip(cif.box, cif2.box):
            self.assertEqual(x, y)

        # Now check CIF writing without anisotropic B-factors and with
        # renumbering
        io = StringIO()
        cif.write_cif(io, coordinates=cif.get_coordinates('all'))
        io.seek(0)
        cif3 = read_CIF(io)
        # cif and cif3 should have equivalent atom properties (basically,
        # everything should be the same except the metadata)
        self.assertEqual(len(cif.atoms), len(cif3.atoms))
        self.assertEqual(len(cif.residues), len(cif3.residues))

        # Check residue properties
        i = 1
        for res1, res2 in zip(cif.residues, cif3.residues):
            self.assertEqual(len(res1), len(res2))
            self.assertEqual(res1.name, res2.name)
            self.assertEqual(res1.chain, res2.chain)
            self.assertEqual(res1.insertion_code, res2.insertion_code)
            self.assertEqual(res2.number, i)
            i += 1

        # Check atom properties
        i = 1
        for a1, a2 in zip(cif.atoms, cif3.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a2.number, i)
            self.assertEqual(a1.element, a2.element)
            self.assertEqual(a1.altloc, a2.altloc)
            self.assertEqual(a1.xx, a2.xx)
            self.assertEqual(a1.xy, a2.xy)
            self.assertEqual(a1.xz, a2.xz)
            self.assertIs(a2.anisou, None)
            i += 1 + len(a2.other_locations)

        # Check box
        self.assertEqual(len(cif.box), len(cif3.box))
        for x, y in zip(cif.box, cif3.box):
            self.assertEqual(x, y)

    def test_cif_detection(self):
        """ Tests CIF file auto-detection """
        fn = self.get_fn('test.cif', written=True)
        with open(fn, 'w') as f:
            pass
        self.assertFalse(formats.CIFFile.id_format(fn))

    def test_4lzt(self):
        """ Test CIF parsing on 4LZT (w/ ANISOU, altlocs, etc.) """
        self._check4lzt(read_CIF(self.lzt))

    @unittest.skipUnless(is_jenkins(), 'PDB blocks Travis from downloading files')
    def test_download(self):
        """ Test CIF downloading on 4LZT """
        fn = self.get_fn('4lzt.cif', written=True)
        self._check4lzt(download_CIF('4lzt', saveto=fn))
        self._check4lzt(read_CIF(fn))
        self.assertRaises(ValueError, lambda: download_CIF('illegal'))
        self.assertRaises(IOError, lambda: download_CIF('#@#%'))

    def test_cif_symmetry(self):
        """ Tests that symmetry is parsed from mmCIF files correctly """
        self.assertEqual(read_CIF(get_fn('1aki.cif')).space_group, 'P 21 21 21')

    def test_cif_space_group_written_from_structure(self):
        """ Tests CIF file writing with space groups """
        parm = pmd.load_file(get_fn('SCM_A.pdb'))
        self.assertEqual(parm.space_group, 'P 1 21 1')
        written = self.get_fn('test.cif', written=True)
        parm.write_cif(written)
        parm2 = pmd.load_file(written)
        self.assertEqual(parm2.space_group, 'P 1 21 1')

    def test_cif_models(self):
        """ Test CIF parsing/writing NMR structure with 20 models (2koc) """
        cif = read_CIF(get_fn('2koc.cif'))
        self.assertEqual(cif.get_coordinates('all').shape, (20, 451, 3))
        self.assertEqual(len(cif.atoms), 451)
        output = StringIO()
        cif.write_cif(output)
        output.seek(0)
        pdbfile2 = read_CIF(output)
        self.assertEqual(len(pdbfile2.atoms), 451)
        self.assertEqual(pdbfile2.get_coordinates('all').shape, (20, 451, 3))
        np.testing.assert_allclose(pdbfile2.get_coordinates('all'),
                                   cif.get_coordinates('all'))
        # Now check parsing and error handling for sample CIF files with
        # multiple models
        cif = formats.CIFFile.parse(get_fn('model.cif'))
        self.assertEqual(len(cif.atoms), 23)
        self.assertEqual(cif.get_coordinates().shape, (3, 23, 3))
        self.assertRaises(ValueError, lambda:
                formats.CIFFile.parse(get_fn('model_error1.cif'))
        )
        self.assertRaises(ValueError, lambda:
                formats.CIFFile.parse(get_fn('model_error2.cif'))
        )
        self.assertRaises(ValueError, lambda:
                formats.CIFFile.parse(get_fn('model_error3.cif'))
        )

    def test_cif_multiple_molecules(self):
        """ Test parsing CIF files with multiple molecules defined """
        # Create a composite CIF file from sample.cif and models.cif (both small
        # files). sample.cif has an extra anisotropic B-factor that is used for
        # error detection. It is the last line of the file, so discard it.
        fn = self.get_fn('test.cif', written=True)
        with open(get_fn('sample.cif'), 'r') as sf, \
                open(get_fn('model.cif'), 'r') as mf, open(fn, 'w') as f:
            nlines = sum(1 for line in sf)
            sf.seek(0)
            for line in range(nlines-1):
                f.write(sf.readline())
            f.write('\n\n')
            f.write(mf.read())
        sample, models = formats.CIFFile.parse(fn)

    def test_cif_write_standard_names(self):
        """ Test PDBx/mmCIF file writing converting to standard names """
        parm = formats.load_file(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        output = StringIO()
        parm.write_cif(output, standard_resnames=True)
        output.seek(0)
        pdb = read_CIF(output)
        for res in pdb.residues:
            self.assertEqual(
                    residue.AminoAcidResidue.get(res.name).abbr, res.name
            )

    def test_parse_cif_element_determination(self):
        """ Test element assignment for CIF files with bad element symbols """
        with self.assertWarns(exceptions.PDBWarning):
            formats.CIFFile.parse(get_fn('sample.cif'))
        cif = formats.CIFFile.parse(get_fn('sample.cif'))
        self.assertEqual(cif[0].atomic_number, 30) # element XX, atom name ZN
        self.assertEqual(cif[1].atomic_number, 6)  # element XX, atom name CA
        self.assertEqual(cif[2].atomic_number, 0)  # element XX, atom name ZZ
        self.assertIsInstance(cif[3], topologyobjects.ExtraPoint)

    def test_cif_altloc_writing(self):
        """ Tests alternate location handling in CIF files upon writing """
        struct = Structure()
        a1 = Atom(name='CA', atomic_number=6, altloc='A', occupancy=0.2)
        a2 = Atom(name='CA', atomic_number=6, altloc='B', occupancy=0.3)
        a3 = Atom(name='CA', atomic_number=6, altloc='C', occupancy=0.5)
        a1.other_locations['B'] = a2
        a1.other_locations['C'] = a3
        a1.xx, a1.xy, a1.xz = 1, 1, 1
        a2.xx, a2.xy, a2.xz = 2, 2, 2
        a3.xx, a3.xy, a3.xz = 3, 3, 3
        struct.add_atom(a1, 'ALA', 'A')
        fobj = StringIO()
        struct.write_cif(fobj, altlocs='all')
        fobj.seek(0)
        pdb = formats.CIFFile.parse(fobj)
        self.assertEqual(len(pdb.atoms), 1)
        self.assertEqual(pdb.atoms[0].occupancy, 0.2)
        self.assertEqual(pdb.atoms[0].other_locations['B'].occupancy, 0.3)
        self.assertEqual(pdb.atoms[0].other_locations['C'].occupancy, 0.5)
        fobj = StringIO()
        struct.write_cif(fobj, altlocs='first')
        fobj.seek(0)
        pdb = formats.CIFFile.parse(fobj)
        self.assertEqual(len(pdb.atoms), 1)
        self.assertEqual(pdb.atoms[0].occupancy, 0.2)
        self.assertEqual(pdb.atoms[0].altloc, 'A')
        self.assertEqual(pdb.atoms[0].xx, 1)
        self.assertEqual(pdb.atoms[0].xy, 1)
        self.assertEqual(pdb.atoms[0].xz, 1)
        fobj = StringIO()
        struct.write_cif(fobj, altlocs='occupancy')
        fobj.seek(0)
        pdb = formats.CIFFile.parse(fobj)
        self.assertEqual(len(pdb.atoms), 1)
        self.assertEqual(pdb.atoms[0].occupancy, 0.5)
        self.assertEqual(pdb.atoms[0].altloc, 'C')
        self.assertEqual(pdb.atoms[0].xx, 3)
        self.assertEqual(pdb.atoms[0].xy, 3)
        self.assertEqual(pdb.atoms[0].xz, 3)
        # Bad input
        with self.assertRaises(ValueError):
            struct.write_cif(self.get_fn('test.cif', written=True), altlocs='bad')

    def _check4lzt(self, cif):
        pdb = read_PDB(self.lztpdb)
        self.assertEqual(len(cif.atoms), len(pdb.atoms))
        nextra = 0
        for a1, a2 in zip(cif.atoms, pdb.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.number + nextra, a2.number)
            self.assertEqual(len(a1.anisou), len(a2.anisou))
            for x, y in zip(a1.anisou, a2.anisou):
                self.assertEqual(x, y)
            self.assertEqual(a1.altloc, a2.altloc)
            self.assertEqual(len(a1.other_locations), len(a2.other_locations))
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.number, a2.residue.number)
            # TER cards consume a serial number in the PDB file, but *not* in a
            # CIF file.
            if a2.residue.ter and a2 is a2.residue.atoms[-1]:
                nextra += 1
        # Check the unit cell info
        self.assertEqual(cif.box[0], 27.240)
        self.assertEqual(cif.box[1], 31.870)
        self.assertEqual(cif.box[2], 34.230)
        self.assertEqual(cif.box[3], 88.520)
        self.assertEqual(cif.box[4], 108.53)
        self.assertEqual(cif.box[5], 111.89)
        # Check the metadata now
        self.assertEqual(cif.experimental, 'X-RAY DIFFRACTION')
        self.assertEqual(cif.authors,
                'Walsh, M.A., Schneider, T., Sieker, L.C., Dauter, Z., '
                'Lamzin, V., Wilson, K.S.')
        self.assertEqual(cif.title,
                'Refinement of triclinic hen egg-white lysozyme at atomic '
                'resolution.; Refinement of Triclinic Lysozyme: I. Fourier '
                'and Least-Squares Methods; Refinement of Triclinic Lysozyme: '
                'II. The Method of Stereochemically Restrained Least Squares')
        self.assertEqual(cif.journal,
                'Acta Crystallogr.,Sect.D; Acta Crystallogr.,Sect.B; '
                'Acta Crystallogr.,Sect.B')
        self.assertEqual(cif.journal_authors,
                'Walsh, M.A., Schneider, T.R., Sieker, L.C., Dauter, Z., '
                'Lamzin, V.S., Wilson, K.S., Hodsdon, J.M., Brown, G.M., '
                'Jensen, L.H., Ramanadham, M.')
        self.assertEqual(cif.year, '1998, 1990, 1990')
        self.assertEqual(cif.page, '522, 54, 63')
        self.assertEqual(cif.keywords, ['HYDROLASE', 'O-GLYCOSYL',
                                        'GLYCOSIDASE'])
        self.assertEqual(cif.volume, '54, 46, 46')
        self.assertEqual(cif.doi, '10.1107/S0907444997013656')
        self.assertEqual(cif.pmid, '9761848')
        self.assertEqual(cif.resolution, 0.95)

class TestMol2File(FileIOTestCase):
    """ Tests the correct parsing and processing of mol2 files """

    def test_multi_mol2(self):
        """ Tests the parsing of a mol2 file with multiple residues """
        cont = formats.Mol2File.parse(get_fn('test_multi.mol2'))
        self.assertIsInstance(cont, ResidueTemplateContainer)
        self.assertEqual(len(cont), 20)
        for i, res in enumerate(cont):
            if i == 0:
                self.assertEqual(res.name, 'DA5')
                self.assertIs(res.head, None)
                self.assertIs(res.tail, [a for a in res if a.name == "O3'"][0])
                self.assertEqual(len(res), 30)
            elif i < 9:
                self.assertEqual(res.name, 'DA')
                self.assertIs(res.head, [a for a in res if a.name == "P"][0])
                self.assertIs(res.tail, [a for a in res if a.name == "O3'"][0])
                self.assertEqual(len(res), 32)
            elif i == 9:
                self.assertEqual(res.name, 'DA3')
                self.assertIs(res.head, [a for a in res if a.name == "P"][0])
                self.assertIs(res.tail, None)
                self.assertEqual(len(res), 33)
            elif i == 10:
                self.assertEqual(res.name, 'DT5')
                self.assertIs(res.head, None)
                self.assertIs(res.tail, [a for a in res if a.name == "O3'"][0])
                self.assertEqual(len(res), 30)
            elif i < 19:
                self.assertEqual(res.name, 'DT')
                self.assertIs(res.head, [a for a in res if a.name == "P"][0])
                self.assertIs(res.tail, [a for a in res if a.name == "O3'"][0])
                self.assertEqual(len(res), 32)
            elif i == 19:
                self.assertEqual(res.name, 'DT3')
                self.assertIs(res.head, [a for a in res if a.name == "P"][0])
                self.assertIs(res.tail, None)
                self.assertEqual(len(res), 33)
        # There are 686 total bonds defined in the mol2 file. However, there are
        # 20 residues, of which only 2 of the termini are *not* bonded to the
        # next residue. Ergo, 18 of those 686 bonds are intra-residue, and so
        # are not added to any of the ResidueTemplate instances. As a result, we
        # should have 686 - 18 = 668 total bonds if we add up the lengths of the
        # bond arrays for every residue
        self.assertEqual(sum([len(x.bonds) for x in cont]), 668)

    def test_multiple_mol2_entries(self):
        """ Tests a mol2 file with multiple @<MOLECULE> sections """
        cont = formats.Mol2File.parse(get_fn('multimol.mol2'))
        self.assertIsInstance(cont, ResidueTemplateContainer)
        self.assertEqual(len(cont), 200)
        for i, res in enumerate(cont):
            self.assertEqual(res.name, 'ZINC00000016_%d' % (i+1))
            self.assertEqual(len(res.atoms), 37)
            self.assertEqual(len(res.bonds), 38)
            self.assertIs(res.head, None)
            self.assertIs(res.tail, None)

    def test_multi_mol2_structure(self):
        """ Tests parsing a multi-residue mol2 into a Structure """
        struct = formats.Mol2File.parse(get_fn('test_multi.mol2'),
                                        structure=True)
        self.assertIsInstance(struct, Structure)
        self.assertEqual(len(struct.residues), 20)
        for i, res in enumerate(struct.residues):
            if i == 0:
                self.assertEqual(res.name, 'DA5')
                self.assertEqual(len(res), 30)
            elif i < 9:
                self.assertEqual(res.name, 'DA')
                self.assertEqual(len(res), 32)
            elif i == 9:
                self.assertEqual(res.name, 'DA3')
                self.assertEqual(len(res), 33)
            elif i == 10:
                self.assertEqual(res.name, 'DT5')
                self.assertEqual(len(res), 30)
            elif i < 19:
                self.assertEqual(res.name, 'DT')
                self.assertEqual(len(res), 32)
            elif i == 19:
                self.assertEqual(res.name, 'DT3')
                self.assertEqual(len(res), 33)
        # This should have the total number of bonds, 686
        self.assertEqual(len(struct.bonds), 686)

    def test_mol3_file(self):
        """ Tests parsing a Mol3 file into a ResidueTemplate """
        mol3 = formats.Mol2File.parse(get_fn('tripos9.mol2'))
        self.assertIsInstance(mol3, ResidueTemplate)
        self.assertEqual(mol3.name, 'GPN')
        self.assertEqual(len(mol3), 34)
        self.assertEqual(len(mol3.bonds), 35)
        self.assertIs(mol3.head, [a for a in mol3 if a.name == "N1'"][0])
        self.assertIs(mol3.tail, [a for a in mol3 if a.name == "C'"][0])
        # Test bad mol3 file
        self.assertRaises(exceptions.Mol2Error, lambda:
                formats.Mol2File.parse(get_fn('error.mol3'))
        )
        self.assertRaises(exceptions.Mol2Error, lambda:
                formats.Mol2File.parse(get_fn('error2.mol3'))
        )
        self.assertRaises(exceptions.Mol2Error, lambda:
                formats.Mol2File.parse(get_fn('error3.mol3'))
        )

    def test_mol3_file2(self):
        """ Tests parsing a Mol3 file with RESIDUECONNECT atoms """
        mol3 = formats.Mol2File.parse(get_fn('m2-c1_f3.mol2'))
        self.assertEqual(len(mol3.atoms), 27)
        self.assertEqual(len(mol3.bonds), 29)
        self.assertIs(mol3.head, None)
        self.assertIs(mol3.tail, None)
        self.assertIs(mol3.connections[0], mol3[5])
        self.assertIs(mol3.connections[1], mol3[9])

    def test_mol2_file_with_blank_lines(self):
        """ Tests parsing a Mol2 file with blank lines at the end """
        mol2 = formats.Mol2File.parse(get_fn('tripos1.mol2'))
        self.assertIsInstance(mol2, ResidueTemplate)
        self.assertEqual(mol2.name, 'DAN')
        self.assertEqual(len(mol2), 31)
        self.assertEqual(len(mol2.bonds), 33)

    def test_mol2_file_with_no_type_names(self):
        """ Tests writing a Mol2 without types uses names instead """
        struct = read_PDB(get_fn('2koc.pdb'))
        output = StringIO()
        formats.Mol2File.write(struct, output)
        output.seek(0)
        mol2 = formats.Mol2File.parse(output, structure=True)

    def test_mol3_structure(self):
        """ Tests parsing a Mol3 file with 1 residue into a Structure """
        mol3 = formats.Mol2File.parse(get_fn('tripos9.mol2'), structure=True)
        self.assertIsInstance(mol3, Structure)
        self.assertEqual(mol3.residues[0].name, 'GPN')
        self.assertEqual(len(mol3.atoms), 34)
        self.assertEqual(len(mol3.bonds), 35)

    def test_mol2_multi_write(self):
        """
        Tests writing mol2 file of multi residues from ResidueTemplateContainer
        """
        mol2 = formats.Mol2File.parse(get_fn('test_multi.mol2'))
        formats.Mol2File.write(mol2, self.get_fn('test_multi.mol2', written=True))
        fn = self.get_fn('test_multi_sep.mol2', written=True)
        formats.Mol2File.write(mol2, fn, split=True)
        fnsq = self.get_fn('test_multi_sep_squashed.mol2', written=True)
        formats.Mol2File.write(mol2, fnsq, split=True, compress_whitespace=True)
        self.assertTrue(diff_files(self.get_fn('test_multi.mol2', saved=True),
                                   self.get_fn('test_multi.mol2', written=True)))
        self.assertTrue(diff_files(self.get_fn('test_multi_sep.mol2', saved=True), fn))
        # Make sure the squashed lines all fall below 80 characters
        with open(fnsq) as f, open(fn) as f2:
            for line1, line2 in zip(f, f2):
                self.assertLessEqual(len(line1), 80)
                self.assertEqual(' '.join(line2.split()), line1.strip())
        mol22 = formats.Mol2File.parse(fn)
        self.assertEqual(len(mol2), len(mol22))
        self.assertEqual([r.name for r in mol2], [r.name for r in mol22])
        for r1, r2 in zip(mol2, mol22):
            self.assertEqual(len(r1.bonds), len(r2.bonds))
            self.assertEqual(len(r1.atoms), len(r2.atoms))
            self.assertFalse(r1.head is None and r1.tail is None)
            self.assertTrue(r2.head is None and r2.tail is None)
        f = StringIO()
        formats.Mol2File.write(mol2, f, mol3=True, split=True)
        f.seek(0)
        mol3 = formats.Mol2File.parse(f)
        self.assertEqual(len(mol2), len(mol3))
        self.assertEqual([r.name for r in mol2], [r.name for r in mol3])
        for r1, r2 in zip(mol2, mol3):
            self.assertEqual(len(r1.bonds), len(r2.bonds))
            self.assertEqual(len(r1.atoms), len(r2.atoms))
            self.assertFalse(r1.head is None and r1.tail is None)
            self.assertFalse(r2.head is None and r2.tail is None)
            if r1.head is None:
                self.assertIs(r2.head, None)
            else:
                self.assertEqual(r1.head.name, r2.head.name)
            if r1.tail is None:
                self.assertIs(r2.tail, None)
            else:
                self.assertEqual(r1.tail.name, r2.tail.name)

    def test_mol2_multi_write_from_structure(self):
        """ Tests writing mol2 file of multi residues from Structure """
        mol2 = formats.Mol2File.parse(get_fn('test_multi.mol2'), structure=True)
        formats.Mol2File.write(mol2, self.get_fn('test_multistruct.mol2', written=True))
        self.assertTrue(diff_files(self.get_fn('test_multistruct.mol2', written=True),
                                   self.get_fn('test_multistruct.mol2', saved=True)))
        fn = self.get_fn('test_splitmultistruct.mol2', written=True)
        formats.Mol2File.write(mol2, fn, split=True)
        self.assertTrue(diff_files(fn, self.get_fn('test_splitmultistruct.mol2', saved=True)))
        # Make residue names unrecognizable as amino or nucleic acids
        for i, res in enumerate(mol2.residues):
            res.name = '%03d' % i
        fobj = StringIO()
        formats.Mol2File.write(mol2, fobj)
        fobj.seek(0)
        self.assertIn('BIOPOLYMER', fobj.read())
        # Now delete *some* of the coordinates -- mol2 should fill with 0s
        fobj = StringIO()
        del mol2.atoms[0].xx, mol2.atoms[1].xy, mol2.atoms[2].xz
        formats.Mol2File.write(mol2, fobj)
        fobj.seek(0)
        read = formats.Mol2File.parse(fobj, structure=True)
        self.assertEqual(read.atoms[0].xx, 0)
        self.assertEqual(read.atoms[1].xy, 0)
        self.assertEqual(read.atoms[2].xz, 0)

    def test_mol3_multi_write(self):
        """
        Tests writing mol3 file of multi residues from ResidueTemplateContainer
        """
        mol2 = formats.Mol2File.parse(get_fn('test_multi.mol2'))
        formats.Mol2File.write(mol2, self.get_fn('test_multi.mol3', written=True), mol3=True)
        self.assertTrue(diff_files(self.get_fn('test_multi.mol3', written=True),
                                   self.get_fn('test_multi.mol3', saved=True)))

    def test_mol3_multi_write_from_structure(self):
        """ Tests writing mol3 file of multi residues from Structure """
        mol2 = formats.Mol2File.parse(get_fn('test_multi.mol2'), structure=True)
        formats.Mol2File.write(mol2, self.get_fn('test_multistruct.mol3', written=True), mol3=True)
        self.assertTrue(diff_files(self.get_fn('test_multistruct.mol3', written=True),
                                   self.get_fn('test_multistruct.mol3', saved=True)))

    def test_mol2_single_write(self):
        """ Tests writing mol2 file of single ResidueTemplate """
        mol2 = formats.Mol2File.parse(get_fn('tripos9.mol2'))
        formats.Mol2File.write(mol2, self.get_fn('tripos9.mol2', written=True))
        self.assertTrue(diff_files(self.get_fn('tripos9.mol2', written=True),
                                   self.get_fn('tripos9.mol2', saved=True)))

    def test_mol2_single_write_struct(self):
        """ Tests writing mol2 file of single-residue Structure """
        mol2 = formats.Mol2File.parse(get_fn('tripos9.mol2'), structure=True)
        self.assertIs(mol2.box, None)
        fn = self.get_fn('tripos9struct.mol2', written=True)
        formats.Mol2File.write(mol2, fn)
        self.assertTrue(diff_files(fn, self.get_fn('tripos9struct.mol2', saved=True)))
        self.assertIs(formats.Mol2File.parse(fn, structure=True).box, None)
        # Now add a box
        mol2.box = [10, 10, 10, 90, 90, 90]
        formats.Mol2File.write(mol2, fn)
        np.testing.assert_equal(formats.Mol2File.parse(fn, structure=True).box, [10, 10, 10, 90, 90, 90])

    def test_mol3_single_write(self):
        """ Tests writing mol3 file of single ResidueTemplate """
        mol2 = formats.Mol2File.parse(get_fn('tripos9.mol2'))
        formats.Mol2File.write(mol2, self.get_fn('tripos9.mol3', written=True), mol3=True)
        self.assertTrue(diff_files(self.get_fn('tripos9.mol3', written=True),
                                   self.get_fn('tripos9.mol3', saved=True)))
        # Now make sure it can write a ResidueTemplate with a connection
        mol2.connections.append(mol2.atoms[2])
        fobj = StringIO()
        formats.Mol2File.write(mol2, fobj, mol3=True)
        fobj.seek(0)
        new = formats.Mol2File.parse(fobj)
        self.assertEqual(len(new.connections), 1)
        self.assertIs(new.atoms[2], new.connections[0])

    def test_mol3_single_write_struct(self):
        """ Tests writing mol3 file of single-residue Structure """
        mol2 = formats.Mol2File.parse(get_fn('tripos9.mol2'), structure=True)
        formats.Mol2File.write(mol2, self.get_fn('tripos9struct.mol3', written=True), mol3=True)
        self.assertTrue(diff_files(self.get_fn('tripos9struct.mol3', written=True),
                                   self.get_fn('tripos9struct.mol3', written=True)))

    def test_mol2_atomic_number_assignment(self):
        """ Tests assignment of atomic numbers for mol2 files """
        mol2 = formats.Mol2File.parse(get_fn('tripos9.mol2'), structure=True)
        templ = formats.Mol2File.parse(get_fn('tripos9.mol2'))
        self.assertEqual(len(templ.atoms), len(mol2.atoms))
        for a1, a2 in zip(mol2.atoms, templ.atoms):
            self.assertEqual(a1.atomic_number, a2.atomic_number)
        # Now check that element assignment from GRO files (which has good
        # element assignment routines) matches what the mol2 does
        fn = self.get_fn('test.gro', written=True)
        mol2.save(fn, overwrite=True)
        for a1, a2 in zip(formats.load_file(fn).atoms, mol2.atoms):
            self.assertEqual(a1.atomic_number, a2.atomic_number)

    def test_mol2_box(self):
        """ Tests parsing Mol2 file with CRYSIN section """
        mol2 = formats.load_file(get_fn('tripos3.mol2'), structure=True)
        np.testing.assert_equal(mol2.box, [20, 20, 20, 90, 90, 90])

    def test_mol2_duplicate_atoms(self):
        """ Tests parsing Mol2 files with duplicate atoms """
        self.assertRaises(exceptions.Mol2Error, lambda:
                formats.Mol2File.parse(get_fn('duplicate_names.mol2')))
        mol2 = formats.Mol2File.parse(get_fn('duplicate_names.mol2'), structure=True)
        self.assertEqual(len(mol2.atoms), 89)
        self.assertEqual(len(mol2.bonds), 89)
        # Make sure that atom types are used to guess element if atom name is
        # ambiguous
        for atom in mol2.atoms:
            if atom.name == '****':
                self.assertEqual(atom.atomic_number, 1)

    def test_mol2_bond_order(self):
        """ Tests that mol2 file parsing remembers bond order/type """
        mol2 = formats.Mol2File.parse(get_fn('multimol.mol2'))[0]
        fn = self.get_fn('test.mol2', written=True)
        mol2.save(fn)
        with open(fn, 'r') as f:
            for line in f:
                if line.startswith('@<TRIPOS>BOND'):
                    break
            # Collect all bond orders
            bos = set()
            for line in f:
                if line.startswith('@<TRIPOS>'):
                    break
                bos.add(line.split()[3])
        # This structure has bond orders 1, 2, am, and ar
        self.assertEqual(bos, {'1', '2', 'am', 'ar'})
        mol2_2 = formats.Mol2File.parse(fn)
        for b1, b2 in zip(mol2.bonds, mol2_2.bonds):
            self.assertEqual(b1.order, b2.order)

    @unittest.skipUnless(HAS_GROMACS, 'Cannot test without gromacs')
    def test_mol3_disulfide(self):
        """ Tests writing mol3 file w/ disulfide (for RESIDUECONNECT) """
        top = formats.load_file(get_fn('1aki.ff99sbildn.top'))['!:SOL']
        top.coordinates = formats.load_file(get_fn('1aki.ff99sbildn.gro')).coordinates
        fn = self.get_fn('hewl.mol3', written=True)
        formats.Mol2File.write(top, fn, mol3=True)
        with open(fn, 'r') as f:
            for line in f:
                if ' '.join(line.split()) == '64 N C SG 0 0 0':
                    break
            else:
                assert False, 'Expected line not found'

class TestRegistry(FileIOTestCase):
    """ Tests properties of the FileFormatType registry """

    def test_file_format_type(self):
        """ Tests the FileFormatType metaclass """
        def create_metaclass():
            @add_metaclass(formats.registry.FileFormatType)
            class PDBFile(object):
                def id_format(fname):
                    return False
                def parse(fname):
                    raise NotImplementedError('Not implemented!')
            return PDBFile

        self.assertRaises(ValueError, create_metaclass)

    def test_load_file_errors(self):
        """ Test error handling in load_file """
        fn = self.get_fn('test.file', written=True)
        with open(fn, 'w') as f:
            pass
        os.chmod(fn, int('311', 8))
        self.assertRaises(IOError, lambda: formats.load_file(fn))

class TestFileDownloader(unittest.TestCase):
    """ Tests load_file with URLs for each format """

    def setUp(self):
        self.url = 'https://github.com/ParmEd/ParmEd/raw/master/test/files/'

    def test_download_off(self):
        """ Tests automatic loading of downloaded OFF files """
        self.assertTrue(amber.AmberOFFLibrary.id_format(self.url + 'amino12.lib'))
        off = amber.AmberOFFLibrary.parse(self.url + 'amino12.lib')
        self.assertIsInstance(off, dict)
        for key, item in iteritems(off):
            self.assertIsInstance(item, ResidueTemplate)

    def test_download_amber_prmtop(self):
        """ Tests automatic loading of downloaded AmberParm object """
        self.assertTrue(amber.AmberFormat.id_format(self.url + 'tip4p.parm7'))
        parm = amber.AmberFormat.parse(self.url + 'tip4p.parm7')
        self.assertIsInstance(parm, amber.AmberParm)

    def test_download_amoeba_prmtop(self):
        """ Tests automatic loading of downloaded AmoebaParm object """
        self.assertTrue(amber.AmberFormat.id_format(self.url + 'nma.parm7'))
        parm = amber.AmberFormat.parse(self.url + 'nma.parm7')
        self.assertIsInstance(parm, amber.AmoebaParm)

    def test_download_chamber_prmtop(self):
        """ Tests automatic loading of downloaded ChamberParm object """
        self.assertTrue(amber.AmberFormat.id_format(self.url + 'ala_ala_ala.parm7'))
        parm = amber.AmberFormat.parse(self.url + 'ala_ala_ala.parm7')
        self.assertIsInstance(parm, amber.ChamberParm)

    def test_download_raw_amber_format(self):
        """ Tests automatic loading of downloaded AmberFormat object """
        self.assertTrue(amber.AmberFormat.id_format(self.url + 'cSPCE.mdl'))
        parm = amber.AmberFormat.parse(self.url + 'cSPCE.mdl')
        self.assertIsInstance(parm, amber.AmberFormat)
        self.assertNotIsInstance(parm, amber.AmberParm)

    def test_download_amber_restart_ascii(self):
        """ Tests automatic loading of downloaded Amber ASCII restart file """
        self.assertTrue(amber.AmberAsciiRestart.id_format(self.url + 'trx.inpcrd'))
        parm = amber.AmberAsciiRestart(self.url + 'trx.inpcrd')
        self.assertIsInstance(parm, amber.AmberAsciiRestart)

    def test_download_amber_traj_ascii(self):
        """ Tests automatic loading of downloaded Amber mdcrd file """
        self.assertTrue(amber.AmberMdcrd.id_format(self.url + 'tz2.truncoct.crd'))
        crd = amber.AmberMdcrd(self.url + 'tz2.truncoct.crd', natom=5827, hasbox=True)
        self.assertIsInstance(crd, amber.AmberMdcrd)

    def test_download_charmm_psf(self):
        """ Tests automatic loading of downloaded CHARMM PSF file """
        self.assertTrue(formats.PSFFile.id_format(self.url + 'ala_ala_ala.psf'))
        parm = formats.PSFFile.parse(self.url + 'ala_ala_ala.psf')
        self.assertIsInstance(parm, charmm.CharmmPsfFile)

    def test_download_charmm_crd(self):
        """ Tests automatic loading of downloaded CHARMM crd file """
        self.assertTrue(charmm.CharmmCrdFile.id_format(self.url + 'dhfr_min_charmm.crd'))
        crd = charmm.CharmmCrdFile(self.url + 'dhfr_min_charmm.crd')
        self.assertIsInstance(crd, charmm.CharmmCrdFile)

    def test_download_charmm_restart(self):
        """ Tests automatic loading of downloaded CHARMM restart file """
        self.assertTrue(charmm.CharmmRstFile.id_format(self.url + 'sample-charmm.rst'))
        crd = charmm.CharmmRstFile(self.url + 'sample-charmm.rst')
        self.assertIsInstance(crd, charmm.CharmmRstFile)

    def test_download_pdb(self):
        """ Tests automatic loading of downloaded PDB files """
        self.assertTrue(formats.PDBFile.id_format(self.url + '4lzt.pdb'))
        pdb = formats.PDBFile.parse(self.url + '4lzt.pdb')
        self.assertIsInstance(pdb, Structure)
        self.assertEqual(len(pdb.atoms), 1164)

    def test_download_cif(self):
        """ Tests automatic loading of downloaded PDBx/mmCIF files """
        self.assertTrue(formats.CIFFile.id_format(self.url + '4LZT.cif'))
        cif = formats.CIFFile.parse(self.url + '4LZT.cif')
        self.assertIsInstance(cif, Structure)
        self.assertEqual(len(cif.atoms), 1164)

    def test_download_mol2(self):
        """ Tests automatic loading of downloaded mol2 and mol3 files """
        self.assertTrue(formats.Mol2File.id_format(self.url + 'test_multi.mol2'))
        self.assertTrue(formats.Mol2File.id_format(self.url + 'tripos9.mol2'))
        mol2 = formats.Mol2File.parse(self.url + 'test_multi.mol2')
        self.assertIsInstance(mol2, ResidueTemplateContainer)
        mol3 = formats.Mol2File.parse(self.url + 'tripos9.mol2')
        self.assertIsInstance(mol3, ResidueTemplate)

    @unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without GROMACS")
    def test_download_gromacs_topology(self):
        """ Tests automatic loading of downloaded Gromacs topology file """
        self.assertTrue(gromacs.GromacsTopologyFile.id_format(self.url + '1aki.charmm27.top'))
        top = gromacs.GromacsTopologyFile(self.url + '1aki.charmm27.top')
        self.assertIsInstance(top, gromacs.GromacsTopologyFile)

    def test_download_gromacs_gro(self):
        """ Tests automatic loading of downloaded Gromacs GRO file """
        self.assertTrue(gromacs.GromacsGroFile.id_format(self.url + '1aki.ff99sbildn.gro'))
        gro = gromacs.GromacsGroFile.parse(self.url + '1aki.ff99sbildn.gro')
        self.assertIsInstance(gro, Structure)

    def test_download_netcdf(self):
        """ Tests that NetCDF files always fail when trying to download them """
        self.assertFalse(amber.NetCDFRestart.id_format(self.url + 'ncinpcrd.rst7'))
        self.assertFalse(amber.NetCDFTraj.id_format(self.url + 'tz2.truncoct.nc'))
