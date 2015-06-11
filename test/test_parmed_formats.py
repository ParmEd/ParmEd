"""
Tests parmed.formats package
"""
from __future__ import division
import utils

from parmed import amber, charmm, exceptions, formats, gromacs
from parmed import (Structure, read_PDB, write_PDB, read_CIF,
                    download_PDB, download_CIF)
from parmed.modeller import ResidueTemplate, ResidueTemplateContainer
from parmed.utils.six import iteritems
try:
    import cStringIO as StringIO
    from itertools import izip as zip
except ImportError:
    # Must be Python 3
    import io as StringIO
import os
import unittest
from utils import get_fn, has_numpy, diff_files, get_saved_fn, skip_big_tests

def reset_stringio(io):
    """ Resets a StringIO instance to "empty-file" state """
    io.seek(0)
    io.truncate()
    return io

class TestFileLoader(unittest.TestCase):
    """ Tests the automatic file loader """

    def testLoadOFF(self):
        """ Tests automatic loading of OFF files """
        off = formats.load_file(get_fn('amino12.lib'))
        self.assertIsInstance(off, dict)
        for key, item in iteritems(off):
            self.assertIsInstance(item, ResidueTemplate)

    def testLoadAmberParm(self):
        """ Tests automatic loading of AmberParm object """
        parm = formats.load_file(get_fn('trx.prmtop'))
        self.assertIsInstance(parm, amber.AmberParm)

    def testLoadAmoebaParm(self):
        """ Tests automatic loading of AmoebaParm object """
        parm = formats.load_file(get_fn('amoeba.parm7'))
        self.assertIsInstance(parm, amber.AmoebaParm)

    def testLoadChamberParm(self):
        """ Tests automatic loading of ChamberParm object """
        parm = formats.load_file(get_fn('ala_ala_ala.parm7'))
        self.assertIsInstance(parm, amber.ChamberParm)

    def testLoadAmberFormat(self):
        """ Tests automatic loading of RISM mdl (AmberFormat) object """
        parm = formats.load_file(get_fn('cSPCE.mdl'))
        self.assertIsInstance(parm, amber.AmberFormat)
        self.assertNotIsInstance(parm, amber.AmberParm)

    def testLoadAmberRestart(self):
        """ Tests automatic loading of Amber ASCII restart file """
        parm = formats.load_file(get_fn('trx.inpcrd'))
        self.assertIsInstance(parm, amber.AmberAsciiRestart)

    def testLoadAmberMdcrd(self):
        """ Tests automatic loading of Amber mdcrd file """
        crd = formats.load_file(get_fn('tz2.truncoct.crd'), natom=5827,
                                hasbox=True)
        self.assertIsInstance(crd, amber.AmberMdcrd)
        self.assertRaises(TypeError, lambda:
                formats.load_file(get_fn('tz2.truncoct.crd')))

    def testLoadCharmmPsfFile(self):
        """ Tests automatic loading of CHARMM PSF file """
        parm = formats.load_file(get_fn('ala_ala_ala.psf'))
        self.assertIsInstance(parm, charmm.CharmmPsfFile)

    def testLoadCharmmCrdFile(self):
        """ Tests automatic loading of CHARMM crd file """
        crd = formats.load_file(get_fn('dhfr_min_charmm.crd'))
        self.assertIsInstance(crd, charmm.CharmmCrdFile)

    def testLoadCharmmRestart(self):
        """ Tests automatic loading of CHARMM restart file """
        crd = formats.load_file(get_fn('sample-charmm.rst'))
        self.assertIsInstance(crd, charmm.CharmmRstFile)

    def testLoadNetCDFRestart(self):
        """ Tests automatic loading of Amber NetCDF restart file """
        if amber.HAS_NETCDF:
            crd = formats.load_file(get_fn('ncinpcrd.rst7'))
            self.assertIsInstance(crd, amber.NetCDFRestart)
        else:
            self.assertRaises(exceptions.FormatNotFound, lambda:
                    formats.load_file(get_fn('ncinpcrd.rst7')))

    def testLoadNetCDFTraj(self):
        """ Tests automatic loading of Amber NetCDF trajectory file """
        if amber.HAS_NETCDF:
            crd = formats.load_file(get_fn('tz2.truncoct.nc'))
            self.assertIsInstance(crd, amber.NetCDFTraj)
        else:
            self.assertRaises(exceptions.FormatNotFound, lambda:
                    formats.load_file(get_fn('ncinpcrd.rst7')))

    def testLoadPDB(self):
        """ Tests automatic loading of PDB files """
        pdb = formats.load_file(get_fn('4lzt.pdb'))
        self.assertIsInstance(pdb, Structure)
        self.assertEqual(len(pdb.atoms), 1164)

    def testLoadCIF(self):
        """ Tests automatic loading of PDBx/mmCIF files """
        cif = formats.load_file(get_fn('4LZT.cif'))
        self.assertIsInstance(cif, Structure)
        self.assertEqual(len(cif.atoms), 1164)

    def testLoadMol2(self):
        """ Tests automatic loading of mol2 and mol3 files """
        mol2 = formats.load_file(get_fn('test_multi.mol2'))
        self.assertIsInstance(mol2, ResidueTemplateContainer)
        mol3 = formats.load_file(get_fn('tripos9.mol2'))
        self.assertIsInstance(mol3, ResidueTemplate)

    def testLoadGromacsTop(self):
        """ Tests automatic loading of Gromacs topology file """
        top = formats.load_file(get_fn('1aki.charmm27.top'))
        self.assertIsInstance(top, gromacs.GromacsTopologyFile)

    def testLoadGromacsGro(self):
        """ Tests automatic loading of Gromacs GRO file """
        gro = formats.load_file(get_fn('1aki.ff99sbildn.gro'))
        self.assertIsInstance(gro, Structure)

    def testBadLoads(self):
        """ Test exception handling when non-recognized files are loaded """
        self.assertRaises(exceptions.FormatNotFound, lambda:
                formats.load_file(get_fn('../test_parmed_formats.py')))
        self.assertRaises(IOError, lambda: formats.load_file('no_file'))

class TestChemistryPDBStructure(unittest.TestCase):
    
    def setUp(self):
        self.pdb = get_fn('4lzt.pdb')
        self.pdbgz = get_fn('4lzt.pdb.gz')
        self.pdbbz2 = get_fn('4lzt.pdb.bz2')
        self.models = get_fn('2koc.pdb')
        self.overflow = get_fn('4lyt_vmd.pdb')
        self.simple = get_fn('ala_ala_ala.pdb')
        self.format_test = get_fn('SCM_A.pdb')
        self.overflow2 = get_fn('overflow.pdb')
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

    def testAscii(self):
        """ Test PDB file parsing """
        self._check4lzt(read_PDB(self.pdb))
        # The PDB file with multiple models
        pdbfile = read_PDB(open(self.models))
        self.assertEqual(len(pdbfile.pdbxyz), 20)
        self.assertEqual(pdbfile.pdbxyz[0][:3], [-8.886, -5.163, 9.647])
        self.assertEqual(pdbfile.pdbxyz[19][-3:], [-12.051, 5.205, -2.146])

    def testDownload(self):
        """ Tests downloading PDB files """
        self._check4lzt(download_PDB('4lzt'))

    def testPositions(self):
        """ Tests that positions are Vec3's with units """
        from parmed import unit as u
        from parmed import Vec3
        pdbfile = read_PDB(open(self.models))
        self.assertIsInstance(pdbfile.positions[0], u.Quantity)
        self.assertIsInstance(pdbfile.positions[0].value_in_unit(u.angstroms), Vec3)

    def testGzip(self):
        """ Test Gzipped-PDB file parsing """
        self._check4lzt(read_PDB(self.pdbgz))

    def testBzip(self):
        """ Test Bzipped-PDB file parsing """
        self._check4lzt(read_PDB(self.pdbbz2))

    @unittest.skipIf(skip_big_tests(), 'Skipping large tests')
    def testVmdOverflow(self):
        """ Test PDB file where atom and residue numbers overflow """
        pdbfile = read_PDB(self.overflow)
        self.assertEqual(len(pdbfile.atoms), 110237)
        self.assertEqual(len(pdbfile.residues), 35697)
        self.assertEqual(pdbfile.box, [0, 0, 0, 90, 90, 90])

    @unittest.skipIf(skip_big_tests(), 'Skipping large tests')
    def testRegularOverflow(self):
        """ Test PDB file where atom number goes to ***** after 99999 """
        pdbfile = read_PDB(self.overflow2)
        self.assertEqual(len(pdbfile.atoms), 114277)
        self.assertEqual(len(pdbfile.residues), 25042)
        for i, atom in enumerate(pdbfile.atoms):
            self.assertEqual(atom.number, i+1)
            self.assertEqual(atom.idx, i)

    def testPdbWriteSimple(self):
        """ Test PDB file writing on a very simple input structure """
        pdbfile = read_PDB(self.simple)
        self.assertEqual(len(pdbfile.atoms), 33)
        self.assertEqual(len(pdbfile.residues), 3)
        output = StringIO.StringIO()
        pdbfile.write_pdb(output)
        output.seek(0)
        pdbfile2 = read_PDB(output)
        self.assertEqual(len(pdbfile2.atoms), 33)
        self.assertEqual(len(pdbfile2.residues), 3)
        self._compareInputOutputPDBs(pdbfile, pdbfile2)

    def testPdbWriteModels(self):
        """ Test PDB file writing from NMR structure with models """
        pdbfile = read_PDB(self.models)
        self.assertEqual(len(pdbfile.pdbxyz), 20)
        self.assertEqual(len(pdbfile.atoms), 451)
        output = StringIO.StringIO()
        write_PDB(pdbfile, output)
        output.seek(0)
        pdbfile2 = read_PDB(output)
        self.assertEqual(len(pdbfile2.atoms), 451)
        self._compareInputOutputPDBs(pdbfile, pdbfile2)

    def testPdbWriteXtal(self):
        """ Test PDB file writing from a Xtal structure """
        pdbfile = read_PDB(self.pdb)
        self._check4lzt(pdbfile)
        output = StringIO.StringIO()
        pdbfile.write_pdb(output, renumber=False)
        output.seek(0)
        pdbfile2 = read_PDB(output)
        self._check4lzt(pdbfile2, check_meta=False)
        self._compareInputOutputPDBs(pdbfile, pdbfile2)
        output = reset_stringio(output)
        write_PDB(pdbfile, output)
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

    def testPdbWriteAltlocOptions(self):
        """ Test PDB file writing with different altloc options """
        pdbfile = read_PDB(self.pdb)
        self._check4lzt(pdbfile)
        output = StringIO.StringIO()
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
        write_PDB(pdbfile, output, renumber=False, altlocs='occupancy')
        output.seek(0)
        pdbfile4 = read_PDB(output)
        self._check4lzt(pdbfile4, check_meta=False, has_altloc=False)
        self._compareInputOutputPDBs(pdbfile, pdbfile4, altloc_option='occupancy')
        # Double-check 'first' vs. 'occupancy'. Residue 85 (SER) has a conformer
        # A that has an occupancy of 0.37 and conformer B with occupancy 0.63
        self.assertEqual(pdbfile3.residues[84][4].xx, -4.162)
        self.assertEqual(pdbfile4.residues[84][4].xx, -4.157)

    def testAnisouRead(self):
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

    def testAnisouWrite(self):
        """ Tests that write_PDB properly writes ANISOU records """
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
        output = StringIO.StringIO()
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
            if has_numpy():
                self.assertEqual(a1.anisou.shape, a2.anisou.shape)
            else:
                self.assertEqual(len(a1.anisou), len(a2.anisou))
            for x, y in zip(a1.anisou, a2.anisou):
                self.assertAlmostEqual(x, y, delta=1e-4)
            self.assertEqual(len(a1.other_locations), len(a2.other_locations))
            for key in sorted(a1.other_locations.keys()):
                oa1 = a1.other_locations[key]
                oa2 = a2.other_locations[key]
                if has_numpy():
                    self.assertEqual(oa1.anisou.shape, oa2.anisou.shape)
                else:
                    self.assertEqual(len(oa1.anisou), len(oa2.anisou))
                for x, y in zip(oa1.anisou, oa2.anisou):
                    self.assertAlmostEqual(x, y, delta=1e-4)

    def testPDBWriteFormat(self):
        """ Test PDB atom names are properly justified per PDB standard """
        pdbfile = read_PDB(self.format_test)
        f = get_fn('pdb_format_test.pdb', written=True)
        pdbfile.write_pdb(f, write_anisou=True)
        self.assertTrue(diff_files(get_saved_fn('SCM_A_formatted.pdb'), f))

    @unittest.skipIf(skip_big_tests(), 'Skipping large tests')
    def testSegidHandling(self):
        """ Test handling of CHARMM-specific SEGID identifier (r/w) """
        pdbfile = read_PDB(self.overflow2)
        allsegids = set(['PROA', 'PROB', 'CARA', 'CARE', 'CARC', 'CARD', 'CARB',
                         'MEMB', 'TIP3', 'POT', 'CLA'])
        foundsegids = set()
        for atom in pdbfile.atoms:
            self.assertTrue(hasattr(atom, 'segid'))
            foundsegids.add(atom.segid)
        self.assertEqual(foundsegids, allsegids)
        self.assertEqual(pdbfile.atoms[0].segid, 'PROA')
        self.assertEqual(pdbfile.atoms[5161].segid, 'PROA')
        self.assertEqual(pdbfile.atoms[5162].segid, 'PROB')
        self.assertEqual(pdbfile.atoms[-1].segid, 'CLA')
        f = get_fn('pdb_segid_test1.pdb', written=True)
        f2 = get_fn('pdb_segid_test2.pdb', written=True)
        pdbfile.write_pdb(f)
        pdbfile2 = read_PDB(f)
        for atom in pdbfile2.atoms:
            self.assertFalse(hasattr(atom, 'segid'))
        pdbfile.write_pdb(f2, charmm=True)
        pdbfile3 = read_PDB(f2)
        for atom in pdbfile3.atoms:
            self.assertTrue(hasattr(atom, 'segid'))
            self.assertEqual(atom.segid, pdbfile.atoms[atom.idx].segid)

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

    # Private helper test functions
    def _check4lzt(self, obj, check_meta=True, has_altloc=True):
        self.assertEqual(len(obj.pdbxyz), 1)
        self.assertEqual(obj.box,
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
        # Check the TER card is picked up
        for i, residue in enumerate(obj.residues):
            if i == 128:
                self.assertTrue(residue.ter)
            else:
                self.assertFalse(residue.ter)

class TestChemistryCIFStructure(unittest.TestCase):

    def setUp(self):
        self.lztpdb = get_fn('4lzt.pdb')
        self.lzt = get_fn('4LZT.cif')
        self.largecif = get_fn('1ffk.cif')
        try:
            os.makedirs(get_fn('writes'))
        except OSError:
            pass # Already exists

    def tearDown(self):
        try:
            for f in os.listdir(get_fn('writes')):
                os.unlink(get_fn(f, written=True))
            os.rmdir(get_fn('writes'))
        except OSError:
            pass

    def testWriteCIF(self):
        """ Test CIF writing capabilities """
        cif = read_CIF(self.lzt)
        written = get_fn('test.cif', written=True)
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
        io = StringIO.StringIO()
        cif.write_cif(io)
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

    def test4LZT(self):
        """ Test CIF parsing on 4LZT (w/ ANISOU, altlocs, etc.) """
        self._check4lzt(read_CIF(self.lzt))

    def testDownload(self):
        """ Test CIF downloading on 4LZT """
        self._check4lzt(download_CIF('4lzt'))

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

class TestMol2File(unittest.TestCase):
    """ Tests the correct parsing and processing of mol2 files """

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

    def testMultiMol2(self):
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

    def testMultiMol2Structure(self):
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

    def testMol3File(self):
        """ Tests parsing a Mol3 file into a ResidueTemplate """
        mol3 = formats.Mol2File.parse(get_fn('tripos9.mol2'))
        self.assertIsInstance(mol3, ResidueTemplate)
        self.assertEqual(mol3.name, 'GPN')
        self.assertEqual(len(mol3), 34)
        self.assertEqual(len(mol3.bonds), 35)
        self.assertIs(mol3.head, [a for a in mol3 if a.name == "N1'"][0])
        self.assertIs(mol3.tail, [a for a in mol3 if a.name == "C'"][0])

    def testMol3Structure(self):
        """ Tests parsing a Mol3 file with 1 residue into a Structure """
        mol3 = formats.Mol2File.parse(get_fn('tripos9.mol2'), structure=True)
        self.assertIsInstance(mol3, Structure)
        self.assertEqual(mol3.residues[0].name, 'GPN')
        self.assertEqual(len(mol3.atoms), 34)
        self.assertEqual(len(mol3.bonds), 35)

    def testMol2MultiWrite(self):
        """
        Tests writing mol2 file of multi residues from ResidueTemplateContainer
        """
        mol2 = formats.Mol2File.parse(get_fn('test_multi.mol2'))
        formats.Mol2File.write(mol2, get_fn('test_multi.mol2', written=True))
        self.assertTrue(diff_files(get_fn('test_multi.mol2', written=True),
                                   get_saved_fn('test_multi.mol2')))

    def testMol2MultiWriteFromStructure(self):
        """ Tests writing mol2 file of multi residues from Structure """
        mol2 = formats.Mol2File.parse(get_fn('test_multi.mol2'), structure=True)
        formats.Mol2File.write(mol2, get_fn('test_multistruct.mol2', written=True))
        self.assertTrue(diff_files(get_fn('test_multistruct.mol2', written=True),
                                   get_saved_fn('test_multistruct.mol2')))

    def testMol3MultiWrite(self):
        """
        Tests writing mol3 file of multi residues from ResidueTemplateContainer
        """
        mol2 = formats.Mol2File.parse(get_fn('test_multi.mol2'))
        formats.Mol2File.write(mol2, get_fn('test_multi.mol3', written=True),
                               mol3=True)
        self.assertTrue(diff_files(get_fn('test_multi.mol3', written=True),
                                   get_saved_fn('test_multi.mol3')))

    def testMol3MultiWriteFromStructure(self):
        """ Tests writing mol3 file of multi residues from Structure """
        mol2 = formats.Mol2File.parse(get_fn('test_multi.mol2'), structure=True)
        formats.Mol2File.write(mol2, get_fn('test_multistruct.mol3', written=True),
                               mol3=True)
        self.assertTrue(diff_files(get_fn('test_multistruct.mol3', written=True),
                                   get_saved_fn('test_multistruct.mol3')))

    def testMol2SingleWrite(self):
        """ Tests writing mol2 file of single ResidueTemplate """
        mol2 = formats.Mol2File.parse(get_fn('tripos9.mol2'))
        formats.Mol2File.write(mol2, get_fn('tripos9.mol2', written=True))
        self.assertTrue(diff_files(get_fn('tripos9.mol2', written=True),
                                   get_saved_fn('tripos9.mol2')))

    def testMol2SingleWriteStruct(self):
        """ Tests writing mol2 file of single-residue Structure """
        mol2 = formats.Mol2File.parse(get_fn('tripos9.mol2'), structure=True)
        formats.Mol2File.write(mol2, get_fn('tripos9struct.mol2', written=True))
        self.assertTrue(diff_files(get_fn('tripos9struct.mol2', written=True),
                                   get_saved_fn('tripos9struct.mol2')))

    def testMol3SingleWrite(self):
        """ Tests writing mol3 file of single ResidueTemplate """
        mol2 = formats.Mol2File.parse(get_fn('tripos9.mol2'))
        formats.Mol2File.write(mol2, get_fn('tripos9.mol3', written=True),
                               mol3=True)
        self.assertTrue(diff_files(get_fn('tripos9.mol3', written=True),
                                   get_saved_fn('tripos9.mol3')))

    def testMol3SingleWriteStruct(self):
        """ Tests writing mol3 file of single-residue Structure """
        mol2 = formats.Mol2File.parse(get_fn('tripos9.mol2'), structure=True)
        formats.Mol2File.write(mol2, get_fn('tripos9struct.mol3', written=True),
                               mol3=True)
        self.assertTrue(diff_files(get_fn('tripos9struct.mol3', written=True),
                                   get_saved_fn('tripos9struct.mol3')))

if __name__ == '__main__':
    unittest.main()
