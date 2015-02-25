"""
Tests chemistry.formats package
"""

from chemistry import amber
from chemistry import charmm
from chemistry import exceptions
from chemistry import formats
from chemistry import Structure
from chemistry.modeller import ResidueTemplate
import unittest
import utils
get_fn = utils.get_fn

class TestFileLoader(unittest.TestCase):
    """ Tests the automatic file loader """

    def testLoadOFF(self):
        """ Tests automatic loading of OFF files """
        off = formats.load_file(get_fn('amino12.lib'))
        self.assertIsInstance(off, dict)
        for key, item in off.iteritems():
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
        crd = formats.load_file(get_fn('ncinpcrd.rst7'))
        self.assertIsInstance(crd, amber.NetCDFRestart)

    def testLoadNetCDFTraj(self):
        """ Tests automatic loading of Amber NetCDF trajectory file """
        crd = formats.load_file(get_fn('tz2.truncoct.nc'))
        self.assertIsInstance(crd, amber.NetCDFTraj)

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

    def testBadLoads(self):
        """ Test exception handling when non-recognized files are loaded """
        self.assertRaises(exceptions.FormatNotFound, lambda:
                formats.load_file(get_fn('../test_chemistry_formats.py')))
        self.assertRaises(IOError, lambda: formats.load_file('no_file'))
