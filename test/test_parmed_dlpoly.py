"""
Tests the functionality in the parmed.dlpoly package
"""
import parmed
from utils import get_fn, get_saved_fn, diff_files, FileIOTestCase

class TestDlpoly(FileIOTestCase):
    """ Tests the generation of DLPOLY input files """

    def test_ala_gas(self):
        """ Test for alanine dipeptide in gas-phase """
        parm = get_fn('ala_gas.prmtop')
        rst7 = get_fn('ala_gas.rst7')
        output1 = self.get_fn('ala_gas.field', written=True)
        output2 = self.get_fn('ala_gas.config', written=True)
        amber = parmed.load_file(parm,rst7)
        amber.save(output1)
        self.assertTrue(
            diff_files(self.get_fn('ala_gas.field', saved=True), output1, absolute_error=1e-6)
        )
        amber.save(output2)
        self.assertTrue(
            diff_files(self.get_fn('ala_gas.config', saved=True), output2, absolute_error=1e-6)
        )

    def test_ksi_gas(self):
      """ Test for ksi protein in gas-phase """
      parm = get_fn('ksi_gas.prmtop')
      rst7 = get_fn('ksi_gas.rst7')
      output1 = self.get_fn('ksi_gas.field', written=True)
      output2 = self.get_fn('ksi_gas.config', written=True)
      amber = parmed.load_file(parm,rst7)
      amber.save(output1)
      self.assertTrue(
          diff_files(self.get_fn('ksi_gas.field', saved=True), output1, absolute_error=1e-6)
      )
      amber.save(output2)
      self.assertTrue(
          diff_files(self.get_fn('ksi_gas.config', saved=True), output2, absolute_error=1e-6)
      )

    def test_ksi_solv(self):
      """ Test for ksi protein in solution-phase """
      parm = get_fn('ksi_solv.prmtop')
      rst7 = get_fn('ksi_solv.rst7')
      output1 = self.get_fn('ksi_solv.field', written=True)
      output2 = self.get_fn('ksi_solv.config', written=True)
      amber = parmed.load_file(parm,rst7)
      amber.save(output1)
      self.assertTrue(
          diff_files(self.get_fn('ksi_solv.field', saved=True), output1, absolute_error=1e-6)
      )
      amber.save(output2)
      self.assertTrue(
          diff_files(self.get_fn('ksi_solv.config', saved=True), output2, absolute_error=1e-6)
      )
