""" Tests some OpenMM-specific functionality """
from __future__ import print_function, division, absolute_import

import numpy as np
from parmed import openmm, load_file, exceptions, amber, charmm
from parmed.utils.six.moves import StringIO
import os
import unittest
from utils import get_fn, mm, app, has_openmm, FileIOTestCase

@unittest.skipIf(not has_openmm, "Cannot test without OpenMM")
class TestOpenMM(FileIOTestCase):

    def setUp(self):
        super(TestOpenMM, self).setUp()
        # Take one of the distributed OpenMM FF XML files as a test
        self.ffxml = os.path.join(os.path.split(app.__file__)[0], 'data',
                                  'amber99sbildn.xml')

    def testFormatID(self):
        """ Tests automatic format determination of OpenMM XML files """
        self.assertTrue(openmm.XmlFile.id_format(get_fn('system_974wat.xml')))
        self.assertTrue(openmm.XmlFile.id_format(get_fn('state_974wat.xml')))
        self.assertTrue(openmm.XmlFile.id_format(get_fn('integrator.xml')))
        self.assertTrue(openmm.XmlFile.id_format(self.ffxml))

    def testDeserializeSystem(self):
        """ Tests automatic deserialization of a System XML file """
        system = openmm.XmlFile.parse(get_fn('system_974wat.xml'))
        self.assertIsInstance(system, mm.System)
        self.assertEqual(system.getNumParticles(), 6638)

    def testDeserializeState(self):
        """ Tests automatic deserialization of a State XML file """
        state = openmm.XmlFile.parse(get_fn('state_974wat.xml'))
        self.assertEqual(state.coordinates.shape, (1, 6638, 3))
        self.assertEqual(state.velocities.shape, (1, 6638, 3))
        self.assertEqual(state.forces.shape, (1, 6638, 3))
        np.testing.assert_allclose(state.box, [35.05, 40.5, 42.37, 90, 90, 90])
        self.assertAlmostEqual(state.time, 20000.000003615783)
        self.assertIs(state.energy, None)

    def testDeserializeIntegrator(self):
        """ Tests automatic deserialization of an Integrator XML file """
        integrator = openmm.XmlFile.parse(get_fn('integrator.xml'))
        self.assertIsInstance(integrator, mm.Integrator)
        self.assertIsInstance(integrator, mm.LangevinIntegrator)

    def testDeserializeForceField(self):
        """ Tests automatic deserialization of an OpenMM ForceField XML file """
        ff = openmm.XmlFile.parse(self.ffxml)
        self.assertIsInstance(ff, app.ForceField)

    def testLoadTopology(self):
        """ Tests loading an OpenMM Topology and System instance """
        import warnings
        warnings.filterwarnings('error', category=exceptions.OpenMMWarning)
        ommparm = app.AmberPrmtopFile(get_fn('complex.prmtop'))
        parm = load_file(get_fn('complex.prmtop'))
        system = ommparm.createSystem(implicitSolvent=app.OBC1)
        structure = openmm.load_topology(ommparm.topology, system)
        self.assertEqual(len(parm.atoms), len(structure.atoms))
        self.assertEqual(len(parm.residues), len(structure.residues))
        self.assertEqual(len(parm.bonds), len(structure.bonds))
        warnings.filterwarnings('always', category=exceptions.OpenMMWarning)

class TestWriteParameters(FileIOTestCase):

    @unittest.skipIf(os.getenv('AMBERHOME') is None, 'Cannot test w/out Amber')
    def testWriteXMLParameters(self):
        """ Test writing XML parameters loaded from Amber files """
        leaprc = StringIO("""\
logFile leap.log
#
# ----- leaprc for loading the ff14SB force field
# ----- NOTE: this is designed for PDB format 3!
#    Uses frcmod.ff14SB for proteins; ff99bsc0 for DNA; ff99bsc0_chiOL3 for RNA
#
#	load atom type hybridizations
#
addAtomTypes {
	{ "H"   "H" "sp3" }
	{ "HO"  "H" "sp3" }
	{ "HS"  "H" "sp3" }
	{ "H1"  "H" "sp3" }
	{ "H2"  "H" "sp3" }
	{ "H3"  "H" "sp3" }
	{ "H4"  "H" "sp3" }
	{ "H5"  "H" "sp3" }
	{ "HW"  "H" "sp3" }
	{ "HC"  "H" "sp3" }
	{ "HA"  "H" "sp3" }
	{ "HP"  "H" "sp3" }
	{ "HZ"  "H" "sp3" }
	{ "OH"  "O" "sp3" }
	{ "OS"  "O" "sp3" }
	{ "O"   "O" "sp2" }
	{ "O2"  "O" "sp2" }
	{ "OP"  "O" "sp2" }
	{ "OW"  "O" "sp3" }
	{ "CT"  "C" "sp3" }
	{ "CX"  "C" "sp3" }
	{ "C8"  "C" "sp3" }
	{ "2C"  "C" "sp3" }
	{ "3C"  "C" "sp3" }
	{ "CH"  "C" "sp3" }
	{ "CS"  "C" "sp2" }
	{ "C"   "C" "sp2" }
	{ "CO"   "C" "sp2" }
	{ "C*"  "C" "sp2" }
	{ "CA"  "C" "sp2" }
	{ "CB"  "C" "sp2" }
	{ "CC"  "C" "sp2" }
	{ "CN"  "C" "sp2" }
	{ "CM"  "C" "sp2" }
	{ "CK"  "C" "sp2" }
	{ "CQ"  "C" "sp2" }
	{ "CD"  "C" "sp2" }
	{ "C5"  "C" "sp2" }
	{ "C4"  "C" "sp2" }
	{ "CP"  "C" "sp2" }
	{ "CI"  "C" "sp3" }
	{ "CJ"  "C" "sp2" }
	{ "CW"  "C" "sp2" }
	{ "CV"  "C" "sp2" }
	{ "CR"  "C" "sp2" }
	{ "CA"  "C" "sp2" }
	{ "CY"  "C" "sp2" }
	{ "C0"  "Ca" "sp3" }
	{ "MG"  "Mg" "sp3" }
	{ "N"   "N" "sp2" }
	{ "NA"  "N" "sp2" }
	{ "N2"  "N" "sp2" }
	{ "N*"  "N" "sp2" }
	{ "NP"  "N" "sp2" }
	{ "NQ"  "N" "sp2" }
	{ "NB"  "N" "sp2" }
	{ "NC"  "N" "sp2" }
	{ "NT"  "N" "sp3" }
	{ "NY"  "N" "sp2" }
	{ "N3"  "N" "sp3" }
	{ "S"   "S" "sp3" }
	{ "SH"  "S" "sp3" }
	{ "P"   "P" "sp3" }
	{ "LP"  ""  "sp3" }
	{ "EP"  ""  "sp3" }
	{ "F"   "F" "sp3" }
	{ "Cl"  "Cl" "sp3" }
	{ "Br"  "Br" "sp3" }
	{ "I"   "I"  "sp3" }
	{ "F-"   "F" "sp3" }
	{ "Cl-"  "Cl" "sp3" }
	{ "Br-"  "Br" "sp3" }
	{ "I-"   "I"  "sp3" }
	{ "Li+"  "Li"  "sp3" }
	{ "Na+"  "Na"  "sp3" }
	{ "K+"  "K"  "sp3" }
	{ "Rb+"  "Rb"  "sp3" }
	{ "Cs+"  "Cs"  "sp3" }
	{ "Mg+"  "Mg"  "sp3" }
# glycam
	{ "OG"  "O" "sp3" }
	{ "OL"  "O" "sp3" }
	{ "AC"  "C" "sp3" }
	{ "EC"  "C" "sp3" }
}
#
#	Load the main parameter set.
#
parm10 = loadamberparams parm10.dat
frcmod14SB = loadamberparams frcmod.ff14SB
#
#	Load main chain and terminating amino acid libraries, nucleic acids
#
loadOff amino12.lib
loadOff aminoct12.lib
loadOff aminont12.lib
loadOff nucleic12.lib
#
#       Load water and ions
# 
#loadOff atomic_ions.lib
#loadOff solvents.lib
#HOH = TP3
#WAT = TP3

#
#	Define the PDB name map for the amino acids and nucleic acids
#
addPdbResMap {
  { 0 "HYP" "NHYP" } { 1 "HYP" "CHYP" }
  { 0 "ALA" "NALA" } { 1 "ALA" "CALA" }
  { 0 "ARG" "NARG" } { 1 "ARG" "CARG" }
  { 0 "ASN" "NASN" } { 1 "ASN" "CASN" }
  { 0 "ASP" "NASP" } { 1 "ASP" "CASP" }
  { 0 "CYS" "NCYS" } { 1 "CYS" "CCYS" }
  { 0 "CYX" "NCYX" } { 1 "CYX" "CCYX" }
  { 0 "GLN" "NGLN" } { 1 "GLN" "CGLN" }
  { 0 "GLU" "NGLU" } { 1 "GLU" "CGLU" }
  { 0 "GLY" "NGLY" } { 1 "GLY" "CGLY" }
  { 0 "HID" "NHID" } { 1 "HID" "CHID" }
  { 0 "HIE" "NHIE" } { 1 "HIE" "CHIE" }
  { 0 "HIP" "NHIP" } { 1 "HIP" "CHIP" }
  { 0 "ILE" "NILE" } { 1 "ILE" "CILE" }
  { 0 "LEU" "NLEU" } { 1 "LEU" "CLEU" }
  { 0 "LYS" "NLYS" } { 1 "LYS" "CLYS" }
  { 0 "MET" "NMET" } { 1 "MET" "CMET" }
  { 0 "PHE" "NPHE" } { 1 "PHE" "CPHE" }
  { 0 "PRO" "NPRO" } { 1 "PRO" "CPRO" }
  { 0 "SER" "NSER" } { 1 "SER" "CSER" }
  { 0 "THR" "NTHR" } { 1 "THR" "CTHR" }
  { 0 "TRP" "NTRP" } { 1 "TRP" "CTRP" }
  { 0 "TYR" "NTYR" } { 1 "TYR" "CTYR" }
  { 0 "VAL" "NVAL" } { 1 "VAL" "CVAL" }
  { 0 "HIS" "NHIS" } { 1 "HIS" "CHIS" }
  { 0 "G" "G5"  } { 1 "G" "G3"  } 
  { 0 "A" "A5"  } { 1 "A" "A3"  } 
  { 0 "C" "C5"  } { 1 "C" "C3"  } 
  { 0 "U" "U5"  } { 1 "U" "U3"  } 
  { 0 "DG" "DG5"  } { 1 "DG" "DG3"  }  
  { 0 "DA" "DA5"  } { 1 "DA" "DA3"  }  
  { 0 "DC" "DC5"  } { 1 "DC" "DC3"  }  
  { 0 "DT" "DT5"  } { 1 "DT" "DT3"  }  
#  some old Amber residue names for RNA:
  { 0  "RA5" "A5" } { 1 "RA3" "A3"} {"RA" "A" }
  { 0  "RC5" "C5" } { 1 "RC3" "C3"} {"RC" "C" }
  { 0  "RG5" "G5" } { 1 "RG3" "G3"} {"RG" "G" }
  { 0  "RU5" "U5" } { 1 "RU3" "U3"} {"RU" "U" }
#  some really old Amber residue names, assuming DNA:
  { 0 "GUA" "DG5"  } { 1 "GUA" "DG3"  } { "GUA" "DG" }
  { 0 "ADE" "DA5"  } { 1 "ADE" "DA3"  } { "ADE" "DA" }
  { 0 "CYT" "DC5"  } { 1 "CYT" "DC3"  } { "CYT" "DC" }
  { 0 "THY" "DT5"  } { 1 "THY" "DT3"  } { "THY" "DT" }
#  uncomment out the following if you have this old style RNA files:
# { 0 "GUA" "G5"  } { 1 "GUA" "G3"  } { "GUA" "G" }
# { 0 "ADE" "A5"  } { 1 "ADE" "A3"  } { "ADE" "A" }
# { 0 "CYT" "C5"  } { 1 "CYT" "C3"  } { "CYT" "C" }
# { 0 "URA" "R5"  } { 1 "URA" "R3"  } { "URA" "R" }

}

#  try to be good about reading in really old atom names as well:
addPdbAtomMap {
  { "O5*" "O5'" }
  { "C5*" "C5'" }
  { "C4*" "C4'" }
  { "O4*" "O4'" }
  { "C3*" "C3'" }
  { "O3*" "O3'" }
  { "C2*" "C2'" }
  { "O2*" "O2'" }
  { "C1*" "C1'" }
  { "C5M" "C7"  }
  { "H1*" "H1'" }
  { "H2*1" "H2'" }
  { "H2*2" "H2''" }
  { "H2'1" "H2'" }
  { "H2'2" "H2''" }
  { "H3*" "H3'" }
  { "H4*" "H4'" }
  { "H5*1" "H5'" }
  { "H5*2" "H5''" }
  { "H5'1" "H5'" }
  { "H5'2" "H5''" }
  { "HO'2" "HO2'" }
  { "H5T"  "HO5'" }
  { "H3T"  "HO3'" }
  { "O1'" "O4'" }
  { "OA"  "OP1" }
  { "OB"  "OP2" }
  { "O1P" "OP1" }
  { "O2P" "OP2" }
}

#
# assume that most often proteins use HIE
#
NHIS = NHIE
HIS = HIE
CHIS = CHIE
""")
        params = openmm.OpenMMParameterSet.from_parameterset(
                amber.AmberParameterSet.from_leaprc(leaprc)
        )
        params.write(get_fn('amber_conv.xml', written=True),
                     provenance=dict(OriginalFile='leaprc.ff14SB',
                     Reference='Maier & Simmerling')
        )


    @unittest.skipIf(os.getenv('AMBERHOME') is None, 'Cannot test w/out Amber')
    def testWriteXMLParametersGAFF(self):
        """ Test writing XML parameters loaded from Amber GAFF parameter files """
        leaprc = StringIO("""\
parm10 = loadamberparams gaff.dat
""")
        params = openmm.OpenMMParameterSet.from_parameterset(
                amber.AmberParameterSet.from_leaprc(leaprc)
        )
        citations = """\
Wang, J., Wang, W., Kollman P. A.; Case, D. A. "Automatic atom type and bond type perception in molecular mechanical calculations". Journal of Molecular Graphics and Modelling , 25, 2006, 247260.
Wang, J., Wolf, R. M.; Caldwell, J. W.;Kollman, P. A.; Case, D. A. "Development and testing of a general AMBER force field". Journal of Computational Chemistry, 25, 2004, 1157-1174.
"""
        params.write(get_fn('gaff.xml', written=True),
                     provenance=dict(OriginalFile='gaff.dat',
                     Reference=citations)
        )

    def testWriteXMLParametersCharmm(self):
        """ Test writing XML parameter files from Charmm parameter files"""

        params = openmm.OpenMMParameterSet.from_parameterset(
                charmm.CharmmParameterSet(get_fn('par_all36_prot.prm'),
                                          get_fn('top_all36_prot.rtf'))
        )
        params.write(get_fn('charmm_conv.xml', written=True),
                     provenance=dict(
                         OriginalFile='par_all36_prot.prm & top_all36_prot.rtf',
                         Reference='MacKerrell'
                     )
        )

