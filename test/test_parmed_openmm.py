""" Tests some OpenMM-specific functionality """
from __future__ import print_function, division, absolute_import

from collections import OrderedDict
import math
import os
import unittest
import warnings
from unittest import skipIf

import numpy as np

import parmed as pmd
from parmed.utils.six import StringIO
from io import TextIOWrapper
from parmed.charmm import (CharmmPsfFile, CharmmCrdFile, CharmmRstFile,
                           CharmmParameterSet)
from parmed import openmm, load_file, exceptions, ExtraPoint, unit as u
from utils import (get_fn, mm, app, has_openmm, has_networkx, has_lxml,
                   FileIOTestCase, CPU, TestCaseRelative, EnergyTestCase)
from parmed.exceptions import ParameterWarning

@unittest.skipUnless(has_openmm, "Cannot test without OpenMM")
class TestOpenMM(FileIOTestCase, EnergyTestCase):

    def setUp(self):
        super(TestOpenMM, self).setUp()
        # Take one of the distributed OpenMM FF XML files as a test
        self.ffxml = os.path.join(os.path.split(app.__file__)[0], 'data',
                                  'amber99sbildn.xml')
        super(TestOpenMM, self).setUp()

    def test_format_id(self):
        """ Tests automatic format determination of OpenMM XML files """
        self.assertTrue(openmm.XmlFile.id_format(get_fn('system_974wat.xml')))
        self.assertTrue(openmm.XmlFile.id_format(get_fn('state_974wat.xml')))
        self.assertTrue(openmm.XmlFile.id_format(get_fn('integrator.xml')))
        self.assertTrue(openmm.XmlFile.id_format(self.ffxml))
        fn = self.get_fn('test.xml', written=True)
        with open(fn, 'w') as f:
            f.write('<NoRecognizedObject>\n <SomeElement attr="yes" '
                    '/>\n</NoRecognizedObject>\n')
        self.assertFalse(openmm.XmlFile.id_format(fn))

    def test_deserialize_system(self):
        """ Tests automatic deserialization of a System XML file """
        system = openmm.XmlFile.parse(get_fn('system_974wat.xml'))
        self.assertIsInstance(system, mm.System)
        self.assertEqual(system.getNumParticles(), 6638)

    def test_deserialize_state(self):
        """ Tests automatic deserialization of a State XML file """
        state = openmm.XmlFile.parse(get_fn('state_974wat.xml'))
        self.assertEqual(state.coordinates.shape, (1, 6638, 3))
        self.assertEqual(state.velocities.shape, (1, 6638, 3))
        self.assertEqual(state.forces.shape, (1, 6638, 3))
        np.testing.assert_allclose(state.box, [35.05, 40.5, 42.37, 90, 90, 90])
        self.assertAlmostEqual(state.time, 20000.000003615783)
        self.assertIs(state.energy, None)
        # Test with open file handle
        with open(get_fn('state_974wat.xml'), 'r') as f:
            state = openmm.XmlFile.parse(f)
        self.assertEqual(state.coordinates.shape, (1, 6638, 3))
        self.assertEqual(state.velocities.shape, (1, 6638, 3))
        self.assertEqual(state.forces.shape, (1, 6638, 3))
        np.testing.assert_allclose(state.box, [35.05, 40.5, 42.37, 90, 90, 90])
        self.assertAlmostEqual(state.time, 20000.000003615783)
        self.assertIs(state.energy, None)

    def test_deserialize_integrator(self):
        """ Tests automatic deserialization of an Integrator XML file """
        integrator = openmm.XmlFile.parse(get_fn('integrator.xml'))
        self.assertIsInstance(integrator, mm.Integrator)
        self.assertIsInstance(integrator, mm.LangevinIntegrator)

    def test_deserialize_force_field(self):
        """ Tests automatic deserialization of an OpenMM ForceField XML file """
        ff = openmm.XmlFile.parse(self.ffxml)
        self.assertIsInstance(ff, app.ForceField)

    def test_load_topology(self):
        """ Tests loading an OpenMM Topology and System instance """
        ommparm = app.AmberPrmtopFile(get_fn('complex.prmtop'))
        parm = load_file(get_fn('complex.prmtop'))
        system = ommparm.createSystem(implicitSolvent=app.OBC1)
        structure = openmm.load_topology(ommparm.topology, system)
        self.assertEqual(len(parm.atoms), len(structure.atoms))
        self.assertEqual(len(parm.residues), len(structure.residues))
        self.assertEqual(len(parm.bonds), len(structure.bonds))

    def test_load_topology_uncondensed_atom_types(self):
        """ Tests condense_atom_types arg """
        ommparm = app.AmberPrmtopFile(get_fn('complex.prmtop'))
        parm = load_file(get_fn('complex.prmtop'))
        system = ommparm.createSystem(implicitSolvent=app.OBC1)
        structure_condensed = openmm.load_topology(ommparm.topology, system, condense_atom_types=True)
        structure_uncondensed = openmm.load_topology(ommparm.topology, system, condense_atom_types=False)

        num_unique_atom_types_condensed = len(set([a.atom_type for a in structure_condensed]))
        num_unique_atom_types_uncondensed = len(set([a.atom_type for a in structure_uncondensed]))

        assert num_unique_atom_types_uncondensed > num_unique_atom_types_condensed

        # This may change if this or a similar flag does a partial condensing, see PR #1087
        assert num_unique_atom_types_uncondensed == len(structure_uncondensed.atoms)

    def test_load_topology_use_atom_id_as_typename(self):
        """ Tests loading an OpenMM Topology and using Atom.id to name types """
        ommparm = load_file(get_fn('ash.parm7'), get_fn('ash.rst7'))
        parm = load_file(get_fn('ash.parm7'), get_fn('ash.rst7'))
        system = ommparm.createSystem(implicitSolvent=app.OBC1)

        for pmd_atom, omm_atom in zip(parm.atoms, ommparm.topology.atoms()):
            omm_atom.id = pmd_atom.type
        structure = openmm.load_topology(ommparm.topology, system,
                                         xyz=parm.positions)

        self.assertEqual(len(parm.atoms), len(structure.atoms))
        self.assertEqual([a.type for a in parm.atoms],
                         [a.type for a in structure.atoms])
        self.assertEqual(len(parm.residues), len(structure.residues))
        self.assertEqual(len(parm.bonds), len(structure.bonds))

        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(parm.positions)
        con2.setPositions(structure.positions)

        self.check_energies(parm, con1, structure, con2)

    def test_load_topology_extra_bonds(self):
        """ Test loading extra bonds not in Topology """
        parm = load_file(get_fn('ash.parm7'))
        system = parm.createSystem()
        for force in system.getForces():
            if isinstance(force, mm.HarmonicBondForce):
                force.addBond(0, parm[-1].idx, 1, 500)
        with self.assertWarns(exceptions.OpenMMWarning):
            top = openmm.load_topology(parm.topology, system)
        self.assertIn(top[-1], top[0].bond_partners)
        self.assertEqual(len(top.bonds), len(parm.bonds)+1)

    def test_load_topology_eps(self):
        """ Tests loading an OpenMM Topology with Extra Points """
        parm = load_file(get_fn('tip4p.parm7'), get_fn('tip4p.rst7'))
        struct = openmm.load_topology(parm.topology, xyz=get_fn('tip4p.rst7'))
        has_eps = False
        self.assertEqual(len(parm.atoms), len(struct.atoms))
        for a1, a2 in zip(parm.atoms, struct.atoms):
            self.assertIs(type(a1), type(a2))
            has_eps = isinstance(a1, ExtraPoint) or has_eps
        self.assertTrue(has_eps)
        np.testing.assert_equal(parm.coordinates, struct.coordinates)
        np.testing.assert_equal(parm.box, struct.box)
        # Now try passing in coordinates
        struct = openmm.load_topology(parm.topology, xyz=parm.coordinates)
        np.testing.assert_equal(parm.coordinates, struct.coordinates)
        np.testing.assert_allclose(parm.box, struct.box) # taken from Topology
        struct = openmm.load_topology(parm.topology, xyz=parm.coordinates)
        np.testing.assert_equal(parm.coordinates, struct.coordinates)
        struct = openmm.load_topology(parm.topology, xyz=parm.coordinates,
                                      box=[10, 10, 10, 90, 90, 90])
        np.testing.assert_equal(parm.coordinates, struct.coordinates)
        np.testing.assert_equal(struct.box, [10, 10, 10, 90, 90, 90])

    def test_load_topology_error(self):
        """ Test error handling in load_topology """
        parm = load_file(get_fn('ash.parm7'), get_fn('ash.rst7'))
        parm2 = load_file(get_fn('solv2.parm7'), get_fn('solv2.rst7'))
        self.assertRaises(TypeError, lambda:
                openmm.load_topology(parm.topology, system=get_fn('integrator.xml'))
        )
        self.assertRaises(TypeError, lambda:
                openmm.load_topology(parm2.topology, system=parm.createSystem())
        )
        system = parm.createSystem()
        system.addForce(mm.CustomTorsionForce('theta**2'))
        with self.assertWarns(exceptions.OpenMMWarning):
            openmm.load_topology(parm.topology, system)

    def test_box_from_system(self):
        """ Tests loading box from System """
        parm = load_file(get_fn('solv2.parm7'), get_fn('solv2.rst7'))
        system = parm.createSystem(nonbondedMethod=app.PME,
                                   nonbondedCutoff=8*u.angstroms)
        top = openmm.load_topology(parm.topology, system)
        np.testing.assert_allclose(parm.box, top.box)

@unittest.skipUnless(has_openmm, 'Cannot test without OpenMM')
class TestSystemCreation(unittest.TestCase):
    """ Test various aspect of System creation """

    def test_two_particle_vsite(self):
        """ Tests assignment of 2-particle virtual site """
        struct = pmd.Structure()
        struct.add_atom(pmd.Atom(name='C', atomic_number=6), 'RES', 1)
        struct.add_atom(pmd.Atom(name='C', atomic_number=6), 'RES', 1)
        struct.add_atom(pmd.ExtraPoint(name='EP', atomic_number=0), 'RES', 1)
        struct.bond_types.append(pmd.BondType(10, 1.0))
        struct.bond_types.append(pmd.BondType(10, 0.5))
        struct.bond_types.claim()
        struct.bonds.append(pmd.Bond(struct[0], struct[1],
                                     type=struct.bond_types[0])
        )
        struct.bonds.append(pmd.Bond(struct[1], struct[2],
                                     type=struct.bond_types[1])
        )
        # This should be a two-particle virtual site
        struct.coordinates = [[0, 0, 0], [0, 0, 1], [0, 0, 1.5]]
        system = mm.System()
        system.addParticle(struct[0].mass)
        system.addParticle(struct[1].mass)
        system.addParticle(struct[2].mass)
        struct.omm_set_virtual_sites(system)
        # Make sure the third atom is a virtual site
        self.assertTrue(system.isVirtualSite(2))
        self.assertIsInstance(system.getVirtualSite(2), mm.TwoParticleAverageSite)

    def test_ep_exceptions(self):
        """ Test Nonbonded exception handling with virtual sites """
        # Analyze the exception parameters for bonding pattern
        #
        # E1 -- A1 -- A2 -- A3 -- A4 -- A5 -- E5
        #             |     |     |
        #             E2    E3    E4
        struct = pmd.Structure()
        ep1 = ExtraPoint(name='E1', type='EP', atomic_number=0, weights=[1, 2])
        ep2 = ExtraPoint(name='E2', type='EP', atomic_number=0)
        ep3 = ExtraPoint(name='E3', type='EP', atomic_number=0)
        ep4 = ExtraPoint(name='E4', type='EP', atomic_number=0)
        ep5 = ExtraPoint(name='E5', type='EP', atomic_number=0)
        self.assertIs(ep1.parent, None)
        self.assertEqual(ep1.bond_partners, [])
        self.assertEqual(ep1.angle_partners, [])
        self.assertEqual(ep1.dihedral_partners, [])
        self.assertEqual(ep1.tortor_partners, [])
        self.assertEqual(ep1.exclusion_partners, [])
        a1 = pmd.Atom(name='A1', type='AX', charge=0.1, atomic_number=6)
        a2 = pmd.Atom(name='A2', type='AY', charge=0.1, atomic_number=6)
        a3 = pmd.Atom(name='A3', type='AZ', charge=0.1, atomic_number=7)
        a4 = pmd.Atom(name='A4', type='AX', charge=0.1, atomic_number=6)
        a5 = pmd.Atom(name='A5', type='AY', charge=0.1, atomic_number=6)
        a1.rmin = a2.rmin = a3.rmin = a4.rmin = a5.rmin = 0.5
        a1.epsilon = a2.epsilon = a3.epsilon = a4.epsilon = a5.epsilon = 1.0
        bond_type = pmd.BondType(10.0, 1.0)
        bond_type2 = pmd.BondType(10.0, 2.0)
        bond_type3 = pmd.BondType(10.0, 0.5)
        bond_type4 = pmd.BondType(10.0, math.sqrt(2))
        angle_type = pmd.AngleType(10.0, 90)
        dihedral_type = pmd.DihedralType(10.0, 2, 0)
        struct.add_atom(a1, 'RES', 1)
        struct.add_atom(a2, 'RES', 1)
        struct.add_atom(a3, 'RES', 1)
        struct.add_atom(a4, 'RES', 1)
        struct.add_atom(a5, 'RES', 1)
        struct.add_atom(ep1, 'RES', 1)
        struct.add_atom(ep2, 'RES', 1)
        struct.add_atom(ep3, 'RES', 1)
        struct.add_atom(ep4, 'RES', 1)
        struct.add_atom(ep5, 'RES', 1)
        struct.bonds.extend(
                [pmd.Bond(a1, ep1, type=bond_type),
                 pmd.Bond(ep2, a2, type=bond_type),
                 pmd.Bond(a3, ep3, type=bond_type3),
                 pmd.Bond(a4, ep4, type=bond_type)]
        )
        struct.bonds.extend(
                [pmd.Bond(a1, a2, type=bond_type),
                 pmd.Bond(a4, a3, type=bond_type4),
                 pmd.Bond(a3, a2, type=bond_type4),
                 pmd.Bond(a4, a5, type=bond_type2),
                 pmd.Bond(a5, ep5, type=bond_type)]
        )
        struct.angles.extend(
                [pmd.Angle(a1, a2, a3, type=angle_type),
                 pmd.Angle(a2, a3, a4, type=angle_type),
                 pmd.Angle(a3, a4, a5, type=angle_type)]
        )
        struct.dihedrals.extend(
                [pmd.Dihedral(a1, a2, a3, a4, type=dihedral_type),
                 pmd.Dihedral(a2, a3, a4, a5, type=dihedral_type)]
        )
        struct.bond_types.extend([bond_type, bond_type3, bond_type2, bond_type4])
        struct.angle_types.append(angle_type)
        struct.dihedral_types.append(dihedral_type)
        # Test exclusions now
        a1.exclude(a5)
        system = struct.createSystem()

@unittest.skipUnless(has_openmm, 'Cannot test without OpenMM')
@unittest.skipUnless(has_lxml, 'Cannot test without lxml')
@unittest.skipUnless(has_networkx, 'Cannot test without networkx')
@unittest.skipUnless(os.getenv('AMBERHOME'), 'Cannot test without AMBERHOME')
class TestWriteAmberParameters(FileIOTestCase):

    def test_write_xml_parameters(self):
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
                pmd.amber.AmberParameterSet.from_leaprc(leaprc)
        )
        ffxml_filename = self.get_fn('amber_conv.xml', written=True)
        params.write(ffxml_filename,
                     provenance=dict(OriginalFile='leaprc.ff14SB',
                     Reference=['Maier and Simmerling', 'Simmerling and Maier'],
                     Source=dict(Source='leaprc.ff14SB',
                     sourcePackage='AmberTools', sourcePackageVersion='15'))
        )
        forcefield = app.ForceField(ffxml_filename)
        # Make sure the forcefield can handle proteins with disulfide bonds
        pdbfile = app.PDBFile(get_fn('3qyt_fix.pdb'))
        forcefield.createSystem(pdbfile.topology, nonbondedMethod=app.NoCutoff)


    def test_write_xml_parameters_gaff(self):
        """ Test writing XML parameters loaded from Amber GAFF parameter files """
        leaprc = StringIO("""\
parm10 = loadamberparams gaff.dat
""")
        params = openmm.OpenMMParameterSet.from_parameterset(
                pmd.amber.AmberParameterSet.from_leaprc(leaprc)
        )
        citations = """\
Wang, J., Wang, W., Kollman P. A.; Case, D. A. "Automatic atom type and bond type perception in molecular mechanical calculations". Journal of Molecular Graphics and Modelling , 25, 2006, 247260.
Wang, J., Wolf, R. M.; Caldwell, J. W.;Kollman, P. A.; Case, D. A. "Development and testing of a general AMBER force field". Journal of Computational Chemistry, 25, 2004, 1157-1174.
"""
        ffxml_filename = self.get_fn('gaff.xml', written=True)
        params.write(ffxml_filename,
                     provenance=dict(OriginalFile='gaff.dat',
                     Reference=citations)
        )
        forcefield = app.ForceField(ffxml_filename)

    def test_write_xml_parameters_antechamber(self):
        """ Test writing XML residue definition from Antechamber mol2 """
        leaprc = "molecule = loadmol2 %s\n" % get_fn('molecule.mol2')
        leaprc += "loadamberparams %s\n" % get_fn('molecule.frcmod')
        leaprc = StringIO(leaprc)
        params = openmm.OpenMMParameterSet.from_parameterset(
                pmd.amber.AmberParameterSet.from_leaprc(leaprc),
                remediate_residues=False
        )
        ffxml_filename = self.get_fn('residue.xml', written=True)
        params.write(ffxml_filename)
        try:
            forcefield = app.ForceField(ffxml_filename)
        except KeyError:
            # A KeyError is expected
            pass
        else:
            assert False, "app.ForceField() should fail with a KeyError when residue templates without parameters are not removed, but it did not."

    def test_write_xml_parameters_amber_write_unused(self):
        """Test the write_unused argument in writing XML files"""
        params = openmm.OpenMMParameterSet.from_parameterset(
                pmd.amber.AmberParameterSet(get_fn('amino12.lib'),
                os.path.join(get_fn('parm'), 'parm10.dat'),
                os.path.join(get_fn('parm'), 'frcmod.ff14SB'))
        )
        ffxml = StringIO()
        params.write(ffxml)
        ffxml.seek(0)
        self.assertEqual(len(ffxml.readlines()), 2179)
        ffxml = StringIO()
        params.write(ffxml, write_unused=False)
        ffxml.seek(0)
        self.assertEqual(len(ffxml.readlines()), 1647)
        ffxml.seek(0)
        forcefield = app.ForceField(ffxml)

        params = openmm.OpenMMParameterSet.from_parameterset(
                  pmd.amber.AmberParameterSet(get_fn('atomic_ions.lib'),
                  os.path.join(get_fn('parm'), 'frcmod.ionsjc_tip3p'))
        )
        ffxml = StringIO()
        params.write(ffxml, write_unused=True)
        ffxml.seek(0)
        self.assertEqual(len(ffxml.readlines()), 222)
        ffxml = StringIO()
        params.write(ffxml, write_unused=False)
        ffxml.seek(0)
        self.assertEqual(len(ffxml.readlines()), 57)
        ffxml.seek(0)

        forcefield = app.ForceField(ffxml)

        # Load TIP3P water box to ensure there are no duplicate ion parameters
        pdbfile = app.PDBFile(get_fn('ionsjc.pdb'))
        system = forcefield.createSystem(pdbfile.topology)

    def test_write_xml_small_amber(self):
        """ Test writing small XML modifications """
        params = openmm.OpenMMParameterSet.from_parameterset(
                load_file(os.path.join(get_fn('parm'), 'frcmod.constph'))
        )
        ffxml_filename = self.get_fn('test.xml', written=True)
        params.write(ffxml_filename)
        forcefield = app.ForceField(ffxml_filename)

    def test_not_write_residues_with_same_templhash(self):
        """Test that no identical residues are written to XML, using the templhasher function."""
        # TODO add testing for multiatomic residues when support for those added
        params = openmm.OpenMMParameterSet.from_parameterset(
                  pmd.amber.AmberParameterSet(get_fn('atomic_ions.lib'),
                  os.path.join(get_fn('parm'), 'frcmod.ionsjc_tip3p'))
        )
        new_residues = OrderedDict()
        for name in ('K', 'K+', 'NA', 'Na+', 'CL', 'Cl-'):
            new_residues[name] = params.residues[name]
        params.residues = new_residues
        ffxml = StringIO()
        params.write(ffxml)
        # TODO: Overhaul tests with lxml (xpath?) queries to eliminate dependence on file order
        ffxml.seek(0)
        nlines = len(ffxml.readlines())
        self.assertEqual(nlines, 39, 'File contents:\n{}\nActual length: {} lines.'.format(ffxml.getvalue(), nlines))

    def test_override_level(self):
        """Test correct support for the override_level attribute of ResidueTemplates and correct writing to XML tag"""
        params = openmm.OpenMMParameterSet.from_parameterset(
                  pmd.amber.AmberParameterSet(get_fn('atomic_ions.lib'),
                  os.path.join(get_fn('parm'), 'frcmod.ionsjc_tip3p'))
        )
        new_residues = OrderedDict()
        new_residues['K'] = params.residues['K']
        new_residues['NA'] = params.residues['NA']
        new_residues['K'].override_level = 1
        params.residues = new_residues
        ffxml = StringIO()
        params.write(ffxml)
        ffxml.seek(0)
        # TODO: Overhaul tests with lxml (xpath?) queries to eliminate dependence on file order
        output_lines = ffxml.readlines()
        control_line1 = '  <Residue name="K" override="1">\n'
        control_line2 = '  <Residue name="NA">\n'
        self.assertEqual(output_lines[16].strip(), control_line1.strip(), 'File contents:\n{}'.format(ffxml.getvalue()))
        self.assertEqual(output_lines[19].strip(), control_line2.strip(), 'File contents:\n{}'.format(ffxml.getvalue()))

@unittest.skipUnless(has_openmm, 'Cannot test without OpenMM')
@unittest.skipUnless(has_lxml, 'Cannot test without lxml')
@unittest.skipUnless(has_networkx, 'Cannot test without networkx')
class TestWriteCHARMMParameters(FileIOTestCase):

    def test_write_xml_parameters_charmm_multisite_waters(self):
        """ Test writing XML parameter files from Charmm multisite water parameter files and reading them back into OpenMM ForceField """

        params = openmm.OpenMMParameterSet.from_parameterset(
                pmd.charmm.CharmmParameterSet(get_fn('toppar_water_ions_tip5p.str'))
        )
        ffxml_filename = self.get_fn('charmm_conv.xml', written=True)
        params.write(ffxml_filename,
                     provenance=dict(
                         OriginalFile='toppar_water_ions_tip5p.str',
                         Reference='MacKerrell'
                     )
        )
        forcefield = app.ForceField(ffxml_filename)

        # Check that water has the right number of bonds
        assert len(params.residues['TIP5'].bonds) == 2, "TIP5P should only have two bonds, but instead has {}".format(params.residues['TIP5'].bonds)

        # Parameterize water box
        pdbfile = app.PDBFile(get_fn('waterbox.pdb'))
        modeller = app.Modeller(pdbfile.topology, pdbfile.positions)
        modeller.addExtraParticles(forcefield)
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff)

        # Parameterize water box
        pdbfile = app.PDBFile(get_fn('waterbox-tip3p.pdb'))
        modeller = app.Modeller(pdbfile.topology, pdbfile.positions)
        modeller.addExtraParticles(forcefield)
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME)

    def test_write_xml_parameters_charmm(self):
        """ Test writing XML parameter files from Charmm parameter files and reading them back into OpenMM ForceField """

        params = openmm.OpenMMParameterSet.from_parameterset(
                pmd.charmm.CharmmParameterSet(get_fn('par_all36_prot.prm'),
                                              get_fn('top_all36_prot.rtf'),
                                              get_fn('toppar_water_ions.str')) # WARNING: contains duplicate water templates
        )

        # Check that patch ACE remains but duplicate patches ACED, ACP, ACPD
        assert 'ACE' in params.patches, "Patch ACE was expected to be retained, but is missing"
        for name in ['ACED', 'ACP', 'ACPD']:
            assert name not in params.patches, "Duplicate patch {} was expected to be pruned, but is present".format(name)

        del params.residues['TP3M'] # Delete to avoid duplicate water template topologies
        ffxml_filename = self.get_fn('charmm_conv.xml', written=True)
        params.write(ffxml_filename,
                     provenance=dict(
                         OriginalFile='par_all36_prot.prm, top_all36_prot.rtf',
                         Reference='MacKerrell'
                     )
        )
        forcefield = app.ForceField(ffxml_filename)

        # Check that water has the right number of bonds
        assert len(params.residues['TIP3'].bonds) == 2, "TIP3P should only have two bonds"

        # Parameterize alanine tripeptide in vacuum
        pdbfile = app.PDBFile(get_fn('ala_ala_ala.pdb'))
        system = forcefield.createSystem(pdbfile.topology, nonbondedMethod=app.NoCutoff)
        # Parameterize ACE-NME in water
        pdbfile = app.PDBFile(get_fn('2igd_924wat.pdb'))
        system = forcefield.createSystem(pdbfile.topology, nonbondedMethod=app.PME)

    def test_write_xml_parameters_methanol_ions_energy(self):
        """ Test writing XML parameter files from Charmm parameter files, reading them back into OpenMM ForceField, and computing energy of methanol and NaCl """

        params = openmm.OpenMMParameterSet.from_parameterset(
                pmd.charmm.CharmmParameterSet(get_fn('par_all36_cgenff.prm'),
                                              get_fn('top_all36_cgenff.rtf'),
                                              get_fn('toppar_water_ions.str')) # WARNING: contains duplicate water templates
        )
        del params.residues['TP3M'] # Delete to avoid duplicate water template topologies
        ffxml_filename = self.get_fn('charmm_conv.xml', written=True)
        params.write(ffxml_filename,
                     provenance=dict(
                         OriginalFile='par_all36_cgenff.prm, top_all36_cgenff.rtf, toppar_water_ions.str',
                         Reference='MacKerrell'
                     )
        )
        forcefield = app.ForceField(ffxml_filename)
        # Parameterize methanol and ions in vacuum
        pdbfile = app.PDBFile(get_fn('methanol_ions.pdb'))
        system = forcefield.createSystem(pdbfile.topology, nonbondedMethod=app.NoCutoff)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        context = mm.Context(system, integrator, CPU)
        context.setPositions(pdbfile.positions)
        ffxml_potential = context.getState(getEnergy=True).getPotentialEnergy() / u.kilocalories_per_mole
        del context, integrator
        # Compute energy via ParmEd reader
        psf = CharmmPsfFile(get_fn('methanol_ions.psf'))
        system = psf.createSystem(params)
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        context = mm.Context(system, integrator, CPU)
        context.setPositions(pdbfile.positions)
        parmed_potential = context.getState(getEnergy=True).getPotentialEnergy() / u.kilocalories_per_mole
        del context, integrator
        # Ensure potentials are almost equal
        self.assertAlmostEqual(ffxml_potential, parmed_potential)

    def test_duplicate_patch_removal(self):
        """ Test importing CHARMM36 protein forcefield correctly removes duplicate patches """

        charmm_params = pmd.charmm.CharmmParameterSet(get_fn('par_all36_prot.prm'),
                                                      get_fn('top_all36_prot.rtf'))

        openmm_params = openmm.OpenMMParameterSet.from_parameterset(charmm_params)

        assert 'ACE' in openmm_params.patches # appears first, should be retained
        assert 'ACED' not in openmm_params.patches # duplicate for OpenMM, should be removed
        assert 'ACP' not in openmm_params.patches # duplicate for OpenMM, should be removed
        assert 'ACPD' not in openmm_params.patches # duplicate for OpenMM, should be removed

    def test_ljforce_charmm(self):
        """ Test writing LennardJonesForce without NBFIX from Charmm parameter files and reading them back into OpenMM ForceField """

        charmm_params = pmd.charmm.CharmmParameterSet(get_fn('par_all36_prot.prm'),
                                                      get_fn('top_all36_prot.rtf'))

        openmm_params = openmm.OpenMMParameterSet.from_parameterset(charmm_params)

        #openmm_params.write(self.get_fn('charmm.xml', written=True),
        ffxml_filename = get_fn('charmm36.xml')
        openmm_params.write(ffxml_filename,
                            provenance=dict(
                                OriginalFile='par_all36_prot.prm & top_all36_prot.rtf',
                                Reference='MacKerrell'
                                ),
                            separate_ljforce=True
                            )
        forcefield = app.ForceField(ffxml_filename)
        # Parameterize alanine tripeptide in vacuum
        pdbfile = app.PDBFile(get_fn('ala_ala_ala.pdb'))
        system = forcefield.createSystem(pdbfile.topology, nonbondedMethod=app.NoCutoff)

    def test_explicit_improper(self):
        """ Test writing out the improper explicitly and reading it back into OpenMM ForceField """

        params = openmm.OpenMMParameterSet.from_parameterset(
                pmd.charmm.CharmmParameterSet(get_fn('par_all36_prot.prm'),
                                              get_fn('top_all36_prot.rtf'))
        )
        ffxml_filename = self.get_fn('charmm.xml', written=True)
        params.write(ffxml_filename,
                     provenance=dict(
                         OriginalFiles='par_all36_prot.prm & top_all36_prot.rtf',
                         Reference='MacKerrel'
                     ),
                     charmm_imp=True,
                     separate_ljforce=True,
                     )
        forcefield = app.ForceField(ffxml_filename)
