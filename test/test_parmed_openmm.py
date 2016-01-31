""" Tests some OpenMM-specific functionality """
from __future__ import print_function, division, absolute_import

import math
import numpy as np
import os
import parmed as pmd
from parmed import openmm, load_file, exceptions, ExtraPoint, unit as u
import unittest
from utils import get_fn, mm, app, has_openmm, FileIOTestCase
import warnings

@unittest.skipUnless(has_openmm, "Cannot test without OpenMM")
class TestOpenMM(FileIOTestCase):

    def setUp(self):
        # Take one of the distributed OpenMM FF XML files as a test
        self.ffxml = os.path.join(os.path.split(app.__file__)[0], 'data',
                                  'amber99sbildn.xml')
        warnings.filterwarnings('error', category=exceptions.OpenMMWarning)
        super(TestOpenMM, self).setUp()

    def tearDown(self):
        warnings.filterwarnings('always', category=exceptions.OpenMMWarning)
        super(TestOpenMM, self).tearDown()

    def test_format_id(self):
        """ Tests automatic format determination of OpenMM XML files """
        self.assertTrue(openmm.XmlFile.id_format(get_fn('system_974wat.xml')))
        self.assertTrue(openmm.XmlFile.id_format(get_fn('state_974wat.xml')))
        self.assertTrue(openmm.XmlFile.id_format(get_fn('integrator.xml')))
        self.assertTrue(openmm.XmlFile.id_format(self.ffxml))
        fn = get_fn('test.xml', written=True)
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
        import warnings
        ommparm = app.AmberPrmtopFile(get_fn('complex.prmtop'))
        parm = load_file(get_fn('complex.prmtop'))
        system = ommparm.createSystem(implicitSolvent=app.OBC1)
        structure = openmm.load_topology(ommparm.topology, system)
        self.assertEqual(len(parm.atoms), len(structure.atoms))
        self.assertEqual(len(parm.residues), len(structure.residues))
        self.assertEqual(len(parm.bonds), len(structure.bonds))

    def test_load_topology_extra_bonds(self):
        """ Test loading extra bonds not in Topology """
        parm = load_file(get_fn('ash.parm7'))
        system = parm.createSystem()
        for force in system.getForces():
            if isinstance(force, mm.HarmonicBondForce):
                force.addBond(0, parm[-1].idx, 1, 500)
        self.assertRaises(exceptions.OpenMMWarning, lambda:
                openmm.load_topology(parm.topology, system))
        warnings.filterwarnings('ignore', category=exceptions.OpenMMWarning)
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
        self.assertRaises(exceptions.OpenMMWarning, lambda:
                openmm.load_topology(parm.topology, system))

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
