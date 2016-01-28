""" Tests some OpenMM-specific functionality """
from __future__ import print_function, division, absolute_import

import numpy as np
import os
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
