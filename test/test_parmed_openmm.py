""" Tests some OpenMM-specific functionality """
from __future__ import print_function, division, absolute_import

import numpy as np
from parmed import openmm, load_file, exceptions
import os
import unittest
from utils import get_fn, mm, app, has_openmm

@unittest.skipIf(not has_openmm, "Cannot test without OpenMM")
class TestOpenMM(unittest.TestCase):

    def setUp(self):
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
