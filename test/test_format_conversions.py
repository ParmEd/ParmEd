""" Test various topology format conversions """
from __future__ import print_function, division, absolute_import

import parmed as chem
from parmed.exceptions import GromacsTopologyWarning
from parmed.gromacs._gromacsfile import GromacsFile
from parmed import unit as u
import os
try:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    CPU = mm.Platform.getPlatformByName('CPU')
    has_openmm = True
except ImportError:
    has_openmm = False
import unittest
from utils import get_fn, get_saved_fn, diff_files, TestCaseRelative
import warnings

class TestCase(TestCaseRelative):
    def setUp(self):
        warnings.filterwarnings('default', category=GromacsTopologyWarning)
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

class TestAmberToGromacs(TestCase):
    """ Tests converting Amber prmtop files to Gromacs topologies """

    def testBenzeneCyclohexane(self):
        """ Test converting binary liquid from Amber prmtop to Gromacs top """
        parm = chem.load_file(get_fn('benzene_cyclohexane_10_500.prmtop'),
                              get_fn('benzene_cyclohexane_10_500.inpcrd'))
        top = chem.gromacs.GromacsTopologyFile.from_structure(parm)
        groname = get_fn('benzene_cyclohexane_10_500.gro', written=True)
        chem.gromacs.GromacsGroFile.write(parm, groname, precision=8)
        gro = chem.gromacs.GromacsGroFile.parse(groname)
        self.assertEqual(len(gro.atoms), len(parm.atoms))
        for a1, a2 in zip(gro.atoms, parm.atoms):
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.idx, a2.residue.idx)
            self.assertEqual(a1.name, a2.name)
            self.assertAlmostEqual(a1.xx, a2.xx)
            self.assertAlmostEqual(a1.xy, a2.xy)
            self.assertAlmostEqual(a1.xz, a2.xz)
        top.write(get_fn('benzene_cyclohexane_10_500.top', written=True))
        saved = GromacsFile(get_saved_fn('benzene_cyclohexane_10_500.top'))
        written = GromacsFile(get_fn('benzene_cyclohexane_10_500.top', written=True))
        self.assertTrue(diff_files(saved, written))

class TestGromacsToAmber(TestCase):
    """ Tests converting Gromacs top/gro files to Amber """

    def testSimple(self):
        """ Tests converting standard Gromacs system into Amber prmtop """
        top = chem.load_file(get_fn(os.path.join('03.AlaGlu', 'topol.top')))
        parm = chem.amber.AmberParm.from_structure(top)
        parm.write_parm(get_fn('ala_glu.parm7', written=True))
        parm = chem.load_file(get_fn('ala_glu.parm7', written=True))
        self.assertIsInstance(parm, chem.amber.AmberParm)
        self.assertEqual(len(top.atoms), len(parm.atoms))
        self.assertEqual(len(top.bonds), len(parm.bonds))
        self.assertEqual(len(top.angles), len(parm.angles))
        self.assertEqual(len(top.residues), len(parm.residues))
        for a1, a2 in zip(top.atoms, parm.atoms):
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.idx, a2.residue.idx)
            self.assertAlmostEqual(a1.mass, a2.mass)
            self.assertIs(type(a1), type(a2))
            self.assertAlmostEqual(a1.charge, a2.charge)
            self.assertEqual(set([a.idx for a in a1.bond_partners]),
                             set([a.idx for a in a2.bond_partners]))
            self.assertEqual(set([a.idx for a in a1.angle_partners]),
                             set([a.idx for a in a2.angle_partners]))
            self.assertEqual(set([a.idx for a in a1.dihedral_partners]),
                             set([a.idx for a in a2.dihedral_partners]))

    @unittest.skipIf(not has_openmm, "Cannot test without OpenMM")
    def testEnergySimple(self):
        """ Check equal energies for Gromacs -> Amber conversion of Amber FF """
        top = chem.load_file(get_fn(os.path.join('03.AlaGlu', 'topol.top')))
        gro = chem.load_file(get_fn(os.path.join('03.AlaGlu', 'conf.gro')))
        parm = chem.amber.AmberParm.from_structure(top)
        parm.write_parm(get_fn('ala_glu.parm7', written=True))
        parm = chem.load_file(get_fn('ala_glu.parm7', written=True))

        sysg = top.createSystem()
        sysa = parm.createSystem()

        cong = mm.Context(sysg, mm.VerletIntegrator(0.001), CPU)
        cona = mm.Context(sysa, mm.VerletIntegrator(0.001), CPU)

        cong.setPositions(gro.positions)
        cona.setPositions(gro.positions)

        self._check_energies(top, cong, parm, cona)

    @unittest.skipIf(not has_openmm, "Cannot test without OpenMM")
    def testEnergyComplicated(self):
        """ Check equal energies for Gmx -> Amber conversion of complex FF """
        top = chem.load_file(get_fn(os.path.join('12.DPPC', 'topol2.top')))
        gro = chem.load_file(get_fn(os.path.join('12.DPPC', 'conf.gro')))
        parm = chem.amber.AmberParm.from_structure(top)
        parm.write_parm(get_fn('dppc.parm7', written=True))
        parm = chem.load_file(get_fn('dppc.parm7', written=True))

        sysg = top.createSystem()
        sysa = parm.createSystem()

        cong = mm.Context(sysg, mm.VerletIntegrator(0.001), CPU)
        cona = mm.Context(sysa, mm.VerletIntegrator(0.001), CPU)

        cong.setPositions(gro.positions)
        cona.setPositions(gro.positions)

        self._check_energies(top, cong, parm, cona)


    def _check_energies(self, parm1, con1, parm2, con2):
        ene1 = chem.openmm.utils.energy_decomposition(parm1, con1)
        ene2 = chem.openmm.utils.energy_decomposition(parm2, con2)

        all_terms = set(ene1.keys()) | set(ene2.keys())

        for term in all_terms:
            if term not in ene1:
                self.assertAlmostEqual(ene2[term], 0)
            elif term not in ene2:
                self.assertAlmostEqual(ene1[term], 0)
            else:
                self.assertRelativeEqual(ene2[term], ene1[term], places=5)

@unittest.skipIf(not has_openmm, "Cannot test without OpenMM")
class TestOpenMMToAmber(TestCase):
    """
    Tests that OpenMM system/topology combo can be translated to other formats
    """

    def testSimple(self):
        """ Test OpenMM System/Topology -> Amber prmtop conversion """
        parm = chem.load_file(get_fn('ash.parm7'), get_fn('ash.rst7'))
        system = parm.createSystem()
        chem.amber.AmberParm.from_structure(
                chem.openmm.load_topology(parm.topology, system)
        ).write_parm(get_fn('ash_from_omm.parm7', written=True))
        parm2 = chem.load_file(get_fn('ash_from_omm.parm7', written=True))
        system2 = parm2.createSystem()
        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(parm.positions)
        con2.setPositions(parm.positions)

        self._check_energies(parm, con1, parm2, con2)

    def _check_energies(self, parm1, con1, parm2, con2):
        ene1 = chem.openmm.utils.energy_decomposition(parm1, con1)
        ene2 = chem.openmm.utils.energy_decomposition(parm2, con2)

        all_terms = set(ene1.keys()) | set(ene2.keys())

        for term in all_terms:
            if term not in ene1:
                self.assertAlmostEqual(ene2[term], 0)
            elif term not in ene2:
                self.assertAlmostEqual(ene1[term], 0)
            else:
                self.assertRelativeEqual(ene2[term], ene1[term], places=5)

@unittest.skipIf(not has_openmm, "Cannot test without OpenMM")
class TestOpenMMToGromacs(TestCase):
    """
    Tests that OpenMM system/topology combo can be translated to other formats
    """

    def testSimple(self):
        """ Test OpenMM System/Topology -> Gromacs topology conversion """
        parm = chem.load_file(get_fn('ash.parm7'), get_fn('ash.rst7'))
        system = parm.createSystem()
        chem.gromacs.GromacsTopologyFile.from_structure(
                chem.openmm.load_topology(parm.topology, system)
        ).write(get_fn('ash_from_omm.top', written=True))
        parm2 = chem.gromacs.GromacsTopologyFile(get_fn('ash_from_omm.top', written=True))
        system2 = parm2.createSystem()
        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(parm.positions)
        con2.setPositions(parm.positions)

        self._check_energies(parm, con1, parm2, con2)

    def _check_energies(self, parm1, con1, parm2, con2):
        ene1 = chem.openmm.utils.energy_decomposition(parm1, con1)
        ene2 = chem.openmm.utils.energy_decomposition(parm2, con2)

        all_terms = set(ene1.keys()) | set(ene2.keys())

        for term in all_terms:
            if term not in ene1:
                self.assertAlmostEqual(ene2[term], 0)
            elif term not in ene2:
                self.assertAlmostEqual(ene1[term], 0)
            else:
                self.assertRelativeEqual(ene2[term], ene1[term], places=5)

