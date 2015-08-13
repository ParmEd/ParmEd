""" Test various topology format conversions """
from __future__ import print_function, division, absolute_import

import os
import numpy as np
from parmed import load_file, gromacs, amber, openmm, charmm
from parmed.exceptions import GromacsWarning
from parmed.gromacs._gromacsfile import GromacsFile
from parmed import unit as u
try:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    CPU = mm.Platform.getPlatformByName('CPU')
    has_openmm = True
except ImportError:
    has_openmm = False
import unittest
from utils import get_fn, get_saved_fn, diff_files, TestCaseRelative, FileIOTestCase
import warnings

class TestAmberToGromacs(FileIOTestCase, TestCaseRelative):
    """ Tests converting Amber prmtop files to Gromacs topologies """

    def testBenzeneCyclohexane(self):
        """ Test converting binary liquid from Amber prmtop to Gromacs top """
        parm = load_file(get_fn('benzene_cyclohexane_10_500.prmtop'),
                              get_fn('benzene_cyclohexane_10_500.inpcrd'))
        top = gromacs.GromacsTopologyFile.from_structure(parm)
        self.assertEqual(top.combining_rule, 'lorentz')
        groname = get_fn('benzene_cyclohexane_10_500.gro', written=True)
        gromacs.GromacsGroFile.write(parm, groname, precision=8)
        gro = gromacs.GromacsGroFile.parse(groname)
        self.assertEqual(len(gro.atoms), len(parm.atoms))
        np.testing.assert_allclose(gro.box, parm.box)
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
        # Check that Gromacs topology is given the correct box information when
        # generated from a Structure
        gromacs.GromacsGroFile.write(top, groname, precision=8)
        gro = gromacs.GromacsGroFile.parse(groname)
        self.assertEqual(len(gro.atoms), len(parm.atoms))
        for a1, a2 in zip(gro.atoms, parm.atoms):
            self.assertEqual(a1.residue.name, a2.residue.name)
            self.assertEqual(a1.residue.idx, a2.residue.idx)
            self.assertEqual(a1.name, a2.name)
            self.assertAlmostEqual(a1.xx, a2.xx)
            self.assertAlmostEqual(a1.xy, a2.xy)
            self.assertAlmostEqual(a1.xz, a2.xz)
        np.testing.assert_allclose(gro.box, parm.box)
        np.testing.assert_allclose(top.box, parm.box)

class TestGromacsToAmber(FileIOTestCase, TestCaseRelative):
    """ Tests converting Gromacs top/gro files to Amber """

    def testSimple(self):
        """ Tests converting standard Gromacs system into Amber prmtop """
        top = load_file(get_fn(os.path.join('03.AlaGlu', 'topol.top')))
        self.assertEqual(top.combining_rule, 'lorentz')
        parm = amber.AmberParm.from_structure(top)
        parm.write_parm(get_fn('ala_glu.parm7', written=True))
        parm = load_file(get_fn('ala_glu.parm7', written=True))
        self.assertIsInstance(parm, amber.AmberParm)
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

    def testChamber(self):
        """ Tests converting standard Gromacs system into Chamber prmtop """
        top = load_file(get_fn('1aki.charmm27.solv.top'),
                        xyz=get_fn('1aki.charmm27.solv.gro'))
        parm = amber.ChamberParm.from_structure(top)
        parm.write_parm(get_fn('1aki.charmm27_fromgmx.parm7', written=True))
        self.assertTrue(diff_files(get_fn('1aki.charmm27_fromgmx.parm7',
                                          written=True),
                                   get_saved_fn('1aki.charmm27_fromgmx.parm7'),
                                   relative_error=1e-8)
        )

    def testGeometricCombiningRule(self):
        """ Tests converting geom. comb. rule from Gromacs to Amber """
        top = load_file(os.path.join(get_fn('05.OPLS'), 'topol.top'),
                        xyz=os.path.join(get_fn('05.OPLS'), 'conf.gro'))
        self.assertEqual(top.combining_rule, 'geometric')
        del top.rb_torsions[:]
        parm = amber.AmberParm.from_structure(top)
        self.assertEqual(parm.combining_rule, 'geometric')
        parm.write_parm(get_fn('opls.parm7', written=True))
        self.assertTrue(diff_files(get_fn('opls.parm7', written=True),
                                   get_saved_fn('opls.parm7'))
        )

    @unittest.skipIf(not has_openmm, "Cannot test without OpenMM")
    def testGeometricCombiningRuleEnergy(self):
        """ Tests converting geom. comb. rule energy from Gromacs to Amber """
        top = load_file(os.path.join(get_fn('05.OPLS'), 'topol.top'),
                        xyz=os.path.join(get_fn('05.OPLS'), 'conf.gro'))
        self.assertEqual(top.combining_rule, 'geometric')
        del top.rb_torsions[:]
        parm = load_file(get_saved_fn('opls.parm7'),
                         xyz=os.path.join(get_fn('05.OPLS'), 'conf.gro'))
        self.assertEqual(parm.combining_rule, 'geometric')
        self.assertFalse(parm.has_NBFIX())

        sysg = top.createSystem()
        sysa = parm.createSystem()

        cong = mm.Context(sysg, mm.VerletIntegrator(0.001), CPU)
        cona = mm.Context(sysa, mm.VerletIntegrator(0.001), CPU)

        cong.setPositions(top.positions)
        cona.setPositions(top.positions)
        
        self._check_energies(top, cong, parm, cona)

    @unittest.skipIf(not has_openmm, "Cannot test without OpenMM")
    def testEnergySimple(self):
        """ Check equal energies for Gromacs -> Amber conversion of Amber FF """
        top = load_file(get_fn(os.path.join('03.AlaGlu', 'topol.top')))
        gro = load_file(get_fn(os.path.join('03.AlaGlu', 'conf.gro')))
        parm = amber.AmberParm.from_structure(top)
        parm.write_parm(get_fn('ala_glu.parm7', written=True))
        parm = load_file(get_fn('ala_glu.parm7', written=True))

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
        top = load_file(get_fn(os.path.join('12.DPPC', 'topol2.top')))
        gro = load_file(get_fn(os.path.join('12.DPPC', 'conf.gro')))
        parm = amber.AmberParm.from_structure(top)
        parm.write_parm(get_fn('dppc.parm7', written=True))
        parm = load_file(get_fn('dppc.parm7', written=True))

        sysg = top.createSystem()
        sysa = parm.createSystem()

        cong = mm.Context(sysg, mm.VerletIntegrator(0.001), CPU)
        cona = mm.Context(sysa, mm.VerletIntegrator(0.001), CPU)

        cong.setPositions(gro.positions)
        cona.setPositions(gro.positions)

        self._check_energies(top, cong, parm, cona)


    def _check_energies(self, parm1, con1, parm2, con2):
        ene1 = openmm.utils.energy_decomposition(parm1, con1)
        ene2 = openmm.utils.energy_decomposition(parm2, con2)

        all_terms = set(ene1.keys()) | set(ene2.keys())

        for term in all_terms:
            if term not in ene1:
                self.assertAlmostEqual(ene2[term], 0)
            elif term not in ene2:
                self.assertAlmostEqual(ene1[term], 0)
            else:
                self.assertRelativeEqual(ene2[term], ene1[term], places=5)

class TestAmberToCharmm(FileIOTestCase, TestCaseRelative):
    """ Tests converting Amber files to CHARMM """

    def testSimple(self):
        """ Tests converting simple Amber system to CHARMM PSF/parameters """
        parm = load_file(get_fn('trx.prmtop'), get_fn('trx.inpcrd'))
        parm.save(get_fn('amber_to_charmm.psf', written=True))
        params = charmm.CharmmParameterSet.from_structure(parm)
        params.write(str=get_fn('amber_to_charmm.str', written=True))

        self.assertTrue(
                diff_files(get_saved_fn('amber_to_charmm.psf'),
                           get_fn('amber_to_charmm.psf', written=True)
                )
        )
        self.assertTrue(
                diff_files(get_saved_fn('amber_to_charmm.str'),
                           get_fn('amber_to_charmm.str', written=True)
                )
        )

@unittest.skipIf(not has_openmm, "Cannot test without OpenMM")
class TestOpenMMToAmber(FileIOTestCase, TestCaseRelative):
    """
    Tests that OpenMM system/topology combo can be translated to other formats
    """

    def testSimple(self):
        """ Test OpenMM System/Topology -> Amber prmtop conversion """
        parm = load_file(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem()
        amber.AmberParm.from_structure(
                openmm.load_topology(parm.topology, system)
        ).write_parm(get_fn('ash_from_omm.parm7', written=True))
        parm2 = load_file(get_fn('ash_from_omm.parm7', written=True))
        system2 = parm2.createSystem()
        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(parm.positions)
        con2.setPositions(parm.positions)

        self._check_energies(parm, con1, parm2, con2)

    def _check_energies(self, parm1, con1, parm2, con2):
        ene1 = openmm.utils.energy_decomposition(parm1, con1)
        ene2 = openmm.utils.energy_decomposition(parm2, con2)

        all_terms = set(ene1.keys()) | set(ene2.keys())

        for term in all_terms:
            if term not in ene1:
                self.assertAlmostEqual(ene2[term], 0)
            elif term not in ene2:
                self.assertAlmostEqual(ene1[term], 0)
            else:
                self.assertRelativeEqual(ene2[term], ene1[term], places=5)

@unittest.skipIf(not has_openmm, "Cannot test without OpenMM")
class TestOpenMMToGromacs(FileIOTestCase, TestCaseRelative):
    """
    Tests that OpenMM system/topology combo can be translated to other formats
    """

    def testSimple(self):
        """ Test OpenMM System/Topology -> Gromacs topology conversion """
        parm = load_file(get_fn('ash.parm7'), get_fn('ash.rst7'))
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem()
        gromacs.GromacsTopologyFile.from_structure(
                openmm.load_topology(parm.topology, system)
        ).write(get_fn('ash_from_omm.top', written=True))
        parm2 = gromacs.GromacsTopologyFile(get_fn('ash_from_omm.top', written=True))
        system2 = parm2.createSystem()
        con1 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con2 = mm.Context(system, mm.VerletIntegrator(0.001), CPU)
        con1.setPositions(parm.positions)
        con2.setPositions(parm.positions)

        self._check_energies(parm, con1, parm2, con2)

    def _check_energies(self, parm1, con1, parm2, con2):
        ene1 = openmm.utils.energy_decomposition(parm1, con1)
        ene2 = openmm.utils.energy_decomposition(parm2, con2)

        all_terms = set(ene1.keys()) | set(ene2.keys())

        for term in all_terms:
            if term not in ene1:
                self.assertAlmostEqual(ene2[term], 0)
            elif term not in ene2:
                self.assertAlmostEqual(ene1[term], 0)
            else:
                self.assertRelativeEqual(ene2[term], ene1[term], places=5)

