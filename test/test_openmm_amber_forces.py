"""
Contains unittests for running OpenMM calculations using the Amber file parsers
"""
from __future__ import division, print_function, absolute_import

from collections import defaultdict
from copy import copy
from math import sqrt
import os
from parmed import Structure, topologyobjects, load_file
from parmed.exceptions import ParameterError
from parmed.amber import AmberParm, ChamberParm, Rst7
from parmed.openmm import load_topology, energy_decomposition_system
import parmed.unit as u
from parmed.utils.six.moves import range, zip
import parmed.tools as PT
import sys
import unittest
from utils import get_fn, CPU, mm, app, has_openmm, TestCaseRelative, QuantityTestCase

# OpenMM NonbondedForce methods are enumerated values. From NonbondedForce.h,
# they are:
#   0 - NoCutoff
#   1 - CutoffNonPeriodic
#   2 - CutoffPeriodic
#   3 - Ewald
#   4 - PME

def energy_decomposition(parm, context):
    from parmed.openmm.utils import energy_decomposition
    ret = defaultdict(float)
    for key, val in energy_decomposition(parm, context).items():
        ret[key] = val
    return ret

@unittest.skipUnless(has_openmm, "Cannot test without OpenMM")
class TestAmberParm(TestCaseRelative, QuantityTestCase):

    def test_gas_energy_conf_1(self):
        """ Compare Amber and OpenMM gas phase energies and forces (topology 1) """
        parm = AmberParm(
            get_fn('hydrogen-peroxide_T0_1.prmtop'), get_fn('hydrogen-peroxide_T0_1.rst7')
        )
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem() # Default, no cutoff
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#BOND    =        0.0319  ANGLE   =        2.1690  DIHED      =        0.2218
#VDWAALS =        0.0000  EEL     =        0.0000  HBOND      =        0.0000
#1-4 VDW =        0.0000  1-4 EEL =        0.0000  RESTRAINT  =        0.0000
        # Compare OpenMM energies with the Amber energies (above)
        self.assertAlmostEqual(energies['bond'], 0.0319, places=4)
        self.assertAlmostEqual(energies['angle'], 2.1690, places=4)
        self.assertAlmostEqual(energies['dihedral'], 0.2218, places=4)
        self.assertRelativeEqual(energies['nonbonded'], 0.0, places=3)

        # Now test the forces to make sure that they are computed correctly in
        # the presence of extra points
        pstate = sim.context.getState(getForces=True)
        pf = pstate.getForces().value_in_unit(u.kilojoule_per_mole/u.nanometer)
# gmx forces:
#     f[    0]={ 5.37936e+02, -5.90363e+02, -6.71315e+01}
#     f[    1]={-5.37936e+02,  5.90363e+02, -6.71315e+01}
#     f[    2]={ 9.92309e+01,  5.90363e+02,  6.71315e+01}
#     f[    3]={-9.92309e+01, -5.90363e+02,  6.71315e+01}
        gf = [  ( 5.37936e+02, -5.90363e+02, -6.71315e+01),
                (-5.37936e+02,  5.90363e+02, -6.71315e+01),
                ( 9.92309e+01,  5.90363e+02,  6.71315e+01),
                (-9.92309e+01, -5.90363e+02,  6.71315e+01)]
        for p, s in zip(pf, gf):
            for x1, x2 in zip(p, s):
                self.assertAlmostEqual(x1, x2, places=3)

    def test_gas_energy_conf_2(self):
        """ Compare Amber and OpenMM gas phase energies and forces (topology 2) """
        parm = AmberParm(
            get_fn('hydrogen-peroxide_T0_2.prmtop'), get_fn('hydrogen-peroxide_T0_2.rst7')
        )
        self.assertEqual(parm.combining_rule, 'lorentz')
        system = parm.createSystem() # Default, no cutoff
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#BOND    =        0.0319  ANGLE   =        2.1690  DIHED      =       -2.7931
#VDWAALS =        0.0000  EEL     =        0.0000  HBOND      =        0.0000
#1-4 VDW =        0.0000  1-4 EEL =        0.0000  RESTRAINT  =        0.0000
        # Compare OpenMM energies with the Amber energies (above)
        self.assertAlmostEqual(energies['bond'], 0.0319, places=4)
        self.assertAlmostEqual(energies['angle'], 2.1690, places=4)
        self.assertAlmostEqual(energies['dihedral'],-2.7931, places=4)
        self.assertRelativeEqual(energies['nonbonded'], 0.0, places=3)

        # Now test the forces to make sure that they are computed correctly in
        # the presence of extra points
        pstate = sim.context.getState(getForces=True)
        pf = pstate.getForces().value_in_unit(u.kilojoule_per_mole/u.nanometer)
# gmx forces:
#     f[    0]={ 3.14357e+02, -5.90363e+02,  2.11410e+02}
#     f[    1]={-3.14357e+02,  5.90363e+02,  2.11410e+02}
#     f[    2]={ 2.70641e+02,  5.90363e+02, -2.11410e+02}
#     f[    3]={-2.70641e+02, -5.90363e+02, -2.11410e+02}
        gf = [  ( 3.14357e+02, -5.90363e+02,  2.11410e+02), 
                (-3.14357e+02,  5.90363e+02,  2.11410e+02), 
                ( 2.70641e+02,  5.90363e+02, -2.11410e+02), 
                (-2.70641e+02, -5.90363e+02, -2.11410e+02)] 
        for p, s in zip(pf, gf):
            for x1, x2 in zip(p, s):
                self.assertAlmostEqual(x1, x2, places=3)

    def test_gas_energy_conf_3(self):
        """ Compare Amber and OpenMM gas phase energies and forces (topology 3) """
        parm = AmberParm(get_fn('methanol_T0.prmtop'), get_fn('methanol_T0.rst7'))
        self.assertEqual(parm.combining_rule, 'geometric')
        system = parm.createSystem() # Default, no cutoff
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#BOND    =        0.0224  ANGLE   =        0.0469  DIHED      =        0.0039
#VDWAALS =        0.0000  EEL     =        0.0000  HBOND      =        0.0000
#1-4 VDW =        0.0000  1-4 EEL =        3.3789  RESTRAINT  =        0.0000
        # Compare OpenMM energies with the Amber energies (above)
        self.assertAlmostEqual(energies['bond'], 0.0224, places=4)
        self.assertAlmostEqual(energies['angle'], 0.0469, places=4)
        self.assertAlmostEqual(energies['dihedral'], 0.0039, places=4)
        self.assertRelativeEqual(energies['nonbonded'], 3.3789, places=3)

        # Now test the forces to make sure that they are computed correctly in
        # the presence of extra points
        pstate = sim.context.getState(getForces=True)
        pf = pstate.getForces().value_in_unit(u.kilojoule_per_mole/u.nanometer)
# gmx forces:
#     f[    0]={ 2.19935e+01,  8.60143e+01,  1.77895e+02}
#     f[    1]={ 5.25347e+01,  3.46119e+01, -9.73738e+01}
#     f[    2]={ 1.76138e+01, -7.13778e+01, -1.08807e+00}
#     f[    3]={-3.74994e+01, -5.39602e+01, -6.75679e+01}
#     f[    4]={-8.35142e+01, -7.86616e+01, -1.60611e+02}
#     f[    5]={ 2.88716e+01,  8.33735e+01,  1.48746e+02}

        gf = [  ( 2.19935e+01,  8.60143e+01,  1.77895e+02),
                ( 5.25347e+01,  3.46119e+01, -9.73738e+01),
                ( 1.76138e+01, -7.13778e+01, -1.08807e+00),
                (-3.74994e+01, -5.39602e+01, -6.75679e+01),
                (-8.35142e+01, -7.86616e+01, -1.60611e+02),
                ( 2.88716e+01,  8.33735e+01,  1.48746e+02)]
        for p, s in zip(pf, gf):
            for x1, x2 in zip(p, s):
                self.assertAlmostEqual(x1, x2, places=3)


    def test_gas_energy_conf_4(self):
        """ Compare Amber and OpenMM gas phase energies and forces (topology 4) """
        parm = AmberParm(get_fn('ethanol_T0.prmtop'), get_fn('ethanol_T0.rst7'))
        self.assertEqual(parm.combining_rule, 'geometric')
        system = parm.createSystem() # Default, no cutoff
        integrator = mm.VerletIntegrator(1.0*u.femtoseconds)
        sim = app.Simulation(parm.topology, system, integrator, platform=CPU)
        sim.context.setPositions(parm.positions)
        energies = energy_decomposition(parm, sim.context)
#BOND    =        0.0239  ANGLE   =        0.0298  DIHED      =        0.0093
#VDWAALS =        0.0000  EEL     =        6.7526  HBOND      =        0.0000
#1-4 VDW =        0.0492  1-4 EEL =       -6.0430  RESTRAINT  =        0.0000
        # Compare OpenMM energies with the Amber energies (above)
        self.assertAlmostEqual(energies['bond'], 0.0239, places=4)
        self.assertAlmostEqual(energies['angle'], 0.0298, places=4)
        self.assertAlmostEqual(energies['dihedral'], 0.0093, places=4)
        self.assertRelativeEqual(energies['nonbonded'], 0.0000+6.7526+0.0492-6.0430, places=3)

        # Now test the forces to make sure that they are computed correctly in
        # the presence of extra points
        pstate = sim.context.getState(getForces=True)
        pf = pstate.getForces().value_in_unit(u.kilojoule_per_mole/u.nanometer)
# gmx forces:
#      f[    0]={ 1.73101e+02,  5.67250e+01, -4.02950e+01}
#      f[    1]={-8.13704e+00, -1.79612e+01,  9.88537e+01}
#      f[    2]={-2.83120e+01,  6.23352e+00, -5.47393e+00}
#      f[    3]={-1.03312e+01, -1.00966e+01, -5.12129e+00}
#      f[    4]={-1.69636e+02,  3.34850e+01, -1.73612e+02}
#      f[    5]={ 4.19932e+00, -2.58283e+00,  1.29999e+02}
#      f[    6]={ 3.02865e+01, -6.68331e+00, -2.99153e+01}
#      f[    7]={ 1.00113e+02, -5.22480e+01,  4.80526e+01}
#      f[    8]={-9.12828e+01, -6.87157e+00, -2.24877e+01}

        gf = [  ( 1.73101e+02,  5.67250e+01, -4.02950e+01),
                (-8.13704e+00, -1.79612e+01,  9.88537e+01),
                (-2.83120e+01,  6.23352e+00, -5.47393e+00),
                (-1.03312e+01, -1.00966e+01, -5.12129e+00),
                (-1.69636e+02,  3.34850e+01, -1.73612e+02),
                ( 4.19932e+00, -2.58283e+00,  1.29999e+02),
                ( 3.02865e+01, -6.68331e+00, -2.99153e+01),
                ( 1.00113e+02, -5.22480e+01,  4.80526e+01),
                (-9.12828e+01, -6.87157e+00, -2.24877e+01)]
        for p, s in zip(pf, gf):
            for x1, x2 in zip(p, s):
                self.assertAlmostEqual(x1, x2, places=3)


