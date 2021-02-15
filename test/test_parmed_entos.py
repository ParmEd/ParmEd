""" Test cases for integrating with Entos """
from unittest import skipIf, TestCase
import numpy as np

import parmed as pmd
from parmed.entos import HAS_ENTOS, HAS_MISTRAL, to_entos_molecule, to_entos_qmmm_system
from parmed.entos.imports import constants, QMMMSystem
from parmed.entos.imports import Molecule
from utils import get_fn

@skipIf(not HAS_ENTOS, "Entos test cases require Entos packages to be installed")
class TestEntos(TestCase):

    def test_to_entos_molecule(self):
        """ Tests creation of ParmEd molecule """
        struct = pmd.load_file(get_fn("ash.parm7"), get_fn("ash.rst7"))
        molecule = to_entos_molecule(struct)
        self.assertIsInstance(molecule, Molecule)
        self.assertEqual(molecule.atomic_numbers.size, len(struct.atoms))
        for a1, a2 in zip(molecule.atomic_numbers, struct.atoms):
            self.assertEqual(a1, a2.atomic_number)
        np.testing.assert_allclose(
            molecule.geometry * constants.cf("bohr", "angstrom"),
            struct.coordinates,
            atol=1e-6,
        )

    def test_to_entos_molecule_requires_coordinates(self):
        """ Check that to_entos_molecule requires coordinates """
        with self.assertRaises(ValueError):
            to_entos_molecule(pmd.load_file(get_fn("ash.parm7")))

    def test_to_entos_qmmm_system(self):
        struct = pmd.load_file(get_fn("ala3_solv.parm7"), get_fn("ala3_solv.rst7"))
        qmmm_system = to_entos_qmmm_system(struct, ':1-3')
        self.assertIsInstance(qmmm_system, QMMMSystem)
