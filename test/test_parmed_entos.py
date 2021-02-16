""" Test cases for integrating with Entos """
from unittest import skipUnless, TestCase
import numpy as np

import parmed as pmd
from parmed.entos import HAS_ENTOS, HAS_MISTRAL, to_entos_molecule, to_entos_qmmm_system
from parmed.entos.imports import constants, QMMMSystem
from parmed.entos.imports import Molecule
from utils import get_fn

@skipUnless(HAS_ENTOS, "Entos test cases require Entos packages to be installed")
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

    @skipUnless(HAS_MISTRAL, "QM/MM system creation tests require mistral")
    def test_to_entos_qmmm_system(self):
        """ Test the creation of an Entos QMMMSystem from Structure """
        struct = pmd.load_file(get_fn("ala3_solv.parm7"), get_fn("ala3_solv.rst7"))
        qmmm_system = to_entos_qmmm_system(struct, ':1-3')
        self.assertIsInstance(qmmm_system, QMMMSystem)

    @skipUnless(HAS_MISTRAL, "QM/MM system creation tests require mistral")
    def test_to_entos_qmmm_system_requires_coordinates(self):
        """ Check that to_entos_qmmm_system requires coordinates """
        with self.assertRaises(ValueError):
            to_entos_qmmm_system(pmd.load_file(get_fn("ala3_solv.parm7")), ":1-3")
