"""
Useful functions for the test cases
"""
import os
from os.path import join, split, abspath
import random
import tempfile
import unittest
import warnings
from pathlib import Path
from typing import Union

import numpy as np

from parmed import gromacs, openmm
from parmed.utils.six import string_types
from parmed.utils.six.moves import zip

try:
    from simtk import openmm as mm
    from simtk.openmm import app
    openmm_version = tuple([int(x) for x in mm.__version__.split('.')])
    CPU = mm.Platform.getPlatformByName('CPU')
    Reference = mm.Platform.getPlatformByName('Reference')
    has_openmm = True
except ImportError:
    has_openmm = False
    app = openmm_version = CPU = Reference = mm = None

try:
    import networkx as nx
    has_networkx = True
except ImportError:
    has_networkx = False

try:
    from lxml import etree
    has_lxml = True
except ImportError:
    has_lxml = False

try:
    from string import uppercase
except ImportError:
    from string import ascii_uppercase as uppercase

run_all_tests = os.getenv('PARMED_RUN_ALL_TESTS') is not None

HAS_GROMACS = os.path.isdir(gromacs.GROMACS_TOPDIR)

class QuantityTestCase(unittest.TestCase):

    def assertAlmostEqualQuantities(self, item1, item2, places=6):
        try:
            val1 = item1.value_in_unit(item1.unit)
            val2 = item2.value_in_unit(item1.unit)
        except TypeError:
            raise self.failureException('Incompatible units %s and %s' %
                                        (item1.unit, item2.unit))
        try:
            if len(val1) != len(val2):
                raise self.failureException('collections are different lengths')
            for x, y in zip(val1, val2):
                self.assertAlmostEqual(x, y, places=places)
        except TypeError:
            self.assertAlmostEqual(val1, val2, places=places)

class TestCaseRelative(unittest.TestCase):

    def assertRelativeEqual(self, val1, val2, places=7, delta=None):
        if val1 == val2: return
        try:
            ratio = val1 / val2
        except ZeroDivisionError:
            return self.assertAlmostEqual(val1, val2, places=places)
        else:
            if delta is None:
                if abs(round(ratio - 1, places)) == 0:
                    return
                raise self.failureException('%s != %s with relative tolerance %g (%f)' %
                                            (val1, val2, 10**-places, ratio))
            else:
                if abs(ratio - 1) < delta:
                    return
                raise self.failureException('%s != %s with relative tolerance %g (%f)' %
                                            (val1, val2, delta, ratio))

class EnergyTestCase(TestCaseRelative):

    def check_energies(self, parm1, con1, parm2, con2):
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

class FileIOTestCase(unittest.TestCase):

    def setUp(self):
        self._temporary_directory = tempfile.TemporaryDirectory(suffix="pmdtest")

    def tearDown(self):
        self._temporary_directory.cleanup()

    def get_fn(self, filename: str, written: bool = False, saved: bool = False):
        assert not (written and saved), "Cannot get saved and written"
        if written:
            return join(self._temporary_directory.name, filename)
        elif saved:
            return str(Path(__file__).parent / "files" / "saved" / filename)
        return get_fn(filename)

    def _empty_writes(self):
        for fname in os.listdir(self._temporary_directory.name):
            os.unlink(os.path.join(self._temporary_directory.name, fname))

def get_fn(filename):
    """
    Gets the full path of the file name for a particular test file

    Parameters
    ----------
    filename : str
        Name of the file to get
    written : bool=False
        Was this file written by a test? (If so, it is put in a different
        location)

    Returns
    -------
    str
        Name of the test file with the full path location
    """
    return str(Path(__file__).parent / "files" / filename)

def get_saved_fn(filename):
    """
    Gets the full path of a file name of a saved test case that is used for
    comparison with a generated file

    Parameters
    ----------
    filename : str
        Name of the file to get

    Returns
    -------
    str
        Name of the test file with the full path location
    """
    return str(Path(__file__).parent / "files" / "saved" / filename)

def diff_files(file1, file2, ignore_whitespace=True,
               absolute_error=None, relative_error=None,
               comment=None, spacechar=None):
    """
    Compares 2 files line-by-line

    Parameters
    ----------
    file1 : str or file-like
        Name of the first file to compare or first file object to compare
    file2 : str or file-like
        Name of the second file to compare or second file object to compare
    ignore_whitespace : bool=True
        If true, ignores differences in leading and trailing whitespace
    absolute_error : float=None
        If set, only differences greater than the absolute error will trigger
        failures. Cannot be used with relative_error
    relative_error : float=None
        If set, only relative differences greater than the relative error will
        trigger failures. Cannot be used with absolute_error
    comment : str or None
        The character to identify comments in this file
    spacechar : str or None
        A collection of characters to turn into spaces (to facilitate proper
        tokenization)

    Returns
    -------
    bool
        True if files match. False if they do not or one file does not exist

    Notes
    -----
    This routine is not protected against bad types of input. AttributeError may
    be raised if readline() is not implemented on any file-like object passed in
    """
    if absolute_error is not None and relative_error is not None:
        raise ValueError('Cannot specify absolute_error AND relative_error')
    if absolute_error is not None: absolute_error = float(absolute_error)
    if relative_error is not None: relative_error = float(relative_error)
    if isinstance(file1, string_types):
        try:
            f1 = open(file1, 'r')
        except IOError:
            print('Could not find %s' % file1)
            return False
    else:
        f1 = file1
        file1 = str(file1)
    if isinstance(file2, string_types):
        try:
            f2 = open(file2, 'r')
        except IOError:
            print('Could not find %s' % file2)
            return False
    else:
        f2 = file2
        file2 = str(file2)
    try:
        l1 = f1.readline()
        l2 = f2.readline()
        i = 1
        same = True
        if ignore_whitespace:
            while l1 or l2:
                while l1 and l1[0] == comment:
                    l1 = f1.readline()
                while l2 and l2[0] == comment:
                    l2 = f2.readline()
                if l1.strip() != l2.strip():
                    if l1.startswith('%VERSION') and l2.startswith('%VERSION'):
                        l1 = f1.readline()
                        l2 = f2.readline()
                        continue
                    if not detailed_diff(l1,l2,absolute_error,relative_error,spacechar):
                        same = False
                        record_diffs(i, file1, file2, l1, l2)
                l1 = f1.readline()
                l2 = f2.readline()
                i += 1
        else:
            while l1 and l2:
                if l1 != l2:
                    if 'At date:' in l1 and 'At date:' in l2:
                        l1 = f1.readline()
                        l2 = f2.readline()
                        continue
                    if l1.startswith('%VERSION') and l2.startswith('%VERSION'):
                        l1 = f1.readline()
                        l2 = f2.readline()
                        continue
                    if not detailed_diff(l1,l2,absolute_error,relative_error,spacechar):
                        same = False
                        record_diffs(i, file1, file2, l1, l2)
                l1 = f1.readline()
                l2 = f2.readline()
                i += 1
    finally:
        f1.close()
        f2.close()

    return same

def record_diffs(i, f1, f2, l1, l2):
    if not os.path.isdir(get_fn('diffs')):
        os.makedirs(get_fn('diffs'))
    f = open(os.path.join(get_fn('diffs'), 'TEST_FAILURES.diff'), 'a')
    f.write('# diff %s %s [line %d]\n' % (f1, f2, i))
    f.write('< %s> %s' % (l1, l2))
    f.close()

def detailed_diff(l1, l2, absolute_error=None, relative_error=None, spacechar=None):
    """
    Check individual fields to make sure numbers are numerically equal if the
    lines differ. Also ignore fields that appear to be a file name, since those
    will be system-dependent
    """
    fdir = os.path.split(get_fn('writes'))[0]
    if spacechar is not None:
        for char in spacechar:
            l1 = l1.replace(char, ' ')
            l2 = l2.replace(char, ' ')
    w1 = l1.split()
    w2 = l2.split()
    if len(w1) != len(w2):
        return False
    for wx, wy in zip(w1, w2):
        try:
            wx = float(wx)
            wy = float(wy)
        except ValueError:
            y_is_filename = wy.startswith(fdir) or wy.startswith(tempfile.tempdir)
            x_is_filename = isinstance(wx, str) and (wx.startswith(fdir) or wx.startswith(tempfile.tempdir))
            if isinstance(wx, float) or (wx != wy and not (y_is_filename and x_is_filename)):
                return False
        else:
            if wx != wy:
                if absolute_error is not None and abs(wx - wy) > absolute_error:
                    return False
                elif relative_error is not None:
                    if wx == 0 or wy == 0 and abs(wx - wy) > relative_error:
                        return False
                    if abs((wx / wy) - 1) > relative_error:
                        return False
                elif absolute_error is None and relative_error is None:
                    return False
    return True

def which(prog):
    """ Like the Unix program ``which``

    Parameters
    ----------
    prog : str
        Name of the program to look for in PATH

    Returns
    -------
    path
        The full path of the program, or None if it cannot be found
    """
    def is_exe(filename):
        return os.path.exists(filename) and os.access(filename, os.X_OK)

    # See if full path is provided
    fpath, fname = os.path.split(prog)
    if fpath:
        if is_exe(prog):
            return prog
        return None

    # If not, see if the program exists anywhere in PATH
    pathparts = os.environ['PATH'].split(os.path.pathsep)
    for part in pathparts:
        trial = os.path.join(part, prog)
        if is_exe(trial):
            return trial
    return None

def create_random_structure(parametrized, novalence=False):
    """ Create a random Structure with random attributes

    Parameters
    ----------
    parametrized : bool
        If True, add at least two of all kinds of parameters to the
        generated random structure. If False, just fill in the atoms and
        residues and some random valence terms, but no "types"
    novalence : bool, optional
        If True, no valence terms will be added. Default is False. This is
        set to False if parametrized is True
    """
    from parmed.topologyobjects import (Atom, Bond, AtomType, BondType,
            AngleType, DihedralType, ImproperType, CmapType, OutOfPlaneBendType,
            StretchBendType, TorsionTorsionType, AmoebaNonbondedExceptionType,
            Angle, UreyBradley, Dihedral, Improper, Cmap, TrigonalAngle,
            OutOfPlaneBend, StretchBend, PiTorsion, TorsionTorsion,
            AcceptorDonor, Group, ChiralFrame, MultipoleFrame,
            NonbondedException, RBTorsionType)
    from parmed import structure
    from copy import copy
    if parametrized: novalence = False
    # Generate random atom and parameter types
    atom_types = [AtomType(''.join(random.sample(uppercase, 3)),
                           i, random.random()*16+1, random.randint(1, 8))
                  for i in range(random.randint(8, 20))]
    bond_types = [BondType(random.random()*2, random.random()*100)
                  for i in range(random.randint(10, 20))]
    angle_types = [AngleType(random.random()*50, random.random()*120)
                   for i in range(random.randint(10, 20))]
    dihed_types = [DihedralType(random.random()*10, random.randint(1, 6),
                                random.choice([0, 180]))
                   for i in range(random.randint(10, 20))]
    rb_types = [RBTorsionType(*[random.random()*10 for i in range(6)])]
    imp_types = [ImproperType(random.random()*100, random.choice([0, 180]))
                 for i in range(random.randint(10, 20))]
    cmap_types = [CmapType(24, [random.random()*5 for i in range(24*24)])
                  for i in range(random.randint(5, 10))]
    oop_types = [OutOfPlaneBendType(random.random()*100)
                 for i in range(random.randint(10, 20))]
    strbnd_types = [StretchBendType(random.random()*10, random.random()*10,
                                    random.random()*2, random.random()*2,
                                    random.random()*120)
                    for i in range(random.randint(10, 20))]
    ang1, ang2 = list(range(-180,180,36)), list(range(-180,180,18))
    tortor_types = [TorsionTorsionType((10, 20), ang1[:], ang2[:],
                            [random.random()*10 for j in range(200)])
                    for i in range(random.randint(5, 10))]
    for typ in atom_types:
        typ.set_lj_params(random.random()*2, random.random()*2)

    struct = structure.Structure()
    # Add atoms in residues
    for res in range(random.randint(20, 30)):
        resname = ''.join(random.sample(uppercase, 3))
        resid = res + 1
        for i in range(random.randint(10, 25)):
            name = ''.join(random.sample(uppercase, 4))
            if parametrized:
                typ = random.choice(atom_types)
                type = str(typ)
                mass = typ.mass
                atomic_number = typ.atomic_number
            else:
                type = ''.join(random.sample(uppercase, 3))
                mass = random.random() * 16 + 1
                atomic_number = random.randint(1, 8)
            charge = random.random() * 2 - 1
            solvent_radius = random.random() * 2
            screen = random.random() * 2
            atom = Atom(atomic_number=atomic_number, type=type, charge=charge,
                        mass=mass, solvent_radius=solvent_radius,
                        screen=screen, name=name)
            if parametrized:
                atom.atom_type = typ
            struct.add_atom(atom, resname, resid)
    if novalence:
        return struct
    # Possibly add parameter type lists
    if parametrized:
        struct.bond_types.extend([copy(x) for x in bond_types])
        struct.bond_types.claim()
        struct.angle_types.extend([copy(x) for x in angle_types])
        struct.angle_types.claim()
        struct.dihedral_types.extend([copy(x) for x in dihed_types])
        struct.dihedral_types.claim()
        struct.rb_torsion_types.extend([copy(x) for x in rb_types])
        struct.rb_torsion_types.claim()
        struct.urey_bradley_types.extend([copy(x) for x in bond_types])
        struct.urey_bradley_types.claim()
        struct.improper_types.extend([copy(x) for x in imp_types])
        struct.improper_types.claim()
        struct.cmap_types.extend([copy(x) for x in cmap_types])
        struct.cmap_types.claim()
        struct.trigonal_angle_types.extend([copy(x) for x in angle_types])
        struct.trigonal_angle_types.claim()
        struct.out_of_plane_bend_types.extend([copy(x) for x in oop_types])
        struct.out_of_plane_bend_types.claim()
        struct.pi_torsion_types.extend([copy(x) for x in dihed_types])
        struct.pi_torsion_types.claim()
        struct.stretch_bend_types.extend([copy(x) for x in strbnd_types])
        struct.stretch_bend_types.claim()
        struct.torsion_torsion_types.extend([copy(x) for x in tortor_types])
        struct.torsion_torsion_types.claim()
        struct.adjust_types.extend([AmoebaNonbondedExceptionType(0.5, 0.5, 0.6, 0.6, 0.7)
                                    for i in range(random.randint(10, 20))])
        struct.adjust_types.claim()
    # Add valence terms with optional
    for i in range(random.randint(40, 50)):
        struct.bonds.append(Bond(*random.sample(struct.atoms, 2)))
        if parametrized:
            struct.bonds[-1].type = random.choice(struct.bond_types)
    for i in range(random.randint(35, 45)):
        struct.angles.append(Angle(*random.sample(struct.atoms, 3)))
        if parametrized:
            struct.angles[-1].type = random.choice(struct.angle_types)
    for i in range(random.randint(35, 45)):
        struct.urey_bradleys.append(UreyBradley(*random.sample(struct.atoms, 2)))
        if parametrized:
            struct.urey_bradleys[-1].type = random.choice(struct.urey_bradley_types)
    for i in range(random.randint(30, 40)):
        struct.dihedrals.append(Dihedral(*random.sample(struct.atoms, 4),
                                         improper=random.choice([True, False])))
        if parametrized:
            struct.dihedrals[-1].type = random.choice(struct.dihedral_types)
    for i in range(random.randint(30, 40)):
        struct.rb_torsions.append(Dihedral(*random.sample(struct.atoms, 4)))
        if parametrized:
            struct.rb_torsions[-1].type = random.choice(struct.rb_torsion_types)
    for i in range(random.randint(10, 20)):
        struct.impropers.append(Improper(*random.sample(struct.atoms, 4)))
        if parametrized:
            struct.impropers[-1].type = random.choice(struct.improper_types)
    for i in range(random.randint(25, 35)):
        struct.cmaps.append(Cmap(*random.sample(struct.atoms, 5)))
        if parametrized:
            struct.cmaps[-1].type = random.choice(struct.cmap_types)
    for i in range(random.randint(30, 40)):
        struct.trigonal_angles.append(TrigonalAngle(*random.sample(struct.atoms, 4)))
        if parametrized:
            struct.trigonal_angles[-1].type = random.choice(struct.trigonal_angle_types)
    for i in range(random.randint(30, 40)):
        struct.out_of_plane_bends.append(OutOfPlaneBend(*random.sample(struct.atoms, 4)))
        if parametrized:
            struct.out_of_plane_bends[-1].type = random.choice(struct.out_of_plane_bend_types)
    for i in range(random.randint(30, 40)):
        struct.stretch_bends.append(StretchBend(*random.sample(struct.atoms, 3)))
        if parametrized:
            struct.stretch_bends[-1].type = random.choice(struct.stretch_bend_types)
    for i in range(random.randint(20, 30)):
        struct.pi_torsions.append(PiTorsion(*random.sample(struct.atoms, 6)))
        if parametrized:
            struct.pi_torsions[-1].type = random.choice(struct.pi_torsion_types)
    for i in range(random.randint(10, 20)):
        struct.torsion_torsions.append(TorsionTorsion(*random.sample(struct.atoms, 5)))
        if parametrized:
            struct.torsion_torsions[-1].type = random.choice(struct.torsion_torsion_types)
    # Now use some lesser-used features
    for i in range(random.randint(5, 10)):
        struct.acceptors.append(AcceptorDonor(*random.sample(struct.atoms, 2)))
        struct.donors.append(AcceptorDonor(*random.sample(struct.atoms, 2)))
        struct.groups.append(Group(random.choice(struct.atoms), 2, 0))
        struct.chiral_frames.append(ChiralFrame(*random.sample(struct.atoms, 2),
                                                chirality=random.choice([-1, 1])))
        struct.multipole_frames.append(MultipoleFrame(random.choice(struct.atoms),
                                                      0, 1, 2, 3))
    for i in range(random.randint(20, 30)):
        struct.adjusts.append(NonbondedException(*random.sample(struct.atoms, 2)))
        if parametrized:
            struct.adjusts[-1].type = random.choice(struct.adjust_types)
    struct.prune_empty_terms()
    struct.unchange()
    struct.update_dihedral_exclusions()
    return struct

def equal_atoms(tester, a1, a2):
    """ Tests equality of two atoms based on properties

    Parameters
    ----------
    tester : unittest.TestCase
        TestCase instance
    a1 : Atom
        First atom to compare
    a2 : Atom
        Second atom to compare
    """
    tester.assertEqual(a1.atomic_number, a2.atomic_number)
    tester.assertEqual(a1.screen, a2.screen)
    tester.assertEqual(a1.name, a2.name)
    tester.assertEqual(a1.type, a2.type)
    tester.assertEqual(a1.atom_type, a2.atom_type)
    tester.assertEqual(a1.charge, a2.charge)
    tester.assertEqual(a1.mass, a2.mass)
    tester.assertEqual(a1.nb_idx, a2.nb_idx)
    tester.assertEqual(a1.solvent_radius, a2.solvent_radius)
    tester.assertEqual(a1.tree, a2.tree)
    tester.assertEqual(a1.join, a2.join)
    tester.assertEqual(a1.irotat, a2.irotat)
    tester.assertEqual(a1.occupancy, a2.occupancy)
    tester.assertEqual(a1.bfactor, a2.bfactor)
    tester.assertEqual(a1.rmin, a2.rmin)
    tester.assertEqual(a1.epsilon, a2.epsilon)
    tester.assertEqual(a1.rmin_14, a2.rmin_14)
    tester.assertEqual(a1.epsilon_14, a2.epsilon_14)
    for key in ('xx', 'xy', 'xz', 'vx', 'vy', 'vz', 'multipoles',
                'type_idx', 'class_idx', 'polarizability', 'vdw_weight'):
        if hasattr(a2, key):
            if isinstance(getattr(a2, key), np.ndarray):
                np.testing.assert_equal(
                        getattr(a1, key), getattr(a2, key)
                )
            else:
                tester.assertEqual(getattr(a1, key), getattr(a2, key))
        else:
            tester.assertFalse(hasattr(a1, key))

def is_jenkins():
    return 'JENKINS_URL' in os.environ

def has_old_vec3():
    from parmed.vec3 import Vec3
    v1 = Vec3(1, 2, 3)
    if not hasattr(v1, 'x') or v1.x != 1:
        return True
    try:
        v2 = -v1
    except TypeError:
        return True
    return False
