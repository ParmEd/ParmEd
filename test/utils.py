"""
Useful functions for the test cases
"""
from parmed.utils.six import string_types
from parmed.utils.six.moves import zip
import os
from os.path import join, split, abspath
import random
import unittest
import warnings
warnings.filterwarnings('error', category=DeprecationWarning)

try:
    from simtk import openmm
    openmm_version = tuple([int(x) for x in openmm.__version__.split('.')])
except ImportError:
    openmm_version = None

try:
    from string import uppercase
except ImportError:
    from string import ascii_uppercase as uppercase

def skip_big_tests():
    return os.getenv('PARMED_SKIP_BIG_TESTS') is not None

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
                raise self.failureException(
                            '%s != %s with relative tolerance %g (%f)' %
                            (val1, val2, 10**-places, ratio)
                )
            else:
                if abs(ratio - 1) < delta:
                    return
                raise self.failureException(
                            '%s != %s with relative tolerance %g (%f)' %
                            (val1, val2, delta, ratio))

class FileIOTestCase(unittest.TestCase):

    def setUp(self):
        try:
            os.makedirs(get_fn('writes'))
        except OSError:
            pass

    def tearDown(self):
        self._empty_writes()
        try:
            os.rmdir(get_fn('writes'))
        except OSError:
            pass

    def _empty_writes(self):
        """ Empty the "writes" directory """
        try:
            for f in os.listdir(get_fn('writes')):
                os.unlink(get_fn(f, written=True))
        except OSError:
            pass

def get_fn(filename, written=False):
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
    if written:
        return join(split(abspath(__file__))[0], 'files', 'writes', filename)
    else:
        return join(split(abspath(__file__))[0], 'files', filename)

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
    return join(split(abspath(__file__))[0], 'files', 'saved', filename)

def has_scipy():
    try:
        import scipy.io.netcdf as nc
        return True
    except ImportError:
        return False

def has_netcdf4():
    try:
        import netCDF4
        return True
    except ImportError:
        return False

def has_scientific():
    try:
        from Scientific.IO.NetCDF import NetCDFFile
        return True
    except ImportError:
        return False

def has_pynetcdf():
    try:
        import pynetcdf
        return True
    except ImportError:
        return False

def has_numpy():
    return True

def diff_files(file1, file2, ignore_whitespace=True,
               absolute_error=None, relative_error=None,
               comment=None):
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
                    if not detailed_diff(l1,l2,absolute_error,relative_error):
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
                    if not detailed_diff(l1,l2,absolute_error,relative_error):
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

def detailed_diff(l1, l2, absolute_error=None, relative_error=None):
    """
    Check individual fields to make sure numbers are numerically equal if the
    lines differ. Also ignore fields that appear to be a file name, since those
    will be system-dependent
    """
    fdir = os.path.split(get_fn('writes'))[0]
    w1 = l1.split()
    w2 = l2.split()
    if len(w1) != len(w2): return False
    for wx, wy in zip(w1, w2):
        try:
            wx = float(wx)
            wy = float(wy)
        except ValueError:
            if isinstance(wx, float) or (wx != wy and not
                    (wx.startswith(fdir) or wy.startswith(fdir))):
                return False
        else:
            if wx != wy:
                if absolute_error is not None and abs(wx-wy) > absolute_error:
                    return False
                elif relative_error is not None:
                    if wx == 0 or wy == 0 and abs(wx-wy) > relative_error:
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
            radii = random.random() * 2
            screen = random.random() * 2
            atom = Atom(atomic_number=atomic_number, type=type,
                        charge=charge, mass=mass, radii=radii,
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
        struct.groups.append(Group(*random.sample(range(1, 11), 3)))
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
