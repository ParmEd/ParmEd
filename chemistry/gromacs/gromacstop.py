"""
This module contains functionality relevant to loading a GROMACS topology file
and building a Structure from it
"""
from chemistry.exceptions import GromacsTopologyError, GromacsTopologyWarning
from chemistry.formats.registry import FileFormatType
from chemistry import gromacs as gmx, ParameterSet
from chemistry.gromacs._gromacsfile import GromacsFile
from chemistry.structure import Structure
from chemistry.topologyobjects import (Atom, Bond, Angle, Dihedral, Improper,
            NonbondedException, ExtraPoint, BondType)
from chemistry.periodic_table import element_by_mass, AtomicNum
from chemistry import unit as u
from chemistry.utils.io import genopen
from chemistry.utils.six import add_metaclass
from chemistry.utils.six.moves import range
from contextlib import closing
import warnings

# Gromacs uses "funct" flags in its parameter files to indicate what kind of
# functional form is used for each of its different parameter types. This is
# taken from the topdirs.c source code file along with a table in the Gromacs
# user manual. The table below summarizes my findings, for reference:

# Bonds
# -----
#  1 - F_BONDS : simple harmonic potential
#  2 - F_G96BONDS : fourth-power potential
#  3 - F_MORSE : morse potential
#  4 - F_CUBICBONDS : cubic potential
#  5 - F_CONNBONDS : not even implemented in GROMACS
#  6 - F_HARMONIC : seems to be the same as (1) ??
#  7 - F_FENEBONDS : finietely-extensible-nonlinear-elastic (FENE) potential
#  8 - F_TABBONDS : bond function from tabulated function
#  9 - F_TABBONDSNC : bond function from tabulated function (no exclusions)
# 10 - F_RESTRBONDS : restraint bonds

# Angles
# ------
#  1 - F_ANGLES : simple harmonic potential
#  2 - F_G96ANGLES : cosine-based angle potential
#  3 - F_CROSS_BOND_BONDS : bond-bond cross term potential
#  4 - F_CROSS_BOND_ANGLES : bond-angle cross term potential
#  5 - F_UREY_BRADLEY : Urey-Bradley angle-bond potential
#  6 - F_QUARTIC_ANGLES : 4th-order polynomial potential
#  7 - F_TABANGLES : angle function from tabulated function
#  8 - F_LINEAR_ANGLES : angle function from tabulated function
#  9 - F_RESTRANGLES : restricted bending potential

# Dihedrals
# ---------
#  1 - F_PDIHS : periodic proper torsion potential [ k(1+cos(n*phi-phase)) ]
#  2 - F_IDIHS : harmonic improper torsion potential
#  3 - F_RBDIHS : Ryckaert-Bellemans torsion potential
#  4 - F_PIDIHS : periodic harmonic improper torsion potential (same as 1)
#  5 - F_FOURDIHS : Fourier dihedral torsion potential
#  8 - F_TABDIHS : dihedral potential from tabulated function
#  9 - F_PDIHS : Same as 1, but can be multi-term
# 10 - F_RESTRDIHS : Restricted torsion potential
# 11 - F_CBTDIHS : combined bending-torsion potential

@add_metaclass(FileFormatType)
class GromacsTopologyFile(Structure):
    """
    Loads a GROMACS topology file

    Parameters
    ----------
    fname : str=None
        Name of the GROMACS topology file to parse, if any

    Attributes
    ----------
    parameterset : ParameterSet
        The set of parameters defining a force field
    """
    #===================================================

    @staticmethod
    def id_format(filename):
        """ Identifies the file as a GROMACS topology file

        Parameters
        ----------
        filename : str
            Name of the file to check if it is a gromacs topology file

        Returns
        -------
        is_fmt : bool
            If it is identified as a gromacs topology, return True. False
            otherwise
        """
        with closing(genopen(filename)) as f:
            for line in f:
                if line.startswith(';'): continue
                if line.startswith('#'):
                    if line.startswith('#if'): continue
                    if line.startswith('#define'): continue
                    if line.startswith('#include'): continue
                    if line.startswith('#undef'): continue
                    if line.startswith('#endif'): continue
                    return False
                if line.strip() == '[ moleculetype ]': return True
                return False
            return False

    #===================================================

    def __init__(self, fname=None, defines=None):
        super(GromacsTopologyFile, self).__init__()
        self.name = fname
        self.parameterset = ParameterSet()
        # This protects us agaself using topologies defining a functional form
        # that I have no idea how to deal with
        self.unknown_functional = False
        if fname is not None:
            self.rdparm(fname, defines)

    #===================================================

    def rdparm(self, fname, defines=None):
        """
        Reads the GROMACS topology file

        Parameters
        ----------
        fname : str
            The name of the file to read
        defines : list of str=None
            If specified, this is the set of defines to use when parsing the
            topology file
        """
        params = self.parameterset
        molecules = dict()
        structure_contents = []
        with closing(GromacsFile(fname, includes=[gmx.GROMACS_TOPDIR],
                                 defines=defines)) as f:
            current_section = None
            current_define = None
            for line in f:
                line = line.strip()
                if not line: continue

                if line[0] == '[':
                    current_section = line[1:-1].strip()
                elif current_section == 'moleculetype':
                    molname, nrexcl = line.split()
                    nrexcl = int(nrexcl)
                    if molname in molecules:
                        raise GromacsTopologyError('Duplicate definition of '
                                                   'molecule %s' % molname)
                    molecule = Structure()
                    molecules[molname] = (molecule, nrexcl)
                elif current_section == 'atoms':
                    words = line.split()
                    if len(words) < 8:
                        mass = -1
                        atomic_number = -1
                    else:
                        mass = float(words[7])
                        atomic_number = AtomicNum[element_by_mass(mass)]
                    if len(words) < 7:
                        charge = None
                    else:
                        charge = float(words[6])
                    if atomic_number == 0:
                        atom = ExtraPoint(name=words[4], type=words[1],
                                          charge=charge)
                    else:
                        atom = Atom(atomic_number=atomic_number, name=words[4],
                                    type=words[1], charge=charge)
                    molecule.add_atom(atom, words[3], int(words[2]))
                elif current_section == 'bonds':
                    words = line.split()
                    i, j = int(words[0]), int(words[1])
                    if words[2] != '1':
                        warnings.warn('bond funct != 1; unknown functional',
                                      GromacsTopologyWarning)
                        self.unknown_functional = True
                    molecule.bonds.append(Bond(molecule.atoms[i-1],
                                               molecule.atoms[j-1]))
                elif current_section == 'pairs':
                    words = line.split()
                    i, j = int(words[0]), int(words[1])
                    if words[2] != '1':
                        # This is not even supported in Gromacs
                        warnings.warn('pairs funct != 1; unknown functional',
                                      GromacsTopologyWarning)
                        self.unknown_functional = True
                    molecule.adjusts.append(
                            NonbondedException(molecule.atoms[i-1],
                                               molecule.atoms[j-1])
                    )
                elif current_section == 'angles':
                    words = line.split()
                    i, j, k = int(words[0]), int(words[1]), int(words[2])
                    if words[3] != '1':
                        warnings.warn('angles funct != 1; unknown functional',
                                      GromacsTopologyWarning)
                        self.unknown_functional = True
                    molecule.angles.append(
                            Angle(molecule.atoms[i-1], molecule.atoms[j-1],
                                  molecule.atoms[k-1])
                    )
                elif current_section == 'dihedrals':
                    words = line.split()
                    i, j, k, l = [int(x) for x in words[:4]]
                    if words[4] in ('1', '9', '4'):
                        # Normal dihedral
                        dih = Dihedral(molecule.atoms[i-1], molecule.atoms[j-1],
                                       molecule.atoms[k-1], molecule.atoms[l-1])
                        molecule.dihedrals.append(dih)
                    elif words[4] == '2':
                        # Improper
                        imp = Improper(molecule.atoms[i-1], molecule.atoms[j-1],
                                       molecule.atoms[k-1], molecule.atoms[l-1])
                        self.impropers.append(imp)
                    else:
                        # ??? unknown
                        warnings.warn('torsions funct != 1, 2, 4, or 9; unknown'
                                      ' functional', GromacsTopologyWarning)
                        dih = Dihedral(molecule.atoms[i-1], molecule.atoms[j-1],
                                       molecule.atoms[k-1], molecule.atoms[l-1])
                        self.dihedrals.append(dih)
                elif current_section == 'system':
                    self.title = line
                elif current_section == 'defaults':
                    words = line.split()
                    if len(words) < 4:
                        raise GromacsTopologyError('Too few fields in '
                                                   '[ defaults ]')
                    if words[0] != '1':
                        warnings.warn('Unsupported nonbonded type; unknown '
                                      'functional', GromacsTopologyWarning)
                        self.unknown_functional = True
                    if words[1] != '2':
                        warnings.warn('Unsupported combining rule',
                                      GromacsTopologyWarning)
                        self.unknown_functional = True
                    if words[2].lower() == 'no':
                        warnings.warn('gen_pairs=no is not supported')
                        self.unknown_functional = True
                    self._fudgeLJ = float(words[3])
                    self._fudgeQQ = float(words[4])
                elif current_section == 'molecules':
                    name, num = line.split()
                    num = int(num)
                    structure_contents.append((name, num))
                elif current_section == 'settles':
                    # Instead of adding bonds that get constrained for waters
                    # (or other 3-atom molecules), GROMACS uses a "settles"
                    # section to specify the constraint geometry. We have to
                    # translate that into bonds.
                    natoms = len([a for a in molecule.atoms
                                    if not isinstance(a, ExtraPoint)])
                    if natoms != 3:
                        raise GromacsTopologyError("Cannot SETTLE a %d-atom "
                                                   "molecule" % natoms)
                    try:
                        oxy, = [atom for atom in molecule.atoms
                                    if atom.atomic_number == 8]
                        hyd1, hyd2 = [atom for atom in molecule.atoms
                                        if atom.atomic_number == 1]
                    except ValueError:
                        raise GromacsTopologyError("Can only SETTLE water; "
                                    "Could not detect 2 hydrogens and 1 oxygen")
                    #TODO see if there's a bond_type entry in the parameter set
                    #     that we can fill in
                    try:
                        i, funct, doh, dhh = line.split()
                        doh, dhh = float(doh), float(dhh)
                    except ValueError:
                        raise GromacsTopologyError('Bad [ settles ] line')
                    bt_oh = BondType(1000*u.kilojoules_per_mole/u.nanometers**2,
                                     doh*u.nanometers, list=molecule.bond_types)
                    bt_hh = BondType(1000*u.kilojoules_per_mole/u.nanometers**2,
                                     dhh*u.nanometers, list=molecule.bond_types)
                    molecule.bond_types.extend([bt_oh, bt_hh])
                    molecule.bonds.append(Bond(oxy, hyd1, bt_oh))
                    molecule.bonds.append(Bond(oxy, hyd2, bt_oh))
                    molecule.bonds.append(Bond(hyd1, hyd2, bt_hh))
                elif current_section in ('virtual_sites3', 'dummies3'):
                    words = line.split()
                    vsite = molecule.atoms[int(words[0])-1]
                    atoms = [molecule.atoms[int(i)-1] for i in words[1:4]]
                    funct = int(words[4])
                    if funct == 1:
                        a, b = float(words[5]), float(words[6])
                        if abs(a - b) > TINY:
                            raise GromacsTopologyError('Cannot handle virtual '
                                    'site frames with different weights')
                    else:
                        raise GromacsTopologyError('Only 3-point virtual site '
                                                   'type "1" is supported')
                    parent = atoms[0]
                    bondlen = ThreeParticleVirtualSite.from_weights(parent,
                            atoms[1], atoms[2], a, b)
                    bt_vs = BondType(0, bondlen*u.nanometers,
                                     list=molecule.bond_types)
                    if vsite in parent.bond_partners:
                        raise GromacsTopologyError('Unexpected bond b/w '
                                    'virtual site and its parent')
                    molecule.bonds.append(Bond(vsite, parent, bt_vs))
                elif current_section == 'exclusions':
                    atoms = [molecule.atoms[int(w)-1] for w in line.split()]
                    for a in atoms[1:]:
                        atom.exclude(a)
#               elif current_section == 'bondtypes':
#                   words = line.split()
#                   r = float(words[3]) * u.nanometers
#                   k = (float(words[4]) / 2) * (
#                           u.kilojoules_per_mole / u.nanometers**2)
#                   if words[2] != '1':
#                       warnings.warn('bondtypes funct != 1; unknown '
#                                     'functional', GromacsTopologyWarning)
#                       unst.unknown_functional = True
#                   ptype = BondType(k, r)
#                   params.bond_types[(words[0], words[1])] = ptype
#                   params.bond_types[(words[1], words[0])] = ptype
#               elif current_section == 'angletypes':
#                   words = line.split()
#                   theta = float(words[4]) * u.degrees
#                   k = (float(words[5]) / 2) * (
#                           u.kilojoules_per_mole / u.radians**2)
#                   if words[2] != '1' and words[2] != '5':
#                       warnings.warn('angletypes funct != 1; unknown '
#                                     'functional', GromacsTopologyWarning)
#                       self.unknown_functional = True
#                   if words[2] == '5':
#                       # Contains the angle with urey-bradley
#                       ub0 = float(words[6])
#                       cub = float(words[7])
#                       if cub == 0:
#                           ub = NoUreyBradley
#                       else:
#                           ub0 *= u.nanometers
#                           cub *= u.kilojoules_per_mole / u.nanometers**2
#                           ub = BondType(cub, ub0)
#                           params.urey_bradley_types[(words[0], words[2])] = ub
#                           params.urey_bradley_types[(words[2], words[0])] = ub
#                   ptype = AngleType(k, theta)
#                   params.angle_types[(words[0], words[1], words[2])] = ptype
#                   params.angle_types[(words[2], words[1], words[0])] = ptype
#               elif current_section == 'dihedraltypes':
#                   words = line.split()
#                   replace = False
#                   dtype = 'normal'
#                   a1, a2, a3, a4 = words[:4]
#                   if words[4] == '1':
#                       pass
#                   if words[4] == '4':
#                       replace = True
#                   elif words[4] == '9':
#                       pass
#                   elif words[4] == '2':
#                       replace = True
#                       dtype = 'improper'
#                   elif words[4] == '5':
#                       dtype = 'rbtorsion'
#                   else:
#                       warnings.warn('dihedraltypes funct not supported',
#                                     GromacsTopologyWarning)
#                       self.unknown_functional = True
#                   # Do the proper types
#                   if dtype == 'normal':
#                       phase = float(words[5]) * u.degrees
#                       phi_k = float(words[6]) * u.kilojoules_per_mole
#                       per = int(words[7])
#                       dt = DihedralType(phi_k, per, phase)
#                       key = (words[0], words[1], words[2], words[3])
#                       rkey = (words[3], words[2], words[1], words[0])
#                       if replace or not key in params.dihedral_types:
#                           dtl = DihedralTypeList()
#                           dtl.append(dt)
#                           params.dihedral_types[key] = dtl
#                           params.dihedral_types[rkey] = dtl
#                       else:
#                           params.dihedral_types[key].append(dt)
#                   elif dtype == 'improper':
#                       theta = float(words[5])*u.degrees
#                       k = float(words[6])*u.kilojoules_per_mole/u.radians**2
#                       a1, a2, a3, a4 = words[:4]
#                       ptype = ImproperType(k, theta)
#                       params.improper_types[(a1, a2, a3, a4)] = ptype
#                   elif dtype == 'rbtorsion':
#                       a1, a2, a3, a4 = words[:4]
#                       c0, c1, c2, c3, c4, c5 = [float(x) for x in words[5:11]]
#                       ptype = RBTorsionType(c0, c1, c2, c3, c4, c5)
#                       params.rb_torsion_types[(a1, a2, a3, a4)] = ptype
#                       params.rb_torsion_types[(a4, a3, a2, a1)] = ptype
        # TODO: What should come first? Combining, or parametrization?
        for molname, num in structure_contents:
            if molname not in molecules:
                raise GromacsTopologyError('Structure contains %s molecules, '
                                           'but no template defined' % molname)
            molecule, nrexcl = molecules[molname]
            if nrexcl < 3 and _any_atoms_farther_than(molecule, nrexcl):
                warnings.warn('nrexcl %d not currently supported' % nrexcl,
                              GromacsTopologyWarning)
            elif nrexcl > 3 and _any_atoms_farther_than(molecule, 3):
                warnings.warn('nrexcl %d not currently supported' % nrexcl,
                              GromacsTopologyWarning)
            if num == 0:
                warnings.warn('Detected addition of 0 %s molecules in topology '
                              'file' % molname, GromacsTopologyWarning)
            if num == 1:
                self += molecules[molname][0]
            elif num > 1:
                self += molecules[molname][0] * num
            else:
                raise GromacsTopologyError('Cannot add %d %s molecules' %
                                           (num, molname))
        self.molecules = molecules

def _any_atoms_farther_than(structure, limit=3):
    """
    This function checks to see if there are any atom pairs farther away in the
    bond graph than the desired limit

    Parameters
    ----------
    structure : :class:`Structure`
        The structure to search through
    limit : int, optional
        The most number of bonds away to check for. Default is 3

    Returns
    -------
    within : bool
        True if any atoms are *more* than ``limit`` bonds away from any other
        atom
    """
    import sys
    if len(structure.atoms) <= limit + 1: return False
    sys.setrecursionlimit(max(sys.getrecursionlimit(), limit+1))
    for atom in structure.atoms:
        for atom in structure.atoms: atom.marked = limit + 1
        _mark_graph(atom, 0)
        if any((atom.marked > limit for atom in structure.atoms)):
            return True
    return False

def _mark_graph(atom, num):
    """ Marks all atoms in the graph listing the minimum number of bonds each
    atom is away from the current atom

    Parameters
    ----------
    atom : :class:`Atom`
        The current atom to evaluate in the bond graph
    num : int
        The current depth in our search
    limit : int
        The maximum depth we want to search
    """
    atom.marked = num
    for a in atom.bond_partners:
        if a.marked <= num: continue
        _mark_graph(a, num+1)
