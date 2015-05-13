"""
This module contains functionality relevant to loading a GROMACS topology file
and building a Structure from it
"""
from __future__ import print_function, division, absolute_import

from chemistry.constants import TINY
from chemistry.exceptions import GromacsTopologyError, GromacsTopologyWarning
from chemistry.formats.registry import FileFormatType
from chemistry.parameters import ParameterSet
from chemistry.gromacs._gromacsfile import GromacsFile
from chemistry.structure import Structure
from chemistry.topologyobjects import (Atom, Bond, Angle, Dihedral, Improper,
            NonbondedException, ExtraPoint, BondType, Cmap, NoUreyBradley,
            AngleType, DihedralType, DihedralTypeList, ImproperType,
            RBTorsionType, ThreeParticleExtraPointFrame, AtomType,
            TwoParticleExtraPointFrame, OutOfPlaneExtraPointFrame)
from chemistry.periodic_table import element_by_mass, AtomicNum
from chemistry import unit as u
from chemistry.utils.io import genopen
from chemistry.utils.six import add_metaclass, string_types
from chemistry.utils.six.moves import range
from contextlib import closing
from datetime import datetime
import os
import re
import sys
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

_sectionre = re.compile(r'\[ (\w+) \]\s*$')

class _Defaults(object):
    """ Global properties of force fields as implemented in GROMACS """
    def __init__(self, nbfunc=1, comb_rule=2, gen_pairs='yes',
                 fudgeLJ=0.5, fudgeQQ=0.83333333):
        if int(nbfunc) not in (1, 2):
            raise ValueError('nbfunc must be 1 (L-J) or 2 (Buckingham)')
        if int(comb_rule) not in (1, 2, 3):
            raise ValueError('comb_rule must be 1, 2, or 3')
        if gen_pairs not in ('yes', 'no'):
            raise ValueError("gen_pairs must be 'yes' or 'no'")
        if float(fudgeLJ) < 0:
            raise ValueError('fudgeLJ must be non-negative')
        if float(fudgeQQ) < 0:
            raise ValueError('fudgeQQ must be non-negative')
        self.nbfunc = int(nbfunc)
        self.comb_rule = int(comb_rule)
        self.gen_pairs = gen_pairs
        self.fudgeLJ = float(fudgeLJ)
        self.fudgeQQ = float(fudgeQQ)

    def __getitem__(self, idx):
        # Treat it like the array that it is in the topology file
        if idx < 0: idx += 5
        if idx == 0: return self.nbfunc
        if idx == 1: return self.comb_rule
        if idx == 2: return self.gen_pairs
        if idx == 3: return self.fudgeLJ
        if idx == 4: return self.fudgeQQ
        raise IndexError('Index %d out of range' % idx)

    def __setitem__(self, idx, value):
        if idx < 0: idx += 5
        if idx == 0:
            if int(value) not in (1, 2):
                raise ValueError('nbfunc must be 1 or 2')
            self.nbfunc = int(value)
        elif idx == 1:
            if int(value) not in (1, 2, 3):
                raise ValueError('comb_rule must be 1, 2, or 3')
            self.comb_rule = int(value)
        elif idx == 2:
            if value not in ('yes', 'no'):
                raise ValueError('gen_pairs must be "yes" or "no"')
            self.gen_pairs = value
        elif idx == 3:
            if float(value) < 0:
                raise ValueError('fudgeLJ must be non-negative')
            self.fudgeLJ = float(value)
        elif idx == 4:
            if float(value) < 0:
                raise ValueError('fudgeQQ must be non-negative')
            self.fudgeLJ = value

@add_metaclass(FileFormatType)
class GromacsTopologyFile(Structure):
    """ Class providing a parser and writer for a GROMACS topology file

    Parameters
    ----------
    fname : str
        The name of the file to read
    defines : list of str=None
        If specified, this is the set of defines to use when parsing the
        topology file
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
                if line.startswith(';'):
                    line = line[:line.index(';')]
                if not line.strip(): continue
                if line.startswith('#'):
                    if line.startswith('#if'): continue
                    if line.startswith('#define'): continue
                    if line.startswith('#include'): continue
                    if line.startswith('#undef'): continue
                    if line.startswith('#endif'): continue
                    return False
                rematch = _sectionre.match(line)
                if not rematch:
                    return False
                sec, = rematch.groups()
                return sec in ('atoms', 'atomtypes', 'defaults', 'moleculetype')
            return False

    #===================================================

    def __init__(self, fname=None, defines=None):
        super(GromacsTopologyFile, self).__init__()
        self.parameterset = ParameterSet()
        if fname is not None:
            self.read(fname, defines)
        self.defaults = _Defaults()

    #===================================================

    def read(self, fname, defines=None):
        """ Reads the topology file into the current instance """
        from chemistry import gromacs as gmx
        params = self.parameterset
        molecules = dict()
        structure_contents = []
        with closing(GromacsFile(fname, includes=[gmx.GROMACS_TOPDIR],
                                 defines=defines)) as f:
            current_section = None
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
                    molecule.nrexcl = nrexcl
                elif current_section == 'atoms':
                    words = line.split()
                    try:
                        attype = params.atom_types[words[1]]
                    except KeyError:
                        attype = None
                    if len(words) < 8:
                        if attype is not None:
                            mass = attype.mass
                            atomic_number = attype.atomic_number
                        else:
                            mass = -1
                            atomic_number = -1
                    else:
                        mass = float(words[7])
                        if attype is not None and attype.atomic_number >= 0:
                            atomic_number = attype.atomic_number
                        else:
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
                                    type=words[1], charge=charge, mass=mass)
                    molecule.add_atom(atom, words[3], int(words[2]))
                elif current_section == 'bonds':
                    words = line.split()
                    i, j = int(words[0])-1, int(words[1])-1
                    funct = int(words[2])
                    if funct != 1:
                        warnings.warn('bond funct != 1; unknown functional',
                                      GromacsTopologyWarning)
                        self.unknown_functional = True
                    molecule.bonds.append(Bond(molecule.atoms[i],
                                               molecule.atoms[j]))
                    molecule.bonds[-1].funct = funct
                elif current_section == 'pairs':
                    words = line.split()
                    i, j = int(words[0])-1, int(words[1])-1
                    funct = int(words[2])
                    if funct != 1:
                        # This is not even supported in Gromacs
                        warnings.warn('pairs funct != 1; unknown functional',
                                      GromacsTopologyWarning)
                        self.unknown_functional = True
                    molecule.adjusts.append(
                            NonbondedException(molecule.atoms[i],
                                               molecule.atoms[j])
                    )
                    molecule.adjusts[-1].funct = funct
                elif current_section == 'angles':
                    words = line.split()
                    i, j, k = [int(w)-1 for w in words[:3]]
                    funct = int(words[3])
                    if funct not in (1, 5):
                        warnings.warn('angles funct != 1 or 5; unknown '
                                      'functional', GromacsTopologyWarning)
                        self.unknown_functional = True
                    molecule.angles.append(
                            Angle(molecule.atoms[i], molecule.atoms[j],
                                  molecule.atoms[k])
                    )
                    molecule.angles[-1].funct = funct
                elif current_section == 'dihedrals':
                    words = line.split()
                    i, j, k, l = [int(x)-1 for x in words[:4]]
                    funct = int(words[4])
                    if funct in (1, 9, 4):
                        # Normal dihedral
                        dih = Dihedral(molecule.atoms[i], molecule.atoms[j],
                                       molecule.atoms[k], molecule.atoms[l])
                        molecule.dihedrals.append(dih)
                        molecule.dihedrals[-1].funct = funct
                    elif funct == 2:
                        # Improper
                        imp = Improper(molecule.atoms[i], molecule.atoms[j],
                                       molecule.atoms[k], molecule.atoms[l])
                        molecule.impropers.append(imp)
                        molecule.impropers[-1].funct = funct
                    else:
                        # ??? unknown
                        warnings.warn('torsions funct != 1, 2, 4, or 9; unknown'
                                      ' functional', GromacsTopologyWarning)
                        dih = Dihedral(molecule.atoms[i], molecule.atoms[j],
                                       molecule.atoms[k], molecule.atoms[l])
                        molecule.dihedrals.append(dih)
                        molecule.dihedrals[-1].funct == funct
                elif current_section == 'cmap':
                    words = line.split()
                    i, j, k, l, m = [int(w)-1 for w in words[:5]]
                    funct = int(words[5])
                    if funct != 1:
                        warnings.warn('cmap funct != 1; unknown functional',
                                      GromacsTopologyWarning)
                    cmap = Cmap(molecule.atoms[i], molecule.atoms[j],
                                molecule.atoms[k], molecule.atoms[l],
                                molecule.atoms[m])
                    molecule.cmaps.append(cmap)
                    molecule.cmaps[-1].funct = funct
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
                        warnings.warn('gen_pairs=no is not supported',
                                      GromacsTopologyWarning)
                        self.unknown_functional = True
                    self.defaults = _Defaults(*words)
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
                    bt_oh = BondType(5e5*u.kilojoules_per_mole/u.nanometers**2,
                                     doh*u.nanometers, list=molecule.bond_types)
                    bt_hh = BondType(5e5*u.kilojoules_per_mole/u.nanometers**2,
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
                    bondlen = ThreeParticleExtraPointFrame.from_weights(parent,
                            atoms[1], atoms[2], a, b)
                    bt_vs = BondType(0, bondlen*u.angstroms,
                                     list=molecule.bond_types)
                    if vsite in parent.bond_partners:
                        raise GromacsTopologyError('Unexpected bond b/w '
                                    'virtual site and its parent')
                    molecule.bonds.append(Bond(vsite, parent, bt_vs))
                    molecule.bond_types.append(bt_vs)
                elif current_section == 'exclusions':
                    atoms = [molecule.atoms[int(w)-1] for w in line.split()]
                    for a in atoms[1:]:
                        atoms[0].exclude(a)
                elif current_section == 'atomtypes':
                    words = line.split()
                    attype = words[0]
                    try:
                        atnum = int(words[1])
                    except ValueError:
                        atnum = -1
                    mass = float(words[2])
                    if mass > 0 and atnum == -1:
                        atnum = AtomicNum[element_by_mass(mass)]
                    chg = float(words[3])
                    ptype = words[4]
                    sig = float(words[5]) * u.nanometers
                    eps = float(words[6]) * u.kilojoules_per_mole
                    typ = AtomType(attype, None, mass, atnum)
                    typ.set_lj_params(eps, sig*2**(1/6))
                    params.atom_types[attype] = typ
#               elif current_section == 'bondtypes':
#                   words = line.split()
#                   r = float(words[3]) * u.nanometers
#                   k = (float(words[4]) / 2) * (
#                           u.kilojoules_per_mole / u.nanometers**2)
#                   if words[2] != '1':
#                       warnings.warn('bondtypes funct != 1; unknown '
#                                     'functional', GromacsTopologyWarning)
#                       self.unknown_functional = True
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
#                       c0, c1, c2, c3, c4, c5 = [float(x)*u.kilojoules_per_mole
#                                                   for x in words[5:11]]
#                       ptype = RBTorsionType(c0, c1, c2, c3, c4, c5)
#                       params.rb_torsion_types[(a1, a2, a3, a4)] = ptype
#                       params.rb_torsion_types[(a4, a3, a2, a1)] = ptype
            itplist = f.included_files

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
        self.itps = itplist

    #===================================================

    @classmethod
    def from_structure(cls, struct, copy=False):
        """ Instantiates a GromacsTopologyFile instance from a Structure

        Parameters
        ----------
        struct : :class:`chemistry.Structure`
            The input structure to generate from
        copy : bool, optional
            If True, assign from a *copy* of ``struct`` (this is a lot slower).
            Default is False

        Returns
        -------
        gmxtop : :class:`GromacsTopologyFile`
            The topology file defined by the given struct
        """
        from copy import copy as _copy
        gmxtop = cls()
        if copy:
            struct = _copy(struct)
        gmxtop.atoms = struct.atoms
        gmxtop.bonds = struct.bonds
        gmxtop.angles = struct.angles
        gmxtop.dihedrals = struct.dihedrals
        gmxtop.impropers = struct.impropers
        gmxtop.cmaps = struct.cmaps
        gmxtop.rb_torsions = struct.rb_torsions
        gmxtop.urey_bradleys = struct.urey_bradleys
        gmxtop.adjusts = struct.adjusts
        gmxtop.bond_types = struct.bond_types
        gmxtop.angle_types = struct.angle_types
        gmxtop.dihedral_types = struct.dihedral_types
        gmxtop.improper_types = struct.improper_types
        gmxtop.cmap_types = struct.cmap_types
        gmxtop.rb_torsion_types = struct.rb_torsion_types
        if (struct.trigonal_angles or
                struct.out_of_plane_bends or
                struct.pi_torsions or
                struct.stretch_bends or
                struct.torsion_torsions or
                struct.chiral_frames or
                struct.multipole_frames):
            raise TypeError('GromacsTopologyFile does not support Amoeba '
                            'potential terms')
        return gmxtop

    #===================================================

    def write(self, dest, include_itps=None, combine=None):
        """ Write a Gromacs Topology File from a Structure

        Parameters
        ----------
        dest : str or file-like
            The name of a file or a file object to write the Gromacs topology to
        include_itps : list of str
            This keyword-only parameter contains a list of ITP files to include
            in the beginning of the generated Gromacs topology file
        combine : 'all', None, or list of iterables
            If None, no molecules are combined into a single moleculetype. If
            'all', all molecules are combined into a single moleculetype.
            Otherwise, the list of molecule indices will control which atoms are
            combined into single moleculetype's. Each index can appear *only*
            once, and start from 0. The same molecule number cannot appear in
            two lists

        Raises
        ------
        ValueError if the same molecule number appears in multiple combine lists
        """
        import chemistry.gromacs as gmx
        from chemistry import __version__
        own_handle = False
        fname = ''
        if isinstance(dest, string_types):
            fname = '%s ' % dest
            dest = genopen(dest, 'w')
            own_handle = True
        elif not hasattr(dest, 'write'):
            raise TypeError('dest must be a file name or file-like object')

        try:
            # Write the header
            now = datetime.now()
            dest.write('''\
;
;   File %swas generated
;   By user: %s (%d)
;   On host: %s
;   At date: %s
;
;   This is a standalone topology file
;
;   Created by:
;   ParmEd:       %s, VERSION %s
;   Executable:   %s
;   Library dir:  %s
;   Command line:
;     %s
;   Force field was read from the standard Gromacs share directory
;
''' % (fname, os.getlogin(), os.getuid(), os.uname()[1],
       now.strftime('%a. %B  %w %X %Y'), os.path.split(sys.argv[0])[1],
       __version__, os.path.split(sys.argv[0])[1], gmx.GROMACS_TOPDIR,
       ' '.join(sys.argv)))
            if include_itps is not None:
                dest.write('; Include forcefield parameters\n')
                if isinstance(include_itps, string_types):
                    dest.write('#include "%s"\n' % include_itps)
                else:
                    for include in include_itps:
                        dest.write('#include "%s"\n' % include)
            if combine is None:
                molecules = self.split()
                sysnum = 1
                names = []
                for molecule, num in molecules:
                    if len(molecule.residues) == 1:
                        title = molecule.residues[0].name
                    else:
                        title = 'system%d' % sysnum
                        sysnum += 1
                    names.append(title)
                    GromacsTopologyFile._write_molecule(molecule, dest, title)
                # System
                dest.write('[ system ]\n; Name\n')
                if self.title:
                    dest.write(self.title)
                else:
                    dest.write('Generic title')
                dest.write('\n\n')
                # Molecules
                dest.write('[ molecules ]\n; Compound       #mols\n')
                for i, (molecule, num) in enumerate(molecules):
                    dest.write('%-15s %6d\n' % (names[i], num))
            elif isinstance(combine, string_types) and combine.lower() == 'all':
                GromacsTopologyFile._write_molecule(self, dest, 'system')
            else:
                raise NotImplementedError('Specialized molecule splitting is '
                                          'not yet supported')
        finally:
            if own_handle:
                dest.close()

    #===================================================

    @staticmethod
    def _write_molecule(struct, dest, title):
        dest.write('\n[ moleculetype ]\n; Name            nrexcl\n')
        dest.write('%s          %d\n\n' % (title, struct.nrexcl))
        dest.write('[ atoms ]\n')
        dest.write(';   nr       type  resnr residue  atom   cgnr    '
                   'charge       mass  typeB    chargeB      massB\n')
        runchg = 0
        for residue in struct.residues:
            dest.write('; residue %4d %s rtp %s q %.1f\n' %
                       (residue.idx+1, residue.name, residue.name,
                        sum(a.charge for a in residue)))
            for atom in residue:
                runchg += atom.charge
                dest.write('%5d %10s %6d %6s %6s %6d %8.4f %10.3f   ; '
                           'qtot %.4f\n' % (atom.idx+1, atom.type,
                            residue.idx+1, residue.name, atom.name,
                            atom.idx+1, atom.charge, atom.mass, runchg))
        dest.write('\n')
        # Do valence terms now
        EPs = [a for a in struct.atoms if isinstance(a, ExtraPoint)]
        settle = False
        if len(struct.atoms) - len(EPs) == 3:
            try:
                oxy, = (a for a in struct.atoms if a.atomic_number == 8)
                hyd1, hyd2 = (a for a in struct.atoms if a.atomic_number == 1)
                settle = True
            except ValueError:
                pass
        if struct.bonds and not settle:
            dest.write('[ bonds ]\n')
            dest.write(';%6s %6s %5s %10s %10s %10s %10s\n' % ('ai', 'aj',
                       'funct', 'c0', 'c1', 'c2', 'c3'))
            for bond in struct.bonds:
                dest.write('%7d %6d %5d\n' % (bond.atom1.idx+1,
                           bond.atom2.idx+1, bond.funct))
            dest.write('\n')
        # Do the pair-exceptions
        if struct.adjusts:
            dest.write('[ pairs ]\n')
            dest.write(';%6s %6s %5s %10s %10s %10s %10s\n' % ('ai', 'aj',
                       'funct', 'c0', 'c1', 'c2', 'c3'))
            for adjust in struct.adjusts:
                dest.write('%7d %6d %5d\n' % (adjust.atom1.idx+1,
                           adjust.atom2.idx+1, adjust.funct))
            dest.write('\n')
        elif struct.dihedrals:
            dest.write('[ pairs ]\n')
            dest.write(';%6s %6s %5s %10s %10s %10s %10s\n' % ('ai', 'aj',
                       'funct', 'c0', 'c1', 'c2', 'c3'))
            # Get the 1-4 pairs from the dihedral list
            for dihed in struct.dihedrals:
                if dihed.ignore_end or dihed.improper: continue
                a1, a2 = dihed.atom1, dihed.atom4
                if a1 in a2.bond_partners or a1 in a2.angle_partners:
                    continue
                dest.write('%7d %6d %5d\n' % (a1.idx+1, a2.idx+1, 1))
            dest.write('\n')
        # Angles
        if struct.angles:
            dest.write('[ angles ]\n')
            dest.write(';%6s %6s %6s %5s %10s %10s %10s %10s\n' %
                       ('ai', 'aj', 'ak', 'funct', 'c0', 'c1', 'c2', 'c3'))
            for angle in struct.angles:
                dest.write('%7d %6d %6d %5d\n' % (angle.atom1.idx+1,
                           angle.atom2.idx+1, angle.atom3.idx+1,
                           angle.funct))
            dest.write('\n')
        # Dihedrals
        if struct.dihedrals:
            dest.write('[ dihedrals ]\n')
            dest.write((';%6s %6s %6s %6s %5s'+' %10s'*6) % ('ai', 'aj',
                       'ak', 'al', 'funct', 'c0', 'c1', 'c2', 'c3',
                       'c4', 'c5'))
            dest.write('\n')
            for dihed in struct.dihedrals:
                dest.write('%7d %6d %6d %6d %5d\n' % (dihed.atom1.idx+1,
                           dihed.atom2.idx+1, dihed.atom3.idx+1,
                           dihed.atom4.idx+1, dihed.funct))
            dest.write('\n')
        # RB-torsions
        if struct.rb_torsions:
            dest.write('[ dihedrals ]\n')
            dest.write((';%6s %6s %6s %6s %5s'+' %10s'*6) % ('ai', 'aj',
                       'ak', 'al', 'funct', 'c0', 'c1', 'c2', 'c3',
                       'c4', 'c5'))
            dest.write('\n')
            for dihed in struct.rb_torsions:
                dest.write('%7d %6d %6d %6d %5d\n' % (dihed.atom1.idx+1,
                           dihed.atom2.idx+1, dihed.atom3.idx+1,
                           dihed.atom4.idx+1, dihed.funct))
            dest.write('\n')
        # Impropers
        if struct.impropers:
            dest.write('[ dihedrals ]\n')
            dest.write((';%6s %6s %6s %6s %5s'+' %10s'*4) % ('ai', 'aj',
                       'ak', 'al', 'funct', 'c0', 'c1', 'c2', 'c3'))
            dest.write('\n')
            for dihed in struct.impropers:
                dest.write('%7d %6d %6d %6d %5d\n' % (dihed.atom1.idx+1,
                           dihed.atom2.idx+1, dihed.atom3.idx+1,
                           dihed.atom4.idx+1, dihed.funct))
            dest.write('\n')
        # Cmaps
        if struct.cmaps:
            dest.write('[ cmap ]\n')
            dest.write(';%6s %6s %6s %6s %6s %5s\n' % ('ai', 'aj', 'ak',
                       'al', 'am', 'funct'))
            for cmap in struct.cmaps:
                dest.write('%7d %6d %6d %6d %6d %5d\n' % (cmap.atom1.idx+1,
                           cmap.atom2.idx+1, cmap.atom3.idx+1,
                           cmap.atom4.idx+1, cmap.atom5.idx+1, cmap.funct))
        # See if this is a solvent molecule with 3 or fewer particles that can
        # be SETTLEd
        if settle:
            dest.write('[ settles ]\n')
            dest.write('; i	funct	doh	dhh\n')
            for b in oxy.bonds:
                if hyd1 in b:
                    if b.type is None:
                        # Use default values for TIPnP
                        doh = 0.09572
                    else:
                        doh = b.type.req / 10
                    break
            for b in hyd1.bonds:
                if hyd2 in b:
                    if b.type is None:
                        # Use default values for TIPnP
                        dhh = 0.15139
                    else:
                        dhh = b.type.req / 10
                    break
            dest.write('1     1   %.5f   %.5f\n\n' % (doh, dhh))
        # Virtual sites
        if EPs:
            ftypes = set(type(a.frame_type) for a in EPs)
            for ftype in ftypes:
                if ftype is TwoParticleExtraPointFrame:
                    dest.write('[ virtual_sites2 ]\n')
                    dest.write('; Site  from    funct  a\n')
                    for EP in EPs:
                        if not isinstance(EP.frame_type, ftype): continue
                        a1, a2 = EP.frame_type.get_atoms()
                        dest.write('%-5d %-4d %-4d %-4d   %.6f\n' %
                                   (EP.idx+1, a1.idx+1, a2.idx+1, 1,
                                    EP.get_weights()[0]))
                    dest.write('\n')
                elif ftype in (ThreeParticleExtraPointFrame,
                               OutOfPlaneExtraPointFrame):
                    dest.write('[ virtual_sites3 ]\n')
                    dest.write('; Site  from                   funct\n')
                    for EP in EPs:
                        if isinstance(EP.frame_type,
                                ThreeParticleExtraPointFrame):
                            a1, a2, a3 = EP.frame_type.get_atoms()
                            junk, w1, w2 = EP.frame_type.get_weights()
                            dest.write('%-5d %-4d %-4d %-4d %-4d   %.6f  %.6f\n'
                                    % (EP.idx+1, a1.idx+1, a2.idx+1, a3.idx+1,
                                       1, w1, w2))
                        elif isinstance(EP.frame_type,
                                OutOfPlaneExtraPointFrame):
                            a1, a2, a3 = EP.frame_type.get_atoms()
                            w1, w2, w3 = EP.frame_type.get_weights()
                            dest.write('%-5d %-4d %-4d %-4d %-4d   %.6f  %.6f  '
                                       '%.6f\n' % (EP.idx+1, a1.idx+1, a2.idx+1,
                                       a3.idx+1, w1, w2, w3))
                    dest.write('\n')
        # Do we need to list exclusions for systems with EPs?
        if EPs or settle:
            dest.write('[ exclusions ]\n')
            for i, atom in enumerate(struct.atoms):
                dest.write('%d' % (i+1))
                for a in atom.bond_partners:
                    dest.write(' %d' % (a.idx+1))
                for a in atom.angle_partners:
                    dest.write(' %d' % (a.idx+1))
                dest.write('\n')
            dest.write('\n')

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
