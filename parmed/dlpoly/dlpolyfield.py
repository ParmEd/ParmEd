"""
This module contains functionality relevant to building a DLPOLY topology file
"""
from __future__ import print_function, division, absolute_import

from collections import defaultdict
import copy
import math

from ..exceptions import DlpolyError
from ..formats.registry import FileFormatType
from ..parameters import ParameterSet
from ..structure import Structure
from ..topologyobjects import (ExtraPoint, DihedralType, DihedralTypeList,
                                    UnassignedAtomType)
from ..utils.io import genopen
from ..utils.six import add_metaclass, string_types, iteritems
from ..utils.six.moves import range


# Dlpoly uses "funct" flags in its parameter files to indicate what kind of
# functional form is used for each of its different parameter types. This is
# taken from the topdirs.c source code file along with a table in the Dlpoly
# user manual. The table below summarizes my findings, for reference:

# Bonds
# -----
#  1 - F_BONDS : simple harmonic potential
#  2 - F_G96BONDS : fourth-power potential
#  3 - F_MORSE : morse potential
#  4 - F_CUBICBONDS : cubic potential
#  5 - F_CONNBONDS : not even implemented in DLPOLY
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

class _Defaults(object):
    """ Global properties of force fields as implemented in DLPOLY """
    def __init__(self, nbfunc=1, comb_rule=2, gen_pairs='no',
                 fudgeLJ=1.0, fudgeQQ=1.0):
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

    def __repr__(self):
        return ('<_Defaults: nbfunc=%d, comb-rule=%d, gen-pairs="%s", '
                'fudgeLJ=%g, fudgeQQ=%g>' % (self.nbfunc, self.comb_rule,
                    self.gen_pairs, self.fudgeLJ, self.fudgeQQ))

    def __eq__(self, other):
        return (self.nbfunc == other.nbfunc and
                self.comb_rule == other.comb_rule and
                self.gen_pairs == other.gen_pairs and
                self.fudgeLJ == other.fudgeLJ and
                self.fudgeQQ == other.fudgeQQ)

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
            self.fudgeQQ = value
        else:
            raise IndexError('Index %d out of range' % idx)

@add_metaclass(FileFormatType)
class DlpolyFieldFile(Structure):
    """ Class providing a writer for a DLPOLY topology file
    """

    #===================================================

    def __init__(self):
        super(DlpolyFieldFile, self).__init__()
        self.parameterset = None
        self.defaults = _Defaults(gen_pairs='yes') # make ParmEd's default yes

    #===================================================

    @classmethod
    def from_structure(cls, struct, copy=False):
        """ Instantiates a DlpolyFieldFile instance from a Structure

        Parameters
        ----------
        struct : :class:`parmed.Structure`
            The input structure to generate from
        copy : bool, optional
            If True, assign from a *copy* of ``struct`` (this is a lot slower).
            Default is False

        Returns
        -------
        dlpolyfield : :class:`DlpolyFieldFile`
            The topology file defined by the given struct
        """
        from copy import copy as _copy
        dlpolyfield = cls()
        if copy:
            struct = _copy(struct)
            struct.join_dihedrals()
        dlpolyfield.atoms = struct.atoms
        dlpolyfield.residues = struct.residues
        dlpolyfield.bonds = struct.bonds
        dlpolyfield.angles = struct.angles
        dlpolyfield.dihedrals = struct.dihedrals
        dlpolyfield.impropers = struct.impropers
        dlpolyfield.cmaps = struct.cmaps
        dlpolyfield.rb_torsions = struct.rb_torsions
        dlpolyfield.urey_bradleys = struct.urey_bradleys
        dlpolyfield.adjusts = struct.adjusts
        dlpolyfield.bond_types = struct.bond_types
        dlpolyfield.angle_types = struct.angle_types
        dlpolyfield.dihedral_types = struct.dihedral_types
        dlpolyfield.improper_types = struct.improper_types
        dlpolyfield.cmap_types = struct.cmap_types
        dlpolyfield.rb_torsion_types = struct.rb_torsion_types
        dlpolyfield.urey_bradley_types = struct.urey_bradley_types
        dlpolyfield.adjust_types = struct.adjust_types
        dlpolyfield.combining_rule = struct.combining_rule
        dlpolyfield.box = struct.box
        dlpolyfield.nrexcl = struct.nrexcl
        if (struct.trigonal_angles or
                struct.out_of_plane_bends or
                struct.pi_torsions or
                struct.stretch_bends or
                struct.torsion_torsions or
                struct.chiral_frames or
                struct.multipole_frames):
            raise TypeError('DlpolyFieldFile does not support Amoeba FF')
        # Now check what the 1-4 scaling factors should be
        if hasattr(struct, 'defaults') and isinstance(struct.defaults,
                                                      _Defaults):
            dlpolyfield.defaults = struct.defaults
        else:
            scee_values = set()
            scnb_values = set()
            if struct.adjusts:
                for adjust in struct.adjusts:
                    if adjust.type is None: continue
                    scee_values.add(1/adjust.type.chgscale)
                    # Do not add scnb_values, since we can just set explicit
                    # exception pair parameters in DLPOLY (which this structure
                    # already has)
                # In order to specify specific pair parameters, we need to set
                # gen_pairs to 'no' so that the pair-specific L-J parameters are
                # printed to the topology file (rather than being auto-created)
                dlpolyfield.defaults.gen_pairs = 'no'
            else:
                for dihedral in struct.dihedrals:
                    if dihedral.type is None or dihedral.ignore_end: continue
                    if isinstance(dihedral.type, DihedralTypeList):
                        for dt in dihedral.type:
                            if dt.scee:
                                scee_values.add(dt.scee)
                            if dt.scnb:
                                scnb_values.add(dt.scnb)
                    else:
                        if dihedral.type.scee:
                            scee_values.add(dihedral.type.scee)
                        if dihedral.type.scnb:
                            scnb_values.add(dihedral.type.scnb)
            if len(set('%.5f' % x for x in scee_values)) > 1:
                raise DlpolyError('Structure has mixed 1-4 scaling which is '
                                   'not supported by Dlpoly')
            scee_values = list(scee_values)
            scnb_values = list(scnb_values)
            if len(set('%.5f' % x for x in scee_values)) == 1:
                dlpolyfield.defaults.fudgeQQ = 1/scee_values[0]
            else:
                dlpolyfield.defaults.fudgeQQ = 1.0
            if len(set('%.5f' % x for x in scnb_values)) == 1:
                dlpolyfield.defaults.fudgeLJ = 1/scnb_values[0]
            else:
                dlpolyfield.defaults.fudgeLJ = 1.0
        if dlpolyfield.combining_rule == 'geometric':
            dlpolyfield.defaults.comb_rule = 3

        dlpolyfield.parameterset = ParameterSet.from_structure(struct,
                                            allow_unequal_duplicates=True)
        return dlpolyfield

    #===================================================

    def write(self, dest, combine=None, parameters='inline'):
        """ Write a Dlpoly Topology File from a Structure

        Parameters
        ----------
        dest : str or file-like
            The name of a file or a file object to write the Dlpoly topology to
        combine : 'all', None, or list of iterables, optional
            If None, no molecules are combined into a single moleculetype. If
            'all', all molecules are combined into a single moleculetype.
            Otherwise, the list of molecule indices (start from 0) will control
            which atoms are combined into single moleculetype's. Default is None
        parameters : 'inline' or str or file-like object, optional
            This specifies where parameters should be printed. If 'inline'
            (default), the parameters are written on the same lines as the
            valence terms are defined on. Any other string is interpreted as a
            filename for an ITP that will be written to and then included at the
            top of `dest`. If it is a file-like object, parameters will be
            written there.  If parameters is the same as ``dest``, then the
            parameter types will be written to the same topologyfile.

        Raises
        ------
        ValueError if the same molecule number appears in multiple combine lists
        TypeError if the dest input cannot be parsed
        ValueError if the combine, or parameters input cannot be parsed
        """
        own_handle = False
        fname = ''
        params = ParameterSet.from_structure(self, allow_unequal_duplicates=True)
        if isinstance(dest, string_types):
            fname = '%s ' % dest
            dest = genopen(dest, 'w')
            own_handle = True
        elif not hasattr(dest, 'write'):
            raise TypeError('dest must be a file name or file-like object')

        # Determine where to write the parameters
        own_parfile_handle = False
        if parameters == 'inline':
            parfile = dest
        elif isinstance(parameters, string_types):
            if parameters == fname.strip():
                parfile = dest
            else:
                own_parfile_handle = True
                parfile = genopen(parameters, 'w')
        elif hasattr(parameters, 'write'):
            parfile = parameters
        else:
            raise ValueError('parameters must be "inline", a file name, or '
                             'a file-like object')

        # Determine where to write the molecules
        own_molfile_handle = False
        _molfile = dest

        # Error-checking for combine
        if combine is not None:
            if isinstance(combine, string_types):
                if combine.lower() != 'all':
                    raise ValueError('combine must be None, list of indices, '
                                     'or "all"')
            else:
                combine_lists = []
                for indices in combine:
                    indices = sorted(set(indices))
                    if any((indices[i+1] - indices[i]) != 1
                                for i in range(len(indices)-1)):
                        raise ValueError('Can only combine adjacent molecules')
                    combine_lists.append(indices)
        try:
            # Write the header
            if self.title:
                dest.write(self.title+"\n")
            else:
                dest.write('AMBER_SYSTEM\n') # pragma: no cover
            dest.write('UNITS kcal\n')
            if combine is None:
                molecules = self.split()
                dest.write('MOLECULAR types %d\n'%(len(molecules)))
                sysnum = 1
                names = []
                nameset = set()
                for molecule, num in molecules:
                    if len(molecule.residues) == 1:
                        title = molecule.residues[0].name
                        if title in nameset:
                            orig = title
                            sfx = 2
                            while title in nameset:
                                title = '%s%d' % (orig, sfx)
                                sfx += 1
                    else:
                        title = 'system%d' % sysnum
                        sysnum += 1
                    names.append(title)
                    nameset.add(title)
                    DlpolyFieldFile._write_molecule(molecule, _molfile,
                                                        title, len(num), params,
                                                        parameters == 'inline')
            elif isinstance(combine, string_types) and combine.lower() == 'all':
                dest.write('MOLECULAR types 1\n')
                DlpolyFieldFile._write_molecule(self, _molfile, 'system', 1,
                                                    params,
                                                    parameters == 'inline')
                # Molecules
                dest.write('[ molecules 2 ]\n; Compound       #mols\n')
                dest.write('%-15s %6d\n' % ('system', 1))
            else:
                molecules = self.split()
                nmols = sum(len(m[1]) for m in molecules)
                moleculedict = dict()
                # Hash our molecules by indices
                for m, num in molecules:
                    for i in num:
                        moleculedict[i] = m
                combined_molecules = []
                for cl in combine_lists:
                    counts = defaultdict(int)
                    mols_in_mol = []
                    for molid in cl:
                        try:
                            mol = moleculedict[molid]
                        except KeyError:
                            raise IndexError('Molecule ID out of range')
                        counts[id(moleculedict[molid])] += 1
                        if counts[id(moleculedict[molid])] == 1:
                            mols_in_mol.append(mol)
                    if counts[id(mols_in_mol[0])] > 1:
                        combmol = mols_in_mol[0] * counts[id(mols_in_mol[0])]
                    else:
                        combmol = copy.copy(mols_in_mol[0])
                    for i, mol in enumerate(mols_in_mol):
                        if i == 0: continue
                        assert id(mol) in counts and counts[id(mol)] > 0
                        if counts[id(mol)] > 1:
                            combmol += mol * counts[id(mol)]
                        else:
                            combmol += mol
                    combined_molecules.append((combmol, cl[0], len(cl)))
                    nmols -= (len(cl) - 1)
                # combined_molecules now contains a list of tuples, and that
                # tuple stores the combined molecule, first molecule index of
                # the pre-combined molecule, and how many molecules were
                # combined

                # Sort combined molecules by starting location
                combined_molecules.sort(key=lambda x: x[1])
                new_molecules = []
                counts = defaultdict(set)
                cmc = 0 # Combined Molecule Counter
                add = 0 # How many molecules to "skip" due to combining
                for i in range(nmols):
                    ii = i + add
                    if (cmc < len(combined_molecules) and
                            combined_molecules[cmc][1] == ii):
                        new_molecules.append([combined_molecules[cmc][0],
                                              set([i])])
                        add += combined_molecules[cmc][2] - 1
                        cmc += 1
                    elif len(counts[id(moleculedict[ii])]) == 0:
                        counts[id(moleculedict[ii])].add(i)
                        new_molecules.append([moleculedict[ii],
                                              counts[id(moleculedict[ii])]])
                    else:
                        counts[id(moleculedict[ii])].add(i)
                sysnum = 1
                names = []
                nameset = set()
                dest.write('MOLECULAR types %d\n'%(len(new_molecules)))
                for molecule, num in new_molecules:
                    if len(molecule.residues) == 1:
                        title = molecule.residues[0].name
                        if title in nameset:
                            orig = title
                            sfx = 2
                            while title in nameset:
                                title = '%s%d' % (orig, sfx)
                                sfx += 1
                    else:
                        title = 'system%d' % sysnum
                        sysnum += 1
                    names.append(title)
                    nameset.add(title)
                    DlpolyFieldFile._write_molecule(molecule, _molfile,
                                                        title, len(num), params,
                                                        parameters == 'inline')
        finally:
            # Print all VdW parameters for the full system
            nattyps = len(params.atom_types)
            parfile.write('VDW %d\n'%(int(nattyps*(nattyps+1)/2)))
            cnt1 = 0
            for key1, atom_type1 in iteritems(params.atom_types):
                cnt1 += 1
                cnt2 = 0
                for key2, atom_type2 in iteritems(params.atom_types):
                    cnt2 += 1
                    if (cnt2 < cnt1): continue
                    eps = math.sqrt(atom_type1.epsilon*atom_type2.epsilon)
                    sig = (atom_type1.sigma+atom_type2.sigma)/2
                    parfile.write('%-8s%-8s%4s %14.8f %14.8f\n'%
                          (atom_type1, atom_type2, 'lj', eps, sig))
            # Close statement
            dest.write('CLOSE\n')
            # Close handles
            if own_handle:
                dest.close()
            if own_parfile_handle:
                parfile.close()
            if own_molfile_handle:
                _molfile.close()

    #===================================================

    @staticmethod
    def _write_molecule(struct, dest, title, nmols, params, writeparams):
        # Printing molecule header
        dest.write('%s\n' % (title))
        dest.write('NUMMOLS %d\n'%(nmols))
        # Getting number of atoms
        natoms = 0
        for residue in struct.residues:
            for atom in residue:
                natoms += 1
        # Print masses and charges
        dest.write('ATOMS %d\n'%(natoms))
        for residue in struct.residues:
            for atom in residue:
                dest.write('%-8s %10.6f %12.8f 1\n' %
                           (atom.type, atom.mass, atom.charge))
        # Do valence terms now
        EPs = [a for a in struct.atoms if isinstance(a, ExtraPoint)]
        if len(struct.atoms) - len(EPs) == 3:
            try:
                oxy, = (a for a in struct.atoms if a.atomic_number == 8)
                hyd1, hyd2 = (a for a in struct.atoms if a.atomic_number == 1)
            except ValueError:
                pass
        # Print bonds
        if struct.bonds:
            conv = 2.0
            dest.write('BONDS %d\n'%(len(struct.bonds)))
            for bond in struct.bonds:
                if (isinstance(bond.atom1, ExtraPoint) or isinstance(bond.atom2, ExtraPoint)):
                    continue # pragma: no cover
                # bond.funct ==1 is for a simple harmonic potential
                if (bond.funct != 1):
                    raise DlpolyError('Bond between atoms %d and %d is of an invalid type!'
                                       % (bond.atom1.idx+1, bond.atom2.idx+1))
                dest.write('%4s %6d %6d' % ('harm', bond.atom1.idx+1, bond.atom2.idx+1))
                if bond.type is None:
                    dest.write('\n')
                    continue # pragma: no cover
                key = (_gettype(bond.atom1), _gettype(bond.atom2))
                if writeparams or key not in params.bond_types or \
                        bond.type != params.bond_types[key]:
                    dest.write('   %12.6f %12.6f %12.6f %12.6f' % (bond.type.k*conv, bond.type.req, 0.0, 0.0))
                dest.write('\n')
        # Angles
        if struct.angles:
            conv = 2.0
            dest.write('ANGLES %d\n'%(len(struct.angles)))
            for angle in struct.angles:
                # angle.funct ==1 is for a simple harmonic potential
                if (angle.funct != 1):
                    raise DlpolyError('Angle between atoms %d, %d and %d is of an invalid type!'
                                       % (angle.atom1.idx+1, angle.atom2.idx+1,
                                          angle.atom3.idx+1))
                dest.write('%4s %6d %6d %6d' % ('harm', angle.atom1.idx+1, angle.atom2.idx+1,
                           angle.atom3.idx+1))
                if angle.type is None:
                    dest.write('\n')
                    continue
                key = (_gettype(angle.atom1), _gettype(angle.atom2), _gettype(angle.atom3))
                param_equal = params.angle_types.get(key) == angle.type
                if writeparams or not param_equal:
                    dest.write('   %12.7f %12.7f %12.7f %12.7f' % (angle.type.k*conv,
                               angle.type.theteq, 0.0, 0.0))
                dest.write('\n')
        # Dihedrals
        if struct.dihedrals:
            dest.write('DIHEDRALS %d\n'%(len(struct.dihedrals)))
            conv = 1.0
            for dihed in struct.dihedrals:
                # dihed.funct ==1 or 4 is for a simple harmonic potential
                if (dihed.funct != 1 and dihed.funct != 4):
                    raise DlpolyError('Dihedral between atoms %d, %d, %d and %d is of an invalid type!'
                                       % (dihed.atom1.idx+1, dihed.atom2.idx+1,
                                          dihed.atom3.idx+1, dihed.atom4.idx+1))
                dest.write('%-4s %6d %6d %6d %6d' % ('cos', dihed.atom1.idx+1, dihed.atom2.idx+1,
                           dihed.atom3.idx+1, dihed.atom4.idx+1))
                if dihed.type is None:
                    dest.write('\n')
                    continue
                if dihed.improper:
                    typedict = params.improper_periodic_types
                else:
                    typedict = params.dihedral_types
                key = (_gettype(dihed.atom1), _gettype(dihed.atom2),
                        _gettype(dihed.atom3), _gettype(dihed.atom4))
                if writeparams or key not in typedict or \
                        _diff_diheds(dihed.type, typedict[key]):
                    scee = 0.0
                    if (dihed.type.scee > 0.0): scee = 1.0/dihed.type.scee
                    scnb = 0.0
                    if (dihed.type.scnb > 0.0): scnb = 1.0/dihed.type.scnb
                    dest.write('  %12.7f  %12.7f  %4d  %10.5f  %10.5f' % (dihed.type.phi_k*conv,
                               dihed.type.phase, int(dihed.type.per), scee,
                               scnb))
                dest.write('\n')
        # Impropers
        if struct.impropers:
            dest.write('INVERSIONS %d\n'%(len(struct.impropers)))
            conv = 2.0
            for dihed in struct.impropers:
                # dihed.funct ==1 is for a simple harmonic potential
                if (dihed.funct != 1):
                    raise DlpolyError('Dihedral between atoms %d, %d, %d and %d is of an invalid type!'
                                       % (dihed.atom1.idx+1, dihed.atom2.idx+1,
                                          dihed.atom3.idx+1, dihed.atom4.idx+1))
                dest.write('%-4s %6d %6d %6d %6d' % ('harm', dihed.atom1.idx+1,
                           dihed.atom2.idx+1, dihed.atom3.idx+1,
                           dihed.atom4.idx+1))
                if dihed.type is None:
                    dest.write('\n')
                    continue
                # BUGBUG: We always write improper types since we don't
                # currently store the correct ordering of the types in the
                # improper section
                dest.write('  %12.7f  %12.7f\n' % (dihed.type.psi_k*conv, dihed.type.psi_eq))
        # Finish
        dest.write('FINISH\n')

def _diff_diheds(dt1, dt2):
    """ Determine if 2 dihedrals are *really* different. dt1 can either be a
    DihedralType or a DihedralTypeList or dt1 can be a DihedralType and dt2 can
    be a DihedralTypeList.  This returns True if dt1 == dt2 *or* dt1 is equal to
    the only element of dt2
    """
    if type(dt1) is type(dt2) and dt1 == dt2:
        return False
    if isinstance(dt2, DihedralTypeList) and isinstance(dt1, DihedralType):
        if len(dt2) == 1 and dt2[0] == dt1: return False
    return True

def _gettype(atom):
    if atom.atom_type not in (None, UnassignedAtomType):
        return atom.atom_type.bond_type
    return atom.type
