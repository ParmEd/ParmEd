"""
Classes helpful for reading/storing Amber parameters
"""
from __future__ import division
from chemistry.topologyobjects import BondType, AngleType, DihedralType
from chemistry.amber.readparm import AmberParm
from chemistry.constants import RAD_TO_DEG, SMALL
from chemistry.exceptions import AmberParameterWarning
import warnings

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondParam(object):
   
    def __init__(self, atype1, atype2, rk=None, req=None, bondtype=None):
        """ Initializes a Bond parameter based on two atom types """
        self.atype1, self.atype2 = atype1, atype2
        if bondtype is not None:
            if not isinstance(bondtype, BondType):
                raise TypeError('bondtype expected to be BondType instance')
            self.type = bondtype
        else:
            self.type = BondType(float(rk), float(req))

    def __str__(self):
        return '%s-%s   %8.3f  %6.3f' % (self.atype1.ljust(2),
                self.atype2.ljust(2), self.type.k, self.type.req)

    def __eq__(self, other):
        """ Two bonds are equal if their atom types and bond type is equal """
        if self.type != other.type: return False
        return (
            (self.atype1 == other.atype1 and self.atype2 == other.atype2) or
            (self.atype2 == other.atype1 and self.atype1 == other.atype2)
        )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleParam(object):
   
    def __init__(self, atype1, atype2, atype3, thetk=None, theteq=None,
                 angletype=None):
        self.atype1, self.atype2, self.atype3 = atype1, atype2, atype3
        if angletype is not None:
            if not isinstance(angletype, AngleType):
                raise TypeError('angletype expected to be AngleType instance')
            self.type = angletype
        else:
            self.type = AngleType(float(thetk), float(theteq))

    def __str__(self):
        return '%s-%s-%s   %8.3f  %6.3f' % (self.atype1.ljust(2),
                self.atype2.ljust(2), self.atype3.ljust(2), self.type.k,
                self.type.theteq * RAD_TO_DEG)

    def __eq__(self, other):
        """ Two angles are equal if their atom types and bond type is equal """
        if self.type != other.type:
            return False
        if self.atype2 != other.atype2:
            return False
        return ((self.atype1 == other.atype1 and self.atype3 == other.atype3) or
                (self.atype1 == other.atype3 and self.atype3 == other.atype1))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _DihedralTerm(object):
    """ A single term in a (potentially multiterm) dihedral """
    def __init__(self, idivf=1, pk=None, phase=None, periodicity=None,
                 dihedraltype=None, scee=1.2, scnb=2.0, dihtype='normal'):
        self.idivf = int(idivf)
        if not dihtype in ('normal', 'improper'):
            raise ValueError('dihtype must be normal or improper')
        improper = dihtype == 'improper'
        self.dihtype = dihtype
        if dihedraltype is not None:
            if not isinstance(dihedraltype, DihedralType):
                raise TypeError('dihedraltype expected to be '
                                'DihedralType instance')
            self.type = dihedraltype
        else:
            self.type = DihedralType(pk / self.idivf, float(periodicity),
                                    float(phase), float(scee), float(scnb),
                                    improper=improper)

    def parmline(self, multiterm=False):
        if self.dihtype == 'improper' or not multiterm:
            return str(self)
        else:
            return '%4i %8.3f %8.3f %5.1f    SCEE=%s SCNB=%s' % (self.idivf,
                self.type.phi_k * self.idivf, self.type.phase * RAD_TO_DEG,
                -self.type.per, self.type.scee, self.type.scnb)
         
    def __str__(self):
        if self.dihtype == 'improper':
            return '%8.3f %8.3f %5.1f' % (self.type.phi_k, 
                    self.type.phase * RAD_TO_DEG, self.type.per)
        else:
            return '%4i %8.3f %8.3f %5.1f    SCEE=%s SCNB=%s' % (self.idivf,
                    self.type.phi_k * self.idivf, self.type.phase * RAD_TO_DEG,
                    self.type.per, self.type.scee, self.type.scnb)

    def __eq__(self, other):
        if self.dihtype != other.dihtype and 'improper' in (self.dihtype,
                                                            other.dihtype):
            return False
        return (abs(self.type.phi_k - other.type.phi_k) < SMALL and
                abs(self.type.phase - other.type.phase) < SMALL and
                abs(self.type.per - other.type.per) < SMALL)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralParam(list):
   
    def __init__(self, atype1, atype2, atype3, atype4):
        list.__init__(self)
        self.atype1, self.atype2 = atype1, atype2
        self.atype3, self.atype4 = atype3, atype4
   
    def add_term_from_args(self, *args, **kwargs):
        """ Convert the args to a _DihedralTerm first """
        term = _DihedralTerm(*args, **kwargs)
        self.add_term(term)

    def add_term(self, term):
        """ Add this to the current dihedral IFF it is not already there """
        if not isinstance(term, _DihedralTerm):
            raise TypeError('DihedralParam must add a _DihedralTerm instance')
        for oterm in self:
            if term == oterm:
                return
        if len(self) == 0:
            return list.append(self, term)
        elif term.dihtype == 'improper' or self[0].dihtype == 'improper':
            raise ValueError('Improper dihedrals cannot have multiple terms')
        # Add dihedrals in terms of decreasing periodicity
        for i, oterm in enumerate(self):
            if oterm.type.per > term.type.per: continue
            return list.insert(self, i, term)
        return list.append(self, term)

    def __str__(self):
        if len(self) == 1:
            return '%s-%s-%s-%s %s' % (self.atype1.ljust(2),
                    self.atype2.ljust(2), self.atype3.ljust(2),
                    self.atype4.ljust(2), self[0].parmline())
        retstr = ''
        for i in xrange(len(self)-1):
            retstr += '%s-%s-%s-%s %s\n' % (self.atype1.ljust(2),
                    self.atype2.ljust(2), self.atype3.ljust(2),
                    self.atype4.ljust(2), self[i].parmline(multiterm=True))
        i = len(self) - 1
        retstr += '%s-%s-%s-%s %s' % (self.atype1.ljust(2),
                    self.atype2.ljust(2), self.atype3.ljust(2),
                    self.atype4.ljust(2), self[i].parmline(multiterm=False))
        return retstr

    def __eq__(self, other):
        sameatoms, sameorder = self.same_atoms((other.atype1, other.atype2,
                                                other.atype3, other.atype4))
        if not sameatoms:
            return False
        if len(self) != len(other):
            return False
        for i in xrange(len(self)):
            if self[i] != other[i]:
                return False
        return True

    def same_atoms(self, atomlist):
        """
        Determine if two dihedrals are assigned to the same sets of atom types

        Parameters
        ----------
        atomlist : list
            4-element list of atom types to compare against this DihedralParam

        Returns
        -------
        bool, bool
            First bool is True if all atoms are the same; False otherwise
            Second bool is True if atoms are in the same order; False otherwise

        Notes
        -----
        If this torsion is an improper, the first atom is fixed and the other 3
        atoms must be the same (but in any order). If this torsion is a proper,
        then the torsions must match in either the forward or reverse
        directions, only.
        """
        if self[0].dihtype == 'improper':
            # For impropers, the third atom is the central atom and the other 3
            # can be in any order
            if self.atype3 != atomlist[2]: return False, False
            if (self.atype1 == atomlist[0] and self.atype2 == atomlist[1] and
                self.atype4 == atomlist[3]):
                return True, True
            # Make sure every atom type is unique so every atom type is added to
            # the set (but we know it is NOT the same order now)
            set1, set2 = set() , set()
            for x, y in zip([atomlist[0], atomlist[1], atomlist[3]],
                            [self.atype1, self.atype2, self.atype4]):
                if x in set1:
                    i = 0
                    while '%s%d%d' % (x, i, i) in set1: i += 1
                    set1.add('%s%d%d' % (x, i, i))
                else:
                    set1.add(x)
                if y in set2:
                    i = 0
                    while '%s%d%d' % (y, i, i) in set2: i += 1
                    set2.add('%s%d%d' % (y, i, i))
                else:
                    set2.add(y)
            return set1 == set2, False
        # If we reached here, this is a proper dihedral
        if self.atype2 == atomlist[1] and self.atype3 == atomlist[2]:
            same = self.atype1 == atomlist[0] and self.atype4 == atomlist[3]
            if not same:
                # Check case when atype2 == atype3 and we're not the same --
                # check for reversed order
                eq = self.atype1 == atomlist[3] and self.atype4 == atomlist[0]
                return eq, False
            return same, same
        if self.atype3 == atomlist[1] and self.atype2 == atomlist[2]:
            return self.atype1==atomlist[3] and self.atype4==atomlist[0], False
        return False, False

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AtomType(object):
   
    def __init__(self, atype, mass, rmin, epsilon):
        self.type = atype
        self.mass = float(mass)
        self.rmin = float(rmin)
        self.epsilon = float(epsilon)

    def __str__(self):
        return '%s%6.3f' % (self.type.ljust(6), self.mass)

    def __eq__(self, other):
        return (self.type == other.type and self.mass == other.mass and
                self.rmin == other.rmin and self.epsilon == other.epsilon)

    def lennard_jones(self):
        return '%s  %8.4f %8.4f' % (self.type.ljust(2), self.rmin, self.epsilon)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ParamList(list):
   
    type = None

    def __init__(self, unique_only=True):
        self.unique_only = unique_only

    def add_param(self, param):
        if not isinstance(param, self.type):
            raise TypeError('Added parameter should be type %s' %
                            self.type.__name__)
        if not self.unique_only or not param in self:
            list.append(self, param)

    def __contains__(self, thing):
        for param in self:
            if param == thing:
                return True
        return False

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondParamList(ParamList):
   
    type = BondParam

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleParamList(ParamList):

    type = AngleParam

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralParamList(ParamList):
   
    type = _DihedralTerm

    def __init__(self, unique_only=True):
        if not unique_only:
            warnings.warn('DihedralParamList can only store unique dihedrals',
                          AmberParameterWarning)
        self.unique_only = True

    def add_param(self, atype1, atype2, atype3, atype4, param):
        """
        Override this since dihedral parameters are actually lists of dihedral
        terms, so adding a DihedralParam to a DihedralParamList is different
        than the other parameter lists
        """
        if not isinstance(param, self.type):
            raise TypeError('Added parameter should be type %s' %
                            self.type.__name__)
        # Loop through all DihedralParam objects in this list. Add it to any
        # DihedralParam that has the same atoms. If it's a duplicate term, it
        # will be discarded.
        added = False
        for oparam in self:
            # if one is an improper and the other is not, skip over this. If
            # both are impropers, let it through to be screened for equality
            # If both are propers, also let them through.
            if ((oparam[0].dihtype == 'improper' and
                 param.dihtype == 'normal') or
                (oparam[0].dihtype == 'normal' and
                 param.dihtype == 'improper')):
                continue
            if oparam.same_atoms((atype1, atype2, atype3, atype4))[0]:
                added = True
                oparam.add_term(param)
                break
        if not added:
            # We need to make a new dihedral parameter now
            p = DihedralParam(atype1, atype2, atype3, atype4)
            p.add_term(param)
            list.append(self, p)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AtomTypeList(ParamList):
   
    type = AtomType

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ParameterSet(object):
   
    def __init__(self):
        self.atoms = AtomTypeList(unique_only=True)
        self.bonds = BondParamList(unique_only=True)
        self.angles = AngleParamList(unique_only=True)
        self.dihedrals = DihedralParamList(unique_only=True)

    def _add_atom(self, *args, **kwargs):
        self.atoms.add_param(AtomType(*args, **kwargs))

    def _add_bond(self, *args, **kwargs):
        self.bonds.add_param(BondParam(*args, **kwargs))

    def _add_angle(self, *args, **kwargs):
        self.angles.add_param(AngleParam(*args, **kwargs))

    def _add_dihedral(self, atype1, atype2, atype3, atype4, *args, **kwargs):
        self.dihedrals.add_param(atype1, atype2, atype3, atype4,
                                 _DihedralTerm(*args, **kwargs))

    def load_from_parm(self, parm):
        if not isinstance(parm, AmberParm):
            raise TypeError('ParameterSet.load_from_parm() expects AmberParm')
        # Loop through all atoms, adding it to the list of atoms
        for atom in parm.atoms:
            rmin = parm.LJ_radius[atom.nb_idx-1]
            eps = parm.LJ_depth[atom.nb_idx-1]
            self._add_atom(atom.type, atom.mass, rmin, eps)
        # Loop through all bonds
        for bond in parm.bonds_without_h:
            self._add_bond(bond.atom1.type, bond.atom2.type, 
                           bondtype=bond.type)
        for bond in parm.bonds_inc_h:
            self._add_bond(bond.atom1.type, bond.atom2.type,
                           bondtype=bond.type)
        # Loop through all angles
        for angle in parm.angles_without_h:
            self._add_angle(angle.atom1.type, angle.atom2.type,
                            angle.atom3.type, angletype=angle.type)
        for angle in parm.angles_inc_h:
            self._add_angle(angle.atom1.type, angle.atom2.type,
                            angle.atom3.type, angletype=angle.type)
        # Loop through all dihedrals
        for i, dihedral in enumerate(parm.dihedrals_without_h):
            if dihedral.improper:
                term = 'improper'
            else:
                term = 'normal'
            self._add_dihedral(dihedral.atom1.type, dihedral.atom2.type,
                               dihedral.atom3.type, dihedral.atom4.type,
                               idivf=1, dihedraltype=dihedral.type,
                               dihtype=term)
        for dihedral in parm.dihedrals_inc_h:
            if dihedral.improper:
                term = 'improper'
            else:
                term = 'normal'
            self._add_dihedral(dihedral.atom1.type, dihedral.atom2.type,
                               dihedral.atom3.type, dihedral.atom4.type,
                               idivf=1, dihedraltype=dihedral.type,
                               dihtype=term)

    def write(self, dest):
        """ Writes a parm.dat file with the current parameters """
        if isinstance(dest, str):
            outfile = open(dest, 'w')
        elif hasattr(dest, 'write'):
            outfile = dest
        else:
            raise TypeError('Cannot write parameter set to type %s' %
                            type(dest).__name__)

        # Write the atom mass
        outfile.write('MASS\n')
        for atom in self.atoms:
            outfile.write('%s%6.3f\n' % (atom.type.ljust(6), atom.mass))
        outfile.write('\n')
        # Write the bonds
        outfile.write('BOND\n')
        for bond in self.bonds:
            outfile.write('%s\n' % bond)
        outfile.write('\n')
        # Write the angles
        outfile.write('ANGLE\n')
        for angle in self.angles:
            outfile.write('%s\n' % angle)
        outfile.write('\n')
        # Write the dihedrals
        outfile.write('DIHE\n')
        for dihedral in self.dihedrals:
            if dihedral[0].dihtype == 'improper': continue
            outfile.write('%s\n' % dihedral)
        outfile.write('\n')
        # Write the impropers
        outfile.write('IMPROPER\n')
        for dihedral in self.dihedrals:
            if dihedral[0].dihtype != 'improper': continue
            outfile.write('%s\n' % dihedral)
        outfile.write('\n')
        # Write the LJ terms
        outfile.write('NONB\n')
        for atom in self.atoms:
            outfile.write('%s\n' % atom.lennard_jones())
        outfile.write('\n')

        if isinstance(dest, str):
            outfile.close()
