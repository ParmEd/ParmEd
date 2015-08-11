"""
This module contains classes for parsing and processing Amber parameter files.

Author: Jason M. Swails
Contributors:
Date: Aug. 11, 2015
"""
from __future__ import division, print_function

from contextlib import closing
from parmed.exceptions import ParameterError
from parmed.formats.registry import FileFormatType
from parmed.parameters import ParameterSet
from parmed.periodic_table import Mass, element_by_mass
from parmed.topologyobjects import (AtomType, BondType, AngleType, DihedralType,
                                    DihedralTypeList)
from parmed.utils.io import genopen
from parmed.utils.six import add_metaclass, string_types
import re

subs = dict(FLOATRE=r'([+-]?(?:\d+(?:\.\d*)?|\.\d+))')
_bondre = re.compile(r'(..?)-(..?) *%(FLOATRE)s *%(FLOATRE)s' % subs)
_anglere = re.compile(r'(..?)-(..?)-(..?) ' '*%(FLOATRE)s *%(FLOATRE)s' % subs)
_dihedre = re.compile(r'(..?)-(..?)-(..?)-(..?) *%(FLOATRE)s '
                      '*%(FLOATRE)s *%(FLOATRE)s *%(FLOATRE)s' % subs)
_sceere = re.compile(r'SCEE=%(FLOATRE)s' % subs)
_scnbre = re.compile(r'SCNB=%(FLOATRE)s' % subs)
_impropre = re.compile(r'(..?)-(..?)-(..?)-(..?) '
                       '*%(FLOATRE)s *%(FLOATRE)s *%(FLOATRE)s' % subs)

@add_metaclass(FileFormatType)
class AmberParameterSet(ParameterSet):
    """ Class storing parameters from an Amber parameter set

    Parameters
    ----------
    filenames : str or list of str
        Either the name of a file or a list of filenames from which parameters
        should be parsed.

    Notes
    -----
    Order is important in the list of files provided. The parameters are loaded
    in the order they are provided, and any parameters that are specified in
    multiple places are overwritten (that is, the *last* occurrence is the
    parameter type that is used)
    """

    @staticmethod
    def id_format(filename):
        """
        Identifies the file type as either an Amber-style frcmod or parm.dat
        file.

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is an Amber-style parameter file. False otherwise.
        """
        with closing(genopen(filename, 'r')) as f:
            f.readline()
            line = f.readline()
            if not line.strip(): # Must be an frcmod file
                while line and not line.strip():
                    line = f.readline()
                if not line:
                    return False
                if line.rstrip() not in ('MASS', 'BOND', 'ANGLE', 'ANGL',
                                         'DIHE', 'DIHED', 'DIHEDRAL', 'IMPR',
                                         'IMPROP', 'IMPROPER', 'NONB', 'NONBON',
                                         'NONBOND', 'NONBONDED'):
                    return False
            if line.rstrip() in ('MASS', 'BOND', 'ANGLE', 'ANGL', 'DIHE',
                                 'DIHED', 'DIHEDRAL', 'IMPR', 'IMPROP',
                                 'IMPROPER', 'NONB', 'NONBON', 'NONBOND',
                                 'NONBONDED'):
                return True # frcmod file
            # This must be an atom definition in the parm.dat file
            words = line.split()
            if len(words) < 2:
                return False
            # The first word is the atom type (must be <= 2 characters, and not
            # an integer)
            if len(words[0]) > 2: return False
            try:
                float(words[0])
                return False
            except ValueError:
                pass
            try:
                float(words[1])
            except ValueError:
                return False
            # Polarizability might not be present...
            try:
                float(words[2])
            except (IndexError, ValueError):
                # UGLY: Check the mass and make sure it matches the element's
                # mass to within 1 amu. Do our best to guess the element. If
                # it's a two-letter element, the element might be a 2-letter
                # element (like Br), or it might be a 1-letter element specified
                # by either the first or second letter. Check all possibilities.
                # Special-case instances like CA, which are just as likely (more
                # likely?) to be a carbon atom as it is to be a calcium (i.e.,
                # check for carbon atoms in anything that starts with C). If at
                # any point, what *should* be the mass doesn't match the mass of
                # the guessed element, tag this as *not* a parameter file
                if len(words[0]) == 2:
                    if words[0][0].isalpha():
                        if words[0][1].isalpha():
                            key = words[0][0].upper() + words[0][1].lower()
                            if (key in Mass and
                                    abs(Mass[key] - float(words[1])) > 1):
                                if (key[0] == 'C' and
                                        abs(Mass['C'] - float(words[1])) > 1):
                                    return False
                                elif key[0] != 'C':
                                    return False
                            elif key not in Mass:
                                if (key[0] in Mass and
                                        abs(Mass[key[0]] - float(words[1])) > 1):
                                    return False
                                else:
                                    return False
                        else:
                            key = words[0][0].upper()
                            if (key in Mass and
                                    abs(Mass[key] - float(words[1])) > 1):
                                return False
                            elif key not in Mass:
                                return False
                    else:
                        key = words[0][1].upper()
                        if key in Mass:
                            if abs(Mass[key] - float(words[1])) > 1:
                                return False
                        else:
                            return False
                else:
                    key = words[0][0].upper()
                    if key in Mass:
                        if abs(Mass[key] - float(words[1])) > 1:
                            return False
                    else:
                        return False
            if len(words) > 3:
                # Heuristic, but anything that comes after the polarizability is
                # a comment, and I have yet to see a leading comment that is a
                # number
                try:
                    float(words[3])
                    return False
                except ValueError:
                    return True
            else:
                return True

    def __init__(self, *filenames):
        super(AmberParameterSet, self).__init__()
        self.titles = []
        for filename in filenames:
            if isinstance(filename, string_types):
                self.load_parameters(filename)
            else:
                for fname in filename:
                    self.load_parameters(fname)

    def load_parameters(self, fname):
        """ Load a set of parameters from a single parameter file

        Parameters
        ----------
        fname : str or file-like
            Parameter file to parse
        """
        if isinstance(fname, string_types):
            f = genopen(fname, 'r')
            own_handle = True
        else:
            f = fname
            own_handle = False
        self.titles.append(f.readline().strip())
        for line in f:
            if not line.strip():
                return self._parse_frcmod(f, line)
            elif line.strip() in ('MASS', 'BOND', 'ANGLE', 'ANGL', 'DIHE',
                                  'DIHED', 'DIHEDRAL', 'IMPR', 'IMPROP',
                                  'IMPROPER', 'NONB', 'NONBON', 'NONBOND',
                                  'NONBONDED'):
                return self._parse_frcmod(f, line)
            else:
                return self._parse_parm_dat(f, line)

    def _parse_frcmod(self, f, line):
        """ Parses an frcmod file from an open file object """
        def fiter():
            yield line
            for l in f: yield l
        section = None
        dihed_type = None
        for line in fiter():
            line = line.rstrip()
            if not line: continue
            if line.startswith('MASS'):
                section = 'MASS'
                continue
            elif line.startswith('BOND'):
                section = 'BOND'
                continue
            elif line.startswith('ANGL'):
                section = 'ANGLE'
                continue
            elif line.startswith('DIHE'):
                section = 'DIHEDRAL'
                continue
            elif line.startswith('IMPR'):
                section = 'IMPROPER'
                continue
            elif line.startswith('NONB'):
                section = 'NONBOND'
                continue
            elif line.startswith('LJEDIT'):
                section = 'NBFIX'
                continue

            if section == 'MASS':
                words = line.split()
                try:
                    mass = float(words[1])
                except ValueError:
                    raise ParameterError('Could not convert mass to float [%s]'
                                         % words[1])
                except IndexError:
                    raise ParameterError('Error parsing MASS line. Not enough '
                                         'tokens')
                if words[0] in self.atom_types:
                    self.atom_types[words[0]].mass = mass
                else:
                    n_types = len(self.atom_types) + 1
                    atype = AtomType(words[0], len(self.atom_types)+1, mass,
                                     element_by_mass(mass))
                    self.atom_types[words[0]] = atype
            elif section == 'BOND':
                rematch = _bondre.match(line)
                if not rematch:
                    raise ParameterError('Could not understand BOND line [%s]' %
                                         line)
                a1, a2, k, eq = rematch.groups()
                a1 = a1.strip(); a2 = a2.strip()
                typ = BondType(float(k), float(eq))
                self.bond_types[(a1, a2)] = typ
                self.bond_types[(a2, a1)] = typ
            elif section == 'ANGLE':
                rematch = _anglere.match(line)
                if not rematch:
                    raise ParameterError('Could not understand ANGLE line [%s]'
                                         % line)
                a1, a2, a3, k, eq = rematch.groups()
                a1 = a1.strip(); a2 = a2.strip(); a3 = a3.strip()
                typ = AngleType(float(k), float(eq))
                self.angle_types[(a1, a2, a3)] = typ
                self.angle_types[(a3, a2, a1)] = typ
            elif section == 'DIHEDRAL':
                rematch = _dihedre.match(line)
                if not rematch:
                    raise ParameterError('Could not understand DIHEDRAL line '
                                         '[%s]' % line)
                a1, a2, a3, a4, div, k, phi, per = rematch.groups()
                scee = [float(x) for x in _sceere.findall(line)] or [1.2]
                scnb = [float(x) for x in _scnbre.findall(line)] or [2.0]
                a1 = a1.strip(); a2 = a2.strip();
                a3 = a3.strip(); a4 = a4.strip()
                per = float(per)
                typ = DihedralType(float(k)/float(div), abs(per), float(phi),
                                   scee[0], scnb[0])
                if per < 0:
                    # Part of a multi-term dihedral definition
                    if dihed_type is not None:
                        # Middle term of a multi-term dihedral
                        self.dihedral_types[dihed_type].append(typ)
                    else:
                        # First term of the multi-term dihedral
                        dihed_type = (a1, a2, a3, a4)
                        typs = DihedralTypeList()
                        typs.append(typ)
                        self.dihedral_types[dihed_type] = typs
                        self.dihedral_types[tuple(reversed(dihed_type))] = typs
                else:
                    if dihed_type is not None:
                        # Finish the existing multi-term dihedral
                        self.dihedral_types[dihed_type].append(typ)
                        dihed_type = None
                    else:
                        typs = DihedralTypeList()
                        typs.append(typ)
                        self.dihedral_types[(a1, a2, a3, a4)] = typs
                        self.dihedral_types[(a4, a3, a2, a1)] = typs
            elif section == 'IMPROPER':
                rematch = _impropre.match(line)
                if not rematch:
                    raise ParameterError('Could not understand IMPROPER line '
                                         '[%s]' % line)
                a1, a2, a3, a4, k, phi, per = rematch.groups()
                a1 = a1.strip(); a2 = a2.strip();
                a3 = a3.strip(); a4 = a4.strip()
                key = tuple(sorted([a1, a2, a3, a4]))
                self.improper_periodic_types[key] = \
                        DihedralType(float(k), float(per), float(phi))
            elif section == 'NONBOND':
                try:
                    atyp, rmin, eps = line.split()[:3]
                except ValueError:
                    raise ParameterError('Could not understand nonbond '
                                         'parameter line [%s]' % line)
                try:
                    self.atom_types[atyp].rmin = float(rmin)
                    self.atom_types[atyp].eps = float(eps)
                except KeyError:
                    raise ParameterError('Atom type %s not present in the '
                                         'database.' % atyp)
                except ValueError:
                    raise ParameterError('Could not convert nonbond parameters '
                                         'to floats [%s, %s]' % (rmin, eps))
            elif section == 'NBFIX':
                try:
                    a1, a2, rmin1, eps1, rmin2, eps2 = line.split()[:6]
                except ValueError:
                    raise ParameterError('Could not understand LJEDIT line [%s]'
                                         % line)
                try:
                    rmin1 = float(rmin1)
                    eps1 = float(eps1)
                    rmin2 = float(rmin2)
                    eps2 = float(eps2)
                except ValueError:
                    raise ParameterError('Could not convert LJEDIT parameters '
                                         'to floats.')
                self.nbfix_types[(min(a1, a2), max(a1, a2))] = \
                        (math.sqrt(eps1*eps2), rmin1+rmin2)

    def _parse_parm_dat(self, f, line):
        """ Internal parser for parm.dat files from open file handle """

    def write(self, dest, style='frcmod'):
        """ Writes a parm.dat file with the current parameters

        Parameters
        ----------
        dest : str or file-like
            The file name or file-like object to write the parameters to
        style : str, optional
            If 'frcmod', the parameters are written in frcmod-format. If 'parm',
            the parameters are written in parm.dat-format. Default is 'frcmod'
        """
        if isinstance(dest, string_types):
            outfile = genopen(dest, 'w')
        elif hasattr(dest, 'write'):
            outfile = dest
        else:
            raise TypeError('Cannot write parameter set to a %s' %
                            type(dest).__name__)

        if style not in ('frcmod', 'parm'):
            raise ValueError('style must be either frcmod or parm, not %s' %
                             style)

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

        if isinstance(dest, string_types):
            outfile.close()
