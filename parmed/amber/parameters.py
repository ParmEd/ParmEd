"""
This module contains classes for parsing and processing Amber parameter files.

Author: Jason M. Swails
Contributors:
Date: Aug. 11, 2015
"""
from __future__ import division, print_function

from collections import defaultdict
from contextlib import closing
import math
import os
import re
import warnings
from .offlib import AmberOFFLibrary
from ..constants import TINY
from ..exceptions import ParameterError, AmberWarning, ParameterWarning
from ..modeller.residue import ResidueTemplateContainer
from ..formats.mol2 import Mol2File
from ..formats.registry import FileFormatType
from ..parameters import ParameterSet
from ..periodic_table import Mass, element_by_mass, AtomicNum
from ..topologyobjects import (AtomType, BondType, AngleType, DihedralType,
                                    DihedralTypeList)
from ..utils.io import genopen
from ..utils.six import add_metaclass, string_types, iteritems, PY2
from ..utils.six.moves import map
if PY2:
    from collections import Sequence
else:
    from collections.abc import Sequence

# parameter file regexes
subs = dict(FLOATRE=r'([+-]?(?:\d+(?:\.\d*)?|\.\d+))',
            FILENAMERE=r'''(".+?"|'.+?'|[\S]*)''')
_bondre = re.compile(r'(..?)-(..?)\s+%(FLOATRE)s\s+%(FLOATRE)s' % subs)
_anglere = re.compile(r'(..?)-(..?)-(..?)\s+%(FLOATRE)s\s+%(FLOATRE)s' % subs)
_dihedre = re.compile(r'(..?)-(..?)-(..?)-(..?)\s+%(FLOATRE)s\s+'
                      '%(FLOATRE)s\s+%(FLOATRE)s\s+%(FLOATRE)s' % subs)
_dihed2re = re.compile(r'\s*%(FLOATRE)s\s+%(FLOATRE)s\s+%(FLOATRE)s\s+'
                       '%(FLOATRE)s' % subs)
_sceere = re.compile(r'SCEE=\s*%(FLOATRE)s' % subs)
_scnbre = re.compile(r'SCNB=\s*%(FLOATRE)s' % subs)
_impropre = re.compile(r'(..?)-(..?)-(..?)-(..?)\s+'
                       '%(FLOATRE)s\s+%(FLOATRE)s\s+%(FLOATRE)s' % subs)
# Leaprc regexes
_atomtypere = re.compile(r"""({\s*["']([\w\+\-]+)["']\s*["'](\w+)["']\s*"""
                         r"""["'](\w+)["']\s*})""")
_loadparamsre = re.compile(r'loadamberparams\s+%(FILENAMERE)s' % subs, re.I)
_loadoffre = re.compile(r'loadoff\s+%(FILENAMERE)s' % subs, re.I)
_loadmol2re = re.compile(r'(\S+)\s*=\s*loadmol[23]\s*%(FILENAMERE)s' % subs, re.I)
del subs

def _find_amber_file(fname, search_oldff):
    """
    Finds an Amber file. Looks in the current directory, then the following
    locations:

    - $AMBERHOME/dat/leap/lib
    - $AMBERHOME/dat/leap/parm
    """
    from parmed.amber import AMBERHOME
    if len(fname) > 1 and {fname[0], fname[-1]} in ({'"'}, {"'"}):
        # Strip quotes
        fname = fname[1:-1]
    leapdir = os.path.join(AMBERHOME, 'dat', 'leap')
    paths = [os.getcwd(), os.path.join(leapdir, 'lib'), os.path.join(leapdir, 'parm')]
    if search_oldff:
        paths.append(os.path.join(leapdir, 'lib', 'oldff'))
    for path in paths:
        if os.path.exists(os.path.join(path, fname)):
            return os.path.join(path, fname)
    raise FileNotFoundError('Cannot find Amber file [%s] in paths %s' % (fname, paths))

@add_metaclass(FileFormatType)
class AmberParameterSet(ParameterSet):
    """ Class storing parameters from an Amber parameter set

    Parameters
    ----------
    filenames : str, list of str, file-like, or list of file-like; optional
        Either the name of a file or a list of filenames from which parameters
        should be parsed.

    Notes
    -----
    Order is important in the list of files provided. The parameters are loaded
    in the order they are provided, and any parameters that are specified in
    multiple places are overwritten (that is, the *last* occurrence is the
    parameter type that is used)

    See Also
    --------
    :class:`parmed.parameters.ParameterSet`
    """

    #===================================================

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
                if line.rstrip() not in ('MASS', 'BOND', 'ANGLE', 'ANGL', 'DIHE', 'DIHED',
                                         'DIHEDRAL', 'IMPR', 'IMPROP', 'IMPROPER', 'NONB', 'NONBON',
                                         'NONBOND', 'NONBONDED'):
                    return False
            if line.rstrip() in ('MASS', 'BOND', 'ANGLE', 'ANGL', 'DIHE', 'DIHED', 'DIHEDRAL',
                                 'IMPR', 'IMPROP', 'IMPROPER', 'NONB', 'NONBON', 'NONBOND',
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
                            if key in Mass and abs(Mass[key] - float(words[1])) > 1:
                                if key[0] == 'C' and abs(Mass['C'] - float(words[1])) > 1:
                                    return False
                                elif key[0] != 'C':
                                    return False
                            elif key not in Mass:
                                if key[0] in Mass and abs(Mass[key[0]] - float(words[1])) > 1:
                                    return False
                                else:
                                    return False
                        else:
                            key = words[0][0].upper()
                            if key in Mass and abs(Mass[key] - float(words[1])) > 1:
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

    #===================================================

    def __init__(self, *filenames):
        super(AmberParameterSet, self).__init__()
        self.default_scee = 1.2
        self.default_scnb = 2.0
        self.titles = []
        for filename in filenames:
            if isinstance(filename, string_types):
                if AmberOFFLibrary.id_format(filename):
                    self.residues.update(AmberOFFLibrary.parse(filename))
                else:
                    self.load_parameters(filename)
            elif isinstance(filename, Sequence):
                for fname in filename:
                    if AmberOFFLibrary.id_format(fname):
                        self.residues.update(AmberOFFLibrary.parse(fname))
                    else:
                        self.load_parameters(fname)
            else:
                # Assume open file object
                self.load_parameters(filename)

    #===================================================

    @classmethod
    def from_leaprc(cls, fname, search_oldff=False):
        """ Load a parameter set from a leaprc file

        Parameters
        ----------
        fname : str or file-like
            Name of the file or open file-object from which a leaprc-style file
            will be read

        search_oldff : bool, optional, default=False
            If True, search the oldff directories in the main Amber leap
            folders. Default is False

        Notes
        -----
        This does not read all parts of a leaprc file -- only those pertinent to
        defining force field information. For instance, the following sections
        and commands are processed:

        - addAtomTypes
        - loadAmberParams
        - loadOFF
        - loadMol2
        - loadMol3
        """
        params = cls()
        if isinstance(fname, string_types):
            f = genopen(fname, 'r')
            own_handle = True
        else:
            f = fname
            own_handle = False
        # To make parsing easier, and because leaprc files are usually quite
        # short, I'll read the whole file into memory
        def joinlines(lines):
            newlines = []
            composite = []
            for line in lines:
                if line.endswith('\\\n'):
                    composite.append(line[:-2])
                    continue
                else:
                    composite.append(line)
                newlines.append(''.join(composite))
                composite = []
            if composite:
                newlines.append(''.join(composite))
            return newlines
        lines = joinlines(map(lambda line: line if '#' not in line else line[:line.index('#')], f))
        text = ''.join(lines)
        if own_handle: f.close()
        lowertext = text.lower() # commands are case-insensitive
        # Now process the parameter files
        def process_fname(fname):
            if fname[0] in ('"', "'"):
                fname = fname[1:-1]
            fname = fname.replace('_BSTOKEN_', r'\ ').replace(r'\ ', ' ')
            return fname
        for line in lines:
            line = line.replace(r'\ ', '_BSTOKEN_')
            if _loadparamsre.findall(line):
                fname = process_fname(_loadparamsre.findall(line)[0])
                params.load_parameters(_find_amber_file(fname, search_oldff))
            elif _loadoffre.findall(line):
                fname = process_fname(_loadoffre.findall(line)[0])
                params.residues.update(AmberOFFLibrary.parse(_find_amber_file(fname, search_oldff)))
            elif _loadmol2re.findall(line):
                (resname, fname), = _loadmol2re.findall(line)
                residue = Mol2File.parse(_find_amber_file(fname, search_oldff))
                if isinstance(residue, ResidueTemplateContainer):
                    warnings.warn('Multi-residue mol2 files not supported by tleap. Loading anyway '
                                  'using names in mol2', AmberWarning)
                    for res in residue:
                        params.residues[res.name] = res
                else:
                    params.residues[resname] = residue
        # Now process the addAtomTypes
        try:
            idx = lowertext.index('addatomtypes')
        except ValueError:
            # Does not exist in this file
            atom_types_str = ''
        else:
            i = idx + len('addatomtypes')
            while i < len(text) and text[i] != '{':
                if text[i] not in '\r\n\t ':
                    raise ParameterError('Unsupported addAtomTypes syntax in leaprc file')
                i += 1
            if i == len(text):
                raise ParameterError('Unsupported addAtomTypes syntax in leaprc file')
            # We are at our first brace
            chars = []
            nopen = 1
            i += 1
            while i < len(text):
                char = text[i]
                if char == '{':
                    nopen += 1
                elif char == '}':
                    nopen -= 1
                    if nopen == 0: break
                elif char == '\n':
                    char = ' '
                chars.append(char)
                i += 1
            atom_types_str = ''.join(chars).strip()
        for _, name, symb, hyb in _atomtypere.findall(atom_types_str):
            if symb not in AtomicNum:
                raise ParameterError('%s is not a recognized element' % symb)
            if name in params.atom_types:
                params.atom_types[name].atomic_number = AtomicNum[symb]
        return params

    #===================================================

    @classmethod
    def from_structure(cls, struct):
        """ Extracts known parameters from a Structure instance

        Parameters
        ----------
        struct : :class:`parmed.structure.Structure`
            The parametrized ``Structure`` instance from which to extract
            parameters into a ParameterSet

        Returns
        -------
        params : :class:`ParameterSet`
            The parameter set with all parameters defined in the Structure
        """
        return super(AmberParameterSet, cls).from_structure(struct, allow_unequal_duplicates=False)

    #===================================================

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
        try:
            for line in f:
                if not line.strip():
                    return self._parse_frcmod(f, line)
                elif line.strip() in ('MASS', 'BOND', 'ANGLE', 'ANGL', 'DIHE', 'DIHED', 'DIHEDRAL',
                                      'IMPR', 'IMPROP', 'IMPROPER', 'NONB', 'NONBON', 'NONBOND',
                                      'NONBONDED'):
                    return self._parse_frcmod(f, line)
                else:
                    return self._parse_parm_dat(f, line)
        finally:
            if own_handle:
                f.close()

    #===================================================

    def _parse_frcmod(self, f, line):
        """ Parses an frcmod file from an open file object """
        def fiter():
            yield line
            for l in f: yield l
        section = None
        finished_diheds = defaultdict(lambda: True)
        key = None
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
                self._process_mass_line(line)
            elif section == 'BOND':
                self._process_bond_line(line)
            elif section == 'ANGLE':
                self._process_angle_line(line)
            elif section == 'DIHEDRAL':
                key = self._process_dihedral_line(line, finished_diheds, key)
            elif section == 'IMPROPER':
                self._process_improper_line(line)
            elif section == 'NONBOND':
                self._process_nonbond_line(line)
            elif section == 'NBFIX':
                self._process_nbfix_line(line)

    #===================================================

    def _parse_parm_dat(self, f, line):
        """ Internal parser for parm.dat files from open file handle """
        def fiter():
            yield line
            for l in f: yield l
            # Keep yielding empty string after file has ended
            yield ''
        fiter = fiter()
        rawline = next(fiter)
        finished_diheds = defaultdict(lambda: True)
        # Parse the masses
        while rawline:
            line = rawline.strip()
            if not line:
                break
            self._process_mass_line(line)
            rawline = next(fiter)
        next(fiter) # Skip the list of hydrophobic atom types
        # Process the bonds
        rawline = next(fiter)
        while rawline:
            line = rawline.strip()
            if not line:
                break
            self._process_bond_line(line)
            rawline = next(fiter)
        # Process the angles
        rawline = next(fiter)
        while rawline:
            line = rawline.strip()
            if not line:
                break
            self._process_angle_line(line)
            rawline = next(fiter)
        # Process the dihedrals
        rawline = next(fiter)
        key = None
        while rawline:
            line = rawline.strip()
            if not line:
                break
            key = self._process_dihedral_line(line, finished_diheds, key)
            rawline = next(fiter)
        # Process the impropers
        rawline = next(fiter)
        while rawline:
            line = rawline.strip()
            if not line:
                break
            self._process_improper_line(line)
            rawline = next(fiter)
        # Process the 10-12 terms
        rawline = next(fiter)
        while rawline:
            line = rawline.strip()
            if not line:
                break
            try:
                a1, a2, acoef, bcoef = line.split()[:4]
                acoef = float(acoef)
                bcoef = float(bcoef)
            except ValueError:
                raise ParameterError('Trouble parsing 10-12 terms')
            if acoef != 0 or bcoef != 0:
                raise ParameterError('10-12 potential not supported in AmberParameterSet currently')
            rawline = next(fiter)
        # Process 12-6 terms. Get Equivalencing first
        rawline = next(fiter)
        equivalent_ljtypes = dict()
        equivalent_types = defaultdict(list)
        while rawline:
            line = rawline.strip()
            if not line:
                break
            words = line.split()
            for typ in words[1:]:
                equivalent_ljtypes[typ] = words[0]
                equivalent_types[words[0]].append(typ)
            rawline = next(fiter)
        words = next(fiter).split()
        if len(words) < 2:
            raise ParameterError('Could not parse the kind of nonbonded '
                                 'parameters in Amber parameter file')
        if words[1].upper() != 'RE':
            raise ParameterError('Only RE nonbonded parameters supported')
        rawline = next(fiter)
        while rawline:
            line = rawline.strip()
            if not line:
                break
            self._process_nonbond_line(line)
            rawline = next(fiter)
        # Now assign all of the equivalenced atoms
        for atyp, otyp in iteritems(equivalent_ljtypes):
            otyp = self.atom_types[otyp]
            if atyp in self.atom_types:
                if (self.atom_types[atyp].rmin is not None and
                        self.atom_types[atyp].epsilon is not None):
                    if (abs(otyp.epsilon-self.atom_types[atyp].epsilon) > TINY
                            or abs(otyp.rmin-self.atom_types[atyp].rmin) > TINY):
                        warnings.warn('Equivalency defined between %s and %s but parameters are '
                                      'not equal' % (otyp.name, atyp), AmberWarning)
                        # Remove from equivalent types
                        equivalent_types[otyp.name].remove(atyp)
                        continue
                self.atom_types[atyp].set_lj_params(otyp.epsilon, otyp.rmin)
        line = next(fiter).strip()
        if line == 'LJEDIT':
            rawline = next(fiter)
            while rawline:
                line = rawline.strip()
                if not line:
                    break
                self._process_nbfix_line(line, equivalent_types)
                rawline = next(fiter)

    #===================================================

    # Private methods for processing parts of the file
    def _process_mass_line(self, line):
        words = line.split()
        try:
            mass = float(words[1])
        except ValueError:
            raise ParameterError('Could not convert mass to float [%s]' % words[1])
        except IndexError:
            raise ParameterError('Error parsing MASS line. Not enough tokens')
        if words[0] in self.atom_types:
            self.atom_types[words[0]].mass = mass
        elif words[0] in ('EP', 'LP'):
            atype = AtomType(words[0], len(self.atom_types)+1, mass, 0)
            self.atom_types[words[0]] = atype
        else:
            atype = AtomType(words[0], len(self.atom_types)+1, mass,
                             AtomicNum[element_by_mass(mass)])
            self.atom_types[words[0]] = atype

    #===================================================

    def _process_bond_line(self, line):
        rematch = _bondre.match(line)
        if not rematch:
            raise ParameterError('Could not understand BOND line [%s]' % line)
        a1, a2, k, eq = rematch.groups()
        a1 = a1.strip(); a2 = a2.strip()
        typ = BondType(float(k), float(eq))
        self.bond_types[(a1, a2)] = typ
        self.bond_types[(a2, a1)] = typ

    def _process_angle_line(self, line):
        rematch = _anglere.match(line)
        if not rematch:
            raise ParameterError('Could not understand ANGLE line [%s]' % line)
        a1, a2, a3, k, eq = rematch.groups()
        a1 = a1.strip(); a2 = a2.strip(); a3 = a3.strip()
        typ = AngleType(float(k), float(eq))
        self.angle_types[(a1, a2, a3)] = typ
        self.angle_types[(a3, a2, a1)] = typ

    def _process_dihedral_line(self, line, finished_diheds, last_key):
        """ Processes a dihedral line, possibly part of a multi-term dihedral

        Parameters
        ----------
        line : str
            Line of the file that contains a dihedral term
        finished_diheds : dict
            Dictionary of dihedral parameters whose final term has been read in
            already (which means additional terms will overwrite, not add)
        last_key : str or None
            If not None, this is the key for the last dihedral type that should
            be implied if the atom types are missing. Atom types seem to only be
            required for the first term in a multi-term torsion definition

        Returns
        -------
        key or None
            If a negative periodicity indicates another term is coming, the
            current key is returned so it can be passed as key to the next
            _process_dihedral_call
        """
        rematch = _dihedre.match(line)
        if not rematch and last_key is None:
            raise ParameterError('Could not understand DIHEDRAL line '
                                 '[%s]' % line)
        elif not rematch:
            rematch = _dihed2re.match(line)
            if not rematch:
                raise ParameterError('Could not understand DIHEDRAL line '
                                     '[%s]' % line)
            div, k, phi, per = rematch.groups()
            key = last_key
            rkey = tuple(reversed(key))
            assert key in finished_diheds
            if finished_diheds[key]:
                raise AssertionError('Cannot have an implied torsion that '
                                     'has already finished!')
        else:
            a1, a2, a3, a4, div, k, phi, per = rematch.groups()
            a1, a2, a3, a4 =  a1.strip(), a2.strip(), a3.strip(), a4.strip()
            key = (a1, a2, a3, a4)
            rkey = (a4, a3, a2, a1)
            if last_key is not None and (last_key != key and last_key != rkey):
                warnings.warn('Expecting next term in dihedral %r, got '
                              'definition for dihedral %r' % (last_key, key),
                              ParameterWarning)
        scee = [float(x) for x in _sceere.findall(line)] or [1.2]
        scnb = [float(x) for x in _scnbre.findall(line)] or [2.0]
        per = float(per)
        typ = DihedralType(float(k)/float(div), abs(per), float(phi),
                           scee[0], scnb[0])
        if finished_diheds[key]:
            # This dihedral is already finished its definition, which means
            # we go ahead and add a new one to override it
            typs = DihedralTypeList()
            typs.append(typ)
            self.dihedral_types[key] = self.dihedral_types[rkey] = typs
        else:
            self.dihedral_types[key].append(typ)
        finished_diheds[key] = finished_diheds[rkey] = per >= 0
        if per < 0:
            return key

    def _process_improper_line(self, line):
        rematch = _impropre.match(line)
        if not rematch:
            raise ParameterError('Could not understand IMPROPER line '
                                    '[%s]' % line)
        a1, a2, a3, a4, k, phi, per = rematch.groups()
        a1 = a1.strip(); a2 = a2.strip();
        a3 = a3.strip(); a4 = a4.strip()
        # Pre-sort the improper types, assuming atom3 is the central atom (which
        # it must be in Amber parameter files!!!!)
        a1, a2, a4 = sorted([a1, a2, a4])
        key = (a1, a2, a3, a4)
        self.improper_periodic_types[key] = \
                DihedralType(float(k), float(per), float(phi))

    def _process_nonbond_line(self, line):
        try:
            atyp, rmin, eps = line.split()[:3]
        except ValueError:
            raise ParameterError('Could not understand nonbond parameter line '
                                 '[%s]' % line)
        try:
            self.atom_types[atyp].set_lj_params(float(eps), float(rmin))
        except KeyError:
            raise ParameterError('Atom type %s not present in the database.' %
                                 atyp)
        except ValueError:
            raise ParameterError('Could not convert nonbond parameters to '
                                 'floats [%s, %s]' % (rmin, eps))

    def _process_nbfix_line(self, line, equivalents=None):
        try:
            a1, a2, rmin1, eps1, rmin2, eps2 = line.split()[:6]
        except ValueError:
            raise ParameterError('Could not understand LJEDIT line [%s]' % line)
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
        if equivalents is not None:
            # We need to add the same nbfixes to all atom types that are
            # equivalent to the atom types defined in the LJEDIT line
            for oa1 in equivalents[a1]:
                self.nbfix_types[(min(oa1, a2), max(oa1, a2))] = \
                        (math.sqrt(eps1*eps2), rmin1+rmin2)
            for oa2 in equivalents[a2]:
                self.nbfix_types[(min(a1, oa2), max(a1, oa2))] = \
                        (math.sqrt(eps1*eps2), rmin1+rmin2)
            # Now do the equivalenced of atom 1 with the equivalenced of atom 2
            for oa1 in equivalents[a1]:
                for oa2 in equivalents[a2]:
                    self.nbfix_types[(min(oa1, oa2), max(oa1, oa2))] = \
                            (math.sqrt(eps1*eps2), rmin1+rmin2)

    #===================================================

    def write(self, dest, title='Created by ParmEd', style='frcmod'):
        """ Writes a parm.dat file with the current parameters

        Parameters
        ----------
        dest : str or file-like
            The file name or file-like object to write the parameters to
        title : str, optional
            The title of the frcmod to write. Default is 'Created by ParmEd'
        style : str, optional
            If 'frcmod', the parameters are written in frcmod-format. If 'parm',
            the parameters are written in parm.dat-format. Default is 'frcmod'
        """
        if isinstance(dest, string_types):
            outfile = genopen(dest, 'w')
            own_handle = True
        else:
            outfile = dest
            own_handle = False

        if style not in ('frcmod', 'parm'):
            raise ValueError('style must be either frcmod or parm, not %s' % style)

        outfile.write(title.rstrip('\r\n'))
        outfile.write('\n')
        # Write the atom mass
        outfile.write('MASS\n')
        for atom, typ in iteritems(self.atom_types):
            outfile.write('%s%6.3f\n' % (atom.ljust(6), typ.mass))
        outfile.write('\n')
        # Write the bonds
        outfile.write('BOND\n')
        done = set()
        for (a1, a2), typ in iteritems(self.bond_types):
            if id(typ) in done: continue
            done.add(id(typ))
            outfile.write('%s-%s   %8.3f  %6.3f\n' % (a1.ljust(2), a2.ljust(2), typ.k, typ.req))
        outfile.write('\n')
        # Write the angles
        outfile.write('ANGLE\n')
        done = set()
        for (a1, a2, a3), typ in iteritems(self.angle_types):
            if id(typ) in done: continue
            done.add(id(typ))
            outfile.write('%s-%s-%s   %8.3f  %6.3f\n' % (a1.ljust(2), a2.ljust(2), a3.ljust(2),
                                                         typ.k, typ.theteq))
        outfile.write('\n')
        # Write the dihedrals
        outfile.write('DIHE\n')
        done = set()
        for (a1, a2, a3, a4), typ in iteritems(self.dihedral_types):
            if id(typ) in done: continue
            done.add(id(typ))
            if isinstance(typ, DihedralType) or len(typ) == 1:
                if not isinstance(typ, DihedralType):
                    typ = typ[0]
                outfile.write('%s-%s-%s-%s %4i %14.8f %8.3f %5.1f    SCEE=%s SCNB=%s\n' %
                              (a1.ljust(2), a2.ljust(2), a3.ljust(2), a4.ljust(2), 1, typ.phi_k,
                               typ.phase, typ.per, typ.scee, typ.scnb))
            else:
                for dtyp in typ[:-1]:
                    outfile.write('%s-%s-%s-%s %4i %14.8f %8.3f %5.1f    SCEE=%s SCNB=%s\n' %
                                  (a1.ljust(2), a2.ljust(2), a3.ljust(2), a4.ljust(2), 1,
                                   dtyp.phi_k, dtyp.phase, -dtyp.per, dtyp.scee, dtyp.scnb))
                dtyp = typ[-1]
                outfile.write('%s-%s-%s-%s %4i %14.8f %8.3f %5.1f    SCEE=%s SCNB=%s\n' %
                              (a1.ljust(2), a2.ljust(2), a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k,
                               dtyp.phase, dtyp.per, dtyp.scee, dtyp.scnb))
        outfile.write('\n')
        # Write the impropers
        outfile.write('IMPROPER\n')
        written_impropers = dict()
        for (a1, a2, a3, a4), typ in iteritems(self.improper_periodic_types):
            # Make sure wild-cards come at the beginning
            if a2 == 'X':
                assert a4 == 'X', 'Malformed generic improper!'
                a1, a2, a3, a4 = a2, a4, a3, a1
            elif a4 == 'X':
                a1, a2, a3, a4 = a4, a1, a3, a2
            a1, a2, a4 = sorted([a1, a2, a4])
            if (a1, a2, a3, a4) in written_impropers:
                if written_impropers[(a1, a2, a3, a4)] != typ:
                    raise ValueError('Multiple impropers with the same atom set not allowed')
                continue
            outfile.write('%s-%s-%s-%s %14.8f %8.3f %5.1f\n' %
                          (a1.ljust(2), a2.ljust(2), a3.ljust(2), a4.ljust(2), typ.phi_k,
                           typ.phase, typ.per))
            written_impropers[(a1, a2, a3, a4)] = typ
        outfile.write('\n')
        # Write the LJ terms
        outfile.write('NONB\n')
        for atom, typ in iteritems(self.atom_types):
            outfile.write('%s  %12.8f %12.8f\n' % (atom.ljust(2), typ.rmin, typ.epsilon))
        outfile.write('\n')
        # Write the NBFIX terms
        if self.nbfix_types:
            outfile.write('LJEDIT\n')
            for (a1, a2), (eps, rmin) in iteritems(self.nbfix_types):
                outfile.write('%s %s %12.8f %12.8f %12.8f %12.8f\n' %
                              (a1.ljust(2), a2.ljust(2), eps, rmin/2, eps, rmin/2))

        if own_handle:
            outfile.close()
