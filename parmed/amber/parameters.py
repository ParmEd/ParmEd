"""
This module contains classes for parsing and processing Amber parameter files.

Author: Jason M. Swails
Contributors:
Date: Aug. 11, 2015
"""
from __future__ import division, print_function

from contextlib import closing
from parmed.formats.registry import FileFormatType
from parmed.parameters import ParameterSet
from parmed.periodic_table import Mass
from parmed.utils.io import genopen
from parmed.utils.six import add_metaclass

@add_metaclass(FileFormatType)
class AmberParameterSet(ParameterSet):

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
