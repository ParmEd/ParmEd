"""
Tool for parsing and writing OFF library files to and from dictionaries of
ResidueTemplate objects
"""
import compat24

from chemistry import Atom, Bond
from chemistry.modeller.residue import ResidueTemplate
from collections import OrderedDict

import re

class AmberOFFLibrary(object):
    """
    Class containing static methods responsible for parsing and writing OFF
    libraries
    """

    # Useful regexes
    _headerre = re.compile(r'!!index *array *str')
    _resre = re.compile(r' *"(\S*)"$')
    _sec1re = re.compile(r'!entry\.(\S*)\.unit\.atoms *table *str *name *str'
                         r' *type *int *typex *int *resx *int *flags *int'
                         r' *seq *int *elmnt *dbl *chg')
    _sec2re = re.compile(r'!entry\.(\S*)\.unit\.atomspertinfo *table *str'
                         r' *pname *str *ptype *int *ptypex *int *pelmnt'
                         r' *dbl *pchg')
    _sec3re = re.compile(r'!entry.(\S*)\.unit\.boundbox *array *dbl')
    _sec4re = re.compile(r'!entry.(\S*)\.unit\.childsequence *single *int')
    _sec5re = re.compile(r'!entry.(\S*)\.unit\.connect *array *int')
    _sec6re = re.compile(r'!entry.(\S*)\.unit\.connectivity *table *int *atom1x'
                         r' *int *atom2x *int *flags')
    _sec7re = re.compile(r'!entry.(\S*)\.unit\.hierarchy *table *str *abovetype'
                         r' *int *abovex *str *belowtype *int *belowx')

    @staticmethod
    def parse(filename):
        """ Parses an Amber OFF library

        Parameters
        ----------
        filename : str or file-like iterable
            The file name or file object to parse. If it is an iterable, it will
            be exhausted

        Returns
        -------
        residues : OrderedDict {str : :class:`ResidueTemplate`}
            Dictionary pairing residue names with their :class:`ResidueTemplate`
            objects

        Raises
        ------
        ValueError if the first line does not match the file format. This line
        will be consumed

        IOError if filename is the name of a file that does not exist

        RuntimeError if EOF is reached prematurely or other formatting issues
        found
        """
        if isinstance(filename, basestring):
            fileobj = open(filename, 'r')
            own_handle = True
        else:
            fileobj = filename
            own_handle = False
        # Now parse the library file
        line = fileobj.readline()
        if not AmberOFFLibrary._headerre.match(line):
            raise ValueError('Unrecognized OFF file format')
        # Build the return value
        residues = OrderedDict()
        # Pull a list of all the residues we expect to find
        line = fileobj.readline()
        rematch = AmberOFFLibrary._resre.match(line)
        while rematch and line:
            name = rematch.groups()[0]
            residues[name] = None
            line = fileobj.readline()
            rematch = AmberOFFLibrary._resre.match(line)
        if not line:
            raise RuntimeError('Unexpected EOF in Amber OFF library')
        # Now make sure we have the next expected line
        rematch = AmberOFFLibrary._sec1re.match(line)
        if not rematch:
            raise RuntimeError('Expected atoms table not found')
        # Now just parse through the rest of the files and m
        
        return residues
