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
    _sec3re = re.compile(r'!entry\.(\S*)\.unit\.boundbox *array *dbl')
    _sec4re = re.compile(r'!entry\.(\S*)\.unit\.childsequence *single *int')
    _sec5re = re.compile(r'!entry\.(\S*)\.unit\.connect *array *int')
    _sec6re = re.compile(r'!entry\.(\S*)\.unit\.connectivity *table *int'
                         r' *atom1x *int *atom2x *int *flags')
    _sec7re = re.compile(r'!entry\.(\S*)\.unit\.hierarchy *table *str'
                         r' *abovetype *int *abovex *str *belowtype *int'
                         r' *belowx')
    _sec8re = re.compile(r'!entry\.(\S*)\.unit\.name *single *str')
    _sec9re = re.compile(r'!entry\.(\S*)\.unit\.positions *table *dbl *x'
                         r' *dbl *y *dbl *z')
    _sec10re = re.compile(r'!entry\.(\S*)\.unit\.residueconnect *table'
                          r' *int *c1x *int *c2x *int *c3x *int *c4x *int'
                          r' *c5x *int *c6x')
    _sec11re = re.compile(r'!entry\.(\S*)\.unit\.residues *table *str *name'
                          r' *int *seq *int *childseq *int *startatomx *str'
                          r' *restype *int *imagingx')
    _sec12re = re.compile(r'!entry\.(\S*)\.unit\.residuesPdbSequenceNumber'
                          r' *array *int')
    _sec13re = re.compile(r'!entry\.(\S*)\.unit\.solventcap *array *dbl')
    _sec14re = re.compile(r'!entry\.(\S*)\.unit\.velocities *table *dbl *x'
                          r' *dbl *y *dbl *z')

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
        while line:
            rematch = AmberOFFLibrary._sec1re.match(line)
            if not rematch:
                raise RuntimeError('Expected atoms table not found')
            name = rematch.groups()[0]
            residues[name] = AmberOFFLibrary._parse_residue(fileobj, name)
            line = fileobj.readline()
        
        return residues

    @staticmethod
    def _parse_residue(fileobj, name):
        """
        Parses the residue information out of the OFF file assuming the file
        is pointed at the first line of an atoms table section of the OFF file

        Parameters
        ----------
        fileobj : file-like
            Assumed to be open for read, this file is parsed until the *next*
            atom table is read
        name : str
            The name of the residue being processed right now
        """
        nres = 1
        templ = ResidueTemplate(name)
        line = fileobj.readline()
        while line[0] != '!':
            nam, typ, typx, resx, flags, seq, elmnt, chg = line.split()
            nam = _strip_enveloping_quotes(nam)
            typ = _strip_enveloping_quotes(typ)
            typx = int(typx)
            resx = int(resx)
            flags = int(flags)
            seq = int(seq)
            elmnt = int(elmnt)
            chg = float(chg)
            atom = Atom(atomic_number=elmnt, type=typ, name=nam, charge=chg)
            templ.add_atom(atom)
            line = fileobj.readline()
        # Make sure we get the next section
        rematch = AmberOFFLibrary._sec2re.match(line)
        if not rematch:
            raise RuntimeError('Expected pertinfo table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        line = fileobj.readline()
        while line[0] != '!':
            if not line:
                raise RuntimeError('Unexpected EOF in Amber OFF library')
            # Not used, just skip
            # TODO sanity check
            line = fileobj.readline()
        rematch = AmberOFFLibrary._sec3re.match(line)
        if not rematch:
            raise RuntimeError('Expected boundbox table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        # Only 5 lines
        try:
            hasbox = float(fileobj.readline().strip())
            angle = float(fileobj.readline().strip())
            a = float(fileobj.readline().strip())
            b = float(fileobj.readline().strip())
            c = float(fileobj.readline().strip())
        except ValueError:
            raise RuntimeError('Error processing boundbox table entries')
        # Get the child sequence entry
        line = fileobj.readline()
        rematch = AmberOFFLibrary._sec4re.match(line)
        if not rematch:
            raise RuntimeError('Expected childsequence table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        n = int(fileobj.readline().strip())
        if isinstance(templ, ResidueTemplate) and n != 2:
            raise RuntimeError('child sequence for single residue must be 2')
        elif not isinstance(templ, ResidueTemplate) and n != len(templ) + 1:
            raise RuntimeError('child sequence must be 1 greater than the '
                               'number of residues in the unit')
        # Get the CONNECT array to set head and tail
        line = fileobj.readline()
        rematch = AmberOFFLibrary._sec5re.match(line)
        if not rematch:
            raise RuntimeError('Expected connect array not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        try:
            head = int(fileobj.readline().strip())
            tail = int(fileobj.readline().strip())
        except ValueError:
            raise RuntimeError('Error processing connect table entries')
        if head > 0:
            templ.head = templ[head-1]
        if tail > 0:
            templ.tail = templ[tail-1]
        # Get the connectivity array to set bonds
        line = fileobj.readline()
        rematch = AmberOFFLibrary._sec6re.match(line)
        if not rematch:
            raise RuntimeError('Expected connectivity table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        line = fileobj.readline()
        while line[0] != '!':
            i, j, flag = line.split()
            line = fileobj.readline()
            templ.add_bond(int(i)-1, int(j)-1)
        # Get the hierarchy table
        rematch = AmberOFFLibrary._sec7re.match(line)
        if not rematch:
            raise RuntimeError('Expected hierarchy table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        line = fileobj.readline()
        while line[0] != '!':
            # Skip this section... not used
            # TODO turn this into a sanity check
            line = fileobj.readline()
        # Get the unit name
        rematch = AmberOFFLibrary._sec8re.match(line)
        if not rematch:
            raise RuntimeError('Expected unit name string not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        fileobj.readline() # Skip this... not used
        line = fileobj.readline()
        # Get the atomic positions
        rematch = AmberOFFLibrary._sec9re.match(line)
        if not rematch:
            raise RuntimeError('Expected unit positions table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        for atom in templ:
            x, y, z = fileobj.readline().split()
            atom.xx, atom.xy, atom.xz = float(x), float(y), float(z)
        line = fileobj.readline()
        # Get the residueconnect table
        rematch = AmberOFFLibrary._sec10re.match(line)
        if not rematch:
            raise RuntimeError('Expected unit residueconnect table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        for i in xrange(nres):
            c1,c2,c3,c4,c5,c6 = [int(x) for x in fileobj.readline().split()]
            if templ.head is not None and templ.head is not templ[c1-1]:
                warnings.warn('HEAD atom is not connect0')
            if templ.tail is not None and templ.tail is not templ[c2-1]:
                warnings.warn('TAIL atom is not connect1')
            for i in (c3, c4, c5, c6):
                if i == 0: continue
                templ.connections.append(templ[i-1])
        # Get the residues table
        line = fileobj.readline()
        rematch = AmberOFFLibrary._sec11re.match(line)
        if not rematch:
            raise RuntimeError('Expected unit residues table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        for i in xrange(nres):
            resname, id, next, start, typ, img = fileobj.readline().split()
            resname = _strip_enveloping_quotes(resname)
            id = int(id)
            start = int(start)
            next = int(next)
            typ = _strip_enveloping_quotes(typ)
            img = int(img)
            if next - start != len(templ):
                raise RuntimeError('residues table predicted %d, not %d atoms' %
                                   (next-start, len(templ)))
        # Get the residues sequence table
        line = fileobj.readline()
        rematch = AmberOFFLibrary._sec12re.match(line)
        if not rematch:
            raise RuntimeError('Expected residue sequence number not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        for i in xrange(nres):
            #TODO sanity check
            fileobj.readline()
        line = fileobj.readline()
        # Get the solventcap array
        rematch = AmberOFFLibrary._sec13re.match(line)
        if not rematch:
            raise RuntimeError('Expected unit solventcap array not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        # Ignore the solvent cap
        fileobj.readline()
        fileobj.readline()
        fileobj.readline()
        fileobj.readline()
        fileobj.readline()
        # Velocities
        line = fileobj.readline()
        rematch = AmberOFFLibrary._sec14re.match(line)
        if not rematch:
            raise RuntimeError('Expected unit solventcap array not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        for atom in templ:
            vx, vy, vz = [float(x) for x in fileobj.readline().split()]
            atom.vx, atom.vy, atom.vz = vx, vy, vz

        return templ

# Helper routines
def _strip_enveloping_quotes(inp):
    """ Strips the quotation marks enveloping a string """
    if inp[0] == inp[-1] == '"' or inp[0] == inp[-1] == "'":
        return inp[1:-1]
    return inp
