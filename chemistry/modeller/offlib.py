"""
Tool for parsing and writing OFF library files to and from dictionaries of
ResidueTemplate objects
"""
import compat24

from chemistry import Atom, Bond
from chemistry.constants import RAD_TO_DEG
from chemistry.exceptions import AmberOFFWarning
from chemistry.modeller.residue import ResidueTemplate, ResidueTemplateContainer
from chemistry.modeller.residue import PROTEIN, NUCLEIC, SOLVENT, UNKNOWN
from collections import OrderedDict

import re
import warnings

class AmberOFFLibrary(object):
    """
    Class containing static methods responsible for parsing and writing OFF
    libraries
    """

    #===================================================

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

    #===================================================

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

        if own_handle: fileobj.close()

        return residues

    #===================================================

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
        container = ResidueTemplateContainer()
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
            if resx == nres + 1:
                container.append(templ)
                nres += 1
                templ = ResidueTemplate(name)
            templ.add_atom(atom)
            line = fileobj.readline()
        container.append(templ)
        if nres > 1:
            start_atoms = []
            runsum = 0
            for res in container:
                start_atoms.append(runsum)
                runsum += len(res)
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
        else:
            if hasbox > 0:
                angle *= RAD_TO_DEG
                container.box = [a, b, c, angle, angle, angle]
        # Get the child sequence entry
        line = fileobj.readline()
        rematch = AmberOFFLibrary._sec4re.match(line)
        if not rematch:
            raise RuntimeError('Expected childsequence table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError('Found residue %s while processing residue %s' %
                               (rematch.groups()[0], name))
        n = int(fileobj.readline().strip())
        if nres + 1 != n:
            warnings.warn('Unexpected childsequence (%d); expected %d for '
                          'residue %s' % (n, nres+1, name), AmberOFFWarning)
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
        if head > 0 and nres == 1:
            templ.head = templ[head-1]
        elif head > 0 and nres > 1:
            if head < sum([len(r) for r in container]):
                raise RuntimeError('HEAD on multi-residue unit not supported')
        if tail > 0 and nres == 1:
            templ.tail = templ[tail-1]
        elif tail > 0 and nres > 1:
            if tail < sum([len(r) for r in container]):
                warnings.warn('TAIL on multi-residue unit not supported (%s). '
                              'Ignored...' % name, AmberOFFWarning)
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
            if nres > 1:
                # Find which residue we belong in
                i = int(i) - 1
                j = int(j) - 1
                for ii, idx in enumerate(start_atoms):
                    if idx > i:
                        ii -= 1
                        break
                start_idx = start_atoms[ii]
                container[ii].add_bond(i-start_idx, j-start_idx)
            else:
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
        for res in container:
            for atom in res:
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
            if next - start != len(container[i]):
                warnings.warn('residue table predicted %d, not %d atoms for '
                              'residue %s' % (next-start, len(container[i]),
                              name), AmberOFFWarning)
            if typ == 'p':
                container[i].type = PROTEIN
            elif typ == 'n':
                container[i].type = NUCLEIC
            elif typ == 'w':
                container[i].type = SOLVENT
            elif typ != '?':
                warnings.warn('Unknown residue type "%s"' % typ,
                              AmberOFFWarning)
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
        for res in container:
            for atom in res:
                vx, vy, vz = [float(x) for x in fileobj.readline().split()]
                atom.vx, atom.vy, atom.vz = vx, vy, vz

        if nres > 1:
            return container
        return templ

    #===================================================

    @staticmethod
    def write(lib, dest):
        """ Writes a dictionary of ResidueTemplate units to a file in OFF format

        Parameters
        ----------
        lib : dict {str : :class:`ResidueTemplate`}
            Items can be either :class:`ResidueTemplate` or
            :class:`ResidueTemplateContainer` instances
        dest : str or file-like
            Either a file name or a file-like object to write the file to
        """
        own_handle = False
        if not hasattr(dest, 'write'):
            dest = open(dest, 'w')
            own_handle = True
        # Write the residues in alphabetical order
        names = sorted(lib.keys())
        dest.write('!!index array str\n')
        for name in names:
            dest.write(' "%s"\n' % name)
        for name in names:
            AmberOFFLibrary._write_residue(dest, lib[name])

        if own_handle: dest.close()

    #===================================================

    @staticmethod
    def _write_residue(dest, res):
        """ Writes a residue to an open file handle

        Parameters
        ----------
        dest : file-like
            File object to write the residue information to
        res : :class:`ResidueTemplate`
            The residue template to write to the file
        """
        dest.write('!entry.%s.unit.atoms table  str name  str type  int typex  '
                   'int resx  int flags  int seq  int elmnt  dbl chg\n' %
                   res.name)
        for atom in res:
            dest.write(' "%s" "%s" 0 1 131072 %d %d %.6f\n' % (atom.name,
                       atom.type, atom.idx+1, atom.atomic_number, atom.charge))
        dest.write('!entry.%s.unit.atomspertinfo table  str pname  str ptype  '
                   'int ptypex  int pelmnt  dbl pchg\n' % res.name)
        for atom in res:
            dest.write(' "%s" "%s" 0 -1 0.0\n' % (atom.name, atom.type))
        dest.write('!entry.%s.unit.boundbox array dbl\n' % res.name)
        dest.write((' -1.000000\n' + ' 0.0\n' * 4))
        dest.write('!entry.%s.unit.childsequence single int\n 2\n' % res.name)
        dest.write('!entry.%s.unit.connect array int\n' % res.name)
        if res.head is not None:
            dest.write(' %d\n' % (res.head.idx + 1))
        else:
            dest.write(' 0\n')
        if res.tail is not None:
            dest.write(' %d\n' % (res.tail.idx + 1))
        else:
            dest.write(' 0\n')
        dest.write('!entry.%s.unit.connectivity table  int atom1x  int atom2x  '
                   'int flags\n' % res.name)
        for bond in res.bonds:
            dest.write(' %d %d 1\n' % (bond.atom1.idx+1, bond.atom2.idx+1))
        dest.write('!entry.%s.unit.hierarchy table  str abovetype  int '
                   'abovex  str belowtype  int belowx\n' % res.name)
        dest.write(' "U" 0 "R" 1\n')
        for atom in res:
            dest.write(' "R" 1 "A" %d\n' % (atom.idx + 1))
        dest.write('!entry.%s.unit.name single str\n' % res.name)
        dest.write(' "%s"\n' % res.name)
        dest.write('!entry.%s.unit.positions table  dbl x  dbl y  dbl z\n' %
                   res.name)
        for atom in res:
            dest.write(' %g %g %g\n' % (atom.xx, atom.xy, atom.xz))
        dest.write('!entry.%s.unit.residueconnect table  int c1x  int c2x  '
                   'int c3x  int c4x  int c5x  int c6x\n' % res.name)
        conn = [0, 0, 0, 0, 0, 0]
        if res.head is not None: conn[0] = res.head.idx + 1
        if res.tail is not None: conn[1] = res.tail.idx + 1
        for i, at in enumerate(res.connections):
            conn[i+2] = at.idx + 1
        dest.write(' %d %d %d %d %d %d\n' % tuple(conn))
        dest.write('!entry.%s.unit.residues table  str name  int seq  int '
                   'childseq  int startatomx  str restype  int imagingx\n' %
                   res.name)
        if res.type is PROTEIN:
            typ = 'p'
        elif res.type is NUCLEIC:
            typ = 'n'
        elif res.type is SOLVENT:
            typ='w'
        elif res.type is UNKNOWN:
            typ='?'
        else:
            warnings.warn('Unrecognized residue type %r' % res.type,
                          AmberOFFWarning)
            typ = '?'
        dest.write(' "%s" 1 %d 1 "%s" 0\n' % (res.name, len(res)+1, typ))
        dest.write('!entry.%s.unit.residuesPdbSequenceNumber array int\n 0\n' %
                   res.name)
        dest.write('!entry.%s.unit.solventcap array dbl\n' % res.name)
        dest.write(' -1.000000\n' + ' 0.0\n' * 4)
        dest.write('!entry.%s.unit.velocities table  dbl x  dbl y  dbl z\n' %
                   res.name)
        for atom in res:
            try:
                s = ' %g %g %g\n' % (atom.vx, atom.vy, atom.vz)
            except AttributeError:
                dest.write(' 0.0 0.0 0.0\n')
            else:
                dest.write(s)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Helper routines
def _strip_enveloping_quotes(inp):
    """ Strips the quotation marks enveloping a string """
    if inp[0] == inp[-1] == '"' or inp[0] == inp[-1] == "'":
        return inp[1:-1]
    return inp
