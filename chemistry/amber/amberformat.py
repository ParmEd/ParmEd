"""
This is a generalization of the readparm.AmberParm class to handle similar
Amber-style files with %FLAG/%FORMAT tags
"""
from __future__ import division

from chemistry.amber.constants import (NATOM, NTYPES, NBONH, NTHETH, NPHIH,
            NEXT, NRES, NBONA, NTHETA, NPHIA, NUMBND, NUMANG, NPTRA, NATYP,
            NPHB, IFBOX, IFCAP, AMBER_ELECTROSTATIC)
from chemistry.exceptions import AmberFormatWarning, FlagError
from copy import copy
import datetime
from math import ceil
import re
from warnings import warn

# Some Py3 compatibility tweaks
if not 'unicode' in dir(__builtins__): unicode = str
if not 'basestring' in dir(__builtins__): basestring = str

class FortranFormat(object):
    """ Handles fortran formats """

    #===================================================

    def __init__(self, format_string, strip_strings=True):
        """
        Sets the format string and determines how we will read and write
        strings using this format
        """
        self.format = format_string
        self.strip_strings = strip_strings # for ease of copying

        # Define a function that processes all arguments prior to adding them to
        # the returned list. By default, do nothing, but this allows us to
        # optionally strip whitespace from strings.
        self.process_method = lambda x: x

        if 'a' in format_string.lower():
            # replace our write() method with write_string to force left-justify
            self.type, self.write = str, self.write_string
            try:
                self.nitems, self.itemlen = format_string.lower().split('a')
                self.nitems, self.itemlen = int(self.nitems), int(self.itemlen)
            except ValueError:
                self.nitems, self.itemlen = 1, 80
            self.fmt = '%s'
            # See if we want to strip the strings
            if strip_strings: self.process_method = lambda x: x.strip()

        elif 'I' in format_string.upper():
            self.type = int
            if len(format_string.upper().split('I')[0]) == 0:
                self.nitems = 1
            else:
                self.nitems = int(format_string.upper().split('I')[0])
            self.itemlen = int(format_string.upper().split('I')[1])
            self.fmt = '%%%dd' % self.itemlen

        elif 'E' in format_string.upper():
            self.type = float
            format_parts = format_string.upper().split('E')
            if len(format_parts[0]) == 0:
                self.nitems = 1
            else:
                self.nitems = int(format_parts[0])
            self.itemlen = int(format_parts[1].split('.')[0])
            self.num_decimals = int(format_parts[1].split('.')[1])
            self.fmt = '%%%s.%sE' % (self.itemlen, self.num_decimals)

        elif 'F' in format_string.upper():
            self.type = float
            # Strip out any parentheses
            format_string = format_string.replace('(', '').replace(')', '')
            format_parts = format_string.upper().split('F')
            if len(format_parts[0].strip()) == 0:
                self.nitems = 1
            else:
                self.nitems = int(format_parts[0])
            self.itemlen = int(format_parts[1].split('.')[0])
            self.num_decimals = int(format_parts[1].split('.')[1])
            self.fmt = '%%%s.%sF' % (self.itemlen, self.num_decimals)

        else:
            # replace our write() method with write_string to force left-justify
            self.type, self.write = str, self.write_string
            warn('Unrecognized format "%s". Assuming string.' % format_string,
                 AmberFormatWarning)
            self.fmt = '%s'
            self.nitems, self.itemlen = 1, 80
            # See if we want to strip the strings
            if strip_strings: self.process_method = lambda x: x.strip()

    #===================================================

    def __copy__(self):
        return type(self)(self.format, self.strip_strings)

    #===================================================

    def __str__(self):
        return self.format

    #===================================================

    def write(self, items, dest):
        """ Writes a list/tuple of data (or a single item) """
        if hasattr(items, '__iter__') and not isinstance(items, basestring):
            mod = self.nitems - 1
            for i, item in enumerate(items):
                dest.write(self.fmt % item)
                if i % self.nitems == mod:
                    dest.write('\n')
            if i % self.nitems != mod:
                dest.write('\n')
        else:
            dest.write(self.fmt % item)
            dest.write('\n')

    #===================================================

    def write_string(self, items, dest):
        """ Writes a list/tuple of strings """
        if hasattr(items, '__iter__') and not isinstance(items, basestring):
            mod = self.nitems - 1
            for i, item in enumerate(items):
                dest.write((self.fmt % item).ljust(self.itemlen))
                if i % self.nitems == mod:
                    dest.write('\n')
            if i % self.nitems != mod:
                dest.write('\n')
        else:
            dest.write((self.fmt % item).ljust(self.itemlen))
            dest.write('\n')

    #===================================================

    def read_nostrip(self, line):
        """
        Reads the line and returns converted data. Special-cased for flags that
        may contain 'blank' data. ugh.
        """
        line = line.rstrip('\n')
        nitems = int(ceil(len(line) / self.itemlen))
        ret = [0 for i in xrange(nitems)]
        start, end = 0, self.itemlen
        for i in xrange(nitems):
            ret[i] = self.process_method(self.type(line[start:end]))
            start = end
            end += self.itemlen
        return ret

    #===================================================

    def read(self, line):
        """ Reads the line and returns the converted data """
        line = line.rstrip()
        nitems = int(ceil(len(line) / self.itemlen))
        ret = [0 for i in xrange(nitems)]
        start, end = 0, self.itemlen
        for i in xrange(nitems):
            ret[i] = self.process_method(self.type(line[start:end]))
            start = end
            end += self.itemlen
        return ret

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AmberFormat(object):
    """ 
    Generalization of the AmberParm class without some of the assumptions made
    about Amber topology files specifically
    """
   
    CHARGE_SCALE = AMBER_ELECTROSTATIC # chamber uses a SLIGHTLY diff value

    #===================================================

    def __init__(self, fname=None):
        """ Constructor.  Read a file if given """
        self._ncopies = 0
        self.parm_data = {}
        self.formats = {}
        self.parm_comments = {}
        self.flag_list = []
        self.version = None
        self.prm_name = fname
        self.charge_flag = 'CHARGE'
        self.valid = True

        if fname is not None:
            self.rdparm(fname)

    #===================================================

    def __copy__(self):
        """ Copy all of the data """
        self._ncopies += 1
        other = type(self)()
        other.flag_list = self.flag_list[:]
        other.version = self.version
        other.prm_name = self.prm_name + '_copy%d' % self._ncopies
        other.charge_flag = self.charge_flag
        other.valid = self.valid
        other.parm_data = {}
        other.parm_comments = {}
        other.formats = {} # formats{} are copied shallow
        for flag in other.flag_list:
            other.parm_data[flag] = self.parm_data[flag][:]
            other.parm_comments[flag] = self.parm_comments[flag][:]
            other.formats[flag] = copy(self.formats[flag])
        return other

    #===================================================

    def view(self, cls):
        """
        Returns a view of the current object as another object.

        Parameters:
            cls Class definition of an AmberParm subclass for the current
                object to be converted into

        Returns:
            instance of cls initialized from data in this object. This is NOT a
            deep copy, so modifying the original object may modify this. The
            copy function will create a deep copy of any AmberFormat-derived
            object
        """
        # If these are the same classes, just return the original instance,
        # since there's nothing to do. Classes are singletons, so use "is"
        if type(self) is cls:
            return self
        if hasattr(cls, 'load_from_rawdata'):
            return cls.load_from_rawdata(self)
        raise ValueError('Cannot instantiate %s from AmberFormat' %
                         cls.__name__)

    #===================================================

    def rdparm(self, fname, slow=False):
        """ Parses the Amber format file """
        # See if we have the optimized parser available

        self.prm_name = fname
        self.version = None # reset all top info each time rdparm is called
        self.formats = {}
        self.parm_data = {}
        self.parm_comments = {}
        self.flag_list = []
        self.valid = False

        try:
            from chemistry.amber import _rdparm
        except ImportError:
            return self.rdparm_slow(fname)

        if slow:
            return self.rdparm_slow(fname)

        # We have the optimized version
        try:
            ret = _rdparm.rdparm(fname)
        except TypeError:
            # This is raised if VERSION is not found
            raise
            return self.rdparm_old(open(fname, 'r').readlines())
        else:
            # Unpack returned contents
            parm_data, parm_comments, formats, unkflg, flag_list, version = ret
            # Now assign them to instance attributes and process where necessary
            self.parm_data = parm_data
            self.parm_comments = parm_comments
            for key in formats:
                self.formats[key] = FortranFormat(formats[key])
            self.flag_list = flag_list
            self.version = version
            # Now we have to process all of those sections that the optimized
            # parser couldn't figure out
            for flag in unkflg:
                rawdata = self.parm_data[flag]
                self.parm_data[flag] = []
                for line in rawdata:
                    self.parm_data[flag].extend(self.formats[flag].read(line))

            try:
                for i, chg in enumerate(self.parm_data[self.charge_flag]):
                    self.parm_data[self.charge_flag][i] = chg / self.CHARGE_SCALE
            except KeyError:
                pass
            self.valid = True

    #===================================================

    def rdparm_slow(self, fname):
        """
        Parses the Amber format file. This parser is written in pure Python and
        is therefore slower than the C++-optimized version
        """

        current_flag = ''
        gathering_data = False
        fmtre = re.compile(r'%FORMAT *\((.+)\)')

        # Open up the file and read the data into memory
        prm = open(self.prm_name, 'r')

        for line in prm:

            if line[0:8] == '%VERSION':
                self.version = line.strip()

            elif line[0:5] == '%FLAG':
                current_flag = line[6:].strip()
                self.formats[current_flag] = ''
                self.parm_data[current_flag] = []
                self.parm_comments[current_flag] = []
                self.flag_list.append(current_flag)
                gathering_data = False

            elif line[0:8] == '%COMMENT':
                self.parm_comments[current_flag].append(line[9:].strip())

            elif line[0:7] == '%FORMAT':
                fmt = FortranFormat(fmtre.match(line).groups()[0])
                # RESIDUE_ICODE can have a lot of blank data...
                if current_flag == 'RESIDUE_ICODE':
                    fmt.read = fmt.read_nostrip
                self.formats[current_flag] = fmt
                gathering_data = True

            elif gathering_data:
                self.parm_data[current_flag].extend(fmt.read(line))

        # convert charges to fraction-electrons
        try:
            for i, chg in enumerate(self.parm_data[self.charge_flag]):
                self.parm_data[self.charge_flag][i] = chg / self.CHARGE_SCALE
        except KeyError:
            pass
        self.valid = True

        prm.close()

        # If we don't have a version, then read in an old-file topology
        if self.version is None:
            self.rdparm_old(open(self.prm_name, 'r').readlines())

        return

    #===================================================

    def rdparm_old(self, prmtop_lines):
        """
        This reads an old-style topology file and stores the results in the
        same data structures as a new-style topology file
        """
        def read_integer(line_idx, lines, num_items):
            # line_idx should be the line _before_ the first line you
            # want data from.
            i, tmp_data = 0, []
            while i < num_items:
                idx = i % 12
                if idx == 0:
                    line_idx += 1
                try:
                    tmp_data.append(int(lines[line_idx][idx*6:idx*6+6]))
                except ValueError:
                    raise ValueError(
                            'Error parsing line %d, token %d [%s]: Problem '
                            'during integer read.' % (line_idx, idx,
                            lines[line_idx][idx*6:idx*6+6])
                    )
                i += 1
            # If we had no items, we need to jump a line:
            if num_items == 0: line_idx += 1
            return tmp_data, line_idx

        def read_string(line_idx, lines, num_items):
            # line_idx should be the line _before_ the first line you
            # want data from.
            i, tmp_data = 0, []
            while i < num_items:
                idx = i % 20
                if idx == 0:
                    line_idx += 1
                tmp_data.append(lines[line_idx][idx*4:idx*4+4])
                i += 1
            # If we had no items, we need to jump a line:
            if num_items == 0: line_idx += 1
            return tmp_data, line_idx

        def read_float(line_idx, lines, num_items):
            # line_idx should be the line _before_ the first line you
            # want data from.
            i, tmp_data = 0, []
            while i < num_items:
                idx = i % 5
                if idx == 0:
                    line_idx += 1
                try:
                    tmp_data.append(float(lines[line_idx][idx*16:idx*16+16]))
                except ValueError:
                    raise ValueError(
                            'Error parsing line %d, token %d [%s]: Problem '
                            'during floating point  read.' % (line_idx, idx,
                            lines[line_idx][idx*16:idx*16+16])
                    )
                i += 1
            # If we had no items, we need to jump a line:
            if num_items == 0: line_idx += 1
            return tmp_data, line_idx

        # First add a title
        self.addFlag('TITLE', '20a4', data=['| Converted old-style topology'])

        # Next, read in the pointers
        line_idx = 0
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 30)
        # Add a final pointer of 0, which corresponds to NUMEXTRA
        tmp_data.append(0)
        self.addFlag('POINTERS', '10I8', data=tmp_data)

        # Set some of the pointers we need
        natom = self.parm_data['POINTERS'][NATOM]
        ntypes = self.parm_data['POINTERS'][NTYPES]
        nres = self.parm_data['POINTERS'][NRES]
        numbnd = self.parm_data['POINTERS'][NUMBND]
        numang = self.parm_data['POINTERS'][NUMANG]
        nptra = self.parm_data['POINTERS'][NPTRA]
        natyp = self.parm_data['POINTERS'][NATYP]
        nbonh = self.parm_data['POINTERS'][NBONH]
        nbona = self.parm_data['POINTERS'][NBONA]
        ntheth = self.parm_data['POINTERS'][NTHETH]
        ntheta = self.parm_data['POINTERS'][NTHETA]
        nex = self.parm_data['POINTERS'][NEXT]
        nphia = self.parm_data['POINTERS'][NPHIA]
        nphb = self.parm_data['POINTERS'][NPHB]
        nphih = self.parm_data['POINTERS'][NPHIH]

        # Next read in the atom names
        tmp_data, line_idx = read_string(line_idx, prmtop_lines, natom)
        self.addFlag('ATOM_NAME', '20a4', data=tmp_data)

        # Next read the charges
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, natom)
        # Divide by the electrostatic constant
        tmp_data = [x / self.CHARGE_SCALE for x in tmp_data]
        self.addFlag('CHARGE', '5E16.8', data=tmp_data)

        # Next read the masses
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, natom)
        self.addFlag('MASS', '5E16.8', data=tmp_data)

        # Next read atom type index
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, natom)
        self.addFlag('ATOM_TYPE_INDEX', '10I8', data=tmp_data)

        # Next read number excluded atoms
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, natom)
        self.addFlag('NUMBER_EXCLUDED_ATOMS', '10I8', data=tmp_data)
      
        # Next read nonbonded parm index
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, ntypes**2)
        self.addFlag('NONBONDED_PARM_INDEX', '10I8', data=tmp_data)

        # Next read residue label
        tmp_data, line_idx = read_string(line_idx, prmtop_lines, nres)
        self.addFlag('RESIDUE_LABEL', '20a4', data=tmp_data)

        # Next read residue pointer
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nres)
        self.addFlag('RESIDUE_POINTER', '10I8', data=tmp_data)
   
        # Next read bond force constant
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, numbnd)
        self.addFlag('BOND_FORCE_CONSTANT', '5E16.8', data=tmp_data)

        # Next read bond equil value
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, numbnd)
        self.addFlag('BOND_EQUIL_VALUE', '5E16.8', data=tmp_data)

        # Next read angle force constant
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, numang)
        self.addFlag('ANGLE_FORCE_CONSTANT', '5E16.8', data=tmp_data)

        # Next read the angle equilibrium value
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, numang)
        self.addFlag('ANGLE_EQUIL_VALUE', '5E16.8', data=tmp_data)

        # Next read the dihedral force constant
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, nptra)
        self.addFlag('DIHEDRAL_FORCE_CONSTANT', '5E16.8', data=tmp_data)

        # Next read dihedral periodicity
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, nptra)
        self.addFlag('DIHEDRAL_PERIODICITY', '5E16.8', data=tmp_data)

        # Next read the dihedral phase
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, nptra)
        self.addFlag('DIHEDRAL_PHASE', '5E16.8', data=tmp_data)

        # Next read SOLTY (?)
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, natyp)
        self.addFlag('SOLTY', '5E16.8', data=tmp_data)

        # Next read lennard jones acoef and bcoef
        numvals = ntypes * (ntypes + 1) / 2
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, numvals)
        self.addFlag('LENNARD_JONES_ACOEF', '5E16.8', data=tmp_data)
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, numvals)
        self.addFlag('LENNARD_JONES_BCOEF', '5E16.8', data=tmp_data)

        # Next read bonds including hydrogen
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nbonh*3)
        self.addFlag('BONDS_INC_HYDROGEN', '10I8', data=tmp_data)

        # Next read bonds without hydrogen
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nbona*3)
        self.addFlag('BONDS_WITHOUT_HYDROGEN', '10I8', data=tmp_data)

        # Next read angles including hydrogen
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, ntheth*4)
        self.addFlag('ANGLES_INC_HYDROGEN', '10I8', data=tmp_data)

        # Next read angles without hydrogen
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, ntheta*4)
        self.addFlag('ANGLES_WITHOUT_HYDROGEN', '10I8', data=tmp_data)

        # Next read dihdrals including hydrogen
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nphih*5)
        self.addFlag('DIHEDRALS_INC_HYDROGEN', '10I8', data=tmp_data)

        # Next read dihedrals without hydrogen
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nphia*5)
        self.addFlag('DIHEDRALS_WITHOUT_HYDROGEN', '10I8', data=tmp_data)

        # Next read the excluded atoms list
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nex)
        self.addFlag('EXCLUDED_ATOMS_LIST', '10I8', data=tmp_data)

        # Next read the hbond terms
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, nphb)
        self.addFlag('HBOND_ACOEF', '5E16.8', data=tmp_data)
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, nphb)
        self.addFlag('HBOND_BCOEF', '5E16.8', data=tmp_data)
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, nphb)
        self.addFlag('HBCUT', '5E16.8', data=tmp_data)

        # Next read amber atom type
        tmp_data, line_idx = read_string(line_idx, prmtop_lines, natom)
        self.addFlag('AMBER_ATOM_TYPE', '20a4', data=tmp_data)

        # Next read tree chain classification
        tmp_data, line_idx = read_string(line_idx, prmtop_lines, natom)
        self.addFlag('TREE_CHAIN_CLASSIFICATION', '20a4', data=tmp_data)

        # Next read the join array
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, natom)
        self.addFlag('JOIN_ARRAY', '10I8', data=tmp_data)

        # Next read the irotat array
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, natom)
        self.addFlag('IROTAT', '10I8', data=tmp_data)

        # Now do PBC stuff
        if self.parm_data['POINTERS'][IFBOX]:
            # Solvent pointers
            tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 3)
            self.addFlag('SOLVENT_POINTERS', '10I8', data=tmp_data)
            nspm = tmp_data[1]

            # Atoms per molecule
            tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nspm)
            self.addFlag('ATOMS_PER_MOLECULE', '10I8', data=tmp_data)

            # Box dimensions
            tmp_data, line_idx = read_float(line_idx, prmtop_lines, 4)
            self.addFlag('BOX_DIMENSIONS', '5E16.8', data=tmp_data)

        # Now do CAP stuff
        if self.parm_data['POINTERS'][IFCAP]:
            # CAP_INFO
            tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 1)
            self.addFlag('CAP_INFO', '10I8', data=tmp_data)
            tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 4)
            self.addFlag('CAP_INFO2', '10I8', data=tmp_data)
        # end if self.parm_data['POINTERS'][IFCAP]

    #===================================================

    def set_version(self):
        """ Sets the version string """
        now = datetime.datetime.now()
        self.version = (
                '%%VERSION  VERSION_STAMP = V0001.000  DATE = %02d/%02d/%02d  '
                '%02d:%02d:%02d' % (now.month, now.day, now.year % 100,
                                    now.hour, now.minute, now.second)
        )

    #===================================================

    def writeParm(self, name):
        """
        Writes the current data in parm_data into a new topology file with
        the given name
        """
        # now that we know we will write the new prmtop file, open the new file
        new_prm = open(name, 'w')

        # get current time to put into new prmtop file if we had a %VERSION
        self.set_version()

        # convert charges back to amber charges...
        if self.charge_flag in self.parm_data.keys():
            for i in xrange(len(self.parm_data[self.charge_flag])):
                self.parm_data[self.charge_flag][i] *= self.CHARGE_SCALE

        # write version to top of prmtop file
        new_prm.write('%s\n' % self.version)

        # write data to prmtop file, inserting blank line if it's an empty field
        for i in xrange(len(self.flag_list)):
            flag = self.flag_list[i]
            new_prm.write('%%FLAG %s\n' % flag)
            # Insert any comments before the %FORMAT specifier
            for comment in self.parm_comments[flag]:
                new_prm.write('%%COMMENT %s\n' % comment)
            new_prm.write('%%FORMAT(%s)\n' % self.formats[flag])
            if len(self.parm_data[flag]) == 0: # empty field...
                new_prm.write('\n')
                continue
            self.formats[flag].write(self.parm_data[flag], new_prm)

        new_prm.close() # close new prmtop

        if self.charge_flag in self.parm_data.keys():
            # Convert charges back to electron-units
            for i in xrange(len(self.parm_data[self.charge_flag])):
                self.parm_data[self.charge_flag][i] /= self.CHARGE_SCALE

    #===================================================

    def addFlag(self, flag_name, flag_format, data=None, num_items=-1,
                comments=[], after=None):
        """
        Adds a new flag with the given flag name and Fortran format string and
        initializes the array with the values given, or as an array of 0s
        of length num_items

        Parameters:
            flag_name (str): Name of the flag to insert. It is converted to all
                             upper case
            flag_format (str): Fortran Format statement (do NOT enclose in ())
            data (list): Sequence with data for the new flag. If None, a list
                         of zeros of length num_items is given as a holder
            num_items (int): Number of items in the section (only used if data
                             is None)
            comments (list): List of comments to put in this section
            after (str): If not None, this flag will be added after the given
                         flag. If the 'after' flag does not exist, IndexError is
                         raised.
        """
        if after is not None:
            after = after.upper()
            if not after in self.flag_list:
                raise IndexError('%s not found in topology flag list' % after)
            # If the 'after' flag is the last one, just append
            if self.flag_list[-1] == after:
                self.flag_list.append(flag_name.upper())
            # Otherwise find the index and add it after
            idx = self.flag_list.index(after) + 1
            self.flag_list.insert(idx, flag_name.upper())
        else:
            self.flag_list.append(flag_name.upper())
        self.formats[flag_name.upper()] = FortranFormat(flag_format)
        if data is not None:
            self.parm_data[flag_name.upper()] = list(data)
        else:
            if num_items < 0:
                raise FlagError("If you do not supply prmtop data, num_items "
                                "must be non-negative!")
            self.parm_data[flag_name.upper()] = [0 for i in xrange(num_items)]
        if comments:
            if isinstance(comments, str) or isinstance(comments, unicode):
                comments = [comments]
            elif isinstance(comments, tuple):
                comments = list(comments)
            elif isinstance(comments, list):
                pass
            else:
                raise TypeError('Comments must be string, list, or tuple')
            self.parm_comments[flag_name.upper()] = comments
        else:
            self.parm_comments[flag_name.upper()] = []

    #===================================================

    def deleteFlag(self, flag_name):
        """ Removes a flag from the topology file """
        flag_name = flag_name.upper()
        if not flag_name in self.flag_list:
            return # already gone
        del self.flag_list[self.flag_list.index(flag_name)]
        del self.parm_comments[flag_name]
        del self.formats[flag_name]
        del self.parm_data[flag_name]

