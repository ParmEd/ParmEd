"""
This module contains classes for reading various TINKER-style files
"""
from compat24 import all
from chemistry.exceptions import TinkerKeyFileError, TinkerDynFileError

class KeywordControlFile(object):
    """ Reads and processes a keyword control file for TINKER simulations """
    _datatypes = {'PARAMETERS' : str, 'A-AXIS' : float, 'B-AXIS' : float,
                  'C-AXIS' : float, 'ALPHA' : float, 'BETA' : float,
                  'GAMMA' : float}

    def __init__(self, fname):
        # The currently recognized/parsed keywords (these are the only ones
        # necessary for running Amoeba in Amber)
        self.keywords = {'PARAMETERS' : None, 'A-AXIS' : None, 'B-AXIS' : None,
                         'C-AXIS' : None, 'ALPHA' : None, 'BETA' : None,
                         'GAMMA' : None}
        # Parse the file
        for line in open(fname, 'r'):
            # Skip over any blank lines
            if not line.strip(): continue
            words = line.split()
            key = words[0].upper()
            # Get rid of the keyword and all whitespace
            result = line.replace(words[0], '').strip()
            try:
                self.keywords[key] = self._datatypes[key](result)
            except KeyError:
                pass # All non-keyword lines are ignored comments
            except ValueError:
                raise TinkerKeyFileError(
                        'Malformed keyword control file! Could not convert the '
                        'value of %s into type %s' % (key, self._datatypes[key])
                )

class XyzFile(object):
    """ Reads and processes a Tinker XYZ file """
    class _Atom(object):
        """
        An atom object that stores the atomic information stored in the XyzFile
        """
        def __init__(self, name, x, y, z, type, bonded_partners):
            self.name = str(name)
            self.position = [float(x), float(y), float(z)]
            self.type = int(type)
            # index from 1
            self.bonded_partners = [int(i) for i in bonded_partners]

    class _AtomList(list):
        " A list of _Atom objects "
        def __init__(self):
            super(XyzFile._AtomList, self).__init__()

        def add(self, name, x, y, z, type, bonded_partners):
            self.append(XyzFile._Atom(name, x, y, z, type, bonded_partners))

        def append(self, thing):
            if not isinstance(thing, XyzFile._Atom):
                raise TypeError('XyzFile._AtomList can only append _Atom\'s')
            super(XyzFile._AtomList, self).append(thing)

        def extend(self, things):
            for thing in things: self.append(thing)

    def __init__(self, fname):
        self.natom = 0
        # Create the list of atoms
        self.atom_list = XyzFile._AtomList()
        f = open(fname, 'r')
        for line in f:
            if self.natom == 0:
                self.natom = int(line.strip())
                # Set up blank positions and connections arrays
                self.connections = [[] for i in xrange(self.natom)]
                continue
            words = line.split()
            if len(words) == 6 and all([is_float(x) for x in words]):
                # The first line after the number of atoms _could_ be the box
                # information. So capture that here and store the info
                self.box = [float(x) for x in words]
                continue
            self.atom_list.add(words[1], words[2], words[3],
                               words[4], words[5], words[6:])
        f.close()

class DynFile(object):
    """ Reads and processes a Tinker DYN file """
    def __init__(self, fname=None):
        if fname is not None:
            self.read(fname)

    def read(self, fname):
        """ Parses the .dyn file """
        f = open(fname, 'r')
        try:
            if f.readline().strip() != 'Number of Atoms and Title :':
                raise TinkerDynFileError('%s is not a recognized TINKER '
                                         '.dyn file' % fname)
            line = f.readline()
            self.natom, self.title = int(line[0:8].strip()), line[8:].strip()
            if f.readline().strip() != 'Periodic Box Dimensions :':
                raise TinkerDynFileError('No periodic box dimension line')
            self.box = [0.0 for i in xrange(6)]
            words = f.readline().upper().split()
            self.box[0] = float(words[0].replace('D', 'E'))
            self.box[1] = float(words[1].replace('D', 'E'))
            self.box[2] = float(words[2].replace('D', 'E'))
            words = f.readline().upper().split()
            self.box[3] = float(words[0].replace('D', 'E'))
            self.box[4] = float(words[1].replace('D', 'E'))
            self.box[5] = float(words[2].replace('D', 'E'))
            if f.readline().strip() != 'Current Atomic Positions :':
                raise TinkerDynFileError('No atomic positions in .dyn file')
            self.positions = [[0.0, 0.0, 0.0] for i in xrange(self.natom)]
            DynFile._read_section(f, self.positions, self.natom)
            line = f.readline().strip()
            if line == 'Current Translational Velocities :':
                self.rigidbody = True
                self.translational_velocities = [[0.0, 0.0, 0.0] for i in
                            xrange(self.natom)]
                DynFile._read_section(f, self.translational_velocities,
                                      self.natom)
                if f.readline().strip() != 'Current Angular Velocities :':
                    raise TinkerDynFileError('Could not find Angular velocity '
                            'section in .dyn file')
                self.angular_velocities = [[0.0, 0.0, 0.0]
                                            for i in xrange(self.natom)]
                DynFile._read_section(f, self.angular_velocities, self.natom)
                if f.readline().strip() != 'Current Angular Momenta :':
                    raise TinkerDynFileError('Could not find angular momenta '
                                             'section in .dyn file')
                self.angular_momenta = [[0.0, 0.0, 0.0]
                                        for i in xrange(self.natom)]
                DynFile._read_section(f, self.angular_momenta, self.natom)
            elif line == 'Current Atomic Velocities :':
                self.rigidbody = False
                self.velocities = [[0.0, 0.0, 0.0] for i in xrange(self.natom)]
                DynFile._read_section(f, self.velocities, self.natom)
                if f.readline().strip() != 'Current Atomic Accelerations :':
                    raise TinkerDynFileError('Could not find accelerations '
                                             'in %s ' % fname)
                self.accelerations = [[0.0, 0.0, 0.0]
                                      for i in xrange(self.natom)]
                DynFile._read_section(f, self.accelerations, self.natom)
                if f.readline().strip() != 'Alternate Atomic Accelerations :':
                    raise TinkerDynFileError('Could not find old '
                                             'accelerations in %s' % fname)
                self.old_accelerations = [[0.0, 0.0, 0.0] 
                                          for i in xrange(self.natom)]
                DynFile._read_section(f, self.old_accelerations, self.natom)
            else:
                raise TinkerDynFileError('No velocities in %s' % fname)
        finally:
            f.close()

    @staticmethod
    def _read_section(f, container, natom):
        """
        Reads a section of an open file into the passed container. The
        container should be a 2-dimensional list of length natom x 3
        """
        for i in xrange(natom):
            words = f.readline().upper().split()
            try:
                container[i][0] = float(words[0].replace('D', 'E'))
                container[i][1] = float(words[1].replace('D', 'E'))
                container[i][2] = float(words[2].replace('D', 'E'))
            except (IndexError, ValueError):
                raise TinkerDynFileError('Could not parse values from dyn file')

def is_float(thing):
    try:
        float(thing)
        return True
    except ValueError:
        return False
