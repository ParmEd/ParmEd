"""
Provides a class for parsing CHARMM-style coordinate files, namely CHARMM .crd
(coordinate) files and CHARMM .rst (restart) file. Uses CharmmFile class in
_charmmfile.py for reading files 

Author: Jason Deckman
Contributors: Jason Swails
Date: Feb. 24, 2015
"""

from chemistry.formats import io
from chemistry.formats.registry import FileFormatType
from chemistry.exceptions import CharmmFileError
from chemistry import unit as u

charmlen = 22
TIMESCALE = 4.888821E-14 * 1e12 # AKMA time units to picoseconds
ONE_TIMESCALE = 1 / TIMESCALE

class CharmmCrdFile(object):
    """
    Reads and parses a CHARMM coordinate file (.crd) into its components,
    namely the coordinates, CHARMM atom types, resid, resname, etc.

    Parameters
    ----------
    fname : str
        Name of the restart file to parse

    Attributes
    ----------
    natom : int
        Number of atoms in the system
    resname : list of str
        List of all residue names in the system
    coords : list of float
        List of all coordinates in the format [x1, y1, z1, x2, y2, z2, ...]
    positions : natom x 3 distance Quantity
        2-D list of all coordinates with the appropriate distance unit attached.
        Has the format [ [x1, y1, z1], [x2, y2, z2], ... ]

    Example:
    >>> chm = CharmmCrdFile('testfiles/1tnm.crd')
    >>> print '%d atoms; %d coords' % (chm.natom, len(chm.coords))
    1414 atoms; 4242 coords
    """
    __metaclass__ = FileFormatType

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as a CHARMM coordinate file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is a CHARMM coordinate file
        """
        f = io.genopen(filename)
        line = f.readline().decode()

        try:
            while len(line.strip()) == 0:   # Skip whitespace, as a precaution
                line = f.readline().decode()

            intitle = True

            while intitle:
                line = f.readline().decode()
                if len(line.strip()) == 0:
                    intitle = False
                elif line.strip()[0] != '*':
                    intitle = False
                else:
                    intitle = True

            while len(line.strip()) == 0:      # Skip whitespace
                line = f.readline().decode()

            try:
                natom = int(line.strip().split()[0])

                for row in xrange(min(natom, 3)):
                    line = f.readline().decode().strip().split()
                    int(line[0])
                    int(line[1])
                    float(line[4])
                    float(line[5])
                    float(line[6])
                    int(line[8])
                    float(line[9])
            except (IndexError, ValueError):
                return False

            return True
        finally:
            f.close()


    def __init__(self, fname):
        self.atomno = []                   # Atom number
        self.resno = []                    # Residue number
        self.resname = []                  # Residue name
        self.resid = []                    # Residue ID
        self.attype = []                   # Atom type
        self.coords = []                   # 3N atomic coordinates
        self.title = []                    # .crd file title block
        self.segid = []                    # Segment ID
        self.weighting = []                # Atom weighting

        self.natom = 0                     # Number of atoms specified in file
        self._parse(fname)

    @property
    def positions(self):
        """
        Atomic coordinates with units attached to them with the shape (natom, 3)
        """
        return ([self.coords[i:i+3] for i in xrange(0, self.natom*3, 3)] *
                        u.angstroms)

    @property
    def coordinates(self):
        return self.coords

    def _parse(self, fname):

        crdfile = open(fname, 'r')
        try:
            line = crdfile.readline()

            while len(line.strip()) == 0:   # Skip whitespace, as a precaution
                line = crdfile.readline()

            intitle = True

            while intitle:
                self.title.append(line.strip())
                line = crdfile.readline()
                if len(line.strip()) == 0:
                    intitle = False
                elif line.strip()[0] != '*':
                    intitle = False
                else: 
                    intitle = True

            while len(line.strip()) == 0:      # Skip whitespace
                line = crdfile.readline()
            
            try:
                self.natom = int(line.strip().split()[0])
                
                for row in xrange(self.natom):
                    line = crdfile.readline().strip().split()
                    self.atomno.append(int(line[0]))
                    self.resno.append(int(line[1]))
                    self.resname.append(line[2])
                    self.attype.append(line[3])
                    self.coords.append(float(line[4]))
                    self.coords.append(float(line[5]))
                    self.coords.append(float(line[6]))
                    self.segid.append(line[7])
                    self.resid.append(int(line[8]))
                    self.weighting.append(float(line[9]))

                if 3*self.natom != len(self.coords):
                    raise CharmmFileError("Error parsing CHARMM .crd file: %d "
                                          "atoms requires %d coords (not %d)" %
                                          (self.natom, 3*self.natom,
                                           len(self.coords))
                    )

            except (ValueError, IndexError):
                raise CharmmFileError('Error parsing CHARMM coordinate file')
        finally:
            crdfile.close()

class CharmmRstFile(object):
    """
    Reads and parses data, velocities and coordinates from a CHARMM restart
    file (.rst) of file name 'fname' into class attributes

    Parameters
    ----------
    fname : str
        Name of the restart file to parse

    Attributes
    ----------
    natom : int
        Number of atoms in the system
    resname : list of str
        Names of all residues in the system
    coords : list of float
        List of all coordinates in the format [x1, y1, z1, x2, y2, z2, ...]
    coordsold : list of float
        List of all old coordinates in the format [x1, y1, z1, x2, y2, z2, ...]
    vels : list of float
        List of all velocities in the format [x1, y1, z1, x2, y2, z2, ...]
    positions : natom x 3 distance Quantity
        2-D list of all coordinates with the appropriate distance unit attached.
        Has the format [ [x1, y1, z1], [x2, y2, z2], ... ]
    positionsold : natom x 3 distance Quantity
        2-D list of all old coordinates with the appropriate distance unit
        attached.  Has the format [ [x1, y1, z1], [x2, y2, z2], ... ]
    velocities : natom x 3 distance/time Quantity
        2-D list of all old coordinates with the appropriate distance unit
        attached. Has the format [ [x1, y1, z1], [x2, y2, z2], ... ]

    Example:
    >>> chm = CharmmRstFile('testfiles/sample-charmm.rst')
    >>> print chm.header[0]
    REST    37     1
    >>> natom, nc, nco = chm.natom, len(chm.coords), len(chm.coordsold)
    >>> nv = len(chm.vels)
    >>> print '%d atoms; %d crds; %d old crds; %d vels' % (natom, nc, nco, nv)
    256 atoms; 768 crds; 768 old crds; 768 vels
    """
    __metaclass__ = FileFormatType

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as a CHARMM restart file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is a CHARMM restart file
        """
        f = io.genopen(filename)
        line = f.readline().decode()
        f.close()
        return line.startswith('REST')

    def __init__(self, fname):
        self.header = []
        self.title = []
        self.enrgstat = []
        self.coordsold = []
        self.coords = []
        self.vels = []
        
        self.ff_version = 0
        self.natom = 0
        self.npriv = 0
        self.nstep = 0
        self.nsavc = 0
        self.nsavv = 0
        self.jhstrt = 0

        self._parse(fname)

    @property
    def coordinates(self):
        return self.coords

    @property
    def positions(self):
        """ Atomic positions with units """
        return ([self.coords[i:i+3] for i in xrange(0, self.natom*3, 3)] *
                        u.angstroms)

    @property
    def positionsold(self):
        """ Old atomic positions with units """
        return ([self.coordsold[i:i+3] for i in xrange(0, self.natom*3, 3)] *
                        u.angstroms)

    @property
    def velocities(self):
        """ Atomic velocities with units """
        return ([self.vels[i:i+3] for i in xrange(0, self.natom*3, 3)] *
                        u.angstroms / u.picoseconds)

    def _parse(self, fname):

        crdfile = open(fname, 'r')
        try:
            readingHeader = True 

            while readingHeader:
                line = crdfile.readline()
                if not len(line):
                    raise CharmmFileError('Premature end of file')
                line = line.strip()
                words = line.split()
                if len(line) != 0:  
                    if words[0] == 'ENERGIES' or words[0] == '!ENERGIES':
                        readingHeader = False
                    else:
                        self.header.append(line.strip())
                else:
                    self.header.append(line.strip())

            for row in xrange(len(self.header)):
                if len(self.header[row].strip()) != 0:
                    line = self.header[row].strip().split()
                    if line[0][0:5] == 'NATOM' or line[0][0:6] == '!NATOM':
                        try:
                            line = self.header[row+1].strip().split()
                            self.natom = int(line[0])     
                            self.npriv = int(line[1])     # num. previous steps
                            self.nstep = int(line[2])     # num. steps in file  
                            self.nsavc = int(line[3])     # coord save frequency 
                            self.nsavv = int(line[4])     # velocities "
                            self.jhstrt = int(line[5])    # Num total steps?
                            break
                       
                        except (ValueError, IndexError):
                            raise CharmmFileError('Problem parsing CHARMM restart')

            self.scan(crdfile, '!XOLD')
            self._get_formatted_crds(crdfile, self.coordsold)

            self.scan(crdfile, '!VX')
            self._get_formatted_crds(crdfile, self.vels)

            self.scan(crdfile, '!X')
            self._get_formatted_crds(crdfile, self.coords)

            # Convert velocities to angstroms/ps
            self.vels = [v * ONE_TIMESCALE for v in self.vels]
        finally:
            crdfile.close()

    def scan(self, handle, str, r=0): # read lines in file till 'str' is found
        scanning = True

        if(r): handle.seek(0)

        while scanning:
            line = handle.readline()
            if not line:
                raise CharmmFileError('Premature end of file')

            if len(line.strip()) != 0:
                if line.strip().split()[0][0:len(str)] == str:
                    scanning = False


    def _get_formatted_crds(self, crdfile, crds):
        for row in xrange(self.natom):
            line = crdfile.readline()

            if not line:
                raise CharmmFileError('Premature end of file')

            if len(line) < 3*charmlen:
                raise CharmmFileError("Less than 3 coordinates present in "
                                      "coordinate row or coords may be "
                                      "truncated.") 

            line = line.replace('D','E')     # CHARMM uses 'D' for exponentials

            # CHARMM uses fixed format (len = charmlen = 22) for crds in .rst's

            crds.append(float(line[0:charmlen]))  
            crds.append(float(line[charmlen:2*charmlen]))
            crds.append(float(line[2*charmlen:3*charmlen]))


    def printcoords(self, crds):
        for crd in xrange(len(crds)):
            print crds[crd],
            if not (crd+1) % 3:
                print '\n',

if __name__ == '__main__':
    import doctest
    doctest.testmod()
