"""
This module contains classes for reading and writing (single frame) NAMD binary
files.  The main functionality is currently aimed at manipulation rather than
analysis.
"""
from struct import unpack, pack

import numpy as np


class NamdBinFile(object):
    """From the NAMD manual:

        NAMD uses a trivial double-precision binary file format for 
        coordinates, velocities, and forces ... The file consists of the atom 
        count as a 32-bit integer followed by all three position or velocity 
        components for each atom as 64-bit double-precision floating point ...

    The main attributes are the number of atom entries (natom) and a (flat)
    numpy array of size 3*natom (values).  The meaning of "values" is
    effectively arbitrary, but for convenience derived classes are provided
    which alias the values to more descriptive names (e.g. "coordinates").

    See NamdBinCoor and NamdBinVel.
    """
    def __init__(self, values=None):
        self.values = np.asarray(values,np.float64)

    @property
    def natom(self):
        """The current number of atom entries."""
        return self.values.size / 3

    @classmethod
    def read(cls, fname):
        """Return an object from the values in a NAMD binary file."""
        infile = open(fname,'r')
        natoms = int(unpack('i',infile.read(4))[0])
        values = [unpack('d',infile.read(8))[0] for n in range(3*natoms)]
        infile.close()
        return cls(values)

    def write(self, fname):
        """Write the current attributes to a file."""
        outfile = open(fname,'w')
        outfile.write(pack('i',self.natom))
        for n in range(3*self.natom):
            outfile.write(pack('d',self.values[n]))
        outfile.close()

    def delatoms(self, indices):
        """Delete entries corresponding to the given atom indices."""
        del_indices = np.atleast_1d(np.asarray(indices,np.int32))
        mask = np.ones(self.natom,np.bool)
        mask[del_indices] = False
        newvalues = np.zeros(3*(self.natom - del_indices.size),np.float64)
        newvalues += self.values.reshape((self.natom,3))[mask].flatten()
        self.values = newvalues

    def insertatoms(self, start_index, natoms, values=None):
        """Insert space for natom entries beginning at start_index. If 
        specified, give them the provided values, otherwise set them to zero.
        """
        if values is not None:
            values = np.asarray(values)
            assert values.size == 3*natoms
        else:
            values = np.zeros(3*natoms)
        hinge = 3*start_index
        newvalues = np.concatenate(
            (self.values[:hinge],values,self.values[hinge:])
        )
        self.values = newvalues

    def copyatoms(self, start_index, natoms):
        """Convenience function, same as insertatoms() but set 'values' to
        be the same as the previous natoms' values (i.e. make a copy of them).
        """
        values = self.values[3*(start_index-natoms):3*start_index]       
        self.insertatoms(start_index,natoms,values)


class NamdBinCoor(NamdBinFile):
    """Class to read or write NAMD "bincoordinates" files."""
    @property
    def coordinates(self):
        return self.values


class NamdBinVel(NamdBinFile):
    """Class to read or write NAMD "binvelocities" files.

    NAMD internal units are assumed. These can be converted to 
    Angstrom / picosecond by multiplying by NamdBinVel.PDBVELFACTOR.
    """
    PDBVELFACTOR = 20.45482706

    @property
    def velocities(self):
        return self.values
