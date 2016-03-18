"""
The parmed package manipulates molecular structures and provides a way to go
between standard and amber file formats, manipulate structures, etc.
"""

# Version format should be "major.minor.patch". For beta releases, attach
# "-beta#" to the end. The beta number will be turned into another number in the
# version tuple
__version__ = '2.4.0'
__author__ = 'Jason Swails'

__all__ = ['exceptions', 'periodic_table', 'residue', 'unit', 'utils',
           'Structure', 'StructureView', 'amber', 'charmm', 'namd', 'gromacs',
           'tinker', 'openmm', 'rosetta', 'rdkit', 'formats', 'Vec3', 'ParameterSet',
           'load_file', 'read_PDB', 'read_CIF', 'write_PDB', 'write_CIF',
           'load_rosetta', 'load_rdkit', 'download_PDB', 'download_CIF', 'tools', 'version']

from parmed import exceptions, periodic_table, residue
from parmed import unit, utils
from parmed.topologyobjects import *
from parmed.structure import Structure, StructureView
from parmed import amber, charmm, gromacs, namd, openmm, rosetta, tinker
from parmed import formats
from parmed.vec3 import Vec3
from parmed.parameters import ParameterSet
from parmed.rdkit import load_rdkit
from parmed import rdkit
from parmed.utils.decorators import deprecated as _deprecated
load_file = formats.load_file
read_PDB = formats.PDBFile.parse
read_CIF = formats.CIFFile.parse
write_PDB = _deprecated(formats.PDBFile.write)
write_CIF = _deprecated(formats.CIFFile.write)
load_rosetta = rosetta.RosettaPose.load

download_PDB = formats.PDBFile.download
download_CIF = formats.CIFFile.download

# The tools package depends on *everything*, so import this at the end to avoid
# circular imports.
from parmed import tools

# Add all of the objects from parmed.topologyobjects to the top-level namespace
__all__ += topologyobjects.__all__

# Build a version tuple from __version__ for easy comparison
import re
from collections import namedtuple
class version(namedtuple('version', ['major', 'minor', 'patchlevel'])):
    # Make sure betas always compare *less* than non-betas, and beta levels
    # compare greater than lower beta levels
    def __lt__(self, other):
        try:
            other = tuple(other)
        except TypeError:
            return NotImplemented
        other_is_beta = False
        if len(other) >= 4:
            other_is_beta = True
        if other_is_beta and self.beta is None:
            return tuple.__lt__(self, other[:3])
        elif other_is_beta and self.beta is not None:
            if tuple.__eq__(self, other[:3]):
                return self.beta < other[3]
            return tuple.__lt__(self, other)
        elif not other_is_beta and self.beta is not None:
            return tuple.__le__(self, other[:3])
        else:
            return tuple.__lt__(self, other)

    def __gt__(self, other):
        try:
            other = tuple(other)
        except TypeError:
            return NotImplemented
        other_is_beta = False
        if len(other) >= 4:
            other_is_beta = True
        if other_is_beta and self.beta is None:
            return tuple.__ge__(self, other[:3])
        elif other_is_beta and self.beta is not None:
            if tuple.__eq__(self, other[:3]):
                return self.beta > other[3]
            return tuple.__gt__(self, other)
        elif not other_is_beta and self.beta is not None:
            return tuple.__gt__(self, other)
        else:
            return tuple.__gt__(self, other)

    def __eq__(self, other):
        try:
            other = tuple(other)
        except TypeError:
            return NotImplemented
        other_is_beta = False
        if len(other) >= 4:
            other_is_beta = True
        if other_is_beta and self.beta is None:
            return False
        elif not other_is_beta and self.beta is not None:
            return False
        if other_is_beta:
            return tuple.__eq__(self, other[:3]) and self.beta == other[3]
        return tuple.__eq__(self, other)

    def __ne__(self, other):
        return not self == other

    def __ge__(self, other):
        return self > other or self == other
    def __le__(self, other):
        return self < other or self == other

    def __repr__(self):
        ret = super(type(self), self).__repr__()
        if self.beta is None:
            return ret
        return ret[:-1] + ', beta=%d)' % self.beta

# Build the version
_betare = re.compile(r'-beta*')
versionlist = list(int(x) for x in _betare.sub('.', __version__).split('.'))
version = version(*versionlist[:3])
if len(versionlist) > 3:
    version.beta = versionlist[3]
else:
    version.beta = None

# Clean up
del namedtuple, re, _betare, versionlist
