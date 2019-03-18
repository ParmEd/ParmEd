"""
The parmed package manipulates molecular structures and provides a way to go
between standard and amber file formats, manipulate structures, etc.
"""

# Version format should be "major.minor.patch". For beta releases, attach
# "-beta#" to the end. The beta number will be turned into another number in the
# version tuple
__author__ = 'Jason Swails'

__all__ = ['exceptions', 'periodic_table', 'residue', 'unit', 'utils',
           'Structure', 'StructureView', 'amber', 'charmm', 'namd', 'gromacs',
           'tinker', 'openmm', 'rosetta', 'rdkit', 'formats', 'Vec3', 'ParameterSet',
           'load_file', 'read_PDB', 'read_CIF', 'write_PDB', 'write_CIF',
           'load_rosetta', 'load_rdkit', 'download_PDB', 'download_CIF', 'tools', 'version']

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

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
read_pdb = read_PDB
read_cif = read_CIF
write_PDB = _deprecated(formats.PDBFile.write)
write_CIF = _deprecated(formats.CIFFile.write)
load_rosetta = rosetta.RosettaPose.load

download_PDB = formats.PDBFile.download
download_CIF = formats.CIFFile.download
download_pdb = download_PDB
download_cif = download_CIF

# The tools package depends on *everything*, so import this at the end to avoid
# circular imports.
from parmed import tools

# Add all of the objects from parmed.topologyobjects to the top-level namespace
__all__ += topologyobjects.__all__

# Build a version tuple from __version__ for easy comparison
import re
from collections import namedtuple
class version(namedtuple('version', ['major', 'minor',
                                     'patchlevel', 'commits_ahead'])):
    def __eq__(self, other):
        try:
            if len(other) == 3:
                return (other == self[:3] and self.commits_ahead == 0
                        and not self.dirty)
        except TypeError:
            return NotImplemented
        if tuple(other) != tuple(self):
            return False
        elif not hasattr(other, 'git_hash') or not hasattr(other, 'dirty'):
            return not self.dirty
        return self.git_hash == other.git_hash and self.dirty is other.dirty

    def __ne__(self, other):
        return not self == other

    def __le__(self, other):
        return self == other or self < other
    def __ge__(self, other):
        return self == other or self > other

# Build the version
_versionre = re.compile(r'(\d+)\.(\d+)\.(\d+)\+?(\d*)\.?g?([\dabcdefABCDEF]*)'
                        r'\.*(dirty)?')
if _versionre.match(__version__):
    versionlist = list(_versionre.match(__version__).groups())
    versionlist[3] = versionlist[3] or 0
    version = version(*[int(v) for v in versionlist[:4]])
    version.git_hash = versionlist[4]
    version.dirty = bool(versionlist[5])
else:
    versionlist = None
    version = version(0, 0, 0, 0)

# Clean up
del namedtuple, re, versionlist
