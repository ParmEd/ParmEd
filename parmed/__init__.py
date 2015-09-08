"""
The parmed package manipulates molecular structures and provides a way to go
between standard and amber file formats, manipulate structures, etc.
"""

__version__ = '2.0'
__author__ = 'Jason Swails'

__all__ = ['exceptions', 'periodic_table', 'residue', 'unit', 'utils',
           'Structure', 'StructureView', 'amber', 'charmm', 'gromacs', 'tinker',
           'openmm', 'rosetta', 'formats', 'Vec3', 'ParameterSet', 'load_file',
           'read_PDB', 'read_CIF', 'write_PDB', 'write_CIF', 'load_rosetta',
           'download_PDB', 'download_CIF', 'tools']

from parmed import exceptions, periodic_table, residue
from parmed import unit, utils
from parmed.topologyobjects import *
from parmed.structure import Structure, StructureView
from parmed import amber, charmm, gromacs, tinker, openmm, rosetta
from parmed import formats
from parmed.vec3 import Vec3
from parmed.parameters import ParameterSet
load_file = formats.load_file
read_PDB = formats.PDBFile.parse
read_CIF = formats.CIFFile.parse
write_PDB = formats.PDBFile.write
write_CIF = formats.CIFFile.write
load_rosetta = rosetta.RosettaPose.load

download_PDB = formats.PDBFile.download
download_CIF = formats.CIFFile.download

# The tools package depends on *everything*, so import this at the end to avoid
# circular imports.
from parmed import tools

# Add all of the objects from parmed.topologyobjects to the top-level namespace
__all__ += topologyobjects.__all__
