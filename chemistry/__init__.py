"""
The chemistry package manipulates molecular structures and provides a way to go
between standard and amber file formats, manipulate structures, etc.
"""

__version__ = '2.0b'
__author__ = 'Jason Swails'

from chemistry import exceptions, periodic_table
from chemistry.structure import Structure
from chemistry.topologyobjects import *
from chemistry import unit
from chemistry.residue import *
from chemistry import amber, charmm, tinker, openmm
from chemistry import formats
load_file = formats.load_file
read_PDB = formats.PDBFile.parse
read_CIF = formats.CIFFile.parse
write_PDB = formats.PDBFile.write
write_CIF = formats.CIFFile.write
