"""
The chemistry package manipulates molecular structures and provides a way to go
between standard and amber file formats, manipulate structures, etc.
"""

__version__ = '15.1'
__author__ = 'Jason Swails'

# Order here matters. Import order goes from fewest dependencies to most
from chemistry import exceptions, periodic_table
from chemistry import unit, utils
from chemistry.topologyobjects import *
from chemistry.residue import *
from chemistry.structure import Structure
from chemistry import amber, charmm, tinker, openmm, rosetta
from chemistry import formats
load_file = formats.load_file
read_PDB = formats.PDBFile.parse
read_CIF = formats.CIFFile.parse
write_PDB = formats.PDBFile.write
write_CIF = formats.CIFFile.write
load_rosetta = rosetta.RosettaPose.load
