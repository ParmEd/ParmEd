"""
The chemistry package manipulates molecular structures and provides a way to go
between standard and amber file formats, manipulate structures, etc.
"""

__version__ = '2.0b'
__author__ = 'Jason Swails'

from chemistry import exceptions, periodic_table
from chemistry.structure import Structure, read_PDB, write_PDB
from chemistry.topologyobjects import *
from chemistry import unit
from chemistry.residue import AminoAcidResidue
