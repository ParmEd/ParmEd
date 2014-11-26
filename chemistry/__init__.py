"""
The chemistry package manipulates molecular structures and provides a way to go
between standard and amber file formats, manipulate structures, etc.
"""

__version__ = '2.0b'
__author__ = 'Jason Swails'

from chemistry import exceptions
from chemistry.structure import Structure, read_PDB
from chemistry.topologyobjects import *
