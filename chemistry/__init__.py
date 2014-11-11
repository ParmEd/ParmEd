"""
The chemistry package manipulates molecular structures and provides a way to go
between standard and amber file formats, manipulate structures, etc.
"""

__all__ = ('amber',
           'exceptions',
           'formats',
           'molecule',
           'periodic_table',
           'system',
)

__version__ = '2.0b'
__author__ = 'Jason Swails'

from chemistry.topologyobjects import *
from chemistry.structure import Structure, read_PDB
