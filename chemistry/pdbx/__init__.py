"""
The PDBx package here was taken from the wwPDB website written by John Westbook
(and modified slightly to fix issues with Python 3 support as well as unified
exception handling)
"""

__author__ = "John Westbrook"
__contributors__ = "Jason Swails"
__version__ = "V0.01"
__license__ = "Creative Commons Attribution 3.0 Unported"
__email__ = "jwest@rcsb.rutgers.edu <or> jason.swails@gmail.com"

from chemistry.pdbx.PdbxReader import PdbxReader
from chemistry.pdbx.PdbxWriter import PdbxWriter
from chemistry.pdbx import PdbxContainers as containers
