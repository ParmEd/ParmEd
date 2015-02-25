"""
A package dealing with different file formats and automatic detection of those
formats
"""

from chemistry.formats.registry import load_file
from chemistry.formats import io
from chemistry.formats.pdb import PDBFile, CIFFile
