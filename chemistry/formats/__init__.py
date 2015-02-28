"""
A package dealing with different file formats and automatic detection of those
formats
"""

__all__ = ['load_file', 'io', 'PDBFile', 'CIFFile']

from chemistry.formats.registry import load_file
from chemistry.formats import io
from chemistry.formats.mol2 import Mol2File
from chemistry.formats.pdb import PDBFile, CIFFile
from chemistry.formats.psf import PSFFile

# Now let's modify structure.Structure and add our write methods from our
# various formats
from chemistry.structure import Structure as _Structure
_Structure.write_pdb = PDBFile.write
_Structure.write_cif = CIFFile.write
_Structure.write_psf = PSFFile.write
