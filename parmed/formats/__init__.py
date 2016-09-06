"""
A package dealing with different file formats and automatic detection of those
formats
"""

__all__ = ['load_file', 'PDBFile', 'CIFFile', 'Mol2File', 'PSFFile', 'PQRFile', 'SDFFile']

from parmed.formats.registry import load_file
from parmed.formats.mol2 import Mol2File
from parmed.formats.pdb import PDBFile, CIFFile
from parmed.formats.pqr import PQRFile
from parmed.formats.psf import PSFFile
from parmed.formats.sdf import SDFFile

# Now let's modify structure.Structure and add our write methods from our
# various formats
from parmed.structure import Structure as _Structure
_Structure.write_pdb = PDBFile.write
_Structure.write_cif = CIFFile.write
_Structure.write_psf = PSFFile.write
