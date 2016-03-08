"""
This package contains classes responsible for loading rdkit objects
"""

from __future__ import print_function
from parmed.formats import PDBFile

try:
    from rdkit import Chem
except ImportError:
    pass

try:
    # py2
    import StringIO
except ImportError:
    import io
    StringIO = io.StringIO

class RDKit(object):

    @staticmethod
    def load(rmol):
        """
        Load a :class:`Mol` object and return a populated :class:`Structure`
        instance

        Parameters
        ----------
        rmol: :class:`Mol`
            RDKit :class:`Mol` object to convert

        Examples
        --------
        >>> from rdkit import Chem
        >>> import parmed as pmd
        >>> mol = Chem.MolFromSmiles('Cc1ccccc1')
        >>> struct = pmd.load_rdkit(mol)
        """
        fh = StringIO(Chem.MolToPDBBlock(rmol))
        return PDBFile.parse(fh)

    @staticmethod
    def from_smile(smile):
        """
        Load smile string to :class:`Structure`
        """
        return RDKit.load(Chem.MolFromSmiles(smile))

