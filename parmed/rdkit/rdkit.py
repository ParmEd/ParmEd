"""
This package contains classes responsible for loading rdkit objects
"""

from __future__ import print_function, absolute_import
from parmed.formats import PDBFile
from parmed.utils.six.moves import StringIO


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
        from rdkit import Chem
        fh = StringIO(Chem.MolToPDBBlock(rmol))
        return PDBFile.parse(fh)

    @staticmethod
    def from_smiles(smiles):
        """
        Load smiles string to :class:`Structure`
        """
        from rdkit import Chem
        return RDKit.load(Chem.MolFromSmiles(smiles))

