"""
This module contains basic information and functionality related to individual
residues in typical biopolymers.
"""

class AminoAcidResidue(object):
    """
    An individual amino acid residue.

    Parameters
    ----------
    name : str
        The name of the residue
    abbr : str
        The 3-letter abbreviation of the amino acid residue
    symbol : str
        The 1-letter symbol of the amino acid
    """
    _all_residues_by_name = dict()
    _all_residues_by_abbr = dict()
    _all_residues_by_symbol = dict()

    def __init__(self, name, abbr, symbol):
        self.name = name
        self.abbr = abbr
        self.symbol = symbol
        AminoAcidResidue._all_residues_by_name[name.upper()] = self
        AminoAcidResidue._all_residues_by_abbr[abbr.upper()] = self
        AminoAcidResidue._all_residues_by_symbol[symbol.upper()] = self

    @classmethod
    def get(cls, key):
        """
        Gets the amino acid corresponding to either the residue name, 3-letter
        abbreviation or 1-letter symbol. It is case-insensitive.

        Parameters
        ----------
        key : str
            1-letter symbol, 3-letter abbreviation, or residue name

        Returns
        -------
        residue : AminoAcidResidue
            The residue corresponding to the given key

        Notes
        -----
        If the symbol is not defined, a KeyError is raised
        """
        if len(key) == 1:
            return cls._all_residues_by_symbol[key.upper()]
        if len(key) == 3:
            return cls._all_residues_by_abbr[key.upper()]
        return cls._all_residues_by_name[key.upper()]

ALA = AminoAcidResidue('Alanine', 'ALA', 'A')
GLY = AminoAcidResidue('Glycine', 'GLY', 'G')
VAL = AminoAcidResidue('Valine', 'VAL', 'V')
LEU = AminoAcidResidue('Leucine', 'LEU', 'L')
ILE = AminoAcidResidue('Isoleucine', 'ILE', 'I')
MET = AminoAcidResidue('Methionine', 'MET', 'M')
PHE = AminoAcidResidue('Phenylalanine', 'PHE', 'F')
TRP = AminoAcidResidue('Tryptophan', 'TRP', 'W')
PRO = AminoAcidResidue('Proline', 'PRO', 'P')
SER = AminoAcidResidue('Serine', 'SER', 'S')
THR = AminoAcidResidue('Threonine', 'THR', 'T')
CYS = AminoAcidResidue('Cysteine', 'CYS', 'C')
TYR = AminoAcidResidue('Tyrosine', 'TYR', 'Y')
ASN = AminoAcidResidue('Asparagine', 'ASN', 'N')
GLN = AminoAcidResidue('Glutamine', 'GLN', 'Q')
ASP = AminoAcidResidue('Aspartate' ,'ASP', 'D')
GLU = AminoAcidResidue('Glutamate', 'GLU', 'E')
LYS = AminoAcidResidue('Lysine', 'LYS', 'K')
ARG = AminoAcidResidue('Arginine', 'ARG', 'R')
HIS = AminoAcidResidue('Histidine', 'HIS', 'H')
