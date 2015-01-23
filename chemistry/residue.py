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
    all_residues = []

    def __init__(self, name, abbr, symbol):
        self.name = name
        self.abbr = abbr
        self.symbol = symbol
        AminoAcidResidue._all_residues_by_name[name.upper()] = self
        AminoAcidResidue._all_residues_by_abbr[abbr.upper()] = self
        AminoAcidResidue._all_residues_by_symbol[symbol.upper()] = self
        AminoAcidResidue.all_residues.append(self)

    def __repr__(self):
        return '<Amino Acid Residue %s: %s [%s]>' % (self.name, self.abbr,
                self.symbol)

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
ARG = AminoAcidResidue('Arginine', 'ARG', 'R')
ASN = AminoAcidResidue('Asparagine', 'ASN', 'N')
ASP = AminoAcidResidue('Aspartate' ,'ASP', 'D')
CYS = AminoAcidResidue('Cysteine', 'CYS', 'C')
GLU = AminoAcidResidue('Glutamate', 'GLU', 'E')
GLN = AminoAcidResidue('Glutamine', 'GLN', 'Q')
GLY = AminoAcidResidue('Glycine', 'GLY', 'G')
HIS = AminoAcidResidue('Histidine', 'HIS', 'H')
ILE = AminoAcidResidue('Isoleucine', 'ILE', 'I')
LEU = AminoAcidResidue('Leucine', 'LEU', 'L')
LYS = AminoAcidResidue('Lysine', 'LYS', 'K')
MET = AminoAcidResidue('Methionine', 'MET', 'M')
PHE = AminoAcidResidue('Phenylalanine', 'PHE', 'F')
PRO = AminoAcidResidue('Proline', 'PRO', 'P')
SER = AminoAcidResidue('Serine', 'SER', 'S')
THR = AminoAcidResidue('Threonine', 'THR', 'T')
TRP = AminoAcidResidue('Tryptophan', 'TRP', 'W')
TYR = AminoAcidResidue('Tyrosine', 'TYR', 'Y')
VAL = AminoAcidResidue('Valine', 'VAL', 'V')

WATER_NAMES = ['WAT', 'HOH', 'TIP3']
EXTRA_POINT_NAMES = ['EP', 'LP']
