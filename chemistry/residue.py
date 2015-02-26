"""
This module contains basic information and functionality related to individual
residues in typical biopolymers.
"""

__all__ = ['AminoAcidResidue', 'RNAResidue', 'DNAResidue', 'ALA', 'ARG', 'ASN',
           'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
           'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'DA', 'DT', 'DG',
           'DC', 'A', 'U', 'G', 'C', 'WATER_NAMES', 'EXTRA_POINT_NAMES']

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BiomolecularResidue(object):
    """ Base class for different classes of biopolymer residues """
    _all_residues_by_name = dict()
    _all_residues_by_abbr = dict()
    _all_residues_by_symbol = dict()
    all_residues = []

    def __init_(self, *args, **kwargs):
        raise NotImplementedError('BiomolecularResidue must be subclassed')

    @classmethod
    def get(cls, key):
        raise NotImplementedError('BiomolecularResidue must be subclassed')

    def __str__(self):
        return self.name

    @classmethod
    def has(cls, thing):
        """
        Determines if a particular BiomolecularResidue or residue name is
        present in this classification of biomolecular residues

        Parameters
        ----------
        thing : str or :class:`BiomolecularResidue`

        Returns
        -------
        contains : bool
            If the residue or residue name *is* of this type, True. Otherwise,
            False.
        """
        if isinstance(thing, BiomolecularResidue):
            return thing in cls.all_residues
        try:
            cls.get(thing)
        except KeyError:
            return False
        else:
            return True

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AminoAcidResidue(BiomolecularResidue):
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
    aliases : list of str, optional
        A list of other abbreviations that *also* refer to this residue

    Raises
    ------
    ValueError
        If any aliases have the same abbreviation as *other* 
    """
    _all_residues_by_name = dict()
    _all_residues_by_abbr = dict()
    _all_residues_by_symbol = dict()
    all_residues = []

    def __init__(self, name, abbr, symbol, aliases=None):
        self.name = name
        self.abbr = abbr
        self.symbol = symbol
        type(self)._all_residues_by_name[name.upper()] = self
        type(self)._all_residues_by_abbr[abbr.upper()] = self
        type(self)._all_residues_by_symbol[symbol.upper()] = self
        type(self).all_residues.append(self)
        if aliases is not None:
            for alias in aliases:
                alias = alias.upper()
                if alias in type(self)._all_residues_by_abbr:
                    raise ValueError('%s is already an abbreviation' % alias)
                type(self)._all_residues_by_abbr[alias] = self

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
        residue : :class:`AminoAcidResidue`
            The residue corresponding to the given key

        Notes
        -----
        If the symbol is not defined, a KeyError is raised
        """
        if len(key) == 1:
            return cls._all_residues_by_symbol[key.upper()]
        if len(key) == 3:
            return cls._all_residues_by_abbr[key.upper()]
        # Handle C- and N-termini that may be prepended with C or N
        if len(key) == 4 and key[0].upper() in 'CN':
            return cls._all_residues_by_abbr[key[1:].upper()]
        return cls._all_residues_by_name[key.upper()]

ALA = AminoAcidResidue('Alanine', 'ALA', 'A')
ARG = AminoAcidResidue('Arginine', 'ARG', 'R')
ASN = AminoAcidResidue('Asparagine', 'ASN', 'N')
ASP = AminoAcidResidue('Aspartate' ,'ASP', 'D', ['ASH', 'AS4'])
CYS = AminoAcidResidue('Cysteine', 'CYS', 'C', ['CYM', 'CYX'])
GLU = AminoAcidResidue('Glutamate', 'GLU', 'E', ['GLH', 'GL4'])
GLN = AminoAcidResidue('Glutamine', 'GLN', 'Q')
GLY = AminoAcidResidue('Glycine', 'GLY', 'G')
HIS = AminoAcidResidue('Histidine', 'HIS', 'H', ['HIP', 'HIE', 'HID'])
ILE = AminoAcidResidue('Isoleucine', 'ILE', 'I')
LEU = AminoAcidResidue('Leucine', 'LEU', 'L')
LYS = AminoAcidResidue('Lysine', 'LYS', 'K', ['LYN'])
MET = AminoAcidResidue('Methionine', 'MET', 'M')
PHE = AminoAcidResidue('Phenylalanine', 'PHE', 'F')
PRO = AminoAcidResidue('Proline', 'PRO', 'P')
SER = AminoAcidResidue('Serine', 'SER', 'S')
THR = AminoAcidResidue('Threonine', 'THR', 'T')
TRP = AminoAcidResidue('Tryptophan', 'TRP', 'W')
TYR = AminoAcidResidue('Tyrosine', 'TYR', 'Y')
VAL = AminoAcidResidue('Valine', 'VAL', 'V')

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DNAResidue(BiomolecularResidue):
    """ An individual DNA residue

    Parameters
    ----------
    name : str
        The name of the residue
    abbr : str
        The abbreviation of the nucleic acid residue
    aliases : list of str, optional
        A list of other abbreviations that *also* refer to this residue
    """
    _all_residues_by_name = dict()
    _all_residues_by_abbr = dict()
    _all_residues_by_symbol = dict()
    all_residues = []

    def __init__(self, name, abbr, aliases=None):
        self.name = name
        self.abbr = abbr
        type(self)._all_residues_by_name[name.upper()] = self
        type(self)._all_residues_by_abbr[abbr.upper()] = self
        type(self).all_residues.append(self)
        if aliases is not None:
            for alias in aliases:
                alias = alias.upper()
                if alias in type(self)._all_residues_by_abbr:
                    raise ValueError('%s is already an abbreviation' % alias)
                type(self)._all_residues_by_abbr[alias] = self

    def __repr__(self):
        return '<DNA Residue %s: %s>' % (self.name, self.abbr)

    @classmethod
    def get(cls, key):
        """
        Gets the nucleic acid corresponding to either the residue name or
        abbreviation. It is case-insensitive. 

        Parameters
        ----------
        key : str
            abbreviation or residue name

        Returns
        -------
        residue : :class:`DNAResidue`
            The residue corresponding to the given key

        Notes
        -----
        If the symbol is not defined, a KeyError is raised
        """
        try:
            if key[-1] in '35':
                return cls._all_residues_by_abbr[key[:-1].upper()]
            return cls._all_residues_by_abbr[key.upper()]
        except KeyError:
            return cls._all_residues_by_name[key.upper()]

class RNAResidue(DNAResidue):
    """ An individual RNA residue

    Parameters
    ----------
    name : str
        The name of the residue
    abbr : str
        The abbreviation of the nucleic acid residue
    aliases : list of str, optional
        A list of other abbreviations that *also* refer to this residue
    """
    _all_residues_by_name = dict()
    _all_residues_by_abbr = dict()
    _all_residues_by_symbol = dict()
    all_residues = []

    def __repr__(self):
        return '<RNA Residue %s: %s>' % (self.name, self.abbr)

    @classmethod
    def get(cls, key):
        """
        Gets the nucleic acid corresponding to either the residue name or
        abbreviation. It is case-insensitive. 

        Parameters
        ----------
        key : str
            abbreviation or residue name

        Returns
        -------
        residue : :class:`RNAResidue`
            The residue corresponding to the given key

        Notes
        -----
        If the symbol is not defined, a KeyError is raised
        """
        try:
            if key[-1] in '35':
                return cls._all_residues_by_abbr[key[:-1].upper()]
            return cls._all_residues_by_abbr[key.upper()]
        except KeyError:
            return cls._all_residues_by_name[key.upper()]

DG = DNAResidue('Guanine', 'DG', ['GUA'])
DC = DNAResidue('Cytosine', 'DC', ['CYT'])
DA = DNAResidue('Adenine', 'DA', ['ADE'])
DT = DNAResidue('Thymine', 'DT', ['THY'])
G = RNAResidue('Guanine', 'G', ['GUA', 'RG'])
C = RNAResidue('Cytosine', 'C', ['CYT', 'RC'])
A = RNAResidue('Adenine', 'A', ['ADE', 'RA'])
U = RNAResidue('Uracil', 'U', ['URA', 'RU'])

WATER_NAMES = ['WAT', 'HOH', 'TIP3']
EXTRA_POINT_NAMES = ['EP', 'LP']
