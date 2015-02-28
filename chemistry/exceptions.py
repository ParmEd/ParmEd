"""
This module contains all of the exceptions that are used in the chemistry
package
"""

class ChemError(Exception):
    def __init__(self, msg='Error'):
        self.msg = msg
    def __str__(self):
        return self.msg

class APIError(Exception):
    """ Raised if the API is used incorrectly """

class ChemWarning(Warning):
    """ Base warning class """

class ParsingError(ChemError):
    """ If there was a problem parsing any kind of file """

class PDBError(ParsingError):
    """ If there was a problem parsing a PDB file """

class Mol2Error(ParsingError):
    """ If there was a problem parsing a Mol2 file """

class PDBWarning(ChemWarning):
    """ A non-fatal error to indicate a problematic PDB file """

class MoleculeError(ChemError):
    """ Error when molecule is illegally defined """

class FileError(ChemError):
    """ Error when something illegal is done with files """

class FlagError(ChemError):
    """ Error when a FLAG is manipulated in readparm.amberParm """

class MaskError(ChemError):
    """ Error when a Mask is poorly formed """

class BondError(ChemError):
    """ This is what happens when you try to bond an atom to itself """

class AmberParmWarning(ChemWarning):
    """ If there is something that is non-fatal """

class AmberFormatWarning(ChemWarning):
    """ If there is something that is non-fatal """

class AmberParmError(ParsingError):
    """ This is a generic AmberParmError """

class AmberParameterError(ChemError):
    """ If there is an error with parameters """

class AmberParameterWarning(ChemWarning):
    """ If there is a warning to be raised with parameters """

class MoleculeWarning(ChemWarning):
    """ This occurs when there's a problem determining molecularity """

class DihedralError(ChemError):
    " This happens when we try to do disallowed things in the _Dihedral class "

class AmoebaParamFileError(ParsingError):
    """ When a parameter file is incorrect """

class AmoebaParamFileWarning(ChemWarning):
    """ When a parameter file is incorrect """

class TinkerFileError(ParsingError):
    """ Raised when one of the TINKER parsing routines hits a bad file """

class TinkerAnaloutError(ParsingError):
    """ When the analout file is not a valid format """

class TinkerKeyFileError(ParsingError):
    """ When a keyword control file is invalid """

class TinkerDynFileError(ParsingError):
    """ When a .dyn file is corrupt or badly formatted """

class CmapError(ChemError):
    """ If there is an error with a correction-map potential """

class OpenMMError(ChemError):
    """ If there's a problem making an OpenMM system """

class FormatError(ChemError):
    """ If there's a problem in formatting """

class AmoebaError(ChemError):
    """ If there's a problem with the AMOEBA force field """

class CreateInputError(ChemError):
    """ If there's a problem making a mdin input file """

class CharmmPSFError(ParsingError):
    """ If there is a problem parsing CHARMM PSF files """

class SplitResidueWarning(ChemWarning):
    """ For if a residue with the same number but different names is split """

class ResidueError(ChemError):
    """ For when there are problems defining a residue """

class CharmmPSFWarning(ChemWarning):
    """ For non-fatal PSF parsing issues """

class CharmmFileError(ParsingError):
    """ If there is a problem parsing CHARMM files """

class CharmmPsfEOF(ChemError):
    """ If we hit an end-of-file in parsing CHARMM files """

class MissingParameter(ChemError):
    """ If a parameter is missing from a database """

class AnisouWarning(ChemWarning):
    """ If there was a problem parsing an ANISOU record """

class MissingParameterWarning(ChemWarning):
    """ If a type of parameter is missing, but you don't want it to be fatal """

class AmberOFFWarning(ChemWarning):
    """ For badly formatted OFF files... ugh """

class PdbxError(ChemError):
    """ Class for catching general errors with PDBx/mmCIF parsing """

class PdbxSyntaxError(ChemError):
    """ Class for catching errors in mmCIF/PDBx syntax """
    def __init__(self, lineNumber, text):
        Exception.__init__(self)
        self.lineNumber = lineNumber
        self.text = text

    def __str__(self):
        return "%%ERROR - [at line: %d] %s" % (self.lineNumber, self.text)

class CpinResidueError(ChemError):
    """ Error adding a residue to the CPIN file """

class CpinChargeWarning(ChemWarning):
    """ Bad charge definitions that are inconsistent with protonation states """

class CpinRefEneWarning(ChemWarning):
    """ If not all reference energies are properly pKa-adjusted """

class CpinInputWarning(ChemWarning):
    """ If there is a non-fatal problem with the input variables """

class CpinInputError(ChemError):
    """ If the user provides bad input """

class FormatNotFound(ChemError):
    """ If the file format does not have a registered parser with it """
