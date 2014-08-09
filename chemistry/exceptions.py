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

class ReadError(ChemError):
    """ Error when files cannot be properly read """

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

class AmberParmError(ChemError):
    """ This is a generic AmberParmError """

class AmberParameterError(ChemError):
    """ If there is an error with parameters """

class AmberParameterWarning(ChemWarning):
    """ If there is a warning to be raised with parameters """

class MoleculeWarning(ChemWarning):
    """ This occurs when there's a problem determining molecularity """

class DihedralError(ChemError):
    " This happens when we try to do disallowed things in the _Dihedral class "

class AmoebaParamFileError(ChemError):
    """ When a parameter file is incorrect """

class AmoebaParamFileWarning(ChemWarning):
    """ When a parameter file is incorrect """

class TinkerFileError(ChemError):
    """ Raised when one of the TINKER parsing routines hits a bad file """

class TinkerAnaloutError(TinkerFileError):
    """ When the analout file is not a valid format """

class TinkerKeyFileError(TinkerFileError):
    """ When a keyword control file is invalid """

class TinkerDynFileError(TinkerFileError):
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

class CharmmPSFError(ChemError):
    """ If there is a problem parsing CHARMM PSF files """

class SplitResidueWarning(ChemWarning):
    """ For if a residue with the same number but different names is split """

class ResidueError(ChemError):
    """ For when there are problems defining a residue """

class CharmmPSFWarning(ChemWarning):
    """ For non-fatal PSF parsing issues """

class CharmmFileError(ChemError):
    """ If there is a problem parsing CHARMM files """

class CharmmPsfEOF(ChemError):
    """ If we hit an end-of-file in parsing CHARMM files """

class MissingParameter(ChemError):
    """ If a parameter is missing from a database """
