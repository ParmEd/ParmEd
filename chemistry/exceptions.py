""" This module contains all of the exceptions that are used in the chemistry
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
   def __init__(self, msg='Read error.'):
      self.msg = msg

class MoleculeError(ChemError):
   """ Error when molecule is illegally defined """
   def __init__(self, msg='Illegal definition of Molecule.'):
      self.msg = msg

class FileError(ChemError):
   """ Error when something illegal is done with files """
   def __init__(self, msg='Illegal file manipulation'):
      self.msg = msg

class FlagError(ChemError):
   """ Error when a FLAG is manipulated in readparm.amberParm """
   def __init__(self, msg='Bad flag'):
      self.msg = msg

class MaskError(ChemError):
   """ Error when a Mask is poorly formed """
   def __init__(self, msg='Bad mask'):
      self.msg = msg

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

class CmapError(ChemError):
   """ If there is an error with a correction-map potential """

class OpenMMError(ChemError):
   """ If there's a problem making an OpenMM system """

class CreateInputError(ChemError):
   """ If there's a problem making a mdin input file """
