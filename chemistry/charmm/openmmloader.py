"""
This module used to contain OpenMM-enabled components of the CHARMM classes.
This functionality has since been incorporated directly into the original
classes. These classes are deprecated and may be phased out in later versions of
ParmEd

Author: Jason M. Swails
Contributors:
Date: April 29, 2015
"""
from chemistry.charmm.charmmcrds import (CharmmCrdFile as OpenMMCharmmCrdFile,
                                         CharmmRstFile as OpenMMCharmmRstFile)
from chemistry.charmm.psf import CharmmPsfFile as OpenMMCharmmPsfFile
import warnings

warnings.warn('chemistry.charmm.openmmloader has been deprecated. Use '
              'CharmmPsfFile, CharmmCrdFile, and CharmmRstFile from '
              'the chemistry.charmm package instead', DeprecationWarning)
