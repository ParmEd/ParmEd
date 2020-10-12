"""
Contains classes for parsing DLPOLY topology and parameter files
"""
import os as _os
from ..utils import which as _which

__all__ = ['DLPOLY_TOPDIR', 'DlpolyFieldFile', 'DlpolyConfigFile']

DLPOLY_TOPDIR = None

try:
    del _testdir
except NameError:
    pass
del _os, _which

from ..dlpoly.dlpolyfield import DlpolyFieldFile
from ..dlpoly.dlpolyconfig import DlpolyConfigFile
