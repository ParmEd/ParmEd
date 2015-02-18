""" Package that provides an API to Amber-specific files """

from chemistry import __version__ as _chemistry_version

__version__ = _chemistry_version
__author__ = "Jason Swails <jason.swails@gmail.com>"

del _chemistry_version

from chemistry.amber.amberformat import AmberFormat, FortranFormat
from chemistry.amber.asciicrd import AmberAsciiRestart, AmberMdcrd
from chemistry.amber.mask import AmberMask
from chemistry.amber.netcdffiles import use, NetCDFTraj, NetCDFRestart, HAS_NETCDF
from chemistry.amber.readparm import (AmberParm, ChamberParm, AmoebaParm,
                Rst7, BeemanRestart, ConvertFromPSF, LoadParm)
from chemistry.modeller import AmberOFFLibrary
