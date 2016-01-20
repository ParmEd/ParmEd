""" Package that provides an API to Amber-specific files """

__author__ = "Jason Swails <jason.swails@gmail.com>"

__all__ = ['AmberFormat', 'FortranFormat', 'AmberAsciiRestart', 'AmberMdcrd',
           'AmberMask', 'NetCDFTraj', 'NetCDFRestart', 'AmberOFFLibrary',
           'AmberParameterSet', 'AmberParm', 'ChamberParm', 'AmoebaParm',
           'Rst7', 'BeemanRestart', 'ConvertFromPSF', 'LoadParm', 'AMBERHOME',
           'titratable_residues']

from parmed.amber.amberformat import AmberFormat, FortranFormat
from parmed.amber.asciicrd import AmberAsciiRestart, AmberMdcrd
from parmed.amber.mask import AmberMask
from parmed.amber.netcdffiles import NetCDFTraj, NetCDFRestart
from parmed.amber.offlib import AmberOFFLibrary
from parmed.amber.parameters import AmberParameterSet
from parmed.amber.readparm import (AmberParm, ChamberParm, AmoebaParm,
                Rst7, BeemanRestart, ConvertFromPSF, LoadParm)
from parmed.amber import titratable_residues

# See if there is an AMBERHOME defined, which we will use by default. Otherwise,
# set it to the empty string
import os as _os
AMBERHOME = _os.getenv('AMBERHOME') or ''
del _os
