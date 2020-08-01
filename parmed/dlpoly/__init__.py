"""
Contains classes for parsing DLPOLY topology and parameter files
"""
import os as _os
from parmed.utils import which as _which

__all__ = ['DLPOLY_TOPDIR', 'DlpolyFieldFile', 'DlpolyConfigFile']

DLPOLY_TOPDIR = None

if _os.getenv('GMXDATA') is not None and _os.path.isdir(
        _os.path.join(_os.getenv('GMXDATA'), 'top')):
    DLPOLY_TOPDIR = _os.path.join(_os.getenv('GMXDATA'), 'top')
elif _os.getenv('GMXBIN') is not None and _os.path.isdir(
        _os.path.join(_os.getenv('GMXBIN'), '..', 'share', 'dlpoly', 'top')):
    DLPOLY_TOPDIR = _os.path.join(_os.getenv('GMXBIN'), '..', 'share',
                                  'dlpoly', 'top')
else:
    for _testdir in ['/usr', '/usr/local', '/opt/local', '/opt']:
        if _os.path.isdir(_os.path.join(_testdir, 'share', 'dlpoly')):
            DLPOLY_TOPDIR = _os.path.join(_testdir, 'share', 'dlpoly', 'top')
            break

if DLPOLY_TOPDIR is None:
    if _which('gmx') is not None:
        DLPOLY_TOPDIR = _os.path.join(_os.path.split(_which('gmx'))[0],
                                       '..', 'share', 'dlpoly', 'top')
    elif _which('pdb2gmx') is not None:
        DLPOLY_TOPDIR = _os.path.join(_os.path.split(_which('pdb2gmx'))[0],
                                       '..', 'share', 'dlpoly', 'top')

if DLPOLY_TOPDIR is not None:
    # Regularize the include path
    DLPOLY_TOPDIR = _os.path.realpath(DLPOLY_TOPDIR)
else:
    # Use the default Dlpoly installation path
    DLPOLY_TOPDIR = '/usr/local/dlpoly/share/dlpoly/top'

try:
    del _testdir
except NameError:
    pass
del _os, _which

from parmed.dlpoly.dlpolyfield import DlpolyFieldFile
from parmed.dlpoly.dlpolyconfig import DlpolyConfigFile
