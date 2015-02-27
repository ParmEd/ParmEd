"""
Contains classes for parsing GROMACS topology and parameter files
"""
import os as _os

from chemistry.gromacs.gromacstop import GromacsTopologyFile
from chemistry.utils import which as _which

GROMACS_TOPDIR = None

if os.getenv('GMXDATA') is not None and os.path.isdir(
        os.path.join(os.getenv('GMXDATA'), 'top')):
    GROMACS_TOPDIR = os.path.join(os.getenv('GMXDATA'), 'top')
elif os.getenv('GMXBIN') is not None and os.path.isdir(
        os.path.join(os.getenv('GMXBIN'), '..', 'share', 'gromacs', 'top')):
    GROMACS_TOPDIR = os.path.join(os.getenv('GMXBIN'), '..', 'share',
                                  'gromacs', 'top'))
else:
    for _testdir in ['/usr', '/usr/local', '/opt/local', '/opt']:
        if os.path.isdir(os.

del _testdir, _os, _which
