""" Package that provides an API to Amber-specific files """

from chemistry import __version__ as _chemistry_version
import warnings as _warnings

__version__ = _chemistry_version
__author__ = "Jason Swails <jason.swails@gmail.com>"

__all__ = ['leaprc', 'mask', 'mdcrd', 'netcdffiles', 'openmmloader',
           'openmmreporters', 'readparm', 'residue',
           # NetCDF objects
           'open_netcdf', 'get_int_dimension', 'get_float']

del _chemistry_version

# This determines which NetCDF package we're going to use...
NETCDF_PACKAGE = None

ALLOWED_NETCDF_PACKAGES = ('netCDF4', 'Scientific', 'pynetcdf', 'scipy')

NETCDF_INITIALIZED = False
SELECTED_NETCDF = ''

_FMT = 'NETCDF3_64BIT'

open_netcdf = get_int_dimension = get_float = None

try:
   from netCDF4 import Dataset as nc4NetCDFFile
   nc4open_netcdf = lambda name, mode: nc4NetCDFFile(name, mode, format=_FMT)
   nc4get_int_dimension = lambda obj, name: len(obj.dimensions[name])
   # Support 1-dimension arrays as scalars (since that's how Python-NetCDF
   # bindings write out scalars in Amber files)
   def nc4get_float(obj, name):
      try:
         return obj.variables[name].getValue()[0]
      except IndexError:
         return obj.variables[name][0]
      raise RuntimeError('Should not be here')
   _HAS_NC4 = True
except ImportError:
   nc4open_netcdf = nc4get_int_dimension = nc4get_float = None
   _HAS_NC4 = False

try:
   from Scientific.IO.NetCDF import NetCDFFile as sciNetCDFFile
   sciopen_netcdf = lambda name, mode: sciNetCDFFile(name, mode)
   sciget_int_dimension = lambda obj, name: obj.dimensions[name]
   sciget_float = lambda obj, name: obj.variables[name].getValue()
   _HAS_SCIENTIFIC_PYTHON = True
except ImportError:
   sciopen_netcdf = sciget_int_dimension = sciget_float = None
   _HAS_SCIENTIFIC_PYTHON = False

try:
   from pynetcdf.NetCDF import NetCDFFile as pynNetCDFFile
   pynopen_netcdf = lambda name, mode: pynNetCDFFile(name, mode)
   pynget_int_dimension = lambda obj, name: obj.dimensions[name]
   pynget_float = lambda obj, name: obj.variables[name].getValue()
   _HAS_PYNETCDF = True
except ImportError:
   pynopen_netcdf = pynget_int_dimension = pynget_float = None
   _HAS_PYNETCDF = False

try:
   from scipy.io.netcdf import netcdf_file as spNetCDFFile
   spopen_netcdf = lambda name, mode: spNetCDFFile(name, mode)
   spget_int_dimension = lambda obj, name: obj.dimensions[name]
   spget_float = lambda obj, name: obj.variables[name].getValue()
   _HAS_SCIPY_NETCDF = True
except ImportError:
   spopen_netcdf = spget_int_dimension = spget_float = None
   _HAS_SCIPY_NETCDF = False

HAS_NETCDF = (_HAS_NC4 or _HAS_SCIENTIFIC_PYTHON or 
              _HAS_PYNETCDF or _HAS_SCIPY_NETCDF)

def use(package=None):
   """
   Selects the NetCDF package to use

   Parameters:
    - package (string): This specifies which package to use, and may be either
         scipy, netCDF4, Scientific/ScientificPython, pynetcdf, or None.  If
         None, it chooses the first available implementation from the above list
         (in that order).
   
   Notes:
    - use() can only be called once, and once it is called there is no changing
      within that Python session or program. The 'netcdffiles' module calls this
      function to get the default NetCDF implementation if none has been
      selected before, so calling this function is unnecessary if the default
      implementation is sufficient. This is mostly useful for development
      testing as the backend NetCDF choice is virtually invisible to the user.

    - The pynetcdf implementation has long since been abandoned (ca.  2006), and
      is not recommended for use. It appears to parse NetCDF files just fine,
      but it does not seem to write them successfully according to my tests.
    
    - The NetCDF files have been tested against netCDF v. 1.0.4,
      Scientific v. 2.9.1, and scipy v. 0.13.1. Later versions are expected to
      work barring backwards-incompatible changes. Earlier versions are expected
      to work barring bugs or backwards-incompatible changes in the current or
      earlier versions.
   """
   global open_netcdf, get_int_dimension, get_float, NETCDF_INITIALIZED
   global HAS_NETCDF, SELECTED_NETCDF, nc4open_netcdf, nc4get_int_dimension
   global nc4get_float, sciopen_netcdf, sciget_int_dimension, sciget_float
   global pynopen_netcdf, pynget_int, pynget_float, ALLOWED_NETCDF_PACKAGES
   global _HAS_NC4, _HAS_SCIENTIFIC_PYTHON, _HAS_PYNETCDF, _HAS_SCIPY_NETCDF
   # IF we have already selected a package, warn and bail
   if NETCDF_INITIALIZED:
      if package is not None and SELECTED_NETCDF != package:
         if SELECTED_NETCDF in ('Scientific', 'ScientificPython') and \
               package in ('Scientific', 'ScientificPython'):
            pass
         else:
            _warnings.warn('Different NetCDF implementation already selected.'
                  '[%s]' % SELECTED_NETCDF)
      return

   if package is None:
      if _HAS_SCIPY_NETCDF:
         open_netcdf = spopen_netcdf
         get_int_dimension = spget_int_dimension
         get_float = spget_float
         SELECTED_NETCDF = 'scipy'
      elif _HAS_NC4:
         open_netcdf = nc4open_netcdf
         get_int_dimension = nc4get_int_dimension
         get_float = nc4get_float
         SELECTED_NETCDF = 'netCDF4'
      elif _HAS_SCIENTIFIC_PYTHON:
         open_netcdf = sciopen_netcdf
         get_int_dimension = sciget_int_dimension
         get_float = sciget_float
         SELECTED_NETCDF = 'ScientificPython'
      elif _HAS_PYNETCDF:
         open_netcdf = pynopen_netcdf
         get_int_dimension = pynget_int_dimension
         get_float = pynget_float
         SELECTED_NETCDF = 'pynetcdf'
   elif package == 'netCDF4':
      if not _HAS_NC4:
         raise ImportError('Could not find netCDF4 package')
      open_netcdf = nc4open_netcdf
      get_int_dimension = nc4get_int_dimension
      get_float = nc4get_float
      SELECTED_NETCDF = 'netCDF4'
   elif package == 'Scientific' or package == 'ScientificPython':
      if not _HAS_SCIENTIFIC_PYTHON:
         raise ImportError('Could not find package ScientificPython')
      open_netcdf = sciopen_netcdf
      get_int_dimension = sciget_int_dimension
      get_float = sciget_float
      SELECTED_NETCDF = 'ScientificPython'
   elif package == 'scipy':
      if not _HAS_SCIPY_NETCDF:
         raise ImportError('Could not find scipy NetCDF package')
      open_netcdf = spopen_netcdf
      get_int_dimension = spget_int_dimension
      get_float = spget_float
   elif package == 'pynetcdf':
      if not _HAS_PYNETCDF:
         raise ImportError('Could not find package pynetcdf')
      _warnings.warn('pynetcdf is no longer maintained and may not work '
            'properly. If you experience problems, try installing a '
            'different NetCDF package like Scientific, scipy, or netCDF4.')
      open_netcdf = pynopen_netcdf
      get_int_dimension = pynget_int_dimension
      get_float = pynget_float
      SELECTED_NETCDF = 'pynetcdf'
   else:
      raise ImportError('%s not a valid NetCDF package. Available options are '
               '%s' % (package, ', '.join(ALLOWED_NETCDF_PACKAGES)))

   NETCDF_INITIALIZED = True # We have now selected a NetCDF implementation
