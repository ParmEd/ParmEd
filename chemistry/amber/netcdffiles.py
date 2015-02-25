"""
This module contains classes for reading and writing Amber NetCDF-style files,
including both restarts and trajectories. These sets of classes abstract the
interaction with the various NetCDF-Python APIs that are available, namely
    -  scipy
    -  netCDF4
    -  ScientificPython
    -  pynetcdf

This module contains objects relevant to Amber NetCDF files. The use() function
is responsible for selecting the API based on a default choice or user-selection
(the latter is really only helpful for development to ensure that all packages
work correctly---there is no difference from a user perspective). ALL
NetCDF-file manipulation that the chemistry/amber package does should be
contained in this module.
"""
from __future__ import division

from chemistry.formats.registry import FileFormatType
from chemistry import unit as u
import warnings
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
            val = obj.variables[name].getValue()
            if hasattr(val, '__iter__'):
                return val[0]
            return val
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

def _coerce_to_string(string, encoding='ascii'):
    """
    Decodes input to a string with the specified encoding if it is a bytes
    object. Otherwise, it just returns the input string.
    """
    try:
        return string.decode(encoding)
    except AttributeError:
        # Assume string
        return string

def use(package=None):
    """
    Selects the NetCDF package to use

    Parameters
    ----------
    package : str
        This specifies which package to use, and may be either scipy, netCDF4,
        Scientific/ScientificPython, pynetcdf, or None.  If None, it chooses the
        first available implementation from the above list (in that order).

    Notes
    -----
    The 'netcdffiles' module calls this function to get the default NetCDF
    implementation if none has been selected before, so calling this function is
    unnecessary if the default implementation is sufficient. This is mostly
    useful for development testing as the backend NetCDF choice is virtually
    invisible to the user.

    The pynetcdf implementation has long since been abandoned (ca. 2006), and is
    not recommended for use. It appears to parse NetCDF files just fine, but it
    does not seem to write them successfully according to my tests.
    
    The NetCDF files have been tested against netCDF v. 1.0.4, Scientific
    v. 2.9.1, and scipy v. 0.13.1. Later versions are expected to work barring
    backwards-incompatible changes. Other versions are expected to work barring
    bugs or backwards-incompatible changes in the current or earlier versions.
    """
    global open_netcdf, get_int_dimension, get_float, NETCDF_INITIALIZED
    global HAS_NETCDF, SELECTED_NETCDF, nc4open_netcdf, nc4get_int_dimension
    global nc4get_float, sciopen_netcdf, sciget_int_dimension, sciget_float
    global pynopen_netcdf, pynget_int, pynget_float, ALLOWED_NETCDF_PACKAGES
    global _HAS_NC4, _HAS_SCIENTIFIC_PYTHON, _HAS_PYNETCDF, _HAS_SCIPY_NETCDF

    if package is None:
        if NETCDF_INITIALIZED:
            return
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
        warnings.warn('pynetcdf is no longer maintained and may not work '
                      'properly. If you experience problems, try installing a '
                      'different NetCDF package like Scientific, scipy, or '
                      'netCDF4.')
        open_netcdf = pynopen_netcdf
        get_int_dimension = pynget_int_dimension
        get_float = pynget_float
        SELECTED_NETCDF = 'pynetcdf'
    else:
        raise ImportError('%s not a valid NetCDF package. Available options '
                          'are %s' % (package, 
                                      ', '.join(ALLOWED_NETCDF_PACKAGES))
        )
    
    NETCDF_INITIALIZED = True # We have now selected a NetCDF implementation

from chemistry import __version__
from compat24 import property, wraps
try:
    import numpy as np
except ImportError:
    # This is just to prevent NetCDF imports from bringing everything down if
    # numpy is not available. np can be used inside any method since it will be
    # required for any NetCDF to work in the first place.
    np = None

def needs_netcdf(fcn):
    """
    Decorator to protect against functions that need NetCDF so we can provide
    a helpful error message
    """
    @wraps(fcn)
    def new_fcn(*args, **kwargs):
        if not HAS_NETCDF:
            raise ImportError('No NetCDF packages are available!')
        if not NETCDF_INITIALIZED:
            use() # Set up a default NetCDF implementation
        return fcn(*args, **kwargs)
    return new_fcn

class NetCDFRestart(object):
    """ Class to read or write NetCDF restart files """
    __metaclass__ = FileFormatType

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as an Amber NetCDF restart file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is an Amber NetCDF restart file. False otherwise
        """
        if not HAS_NETCDF:
            return False # Can't determine...
        if not NETCDF_INITIALIZED:
            use()
        try:
            f = open_netcdf(filename, 'r')
        except: # Bare except... each package raises different exceptions
            return False
        try:
            try:
                if _coerce_to_string(f.Conventions) != 'AMBERRESTART':
                    return False
            except AttributeError:
                return False
            # Passed all our tests
            return True
        finally:
            f.close()

    @needs_netcdf
    def __init__(self, fname, mode='r'):
        """
        Opens a NetCDF File. The main constructor should never be called
        directly.  The alternative "open_old" and "open_new" constructors should
        be called instead for reading existing or writing new NetCDF files,
        respectively.
        """
        self.closed = False
        self._ncfile = open_netcdf(fname, mode)
   
    @classmethod
    @needs_netcdf
    def open_new(cls, fname, natom, box, vels, title='',
                 remd=None, temp=None, remd_indices=None,
                 remd_groups=None, remd_dimtypes=None):
        """
        Opens a new NetCDF file and sets the attributes

        Parameters
        ----------
        fname : str
            Name of the new file to open (overwritten)
        natom : int
            The number of atoms in the system
        box : bool
            Whether unit cell information is written or not
        vels : bool
            Whether velocity information is written or not
        title : str=''
            The title to write to the NetCDF restart file
        remd : str=None
            None -- No REMD information is written
            'T[emperature]' -- target temperature (or pH) will be written
            'M[ulti-D]' -- remd_indices and remd_groups will be written
        remd_dimtypes : iterable of int=None
            Array of exchange types for each group. The length will be the REMD
            dimension (if `remd` above is "M[ulti-D]")

        Notes
        -----
        `remd` is case-insensitive, and done based on first-letter matching
        """
        if remd is not None:
            if remd[0] in 'tT':
                remd_type = 'TEMPERATURE'
                if temp is None:
                    raise ValueError('temp must be specified for '
                                     'T-REMD restarts.')
            elif remd[0] in 'mM':
                remd_type = 'MULTI'
                if (remd_indices is None or remd_groups is None or 
                    len(remd_indices) != len(remd_groups)):
                    raise ValueError('remd_indices and remd_groups must be '
                                     'given for multi-D REMD, and must have '
                                     'the same length.')
                remd_dimension = len(remd_indices)
            else:
                raise ValueError('remd must be None, T[emperature] or M[ulti]')
        else:
            remd_type = None
        inst = cls(fname, 'w')
        ncfile = inst._ncfile
        inst.hasbox = bool(box)
        inst.hasvels = bool(vels)
        # Assign the main attributes
        ncfile.Conventions = 'AMBERRESTART'
        ncfile.ConventionVersion = "1.0"
        ncfile.title = title
        ncfile.application = "AMBER"
        ncfile.program = "ParmEd"
        ncfile.programVersion = str(__version__)
        # Make all of the dimensions
        ncfile.createDimension('spatial', 3)
        ncfile.createDimension('atom', natom)
        inst.spatial = 3
        inst.atom = natom
        inst.title = ncfile.title
        if box:
            ncfile.createDimension('cell_spatial', 3)
            ncfile.createDimension('label', 5)
            ncfile.createDimension('cell_angular', 3)
            inst.cell_spatial = 3
            inst.label = 5
            inst.cell_angular = 3
        if remd_type == 'MULTI':
            ncfile.createDimension('remd_dimension', remd_dimension)
            inst.remd_dimension = remd_dimension
        ncfile.createDimension('time', 1)
        # Now make the variables
        v = ncfile.createVariable('time', 'd', ('time',))
        v.units = 'picosecond'
        v = ncfile.createVariable('spatial', 'c', ('spatial',))
        v[:] = np.asarray(list('xyz'))
        if inst.hasbox:
            v = ncfile.createVariable('cell_angular', 'c',
                                            ('cell_angular', 'label'))
            v[0] = np.asarray(list('alpha'))
            v[1] = np.asarray(list('beta '))
            v[2] = np.asarray(list('gamma'))
            v = ncfile.createVariable('cell_spatial', 'c',
                                            ('cell_spatial',))
            v[0], v[1], v[2] = 'a', 'b', 'c'
            v = ncfile.createVariable('cell_lengths', 'd',
                                            ('cell_spatial',))
            v.units = 'angstrom'
            v = ncfile.createVariable('cell_angles', 'd',
                                            ('cell_angular',))
            v.units = 'degree'
        v = ncfile.createVariable('coordinates', 'd',
                                        ('atom', 'spatial'))
        v.units = 'angstrom'
        if inst.hasvels:
            v = ncfile.createVariable('velocities', 'd',
                                            ('atom', 'spatial'))
            try:
                # Prevent NetCDF4 from trying to autoscale the values. Ugh.
                v.set_auto_maskandscale(False)
            except AttributeError:
                # Other packages do not have this variable.
                pass
            v.units = 'angstrom/picosecond'
            v.scale_factor = np.float32(20.455)
            inst.velocity_scale = 20.455

        if remd_type == 'TEMPERATURE':
            v = ncfile.createVariable('temp0', 'd', ('time',))
            v.units = 'kelvin'
        elif remd_type == 'MULTI':
            ncfile.createVariable('remd_indices', 'i',
                                        ('remd_dimension',))
            ncfile.createVariable('remd_groups', 'i',
                                        ('remd_dimension',))
            ncfile.createVariable('remd_dimtype', 'i',
                                        ('remd_dimension',))

        return inst

    @classmethod
    @needs_netcdf
    def open_old(cls, fname):
        """
        Opens the NetCDF file and sets the global attributes that the file sets

        Parameters
        ----------
        fname : str
            Name of the file to read
        """
        inst = cls(fname, 'r')
        ncfile = inst._ncfile
        inst.Conventions = _coerce_to_string(ncfile.Conventions)
        inst.ConventionVersion = _coerce_to_string(ncfile.ConventionVersion)
        inst.application = _coerce_to_string(ncfile.application)
        inst.program = _coerce_to_string(ncfile.program)
        inst.programVersion = _coerce_to_string(ncfile.programVersion)
        inst.title = _coerce_to_string(ncfile.title)
        # Set up the dimensions as attributes
        for dim in ncfile.dimensions:
            # Exception for ParmEd-created ncrst files
            if dim == 'time': continue
            setattr(inst, dim, get_int_dimension(ncfile, dim))
        inst.hasvels = 'velocities' in ncfile.variables
        inst.hasbox = ('cell_lengths' in ncfile.variables and
                       'cell_angles' in ncfile.variables)
        if inst.hasvels:
            vels = ncfile.variables['velocities']
            inst.velocity_scale = vels.scale_factor
            # Turn off automatic scaling for netCDF4. Ugh.
            try:
                vels.set_auto_maskandscale(False)
            except AttributeError:
                # Other packages do not have this variable
                pass
        return inst

    parse = open_old

    @property
    def coordinates(self):
        return self._ncfile.variables['coordinates'][:].flatten()

    @coordinates.setter
    def coordinates(self, stuff):
        self._ncfile.variables['coordinates'][:] = \
                                np.reshape(stuff, (self.atom, 3))
        self.flush()

    @property
    def velocities(self):
        return (self._ncfile.variables['velocities'][:].flatten() *
                self.velocity_scale)

    @velocities.setter
    def velocities(self, stuff):
        self._ncfile.variables['velocities'][:] = \
                    np.reshape(stuff, (self.atom, 3)) / self.velocity_scale
        self.flush()

    @property
    def cell_lengths(self):
        return self._ncfile.variables['cell_lengths'][:]
   
    @cell_lengths.setter
    def cell_lengths(self, stuff):
        self._ncfile.variables['cell_lengths'][:] = np.asarray(stuff)
        self.flush()

    @property
    def cell_angles(self):
        return self._ncfile.variables['cell_angles'][:]
   
    @cell_angles.setter
    def cell_angles(self, stuff):
        self._ncfile.variables['cell_angles'][:] = np.asarray(stuff)
        self.flush()

    @property
    def box(self):
        leng, ang = self.cell_lengths, self.cell_angles
        return np.concatenate((leng, ang))

    @box.setter
    def box(self, stuff):
        self.cell_lengths = stuff[:3]
        self.cell_angles = stuff[3:]

    @property
    def time(self):
        return get_float(self._ncfile, 'time')

    @time.setter
    def time(self, stuff):
        self._ncfile.variables['time'][0] = float(stuff)
        self.flush()

    @property
    def temp0(self):
        return get_float(self._ncfile, 'temp0')

    @temp0.setter
    def temp0(self, stuff):
        self._ncfile.variables['temp0'][0] = float(stuff)
        self.flush()
   
    @property
    def remd_indices(self):
        return self._ncfile.variables['remd_indices'][:]

    @remd_indices.setter
    def remd_indices(self, stuff):
        self._ncfile.variables['remd_indices'][:] = np.asarray(stuff, dtype='i')
        self.flush()

    @property
    def remd_groups(self):
        return self._ncfile.variables['remd_groups'][:]

    @remd_groups.setter
    def remd_groups(self, stuff):
        self._ncfile.variables['remd_groups'][:] = np.asarray(stuff, dtype='i')
        self.flush()

    @property
    def remd_dimtype(self):
        return self._ncfile.variables['remd_dimtype'][:]

    @remd_dimtype.setter
    def remd_dimtype(self, stuff):
        self._ncfile.variables['remd_dimtype'][:] = np.asarray(stuff, dtype='i')
        self.flush()

    def close(self):
        self.closed = True
        self._ncfile.close()

    def __del__(self):
        self.closed or (hasattr(self, '_ncfile') and self._ncfile.close())

    def flush(self):
        try:
            self._ncfile.flush()
        except AttributeError:
            pass

class NetCDFTraj(object):
    """ Class to read or write NetCDF restart files

    Parameters
    ----------
    fname : str
        Name of the file to open
    mode : str
        Mode to open in:
            - 'w' means write-mode
            - 'r' means read-mode

    Notes
    -----
    You should use the open_new and open_old alternative constructors instead of
    the default constructor
    """
    __metaclass__ = FileFormatType

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as an Amber NetCDF trajectory file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is an Amber NetCDF trajectory file. False otherwise
        """
        if not HAS_NETCDF:
            return False # Can't determine...
        if not NETCDF_INITIALIZED:
            use()
        try:
            f = open_netcdf(filename, 'r')
        except: # Bare except... each package raises different exceptions
            return False
        try:
            try:
                if _coerce_to_string(f.Conventions) != 'AMBER':
                    return False
            except AttributeError:
                return False
            # Passed all our tests
            return True
        finally:
            f.close()

    @needs_netcdf
    def __init__(self, fname, mode='r'):
        """ Opens a NetCDF File """
        self.closed = False
        self._ncfile = open_netcdf(fname, mode)
   
    @classmethod
    @needs_netcdf
    def open_new(cls, fname, natom, box, crds=True, vels=False, frcs=False,
                 remd=None, remd_dimension=None, title=''):
        """
        Opens a new NetCDF file and sets the attributes

        Parameters
        ----------
        fname : str
            Name of the new file to open (overwritten)
        natom : int
            Number of atoms in the restart
        box : bool
            Indicates if cell lengths and angles are written to the NetCDF file
        crds : bool=True
            Indicates if coordinates are written to the NetCDF file
        vels : bool=False
            Indicates if velocities are written to the NetCDF file
        frcs : bool=False
            Indicates if forces are written to the NetCDF file
        remd : str=None
            'T[emperature]' if replica temperature is written
            'M[ulti]' if Multi-D REMD information is written
            None if no REMD information is written
        remd_dimension : int=None
            If remd above is 'M[ulti]', this is how many REMD dimensions exist
        title : str=''
            The title of the NetCDF trajectory file
        """
        inst = cls(fname, 'w')
        ncfile = inst._ncfile
        if remd is not None:
            if remd[0] in 'Tt':
                inst.remd = 'TEMPERATURE'
            elif remd[0] in 'Mm':
                inst.remd = 'MULTI'
                if remd_dimension is None:
                    raise ValueError('remd_dimension must be given '
                                     'for multi-D REMD')
                inst.remd_dimension = int(remd_dimension)
            else:
                raise ValueError('remd must be T[emperature] or M[ultiD]')
        else:
            inst.remd = None
        inst.hasbox = bool(box)
        inst.hasvels = bool(vels)
        inst.hascrds = bool(crds)
        inst.hasfrcs = bool(frcs)
        # Assign the main attributes
        ncfile.Conventions = "AMBER"
        ncfile.ConventionVersion = "1.0"
        ncfile.application = "AmberTools"
        ncfile.program = "ParmEd"
        ncfile.programVersion = __version__
        ncfile.title = "ParmEd-created trajectory"
        inst.Conventions = "AMBER"
        inst.ConventionVersion = "1.0"
        inst.application = "AmberTools"
        inst.program = "ParmEd"
        inst.programVersion = __version__
        inst.title = ncfile.title
        # Create the dimensions
        ncfile.createDimension('frame', None)
        ncfile.createDimension('spatial', 3)
        ncfile.createDimension('atom', natom)
        if inst.remd == 'MULTI':
            ncfile.createDimension('remd_dimension', inst.remd_dimension)
        inst.frame, inst.spatial, inst.atom = None, 3, natom
        if inst.hasbox:
            ncfile.createDimension('cell_spatial', 3)
            ncfile.createDimension('cell_angular', 3)
            ncfile.createDimension('label', 5)
            inst.cell_spatial, inst.cell_angular, inst.label = 3, 3, 5
        # Create the variables and assign units and scaling factors
        v = ncfile.createVariable('spatial', 'c', ('spatial',))
        v[:] = np.asarray(list('xyz'))
        if inst.hasbox:
            v = ncfile.createVariable('cell_spatial', 'c',
                                            ('cell_spatial',))
            v[:] = np.asarray(list('abc'))
            v = ncfile.createVariable('cell_angular', 'c',
                                            ('cell_angular', 'label',))
            v[:] = np.asarray([list('alpha'), list('beta '), list('gamma')])
        v = ncfile.createVariable('time', 'f', ('frame',))
        v.units = 'picosecond'
        if inst.hascrds:
            v = ncfile.createVariable('coordinates', 'f',
                                            ('frame', 'atom', 'spatial'))
            v.units = 'angstrom'
            inst._last_crd_frame = 0
        if inst.hasvels:
            v = ncfile.createVariable('velocities', 'f',
                                            ('frame', 'atom', 'spatial'))
            v.units = 'angstrom/picosecond'
            inst.velocity_scale = v.scale_factor = 20.455
            inst._last_vel_frame = 0
            try:
                # Prevent NetCDF4 from trying to autoscale the values. Ugh.
                v.set_auto_maskandscale(False)
            except AttributeError:
                # Other packages do not have this variable.
                pass
        if inst.hasfrcs:
            v = ncfile.createVariable('forces', 'f',
                                            ('frame', 'atom', 'spatial'))
            v.units = 'kilocalorie/mole/angstrom'
            inst._last_frc_frame = 0
        if inst.hasbox:
            v = ncfile.createVariable('cell_lengths', 'd',
                                            ('frame', 'cell_spatial'))
            v.units = 'angstrom'
            v = ncfile.createVariable('cell_angles', 'd',
                                            ('frame', 'cell_angular'))
            v.units = 'degree'
            inst._last_box_frame = 0
        if inst.remd == 'TEMPERATURE':
            v = ncfile.createVariable('temp0', 'd', ('frame',))
            v.units = 'kelvin'
            inst._last_remd_frame = 0
        elif inst.remd == 'MULTI':
            ncfile.createVariable('remd_indices', 'i',
                                        ('frame', 'remd_dimension'))
            ncfile.createVariable('remd_dimtype', 'i',
                                        ('remd_dimension',))
            inst._last_remd_frame = 0
    
        inst._last_time_frame = 0

        return inst

    @classmethod
    @needs_netcdf
    def open_old(cls, fname):
        """
        Opens the NetCDF file and sets the global attributes that the file sets

        Parameters
        ----------
        fname : str
            File name of the trajectory to open. It must exist
        """
        inst = cls(fname, 'r')
        ncfile = inst._ncfile
        inst.Conventions = _coerce_to_string(ncfile.Conventions)
        inst.ConventionVersion = _coerce_to_string(ncfile.ConventionVersion)
        inst.application = _coerce_to_string(ncfile.application)
        inst.program = _coerce_to_string(ncfile.program)
        inst.programVersion = _coerce_to_string(ncfile.programVersion)
        inst.title = _coerce_to_string(ncfile.title)
        # Set up the dimensions as attributes
        for dim in ncfile.dimensions:
            setattr(inst, dim, get_int_dimension(ncfile, dim))
        inst.hascrds = 'coordinates' in ncfile.variables
        inst.hasvels = 'velocities' in ncfile.variables
        inst.hasfrcs = 'forces' in ncfile.variables
        inst.hasbox = ('cell_lengths' in ncfile.variables and
                       'cell_angles' in ncfile.variables)
        if inst.hasvels:
            inst.velocity_scale = ncfile.variables['velocities'].scale_factor
            try:
                # Prevent NetCDF4 from trying to autoscale the values. Ugh.
                ncfile.variables['velocities'].set_auto_maskandscale(False)
            except AttributeError:
                # Other packages do not have this variable.
                pass
        if inst.frame is None:
            if 'time' in ncfile.variables:
                inst.frame = len(ncfile.variables['time'][:])
            elif inst.hascrds:
                inst.frame = len(ncfile.variables['coordinates'][:])
            elif inst.hasvels:
                inst.frame = len(ncfile.variables['velocities'][:])
            elif inst.hasfrcs:
                inst.frame = len(ncfile.variables['forces'][:])
        return inst

    def coordinates(self, frame):
        """
        Get the coordinates of a particular frame in the trajectory
    
        Parameters
        ----------
        frame : int
            Which snapshot to get (first snapshot is frame 0)
    
        Returns
        -------
        coordinates
            numpy array of length 3*natom with the given coordinates in the
            format [x1, y1, z1, x2, y2, z2, ...] in Angstroms
        """
        return self._ncfile.variables['coordinates'][frame][:].flatten()

    def add_coordinates(self, stuff):
        """
        Adds a new coordinate frame to the end of a NetCDF trajectory. This
        should only be called on objects created with the "open_new"
        constructor.

        Parameters
        ----------
        stuff : iterable of floats or distance Quantity
            This array of floats is converted into a numpy array of shape
            (natom, 3). It can be passed either in the 2-D format of
            [ [x1, y1, z1], [x2, y2, z2], ... ] or in the 1-D format of
            [x1, y1, z1, x2, y2, z2, ... ].
        """
        if u.is_quantity(stuff):
            stuff = stuff.value_in_unit(u.angstroms)
        stuff = np.asarray(stuff)
        self._ncfile.variables['coordinates'][self._last_crd_frame] = \
                np.reshape(stuff, (self.atom, 3))
        self._last_crd_frame += 1
        self.flush()

    def velocities(self, frame):
        """
        Get the velocities of a particular frame in the trajectory

        Parameters
        ----------
        frame : int
            Which snapshot to get (first snapshot is frame 0)

        Returns
        -------
        velocities
            numpy array of length 3*atom with the given velocities properly
            scaled to be of units angstrom/picosecond
        """
        return (self._ncfile.variables['velocities'][frame][:].flatten() * 
                self.velocity_scale)

    def add_velocities(self, stuff):
        """
        Adds a new velocities frame to the end of a NetCDF trajectory. This
        should only be called on objects created with the "open_new"
        constructor.

        Parameters
        ----------
        stuff : iterable of floats or distance/time Quantity
            This array of floats is converted into a numpy array of shape
            (natom, 3). It can be passed either in the 2-D format of
            [ [x1, y1, z1], [x2, y2, z2], ... ] or in the 1-D format of
            [x1, y1, z1, x2, y2, z2, ... ].
        """
        if u.is_quantity(stuff):
            stuff = stuff.value_in_unit(u.angstrom/u.picosecond)
        stuff = np.asarray(stuff)
        self._ncfile.variables['velocities'][self._last_vel_frame] = \
                np.reshape(stuff, (self.atom, 3)) / self.velocity_scale
        self._last_vel_frame += 1
        self.flush()

    def forces(self, frame):
        """
        Get the forces of a particular frame in the trajectory

        Parameters
        ----------
        frame : int
            Which snapshot to get (first snapshot is frame 0)

        Returns
        -------
        forces
            numpy array of length 3*atom with the given forces properly
            scaled to be of kcal/mol/Angstroms
        """
        return (self._ncfile.variables['forces'][frame][:].flatten())

    def add_forces(self, stuff):
        """
        Adds a new coordinate frame to the end of a NetCDF trajectory. This
        should only be called on objects created with the "open_new"
        constructor.

        Parameters
        ----------
        stuff : iterable of floats or energy/distance Quantity
            This array of floats is converted into a numpy array of shape
            (natom, 3). It can be passed either in the 2-D format of
            [ [x1, y1, z1], [x2, y2, z2], ... ] or in the 1-D format of
            [x1, y1, z1, x2, y2, z2, ... ].
        """
        if u.is_quantity(stuff):
            stuff.value_in_unit(u.kilocalories_per_mole/u.angstroms)
        self._ncfile.variables['forces'][self._last_frc_frame] = \
                np.reshape(stuff, (self.atom, 3))
        self._last_frc_frame += 1
        self.flush()

    def cell_lengths_angles(self, frame):
        """
        Get the cell lengths and cell angles of a particular frame in the
        trajectory

        Parameters
        ----------
        frame : int
            Which snapshot to get (first snapshot is frame 0)

        Returns
        -------
        lengths, angles
            2-element tuple: (length-3 numpy array of cell lengths, length-3
            numpy array of cell angles) in Angstroms
        """
        return (self._ncfile.variables['cell_lengths'][frame][:],
                self._ncfile.variables['cell_angles'][frame][:])

    def box(self, frame):
        """
        Get the cell lengths and angles as one array

        Parameters
        ----------
        frame : int
            Frame to get the box from

        Returns
        -------
            6-element array with 3 lengths followed by 3 angles (in degrees)
        """
        return np.concatenate(self.cell_lengths_angles(frame))
   
    def add_cell_lengths_angles(self, lengths, angles=None):
        """
        Adds a new cell length and angle frame to the end of a NetCDF
        trajectory.  This should only be called on objects created with the
        "open_new" constructor.

        Parameters
        ----------
        lengths : array of 3 (or 6) floats (or Quantities)
            This should be a 1-D array of 3 or 6 elements.  If 6 elements,
            `angles` should be None (below) and the first 3 elements are the box
            lengths (angstroms) and the last 3 are the box angles (degrees).
        angles : 3-item iterable = None
            These are the box angles (if lengths contains only 3 elements) in
            degrees. Must be a 1-D array of 3 elements or None if lengths
            includes angles as well.
        """
        def strip_units(x, desired_units):
            if u.is_quantity(x): return x.value_in_unit(desired_units)
            return x
        if len(lengths) == 3 and angles is None:
            raise ValueError('Both lengths and angles are required.')
        if len(lengths) == 6 and angles is not None:
            raise ValueError('Angles can be provided only once.')
        if len(lengths) != 6 and (len(lengths) != 3 or len(angles) != 3):
            raise ValueError('6 numbers expected -- 3 lengths and 3 angles.')
        if angles is None:
            angles = [strip_units(x, u.degrees) for x in lengths[3:]]
        lengths = [strip_units(x, u.angstroms) for x in lengths[:3]]
        self._ncfile.variables['cell_lengths'][self._last_box_frame] = \
                np.asarray(lengths)
        self._ncfile.variables['cell_angles'][self._last_box_frame] = \
                np.asarray(angles)
        self._last_box_frame += 1
        self.flush()

    add_box = add_cell_lengths_angles

    def time(self, frame):
        """
        Get the time of a particular frame in the trajectory

        Parameters
        ----------
        frame : int
            Which snapshot to get (first snapshot is frame 0)

        Returns
        -------
        time : float
            The time of the given frame in ps
        """
        return self._ncfile.variables['time'][frame]

    def add_time(self, stuff):
        """ Adds the time to the current frame of the NetCDF file

        Parameters
        ----------
        stuff : float or time-dimension Quantity
            The time to add to the current frame
        """
        if u.is_quantity(stuff): stuff = stuff.value_in_unit(u.picoseconds)
        self._ncfile.variables['time'][self._last_time_frame] = float(stuff)
        self._last_time_frame += 1
        self.flush()

    def remd_indices(self, frame):
        """ Returns the REMD indices for the desired frame

        Parameters
        ----------
        frame : int
            The frame to get the REMD indices from (0 is the first frame)
        """
        return self._ncfile.variables['remd_indices'][frame][:]

    def add_remd_indices(self, stuff):
        """ Add REMD indices to the current frame of the NetCDF file

        Parameters
        ----------
        stuff : iterable of int
            The indices in each REMD dimension
        """
        self._ncfile.variables['remd_indices'][self._last_remd_frame] = \
                np.asarray(stuff, dtype='i')
        self._last_remd_frame += 1
        self.flush()

    def temp0(self, frame):
        """ Returns the temperature of the current frame in kelvin

        Parameters
        ----------
        frame : int
            The frame to get the temperature from (0 is the first frame)
        """
        return self._ncfile.variables['temp0'][frame]

    def add_temp0(self, stuff):
        """ The temperature to add to the current frame of the NetCDF file

        Parameters
        ----------
        stuff : float or temperature Quantity
            The temperature to add to the current NetCDF file
        """
        if u.is_quantity(stuff): stuff = stuff.value_in_unit(u.kelvin)
        self._ncfile.variables['temp0'][self._last_remd_frame] = float(stuff)
        self._last_remd_frame += 1
        self.flush()

    @property
    def remd_dimtype(self):
        return self._ncfile.variables['remd_dimtype'][:]

    @remd_dimtype.setter
    def remd_dimtype(self, stuff):
        self._ncfile.variables['remd_dimtype'][:] = np.asarray(stuff, dtype='i')
        self.flush()

    def close(self):
        """ Closes the NetCDF file """
        self._ncfile.close()
        self.closed = True

    def __del__(self):
        self.closed or (hasattr(self, '_ncfile') and self._ncfile.close())

    def flush(self):
        try:
            self._ncfile.flush()
        except AttributeError:
            pass
