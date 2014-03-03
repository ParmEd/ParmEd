"""
This module contains classes for reading and writing Amber NetCDF-style files,
including both restarts and trajectories. These sets of classes abstract the
interaction with the various NetCDF-Python APIs that are available, namely
   o  netCDF4
   o  scipy
   o  ScientificPython
   o  pynetcdf

This module contains objects relevant to Amber NetCDF files. The code in
chemistry/amber/__init__.py is responsible for selecting the API based on a
default choice or user-selection (the latter is really only helpful for
development to ensure that all packages work correctly---there is no difference
from a user perspective). ALL NetCDF-file manipulation that the chemistry/amber
package does should be contained in this module.
"""

from chemistry import amber, __version__
from compat24 import property
try:
   import numpy as np
except ImportError:
   # This is just to prevent NetCDF imports from bringing everything down if
   # numpy is not available. np can be used inside any method since it will be
   # required for any NetCDF to work in the first place.
   np = None
# Make sure we have registered a NetCDF implementation. If not, take the default
if not amber.NETCDF_INITIALIZED: amber.use()

def needs_netcdf(fcn):
   """
   Decorator to protect against functions that need NetCDF so we can provide
   a helpful error message
   """
   def new_fcn(*args, **kwargs):
      if not amber.HAS_NETCDF:
         raise ImportError('No NetCDF packages are available!')
      return fcn(*args, **kwargs)
   return new_fcn


class NetCDFRestart(object):
   """ Class to read or write NetCDF restart files """

   @needs_netcdf
   def __init__(self, fname, mode):
      """
      Opens a NetCDF File. The main constructor should never be called directly.
      The alternative "open_old" and "open_new" constructors should be called
      instead for reading existing or writing new NetCDF files, respectively.
      """
      self.closed = False
      self._ncfile = amber.open_netcdf(fname, mode)
   
   @classmethod
   def open_new(cls, fname, natom, box, vels, title='',
         remd=None, temp=None, remd_indices=None, remd_groups=None,
         remd_dimtypes=None):
      """
      Opens a new NetCDF file and sets the attributes

      Parameters:

         fname (string): Name of the new file to open (overwritten)
         natom (int): Number of atoms in the restart
         box (bool): Write box information or not?
         vels (bool): Write velocities or not?
         title (string): title of the NetCDF file
         remd (string): Whether REMD information needs to be written. Can be
            None - No REMD information is written
            T(emperature) - Target temperature (or pH) will be written
            M(ulti-dimensional) - remd_indices and remd_groups will be written
         remd_dimtypes (int array): Array of exchange types for each group.

      remd is case-insensitive, and done based on first-letter matching, only
      """
      if remd is not None:
         if remd[0] in 'tT':
            remd_type = 'TEMPERATURE'
            if temp is None:
               raise ValueError('temp must be specified for T-REMD restarts.')
         elif remd[0] in 'mM':
            remd_type = 'MULTI'
            if (remd_indices is None or remd_groups is None or 
                  len(remd_indices) != len(remd_groups)):
               raise ValueError('remd_indices and remd_groups must be given '
                     'for multi-D REMD, and must have the same length.')
            remd_dimension = len(remd_indices)
         else:
            raise ValueError('remd must be None, T[emperature] or M[ulti]')
      else:
         remd_type = None
      inst = cls(fname, 'w')
      inst.hasbox = bool(box)
      inst.hasvels = bool(vels)
      # Assign the main attributes
      inst._ncfile.Conventions = 'AMBERRESTART'
      inst._ncfile.ConventionVersion = "1.0"
      inst._ncfile.title = title
      inst._ncfile.application = "AMBER"
      inst._ncfile.program = "ParmEd"
      inst._ncfile.programVersion = str(__version__)
      # Make all of the dimensions
      inst._ncfile.createDimension('spatial', 3)
      inst._ncfile.createDimension('atom', natom)
      inst.spatial = 3
      inst.atom = natom
      inst.title = inst._ncfile.title
      if box:
         inst._ncfile.createDimension('cell_spatial', 3)
         inst._ncfile.createDimension('label', 5)
         inst._ncfile.createDimension('cell_angular', 3)
         inst.cell_spatial = 3
         inst.label = 5
         inst.cell_angular = 3
      if remd_type == 'MULTI':
         inst._ncfile.createDimension('remd_dimension', remd_dimension)
         inst.remd_dimension = remd_dimension
      inst._ncfile.createDimension('time', 1)
      # Now make the variables
      v = inst._ncfile.createVariable('time', 'd', ('time',))
      v.units = 'picosecond'
      v = inst._ncfile.createVariable('spatial', 'c', ('spatial',))
      v[:] = np.asarray(list('xyz'))
      if inst.hasbox:
         v = inst._ncfile.createVariable('cell_angular', 'c', ('cell_angular', 'label'))
         v[0] = np.asarray(list('alpha'))
         v[1] = np.asarray(list('beta '))
         v[2] = np.asarray(list('gamma'))
         v = inst._ncfile.createVariable('cell_spatial', 'c', ('cell_spatial',))
         v[0], v[1], v[2] = 'a', 'b', 'c'
         v = inst._ncfile.createVariable('cell_lengths', 'd', ('cell_spatial',))
         v.units = 'angstrom'
         v = inst._ncfile.createVariable('cell_angles', 'd', ('cell_angular',))
         v.units = 'degree'
      v = inst._ncfile.createVariable('coordinates', 'd', ('atom', 'spatial'))
      v.units = 'angstrom'
      if inst.hasvels:
         v = inst._ncfile.createVariable('velocities', 'd', ('atom', 'spatial'))
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
         v = inst._ncfile.createVariable('temp0', 'd', ('time',))
         v.units = 'kelvin'
      elif remd_type == 'MULTI':
         inst._ncfile.createVariable('remd_indices', 'i', ('remd_dimension',))
         inst._ncfile.createVariable('remd_groups', 'i', ('remd_dimension',))
         inst._ncfile.createVariable('remd_dimtype', 'i', ('remd_dimension',))

      return inst

   @classmethod
   def open_old(cls, fname):
      """
      Opens the NetCDF file and sets the global attributes that the file sets
      """
      inst = cls(fname, 'r')
      inst.Conventions = inst._ncfile.Conventions
      inst.ConventionVersion = inst._ncfile.ConventionVersion
      inst.application = inst._ncfile.application
      inst.program = inst._ncfile.program
      inst.programVersion = inst._ncfile.programVersion
      inst.title = inst._ncfile.title
      # Set up the dimensions as attributes
      for dim in inst._ncfile.dimensions:
         if dim == 'time': continue # Exception for ParmEd-created ncrst files
         setattr(inst, dim, amber.get_int_dimension(inst._ncfile, dim))
      inst.hasvels = 'velocities' in inst._ncfile.variables
      inst.hasbox = ('cell_lengths' in inst._ncfile.variables and
                     'cell_angles' in inst._ncfile.variables)
      if inst.hasvels:
         inst.velocity_scale = inst._ncfile.variables['velocities'].scale_factor
         # Turn off automatic scaling for netCDF4. Ugh.
         try:
            inst._ncfile.variables['velocities'].set_auto_maskandscale(False)
         except AttributeError:
            # Other packages do not have this variable
            pass
      return inst

   @property
   def coordinates(self):
      return self._ncfile.variables['coordinates'][:].flatten()

   @coordinates.setter
   def coordinates(self, stuff):
      self._ncfile.variables['coordinates'][:] = np.reshape(stuff, (self.atom, 3))

   @property
   def velocities(self):
      return self._ncfile.variables['velocities'][:].flatten() * \
            self.velocity_scale

   @velocities.setter
   def velocities(self, stuff):
      self._ncfile.variables['velocities'][:] = \
            np.reshape(stuff, (self.atom, 3)) / self.velocity_scale

   @property
   def cell_lengths(self):
      return self._ncfile.variables['cell_lengths'][:]
   
   @cell_lengths.setter
   def cell_lengths(self, stuff):
      self._ncfile.variables['cell_lengths'][:] = np.asarray(stuff)

   @property
   def cell_angles(self):
      return self._ncfile.variables['cell_angles'][:]
   
   @cell_angles.setter
   def cell_angles(self, stuff):
      self._ncfile.variables['cell_angles'][:] = np.asarray(stuff)

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
      return amber.get_float(self._ncfile, 'time')

   @time.setter
   def time(self, stuff):
      self._ncfile.variables['time'][0] = float(stuff)

   @property
   def temp0(self):
      return amber.get_float(self._ncfile, 'temp0')

   @temp0.setter
   def temp0(self, stuff):
      self._ncfile.variables['temp0'][0] = float(stuff)
   
   @property
   def remd_indices(self):
      return self._ncfile.variables['remd_indices'][:]

   @remd_indices.setter
   def remd_indices(self, stuff):
      self._ncfile.variables['remd_indices'][:] = np.asarray(stuff, dtype='i')

   @property
   def remd_groups(self):
      return self._ncfile.variables['remd_groups'][:]

   @remd_groups.setter
   def remd_groups(self, stuff):
      self._ncfile.variables['remd_groups'][:] = np.asarray(stuff, dtype='i')

   @property
   def remd_dimtype(self):
      return self._ncfile.variables['remd_dimtype'][:]

   @remd_dimtype.setter
   def remd_dimtype(self, stuff):
      self._ncfile.variables['remd_dimtype'][:] = np.asarray(stuff, dtype='i')

   def close(self):
      self.closed = True
      self._ncfile.close()

   def __del__(self):
      self.closed or (hasattr(self, '_ncfile') and self._ncfile.close())

class NetCDFTraj(object):
   """ Class to read or write NetCDF restart files """

   @needs_netcdf
   def __init__(self, fname, mode):
      """ Opens a NetCDF File """
      self.closed = False
      self._ncfile = amber.open_netcdf(fname, mode)
   
   @classmethod
   def open_new(cls, fname, natom, box, crds=True, vels=False, frcs=False,
                remd=None, remd_dimension=None, title=''):
      """
      Opens a new NetCDF file and sets the attributes

      Parameters:

         fname (string): Name of the new file to open (overwritten)
         natom (int): Number of atoms in the restart
         box (bool): Write box information or not?
         crds (bool): Write coordinates or not?
         vels (bool): Write velocities or not?
         frcs (bool): Write forces or not?
         remd (string): What kind of REMD are we doing. Allowed values are:
            None - No REMD
            T[emperature] - Temperature or Constant pH REMD
            M[ulti-D] - Multi-dimensional REMD
         remd_dimension (int): Number of dimensions for multi-D REMD. None for
                               non-multi-D REMD
         title (string): title of the NetCDF file
      """
      inst = cls(fname, 'w')
      if remd is not None:
         if remd[0] in 'Tt':
            inst.remd = 'TEMPERATURE'
         elif remd[0] in 'Mm':
            inst.remd = 'MULTI'
            if remd_dimension is None:
               raise ValueError('remd_dimension must be given for multi-D REMD')
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
      inst._ncfile.Conventions = "AMBER"
      inst._ncfile.ConventionVersion = "1.0"
      inst._ncfile.application = "AmberTools"
      inst._ncfile.program = "ParmEd"
      inst._ncfile.programVersion = __version__
      inst._ncfile.title = "ParmEd-created trajectory"
      inst.Conventions = "AMBER"
      inst.ConventionVersion = "1.0"
      inst.application = "AmberTools"
      inst.program = "ParmEd"
      inst.programVersion = __version__
      inst.title = inst._ncfile.title
      # Create the dimensions
      inst._ncfile.createDimension('frame', None)
      inst._ncfile.createDimension('spatial', 3)
      inst._ncfile.createDimension('atom', natom)
      if inst.remd == 'MULTI':
         inst._ncfile.createDimension('remd_dimension', inst.remd_dimension)
      inst.frame, inst.spatial, inst.atom = None, 3, natom
      if inst.hasbox:
         inst._ncfile.createDimension('cell_spatial', 3)
         inst._ncfile.createDimension('cell_angular', 3)
         inst._ncfile.createDimension('label', 5)
         inst.cell_spatial, inst.cell_angular, inst.label = 3, 3, 5
      # Create the variables and assign units and scaling factors
      v = inst._ncfile.createVariable('spatial', 'c', ('spatial',))
      v[:] = np.asarray(list('xyz'))
      if inst.hasbox:
         v = inst._ncfile.createVariable('cell_spatial', 'c', ('cell_spatial',))
         v[:] = np.asarray(list('abc'))
         v = inst._ncfile.createVariable('cell_angular', 'c',
               ('cell_angular', 'label',))
         v[:] = np.asarray([list('alpha'), list('beta '), list('gamma')])
      v = inst._ncfile.createVariable('time', 'f', ('frame',))
      v.units = 'picosecond'
      if inst.hascrds:
         v = inst._ncfile.createVariable('coordinates', 'f',
               ('frame', 'atom', 'spatial'))
         v.units = 'angstrom'
         inst._last_crd_frame = 0
      if inst.hasvels:
         v = inst._ncfile.createVariable('velocities', 'f',
               ('frame', 'atom', 'spatial'))
         v.units = 'angstrom/picosecond'
         inst.velocity_scale = v.scale_factor = 20.455
         inst._last_vel_frame = 0
      if inst.hasfrcs:
         v = inst._ncfile.createVariable('forces', 'f',
               ('frame', 'atom', 'spatial'))
         v.units = 'kilocalorie/mole/angstrom'
         inst._last_frc_frame = 0
      if inst.hasbox:
         v = inst._ncfile.createVariable('cell_lengths', 'd',
               ('frame', 'cell_spatial'))
         v.units = 'angstrom'
         v = inst._ncfile.createVariable('cell_angles', 'd',
               ('frame', 'cell_angular'))
         v.units = 'degree'
         inst._last_box_frame = 0
      if inst.remd == 'TEMPERATURE':
         v = inst._ncfile.createVariable('temp0', 'd', ('frame',))
         v.units = 'kelvin'
         inst._last_remd_frame = 0
      elif inst.remd == 'MULTI':
         inst._ncfile.createVariable('remd_indices', 'i',
               ('frame', 'remd_dimension'))
         inst._ncfile.createVariable('remd_dimtype', 'i', ('remd_dimension',))
         inst._last_remd_frame = 0

      inst._last_time_frame = 0

      return inst

   @classmethod
   def open_old(cls, fname):
      """
      Opens the NetCDF file and sets the global attributes that the file sets

      Parameters:
         -  fname (string): File name of the trajectory to open. It must exist
      """
      inst = cls(fname, 'r')
      inst.Conventions = inst._ncfile.Conventions
      inst.ConventionVersion = inst._ncfile.ConventionVersion
      inst.application = inst._ncfile.application
      inst.program = inst._ncfile.program
      inst.programVersion = inst._ncfile.programVersion
      inst.title = inst._ncfile.title
      # Set up the dimensions as attributes
      for dim in inst._ncfile.dimensions:
         setattr(inst, dim, amber.get_int_dimension(inst._ncfile, dim))
      inst.hascrds = 'coordinates' in inst._ncfile.variables
      inst.hasvels = 'velocities' in inst._ncfile.variables
      inst.hasfrcs = 'forces' in inst._ncfile.variables
      inst.hasbox = ('cell_lengths' in inst._ncfile.variables and
                     'cell_angles' in inst._ncfile.variables)
      if inst.hasvels:
         inst.velocity_scale = inst._ncfile.variables['velocities'].scale_factor
      if inst.frame is None:
         if 'time' in inst._ncfile.variables:
            inst.frame = len(inst._ncfile.variables['time'][:])
         elif inst.hascrds:
            inst.frame = len(inst._ncfile.variables['coordinates'][:])
         elif inst.hasvels:
            inst.frame = len(inst._ncfile.variables['velocities'][:])
         elif inst.hasfrcs:
            inst.frame = len(inst._ncfile.variables['forces'][:])
      return inst

   def coordinates(self, frame):
      """
      Get the coordinates of a particular frame in the trajectory

      Parameters:
         -  frame (int): Which snapshot to get (first snapshot is frame 0)

      Returns:
         numpy array of length 3*natom with the given coordinates
      """
      return self._ncfile.variables['coordinates'][frame][:].flatten()

   def add_coordinates(self, stuff):
      """
      Adds a new coordinate frame to the end of a NetCDF trajectory. This should
      only be called on objects created with the "open_new" constructor.

      Parameters:
         -  stuff (array): This array of floats is converted into a numpy array
            of shape (natom, 3). It can be passed either in the 2-D format of 
            [ [x1, y1, z1], [x2, y2, z2], ... ] or in the 1-D format of
            [x1, y1, z1, x2, y2, z2, ... ].
      """
      if not isinstance(stuff, np.ndarray): stuff = np.asarray(stuff)
      self._ncfile.variables['coordinates'][self._last_crd_frame] = \
            np.reshape(stuff, (self.atom, 3))
      self._last_crd_frame += 1

   def velocities(self, frame):
      """
      Get the velocities of a particular frame in the trajectory

      Parameters:
         -  frame (int): Which snapshot to get (first snapshot is frame 0)

      Returns:
         numpy array of length 3*atom with the given velocities properly
         scaled to be of units angstrom/picosecond
      """
      return (self._ncfile.variables['velocities'][frame][:].flatten() * 
              self.velocity_scale)

   def add_velocities(self, stuff):
      """
      Adds a new velocities frame to the end of a NetCDF trajectory. This should
      only be called on objects created with the "open_new" constructor.

      Parameters:
         -  stuff (array): This array of floats is converted into a numpy array
            of shape (natom, 3). It can be passed either in the 2-D format of 
            [ [x1, y1, z1], [x2, y2, z2], ... ] or in the 1-D format of
            [x1, y1, z1, x2, y2, z2, ... ].
      """
      if not isinstance(stuff, np.ndarray): stuff = np.asarray(stuff)
      self._ncfile.variables['velocities'][self._last_vel_frame] = \
            np.reshape(stuff, (self.atom, 3)) / self.velocity_scale
      self._last_vel_frame += 1

   def forces(self, frame):
      """
      Get the forces of a particular frame in the trajectory

      Parameters:
         -  frame (int): Which snapshot to get (first snapshot is frame 0)

      Returns:
         numpy array of length 3*atom with the given forces properly
         scaled to be of units amu*angstrom/picosecond^2
      """
      return (self._ncfile.variables['forces'][frame][:].flatten())

   def add_forces(self, stuff):
      """
      Adds a new coordinate frame to the end of a NetCDF trajectory. This should
      only be called on objects created with the "open_new" constructor.

      Parameters:
         -  stuff (array): This array of floats is converted into a numpy array
            of shape (natom, 3). It can be passed either in the 2-D format of 
            [ [x1, y1, z1], [x2, y2, z2], ... ] or in the 1-D format of
            [x1, y1, z1, x2, y2, z2, ... ].
      """
      if not isinstance(stuff, np.ndarray): stuff = np.asarray(stuff)
      self._ncfile.variables['forces'][self._last_frc_frame] = \
            np.reshape(stuff, (self.atom, 3))
      self._last_frc_frame += 1

   def cell_lengths_angles(self, frame):
      """
      Get the cell lengths and cell angles of a particular frame in the
      trajectory

      Parameters:
         -  frame (int): Which snapshot to get (first snapshot is frame 0)

      Returns:
         2-element tuple: (length-3 numpy array of cell lengths, length-3 numpy
         array of cell angles)
      """
      return (self._ncfile.variables['cell_lengths'][frame][:],
              self._ncfile.variables['cell_angles'][frame][:])
   
   def add_cell_lengths_angles(self, lengths, angles=None):
      """
      Adds a new cell length and angle frame to the end of a NetCDF trajectory.
      This should only be called on objects created with the "open_new"
      constructor.

      Parameters:
         -  lengths (array): This should be a 1-D array of 3 or 6 elements. If 6
               elements, angles should be None and the first 3 elements are the
               box lengths (angstroms) and the last 3 are the box angles
               (degrees).
         -  angles (array): These are the box angles (if lengths contains only 3
               elements) in degrees. Must be a 1-D array of 3 elements or None
               if lengths includes angles as well.
      """
      if len(lengths) == 3 and angles is None:
         raise ValueError('Both lengths and angles are required.')
      if len(lengths) == 6 and angles is not None:
         raise ValueError('Angles can be provided only once.')
      if len(lengths) != 6 and (len(lengths) != 3 or len(angles) != 3):
         raise ValueError('6 numbers expected -- 3 lengths and 3 angles.')
      if angles is None:
         angles = lengths[3:]
      lengths = lengths[:3]
      self._ncfile.variables['cell_lengths'][self._last_box_frame] = \
            np.asarray(lengths)
      self._ncfile.variables['cell_angles'][self._last_box_frame] = \
            np.asarray(angles)
      self._last_box_frame += 1

   def time(self, frame):
      """
      Get the time of a particular frame in the trajectory

      Parameters:
         -  frame (int): Which snapshot to get (first snapshot is frame 0)

      Returns:
         float: time of the given frame
      """
      return self._ncfile.variables['time'][frame]

   def add_time(self, stuff):
      self._ncfile.variables['time'][self._last_time_frame] = float(stuff)
      self._last_time_frame += 1

   def remd_indices(self, frame):
      return self._ncfile.variables['remd_indices'][frame][:]

   def add_remd_indices(self, stuff):
      self._ncfile.variables['remd_indices'][self._last_remd_frame] = \
            np.asarray(stuff, dtype='i')
      self._last_remd_frame += 1

   def temp0(self, frame):
      return self._ncfile.variables['temp0'][frame]

   def add_temp0(self, stuff):
      self._ncfile.variables['temp0'][self._last_remd_frame] = float(stuff)
      self._last_remd_frame += 1

   @property
   def remd_dimtype(self):
      return self._ncfile.variables['remd_dimtype'][:]

   @remd_dimtype.setter
   def remd_dimtype(self, stuff):
      self._ncfile.variables['remd_dimtype'][:] = np.asarray(stuff, dtype='i')

   def close(self):
      """ Closes the NetCDF file """
      self._ncfile.close()
      self.closed = True

   def __del__(self):
      self.closed or (hasattr(self, '_ncfile') and self._ncfile.close())
