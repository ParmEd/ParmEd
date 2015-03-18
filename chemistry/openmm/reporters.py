from __future__ import division

from chemistry.amber.asciicrd import AmberMdcrd
from chemistry.geometry import box_vectors_to_lengths_and_angles
from chemistry.amber.netcdffiles import NetCDFTraj
from chemistry.amber.readparm import Rst7
from chemistry import unit as u
from functools import wraps
from math import isnan, isinf
try:
    import simtk.openmm as mm
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
from time import time as timeofday

def needs_openmm(fcn):
    global HAS_OPENMM
    @wraps(fcn)
    def new_fcn(*args, **kwargs):
        if not HAS_OPENMM:
            raise ImportError('Could not find or import OpenMM')
        return fcn(*args, **kwargs)

    return new_fcn

VELUNIT = u.angstrom / u.picosecond
FRCUNIT = u.kilocalorie_per_mole / u.angstrom

class StateDataReporter(object):
    """
    This class acts as a state data reporter for OpenMM simulations, but it is a
    little more generalized. Notable differences are:

      -  It allows the units of the output to be specified, with defaults being
         those used in Amber (e.g., kcal/mol, angstroms^3, etc.)
      -  It will write to any object containing a 'write' method; not just
         files. This allows, for instance, writing to a GUI window that
         implements the desired 'write' attribute.

    Most of this code is copied from the OpenMM StateDataReporter class, with
    the above-mentioned changes made.
    
    Parameters
    ----------
    f : str or file-like
        Destination to write the state data (file name or file object)
    reportInterval : int
        Number of steps between state data reports
    step : bool, optional
        Print out the step number (Default True)
    time : bool, optional
        Print out the simulation time (Defaults True)
    potentialEnergy : bool, optional
        Print out the potential energy of the structure (Default True)
    kineticEnergy : bool, optional
        Print out the kinetic energy of the structure (Default True)
    totalEnergy : bool, optional
        Print out the total energy of the system (Default True)
    temperature : bool, optional
        Print out the temperature of the system (Default True)
    volume : bool, optional
        Print out the volume of the unit cell. If the system is not periodic,
        the value is meaningless (Default False)
    density : bool, optional
        Print out the density of the unit cell. If the system is not periodic,
        the value is meaningless (Default False)
    separator : str, optional
        The string to separate data fields (Default ',')
    systemMass : float, optional
        If not None, the density will be computed from this mass, since setting
        a mass to 0 is used to constrain the position of that particle. (Default
        None)
    unit_system : :class:`UnitSystem`, optional
        The unit system to which to convert all of the resulting quantities to.
        Default is ``akma_unit_system``.
    """

    @needs_openmm
    def __init__(self, f, reportInterval, step=True, time=True,
                 potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                 temperature=True, volume=False, density=False, separator=',',
                 systemMass=None, unit_system=u.akma_unit_system):
        """ Create a StateDataReporter.  """

        self._reportInterval = reportInterval
        self._openedFile = not hasattr(f, 'write')
        if self._openedFile:
            self._out = open(f, 'w')
        else:
            self._out = f
        self._step = step
        self._time = time
        self._potentialEnergy = potentialEnergy
        self._kineticEnergy = kineticEnergy
        self._totalEnergy = totalEnergy
        self._temperature = temperature
        self._volume = volume
        self._density = density
        self._separator = separator
        self._totalMass = systemMass
        self._hasInitialized = False
        self._needsPositions = False
        self._needsVelocities = False
        self._needsForces = False
        self._needEnergy = (potentialEnergy or kineticEnergy or
                            totalEnergy or temperature)
        self.unit_system = unit_system

    def describeNextReport(self, simulation):
        """
        Get information about the next report this object will generate.

        Parameters:
            - simulation (Simulation) The Simulation to generate a report for
        Returns:
            A five element tuple.  The first element is the number of steps
            until the next report.  The remaining elements specify whether that
            report will require positions, velocities, forces, and energies
            respectively.
        """
        steps_left = simulation.currentStep % self._reportInterval
        steps = self._reportInterval - steps_left
        return (steps, self._needsPositions, self._needsVelocities,
                self._needsForces, self._needEnergy)

    def report(self, simulation, state):
        """Generate a report.

        Parameters:
            - simulation (Simulation) The Simulation to generate a report for
            - state (State) The current state of the simulation
        """
        if not self._hasInitialized:
            self._initializeConstants(simulation)
            headers = self._constructHeaders()
            self._out.write('#"%s"\n' % ('"'+self._separator+'"').join(headers))
            self._hasInitialized = True

        # Check for errors.
        self._checkForErrors(simulation, state)

        # Query for the values
        values = self._constructReportValues(simulation, state)

        # Write the values.
        self._out.write(self._separator.join(str(v) for v in values) + '\n')
        hasattr(self._out, 'flush') and self._out.flush()

    def _constructReportValues(self, simulation, state):
        """
        Query the simulation for the current state of our observables of
        interest.

        Parameters
        ----------
        simulation : Simulation
            The Simulation object to generate a report for
        state : State
            The current state of the simulation object

        Returns: A list of values summarizing the current state of
        the simulation, to be printed or saved. Each element in the list
        corresponds to one of the columns in the resulting CSV file.
        """
        values = []
        if self._volume or self._density:
            volume = state.getPeriodicBoxVolume()
        if self._density:
            density = self._totalMass / volume
        ke = state.getKineticEnergy()
        pe = state.getPotentialEnergy()
        if self._temperature:
            temp = 2 * ke / (self._dof * u.MOLAR_GAS_CONSTANT_R)
        time = state.getTime()
        if self._step:
            values.append(simulation.currentStep)
        if self._time:
            if self.unit_system is u.akma_unit_system:
                values.append(time.value_in_unit(u.picoseconds))
            else:
                values.append(time.value_in_unit_system(self.unit_system))
        if self._potentialEnergy:
            values.append(pe.value_in_unit_system(self.unit_system))
        if self._kineticEnergy:
            values.append(ke.value_in_unit_system(self.unit_system))
        if self._totalEnergy:
            values.append((pe + ke).value_in_unit_system(self.unit_system))
        if self._temperature:
            values.append(temp.value_in_unit_system(self.unit_system))
        if self._volume:
            values.append(volume.value_in_unit_system(self.unit_system))
        if self._density:
            values.append(density.value_in_unit_system(self.densityUnit))
        return values

    def _initializeConstants(self, simulation):
        """
        Initialize a set of constants required for the reports

        Parameters
        ----------
        simulation : Simulation
            The simulation to generate a report for
        """
        system = simulation.system
        frclist = system.getForces()
        if self._temperature:
            # Compute the number of degrees of freedom.
            dof = 0
            for i in xrange(system.getNumParticles()):
                if system.getParticleMass(i) > 0*u.dalton:
                    dof += 3
            dof -= system.getNumConstraints()
            if any(isinstance(frc, mm.CMMotionRemover) for frc in frclist):
                dof -= 3
            self._dof = dof
        if self._density:
            if self._totalMass is None:
                # Compute the total system mass.
                self._totalMass = 0*u.dalton
                for i in xrange(system.getNumParticles()):
                    self._totalMass += system.getParticleMass(i)
            elif not u.is_quantity(self._totalMass):
                self._totalMass = self._totalMass*u.dalton

    def _constructHeaders(self):
        """Construct the headers for the CSV output

        Returns: a list of strings giving the title of each observable being
                 reported on.
        """
        headers = []
        if self._step:
            headers.append('Step')
        if self._time:
            headers.append('Time (ps)')
        if self._potentialEnergy:
            headers.append('Potential Energy (%s)' % self.energyUnit)
        if self._kineticEnergy:
            headers.append('Kinetic Energy (%s)' % self.energyUnit)
        if self._totalEnergy:
            headers.append('Total Energy (%s)' % self.energyUnit)
        if self._temperature:
            headers.append('Temperature (%s)' % self.temperatureUnit)
        if self._volume:
            headers.append('Box Volume (%s)' %
                           str(self.volumeUnit).replace('**','^'))
        if self._density:
            headers.append('Density (%s)' %
                           str(self.densityUnit).replace('item*', ''))
        return headers

    def _checkForErrors(self, simulation, state):
        """Check for errors in the current state of the simulation

        Parameters
            - simulation (Simulation) The Simulation to generate a report for
            - state (State) The current state of the simulation
        """
        if self._needEnergy:
            energy = state.getKineticEnergy() + state.getPotentialEnergy()
            energy = energy.value_in_unit(u.kilojoules_per_mole)
            if isnan(energy):
                raise ValueError('Energy is NaN')
            if isinf(energy):
                raise ValueError('Energy is infinite')

    def __del__(self):
        if hasattr(self, '_openedFile') and self._openedFile:
            self._out.close()

    def finalize(self):
        """ Closes any open file """
        try:
            if self._out is not None and self._openedFile:
                self._out.close()
        except NameError:
            pass

class NetCDFReporter(object):
    """ NetCDFReporter prints a trajectory in NetCDF format """

    @needs_openmm
    def __init__(self, file, reportInterval, crds=True, vels=False, frcs=False):
        """
        Create a NetCDFReporter instance.

        Parameters
        ----------
        file : str
            Name of the file to write the trajectory to
        reportInterval : int
            How frequently to write a frame to the trajectory
        crds : bool=True
            Should we write coordinates to this trajectory? (Default True)
        vels : bool=False
            Should we write velocities to this trajectory? (Default False)
        frcs : bool=False
            Should we write forces to this trajectory? (Default False)
        """
        if not crds and not vels and not frcs:
            raise ValueError('You must print either coordinates, velocities, '
                             'or forces in a NetCDFReporter')
        # Control flags
        self.crds, self.vels, self.frcs = crds, vels, frcs
        if not (crds or vels or frcs):
            raise ValueError('Cannot write a trajectory without coordinates, '
                             'velocities, or forces! Pick at least one.')
        self._reportInterval = reportInterval
        self._out = None # not written yet
        self.fname = file

    def describeNextReport(self, simulation):
        """
        Get information about the next report this object will generate.

        Parameters:
            - simulation (Simulation) The Simulation to generate a report for

        Returns:
            A five element tuple.  The first element is the number of steps
            until the next report.  The remaining elements specify whether that
            report will require positions, velocities, forces, and energies
            respectively.
        """
        stepsleft = simulation.currentStep % self._reportInterval
        steps = self._reportInterval - stepsleft
        return (steps, self.crds, self.vels, self.frcs, False)


    def report(self, simulation, state):
        """Generate a report.

        Parameters:
            - simulation (Simulation) The Simulation to generate a report for
            - state (State) The current state of the simulation
        """
        global VELUNIT, FRCUNIT
        if self.crds:
            crds = state.getPositions().value_in_unit(u.angstrom)
        if self.vels:
            vels = state.getVelocities().value_in_unit(VELUNIT)
        if self.frcs:
            frcs = state.getForces().value_in_unit(FRCUNIT)
        if self._out is None:
            # This must be the first frame, so set up the trajectory now
            if self.crds:
                atom = len(crds)
            elif self.vels:
                atom = len(vels)
            elif self.frcs:
                atom = len(frcs)
            self.uses_pbc = simulation.topology.getUnitCellDimensions() is not None
            self._out = NetCDFTraj.open_new(
                    self.fname, atom, self.uses_pbc, self.crds, self.vels,
                    self.frcs, title="ParmEd-created trajectory using OpenMM"
            )
        if self.uses_pbc:
            vecs = state.getPeriodicBoxVectors()
            lengths, angles = box_vectors_to_lengths_and_angles(*vecs)
            self._out.add_cell_lengths_angles(lengths.value_in_unit(u.angstrom),
                                              angles.value_in_unit(u.degree))

        # Add the coordinates, velocities, and/or forces as needed
        if self.crds:
            self._out.add_coordinates(crds)
        if self.vels:
            # The velocities get scaled right before writing
            self._out.add_velocities(vels)
        if self.frcs:
            self._out.add_forces(frcs)
        # Now it's time to add the time.
        self._out.add_time(state.getTime().value_in_unit(u.picosecond))

    def __del__(self):
        try:
            if self._out is not None:
                self._out.close()
        except NameError:
            pass

    def finalize(self):
        """ Closes any open file """
        try:
            if self._out is not None:
                self._out.close()
        except NameError:
            pass

class MdcrdReporter(object):
    """
    MdcrdReporter prints a trajectory in ASCII Amber format. This reporter will
    be significantly slower than binary file reporters (like DCDReporter or
    NetCDFReporter).
    """

    @needs_openmm
    def __init__(self, file, reportInterval, crds=True, vels=False, frcs=False):
        """
        Create a MdcrdReporter instance.

        Parameters
        ----------
        file : str
            Name of the file to write the trajectory to
        reportInterval : int
            Number of steps between writing trajectory frames
        crds : bool=True
            Write coordinates to this trajectory file?
        vels : bool=False
            Write velocities to this trajectory file?
        frcs : bool=False
            Write forces to this trajectory file?

        Notes
        -----
        You can only write one of coordinates, forces, or velocities to a mdcrd
        file.
        """
        # ASCII mdcrds can have either coordinates, forces, or velocities
        ntrue = 0
        if crds: ntrue += 1
        if vels: ntrue += 1
        if frcs: ntrue += 1
        if ntrue != 1:
            raise ValueError('MdcrdReporter must print exactly one of either '
                             'coordinates, velocities, or forces.')
        # Control flags
        self.crds, self.vels, self.frcs = crds, vels, frcs
        self._reportInterval = reportInterval
        self._out = None # not written yet
        self.fname = file

    def describeNextReport(self, simulation):
        """
        Get information about the next report this object will generate.

        Parameters:
            - simulation (Simulation) The Simulation to generate a report for

        Returns:
            A five element tuple.  The first element is the number of steps
            until the next report.  The remaining elements specify whether that
            report will require positions, velocities, forces, and energies
            respectively.
        """
        stepsleft = simulation.currentStep % self._reportInterval
        steps = self._reportInterval - stepsleft
        return (steps, self.crds, self.vels, self.frcs, False)

    def report(self, simulation, state):
        """
        Generate a report.

        Parameters:
            - simulation (Simulation) The Simulation to generate a report for
            - state (State) The current state of the simulation
        """
        from chemistry.amber.asciicrd import VELSCALE
        global VELUNIT, FRCUNIT
        if self.crds:
            crds = state.getPositions().value_in_unit(u.angstrom)
        elif self.vels: # crds/vels/frcs are exclusive, elif works
            vels = state.getVelocities().value_in_unit(VELUNIT)
        elif self.frcs:
            frcs = state.getForces().value_in_unit(FRCUNIT)
        if self._out is None:
            # This must be the first frame, so set up the trajectory now
            if self.crds:
                self.atom = len(crds)
            elif self.vels:
                self.atom = len(vels)
            elif self.frcs:
                self.atom = len(frcs)
            self.uses_pbc = simulation.topology.getUnitCellDimensions() is not None
            self._out = AmberMdcrd(self.fname, self.atom, self.uses_pbc,
                    title="ParmEd-created trajectory using OpenMM", mode='w')

        # Add the coordinates, velocities, and/or forces as needed
        if self.crds:
            flatcrd = [0 for i in xrange(self.atom*3)]
            for i in xrange(self.atom):
                i3 = i*3
                flatcrd[i3], flatcrd[i3+1], flatcrd[i3+2] = crds[i]
            self._out.add_coordinates(flatcrd)
        if self.vels:
            # Divide by the scaling factor (works if vels is a list of Vec3's)
            # This is necessary since AmberMdcrd does not scale before writing
            # (since it expects coordinates)
            vels = [v / VELSCALE for v in vels]
            flatvel = [0 for i in xrange(self.atom*3)]
            for i in xrange(self.atom):
                i3 = i*3
                flatvel[i3], flatvel[i3+1], flatvel[i3+2] = vels[i]
            self._out.add_coordinates(flatvel)
        if self.frcs:
            flatfrc = [0 for i in xrange(self.atom*3)]
            for i in xrange(self.atom):
                i3 = i*3
                flatfrc[i3], flatfrc[i3+1], flatfrc[i3+2] = frcs[i]
            self._out.add_coordinates(flatfrc)
        # Now it's time to add the box lengths
        if self.uses_pbc:
            boxvecs = state.getPeriodicBoxVectors()
            lengths, angles = box_vectors_to_lengths_and_angles(*boxvecs)
            self._out.add_box(lengths.value_in_unit(u.angstroms))

    def __del__(self):
        try:
            if self._out is not None:
                self._out.close()
        except NameError:
            pass

    def finalize(self):
        """ Closes any open file """
        try:
            if self._out is not None:
                self._out.close()
        except NameError:
            pass

class RestartReporter(object):
    """
    Use a reporter to handle writing restarts at specified intervals.

    Parameters
    ----------
    file : str
        Name of the file to write the restart to.
    reportInterval : int
        Number of steps between writing restart files
    write_multiple : bool=False
        Either write a separate restart each time (appending the step number in
        the format .# to the file name given above) if True, or overwrite the
        same file each time if False
    netcdf : bool=False
        Use the Amber NetCDF restart file format
    write_velocities : bool=True
        Write velocities to the restart file. You can turn this off for passing
        in, for instance, a minimized structure.
    """
   
    @needs_openmm
    def __init__(self, file, reportInterval, write_multiple=False, netcdf=False,
                 write_velocities=True):
        self.fname = file
        self._reportInterval = reportInterval
        self.write_multiple = write_multiple
        self.netcdf = netcdf
        self.write_velocities = write_velocities
        self.rst7 = None

    def describeNextReport(self, simulation):
        """
        Get information about the next report this object will generate.

        Parameters:
            - simulation (Simulation) The Simulation to generate a report for

        Returns:
            A five element tuple.  The first element is the number of steps
            until the next report.  The remaining elements specify whether that
            report will require positions, velocities, forces, and energies
            respectively.
        """
        stepsleft = simulation.currentStep % self._reportInterval
        steps = self._reportInterval - stepsleft
        return (steps, True, True, False, False)

    def report(self, simulation, state):
        """Generate a report.

        Parameters:
            - simulation (Simulation) The Simulation to generate a report for
            - state (State) The current state of the simulation
        """
        global VELUNIT
        crds = state.getPositions().value_in_unit(u.angstrom)
        if self.rst7 is None:
            self.uses_pbc = simulation.topology.getUnitCellDimensions() is not None
            self.atom = len(crds)
            # First time written
            self.rst7 = Rst7(natom=self.atom,
                             title='Restart file written by ParmEd with OpenMM')
        self.rst7.time = state.getTime().value_in_unit(u.picosecond)
        flatcrd = [0.0 for i in xrange(self.atom*3)]
        for i in xrange(self.atom):
            i3 = i*3
            flatcrd[i3], flatcrd[i3+1], flatcrd[i3+2] = crds[i]
        self.rst7.coordinates = flatcrd

        if self.write_velocities:
            vels = state.getVelocities().value_in_unit(VELUNIT)
            flatvel = [0.0 for i in xrange(self.atom*3)]
            for i in xrange(self.atom):
                i3 = i*3
                flatvel[i3], flatvel[i3+1], flatvel[i3+2] = vels[i]
            self.rst7.vels = flatvel

        if self.uses_pbc:
            boxvecs = state.getPeriodicBoxVectors()
            lengths, angles = box_vectors_to_lengths_and_angles(*boxvecs)
            lengths = lengths.value_in_unit(u.angstrom)
            angles = angles.value_in_unit(u.degree)
            self.rst7.box = [lengths[0], lengths[1], lengths[2],
                             angles[0], angles[1], angles[2]]

        if self.write_multiple:
            fname = self.fname + '.%d' % simulation.currentStep
        else:
            fname = self.fname

        self.rst7.write(fname, self.netcdf)

    def finalize(self):
        """ No-op here """
        pass

class ProgressReporter(StateDataReporter):
    """
    A class that prints out a progress report of how much MD (or minimization)
    has been done, how fast the simulation is running, and how much time is left
    (similar to the mdinfo file in Amber)
    """

    @needs_openmm
    def __init__(self, f, reportInterval, totalSteps, potentialEnergy=True,
                 kineticEnergy=True, totalEnergy=True, temperature=False,
                 volume=False, density=False, systemMass=None, **kwargs):
        """ Create a StateDataReporter.
        Parameters:
        - f (str) The file name of the progress report (overwritten each time)
        - reportInterval (int) The step interval at which to write frames
        - totalSteps (int) The total number of steps that will be run in the
                           simulation (used to estimate the completion time)
        - potentialEnergy (boolean=False) Whether to write the potential energy
        - kineticEnergy (boolean=False) Whether to write the kinetic energy
        - totalEnergy (boolean=False) Whether to write the total energy
        - temperature (boolean=False) Whether to write the instantaneous temp.
        - volume (boolean=False) Whether to write the periodic box volume
        - density (boolean=False) Whether to write the system density
        - systemMass (mass=None) The total mass to use for the system when
            reporting density.  If this is None (the default), the system mass
            is computed by summing the masses of all particles.  This parameter
            is useful when the particle masses do not reflect their actual
            physical mass, such as when some particles have had their masses
            set to 0 to immobilize them.
        - **kwargs (keyword arguments): any other arguments accepted by
            StateDataReporter
        """
        # Make sure we got a file name rather than a file-like object.
        # Immediately close the file after opening. This erases it, which isn't
        # a bad thing, but also prepares it for reports
        kwargs['time'] = kwargs['step'] = True
        super(ProgressReporter, self).__init__(
                f, reportInterval, potentialEnergy=potentialEnergy,
                kineticEnergy=kineticEnergy, totalEnergy=totalEnergy,
                temperature=temperature, volume=volume, density=density,
                systemMass=systemMass, **kwargs
        )
        if not self._openedFile:
            raise ValueError('ProgressReporter requires a file name '
                             '(not file object)')
        self._out.close()
        del self._out
        self.fname = f
        self._totalSteps = totalSteps

        # For timing
        self._startTime = None
        self._firstStep = None
        self._lastReportTime = None
        self._timeStep = None

    def describeNextReport(self, simulation):
        """
        Get information about the next report this object will generate.

        Parameters:
            - simulation (Simulation) The Simulation to generate a report for

        Returns: A five element tuple.  The first element is the number of steps
            until the next report.  The remaining elements specify whether that
            report will require positions, velocities, forces, and energies
            respectively.
        """
        if self._startTime is None:
            # First time this is called, initialize the timers
            self._startTime = timeofday()
            self._lastReportTime = timeofday()
            self._timeStep = simulation.integrator.getStepSize()
            self._timeStep = self._timeStep.value_in_unit(u.nanosecond)
            self._firstStep = simulation.currentStep
        stepsleft = simulation.currentStep % self._reportInterval
        steps = self._reportInterval - stepsleft
        return (steps, False, False, False, self._needEnergy)
   
    def report(self, simulation, state):
        """
        Generate a report and predict the time to completion (and
        current/overall MD performance)
        """
        if not self._hasInitialized:
            self._initializeConstants(simulation)
            self._hasInitialized = True

        # Check for errors.
        self._checkForErrors(simulation, state)

        # Query for the values
        values = self._constructReportValues(simulation, state)

        now = timeofday()

        total_time = now - self._startTime
        partial_time = now - self._lastReportTime
        self._lastReportTime = now

        total_nsperday = ((values['step'] - self._firstStep) * self._timeStep /
                           total_time)
        partial_nsperday = (self._reportInterval*self._timeStep) / partial_time

        # Get total and partial ns/day (currently ns/second)
        total_nsperday *= 3600 * 24
        partial_nsperday *= 3600 * 24

        # Estimated time to completion (based on performance of last N steps
        remaining_steps = self._totalSteps - values['step'] + self._firstStep
        etc = partial_time / self._reportInterval * remaining_steps

        if etc > 3600:
            etc /= 3600
            unitstr = 'hr.'
        elif etc > 60:
            etc /= 60
            unitstr = 'min.'
        else:
            unitstr = 'sec.'

        # Write the values.
        with open(self.fname, 'w') as f:
            f.write('-+' * 39 + '\n')
            f.write('\n')
            f.write('  On step %d of %d\n' % (values['step'] - self._firstStep,
                                              self._totalSteps))
            f.write('\n')
            if self._totalEnergy:
                f.write('  Total Energy     = %12.4f\n' % values['totalEnergy'])
            if self._potentialEnergy:
                f.write('  Potential Energy = %12.4f\n' %
                        values['potentialEnergy'])
            if self._kineticEnergy:
                f.write('  Kinetic Energy   = %12.4f\n' %
                        values['kineticEnergy'])
            if self._volume:
                f.write('  Volume           = %12.4f\n' % values['volume'])
            if self._density:
                f.write('  Density          = %12.4f\n' % values['density'])
            if self._temperature:
                f.write('  Temperature      = %12.4f\n' % values['temperature'])
            f.write('\n')
            f.write(' Time for last %8d steps: %10.4f s. (%.3f ns/day)\n' %
                    (self._reportInterval, partial_time, partial_nsperday))
            f.write(' Time for all %9d steps: %10.4f s. (%.3f ns/day)\n' %
                    (values['step']-self._firstStep, total_time, total_nsperday)
            )
            f.write('\n')
            f.write(' Estimated time to completion: %.3f %s\n' % (etc, unitstr))
            f.write('\n')
            f.write('-+' * 39 + '\n')

        if remaining_steps == 0: self._startTime = None

    def _constructReportValues(self, simulation, state):
        """
        Query the simulation for the current state of our observables of
        interest.

        Parameters:
            - simulation (Simulation) The Simulation to generate a report for
            - state (State) The current state of the simulation

        Returns: A list of values summarizing the current state of the
            simulation, to be printed or saved. Each element in the list
            corresponds to one of the columns in the resulting CSV file.
        """
        values = dict()
        values['step'] = simulation.currentStep
        values['time'] = state.getTime()
        volume = state.getPeriodicBoxVolume()
        pe = state.getPotentialEnergy()
        ke = state.getKineticEnergy()
        if self._temperature:
            temp = 2 * ke / (self._dof * u.MOLAR_GAS_CONSTANT_R)
        if self._potentialEnergy:
            values['potentialEnergy'] = pe.value_in_unit(self.energyUnit)
        if self._kineticEnergy:
            values['kineticEnergy'] = ke.value_in_unit(self.energyUnit)
        if self._totalEnergy:
            values['totalEnergy'] = (pe + ke).value_in_unit(self.energyUnit)
        if self._temperature:
            values['temperature'] = temp.value_in_unit(self.temperatureUnit)
        if self._volume:
            values['volume'] = volume.value_in_unit(self.volumeUnit)
        if self._density:
            dens = self._totalMass / volume
            values['density'] = dens.value_in_unit(self.densityUnit)

        return values
   
    def __del__(self):
        """ We already closed the file. """

class EnergyMinimizerReporter(StateDataReporter):
    """
    This class acts as a simple energy reporter class for minimizations. This
    is not meant to be used as a reporter for OpenMM's molecular dynamics
    routines, but instead passed a simulation object for printing out
    single-point energies.
    """

    def __init__(self, f, volume=False,
                 energyUnit=u.kilocalories_per_mole,
                 volumeUnit=u.angstrom*u.angstrom*u.angstrom,
                 densityUnit=u.grams/u.item/u.milliliter):
        """
        Initializes an EnergyMinimizerReporter

        Parameters:
            - f (string or file object): File to write the energies to
            - volume (bool): Print system volume?
            - energyUnit (unit): Units to print energy in
            - volumeUnit (unit): Units to print volume in
        """
        self._volume = volume
        self.energyUnit = energyUnit
        self.volumeUnit = volumeUnit
        if hasattr(f, 'write'):
            self._openedFile = False
            self._out = f
        else:
            self._openedFile = True
            self._out = open(f, 'w')

    def describeNextReport(self, *args, **kwargs):
        """ Disable this reporter inside MD """
        raise NotImplemented('EnergyMinimizerReporter is not intended to be '
                             'used for reporting on molecular dynamics')

    def report(self, simulation, frame=None):
        """ Print out the current energy """
        has_pbc = simulation.topology.getUnitCellDimensions() is not None
        state = simulation.context.getState(getEnergy=True,
                                            enforcePeriodicBox=has_pbc)
        if frame is not None:
            self._out.write('Frame: %10d\n' % frame)
        self._out.write('   Potential Energy = %12.4f %s\n' %
                (state.getPotentialEnergy().value_in_unit(self.energyUnit),
                self.energyUnit))
        if has_pbc and (self._volume or self._density):
            vol = state.getPeriodicBoxVolume().value_in_unit(self.volumeUnit)
        if has_pbc and self._volume:
            self._out.write('   Volume = %12.4f %s\n' % (vol, self.volumeUnit))
        self._out.write('\n')

    def finalize(self):
        """ Closes any open file """
        try:
            if self._out is not None:
                self._out.close()
        except NameError:
            pass
