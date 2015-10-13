"""
XML file parsing utilities for OpenMM serialized objects
"""
from __future__ import division, print_function, absolute_import

from contextlib import closing
import numpy as np
from parmed.formats.registry import FileFormatType
from parmed.geometry import box_vectors_to_lengths_and_angles
import parmed.unit as u
from parmed.utils.decorators import needs_openmm
from parmed.utils.io import genopen
from parmed.utils.six import add_metaclass, string_types
from parmed.utils.six.moves import StringIO
try:
    from simtk import openmm as mm
    from simtk.openmm import app
except ImportError:
    mm = None
import re

_xmlre = re.compile('<(.*?)/?>')

@add_metaclass(FileFormatType)
class XmlFile(object):
    """
    Wrapper for parsing OpenMM-serialized objects. Supports serialized State,
    System, Integrator, and ForceField objects.
    """

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as an XML file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is an XML format, False otherwise
        """
        with closing(genopen(filename, 'r')) as f:
            for line in f:
                line = line.strip()
                if not line: continue
                rematch = _xmlre.match(line)
                if not rematch:
                    return False
                stuff = rematch.groups()[0]
                if stuff[0] in '?!':
                    continue
                kind = stuff.split()[0]
                if kind in ('System', 'State', 'ForceField', 'Integrator'):
                    return True
                return False

    @staticmethod
    @needs_openmm
    def parse(filename):
        """
        Parses XML file and returns deserialized object. The return value
        depends on the serialized object, summarized below

            - System : returns simtk.openmm.System
            - State : returns simtk.openmm.State
            - Integrator : returns simtk.openmm.Integrator subclass
            - ForceField : returns simtk.openmm.app.ForceField

        Parameters
        ----------
        filename : str or file-like
            The file name or file object containing the XML-serialized object

        Returns
        -------
        obj : System, State, Integrator, or ForceField
            The deserialized object

        Notes
        -----
        OpenMM requires the entire contents of this file read into memory. As a
        result, this function may require a significant amount of memory.
        """
        if isinstance(filename, string_types):
            with closing(genopen(filename, 'r')) as f:
                contents = f.read()
        else:
            contents = f.read()
        # ForceField is not handled by XmlSerializer
        if '<ForceField' in contents:
            obj = StringIO()
            obj.write(contents)
            obj.seek(0)
            return app.ForceField(obj)

        obj = mm.XmlSerializer.deserialize(contents)
        if isinstance(obj, (mm.System, mm.Integrator)):
            return obj
        elif isinstance(obj, mm.State):
            return _OpenMMStateContents(obj)
        return 

class _OpenMMStateContents(object):
    """
    A container that holds all of the information present in the OpenMM State
    object that gets deserialized. This is an internal class that is only
    intended to be *instantiated* by XmlFile.parse.

    Parameters
    ----------
    state : :class:`simtk.openmm.State`
        OpenMM State containing relevant state information

    Attributes
    ----------
    coordinates : np.ndarray with shape(1, natom, 3)
        Atomic coordinates. If not present in the State, just set to None. Units
        are Angstrom
    positions : list of Vec3
        Atomic coordinates of first frame with applied units
    velocities : np.ndarray with shape(1, natom, 3)
        Atomic velocities. If not present in the State, just set to None. Units
        are Angstrom/picosecond
    forces : np.ndarray with shape(1, natom, 3)
        Atomic forces. If not present in the State, just set to None. Units are
        kcal/mol/Angstrom
    energy : float
        Energy of this configuration in kcal/mol. If not specified, it is set to
        None
    time : float
        The time associated with this State in ps. If not specified, set to None
    """

    @staticmethod
    def _get_data(state, getter, unit, shape=None, **kwargs):
        try:
            stuff = getattr(state, getter)(**kwargs)
        except TypeError:
            return None
        else:
            stuff = stuff.value_in_unit(unit)
            if shape is not None:
                stuff = stuff.reshape(shape)
            return stuff

    def __init__(self, state):
        self.coordinates = self._get_data(state, 'getPositions',
                                u.angstrom, (1, -1, 3), asNumpy=True)
        self.velocities = self._get_data(state, 'getVelocities',
                                u.angstrom/u.picosecond, (1, -1, 3),
                                asNumpy=True)
        self.forces = self._get_data(state, 'getForces',
                                u.kilocalorie_per_mole/u.angstrom,
                                (1, -1, 3), asNumpy=True)
        self.energy = self._get_data(state, 'getPotentialEnergy',
                                u.kilocalorie_per_mole)
        self.time = self._get_data(state, 'getTime', u.picosecond)
        box = self._get_data(state, 'getPeriodicBoxVectors', u.angstroms,
                             (3, 3), asNumpy=True)
        if box is not None:
            leng, ang = box_vectors_to_lengths_and_angles(*box)
            leng = leng.value_in_unit(u.angstrom)
            ang = ang.value_in_unit(u.degree)
            self.box = np.array(list(leng) + list(ang))
        else:
            self.box = None

    @property
    def positions(self):
        return [mm.Vec3(*xyz) for xyz in self.coordinates[0]]
