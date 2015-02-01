"""
This module contains classes that can serve as a drop-in replacement for the
Amber file reading classes provided by OpenMM for creating an OpenMM System
for performing simulations.

It also pulls the box information from the restart file instead of the topology
file.
"""
from __future__ import division

from chemistry.constants import TINY, DEG_TO_RAD
from chemistry.amber.readparm import AmberParm, ChamberParm, Rst7
from chemistry.exceptions import OpenMMError
import chemistry.periodic_table as pt
from math import cos, sin, sqrt, pi
import simtk.openmm as mm
from simtk.openmm.vec3 import Vec3
import simtk.unit as u
from simtk.openmm.app import (forcefield as ff, Topology, element)
from simtk.openmm.app.amberprmtopfile import HCT, OBC1, OBC2, GBn, GBn2
from simtk.openmm.app.internal.customgbforces import (GBSAHCTForce,
                GBSAOBC1Force, GBSAOBC2Force, GBSAGBnForce, GBSAGBn2Force,
                convertParameters)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

OpenMMAmberParm = AmberParm
OpenMMChamberParm = ChamberParm

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class OpenMMRst7(Rst7):
    """
    Contains positions and maybe velocities and box information from a restart
    file using OpenMM data structures
    """

    @property
    def positions(self):
        """ Positions as a sequence of Vec3 objects with angstrom units """
        # Return the cached copy of positions
        try:
            return self._positions
        except AttributeError:
            pass

        # The first time it's called, cache the data structure
        self._positions = tuple([Vec3(self.coordinates[3*i],
                self.coordinates[3*i+1], self.coordinates[3*i+2])
                for i in xrange(self.natom)]) * u.angstroms

        return self._positions

    @property
    def velocities(self):
        """ Same as positions, but for velocities """
        try:
            return self._velocities
        except AttributeError:
            pass
      
        self._velocities = tuple([Vec3(self._ambvels[3*i],
            self._ambvels[3*i+1], self._ambvels[3*i+2])
            for i in xrange(self.natom)]) * u.angstroms / u.picoseconds 

        return self._velocities

    @velocities.setter
    def velocities(self, stuff):
        """
        This is a hack to work around the new Amber Rst7 class using the name
        'velocities'. This forces the Amber object to store the velocities in a
        special '_ambvels' attribute instead
        """
        self._ambvels = stuff

    @property
    def box_vectors(self):
        """ Return tuple of box vectors """
        if self.box is None: return None
        box = [x*u.angstrom for x in self.box[:3]]
        ang = [self.box[3], self.box[4], self.box[5]]
        return _box_vectors_from_lengths_angles(box[0], box[1], box[2],
                                                ang[0], ang[1], ang[2])

    @property
    def box_lengths(self):
        """ Return tuple of floats """
        box = [x*u.angstrom for x in self.box[:3]]
        return tuple(box)

def _box_vectors_from_lengths_angles(a, b, c, alpha, beta, gamma):
    """
    This method takes the lengths and angles from a unit cell and creates unit
    cell vectors.

    Parameters:
        - a (unit, dimension length): Length of the first vector
        - b (unit, dimension length): Length of the second vector
        - c (unit, dimension length): Length of the third vector
        - alpha (float): Angle between b and c in degrees
        - beta (float): Angle between a and c in degrees
        - gamma (float): Angle between a and b in degrees

    Returns:
        Tuple of box vectors (as Vec3 instances)
    """
    if not (u.is_quantity(a) and u.is_quantity(b) and u.is_quantity(c)):
        raise TypeError('a, b, and c must be units of dimension length')
    if u.is_quantity(alpha): alpha = alpha.value_in_unit(u.degree)
    if u.is_quantity(beta): beta = beta.value_in_unit(u.degree)
    if u.is_quantity(gamma): gamma = gamma.value_in_unit(u.degree)
    a = a.value_in_unit(u.angstrom)
    b = b.value_in_unit(u.angstrom)
    c = c.value_in_unit(u.angstrom)

    if alpha <= 2 * pi and beta <= 2 * pi and gamma <= 2 * pi:
        raise ValueError('box angles must be given in degrees')

    alpha *= DEG_TO_RAD
    beta *= DEG_TO_RAD
    gamma *= DEG_TO_RAD

    av = Vec3(a, 0.0, 0.0) * u.angstrom
    bx = b * cos(gamma)
    by = b * sin(gamma)
    bz = 0.0
    cx = c * cos(beta)
    cy = c * (cos(alpha) - cos(beta) * cos(gamma))
    cz = sqrt(c * c - cx * cx - cy * cy)
   
    # Make sure any components that are close to zero are set to zero exactly
    if abs(bx) < TINY: bx = 0.0
    if abs(by) < TINY: by = 0.0
    if abs(cx) < TINY: cx = 0.0
    if abs(cy) < TINY: cy = 0.0
    if abs(cz) < TINY: cz = 0.0

    bv = Vec3(bx, by, bz) * u.angstrom
    cv = Vec3(cx, cy, cz) * u.angstrom

    return (av, bv, cv)
