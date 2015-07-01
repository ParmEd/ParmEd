"""
This module contains the functionality for carrying out geometrical calculations
for molecules and molecular systems

Author: Jason Swails
Contributors:

Copyright (C) 2014 - 2015  Jason Swails

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
59 Temple Place - Suite 330
Boston, MA 02111-1307, USA.
"""
from __future__ import division

from parmed import unit as u
from parmed.constants import TINY, DEG_TO_RAD, RAD_TO_DEG
from parmed.vec3 import Vec3
from math import pi, cos, sin, sqrt, acos
import numpy as np
import warnings

def box_lengths_and_angles_to_vectors(a, b, c, alpha, beta, gamma):
    """
    This function takes the lengths of the unit cell vectors and the angles
    between them and returns 3 unit cell vectors satisfying those dimensions

    Parameters
    ----------
    a : double (or length Quantity)
        Length of the first unit cell vector
    b : double (or length Quantity)
        Length of the second unit cell vector
    c : double (or length Quantity)
        Length of the third unit cell vector
    alpha : double (or angle Quantity)
        Angle between vectors b and c
    beta : double (or angle Quantity)
        Angle between vectors a and c
    gamma : double (or angle Quantity)
        Angle between vectors a and b

    Returns
    -------
    list Quantity, list Quantity, list Quantity
        The 3, 3-element vectors as quantities with dimension length

    Notes
    -----
    The unit cell lengths are assumed to be Angstroms if no explicit unit is
    given. The angles are assumed to be degrees
    """
    if u.is_quantity(a): a = a.value_in_unit(u.angstroms)
    if u.is_quantity(b): b = b.value_in_unit(u.angstroms)
    if u.is_quantity(c): c = c.value_in_unit(u.angstroms)
    if u.is_quantity(alpha): alpha = alpha.value_in_unit(u.degrees)
    if u.is_quantity(beta): beta = beta.value_in_unit(u.degrees)
    if u.is_quantity(gamma): gamma = gamma.value_in_unit(u.degrees)

    if alpha <= 2*pi and beta <= 2*pi and gamma <= 2*pi:
        warnings.warn('Strange unit cell vector angles detected. They '
                      'appear to be in radians...')

    alpha *= DEG_TO_RAD
    beta *= DEG_TO_RAD
    gamma *= DEG_TO_RAD

    av = [a, 0.0, 0.0]
    bx = b * cos(gamma)
    by = b * sin(gamma)
    bv = [bx, by, 0.0]
    cx = c * cos(beta)
    cy = c * (cos(alpha) - cos(beta)*cos(gamma)) / sin(gamma)
    cz = sqrt(c*c - cx*cx - cy*cy)
    cv = [cx, cy, cz]

    # Make sure that any tiny components are exactly 0
    if abs(bx) < TINY: bv[0] = 0
    if abs(by) < TINY: bv[1] = 0
    if abs(cx) < TINY: cv[0] = 0
    if abs(cy) < TINY: cv[1] = 0
    if abs(cz) < TINY: cv[2] = 0

    return (av, bv, cv) * u.angstroms

def box_vectors_to_lengths_and_angles(a, b, c):
    """
    This function takes the lengths of the unit cell vectors and the angles
    between them and returns 3 unit cell vectors satisfying those dimensions

    Parameters
    ----------
    a : collection of 3 floats (or length Quantity)
        The first unit cell vector
    b : collection of 3 floats (or length Quantity)
        The second unit cell vector
    c : collection of 3 floats (or length Quantity)
        The third unit cell vector

    Returns
    -------
    (a, b, c), (alpha, beta, gamma)
        Two tuples, the first is the 3 unit cell vector lengths as
        length-dimension Quantity objects and the second is the set of angles
        between the unit cell vectors as angle-dimension Quantity objects

    Notes
    -----
    The unit cell lengths are assumed to be Angstroms if no explicit unit is
    given.
    """
    if u.is_quantity(a): a = a.value_in_unit(u.angstroms)
    if u.is_quantity(b): b = b.value_in_unit(u.angstroms)
    if u.is_quantity(c): c = c.value_in_unit(u.angstroms)
    # Get the lengths
    la = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
    lb = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2])
    lc = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2])
    # Angles
    alpha = acos((b[0]*c[0] + b[1]*c[1] + b[2]*c[2]) / (lb*lc))
    beta = acos((a[0]*c[0] + a[1]*c[1] + a[2]*c[2]) / (la*lc))
    gamma = acos((b[0]*a[0] + b[1]*a[1] + b[2]*a[2]) / (lb*la))
    # Convert to degrees
    alpha *= RAD_TO_DEG
    beta *= RAD_TO_DEG
    gamma *= RAD_TO_DEG

    return (la, lb, lc) * u.angstroms, (alpha, beta, gamma) * u.degrees

def reduce_box_vectors(a, b, c):
    """
    This function puts three unit cell vectors in a reduced form where a is
    "mostly" in x, b is "mostly" in y, and c is "mostly" in z. This form is
    necessary for some programs (notably OpenMM and Gromacs)

    Parameters
    ----------
    a : 3-element collection of float
        First unit cell vector
    b : 3-element collection of float
        Second unit cell vector
    c : 3-element collection of float
        Third unit cell vector

    Returns
    -------
    red_a, red_b, red_c : Vec3, Vec3, Vec3
        The reduced unit cell vectors in units of angstroms

    Notes
    -----
    The implementation here is taken from the OpenMM Python application layer
    written by Peter Eastman
    """
    if u.is_quantity(a):
        a = a.value_in_unit(u.angstroms)
    if u.is_quantity(b):
        b = b.value_in_unit(u.angstroms)
    if u.is_quantity(c):
        c = c.value_in_unit(u.angstroms)

    a = Vec3(*a)
    b = Vec3(*b)
    c = Vec3(*c)

    c = c - b*round(c[1]/b[1])
    c = c - a*round(c[0]/a[0])
    b = b - a*round(b[0]/a[0])

    return a, b, c


def center_of_mass(coordinates, masses):
    """ Compute the center of mass of a group of coordinates.

    Parameters
    ----------
    coordinates : numpy.ndarray
        Coordinate array
    masses : numpy.ndarray
        Array of masses

    Returns
    -------
    COM
        np.ndarray of shape (3,) identifying the cartesian center of mass

    Notes
    -----
    This method *requires* that the parameters be passed in as numpy arrays.
    AttributeError's will ensue if this is not the case. Also, coordinates must
    be able to be reshaped to (len(masses), 3), or ValueError's will ensue
    """
    masses = masses.flatten()
    coordinates = coordinates.reshape((masses.shape[0], 3))
    return np.average(coordinates, weights=masses, axis=0)
