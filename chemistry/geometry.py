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
from chemistry import unit as u
from chemistry.constants import TINY, DEG_TO_RAD, RAD_TO_DEG
import math
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

    if alpha <= 2*math.pi and beta <= 2*math.pi and gamma <= 2*math.pi:
        warnings.warn('Strange unit cell vector angles detected. They '
                      'appear to be in radians...')

    alpha *= DEG_TO_RAD
    beta *= DEG_TO_RAD
    gamma *= DEG_TO_RAD

    av = [a, 0.0, 0.0]
    bx = b * math.cos(gamma)
    by = b * math.sin(gamma)
    bv = [bx, by, 0.0]
    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta)*math.cos(gamma))
    cz = math.sqrt(c*c - cx*cx - cy*cy)
    cv = [cx, cy, cz]

    # Make sure that any tiny components are exactly 0
    if abs(bx) < TINY: bv[0] = 0
    if abs(by) < TINY: bv[1] = 0
    if abs(cx) < TINY: cv[0] = 0
    if abs(cy) < TINY: cv[1] = 0
    if abs(cz) < TINY: cv[2] = 0

    return av*u.angstroms, bv*u.angstroms, cv*u.angstroms

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
    la = math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
    lb = math.sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2])
    lc = math.sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2])
    # Angles
    alpha = math.acos((b[0]*c[0] + b[1]*c[1] + b[2]*c[2]) / (lb*lc))
    beta = math.acos((a[0]*c[0] + a[1]*c[1] + a[2]*c[2]) / (la*lc))
    gamma = math.acos((b[0]*a[0] + b[1]*a[1] + b[2]*a[2]) / (lb*la))

    return (la, lb, lc) * u.angstroms, (alpha, beta, gamma) * u.degrees
