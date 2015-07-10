"""
This module is simply a namespace for all Amber topology-like classes. There are
different variants; for instance the standard Amber topology, the chamber-made
topology, and the tinker_to_amber-made topology. These classes are all defined
in their own private modules, but imported here to simplify the API.

Copyright (C) 2010 - 2014  Jason Swails

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""

from parmed.amber.amberformat import AmberFormat
from parmed.amber._amberparm import AmberParm, Rst7
from parmed.amber._chamberparm import ChamberParm, ConvertFromPSF
from parmed.amber._tinkerparm import AmoebaParm, BeemanRestart
from warnings import warn as _warn

# Define importables via *
__all__ = ['AmberFormat', 'AmberParm', 'ChamberParm', 'LoadParm', 'Rst7']

# Supply a function to load a topology file in the 'correct' format
def LoadParm(parmname, xyz=None, box=None, rst7name=None):
    """
    Loads a topology file using the correct class.

    Parameters
    ----------
    parmname : ``str``
        The name of the topology file to load
    xyz : str or array, optional
        If provided, the coordinates and unit cell dimensions from the provided
        Amber inpcrd/restart file will be loaded into the molecule, or the
        coordinates will be loaded from the coordinate array
    box : array, optional
        If provided, the unit cell information will be set from the provided
        unit cell dimensions (a, b, c, alpha, beta, and gamma, respectively)

    Returns
    -------
    parm : :class:`AmberParm` (or subclass)
        This function parses the topology file, determines if it is an
        Amber-style (i.e., *traditional* Amber force field), Chamber-style
        (i.e., CHARMM force field), or Amoeba-style (i.e., Amoeba force field),
        and then returns an instance of the appropriate type.
    """
    parm = AmberFormat(parmname)
    if 'CTITLE' in parm.flag_list:
        parm = parm.view(ChamberParm)
    elif 'AMOEBA_FORCEFIELD' in parm.flag_list:
        parm = parm.view(AmoebaParm)
    else:
        parm = parm.view(AmberParm)

    # Now read the coordinate file if applicable
    if xyz is None and rst7_name is not None:
        warn('rst7_name keyword is deprecated. Use xyz instead',
             DeprecationWarning)
        xyz = rst7_name
    elif xyz is not None and rst7_name is not None:
        warn('rst7_name keyword is deprecated and ignored in favor of xyz',
             DeprecationWarning)

    if isinstance(xyz, string_types):
        f = load_file(xyz)
        if not hasattr(f, 'coordinates') or f.coordinates is None:
            raise TypeError('%s does not have coordinates' % xyz)
        self.coordinates = f.coordinates
        if hasattr(f, 'box') and f.box is not None and box is None:
            self.box = f.box
    else:
        self.coordinates = xyz
    if box is not None:
        self.box = box

    # If all else fails, set the box from the prmtop file
    if self.parm_data['POINTERS'][IFBOX] > 0 and self.box is None:
        box = self.parm_data['BOX_DIMENSIONS']
        self.box = list(box[1:]) + [box[0], box[0], box[0]]

    self.hasbox = self.box is not None

    return parm
