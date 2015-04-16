"""
This package contains code for reading, writing, and modifying NAMD files 
needed to start and/or restart simulations with NAMD.  The most common formats
are the so-called "namdbin" files heavily utilized by VMD; these are often
given the extensions .coor and .vel (coordinates and velocities, respectively).
"""

__all__ = ['namdbinfiles']
__authors__ = 'Brian Radak'
__contributors__ = ''
__license__ = 'GPL v.3'


from chemistry.namd.namdbinfiles import NamdBinCoor, NamdBinVel
