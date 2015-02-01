"""
This module used to contain OpenMM-enabled components of the Amber classes.
This functionality has since been incorporated directly into the original
classes. These classes are deprecated and may be phased out in later versions of
ParmEd


It also pulls the box information from the restart file instead of the topology
file.
"""

from chemistry.amber.readparm import AmberParm, ChamberParm, Rst7

OpenMMAmberParm = AmberParm
OpenMMChamberParm = ChamberParm
OpenMMRst7 = Rst7
