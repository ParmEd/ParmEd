"""
This module used to contain OpenMM-enabled components of the Amber classes.
This functionality has since been incorporated directly into the original
classes. These classes are deprecated and may be phased out in later versions of
ParmEd
"""

from chemistry.amber.readparm import (AmberParm as OpenMMAmberParm,
                                      ChamberParm as OpenMMChamberParm,
                                      Rst7 as OpenMMRst7)
