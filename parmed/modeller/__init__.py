"""
This package contains functionality necessary to carry out basic molecular
modelling tasks.
"""

__author__ = 'Jason Swails <jason.swails@gmail.com>'
__date__ = '2015'
__all__ = ['AmberOFFLibrary']

from parmed.modeller.residue import *
from parmed.amber.offlib import AmberOFFLibrary

__all__ += residue.__all__
