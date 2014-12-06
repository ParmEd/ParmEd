"""
Physical quantities with units for dimensional analysis and automatic unit
conversion.
"""
__docformat__ = "epytext en"

__author__ = "Christopher M. Bruns"
__copyright__ = "Copyright 2010, Stanford University and Christopher M. Bruns"
__credits__ = []
__license__ = "MIT"
__maintainer__ = "Christopher M. Bruns"
__email__ = "cmbruns@stanford.edu"

# This code is copied from the `simtk.unit` package distributed as a standalone
# package (https://pypi.python.org/pypi/simtk.unit/) and as part of OpenMM
# (https://simtk.org/home/openmm).

# When OpenMM can be imported, the unit package will be taken from there.
# Otherwise, the implementation here will be used. This way, the
# `chemistry.unit` package can be used interchangeably with OpenMM

try:
    import nomod
    from simtk.unit import *
except ImportError:
    from chemistry.unit.unit import Unit, is_unit
    from chemistry.unit.quantity import Quantity, is_quantity
    from chemistry.unit.unit_math import *
    from chemistry.unit.unit_definitions import *
    from chemistry.unit.constants import *
