"""
Physical quantities with units for dimensional analysis and automatic unit
conversion.
"""
from __future__ import division

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
# `parmed.unit` package can be used interchangeably with OpenMM

try:
    from simtk.unit import *
except ImportError:
    from parmed.unit.unit import Unit, is_unit
    from parmed.unit.quantity import Quantity, is_quantity
    from parmed.unit.unit_math import *
    from parmed.unit.unit_definitions import *
    from parmed.unit.constants import *

# Now create the AKMA unit system, which is in common use by a number of
# programs (like Amber and CHARMM). Most complicated thing to get is the time.
# Do it from the kJ->kcal conversion

_time_scale = (daltons * angstroms**2 /
               kilocalories_per_mole).conversion_factor_to(picoseconds**2)
akma_time_base_unit = BaseUnit(time_dimension, "akma_time_unit", "akma-s")
akma_time_base_unit.define_conversion_factor_to(picosecond_base_unit,
            1 / sqrt(_time_scale))
del _time_scale

akma_unit_system = UnitSystem([
        angstrom_base_unit,
        dalton_base_unit,
        akma_time_base_unit,
        elementary_charge_base_unit,
        kelvin_base_unit,
        mole_base_unit,
        radian_base_unit])
