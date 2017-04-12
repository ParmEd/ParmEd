Dimensional Analysis (working with :mod:`parmed.unit`)
=========================================================

The :mod:`unit` package was originally developed by Christopher M. Bruns at
SimTK (Stanford University) for the ``simtk`` Python package in OpenMM_.

.. _OpenMM: http://simtk.org/home/openmm

Different programs use different systems of units internally while performing
calculations on quantities with specific dimensions. For instance, the *Amber*
and *CHARMM* programs deal with energies in *kcal/mol* and lengths in
*Angstroms*, while *GROMACS* and *OpenMM* deal with energies in *kJ/mol* and
distances in *nanometers*. Electronic structure packages, on the other hand,
tend to use so-called *atomic* units, with energies in Hartrees and lengths in
Bohrs.

Comparing the results between these programs requires that the results be
converted to a common set of units beforehand. More
importantly---interconverting these packages, or preparing common input for
both, requires that the input for each program conform to its expected system of
units. Not doing so can result in `a very costly mistake
<http://www.wired.com/2010/11/1110mars-climate-observer-report/>`_

In all of the examples below, I will assume that the :mod:`unit` package has
been imported as follows::

    from parmed import unit as u

What does :mod:`parmed.unit` do?
-----------------------------------

The main classes in the :mod:`unit` module are shown below:

.. currentmodule:: parmed.unit
.. autosummary::
    :toctree: unitobj/
    
    Quantity
    is_quantity
    Unit
    BaseUnit
    UnitSystem

The :mod:`unit` package allows you to turn *numbers* into a :class:`Quantity`,
that has both a value and a unit.  It allows you to do basic math with units and
it automatically performs the requisite dimensional analysis, raising an
exception if you try to illegally combine incompatible units.  For example::

    >>> x = 1.0 * u.nanometers
    >>> x + 1.0 * u.angstroms
    Quantity(value=1.1, unit=nanometer)
    >>> # Find the volume of a cube with side 1 cm
    >>> x = 1 * u.centimeters
    >>> v = x ** 3
    >>> v
    Quantity(value=1, unit=centimeter**3)
    >>> v + 1.0 * u.milliliters
    Quantity(value=1.9999999999999998, unit=centimeter**3)
    >>> # Now try to combine incompatible units
    >>> 1*u.millimeters + 1*u.kelvin
    TypeError: Cannot add two quantities with incompatible units "millimeter" and "kelvin".

Two units can be multiplied together to form a composite unit with the
corresponding dimensionality---in the example above, cubing 1 centimeter
resulted in a volume, which is a cubic length dimension. This works for much
more complex units as well::

    >>> x = 1 * u.kilograms * u.meters**2 / u.seconds**2
    >>> x
    Quantity(value=1, unit=kilogram*meter**2/(second**2))
    >>> x + 1 * u.joule
    Quantity(value=2, unit=kilogram*meter**2/(second**2))

How do I use ``units``?
-----------------------

There are a variety of tasks you may want to do with units -- the more popular
of which is shown in the sections below.

Checking to see if an object is a :class:`Quantity`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There may be times when you want to see if a Python object is a
:class:`Quantity` instance, or if it is just a raw number. Many functions in
ParmEd do this when taking user input, converting to its internal AKMA unit
system (see below). Checking that an object is a :class:`Quantity` is done using
the :func:`is_quantity` function in the ``parmed.unit`` namespace::

    >>> x = 1
    >>> y = 1 * u.nanometers
    >>> u.is_quantity(x)
    False
    >>> u.is_quantity(y)
    True

Creating units
~~~~~~~~~~~~~~

You can create a unit in one of two ways:

    1. Use the :class:`Quantity` constructor, as shown below::

           >>> u.Quantity(1, u.nanometers)
           Quantity(1, nanometer)
           >>> u.Quantity(1, u.meters**-3) # You can raise units to exponents
           Quantity(value=1, unit=/(meter**3))

    2. Multiply (or divide!) by a specific unit, as shown below::

           >>> 1 * u.nanometers
           Quantity(1, nanometer)
           >>> 1 / u.meters**3 # You can raise units to exponents
           Quantity(value=1, unit=/(meter**3))

Seems easy enough, and I usually use method 2 when turning a ``float``, ``int``,
``list``, or ``tuple`` into a :class:`Quantity`.

Converting to another unit
~~~~~~~~~~~~~~~~~~~~~~~~~~

Aside from using units as a type-checker to make sure you are not making an
arithmetic mistake with your transformations, the most common reason to use the
:mod:`unit` package is to handle unit conversions for you. The :class:`Quantity`
class defines two functions for this:

    1. ``in_units_of(unit)``: This method returns a new :class:`Quantity`
       instance with the specified ``unit`` (a ``TypeError`` is raised).
    2. ``value_in_unit(unit)``: This method returns a copy of the underlying
       data *without* an attached ``unit`` (i.e., it is *not* a
       :class:`Quantity` instance). This should be used when a scalar is needed
       from a quantity for use in a particular API, or if the performance hit
       from the dimensional analysis is unacceptable.

The following example shows both of the above in action, demonstrating their
differences::

    >>> x = 1 * u.kilocalorie
    >>> x.in_units_of(u.kilojoule)
    Quantity(value=4.184, unit=kilojoule)
    >>> x.value_in_unit(u.kilojoule)
    4.184

Converting units in batch -- the :class:`UnitSystem`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The previous section described how to convert a :class:`Quantity` into another
:class:`Quantity` in a different, but *compatible* unit. The argument to
``in_units_of`` and ``value_in_unit`` methods *must* be a compatible unit, which
can be limiting if you want to write a function that changes the units of an
argument with unknown dimensionality into a particular system of units.

This is where :class:`UnitSystem` comes in. Each :class:`Quantity` instance also
has a corresponding methods for conversion to units compatible within a
particular :class:`UnitSystem`:

    1. ``in_unit_system(UnitSystem)``: This method returns a new
       :class:`Quantity` instance composed of the base units defined in that
       unit system.
    2. ``value_in_unit_system(UnitSystem)``: This method returns a scalar of the
       value in that particular unit system.

For example::

    >>> x = 1 * u.kilocalorie
    >>> x.in_unit_system(u.md_unit_system)
    Quantity(value=4.184, unit=nanometer**2*mole*dalton/(picosecond**2))
    >>> x.value_in_unit_system(u.md_unit_system)
    4.184

This can be useful when designing a function that needs to handle arbitrary
units, but it comes with some drawbacks (see below), so you should avoid unit
system conversions when you can.

A warning about :class:`UnitSystem`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note that this method does not provide the added security of checking that the
*dimensions* are correct -- it simply converts all base dimensions to the ones
defining the unit system. For instance, consider the two functions below::

    def area(x, y):
        """ Computes the area of a rectangle

        Parameters
        ----------
        x : distance Quantity
            The length of the first side
        y : distance Quantity
            The length of the second side

        Returns
        -------
        area : float
            Area in nanometers^2
        """
        return (x.value_in_unit_system(u.md_unit_system) *
                    y.value_in_unit_system(u.md_unit_system))

    def area2(x, y):
        return (x.value_in_unit(u.nanometers) *
                    y.value_in_unit(u.nanometers))

As long as you provide the correct dimensions in the arguments, both functions
yield the same result::

    >>> area(1*u.angstroms, 1*u.millimeters)
    100000.0
    >>> area2(1*u.angstroms, 1*u.millimeters)
    100000.0

Consider what happens if you make a small mistake, though::

    >>> area(1/u.angstroms, 1*u.millimeters)
    9999999.999999998
    >>> area2(1/u.angstroms, 1*u.millimeters)
    TypeError: Unit "/angstrom" is not compatible with Unit "nanometer".

Converting ``1/angstrom`` to the ``md_unit_system`` converts the value to
``10/nanometer`` (rather than ``1 angstrom`` to ``0.1 nanometer``), leading to a
result that is 2 orders of magnitude larger than expected! The ``area2``
function catches the error and raises an exception, though, thereby making the
mistake easier to find.

Available :class:`UnitSystem` options
-------------------------------------

The :mod:`unit` module contains a number of :class:`UnitSystem` instances
available in the ``parmed.unit`` namespace. They are summarized in the table
below, along with the units defining them (the energy unit is a composite of the
mass unit times the square of the ratio of the length-to-time units):

+------------------------+-------------+-----------+---------------+-----------------+------------------+-------------+-------------+
| :class:`UnitSystem`    | Length unit | Mass unit | Time unit     | Charge unit [1] | Temperature unit | Amount unit | Energy unit |
+========================+=============+===========+===============+=================+==================+=============+=============+
| ``si_unit_system``     |   meters    | kilograms |  second       |    Ampere       |      kelvin      |    mole     |    joule    |
+------------------------+-------------+-----------+---------------+-----------------+------------------+-------------+-------------+
| ``cgs_unit_system``    | centimeters |   gram    |  second       |    Ampere       |      kelvin      |    mole     | 1e-7 joule  |
+------------------------+-------------+-----------+---------------+-----------------+------------------+-------------+-------------+
| ``md_unit_system``     | nanometers  |  daltons  | picosecond    |  q electron     |      kelvin      |    mole     |  kilojoule  |
+------------------------+-------------+-----------+---------------+-----------------+------------------+-------------+-------------+
| ``planck_unit_system`` | pl. length  | pl. mass  | pl. time      | pl. charge      | pl. temperature  |    item     | pl. energy  |
+------------------------+-------------+-----------+---------------+-----------------+------------------+-------------+-------------+
| ``akma_unit_system``   | angstroms   |  daltons  | akma time [2] |  q electron     |      kelvin      |    mole     | kilocalorie |
+------------------------+-------------+-----------+---------------+-----------------+------------------+-------------+-------------+

[1] Some unit systems use current instead of charge as a base unit, while others
use the charge.

[2] The time in the AKMA unit system (Angstroms, Kilocalories per Mole, Atomic
mass units) does not have an official name, and is roughly equal to 20.455
picoseconds.

Different programs tend to use different unit systems. As a general reference,
the table below summarizes a subset of programs and the unit systems they use:

+---------+----------------------+
| Program |      Unit System     |
+=========+======================+
| AMBER   | ``akma_unit_system`` |
+---------+----------------------+
| CHARMM  | ``akma_unit_system`` |
+---------+----------------------+
| Tinker  | ``akma_unit_system`` |
+---------+----------------------+
| Desmond | ``akma_unit_system`` |
+---------+----------------------+
| LAMMPS  | ``akma_unit_system`` |
+---------+----------------------+
|  NAMD   | ``akma_unit_system`` |
+---------+----------------------+
|  AceMD  | ``akma_unit_system`` |
+---------+----------------------+
| OpenMM  | ``md_unit_system``   |
+---------+----------------------+
| Gromacs | ``md_unit_system``   |
+---------+----------------------+
| Gromos  | ``md_unit_system``   |
+---------+----------------------+
