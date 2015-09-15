Units and Physical Constants in ParmEd
======================================

This page details the unit system used by ParmEd to store information
internally, as well as some relevant values of physical constants that it uses
for representing electrostatic interactions.

Units
-----

As discussed `in the section on dimensional analysis
<dimensional_analysis.html>`_, different programs use different sets of internal
units for different reasons. The default set of units that ParmEd uses to store
quantities internally is the *AKMA* unit system, with the exception that angles
are stored in degrees (rather than radians), and time quantities are stored in
picoseconds (this is also true for quantities with a simple inclusion of "time",
like velocities which are stored in Angstroms per picosecond).

If you are at all unsure, you may assign or initialize attributes or objects
using numbers that bear units, and the conversion to the internal unit set used
by ParmEd will be performed automatically. By default, reading files from
specific formats will result in the unit system associated with the parent
program being applied to the parsed information (see the table at the bottom of
`this page <dimensional_analysis.html#available-unitsystem-options>`_).

Electrostatic Constant
----------------------

The most problematic physical constant that complicates the reproducibility of
interconverted files between different programs is the value chosen for the
electrostatic scaling constant. The *exact* value used for the different
supported programs are shown below:

+--------------+----------------------------+-----------------------------------------------+
| Program      | Electrostatic constant     | Units                                         |
+==============+============================+===============================================+
| AMBER        | :math:`18.2223^2`          | :math:`kcal\ Ang/(mol\ elementary\ charge^2)` |
+--------------+----------------------------+-----------------------------------------------+
| CHARMM       | :math:`332.0716`           | :math:`kcal\ Ang/(mol\ elementary\ charge^2)` |
+--------------+----------------------------+-----------------------------------------------+
| GROMACS      | :math:`138.93545783908448` | :math:`kJ\ nm/(mol\ elementary\ charge^2)`    |
+--------------+----------------------------+-----------------------------------------------+
| OpenMM       | :math:`138.935456`         | :math:`kJ\ nm/(mol\ elementary\ charge^2)`    |
+--------------+----------------------------+-----------------------------------------------+
| ParmEd       | :math:`332.0634322112194`  | :math:`kcal\ Ang/(mol\ elementary\ charge^2)` |
+--------------+----------------------------+-----------------------------------------------+


GROMACS computes the electrostatic scaling constant as a composite calculation
from a number of other fundamental physical constants. The result shown here is
the full precision returned by the calculation of the constant inside Python.
For ParmEd, the value of the constant was computed from constants defined in
Table 1 of the 77th edition of the CRC Handbook of Chemistry and Physics.

What ParmEd does to make sure that the representation of the charges is
identical in every program supported by ParmEd, the input charges stored in the
parameter and/or topology files (e.g., the PSF, GROMACS topology, or Amber
prmtop) is multiplied by the ratio of the electrostatic constant used in the
parent program of the file format and the ParmEd internal constant. Upon
writing, the charge is multiplied by the inverse of this ratio.

The advantage of this approach is that it will work for every file reader/writer
that gets added (as long as the scaling is done appropriately), and
interconversions between supported file formats will result in different
programs using the same *effective* electrostatic constant despite the fact that
the exact values differ. This significantly improves agreement in calculated
electrostatic energies between the different programs.

The disadvantage of this approach is that the stored charges in general *no
longer match the exact charges defined in the input file*, since they are
"converted" to ParmEd's internal representation of the electrostatic scaling
constant. In practice, the difference between these electrostatic scaling
constants is small, but I do not know of any rigorous study *demonstrating* that
they are completely unimportant.

In order to be safe with regards to changing charges, ``parmed.tools.change``
will automatically treat the *input* charge as though it should be treated with
the electrostatic constant of the parent program for the current format (see the
table above).
