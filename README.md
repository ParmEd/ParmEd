ParmEd
======

Cross-program parameter and topology file editor and molecular mechanical
simulator engine.

Build Status
============

[![Linux Build Status](https://travis-ci.org/ParmEd/ParmEd.svg?branch=master)](https://travis-ci.org/ParmEd/ParmEd)

Description
===========

ParmEd is a package designed to facilitate creating and easily manipulating
molecular systems that are fully described by a common classical force field.
Supported force fields include Amber, CHARMM, AMOEBA, and several others that
share a similar functional form (e.g., GROMOS).

ParmEd is capable of reading and writing to a wide array of different file
formats, like the Amber topology and coordinate files, CHARMM PSF, parameter,
topology, and coordinate files, Tinker parameter, topology, and coordinate
files, and many others. The expressive central data structure (the ``Structure``
class) makes it easy to quickly and safely manipulate a chemical system, its
underlying topology, and force field parameters describing its potential energy
function.

There are two parts of ParmEd---a documented API that you can incorporate into
your own Python scripts and programs, and a GUI/CLI pair of programs that
provide a means to quickly perform various modifications to chemical systems for
rapid prototyping.

The API also provides bindings to the [OpenMM](https://simtk.org/home/openmm)
library, permitting you to carry out full molecular dynamics investigations
using ParmEd on high-performant hardware, like AMD and NVidia GPUs.

Documentation
=============

Want to learn more?  Visit the ParmEd documentation page at
https://parmed.github.io/ParmEd for examples, descriptions, and API
documentation.

Authors and Contributors
========================

The following people have contributed directly to the coding and validation
efforts in ParmEd (in alphabetical order).  And a special thanks to all of you
who helped improve this project either by providing feedback, bug reports, or
other general comments!

* Jason Swails (principal developer) | jason.swails@gmail.com

* Carlos Hernandez
* David L. Mobley
* Lee-Ping Wang
* Pawel Janowski

License
=======

```
                    LESSER GPL LICENSE INFO

Copyright (C) 2010 - 2014 Jason Swails

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU Lesser General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.

The `fortranformat` package is released under the MIT license, Copyright (C)
Brendan Arnold. See `fortranformat/__init__.py` for more information.
```
