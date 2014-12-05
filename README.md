ParmEd
======

Amber parameter topology file editor

Build Status
============

[![Linux Build Status](https://travis-ci.org/swails/ParmEd.svg?branch=master)](https://travis-ci.org/swails/ParmEd)

License
=======
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

Description
===========

This package is designed to make it safe and easy(er) to rapidly prototype
changes to the underlying Hamiltonian in Amber as defined in the
parameter-topology file. A full description of its capabilities is described in
in the Amber 14 manual (available at http://ambermd.org/doc12/Amber14.pdf).
The program is comprised of a command interpreter (with readline support), a
GUI-frontend based on Tkinter with clickable buttons for each action, and a
flexible Python API whose use follows naturally from the simple command
interpreter syntax.

The underlying AmberParm class that ParmEd was built around was augmented with
OpenMM support so that it is able to directly create an OpenMM System instance
with a call signature almost identical to the one provided in the OpenMM
application layer. Improvements implemented here include:

- Periodic box information is set up from the input coordinate file rather than the parameter-topology file, since that is where the 'preferred' box information is stored in Amber

- The parameter topology file class can read in both a topology file and a coordinate file at the same time (this is required to set the default box information from the inpcrd file)

- The restart file class/parser is capable of reading both Amber restarts as well as Amber NetCDF restart files, with the format detected automatically

- In many cases, System creation is noticeably faster with this class for large systems.

Also included are a handful of OpenMM reporters in the
chemistry.amber.openmmreporters module, including:

- `AmberStateDataReporter` - Prints state data in standard Amber units (and provides options in the call signature to specify the set of units you want to use)
  
- `RestartReporter` - Allows Amber-style ASCII or NetCDF restart files to be written periodically during a simulation. It can either write a series of numbered restarts or continually overwrite the previous one. Fully compatible with all Amber programs.
  
- `MdcrdReporter` - Allows Amber-style ASCII trajectory files to be written
  
- `NetCDFReporter` - Allows Amber-style NetCDF trajectories to be written

The principle difference between the NetCDF support included here and that
included in MDTraj is that the NetCDF parsing here is backend agnostic,
supporting the implementations in the scipy, netCDF4, and ScientificPython
packages (the first one is written in pure Python and does not require linking
to the NetCDF libraries). A default implementation is chosen if no explicit
choice is made in the order listed above.

Notable OpenMM Capabilities
===========================

The OpenMM action in ParmEd will read a pmemd/sander-like command-line and an
Amber-style input file and run the simulation using OpenMM through the Python
application layer with equivalent settings.  It includes options to specify the
platform and precision model, and also has the option to generate a working
Python script that will run the equivalent calculation.

The "energy" action in ParmEd will compute a single-point energy from the loaded
coordinates, optionally decomposing the energy contributions into "bond",
"angle", "dihedral", and "nonbonded" contributions (finer-grained decomposition
of the nonbonded components is not currently possible except for
direct/reciprocal decomposition).
