.. ParmEd documentation master file, created by
   sphinx-quickstart on Wed Feb  4 09:24:43 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ParmEd
======

ParmEd is a general tool for aiding in investigations of biomolecular systems
using popular molecular simulation packages, like Amber, CHARMM, and OpenMM
written in Python.

What is it?
-----------

There are two parts to ParmEd -- the Python API that exposes the core classes
used in its modeling capabilities, and two front-end Python programs
(``parmed`` and its GUI counterpart, ``xparmed``) that make use of the
ParmEd API to allow rapid prototyping and parameter-topology modifications for
use in molecular simulations.

Why use it?
-----------

You can use it for a variety of modeling purposes, like manipulating system
topologies (i.e., the *atoms*, *bonds*, and other bonded terms like valence
angles) and translating between file formats (e.g. PDB, mmCIF/PDBx, Amber
*prmtop*, and CHARMM *psf*) and even other APIs (e.g. PyRosetta).
What sets ParmEd apart from tools like OpenBabel is that it stores and tracks
force field parameters so that the resulting files can be used to carry out
molecular mechanics simulations with tools like Amber, OpenMM, NAMD,
and CHARMM.

ParmEd has sophisticated machinery built into its core classes that
substantially reduces the burden of bookkeeping when manipulating chemical
structures (for instance, adding or deleting an atom from a structure
automatically updates the indices in each of the parameter arrays defining the
system topology so you don't have to worry about it).

What can it do?
---------------

The core ``parmed`` package provided as part of ParmEd is intended to provide
a powerful, but simple, API for carrying out tasks common in the field of
biomolecular simulations.  For example, some of its features include:

    - Simply iterate through atoms as well as other parameters found in common
      molecular mechanical force fields
    - Rapidly modify a system topology and its parameters, and write files that
      can be used with the Amber and CHARMM program suites (as well as other
      programs that support those file types, like NAMD)
    - Carry out investigations using force fields, like molecular dynamics,
      directly on modern, high-performance hardware (like Graphics Processing
      Units) using the *OpenMM* library and Python application layer
    - Dimensional analysis with a complete system of physical units
    - Translate between a variety of file formats in common use by various
      programs including:

        + Amber prmtop, inpcrd, NetCDF trajectory, and NetCDF restart files
        + CHARMM PSF, coordinate, and restart files
        + Gromacs topology and GRO files
        + PDB files, supporting a wide range of dialects that technically
          violate the PDB standard
        + PDBx/mmCIF files -- the new standard for the Protein Data Bank

    - Extract metadata from the PDB and PDBx/mmCIF files, such as citation
      information and related database entries

Roadmap: Main goals and future directions
-----------------------------------------

One of the main goals of ParmEd is to provide a single interface for all of the
various biomolecular simulation programs and file formats out there and provide
a platform upon which transferring data between these different formats and
programs is easy and, most importantly, reliable.

Every program is different--developed by different people for different purposes
to solve perhaps slightly different scientific problems. As such, each program
has something to offer that the others don't, be it a fancy new method, improved
computational performance, better force field for a particular molecule, easier
simulation setup, etc. It is currently very difficult to combine components of
different program suites into a single workflow, however, given the highly
specialized nature of the various file formats for the programs they were
written around (e.g., the Amber topology file in the Amber programs, the GROMACS
top and itp files for GROMACS, etc.).

Some programs out there will allow you to take, for instance, an Amber topology
file and convert it into something GROMACS will understand, but the reverse is
not available in any tool.  ParmEd hopes to bridge this gap in addition to
providing a flexible API that will allow you to go beyond what each programs'
modeling tools allow you to do by themselves.  It is an ambitious goal, to be
sure, but good progress has already been made.

Check out the `Github repository <http://github.com/ParmEd/ParmEd>`_ and its
issue tracker to keep up-to-date with the planned as well as on-going
developments!

Slides and presentations
------------------------

I will post any slides pertaining to ParmEd from talks that I've given here, in
the hopes that they may be helpful or informative.

- `April 10, 2015 at MSKCC <http://parmed.github.io/ParmEd/ParmEd_Slides_08Apr2015.pdf>`_
- `March 19, 2016 at the Amber developer meeting at UCSD <http://parmed.github.io/ParmEd/ParmEd_Amber_Meeting2016.pdf>`_

Getting Started
---------------

When to use ParmEd?
~~~~~~~~~~~~~~~~~~~

- When you want to extract information about a structure or parameters from one
  of the supported file formats.
- When you want to manipulate molecular mechanical (force field) descriptions of
  chemical systems to quickly prototype ideas involving Hamiltonian
  modifications.
- When you want to manipulate chemical systems by selecting a subset of the
  atoms.
- When you want to combine or replicate the contents of one or more chemical
  systems.
- When you want to convert molecular mechanical descriptions of molecules from a
  format that one program supports to another (e.g., converting an Amber
  topology file to a Gromacs topology file or vice-versa).
- When you want to carry out molecular dynamics investigations on
  high-performance computational hardware (like GPUs) using OpenMM

Common examples
~~~~~~~~~~~~~~~

The simplest example is to use ParmEd to download and inspect a PDB file. We can
find out numerous attributes about a structure defined by a PDB file, like the
number of atoms, residues, and even various citation information. For example::

    >>> import parmed as pmd
    >>> pmd.download_PDB('4lzt')
    <Structure 1164 atoms; 274 residues; 1043 bonds; PBC (triclinic); NOT parametrized>
    >>> lysozyme = pmd.download_PDB('4lzt')
    >>> lysozyme.authors
    'M.A.WALSH,T.SCHNEIDER,L.C.SIEKER,Z.DAUTER,V.LAMZIN,K.S.WILSON'
    >>> lysozyme.experimental
    'X-RAY DIFFRACTION'
    >>> len(lysozyme.atoms)
    1164
    >>> len(lysozyme.residues)
    274
    >>> # Get all backbone atoms
    ... lysozyme['@CA,C,N']
    <Structure 393 atoms; 135 residues; 386 bonds; NOT parametrized>

A simple example demonstrating the file conversion capabilities is to convert a
PDBx/mmCIF file into the more commonly supported PDB format::

    >>> import parmed as pmd
    >>> pmd.download_CIF('4lzt').save('4lzt.pdb')
    >>> # Now read in the PDB file we just created
    ... pmd.load_file('4lzt.pdb')
    <Structure 1164 atoms; 274 residues; 1043 bonds; PBC; NOT parametrized>

Program and API Reference
-------------------------

.. toctree::
   :maxdepth: 1

   Core classes used to represent topologies <topologyobjects>
   The core Structure class <structure>
   Working with units <dimensional_analysis>
   Reading and writing files <readwrite>
   The Amber file classes <amber>
   The CHARMM file classes <charmm>
   The GROMACS file classes <gromacs>
   PyRosetta Integration <rosetta>
   Using the ParmEd programs (``parmed`` and ``xparmed``) <parmed>
   The ParmEd API <parmed_api>
   OpenMM Functionality <openmm>
   Visualization Functionality <visualization>

Click the "modules" link at the top of the page for a full API reference.

Contributing
------------

If you would like to contribute, please head over to the Github repository at
https://github.com/ParmEd/ParmEd and feel free to jump right in! New issues
and pull requests are always welcome.

If you want to start contributing to the code base, please check out `the
developer's documentation <devdoc>`_ for general information about style guides
and the like.

Search
------

* :ref:`search`
