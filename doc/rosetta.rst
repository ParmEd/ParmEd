PyRosetta Functionality
====================

Since ParmEd is a library for reading, converting, and modifying full molecular
mechanical descriptions of chemical systems in a wide variety of different
families of force fields, supporting molecular modeling using the
`PyRosetta <http://www.pyrosetta.org/>`_ was a natural extension.

This page is not meant as a tutorial for the PyRosetta library. Instead, you
should visit the `PyRosetta website <http://www.pyrosetta.org/tutorials>`_ for
that. However, this page *will* provide a brief description of PyRosetta, and
how ParmEd can enhance its utility for molecular modeling. It will also present
a few examples using PyRosetta with the tools provided by ParmEd.

What is PyRosetta
--------------

Rosetta is a popular molecular modeling suite that allows scientists to model
biomolecular systems, and has had a lot of `success in predicting experimental
structures <http://en.wikipedia.org/wiki/Rosetta@home#Project_significance>`_.
PyRosetta is a Python-based interface to the Rosetta, and it allows users to
easily model and manipulate biomolecules, enabling custom algorithms for 
protein folding using Rosetta sampling and scoring functions.
Its features include:

    * Generating a nucleic or amino acid structure from a text string.
    * Mutating residues in an existing structure.
    * Manipulating residue positions and torsions in an algorithmic fashion.
    * Monte Carlo sampling using the Rosetta force field.


How does ParmEd enhance PyRosetta?
-------------------------------

.. currentmodule:: chemistry.rosetta
.. autosummary::
    :toctree: rosettaobj/

    RosettaPose

ParmEd provides a simple function, :func:`load_rosetta`, to load a PyRosetta
:class:`Pose` object into ParmEd's :class:`Structure` class. From this
:class:`Structure`, a user can proceed to use the :attr:`positions` and
:attr:`topology` attributes to start a molecular dynamics simulation.

One might imagine that this functionality would be useful to easily seed
simulations  for *ab initio* structure predictions or studies of a mutant
protein using wildtype structures.


Examples
--------

For examples on how to run an simulation starting from a
PyRosetta :class:`Pose` please refer to the following:

.. toctree::
    :maxdepth: 2

    Using PyRosetta to seed an OpenMM simulation <omm_rosetta>
