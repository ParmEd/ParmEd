PyRosetta Integration
====================

ParmEd provides support for the `PyRosetta <http://www.pyrosetta.org/>`_
API. The ability to create a *de novo* molecular system
and then seamlessly parameterize and simulate it, all within the Python
environment, is powerful and fits right within ParmEd's goal of easing
the molecular modeling process.

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

.. currentmodule:: parmed
.. autosummary::
    :toctree: rosettaobj/

    load_rosetta

ParmEd provides a simple function, :func:`load_rosetta`, to load a PyRosetta
:class:`Pose` object into ParmEd's :class:`Structure` class. This
:class:`Structure` can then be used to easily parameterize a molecular dynamics
simulation.

One can imagine that this functionality would be useful to easily seed
simulations  for *ab initio* structure predictions or studies of a mutant
protein using wildtype structures.


Examples
--------

For examples on how to run an simulation starting from a
PyRosetta :class:`Pose` please refer to the following:

.. toctree::
    :maxdepth: 1

    Using PyRosetta and ParmEd to seed an OpenMM simulation <omm_rosetta>
