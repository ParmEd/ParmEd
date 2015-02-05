The :mod:`structure` module
===========================

The :mod:`structure` module contains the core class defining systems of
particles -- the :class:`Structure` class. The API documentation makes this
class look intimidating, but its core features are the ``atoms`` attribute
(which is an :class:`AtomList`) and its ``residues`` attribute (which is a
:class:`ResidueList`). These two attributes are enough to write numerous files
(e.g., PDB and PDBx/mmCIF), and these are the only two attributes currently
populated through the PDB and mmCIF parsing routines.

------

:class:`Structure` class
------------------------

.. autoclass:: chemistry.structure.Structure
    :members: add_atom, add_atom_to_residue, copy, is_changed, unchange,
              prune_empty_terms, update_dihedral_exclusions, strip, write_pdb,
              write_cif, topology, createSystem

------

:func:`read_PDB` and :func:`read_CIF` functions
-----------------------------------------------
.. automodule:: chemistry.structure
    :members: read_PDB, read_CIF

