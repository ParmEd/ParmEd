The :mod:`structure` module
===========================

The :mod:`structure` module contains the core class defining systems of
particles -- the :class:`Structure` class. The API documentation makes this
class look intimidating, but its core features are the ``atoms`` attribute
(which is an :class:`AtomList`) and its ``residues`` attribute (which is a
:class:`ResidueList`). These two attributes are enough to write numerous files
(e.g., PDB and PDBx/mmCIF), and these are the only two attributes currently
populated through the PDB and mmCIF parsing routines.

:class:`Structure` class
------------------------
.. currentmodule:: chemistry.structure
.. autosummary::
    :toctree: structobj/

    Structure

The :class:`Structure` class may be instantiated directly, but is more often
created by one of the parsers (see below for :func:`read_PDB` and
:func:`read_CIF` which both return a :class:`Structure` instance) or extended to
support the structure files of various computational chemistry programs.

:func:`read_PDB` and :func:`read_CIF` functions
-----------------------------------------------
.. currentmodule:: chemistry.structure
.. autosummary::
    :toctree: structobj/

    read_PDB
    read_CIF

PDB files and mmCIF files downloaded from the RCSB Protein Data Bank or the
world wide Protein Data Bank often contain a large amount of metadata describing
the structure, such as the citation information, experimental method (e.g.,
X-ray crystallography or NMR spectroscopy), authors, and related database
entries (such as BMRB entries for NMR-solved structures). This information is
extracted from both PDB and PDBx/mmCIF files when available, along with
anisotropic B-factors.

The following sections will briefly demonstrate parsing a PDB file and a mmCIF
file to a :class:`Structure` instance.

Examples
--------

For the purposes of this example, we will download the 4LZT structure as both a
PDB file and a CIF file. These structures are both used in the ParmEd unittest
suite, so you can get the files from there, or `download them from here.`__

.. _4lzt_structures: http://www.rcsb.org/pdb/explore/explore.do?structureId=4LZT
__ 4lzt_structures_

------

The first example demonstrates reading the PDB file::

    >>> from chemistry import read_PDB
    >>> lzt_pdb = read_PDB('4lzt.pdb')
    >>> len(lzt_pdb.atoms)
    1164
    >>> len(lzt_pdb.residues)
    274
    >>> lzt_pdb.experimental # See the experimental method
    u'X-RAY DIFFRACTION'
    >>> # See how many chains we have
    >>> chain_list = set()
    >>> for residue in lzt_pdb.residues:
    ...     chain_list.add(residue.chain)
    ... 
    >>> chain_list
    set([u'A'])
    >>> # Only one chain. Now see how many atoms have alternate conformations
    >>> atoms_with_altconf = [atom for atom in lzt_pdb.atoms
    ...                           if atom.other_locations]
    >>> len(atoms_with_altconf)
    19
    >>> # Just for sanity's sake, make sure that the sum of all atoms in all
    >>> # residues is equal to the total number of atoms
    >>> sum(len(residue) for residue in lzt_pdb.residues)
    1164
    >>> len(lzt_pdb.atoms)
    1164
    >>> lzt_pdb.atoms[0].anisou # Look, we even get anisotropic B-factors!
    array([ 0.2066,  0.1204,  0.1269,  0.0044,  0.0126,  0.0191])
    >>> 

The second example demonstrates reading the CIF file::

    >>> from chemistry import read_CIF
    >>> lzt_cif = read_CIF('4LZT.cif')
    >>> len(lzt_cif.atoms)
    1164
    >>> len(lzt_cif.residues)
    274
    >>> lzt_cif.experimental # See the experimental method
    'X-RAY DIFFRACTION'
    >>> # See how many chains we have
    >>> chain_list = set()
    >>> for residue in lzt_cif.residues:
    ...     chain_list.add(residue.chain)
    ... 
    >>> chain_list
    set(['A'])
    >>> # Only one chain. Now see how many atoms have alternate conformations
    >>> atoms_with_altconf = [atom for atom in lzt_cif.atoms
    ...                          if atom.other_locations]
    >>> len(atoms_with_altconf)
    19
    >>> # Just for sanity's sake, make sure that the sum of all atoms in all
    >>> # residues is equal to the total number of atoms
    >>> sum(len(residue) for residue in lzt_cif.residues)
    1164
    >>> len(lzt_pdb.atoms)
    1164
    >>> lzt_cif.atoms[0].anisou # Look, we even get anisotropic B-factors!
    array([ 0.2066,  0.1204,  0.1269,  0.0044,  0.0126,  0.0191])
    >>> 

-------

Converting PDBx/mmCIF files to PDB files
----------------------------------------

If you noticed the ``write_cif`` and ``write_pdb`` methods attached to the
:class:`Structure` class, you may have deduced that you can very simply convert
a PDBx/mmCIF file to a PDB file.

This is likely to be increasingly popular, since the PDB is moving to the mmCIF
format, but many programs in the field of computational chemistry and physics
has decades worth of legacy code built around PDB files. Not to worry!  A quick
1-liner will seamlessly convert PDBx to PDB::

    from chemistry import read_CIF
    read_CIF('4LZT.cif').write_pdb('4lzt_converted.pdb', write_anisou=True,
                                   renumber=False)

In this case, the metadata is *not* copied (i.e., the ``EXPTL``, ``JRNL``,
and ``AUTHOR`` records, to name a few). Only the coordinates, unit cell
(``CYRST1`` record), and optionally anisotropic B-factor lines are translated.

The ``renumber`` argument tells :class:`Structure` to use the original PDB
numbers, rather than its internal number scheme that numbers sequentially from 1
to N, where N is the number of residues.

If your system has more than 99,999 atoms (and/or more than 9,999 residues), the
numbering cycles back, such that the atom serial number after 99,999 is 0, and
the numbering starts again.
