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
created by one of the parsers (see below for :func:`load_file
<chemistry.formats.registry.load_file>`) or extended to
support the structure files of various computational chemistry programs.

:func:`load_file <chemistry.formats.registry.load_file>`
--------------------------------------------------------
.. currentmodule:: chemistry.formats.registry
.. autosummary::
    :toctree: structobj/

    load_file

The :func:`load_file` function automatically determines the format of the file
whose name is passed as an argument. The following formats are currently
recognized and result in the instantiation of either a :class:`Structure
<chemistry.structure.Structure>`. or one of its subclasses:

1. PDB
2. PDBx/mmCIF
3. Gromacs GRO
4. Gromacs topology file
5. Amber topology file
6. CHARMM PSF file
7. Mol2 file

Here, we will focus on instantiating a :class:`Structure
<chemistry.structure.Structure>` instance from PDB and
mmCIF files.  PDB files and mmCIF files downloaded from the RCSB Protein Data
Bank or the world wide Protein Data Bank often contain a large amount of
metadata describing the structure, such as the citation information,
experimental method (e.g., X-ray crystallography or NMR spectroscopy), authors,
and related database entries (such as BMRB entries for NMR-solved structures).
This information is extracted from both PDB and PDBx/mmCIF files when available,
along with anisotropic B-factors.

The following sections will briefly demonstrate parsing a PDB file and a mmCIF
file to a :class:`Structure <chemistry.structure.Structure>` instance.

Examples
--------

For the purposes of this example, we will download the 4LZT structure as both a
PDB file and a CIF file. These structures are both used in the ParmEd unittest
suite, so you can get the files from there, or `download them from here.`__

.. _4lzt_structures: http://www.rcsb.org/pdb/explore/explore.do?structureId=4LZT
__ 4lzt_structures_

------

The first example demonstrates reading the PDB file::

    >>> from chemistry import load_file
    >>> lzt_pdb = load_file('4lzt.pdb')
    >>> len(lzt_pdb.atoms)
    1164
    >>> len(lzt_pdb.residues)
    274
    >>> lzt_pdb.experimental # See the experimental method
    'X-RAY DIFFRACTION'
    >>> # See how many chains we have
    >>> chain_list = set()
    >>> for residue in lzt_pdb.residues:
    ...     chain_list.add(residue.chain)
    ... 
    >>> chain_list
    set(['A'])
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

    >>> from chemistry import load_file
    >>> lzt_cif = load_file('4LZT.cif')
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
:class:`Structure <chemistry.structure.Structure>` class, you may have deduced
that you can very simply convert a PDBx/mmCIF file to a PDB file.

This is likely to be increasingly popular, since the PDB is moving to the mmCIF
format, but many programs in the field of computational chemistry and physics
has decades worth of legacy code built around PDB files. Not to worry!  A quick
1-liner will seamlessly convert PDBx to PDB::

    from chemistry import load_file
    load_file('4LZT.cif').write_pdb('4lzt_converted.pdb', write_anisou=True,
                                    renumber=False)

In this case, the metadata is *not* copied (i.e., the ``EXPTL``, ``JRNL``,
and ``AUTHOR`` records, to name a few). Only the coordinates, unit cell
(``CYRST1`` record), and optionally anisotropic B-factor lines are translated.

The ``renumber`` argument tells :class:`Structure
<chemistry.structure.Structure>` to use the original PDB numbers, rather than
its internal number scheme that numbers sequentially from 1 to N, where N is
the number of residues.

If your system has more than 99,999 atoms (and/or more than 9,999 residues), the
numbering cycles back, such that the atom serial number after 99,999 is 0, and
the numbering starts again.

-------

Structure manipulation: slicing, combining, and replicating
-----------------------------------------------------------

This section describes a number of simple, yet powerful manipulations you can do
to :class:`Structure <chemistry.structure.Structure>` instances (and, by
extension, instances of their subclasses).  In order, these are *slicing*,
*combining* (or *merging*), and *replicating*.

:class:`Structure` Slicing and Selections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One of the things that ParmEd tries to make easy to do is to select a subset of
a :class:`Structure <chemistry.structure.Structure>` instance as a new
:class:`Structure <chemistry.structure.Structure>` (complete with all
remaining atoms, residues, and parameters) by using the common Python idiom of
`slicing <http://www.pythoncentral.io/how-to-slice-listsarrays-and-tuples-in-python/>`_.
The *slice*, or *selection* syntax of :class:`Structure
<chemistry.structure.Structure>` is designed to be flexible, expressive, and
intuitive (although accomplishing all of these is a challenge!).

There are three many ways to select from a :class:`Structure
<chemistry.structure.Structure>` instance:

1. By atom index
2. By Amber selection mask (see :class:`AmberMask <chemistry.amber.AmberMask>`)
3. By selection arrays, mask arrays, or slices of chains, residues, and/or
   atoms.

When selecting from a :class:`Structure <chemistry.structure.Structure>`
instance, the return value can be one of three things:

1. ``None`` if no atoms match the selection criteria
2. :class:`Atom <chemistry.topologyobjects.Atom>` instance if the selection
   matches only a single atom. This instance is a reference to the original
   atom, it is *not* a copy.
3. :class:`Structure <chemistry.structure.Structure>` with all of the selected
   atoms and all parameters that were present between the selected atoms. In
   this case, a copy is made of all selected atoms.

**Let's look at the simplest form of the selection syntax -- by atom index**::

    >>> struct = load_file('4LZT.cif') # use the 4LZT.cif file in test/files/
    >>> struct
    <Structure 1164 atoms; 274 residues; 0 bonds; PBC; NOT parametrized>
    >>> struct[0]
    <Atom N [0]; In LYS 0>
    >>> struct[10]
    <Atom CA [10]; In VAL 1>
    >>> struct[1163]
    <Atom O [1163]; In HOH 273>
    >>> struct[-1]  # Negative indices also work!
    <Atom O [1163]; In HOH 273>
    >>> struct[1164]
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/.../structure.py", line 1021, in __getitem__
        return self.atoms[selection]
      File "/.../topologyobjects.py", line 3932, in __getitem__
        retval = list.__getitem__(self, thing)
    IndexError: list index out of range

Oops!  Our original structure only had 1164 atoms, so selecting the 1165th index
(remember, Python indices start from 0) results in an ``IndexError``. As with
Python slicing, though, slice and Amber mask selections simply return ``None``
if no atoms match the selection.

**Now let's look at the slightly more complex Amber mask selection. See `Amber
mask syntax <amber.html#amber-mask-syntax>`_ for more details**.

If you pass a single string that is recognized as valid Amber mask syntax, it
will be interpreted as an Amber mask. The example below continues from the same
``struct`` we were using above::

    >>> struct['@CA']
    <Structure 129 atoms; 129 residues; 0 bonds; NOT parametrized>
    >>> struct[':1-10']
    <Structure 75 atoms; 10 residues; 0 bonds; NOT parametrized>
    >>> struct[':1-10@CA,CB']
    <Structure 19 atoms; 10 residues; 0 bonds; NOT parametrized>
    >>> struct[':10<@5']
    <Structure 49 atoms; 13 residues; 0 bonds; NOT parametrized>
    >>> struct[':10<:5']
    <Structure 88 atoms; 13 residues; 0 bonds; NOT parametrized>

**Now let's look at mask arrays**.

I'll define my terminology.  A *mask array* is a boolean array (or it is
interpreted as a boolean array) with a Truthy-value (e.g., ``True``, ``1``,
etc.) for atoms you want to select and a Falsey-value (e.g., ``False``, ``0``,
etc.) for atoms you do *not* want to select.  A *selection array*, described in
more detail below, is an array of atom indices (first index is 0) of the atoms
you want to select. These arrays can be any subscriptable iterable (e.g.,
``tuple``, ``list``, ``numpy.ndarray``, etc.).

A *mask array* must have the same number of elements as your system has atoms.
So if an array has the same number of elements as you have atoms, it is
interpreted as a mask array.  A common use-case for a mask array is to convert a
:class:`Structure <chemistry.structure.Structure>` to a Pandas ``DataFrame``
(see :func:`Structure.to_dataframe <chemistry.structure.Structure.to_dataframe>`
for more information), and using Pandas/numpy to generate the mask array. An
example of this is shown below::

    >>> df = struct.to_dataframe()
    >>> struct[df.name == 'CA']
    <Structure 129 atoms; 129 residues; 0 bonds; NOT parametrized>
    >>> struct[(df.name == 'CA')&(df.resid < 10)]
    <Structure 10 atoms; 10 residues; 0 bonds; NOT parametrized>

A mask array cannot be used alongside a residue or chain selection. In fact,
ParmEd will simply interpret that array as a selection array, which may lead to
unexpected results, as seen below::

    >>> struct['A',:,df.name == 'CA']
    <Structure 409 atoms; 274 residues; 0 bonds; NOT parametrized>

While you might expect this to select only the ``CA`` atoms (since all residues
are in chain A), the key is that ``df.name == 'CA'`` is interpreted as a
selection array, and it contains 0s and 1s. So it selects the first two atoms of
every residue (not all residues have more than 1 atom, though, which is why the
number of atoms is not double the number of residues).

**Now let's look at selection arrays and slices of chains, residues, and/or
atoms**

:class:`Structure <chemistry.structure.Structure>` instances can take between 1
and 3 "slots" in their indexing scheme, corresponding to atom, residue, and
chain selections.  If one slot is used, that selection applies to the list of
atoms.  If two slots are used, the first slot applies to the list of residues
while the second slot refers to the list of atoms *within the selected
residues*. Note the difference between the atom selections when a residue
selection is given compared to when one is not. If all three slots are used, the
first is interpreted as a chain selection, the second as a selection of residues
*within each selected chain* and the third as a selection of each atom *within
those residues*.

Slots can be assigned either an index, string, slice, or selection array. Each
is defined below:

- An *index* is a number between 0 and the number of elements in either the
  residue or atom list (depending on which slot is found in).
- A *string* is interpreted as a single name that must match every chain,
  residue, or atom name (as determined by the slot it is found in).
- A *slice* is a standard Python slice (e.g., ``iterable[10:20:2]``) and can be
  used for either the residue or atom selections. A raw ``:`` in any slot means
  a "full slice" that selects everything according to its slot.
- A *selection array* is an array of indexes or strings (interpreted as names),
  and can have as many elements as you want (including negative numbers to count
  from the *end* of the atom list). but no atom index can be outside the range
  of the atom list (or you will get an ``IndexError``). So practically speaking,
  there is no use-case for a selection array to have the same length as a mask
  array, since that would either include duplicates *or* select every atom.

Let's have a look at some examples, again continuing with the ``struct`` object
we defined above::

    >>> struct[10:20:2]
    <Structure 5 atoms; 2 residues; 0 bonds; NOT parametrized>
    >>> struct[0,10:20:2] # note, this is the 10th to 20th atom of the 1st residue!
    >>> struct[:,0] # First atom of every residue
    <Structure 274 atoms; 274 residues; 0 bonds; NOT parametrized>
    >>> struct['A',:,:] # All atoms in chain A
    <Structure 1164 atoms; 274 residues; 0 bonds; NOT parametrized>
    >>> struct[:,'CA'] # All atoms named CA in all residues
    <Structure 129 atoms; 129 residues; 0 bonds; NOT parametrized>

There is so much flexiblity in the Atom selection here that we can't possibly
cover everything. You are encouraged to try things out!

Structure Combining
~~~~~~~~~~~~~~~~~~~

Two structures can be *combined* in the sense that a :class:`Structure
<chemistry.structure.Structure>` instance can be created by taking the atoms,
residues, and parameters of one ``Structure`` and tacking it on to the end of
another one.

This can naturally be thought of as a *sum* of two ``Structure`` instances, so
it was implemented via the addition operator.  This can be done both in-place
(in a way that modifies the first ``Structure``) as well as creating a new copy
that is the sum of the originals.
