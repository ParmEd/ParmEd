File input/output with ParmEd
=============================

ParmEd has several convenience functions designed specifically to make file
input/output easier for a wide range of file formats.

The API for parsing different file types is not always consistent due to the
technical needs of each class.  For instance, the Amber topology file class
subclasses from :class:`Structure <parmed.structure.Structure>`, so it can be
instantiated directly from a filename. On the other hand, :class:`PDBFile
<parmed.formats.pdb.PDBFile>` contains a ``parse`` static method that returns a
:class:`Structure <parmed.structure.Structure>` instance directly rather than
subclassing.

Reading files with :func:`load_file <parmed.formats.registry.load_file>`
------------------------------------------------------------------------

.. currentmodule:: parmed.formats.registry
.. autosummary::
    :toctree: fileio/

    load_file

To provide a single interface for parsing *any* type of file, the
:func:`load_file <parmed.formats.registry.load_file>` function takes a filename
and any relevant arguments or keyword arguments for the supported file formats
and calls the appropriate parsing function.  Supported file formats along with
the extra supported keyword arguments (if applicable) are shown below:

* Amber ASCII restart files
* Amber topology/MDL files (``xyz=None`` is the coordinate file name or numpy
  array of coordinates, ``box=None`` is the array of unit cell lengths and
  angles)
* Amber ASCII trajectory (MDCRD) files (``natom`` is the *required* number of
  atoms in the system, ``hasbox`` is the *required* boolean flag specifying
  whether or not unit cell dimensions are present in the mdcrd file)
* Amber NetCDF restart
* Amber NetCDF trajectory file
* Amber OFF residue template library files
* PDBx/mmCIF files
* CHARMM coordinate files
* CHARMM restart files
* Gromacs GRO files
* Gromacs topology files (``defines=None`` is the (ordered) ``dict`` of
  preprocessor defines to apply to the processing of the topology file,
  ``parametrize=True`` is a flag to specify whether or not to apply parameters,
  ``xyz=None`` is the file name of a coordinate file or an array of coordinates,
  and ``box=None`` is the array of unit cell dimensions and angles)
* SYBYL mol2 files and the R.E.D. mol3 extension (``structure=False`` is a flag
  that determines whether or not the returned object is a :class:`Structure
  <parmed.structure.Structure>` or :class:`ResidueTemplate
  <parmed.modeller.residue.ResidueTemplate>`/:class:`ResidueTemplateContainer
  <parmed.modeller.residue.ResidueTemplateContainer>` instance)
* PDB files
* CHARMM/XPLOR PSF files

:func:`load_file` automatically inspects the contents of the file with the given
name to determine what format the file is based on the first couple lines.
Except in rare, pathological cases, the file format detection mechanism is
fairly robust. If any files *fail* this detection, feel free to file an issue on
the Github issue tracker to improve file type prediction.

Writing files with :meth:`Structure.save <parmed.structure.Structure.save>`
---------------------------------------------------------------------------

Many of the file formats supported by ParmEd either parse directly to a
:class:`Structure <parmed.structure.Structure>` instance or subclass, and many
of the desired file type conversions that ParmEd is designed to facilitate are
between these formats (e.g., Amber topology, PDB, CHARMM PSF file, etc.).

To facilitate the required conversion and file writing, the base
:class:`Structure <parmed.structure.Structure>` class has a ``save`` method that
will convert to the requested file format and write the output file. The desired
format is specified either explicitly or by file name extension (with explicit
format specifications taking precedence).  The supported file formats, along
with their supported extra keyword arguments, are detailed in the following
table.

+--------------+-------------------------+----------------+--------------------------+
| File type    | Recognized extension(s) | Format keyword | Supported arguments      |
+==============+=========================+================+==========================+
| PDB          | ``.pdb``                | ``pdb``        |  ``charmm``\*,           |
|              |                         |                |  ``renumber``,           |
+--------------+-------------------------+----------------+  ``coordinates``,        |
| PDBx/mmCIF   | ``.cif``, ``.pdbx``     | ``cif``        |  ``altlocs``,            |
|              |                         |                |  ``write_anisou``        |
+--------------+-------------------------+----------------+--------------------------+
| Amber prmtop | ``.parm7``, ``.prmtop`` | ``amber``      |  None                    |
+--------------+-------------------------+----------------+--------------------------+
| CHARMM PSF   | ``.psf``                | ``charmm``     | ``vmd``                  |
+--------------+-------------------------+----------------+--------------------------+
| Gromacs      | ``.top``                | ``gromacs``    | ``combine``,             |
| topology     |                         |                | ``parameters``           |
+--------------+-------------------------+----------------+--------------------------+
| Gromacs GRO  | ``.gro``                | ``gro``        | ``precision``, ``nobox`` |
+--------------+-------------------------+----------------+--------------------------+
| Mol2         | ``.mol2``               | ``mol2``       | ``split``                |
+--------------+-------------------------+----------------+--------------------------+
| Mol3         | ``.mol3``               | ``mol3``       | ``split``                |
+--------------+-------------------------+----------------+--------------------------+

\* PDB format only

The meanings and default values of each of the keywords is described in the next
subsection.

Keywords
~~~~~~~~

* ``charmm`` -- If ``True``, the SEGID will be printed in columns 73 to 76 of
  the PDB file (default is ``False``)
* ``renumber`` -- If ``True``, atoms and residues will be numbered according to
  their absolute index (starting from 1) in the system. If ``False``, the
  numbers will be retained from their original source (e.g., in the original PDB
  file). Default is ``True``
* ``coordinates`` -- A set of coordinates for one or multiple frames. If more
  than one frame is provided, the resulting PDB or PDBx/mmCIF file will have
  multiple models defined. Default is ``None``, and the generated coordinates
  are the ones stored on the atoms themselves.
* ``altlocs`` -- Allowable values are ``all`` (print all alternate locations for
  all conformers), ``first`` (print only the first alternate conformer), and
  ``occupancy`` (print the conformer with the highest occupancy). Default is
  ``all``
* ``write_anisou`` -- If ``True``, print anistropic B-factors for the various
  atoms (either as a separate CIF section or as ``ANISOU`` records in a PDB
  file). Default is ``False``
* ``vmd`` -- If ``True``, write a VMD-style PSF file. This is very similar to
  XPLOR format PSF files. Default is ``False``.
* ``combine`` -- Can be ``None`` to combine no molecules when writing a GROMACS
  topology file. A value of ``all`` will combine all of the molecules into a
  single moleculetype. Default is ``None``.
* ``parameters`` -- Can be ``inline`` (write the parameters inline in the
  GROMACS topology file). Can also be a string or a file-like object. If it is
  the same as the topology file name, it will be written in the previous
  sections of the GROMACS topology file. Other strings will be interpreted as
  filenames to print the parameters to as an include topology file. Default is
  ``inline``.
* ``precision`` -- The number of decimal places to print coordinates with in GRO
  files. Default is 3.
* ``nobox`` -- If ``True`` and the ``Structure`` does not have a unit cell
  defined, no box is written to the bottom of the GRO file. Otherwise, an
  enclosing box (buffered by 5 angstroms) is written
* ``split`` -- If ``True``, all residues will be split into separate mol2 or
  mol3 entries in the same file (like the ZINC database, for example). If
  ``False``, all residues will be part of the same mol2 or mol3 entry. Default
  is ``False``.
