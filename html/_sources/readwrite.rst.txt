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
the extra supported keyword arguments (if applicable) are shown in the table
below (the filename argument is omitted):

+--------------------------+--------------------------------+
| File type                | Supported arguments            |
+==========================+================================+
| Amber ASCII restart file | None                           |
+--------------------------+--------------------------------+
| Amber prmtop             | ``xyz``, ``box``               |
+--------------------------+--------------------------------+
| Amber MDL (e.g., RISM)   | None                           |
+--------------------------+--------------------------------+
| Amber MDCRD trajectory   | ``natom``\*, ``hasbox``\*      |
+--------------------------+--------------------------------+
| Amber OFF library        | None                           |
+--------------------------+--------------------------------+
| Amber frcmod/parm.dat    | list of filenames              |
+--------------------------+--------------------------------+
| PDBx/mmCIF file          | ``skip_bonds``                 |
+--------------------------+--------------------------------+
| CHARMM coordinate file   | None                           |
+--------------------------+--------------------------------+
| CHARMM restart file      | None                           |
+--------------------------+--------------------------------+
| Gromacs GRO file         | ``skip_bonds``                 |
+--------------------------+--------------------------------+
| Gromacs topology file    | ``defines``, ``parametrize``   |
|                          | ``xyz``, ``box``               |
+--------------------------+--------------------------------+
| Mol2 and Mol3 files      | ``structure``                  |
+--------------------------+--------------------------------+
| NetCDF restart file      | None                           |
+--------------------------+--------------------------------+
| NetCDF trajectory file   | None                           |
+--------------------------+--------------------------------+
| PDB file                 | ``skip_bonds``                 |
+--------------------------+--------------------------------+
| PQR file                 | ``skip_bonds``                 |
+--------------------------+--------------------------------+
| PSF file                 | None                           |
+--------------------------+--------------------------------+
| Serialized System XML    | None                           |
+--------------------------+--------------------------------+
| Serialized State XML     | None                           |
+--------------------------+--------------------------------+
| Serialized Integrator    | None                           |
| XML file                 |                                |
+--------------------------+--------------------------------+
| ForceField XML file      | None                           |
+--------------------------+--------------------------------+
\*These arguments are *required* when parsing the corresponding format file

The optional keyword arguments are described below:

* ``xyz`` -- Either a file name for a coordinate file or an array with the
  coordinates. If the unit cell information is stored in the coordinate file, it
  is used (unless an explicit ``box`` argument is given; see below). If the file
  is a trajectory file, the first frame is used for the coordinates.
* ``box`` -- The unit cell dimensions in Angstroms (and angles in degrees).
* ``natom`` -- The number of atoms from which the file was written
* ``hasbox`` -- Whether unit cell information is stored in this trajectory file
* ``defines`` -- ``dict`` of preprocessor defines (order is respected if given
  via an ``OrderedDict``)
* ``parametrize`` -- If ``True``, parameters are assigned from the parameter
  database.  If ``False``, they are not (default is ``True``).
* ``structure`` -- If ``True``, return the Mol2/Mol3 file as a :class:`Structure
  <parmed.structure.Structure>` instance. Default is ``False``
* ``skip_bonds`` -- If ``True``, no attempt is made to identify covalent bonds
  between any atoms for these formats. If ``False`` (default), bonds are first
  assigned by comparing to standard residue templates, then based on distance
  criteria for any atom not specified within the templates. This has a
  side-effect of improving element identification for PDB files that have no
  element columns specified as well as substantially improving element
  identification for ions in GRO files. You are suggested to keep the default
  (``False``) unless you are *only* using the coordinates, in which case setting
  ``skip_bonds`` to ``True`` will result in significantly improved performance.

:func:`load_file` automatically inspects the contents of the file with the given
name to determine what format the file is based on the first couple lines.
Except in rare, pathological cases, the file format detection mechanism is
fairly robust. If any files *fail* this detection, feel free to file an issue on
the Github issue tracker to improve file type prediction.

:func:`load_file` has a number of helpful features. For instance, files ending
with the ``.gz`` or ``.bz2`` suffix will automatically be decompressed in-memory
using Gzip or Bzip2, respectively (except for some binary file formats, like
NetCDF). Furthermore, URLs beginning with ``http://``, ``https://``, or
``ftp://`` are valid file names and will result in the remote file being
downloaded and processed (again, in-memory).

Finally, to make it so that you can always retrieve a :class:`Structure
<parmed.structure.Structure>` instance from file types that support returning
one by passing the ``structure=True`` keyword to ``load_file``. If this argument
is not supported by the resulting file type, it is simply ignored, as is the
``natom``, ``hasbox``, and ``skip_bonds`` keywords.

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
format specifications taking precedence).  Because :meth:`Structure.save
<parmed.structure.Structure.save>` is a convenience method, it will protect
against accidentally overwriting an existing file. The ``overwrite`` argument,
when set to ``True``, will allow an existing file to be overwritten. If set to
``False`` (or left at its default), ``IOError`` will be raised when attempting
to overwrite an existing file. The supported file formats, along with their
supported extra keyword arguments, are detailed in the following table.

+--------------+-------------------------+----------------+--------------------------+
| File type    | Recognized extension(s) | Format keyword | Supported arguments      |
+==============+=========================+================+==========================+
| PDB          | ``.pdb``                | ``pdb``        |  ``charmm``\*,           |
|              |                         |                |  ``renumber``,           |
+--------------+-------------------------+----------------+  ``coordinates``,        |
| PDBx/mmCIF   | ``.cif``, ``.pdbx``     | ``cif``        |  ``altlocs``,            |
|              |                         |                |  ``write_anisou``,       |
|              |                         |                |  ``standard_resnames``   |
+--------------+-------------------------+----------------+--------------------------+
| PQR          | ``.pqr``                | ``pqr``        | ``renumber``,            |
|              |                         |                | ``coordinates``,         |
|              |                         |                | ``standard_resnames``    |
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
| Amber ASCII  | ``.rst7``, ``.inpcrd``, | ``rst7``       | ``title``, ``time``      |
| coordinates  | ``.restrt``             |                |                          |
+--------------+-------------------------+----------------+--------------------------+
| Amber NetCDF | ``.ncrst``              | ``ncrst``      | ``title``, ``time``      |
| restart      |                         |                |                          |
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
* ``standard_resnames`` -- If ``True``, residue names will be regularlized from
  common alternatives back to the PDB standard. For example, ASH and GLH will be
  translated to ASP and GLU, respectively, as they often refer to different
  protomers of aspartate and glutamate.
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
* ``title`` -- Purely cosmetic, it will specify the title that will be written
  to the coordinate files
* ``time`` -- Also cosmetic, this is the time corresponding to the snapshot that
  will be written to the coordinate files

Examples
--------

The following examples use various files from the ParmEd test suite, which can
be found in the ``test/files/`` directory of the ParmEd source code::

    >>> import parmed as pmd
    >>> # Load a Mol2 file
    ... pmd.load_file('tripos1.mol2')
    <ResidueTemplate DAN: 31 atoms; 33 bonds; head=None; tail=None>
    >>> # Load a Mol2 file as a Structure
    ... pmd.load_file('tripos1.mol2', structure=True)
    <Structure 31 atoms; 1 residues; 33 bonds; NOT parametrized>
    >>> # Load an Amber topology file
    ... parm = pmd.load_file('trx.prmtop', xyz='trx.inpcrd')
    >>> parm
    <AmberParm 1654 atoms; 108 residues; 1670 bonds; parametrized>
    >>> # Load a CHARMM PSF file
    ... psf = pmd.load_file('ala_ala_ala.psf')
    >>> psf
    <CharmmPsfFile 33 atoms; 3 residues; 32 bonds; NOT parametrized>
    >>> # Load a PDB and CIF file
    ... pdb = pmd.load_file('4lzt.pdb')
    >>> cif = pmd.load_file('4LZT.cif')
    >>> pdb
    <Structure 1164 atoms; 274 residues; 1043 bonds; PBC (triclinic); NOT parametrized>
    >>> cif
    <Structure 1164 atoms; 274 residues; 1043 bonds; PBC (triclinic); NOT parametrized>
    >>> # Load a Gromacs topology file -- only works with Gromacs installed
    ... top = pmd.load_file('1aki.ff99sbildn.top')
    >>> top
    <GromacsTopologyFile 40560 atoms [9650 EPs]; 9779 residues; 30934 bonds; parametrized>

Any of the :class:`Structure <parmed.structure.Structure>` subclasses shown
above can be saved as any other kind of :class:`Structure
<parmed.structure.Structure>` class or subclass *for which the conversion is
supported*. For instance, a raw PDB file has no defined parameters, so it cannot
be saved as an Amber topology file. An Amber topology file, on the other hand,
has all of the information required for a PDB, and so that conversion is
supported::

    >>> parm.save('test_parm.pdb')
    >>> # You can also convert CIF to PDB
    ... cif.save('test_cif.pdb')
    >>> # Or you can convert PDB to CIF
    ... pdb.save('test_pdb.cif')
    >>> # Check the resulting saved files
    ... pmd.load_file('test_parm.pdb')
    <Structure 1654 atoms; 108 residues; 1670 bonds; NOT parametrized>
    >>> pmd.load_file('test_cif.pdb')
    <Structure 1164 atoms; 274 residues; 1043 bonds; PBC (triclinic); NOT parametrized>
    >>> pmd.load_file('test_pdb.cif')
    <Structure 1164 atoms; 274 residues; 1043 bonds; PBC (triclinic); NOT parametrized>
