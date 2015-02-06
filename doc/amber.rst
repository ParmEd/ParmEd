The :mod:`chemistry.amber` package
==================================

The :mod:`amber` package contains classes that can parse most of the file
formats used by the Amber molecular dynamics package.  In particular are Amber
parameter and topology files as well as the various formats of coordinate and
trajectory files that Amber supports (and creates).

It is out of the scope of this site to discuss the details of the file formats
themselves, and you are instead forwarded to the `official Amber website
<http://ambermd.org>`_ and `file format specifications
<http://ambermd.org/formats.html>`_\ .

The parameter-topology (*prmtop*) file
--------------------------------------

The primary file in Amber defining the system topology and the parameters
defining the force field for that system is called the *prmtop* file, whose
format is detailed on http://ambermd.org/formats.html. There are several
variants of this extensible file format, and these parsers rigorously implement
the format specification.

.. currentmodule:: chemistry.amber
.. autosummary::
    :toctree: amberobj/

    AmberFormat
    AmberParm
    ChamberParm
    AmoebaParm
    LoadParm
    AmberMask

The :class:`AmberFormat` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class simply parses the raw data from a file with the format of an Amber
*prmtop* file and makes the data available in a ``dict`` whose keys are the
``FLAG_NAME`` in the ``%FLAG <FLAG_NAME>`` lines. Some files---like the MDL
files used in 3D-RISM calculations---use this format to store raw, structured
data, and this class is appropriate for use with those files.

In most cases, however, this simply serves as a base class, along with
:class:`Structure <chemistry.structure.Structure>`, for the *prmtop* parsers.

An example using this class on `one of the test files in the ParmEd repostory
<https://github.com/swails/ParmEd/blob/master/test/files/cSPCE.mdl>`_ is shown
below::

    >>> amb = AmberFormat('cSPCE.mdl')
    >>> amb.parm_data['ATMTYP']
    [1, 2]
    >>> amb.parm_data['COORD']
    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.333314, 0.942816, 0.0]

If you want to, you can make modifications to the data and then write a new
file::

    >>> amb.parm_data['COORD'][0] = 10
    >>> amb.write_parm('new.mdl')
    >>> # Let's check that our new file has the changed COORD
    >>> new_amb = AmberFormat('new.mdl')
    >>> new_amb.parm_data['COORD']
    [10.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.333314, 0.942816, 0.0]

The :class:`AmberParm`, :class:`ChamberParm`, and :class:`AmoebaParm` classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These classes are subclasses of both :class:`Structure
<chemistry.structure.Structure>` and :class:`AmberFormat`. (The special classes
:class:`ChamberParm` and :class:`AmoebaParm`, implementing the CHARMM and Amoeba
force fields, respectively, both inherit from :class:`AmberParm`.) You need to
make sure that you use the correct class for the variant of *prmtop* file you
are using.  If you are unsure, the :func:`LoadParm` function will automatically
determine the correct class and return an instance of it (see below).

The :class:`AmberParm` class parses a *prmtop* object and stores both the raw
data (see the :class:`AmberFormat` description, above), as well as the standard
topological data structures used in the :class:`Structure
<chemistry.structure.Structure>` class. Different subclasses of
:class:`AmberParm` will populate different term and parameter lists depending on
the terms present in the force field.

Since the raw data and topological features are inextricably linked, there are
risks associated with modifying the raw data in :py:attr:`AmberParm.parm_data`
*and* the parameters and properties of the topological data structures in the
lists inherited from :class:`Structure <chemistry.structure.Structure>`.

The required overhead to keep :py:attr:`AmberParm.parm_data` and the topology
structures synchronized at all times introduces too great a cost, so keeping
them synchronized largely falls on the programmer, although there is
functionality introduced to help.

Generalizing with :func:`LoadParm`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This will automatically read and cast to the proper parm type.

The coordinate and trajectory files
-----------------------------------

.. currentmodule:: chemistry.amber
.. autosummary::
    :toctree: amberobj/

    Rst7
    AmberAsciiRestart
    NetCDFRestart
    NetCDFTraj
    AmberMdcrd
    netcdffiles.use
