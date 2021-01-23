"""
Provides a Python class for parsing a PSF file and setting up a system
structure for it

Author: Jason M. Swails
Contributors:
Date: April 20, 2014
"""
from __future__ import division

from copy import copy as _copy
from ..topologyobjects import (Bond, Angle, Dihedral, Improper, AcceptorDonor, Group, Cmap,
                               UreyBradley, NoUreyBradley, Atom, DihedralType, ImproperType,
                               UnassignedAtomType)
from ..exceptions import CharmmError, CharmmWarning, ParameterError
from ..structure import needs_openmm, Structure
from ..utils.io import genopen
from ..utils.six import wraps
from ..utils.six.moves import zip, range
from ..utils.six import string_types
import re
import warnings

def _catchindexerror(func):
    """
    Protects a function from raising an index error, and replace that exception
    with a CharmmError instead
    """
    @wraps(func)
    def newfunc(*args, **kwargs):
        """ Catch the index error """
        try:
            return func(*args, **kwargs)
        except IndexError as e:
            raise CharmmError('Array is too short: %s' % e)

    return newfunc

class _FileEOF(Exception):
    """ For control flow """

class _ZeroDict(dict):
    """
    Contains a dict that returns dummy (zero) arguments when a key is not
    present rather than raising a KeyError.  The return value for non-existent
    items is (0, []). It also special-case sections that have multiple pointers
    to avoid index errors if those are not present in the PSF file
    """
    def __getitem__(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            if key.startswith('NGRP'):
                for k in self:
                    if k.startswith('NGRP'):
                        return dict.__getitem__(self, k)
                return [0, 0], []
            elif key.startswith('NUMLP'):
                for k in self:
                    if k.startswith('NUMLP'):
                        return dict.__getitem__(self, k)
                return [0, 0], []
            return 0, []

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_resre = re.compile(r'(-?\d+)([a-zA-Z]*)')

class CharmmPsfFile(Structure):
    """
    A chemical :class:`Structure` instantiated from CHARMM files.

    Parameters
    ----------
    psf_name : str, optional
        Name of the PSF file (it must exist)

    Raises
    ------
    IOError : If file ``psf_name`` does not exist
    CharmmPsfError : If any parsing errors are encountered
    """
    @staticmethod
    def _convert(string, type, message):
        """
        Converts a string to a specific type, making sure to raise
        CharmmError with the given message in the event of a failure.

        Parameters
        ----------
        string : str
            Input string to process
        type : type
            Type of data to convert to (e.g., ``int``)
        message : str
            Error message to put in exception if failed
        """
        try:
            return type(string)
        except ValueError:
            raise CharmmError('Could not convert %s [%s]' % (message,string))

    #===================================================

    @classmethod
    def _parse_psf_title_line(cls, line):
        words = line[:line.index('!')].split()
        title = line[line.index('!')+1:].strip().upper()
        # Strip out description
        if ':' in title:
            title = title[:title.index(':')]
        if len(words) == 1:
            pointers = cls._convert(words[0], int, 'pointer')
        else:
            pointers = tuple([cls._convert(w, int, 'pointer') for w in words])
        return title, pointers

    #===================================================

    @classmethod
    def _parse_psf_section(cls, psf):
        """
        This method parses a section of the PSF file

        Parameters
        ----------
        psf : file
            Open file that is pointing to the first line of the section that is
            to be parsed

        Returns
        -------
        title : str
            The label of the PSF section we are parsing
        pointers : (int/tuple of ints)
            If one pointer is set, pointers is simply the integer that is value
            of that pointer. Otherwise it is a tuple with every pointer value
            defined in the first line
        data : list
            A list of all data in the parsed section converted to integers
        """
        line = psf.readline()
        while not line.strip():
            if not line:
                raise _FileEOF('Unexpected EOF in PSF file')
            else:
                line = psf.readline()
        if '!' in line:
            title, pointers = cls._parse_psf_title_line(line)
        else:
            raise CharmmError('Could not determine section title') # pragma: no cover
        line = psf.readline().strip()
        if not line and title.startswith('NNB'):
            # This will correctly handle the NNB section (which has a spurious
            # blank line) as well as any sections that have 0 members.
            line = psf.readline().strip()
            # If this line has a title in it, then it's one of those weird PSF files that
            # has basically no NNB section. Just skip over NNB and take the next section
            # instead
            if '!' in line:
                title, pointers = cls._parse_psf_title_line(line)
                line = psf.readline().strip()
        data = []
        if title == 'NATOM' or title == 'NTITLE':
            # Store these two sections as strings (ATOM section we will parse
            # later). The rest of the sections are integer pointers
            while line:
                data.append(line)
                line = psf.readline().strip()
        else:
            while line:
                words = line.split()
                data.extend([cls._convert(w, int, 'PSF data') for w in words])
                line = psf.readline().strip()
        return title, pointers, data

    #===================================================

    @_catchindexerror
    def __init__(self, psf_name=None):
        """
        Opens and parses a PSF file, then instantiates a CharmmPsfFile
        instance from the data.
        """
        global _resre
        Structure.__init__(self)
        # Bail out if we don't have a filename
        if psf_name is None:
            return
        # Open the PSF and read the first line. It must start with "PSF"
        if isinstance(psf_name, string_types):
            fileobj = genopen(psf_name, 'r')
            own_handle = True
        else:
            fileobj = psf_name
            own_handle = False
        try:
            self.name = psf_name if isinstance(psf_name, string_types) else ''
            line = fileobj.readline()
            if not line.startswith('PSF'):
                raise CharmmError('Unrecognized PSF file. First line is %s' % line.strip())
            # Store the flags
            psf_flags = line.split()[1:]
            # Now get all of the sections and store them in a dict
            fileobj.readline()
            # Now get all of the sections
            psfsections = _ZeroDict()
            while True:
                try:
                    sec, ptr, data = CharmmPsfFile._parse_psf_section(fileobj)
                except _FileEOF:
                    break
                psfsections[sec] = (ptr, data)
            # store the title
            self.title = psfsections['NTITLE'][1]
            # Next is the number of atoms
            natom = self._convert(psfsections['NATOM'][0], int, 'natom')
            # Parse all of the atoms
            for i in range(natom):
                words = psfsections['NATOM'][1][i].split()
                atid = int(words[0])
                if atid != i + 1:
                    raise CharmmError('Nonsequential atoms detected!')
                segid = words[1]
                rematch = _resre.match(words[2])
                if not rematch:
                    raise CharmmError('Could not interpret residue number %s' % # pragma: no cover
                                      words[2])
                resid, inscode = rematch.groups()
                resid = self._convert(resid, int, 'residue number')
                resname = words[3]
                name = words[4]
                attype = words[5]
                # Try to convert the atom type to an integer a la CHARMM
                try:
                    attype = int(attype)
                except ValueError:
                    pass
                charge = self._convert(words[6], float, 'partial charge')
                mass = self._convert(words[7], float, 'atomic mass')
                props = words[8:]
                atom = Atom(name=name, type=attype, charge=charge, mass=mass)
                atom.props = props
                self.add_atom(atom, resname, resid, chain=segid,
                              inscode=inscode, segid=segid)
            # Now get the number of bonds
            nbond = self._convert(psfsections['NBOND'][0], int, 'number of bonds')
            if len(psfsections['NBOND'][1]) != nbond * 2:
                raise CharmmError('Got %d indexes for %d bonds' % # pragma: no cover
                                  (len(psfsections['NBOND'][1]), nbond))
            it = iter(psfsections['NBOND'][1])
            for i, j in zip(it, it):
                self.bonds.append(Bond(self.atoms[i-1], self.atoms[j-1]))
            # Now get the number of angles and the angle list
            ntheta = self._convert(psfsections['NTHETA'][0], int, 'number of angles')
            if len(psfsections['NTHETA'][1]) != ntheta * 3:
                raise CharmmError('Got %d indexes for %d angles' % # pragma: no cover
                                  (len(psfsections['NTHETA'][1]), ntheta))
            it = iter(psfsections['NTHETA'][1])
            for i, j, k in zip(it, it, it):
                self.angles.append(
                        Angle(self.atoms[i-1], self.atoms[j-1], self.atoms[k-1])
                )
                self.angles[-1].funct = 5 # urey-bradley
            # Now get the number of torsions and the torsion list
            nphi = self._convert(psfsections['NPHI'][0], int, 'number of torsions')
            if len(psfsections['NPHI'][1]) != nphi * 4:
                raise CharmmError('Got %d indexes for %d torsions' % # pragma: no cover
                                  (len(psfsections['NPHI']), nphi))
            it = iter(psfsections['NPHI'][1])
            for i, j, k, l in zip(it, it, it, it):
                self.dihedrals.append(
                        Dihedral(self.atoms[i-1], self.atoms[j-1],
                                 self.atoms[k-1], self.atoms[l-1])
                )
            self.dihedrals.split = False
            # Now get the number of improper torsions
            nimphi = self._convert(psfsections['NIMPHI'][0], int, 'number of impropers')
            if len(psfsections['NIMPHI'][1]) != nimphi * 4:
                raise CharmmError('Got %d indexes for %d impropers' % # pragma: no cover
                                  (len(psfsections['NIMPHI'][1]), nimphi))
            it = iter(psfsections['NIMPHI'][1])
            for i, j, k, l in zip(it, it, it, it):
                self.impropers.append(
                        Improper(self.atoms[i-1], self.atoms[j-1],
                                 self.atoms[k-1], self.atoms[l-1])
                )
            # Now handle the donors (what is this used for??)
            ndon = self._convert(psfsections['NDON'][0], int, 'number of donors')
            if len(psfsections['NDON'][1]) != ndon * 2:
                raise CharmmError('Got %d indexes for %d donors' % # pragma: no cover
                                  (len(psfsections['NDON'][1]), ndon))
            it = iter(psfsections['NDON'][1])
            for i, j in zip(it, it):
                self.donors.append(
                        AcceptorDonor(self.atoms[i-1], self.atoms[j-1])
                )
            # Now handle the acceptors (what is this used for??)
            nacc = self._convert(psfsections['NACC'][0], int, 'number of acceptors')
            if len(psfsections['NACC'][1]) != nacc * 2:
                raise CharmmError('Got %d indexes for %d acceptors' % # pragma: no cover
                                  (len(psfsections['NACC'][1]), nacc))
            it = iter(psfsections['NACC'][1])
            for i, j in zip(it, it):
                self.acceptors.append(
                        AcceptorDonor(self.atoms[i-1], self.atoms[j-1])
                )
            # Now get the group sections
            try:
                ngrp, nst2 = psfsections['NGRP NST2'][0]
            except ValueError: # pragma: no cover
                raise CharmmError('Could not unpack GROUP pointers') # pragma: no cover
            tmp = psfsections['NGRP NST2'][1]
            self.groups.nst2 = nst2
            # Now handle the groups
            if len(psfsections['NGRP NST2'][1]) != ngrp * 3:
                raise CharmmError('Got %d indexes for %d groups' % # pragma: no cover
                                     (len(tmp), ngrp))
            it = iter(psfsections['NGRP NST2'][1])
            for i, j, k in zip(it, it, it):
                self.groups.append(Group(self.atoms[i], j, k))
            # Assign all of the atoms to molecules recursively
            tmp = psfsections['MOLNT'][1]
            set_molecules(self.atoms)
            molecule_list = [a.marked for a in self.atoms]
            if len(tmp) == len(self.atoms):
                if molecule_list != tmp:
                    warnings.warn('Detected PSF molecule section that is WRONG. '
                                  'Resetting molecularity.', CharmmWarning)
                # We have a CHARMM PSF file; now do NUMLP/NUMLPH sections
                numlp, numlph = psfsections['NUMLP NUMLPH'][0]
                if numlp != 0 or numlph != 0:
                    raise NotImplementedError('Cannot currently handle PSFs with '
                                              'lone pairs defined in the NUMLP/'
                                              'NUMLPH section.')
            # Now do the CMAPs
            ncrterm = self._convert(psfsections['NCRTERM'][0], int, 'Number of cross-terms')
            if len(psfsections['NCRTERM'][1]) != ncrterm * 8:
                raise CharmmError('Got %d CMAP indexes for %d cmap terms' % # pragma: no cover
                                  (len(psfsections['NCRTERM']), ncrterm))
            it = iter(psfsections['NCRTERM'][1])
            for i, j, k, l, m, n, o, p in zip(it, it, it, it, it, it, it, it):
                self.cmaps.append(
                        Cmap.extended(self.atoms[i-1], self.atoms[j-1],
                                      self.atoms[k-1], self.atoms[l-1],
                                      self.atoms[m-1], self.atoms[n-1],
                                      self.atoms[o-1], self.atoms[p-1])
                )
            self.unchange()
            self.flags = psf_flags
        finally:
            if own_handle:
                fileobj.close()

    #===================================================

    @classmethod
    def from_structure(cls, struct, copy=False):
        """
        Instantiates a CharmmPsfFile from an input Structure instance. This
        method makes sure all atom types have uppercase-only names

        Parameters
        ----------
        struct : :class:`parmed.structure.Structure`
            The input structure to convert to a CharmmPsfFile instance
        copy : bool, optional
            If True, a copy of all items are made. Otherwise, the resulting
            CharmmPsfFile is a shallow copy

        Returns
        -------
        psf : :class:`CharmmPsfFile`
            CHARMM PSF file

        Raises
        ------
        ValueError if the functional form is not recognized or cannot be
        implemented through the PSF and parameter/stream files

        Notes
        -----
        If copy is False, the original object may have its atom type names
        changed if any of them have lower-case letters
        """
        from parmed.charmm.parameters import _typeconv as typeconv
        if (struct.rb_torsions or struct.trigonal_angles or
                struct.out_of_plane_bends or struct.pi_torsions or
                struct.stretch_bends or struct.torsion_torsions or
                struct.chiral_frames or struct.multipole_frames or
                struct.nrexcl != 3):
            raise ValueError('Unsupported functional form for CHARMM PSF')
        if copy:
            struct = _copy(struct)
        psf = cls()
        psf.atoms = struct.atoms
        psf.residues = struct.residues
        psf.bonds = struct.bonds
        psf.angles = struct.angles
        psf.urey_bradleys = struct.urey_bradleys
        psf.dihedrals = struct.dihedrals
        psf.impropers = struct.impropers
        psf.acceptors = struct.acceptors
        psf.donors = struct.donors
        psf.groups = struct.groups
        psf.cmaps = struct.cmaps

        psf.bond_types = struct.bond_types
        psf.angle_types = struct.angle_types
        psf.urey_bradley_types = struct.urey_bradley_types
        psf.dihedral_types = struct.dihedral_types
        psf.improper_types = struct.improper_types
        psf.cmap_types = struct.cmap_types

        for atom in psf.atoms:
            atom.type = typeconv(atom.type)
            if atom.atom_type is not UnassignedAtomType:
                atom.atom_type.name = typeconv(atom.atom_type.name)

        # If no groups are defined, make each residue its own group
        if not psf.groups:
            for residue in psf.residues:
                chg = sum(a.charge for a in residue)
                if chg < 1e-4:
                    psf.groups.append(Group(residue[0], 1, 0))
                else:
                    psf.groups.append(Group(residue[0], 2, 0))
            psf.groups.nst2 = 0

        return psf

    #===================================================

    def __str__(self):
        return self.name

    #===================================================

    @needs_openmm
    def createSystem(self, params=None, *args, **kwargs):
        """
        Creates an OpenMM System object from the CHARMM PSF file. This is a
        shortcut for calling `load_parameters` followed by
        Structure.createSystem. If params is not None, `load_parameters` will be
        called on that parameter set, and Structure.createSystem will be called
        with the remaining args and kwargs

        Parameters
        ----------
        params : CharmmParameterSet=None
            If not None, this parameter set will be loaded

        See Also
        --------
        :meth:`parmed.structure.Structure.createSystem`
            In addition to `params`, this method also takes all arguments for
            :meth:`parmed.structure.Structure.createSystem`
        """
        if params is not None: self.load_parameters(params)
        return super(CharmmPsfFile, self).createSystem(*args, **kwargs)

    #===================================================

    def load_parameters(self, parmset, copy_parameters=True):
        """
        Loads parameters from a parameter set that was loaded via CHARMM RTF,
        PAR, and STR files.

        Parameters
        ----------
        parmset : :class:`CharmmParameterSet`
            List of all parameters

        copy_parameters : bool, optional, default=True
            If False, parmset will not be copied.

            WARNING:
            -------
            Not copying parmset will cause ParameterSet and Structure to share
            references to types.  If you modify the original parameter set, the
            references in Structure list_types will be silently modified.
            However, if you change any reference in the parameter set, then that
            reference will no longer be shared with structure.

            Example where the reference in ParameterSet is changed. The
            following will NOT modify the parameters in the psf::

                psf.load_parameters(parmset, copy_parameters=False)
                parmset.angle_types[('a1', 'a2', a3')] = AngleType(1, 2)

            The following WILL change the parameter in the psf because the
            reference has not been changed in ``ParameterSet``::

                psf.load_parameters(parmset, copy_parameters=False)
                a = parmset.angle_types[('a1', 'a2', 'a3')]
                a.k = 10
                a.theteq = 100

            Extra care should be taken when trying this with dihedral_types.
            Since dihedral_type is a Fourier sequence, ParameterSet stores
            DihedralType for every term in DihedralTypeList. Therefore, the
            example below will STILL modify the type in the :class:`Structure`
            list_types::

                parmset.dihedral_types[('a', 'b', 'c', 'd')][0] = DihedralType(1, 2, 3)

            This assigns a new instance of DihedralType to an existing
            DihedralTypeList that ParameterSet and Structure are tracking and
            the shared reference is NOT changed.

            Use with caution!

        Notes
        -----
        - If any dihedral or improper parameters cannot be found, I will try
          inserting wildcards (at either end for dihedrals and as the two
          central atoms in impropers) and see if that matches.  Wild-cards will
          apply ONLY if specific parameters cannot be found.

        - This method will expand the dihedrals attribute by adding a separate
          Dihedral object for each term for types that have a multi-term
          expansion

        Raises
        ------
        ParameterError if any parameters cannot be found
        """
        if copy_parameters:
            parmset = _copy(parmset)
        self.combining_rule = parmset.combining_rule
        # First load the atom types
        for atom in self.atoms:
            try:
                if isinstance(atom.type, int):
                    atype = parmset.atom_types_int[atom.type]
                else:
                    atype = parmset.atom_types_str[atom.type.upper()]
            except KeyError:
                raise ParameterError('Could not find atom type for %s' % atom.type)
            atom.atom_type = atype
            # Change to string type to look up the rest of the parameters. Use
            # upper-case since all parameter sets were read in as upper-case
            atom.type = str(atom.atom_type).upper()
            atom.atomic_number = atype.atomic_number

        # Next load all of the bonds
        for bond in self.bonds:
            # Construct the key
            key = (min(bond.atom1.type, bond.atom2.type),
                   max(bond.atom1.type, bond.atom2.type))
            try:
                bond.type = parmset.bond_types[key]
            except KeyError:
                raise ParameterError('Missing bond type for %r' % bond)
            bond.type.used = False
        # Build the bond_types list
        del self.bond_types[:]
        for bond in self.bonds:
            if bond.type.used: continue
            bond.type.used = True
            self.bond_types.append(bond.type)
            bond.type.list = self.bond_types
        # Next load all of the angles. If a Urey-Bradley term is defined for
        # this angle, also build the urey_bradley and urey_bradley_type lists
        del self.urey_bradleys[:]
        for ang in self.angles:
            # Construct the key
            key = (min(ang.atom1.type, ang.atom3.type), ang.atom2.type,
                   max(ang.atom1.type, ang.atom3.type))
            try:
                ang.type = parmset.angle_types[key]
                ang.type.used = False
                ubt = parmset.urey_bradley_types[key]
                if ubt is not NoUreyBradley:
                    ub = UreyBradley(ang.atom1, ang.atom3, ubt)
                    self.urey_bradleys.append(ub)
                    ubt.used = False
            except KeyError:
                raise ParameterError('Missing angle type for %r' % ang)
        del self.urey_bradley_types[:]
        del self.angle_types[:]
        for ub in self.urey_bradleys:
            if ub.type.used: continue
            ub.type.used = True
            self.urey_bradley_types.append(ub.type)
            ub.type.list = self.urey_bradley_types
        for ang in self.angles:
            if ang.type.used: continue
            ang.type.used = True
            self.angle_types.append(ang.type)
            ang.type.list = self.angle_types
        # Next load all of the dihedrals.
        active_dih_list = set()
        for dih in self.dihedrals:
            # Store the atoms
            a1, a2, a3, a4 = dih.atom1, dih.atom2, dih.atom3, dih.atom4
            key = (a1.type, a2.type, a3.type, a4.type)
            # First see if the exact dihedral is specified
            if not key in parmset.dihedral_types:
                # Check for wild-cards
                key = ('X', a2.type, a3.type, 'X')
                if not key in parmset.dihedral_types:
                    raise ParameterError('No dihedral parameters found for %r' % dih)
            dih.type = parmset.dihedral_types[key]
            dih.type.used = False
            pair = (dih.atom1.idx, dih.atom4.idx) # To determine exclusions
            if (dih.atom1 in dih.atom4.bond_partners or
                dih.atom1 in dih.atom4.angle_partners):
                dih.ignore_end = True
            elif pair in active_dih_list:
                dih.ignore_end = True
            else:
                active_dih_list.add(pair)
                active_dih_list.add((dih.atom4.idx, dih.atom1.idx))
        del self.dihedral_types[:]
        for dihedral in self.dihedrals:
            if dihedral.type.used: continue
            dihedral.type.used = True
            self.dihedral_types.append(dihedral.type)
            dihedral.type.list = self.dihedral_types
        # Now do the impropers
        for imp in self.impropers:
            a1, a2, a3, a4 = imp.atom1.type, imp.atom2.type, imp.atom3.type, imp.atom4.type
            imp.type = parmset.match_improper_type(a1, a2, a3, a4)
            if imp.type is None:
                raise ParameterError('No improper type for %s, %s, %s, and %s', a1, a2, a3, a4)
            imp.type.used = False
        # prepare list of harmonic impropers present in system
        del self.improper_types[:]
        for improper in self.impropers:
            if improper.type.used: continue
            improper.type.used = True
            if isinstance(improper.type, ImproperType):
                self.improper_types.append(improper.type)
                improper.type.list = self.improper_types
            elif isinstance(improper.type, DihedralType):
                self.dihedral_types.append(improper.type)
                improper.type.list = self.dihedral_types
            else:
                assert False, 'Should not be here'
        # Look through the list of impropers -- if there are any periodic
        # impropers, move them over to the dihedrals list
        for i in reversed(range(len(self.impropers))):
            if isinstance(self.impropers[i].type, DihedralType):
                imp = self.impropers.pop(i)
                dih = Dihedral(imp.atom1, imp.atom2, imp.atom3, imp.atom4,
                               improper=True, ignore_end=True, type=imp.type)
                imp.delete()
                self.dihedrals.append(dih)
        # Now do the cmaps. These will not have wild-cards
        for cmap in self.cmaps:
            key = (cmap.atom1.type, cmap.atom2.type, cmap.atom3.type,
                   cmap.atom4.type, cmap.atom2.type, cmap.atom3.type,
                   cmap.atom4.type, cmap.atom5.type)
            try:
                cmap.type = parmset.cmap_types[key]
            except KeyError:
                raise ParameterError('No CMAP parameters found for %r' % cmap)
            cmap.type.used = False
        del self.cmap_types[:]
        for cmap in self.cmaps:
            if cmap.type.used: continue
            cmap.type.used = True
            self.cmap_types.append(cmap.type)
            cmap.type.list = self.cmap_types

    #===================================================

    def clear_cmap(self):
        " Clear the cmap list to prevent any CMAP parameters from being used "
        del self.cmaps[:]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def set_molecules(atoms):
    """
    Correctly sets the molecularity of the system based on connectivity.
    """
    from sys import setrecursionlimit, getrecursionlimit
    # Since we use a recursive function here, we make sure that the recursion
    # limit is large enough to handle the maximum possible recursion depth we'll
    # need (NATOM). We don't want to shrink it, though, since we use list
    # comprehensions in list constructors in some places that have an implicit
    # (shallow) recursion, therefore, reducing the recursion limit too much here
    # could raise a recursion depth exceeded exception during a _Type/Atom/XList
    # creation. Therefore, set the recursion limit to the greater of the current
    # limit or the number of atoms
    setrecursionlimit(max(len(atoms), getrecursionlimit()))

    # Unmark all atoms so we can track which molecule each goes into
    atoms.unmark()

    # The molecule "ownership" list
    owner = []
    # The way I do this is via a recursive algorithm, in which
    # the "set_owner" method is called for each bonded partner an atom
    # has, which in turn calls set_owner for each of its partners and
    # so on until everything has been assigned.
    molecule_number = 1 # which molecule number we are on
    for i in range(len(atoms)):
        # If this atom has not yet been "owned", make it the next molecule
        # However, we only increment which molecule number we're on if
        # we actually assigned a new molecule (obviously)
        if not atoms[i].marked:
            tmp = [i]
            _set_owner(atoms, tmp, i, molecule_number)
            # Make sure the atom indexes are sorted
            tmp.sort()
            owner.append(tmp)
            molecule_number += 1
    return owner

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _set_owner(atoms, owner_array, atm, mol_id):
    """ Recursively sets ownership of given atom and all bonded partners """
    atoms[atm].marked = mol_id
    for partner in atoms[atm].bond_partners:
        if not partner.marked:
            owner_array.append(partner.idx)
            _set_owner(atoms, owner_array, partner.idx, mol_id)
        assert partner.marked == mol_id, 'Atom in multiple molecules!'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
