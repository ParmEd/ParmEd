"""
Provides a Python class for parsing a PSF file and setting up a system
structure for it

Author: Jason M. Swails
Contributors:
Date: April 20, 2014
"""
from __future__ import division

from contextlib import closing
from copy import copy as _copy
from math import sqrt
from parmed import (Bond, Angle, Dihedral, Improper, AcceptorDonor, Group,
                    Cmap, UreyBradley, NoUreyBradley, Structure, Atom,
                    DihedralType, AngleType, ExtraPoint, DihedralTypeList)
from parmed.constants import SMALL
from parmed.exceptions import (CharmmError, MoleculeError, CharmmWarning,
        ParameterError)
from parmed.structure import needs_openmm
from parmed.utils.io import genopen
from parmed.utils.six import wraps
from parmed.utils.six.moves import zip, range
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

_resre = re.compile(r'(\d+)([a-zA-Z]*)')

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

    @staticmethod
    def _parse_psf_section(psf):
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
        conv = CharmmPsfFile._convert
        line = psf.readline()
        while not line.strip():
            if not line:
                raise _FileEOF('Unexpected EOF in PSF file')
            else:
                line = psf.readline()
        if '!' in line:
            words = line[:line.index('!')].split()
            title = line[line.index('!')+1:].strip().upper()
            # Strip out description
            if ':' in title:
                title = title[:title.index(':')]
        else:
            raise CharmmError('Could not determine section title')
        if len(words) == 1:
            pointers = conv(words[0], int, 'pointer')
        else:
            pointers = tuple([conv(w, int, 'pointer') for w in words])
        line = psf.readline().strip()
        if not line and title.startswith('NNB'):
            # This will correctly handle the NNB section (which has a spurious
            # blank line) as well as any sections that have 0 members.
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
                data.extend([conv(w, int, 'PSF data') for w in words])
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
        conv = CharmmPsfFile._convert
        # Open the PSF and read the first line. It must start with "PSF"
        with closing(genopen(psf_name, 'r')) as psf:
            self.name = psf_name
            line = psf.readline()
            if not line.startswith('PSF'):
                raise CharmmError('Unrecognized PSF file. First line is %s' %
                                     line.strip())
            # Store the flags
            psf_flags = line.split()[1:]
            # Now get all of the sections and store them in a dict
            psf.readline()
            # Now get all of the sections
            psfsections = _ZeroDict()
            while True:
                try:
                    sec, ptr, data = CharmmPsfFile._parse_psf_section(psf)
                except _FileEOF:
                    break
                psfsections[sec] = (ptr, data)
            # store the title
            self.title = psfsections['NTITLE'][1]
            # Next is the number of atoms
            natom = conv(psfsections['NATOM'][0], int, 'natom')
            # Parse all of the atoms
            for i in range(natom):
                words = psfsections['NATOM'][1][i].split()
                atid = int(words[0])
                if atid != i + 1:
                    raise CharmmError('Nonsequential atoms detected!')
                segid = words[1]
                rematch = _resre.match(words[2])
                if not rematch:
                    raise RuntimeError('Could not interpret residue number %s' %
                                       words[2])
                resid, inscode = rematch.groups()
                resid = conv(resid, int, 'residue number')
                resname = words[3]
                name = words[4]
                attype = words[5]
                # Try to convert the atom type to an integer a la CHARMM
                try:
                    attype = int(attype)
                except ValueError:
                    pass
                charge = conv(words[6], float, 'partial charge')
                mass = conv(words[7], float, 'atomic mass')
                props = words[8:]
                atom = Atom(name=name, type=attype, charge=charge, mass=mass)
                atom.segid = segid
                atom.props = props
                self.add_atom(atom,resname,resid,chain=segid,inscode=inscode)
            # Now get the number of bonds
            nbond = conv(psfsections['NBOND'][0], int, 'number of bonds')
            if len(psfsections['NBOND'][1]) != nbond * 2:
                raise CharmmError('Got %d indexes for %d bonds' %
                                     (len(psfsections['NBOND'][1]), nbond))
            it = iter(psfsections['NBOND'][1])
            for i, j in zip(it, it):
                self.bonds.append(Bond(self.atoms[i-1], self.atoms[j-1]))
            # Now get the number of angles and the angle list
            ntheta = conv(psfsections['NTHETA'][0], int, 'number of angles')
            if len(psfsections['NTHETA'][1]) != ntheta * 3:
                raise CharmmError('Got %d indexes for %d angles' %
                                     (len(psfsections['NTHETA'][1]), ntheta))
            it = iter(psfsections['NTHETA'][1])
            for i, j, k in zip(it, it, it):
                self.angles.append(
                        Angle(self.atoms[i-1], self.atoms[j-1], self.atoms[k-1])
                )
                self.angles[-1].funct = 5 # urey-bradley
            # Now get the number of torsions and the torsion list
            nphi = conv(psfsections['NPHI'][0], int, 'number of torsions')
            if len(psfsections['NPHI'][1]) != nphi * 4:
                raise CharmmError('Got %d indexes for %d torsions' %
                                     (len(psfsections['NPHI']), nphi))
            it = iter(psfsections['NPHI'][1])
            for i, j, k, l in zip(it, it, it, it):
                self.dihedrals.append(
                        Dihedral(self.atoms[i-1], self.atoms[j-1],
                                 self.atoms[k-1], self.atoms[l-1])
                )
            self.dihedrals.split = False
            # Now get the number of improper torsions
            nimphi = conv(psfsections['NIMPHI'][0], int, 'number of impropers')
            if len(psfsections['NIMPHI'][1]) != nimphi * 4:
                raise CharmmError('Got %d indexes for %d impropers' %
                                     (len(psfsections['NIMPHI'][1]), nimphi))
            it = iter(psfsections['NIMPHI'][1])
            for i, j, k, l in zip(it, it, it, it):
                self.impropers.append(
                        Improper(self.atoms[i-1], self.atoms[j-1],
                                 self.atoms[k-1], self.atoms[l-1])
                )
            # Now handle the donors (what is this used for??)
            ndon = conv(psfsections['NDON'][0], int, 'number of donors')
            if len(psfsections['NDON'][1]) != ndon * 2:
                raise CharmmError('Got %d indexes for %d donors' %
                                     (len(psfsections['NDON'][1]), ndon))
            it = iter(psfsections['NDON'][1])
            for i, j in zip(it, it):
                self.donors.append(
                        AcceptorDonor(self.atoms[i-1], self.atoms[j-1])
                )
            # Now handle the acceptors (what is this used for??)
            nacc = conv(psfsections['NACC'][0], int, 'number of acceptors')
            if len(psfsections['NACC'][1]) != nacc * 2:
                raise CharmmError('Got %d indexes for %d acceptors' %
                                     (len(psfsections['NACC'][1]), nacc))
            it = iter(psfsections['NACC'][1])
            for i, j in zip(it, it):
                self.acceptors.append(
                        AcceptorDonor(self.atoms[i-1], self.atoms[j-1])
                )
            # Now get the group sections
            try:
                ngrp, nst2 = psfsections['NGRP NST2'][0]
            except ValueError:
                raise CharmmError('Could not unpack GROUP pointers')
            tmp = psfsections['NGRP NST2'][1]
            self.groups.nst2 = nst2
            # Now handle the groups
            if len(psfsections['NGRP NST2'][1]) != ngrp * 3:
                raise CharmmError('Got %d indexes for %d groups' %
                                     (len(tmp), ngrp))
            it = iter(psfsections['NGRP NST2'][1])
            for i, j, k in zip(it, it, it):
                self.groups.append(Group(i, j, k))
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
            ncrterm = conv(psfsections['NCRTERM'][0], int, 'Number of cross-terms')
            if len(psfsections['NCRTERM'][1]) != ncrterm * 8:
                raise CharmmError('Got %d CMAP indexes for %d cmap terms' %
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
        if (struct.rb_torsions or struct.trigonal_angles or
                struct.out_of_plane_bends or struct.pi_torsions or
                struct.stretch_bends or struct.torsion_torsions or
                struct.chiral_frames or struct.multipole_frames):
            raise ValueError('Unsupported functional form for CHARMM PSF')
        if copy:
            struct = _copy(struct)
        psf = cls()
        psf.atoms = struct.atoms
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

        # Make all atom type names upper-case
        def typeconv(name):
            if name.upper() == name:
                return name
            # Lowercase letters present -- decorate the type name with LTU --
            # Lower To Upper
            return '%sLTU' % name.upper()
        for atom in psf.atoms:
            atom.type = typeconv(atom.type)
            if atom.atom_type is not None:
                atom.atom_type.name = typeconv(atom.atom_type.name)

        # In case CHARMM force fields define all of their exclusions from the
        # bond, angle, and dihedral lists, go ahead and fill in any blank,
        # zeroed terms that need to exist to fill them out from the bond graph
        # (the same function is used in the AmberParm class for the same reason)
#       psf._add_missing_13_14()
#       del psf.adjusts[:]

        # In some cases, 1-4 interactions are defined in the dihedral list. If
        # 1-4 van der Waals scalings exist, these are implemented in CHARMM
        # through scaled 1-4 epsilon parameters. So first we update our dihedral
        # exclusions to make sure we skip over torsions that have no effect on
        # the 1-4 pairlist, then walk through the pairlist and assign the
        # correct 1-4 epsilon parameter
        psf.update_dihedral_exclusions()
        for dihedral in psf.dihedrals:
            if dihedral.improper or dihedral.ignore_end: continue
            if dihedral.type is None: continue
            a1, a2 = dihedral.atom1, dihedral.atom4
            if a1.epsilon_14 == a1.epsilon and a2.epsilon_14 == a2.epsilon:
                if isinstance(dihedral.type, DihedralType):
                    scnb = dihedral.type.scnb
                elif isinstance(dihedral.type, DihedralTypeList):
                    scnb = dihedral.type[0].scnb
                elif dihedral.type is None:
                    continue
                else:
                    assert False, "Should not be here"
                if scnb != 1:
                    a1.epsilon_14 = a1.epsilon / scnb
                    a2.epsilon_14 = a2.epsilon / scnb
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

    def load_parameters(self, parmset):
        """
        Loads parameters from a parameter set that was loaded via CHARMM RTF,
        PAR, and STR files.

        Parameters
        ----------
        parmset : :class:`CharmmParameterSet`
            List of all parameters

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
        parmset = _copy(parmset)
        self.combining_rule = parmset.combining_rule
        # First load the atom types
        types_are_int = False
        for atom in self.atoms:
            try:
                if isinstance(atom.type, int):
                    atype = parmset.atom_types_int[atom.type]
                    types_are_int = True # if we have to change back
                else:
                    atype = parmset.atom_types_str[atom.type]
            except KeyError:
                raise ParameterError('Could not find atom type for %s' %
                                     atom.type)
            atom.atom_type = atype
            # Change to string type to look up the rest of the parameters
            atom.type = str(atom.atom_type)
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
                    raise ParameterError('No dihedral parameters found for '
                                           '%r' % dih)
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
            # Store the atoms
            a1, a2, a3, a4 = imp.atom1, imp.atom2, imp.atom3, imp.atom4
            at1, at2, at3, at4 = a1.type, a2.type, a3.type, a4.type
            key = tuple(sorted([at1, at2, at3, at4]))
            if not key in parmset.improper_types:
                # Check for wild-cards
                for anchor in (at2, at3, at4):
                    key = tuple(sorted([at1, anchor, 'X', 'X']))
                    if key in parmset.improper_types:
                        break # This is the right key
            try:
                imp.type = parmset.improper_types[key]
            except KeyError:
                raise ParameterError('No improper parameters found for %r' %
                                       imp)
            imp.type.used = False
        del self.improper_types[:]
        for improper in self.impropers:
            if improper.type.used: continue
            improper.type.used = True
            self.improper_types.append(improper.type)
            improper.type.list = self.improper_types
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
        # If the types started out as integers, change them back
        if types_are_int:
            for atom in self.atoms: atom.type = int(atom.atom_type)

    #===================================================

    def clear_cmap(self):
        " Clear the cmap list to prevent any CMAP parameters from being used "
        del self.cmaps[:]

    #===================================================

    def _add_missing_13_14(self, ignore_inconsistent_vdw=False):
        """
        Uses the bond graph to fill in zero-parameter angles and dihedrals. The
        reason this is necessary is that Amber assumes that the list of angles
        and dihedrals encompasses *all* 1-3 and 1-4 pairs as determined by the
        bond graph, respectively. As a result, Amber programs use the angle and
        dihedral lists to set nonbonded exclusions and exceptions.

        Parameters
        ----------
        ignore_inconsistent_vdw : bool, optional
            If True, do not make inconsistent 1-4 vdW parameters fatal. For
            ChamberParm, the 1-4 specific vdW parameters can compensate. For
            AmberParm, the 1-4 scaling factor cannot represent arbitrary
            exceptions. Default is False (should only be True for ChamberParm)
        """
        # We need to figure out what 1-4 scaling term to use if we don't have
        # explicit exceptions
        if not self.adjusts:
            eel_scale = set()
            for dih in self.dihedrals:
                if dih.ignore_end or dih.improper: continue
                if isinstance(dih.type, DihedralType):
                    eel_scale.add(dih.type.scee)
                elif isinstance(dih.type, DihedralTypeList):
                    eel_scale.add(dih.type[0].scee)
                else:
                    assert False, 'Should not be here'
            if len(eel_scale) > 1:
                raise ValueError('Cannot have mixed 1-4 EEL scaling')
            elif len(eel_scale) == 1:
                scee = list(eel_scale)[0]
            else:
                scee = 1e10
            zero_torsion = DihedralType(0, 1, 0, scee, 1.0)
        else:
            # Turn list of exceptions into a dict so we can look it up quickly
            adjust_dict = dict()
            for pair in self.adjusts:
                adjust_dict[tuple(sorted([pair.atom1, pair.atom2]))] = pair
            ignored_torsion = None
            zero_torsion = None
            # Scan through existing dihedrals to make sure the exceptions match
            # the dihedral list
            if self.combining_rule == 'lorentz':
                comb_sig = lambda sig1, sig2: 0.5 * (sig1 + sig2)
            elif self.combining_rule == 'geometric':
                comb_sig = lambda sig1, sig2: sqrt(sig1 * sig2)
            else:
                assert False, "Unrecognized combining rule"
            fac = 2**(1/6)
            for dihedral in self.dihedrals:
                if dihedral.ignore_end: continue
                key = tuple(sorted([dihedral.atom1, dihedral.atom4]))
                eref = sqrt(dihedral.atom1.epsilon_14*dihedral.atom4.epsilon_14)
                rref = comb_sig(dihedral.atom1.sigma_14,
                                dihedral.atom4.sigma_14) * fac
                if key in adjust_dict:
                    pair = adjust_dict[key]
                    if pair.type.epsilon == 0:
                        scnb = 1e10
                    else:
                        scnb = eref / pair.type.epsilon
                    if pair.type.chgscale == 0:
                        scee = 1e10
                    else:
                        scee = 1 / pair.type.chgscale
                    if (abs(rref - pair.type.rmin) > SMALL and
                            pair.type.epsilon != 0):
                        if ignore_inconsistent_vdw:
                            scnb = 1.0
                        else:
                            raise TypeError('Cannot translate exceptions')
                    if (abs(scnb - dihedral.type.scnb) < SMALL and
                            abs(scee - dihedral.type.scee) < SMALL):
                        continue
                else:
                    scee = scnb = 1e10
                newtype = _copy(dihedral.type)
                newtype.scee = scee
                newtype.scnb = scnb
                dihedral.type = newtype
                newtype.list = self.dihedral_types
                self.dihedral_types.append(newtype)

        zero_angle = AngleType(0, 0)

        n13 = n14 = 0
        if self.combining_rule == 'lorentz':
            comb_sig = lambda sig1, sig2: 0.5 * (sig1 + sig2)
        elif self.combining_rule == 'geometric':
            comb_sig = lambda sig1, sig2: sqrt(sig1 * sig2)
        else:
            assert False, 'Unrecognized combining rule. Should not be here'
        fac = 2**(1/6)
        for atom in self.atoms:
            if isinstance(atom, ExtraPoint): continue
            for batom in atom.bond_partners:
                if isinstance(batom, ExtraPoint): continue
                for aatom in batom.bond_partners:
                    if isinstance(aatom, ExtraPoint) or aatom is atom: continue
                    for datom in aatom.bond_partners:
                        if isinstance(datom, ExtraPoint): continue
                        if (datom in atom.angle_partners + atom.bond_partners +
                                atom.dihedral_partners or datom is atom):
                            continue
                        # Add the missing dihedral
                        if not self.adjusts:
                            tortype = zero_torsion
                            if n14 == 0:
                                tortype.list = self.dihedral_types
                                self.dihedral_types.append(tortype)
                        else:
                            # Figure out what the scale factors must be
                            key = tuple(sorted([atom, datom]))
                            if key not in adjust_dict:
                                if ignored_torsion is None:
                                    ignored_torsion = DihedralType(0, 1, 0,
                                                                   1e10, 1e10)
                                    self.dihedral_types.append(ignored_torsion)
                                    ignored_torsion.list = self.dihedral_types
                                tortype = ignored_torsion
                            elif 0 in (adjust_dict[key].type.epsilon,
                                       adjust_dict[key].type.rmin) and \
                                    adjust_dict[key].type.chgscale == 0:
                                if ignored_torsion is None:
                                    ignored_torsion = \
                                            DihedralType(0, 1, 0, 1e10, 1e10,
                                                    list=self.dihedral_types)
                                    self.dihedral_types.append(ignored_torsion)
                                tortype = ignored_torsion
                            else:
                                pair = adjust_dict[key]
                                epsilon = pair.type.epsilon
                                rmin = pair.type.rmin
                                # Compare it to the 1-4 parameters that are
                                # already present
                                eref = sqrt(pair.atom1.epsilon_14*
                                            pair.atom2.epsilon_14)
                                if pair.type.epsilon == 0:
                                    scnb = 1e10
                                else:
                                    scnb = eref / epsilon
                                if pair.type.chgscale == 0:
                                    scee = 1e10
                                else:
                                    scee = 1 / pair.type.chgscale
                                rref = comb_sig(pair.atom1.sigma_14,
                                                pair.atom2.sigma_14) * fac
                                if abs(rmin - rref) > SMALL:
                                    if ignore_inconsistent_vdw:
                                        scnb = 1.0
                                    else:
                                        raise TypeError(
                                                'Cannot translate exceptions'
                                        )
                                tortype = DihedralType(0, 1, 0, scee, scnb,
                                                       list=self.dihedral_types)
                                self.dihedral_types.append(tortype)
                        dihedral = Dihedral(atom, batom, aatom, datom,
                                            ignore_end=False, improper=False,
                                            type=tortype)
                        self.dihedrals.append(dihedral)
                        n14 += 1
                    if aatom in atom.angle_partners + atom.bond_partners:
                        continue
                    # Add the missing angle
                    self.angles.append(Angle(atom, batom, aatom, zero_angle))
                    n13 += 1

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
        elif partner.marked != mol_id:
            raise MoleculeError('Atom %d in multiple molecules' % 
                                partner.idx)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
