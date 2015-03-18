"""
Provides a Python class for parsing a PSF file and setting up a system
structure for it

Author: Jason M. Swails
Contributors:
Date: April 20, 2014
"""
from __future__ import division

from compat24 import wraps
from chemistry import (Bond, Angle, Dihedral, Improper, AcceptorDonor, Group,
                       Cmap, UreyBradley, NoUreyBradley, Structure, Atom)
from chemistry.exceptions import (CharmmPSFError, MoleculeError,
                CharmmPSFWarning, MissingParameter, CharmmPsfEOF)
from chemistry.structure import needs_openmm, app, mm
from chemistry import unit as u
from math import sqrt
import os
import warnings

def _catchindexerror(func):
    """
    Protects a function from raising an index error, and replace that exception
    with a CharmmPSFError instead
    """
    @wraps(func)
    def newfunc(*args, **kwargs):
        """ Catch the index error """
        try:
            return func(*args, **kwargs)
        except IndexError, e:
            raise CharmmPSFError('Array is too short: %s' % e)

    return newfunc

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

class CharmmPsfFile(Structure):
    """
    A chemical :class:`Structure` instantiated from CHARMM files.

    Parameters
    ----------
    psf_name : str
        Name of the PSF file (it must exist)

    Raises
    ------
    IOError : If file ``psf_name`` does not exist
    CharmmPsfError : If any parsing errors are encountered

    Examples
    --------
    >>> cs = CharmmPsfFile("testfiles/test.psf")
    >>> len(cs.atoms)
    33
    >>> len(cs.bonds)
    32
    """
    @staticmethod
    def _convert(string, type, message):
        """
        Converts a string to a specific type, making sure to raise
        CharmmPSFError with the given message in the event of a failure.

        Parameters:
            - string (str) : Input string to process
            - type (type) : Type of data to convert to
            - message (str) : Error message to put in exception if failed
        """
        try:
            return type(string)
        except ValueError:
            raise CharmmPSFError('Could not convert %s [%s]' % (message,string))

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
                raise CharmmPsfEOF('Unexpected EOF in PSF file')
            else:
                line = psf.readline()
        if '!' in line:
            words = line[:line.index('!')].split()
            title = line[line.index('!')+1:].strip().upper()
            # Strip out description
            if ':' in title:
                title = title[:title.index(':')]
        else:
            raise CharmmPSFError('Could not determine section title')
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
    def __init__(self, psf_name):
        """
        Opens and parses a PSF file, then instantiates a CharmmPsfFile
        instance from the data.
        """
        Structure.__init__(self)
        conv = CharmmPsfFile._convert
        # Make sure the file exists
        if not os.path.exists(psf_name):
            raise IOError('Could not find PSF file %s' % psf_name)
        # Open the PSF and read the first line. It must start with "PSF"
        psf = open(psf_name, 'r')
        self.name = psf_name
        line = psf.readline()
        if not line.startswith('PSF'):
            raise CharmmPSFError('Unrecognized PSF file. First line is %s' %
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
            except CharmmPsfEOF:
                break
            psfsections[sec] = (ptr, data)
        # store the title
        self.title = psfsections['NTITLE'][1]
        # Next is the number of atoms
        natom = conv(psfsections['NATOM'][0], int, 'natom')
        # Parse all of the atoms
        for i in xrange(natom):
            words = psfsections['NATOM'][1][i].split()
            atid = int(words[0])
            if atid != i + 1:
                raise CharmmPSFError('Nonsequential atoms detected!')
            segid = words[1]
            resid = conv(words[2], int, 'residue number')
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
            self.add_atom(atom, resname, resid, chain=segid)
        # Now get the number of bonds
        nbond = conv(psfsections['NBOND'][0], int, 'number of bonds')
        tmp = psfsections['NBOND'][1]
        if len(tmp) != nbond * 2:
            raise CharmmPSFError('Got %d indexes for %d bonds' %
                                 (len(tmp), nbond))
        for i in xrange(0, 2*nbond, 2):
            self.bonds.append(
                    Bond(self.atoms[tmp[i]-1], self.atoms[tmp[i+1]-1])
            )
        # Now get the number of angles and the angle list
        ntheta = conv(psfsections['NTHETA'][0], int, 'number of angles')
        tmp = psfsections['NTHETA'][1]
        if len(tmp) != ntheta * 3:
            raise CharmmPSFError('Got %d indexes for %d angles' %
                                 (len(tmp), ntheta))
        for i in xrange(0, 3*ntheta, 3):
            self.angles.append(
                    Angle(self.atoms[tmp[i]-1], self.atoms[tmp[i+1]-1],
                          self.atoms[tmp[i+2]-1])
            )
        # Now get the number of torsions and the torsion list
        nphi = conv(psfsections['NPHI'][0], int, 'number of torsions')
        tmp = psfsections['NPHI'][1]
        if len(tmp) != nphi * 4:
            raise CharmmPSFError('Got %d indexes for %d torsions' %
                                 (len(tmp), nphi))
        for i in xrange(0, 4*nphi, 4):
            self.dihedrals.append(
                    Dihedral(self.atoms[tmp[i  ]-1], self.atoms[tmp[i+1]-1],
                             self.atoms[tmp[i+2]-1], self.atoms[tmp[i+3]-1])
            )
        self.dihedrals.split = False
        # Now get the number of improper torsions
        nimphi = conv(psfsections['NIMPHI'][0], int, 'number of impropers')
        tmp = psfsections['NIMPHI'][1]
        if len(tmp) != nimphi * 4:
            raise CharmmPSFError('Got %d indexes for %d impropers' %
                                 (len(tmp), nimphi))
        for i in xrange(0, 4*nimphi, 4):
            self.impropers.append(
                    Improper(self.atoms[tmp[i  ]-1], self.atoms[tmp[i+1]-1],
                             self.atoms[tmp[i+2]-1], self.atoms[tmp[i+3]-1])
            )
        # Now handle the donors (what is this used for??)
        ndon = conv(psfsections['NDON'][0], int, 'number of donors')
        tmp = psfsections['NDON'][1]
        if len(tmp) != ndon * 2:
            raise CharmmPSFError('Got %d indexes for %d donors' %
                                 (len(tmp), ndon))
        for i in xrange(0, 2*ndon, 2):
            self.donors.append(
                    AcceptorDonor(self.atoms[tmp[i]-1], self.atoms[tmp[i+1]-1])
            )
        # Now handle the acceptors (what is this used for??)
        nacc = conv(psfsections['NACC'][0], int, 'number of acceptors')
        tmp = psfsections['NACC'][1]
        if len(tmp) != nacc * 2:
            raise CharmmPSFError('Got %d indexes for %d acceptors' %
                                 (len(tmp), ndon))
        for i in xrange(0, 2*nacc, 2):
            self.acceptors.append(
                    AcceptorDonor(self.atoms[tmp[i]-1], self.atoms[tmp[i+1]-1])
            )
        # Now get the group sections
        try:
            ngrp, nst2 = psfsections['NGRP NST2'][0]
        except ValueError:
            raise CharmmPSFError('Could not unpack GROUP pointers')
        tmp = psfsections['NGRP NST2'][1]
        self.groups.nst2 = nst2
        # Now handle the groups
        if len(tmp) != ngrp * 3:
            raise CharmmPSFError('Got %d indexes for %d groups' %
                                 (len(tmp), ngrp))
        for i in xrange(0, 3*ngrp, 3):
            self.groups.append(Group(tmp[i], tmp[i+1], tmp[i+2]))
        # Assign all of the atoms to molecules recursively
        tmp = psfsections['MOLNT'][1]
        set_molecules(self.atoms)
        molecule_list = [a.marked for a in self.atoms]
        if len(tmp) == len(self.atoms):
            if molecule_list != tmp:
                warnings.warn('Detected PSF molecule section that is WRONG. '
                              'Resetting molecularity.', CharmmPSFWarning)
            # We have a CHARMM PSF file; now do NUMLP/NUMLPH sections
            numlp, numlph = psfsections['NUMLP NUMLPH'][0]
            if numlp != 0 or numlph != 0:
                raise NotImplementedError('Cannot currently handle PSFs with '
                                          'lone pairs defined in the NUMLP/'
                                          'NUMLPH section.')
        # Now do the CMAPs
        ncrterm = conv(psfsections['NCRTERM'][0], int, 'Number of cross-terms')
        tmp = psfsections['NCRTERM'][1]
        if len(tmp) != ncrterm * 8:
            raise CharmmPSFError('Got %d CMAP indexes for %d cmap terms' %
                                 (len(tmp), ncrterm))
        for i in xrange(0, 8*ncrterm, 8):
            self.cmaps.append(
                    Cmap.extended(self.atoms[tmp[i  ]-1],
                                  self.atoms[tmp[i+1]-1],
                                  self.atoms[tmp[i+2]-1],
                                  self.atoms[tmp[i+3]-1],
                                  self.atoms[tmp[i+4]-1],
                                  self.atoms[tmp[i+5]-1],
                                  self.atoms[tmp[i+6]-1],
                                  self.atoms[tmp[i+7]-1])
            )
        self.unchange()
        self.flags = psf_flags

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
        :meth:`chemistry.structure.Structure.createSystem`
            In addition to `params`, this method also takes all arguments for
            :meth:`chemistry.structure.Structure.createSystem`
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
        - If any parameters that are necessary cannot be found, MissingParameter
          is raised.

        - If any dihedral or improper parameters cannot be found, I will try
          inserting wildcards (at either end for dihedrals and as the two
          central atoms in impropers) and see if that matches.  Wild-cards will
          apply ONLY if specific parameters cannot be found.

        - This method will expand the dihedrals attribute by adding a separate
          Dihedral object for each term for types that have a multi-term
          expansion
        """
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
                raise MissingParameter('Could not find atom type for %s' %
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
                raise MissingParameter('Missing bond type for %r' % bond)
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
                raise MissingParameter('Missing angle type for %r' % ang)
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
                    raise MissingParameter('No dihedral parameters found for '
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
                raise MissingParameter('No improper parameters found for %r' %
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
                raise MissingParameter('No CMAP parameters found for %r' % cmap)
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

    def has_NBFIX(self):
        """
        Returns whether or not any pairs of atom types have their LJ
        interactions modified by an NBFIX definition

        Returns
        -------
        has_nbfix : bool
            If True, at least two atom types have NBFIX mod definitions
        """
        typemap = dict()
        for a in self.atoms:
            typemap[str(a.atom_type)] = a.atom_type
        # Now we have a map of all atom types that we have defined in our
        # system. Look through all of the atom types and see if any of their
        # NBFIX definitions are also keys in typemap
        for key, type in typemap.iteritems():
            for key in type.nbfix:
                if key in typemap:
                    return True
        return False

    #===================================================

    @needs_openmm
    def omm_nonbonded_force(self, nonbondedMethod=None,
                            nonbondedCutoff=8*u.angstroms,
                            switchDistance=0*u.angstroms,
                            ewaldErrorTolerance=0.0005,
                            reactionFieldDielectric=78.5):
        """ Creates the OpenMM NonbondedForce instance

        Parameters
        ----------
        nonbondedMethod : cutoff method
            This is the cutoff method. It can be either the NoCutoff,
            CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald objects from the
            simtk.openmm.app namespace
        nonbondedCutoff : float or distance Quantity
            The nonbonded cutoff must be either a floating point number
            (interpreted as nanometers) or a Quantity with attached units. This
            is ignored if nonbondedMethod is NoCutoff.
        switchDistance : float or distance Quantity
            The distance at which the switching function is turned on for van
            der Waals interactions. This is ignored when no cutoff is used, and
            no switch is used if switchDistance is 0, negative, or greater than
            the cutoff
        ewaldErrorTolerance : float=0.0005
            When using PME or Ewald, the Ewald parameters will be calculated
            from this value
        reactionFieldDielectric : float=78.5
            If the nonbondedMethod is CutoffPeriodic or CutoffNonPeriodic, the
            region beyond the cutoff is treated using a reaction field method
            with this dielectric constant. It should be set to 1 if another
            implicit solvent model is being used (e.g., GB)

        Returns
        -------
        NonbondedForce
            This just implements the very basic NonbondedForce with the typical
            charge-charge and 12-6 Lennard-Jones interactions with the
            Lorentz-Berthelot combining rules.

        Notes
        -----
        Subclasses of Structure for which this nonbonded treatment is inadequate
        should override this method to implement what is needed
        """
        if not self.atoms: return None
        nonbfrc = super(CharmmPsfFile, self).omm_nonbonded_force(
                nonbondedMethod, nonbondedCutoff, switchDistance,
                ewaldErrorTolerance, reactionFieldDielectric,
        )
        hasnbfix = self.has_NBFIX()
        if not hasnbfix:
            return nonbfrc
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        ene_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        # We need a CustomNonbondedForce to implement the NBFIX functionality.
        # First derive the type lookup tables
        lj_idx_list = [0 for atom in self.atoms]
        lj_radii, lj_depths = [], []
        num_lj_types = 0
        lj_type_list = []
        for i, atom in enumerate(self.atoms):
            atype = atom.atom_type
            if lj_idx_list[i]: continue # already assigned
            num_lj_types += 1
            lj_idx_list[i] = num_lj_types
            ljtype = (atype.rmin, abs(atype.epsilon))
            lj_type_list.append(atype)
            lj_radii.append(atype.rmin)
            lj_depths.append(abs(atype.epsilon))
            for j in xrange(i+1, len(self.atoms)):
                atype2 = self.atoms[j].atom_type
                if lj_idx_list[j] > 0: continue # already assigned
                if atype2 is atype:
                    lj_idx_list[j] = num_lj_types
                elif not atype.nbfix:
                    # Only non-NBFIXed atom types can be compressed
                    ljtype2 = (atype2.rmin, abs(atype2.epsilon))
                    if ljtype == ljtype2:
                        lj_idx_list[j] = num_lj_types
        # Now everything is assigned. Create the A- and B-coefficient arrays
        acoef = [0 for i in xrange(num_lj_types*num_lj_types)]
        bcoef = acoef[:]
        for i in xrange(num_lj_types):
            for j in xrange(num_lj_types):
                namej = lj_type_list[j].name
                try:
                    rij, wdij, rij14, wdij14 = lj_type_list[i].nbfix[namej]
                except KeyError:
                    rij = (lj_radii[i] + lj_radii[j]) * length_conv
                    wdij = sqrt(lj_depths[i] * lj_depths[j]) * ene_conv
                else:
                    rij *= length_conv
                    wdij *= ene_conv
                rij6 = rij**6
                acoef[i+num_lj_types*j] = sqrt(wdij) * rij6
                bcoef[i+num_lj_types*j] = 2 * wdij * rij6
        force = mm.CustomNonbondedForce('(a/r6)^2-b/r6; r6=r2*r2*r2; r2=r^2; '
                                        'a=acoef(type1, type2); '
                                        'b=bcoef(type1, type2)')
        force.addTabulatedFunction('acoef',
                mm.Discrete2DFunction(num_lj_types, num_lj_types, acoef))
        force.addTabulatedFunction('bcoef',
                mm.Discrete2DFunction(num_lj_types, num_lj_types, bcoef))
        force.addPerParticleParameter('type')
        force.setForceGroup(self.NONBONDED_FORCE_GROUP)
        if (nonbondedMethod is app.PME or nonbondedMethod is app.Ewald or
                nonbondedMethod is app.CutoffPeriodic):
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        elif nonbondedMethod is app.NoCutoff:
            force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
        elif nonbondedMethod is app.CutoffNonPeriodic:
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        else:
            raise ValueError('Unrecognized nonbonded method [%s]' %
                             nonbondedMethod)
        # Add the particles
        for i in lj_idx_list:
            force.addParticle((i-1,))
        # Now wipe out the L-J parameters in the nonbonded force
        for i in xrange(nonbfrc.getNumParticles()):
            chg, sig, eps = nonbfrc.getParticleParameters(i)
            nonbfrc.setParticleParameters(i, chg, 0.5, 0.0)
        # Now transfer the exclusions
        for ii in xrange(nonbfrc.getNumExceptions()):
            i, j, qq, ss, ee = nonbfrc.getExceptionParameters(ii)
            force.addExclusion(i, j)
        # Now transfer the other properties (cutoff, switching function, etc.)
        force.setUseLongRangeCorrection(True)
        if nonbondedMethod is app.NoCutoff:
            force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
        elif nonbondedMethod is app.CutoffNonPeriodic:
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        elif nonbondedMethod in (app.PME, app.Ewald, app.CutoffPeriodic):
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        else:
            raise ValueError('Unsupported nonbonded method %s' %
                             nonbondedMethod)
        force.setCutoffDistance(nonbfrc.getCutoffDistance())
        if nonbfrc.getUseSwitchingFunction():
            force.setUseSwitchingFunction(True)
            force.setSwitchingDistance(nonbfrc.getSwitchingDistance())

        return nonbfrc, force

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
    for i in xrange(len(atoms)):
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

if __name__ == '__main__':
    import doctest
    doctest.testmod()
