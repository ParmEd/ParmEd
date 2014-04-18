"""
Provides a Python class for parsing a PSF file and setting up a system
structure for it

Author: Jason M. Swails
Contributors:
Date: April 5, 2014
"""
from __future__ import division

from compat24 import wraps
from chemistry.charmm._charmmfile import CharmmFile
from chemistry.charmm.topologyobjects import (ResidueList, AtomList,
                TrackedList, Bond, Angle, Dihedral, Improper, AcceptorDonor,
                Group, Cmap, UreyBradley, NoUreyBradley)
from chemistry.exceptions import (CharmmPSFError, MoleculeError,
                CharmmPSFWarning, MissingParameter)
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CharmmPsfFile(object):
    """
    A chemical structure instantiated from CHARMM files.

    Example:
    >>> cs = CharmmPsfFile("testfiles/test.psf")
    
    This structure has numerous attributes that are lists of the elements of
    this structure, including atoms, bonds, torsions, etc. The attributes are
        - residue_list
        - atom_list
        - bond_list
        - angle_list
        - dihedral_list
        - improper_list
        - cmap_list
        - donor_list    # hbonds donors?
        - acceptor_list # hbond acceptors?
        - group_list    # list of nonbonded interaction groups

    Additional attribute is available if a CharmmParameterSet is loaded into
    this structure.
        
        - urey_bradley_list

    The lengths of each of these lists gives the pointers (e.g., natom, nres,
    etc.)

    Example:
    >>> cs = CharmmPsfFile("testfiles/test.psf")
    >>> len(cs.atom_list)
    33
    >>> len(cs.bond_list)
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
        except ValueError, e:
            print e
            raise CharmmPSFError('Could not convert %s' % message)

    @staticmethod
    def _parse_psf_section(psf, dtype):
        """
        This method parses a section of the PSF file

        Parameters:
            - psf (CharmmFile) : Open file that is pointing to the first line
                                 of the section that is to be parsed
            - dtype (type) : The data type to convert all of the data into
        
        Returns:
            (pointers, data)
            
            - pointers (int/tuple of ints) : If one pointer is set, pointers is
                    simply the integer that is value of that pointer. Otherwise
                    it is a tuple with every pointer value defined in the first
                    line
            - data (list) : A list of all data in the parsed section converted
                    to `dtype'
        """
        conv = CharmmPsfFile._convert
        words = psf.readline().split()
        if len(words) == 1:
            pointers = conv(words[0], int, 'pointer')
        else:
            pointers = tuple([conv(w, int, 'pointer') for w in words])
        line = psf.readline().strip()
        if not line:
            # This will correctly handle the NNB section (which has a spurious
            # blank line) as well as any sections that have 0 members.
            line = psf.readline().strip()
        data = []
        while line:
            words = line.split()
            data.extend([conv(w, dtype, 'PSF data') for w in words])
            line = psf.readline().strip()
        return pointers, data

    @_catchindexerror
    def __init__(self, psf_name):
        """
        Opens and parses a PSF file, then instantiates a CharmmPsfFile
        instance from the data.
            
        Parameters:
            psf_name (str) : Name of the PSF file (it must exist)
        
        Exceptions Raised:
            IOError : If file "psf_name" does not exist
            CharmmPSFError: If any parsing errors are encountered
        """
        conv = CharmmPsfFile._convert
        # Make sure the file exists
        if not os.path.exists(psf_name):
            raise IOError('Could not find PSF file %s' % psf_name)
        # Open the PSF and read the first line. It must start with "PSF"
        psf = CharmmFile(psf_name, 'r')
        line = psf.readline()
        if not line.startswith('PSF'):
            raise CharmmPSFError('Unrecognized PSF file. First line is %s' %
                                 line.strip())
        # Store the flags
        psf_flags = line.split()[1:]
        # Next line is blank
        psf.readline()
        # The next line has one number -- the number of title lines
        ntitle = conv(psf.readline().strip(), int, 'title count')
        # store the title
        title = list()
        for i in range(ntitle):
            title.append(psf.readline().rstrip())
        # Skip the blank line
        psf.readline()
        # Next is the number of atoms
        natom = conv(psf.readline().strip(), int, 'natom')
        # Parse all of the atoms
        residue_list = ResidueList()
        atom_list = AtomList()
        for i in xrange(natom):
            words = psf.readline().split()
            atid = int(words[0])
            if atid != i + 1:
                raise CharmmPSFError('Nonsequential atoms detected!')
            system = words[1]
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
            atom = residue_list.add_atom(system, resid, resname, name,
                                         attype, charge, mass, props)
            atom_list.append(atom)
        atom_list.assign_indexes()
        # Eat the next line
        psf.readline()
        # Now get the number of bonds
        nbond, holder = CharmmPsfFile._parse_psf_section(psf, int)
        bond_list = TrackedList()
        if len(holder) != nbond * 2:
            raise CharmmPSFError('Got %d indexes for %d bonds' %
                                 (len(holder), nbond))
        for i in range(nbond):
            id1 = holder[2*i  ] - 1
            id2 = holder[2*i+1] - 1
            bond_list.append(Bond(atom_list[id1], atom_list[id2]))
        bond_list.changed = False
        # Now get the number of angles and the angle list
        ntheta, holder = CharmmPsfFile._parse_psf_section(psf, int)
        angle_list = TrackedList()
        if len(holder) != ntheta * 3:
            raise CharmmPSFError('Got %d indexes for %d angles' %
                                 (len(holder), ntheta))
        for i in range(ntheta):
            id1 = holder[3*i  ] - 1
            id2 = holder[3*i+1] - 1
            id3 = holder[3*i+2] - 1
            angle_list.append(
                    Angle(atom_list[id1], atom_list[id2], atom_list[id3])
            )
        angle_list.changed = False
        # Now get the number of torsions and the torsion list
        nphi, holder = CharmmPsfFile._parse_psf_section(psf, int)
        dihedral_list = TrackedList()
        if len(holder) != nphi * 4:
            raise CharmmPSFError('Got %d indexes for %d torsions' %
                                 (len(holder), nphi))
        for i in range(nphi):
            id1 = holder[4*i  ] - 1
            id2 = holder[4*i+1] - 1
            id3 = holder[4*i+2] - 1
            id4 = holder[4*i+3] - 1
            dihedral_list.append(
                    Dihedral(atom_list[id1], atom_list[id2], atom_list[id3],
                             atom_list[id4])
            )
        dihedral_list.changed = False
        # Now get the number of improper torsions
        nimphi, holder = CharmmPsfFile._parse_psf_section(psf, int)
        improper_list = TrackedList()
        if len(holder) != nimphi * 4:
            raise CharmmPSFError('Got %d indexes for %d impropers' %
                                 (len(holder), nimphi))
        for i in range(nimphi):
            id1 = holder[4*i  ] - 1
            id2 = holder[4*i+1] - 1
            id3 = holder[4*i+2] - 1
            id4 = holder[4*i+3] - 1
            improper_list.append(
                    Improper(atom_list[id1], atom_list[id2], atom_list[id3],
                             atom_list[id4])
            )
        improper_list.changed = False
        # Now handle the donors (what is this used for??)
        ndon, holder = CharmmPsfFile._parse_psf_section(psf, int)
        donor_list = TrackedList()
        if len(holder) != ndon * 2:
            raise CharmmPSFError('Got %d indexes for %d donors' %
                                 (len(holder), ndon))
        for i in range(ndon):
            id1 = holder[2*i  ] - 1
            id2 = holder[2*i+1] - 1
            donor_list.append(AcceptorDonor(atom_list[id1], atom_list[id2]))
        donor_list.changed = False
        # Now handle the acceptors (what is this used for??)
        nacc, holder = CharmmPsfFile._parse_psf_section(psf, int)
        acceptor_list = TrackedList()
        if len(holder) != nacc * 2:
            raise CharmmPSFError('Got %d indexes for %d acceptors' %
                                 (len(holder), ndon))
        for i in range(nacc):
            id1 = holder[2*i  ] - 1
            id2 = holder[2*i+1] - 1
            acceptor_list.append(AcceptorDonor(atom_list[id1], atom_list[id2]))
        acceptor_list.changed = False
        # Now get the NNB section. Not sure what this section is for or what it
        # does...
        nnb, holder = CharmmPsfFile._parse_psf_section(psf, int)
        # Now get the group sections
        pointers, holder = CharmmPsfFile._parse_psf_section(psf, int)
        group_list = TrackedList()
        try:
            ngrp, nst2 = pointers
        except ValueError:
            raise CharmmPSFError('Could not unpack GROUP pointers')
        group_list.nst2 = nst2
        # Now handle the groups
        if len(holder) != ngrp * 3:
            raise CharmmPSFError('Got %d indexes for %d groups' %
                                 (len(holder), ngrp))
        for i in range(ngrp):
            i1 = holder[3*i  ]
            i2 = holder[3*i+1]
            i3 = holder[3*i+2]
            group_list.append(Group(i1, i2, i3))
        group_list.changed = False
        # The next section might be the number of molecules or it might be the
        # cross-term (cmap) section. The first thing we'll do is determine
        # molecularity based on the atom connectivity. If every PSF file was
        # guaranteed to be "correct", we could just compare the MOLNT
        # section with the one we compute here. However, CHARMM GUI appears
        # to assign MOLNT as a dummy section (with all 1's), so this
        # approach will not work. Instead, look at the value of the pointer
        # and the number of entries in the group. If the # of entries is
        # NATOM, assume we have MOLNT section. Warn if the MOLNT section is
        # 'wrong'...
        pointer, holder = CharmmPsfFile._parse_psf_section(psf, int)

        # Assign all of the atoms to molecules recursively
        set_molecules(atom_list)
        molecule_list = [atom.marked for atom in atom_list]
        if len(holder) == len(atom_list):
            if molecule_list != holder:
                warnings.warn('Detected PSF molecule section that is WRONG. '
                              'Resetting molecularity.', CharmmPSFWarning)
            # We have a CHARMM PSF file; now do NUMLP/NUMLPH sections
            words = psf.readline().split()
            numlp = conv(words[0], int, 'numlp')
            numlph = conv(words[1], int, 'numlph')
            if numlp != 0 or numlph != 0:
                raise NotImplemented('Cannot currently handle PSFs with lone '
                                     'pairs defined in the NUMLP/NUMLPH '
                                     'section.')
            psf.readline() # blank
            # Now we get to the cross-term section
            ncrterm, holder = CharmmPsfFile._parse_psf_section(psf, int)
        else:
            ncrterm = pointer
        # At this point, ncrterm and holder are both set to the CMAP list for
        # VMD and non-VMD PSFs.
        # Now get the cmaps
        cmap_list = TrackedList()
        if len(holder) != ncrterm * 8:
            raise CharmmPSFError('Got %d CMAP indexes for %d cmap terms' %
                                 (len(holder), ncrterm))
        for i in range(ncrterm):
            id1 = holder[8*i  ] - 1
            id2 = holder[8*i+1] - 1
            id3 = holder[8*i+2] - 1
            id4 = holder[8*i+3] - 1
            id5 = holder[8*i+4] - 1
            id6 = holder[8*i+5] - 1
            id7 = holder[8*i+6] - 1
            id8 = holder[8*i+7] - 1
            cmap_list.append(
                    Cmap(atom_list[id1], atom_list[id2], atom_list[id3],
                         atom_list[id4], atom_list[id5], atom_list[id6],
                         atom_list[id7], atom_list[id8])
            )
        cmap_list.changed = False

        self.residue_list = residue_list
        self.atom_list = atom_list
        self.bond_list = bond_list
        self.angle_list = angle_list
        self.dihedral_list = dihedral_list
        self.improper_list = improper_list
        self.cmap_list = cmap_list
        self.donor_list = donor_list
        self.acceptor_list = acceptor_list
        self.group_list = group_list
        self.title = title
        self.flags = psf_flags

    def write_psf(self, dest, vmd=False):
        """
        Writes a PSF file from the stored molecule

        Parameters:
            - dest (str or file-like) : The place to write the output PSF file.
                    If it has a "write" attribute, it will be used to print the
                    PSF file. Otherwise, it will be treated like a string and a
                    file will be opened, printed, then closed
            - vmd (bool) : If True, it will write out a PSF in the format that
                    VMD prints it in (i.e., no NUMLP/NUMLPH or MOLNT sections)
        Example:
            >>> cs = CharmmPsfFile('testfiles/test.psf')
            >>> cs.write_psf('testfiles/test2.psf')
        """
        # See if this is an extended format
        ext = 'EXT' in self.flags
        own_handle = False
        # Index the atoms and residues
        self.residue_list.assign_indexes()
        self.atom_list.assign_indexes()
        if not hasattr(dest, 'write'):
            own_handle = True
            dest = open(dest, 'w')

        # Assign the formats we need to write with
        if ext:
            atmfmt1 = ('%10d %-8s %-8i %-8s %-8s %4d %10.6f %13.4f' + 11*' ')
            atmfmt2 = ('%10d %-8s %-8i %-8s %-8s %-4s %10.6f %13.4f' + 11*' ')
            intfmt = '%10d' # For pointers
        else:
            atmfmt1 = ('%8d %-4s %-4i %-4s %-4s %4d %10.6f %13.4f' + 11*' ')
            atmfmt2 = ('%8d %-4s %-4i %-4s %-4s %-4s %10.6f %13.4f' + 11*' ')
            intfmt = '%8d' # For pointers

        # Now print the header then the title
        dest.write('PSF ' + ' '.join(self.flags) + '\n')
        dest.write('\n')
        dest.write(intfmt % len(self.title) + ' !NTITLE\n')
        dest.write('\n'.join(self.title) + '\n\n')
        # Now time for the atoms
        dest.write(intfmt % len(self.atom_list) + ' !NATOM\n')
        # atmfmt1 is for CHARMM format (i.e., atom types are integers)
        # atmfmt is for XPLOR format (i.e., atom types are strings)
        for i, atom in enumerate(self.atom_list):
            if isinstance(atom.attype, str):
                fmt = atmfmt2
            else:
                fmt = atmfmt1
            atmstr = fmt % (i+1, atom.system, atom.residue.resnum,
                            atom.residue.resname, atom.name, atom.attype,
                            atom.charge, atom.mass)
            dest.write(atmstr + '   '.join(atom.props) + '\n')
        dest.write('\n')
        # Bonds
        dest.write(intfmt % len(self.bond_list) + ' !NBOND: bonds\n')
        for i, bond in enumerate(self.bond_list):
            dest.write((intfmt*2) % (bond.atom1.idx+1, bond.atom2.idx+1))
            if i % 4 == 3: # Write 4 bonds per line
                dest.write('\n')
        # See if we need to terminate
        if len(self.bond_list) % 4 != 0 or len(self.bond_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # Angles
        dest.write(intfmt % len(self.angle_list) + ' !NTHETA: angles\n')
        for i, angle in enumerate(self.angle_list):
            dest.write((intfmt*3) % (angle.atom1.idx+1, angle.atom2.idx+1,
                                     angle.atom3.idx+1)
            )
            if i % 3 == 2: # Write 3 angles per line
                dest.write('\n')
        # See if we need to terminate
        if len(self.angle_list) % 3 != 0 or len(self.angle_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # Dihedrals
        dest.write(intfmt % len(self.dihedral_list) + ' !NPHI: dihedrals\n')
        for i, dih in enumerate(self.dihedral_list):
            dest.write((intfmt*4) % (dih.atom1.idx+1, dih.atom2.idx+1,
                                     dih.atom3.idx+1, dih.atom4.idx+1)
            )
            if i % 2 == 1: # Write 2 dihedrals per line
                dest.write('\n')
        # See if we need to terminate
        if len(self.dihedral_list) % 2 != 0 or len(self.dihedral_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # Impropers
        dest.write(intfmt % len(self.improper_list) + ' !NIMPHI: impropers\n')
        for i, imp in enumerate(self.improper_list):
            dest.write((intfmt*4) % (imp.atom1.idx+1, imp.atom2.idx+1,
                                     imp.atom3.idx+1, imp.atom4.idx+1)
            )
            if i % 2 == 1: # Write 2 dihedrals per line
                dest.write('\n')
        # See if we need to terminate
        if len(self.improper_list) % 2 != 0 or len(self.improper_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # Donor section
        dest.write(intfmt % len(self.donor_list) + ' !NDON: donors\n')
        for i, don in enumerate(self.donor_list):
            dest.write((intfmt*2) % (don.atom1.idx+1, don.atom2.idx+1))
            if i % 4 == 3: # 4 donors per line
                dest.write('\n')
        if len(self.donor_list) % 4 != 0 or len(self.donor_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # Acceptor section
        dest.write(intfmt % len(self.acceptor_list) + ' !NACC: acceptors\n')
        for i, acc in enumerate(self.acceptor_list):
            dest.write((intfmt*2) % (acc.atom1.idx+1, acc.atom2.idx+1))
            if i % 4 == 3: # 4 donors per line
                dest.write('\n')
        if len(self.acceptor_list) % 4 != 0 or len(self.acceptor_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # NNB section ??
        dest.write(intfmt % 0 + ' !NNB\n\n')
        for i in range(len(self.atom_list)):
            dest.write(intfmt % 0)
            if i % 8 == 7: # Write 8 0's per line
                dest.write('\n')
        if len(self.atom_list) % 8 != 0: dest.write('\n')
        dest.write('\n')
        # Group section
        dest.write((intfmt*2) % (len(self.group_list), self.group_list.nst2))
        dest.write(' !NGRP NST2\n')
        for i, gp in enumerate(self.group_list):
            dest.write((intfmt*3) % (gp.bs, gp.type, gp.move))
            if i % 3 == 2: dest.write('\n')
        if len(self.group_list) % 3 != 0 or len(self.group_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # The next two sections are never found in VMD prmtops...
        if not vmd:
            # Molecule section; first set molecularity
            set_molecules(self.atom_list)
            mollist = [a.marked for a in self.atom_list]
            dest.write(intfmt % max(mollist) + ' !MOLNT\n')
            for i, atom in enumerate(self.atom_list):
                dest.write(intfmt % atom.marked)
                if i % 8 == 7: dest.write('\n')
            if len(self.atom_list) % 8 != 0: dest.write('\n')
            dest.write('\n')
            # NUMLP/NUMLPH section
            dest.write((intfmt*2) % (0, 0) + ' !NUMLP NUMLPH\n')
            dest.write('\n')
        # CMAP section
        dest.write(intfmt % len(self.cmap_list) + ' !NCRTERM: cross-terms\n')
        for i, cmap in enumerate(self.cmap_list):
            dest.write((intfmt*4) % (cmap.atom1.idx+1, cmap.atom2.idx+1,
                                     cmap.atom3.idx+1, cmap.atom4.idx+1)
            )
            if cmap.consecutive:
                dest.write((intfmt*4) % (cmap.atom2.idx+1, cmap.atom3.idx+1,
                                         cmap.atom4.idx+1, cmap.atom5.idx+1)
                )
            else:
                dest.write((intfmt*4) % (cmap.atom5.idx+1, cmap.atom6.idx+1,
                                         cmap.atom7.idx+1, cmap.atom8.idx+1)
                )
            dest.write('\n')
        # Done!
        # If we opened our own handle, close it
        if own_handle:
            dest.close()
        
    def load_parameters(self, parmset):
        """
        Loads parameters from a parameter set that was loaded via CHARMM RTF,
        PAR, and STR files.

        Parameters:
            - parmset (ParameterSet) : List of all parameters

        Notes:
            - If any parameters that are necessary cannot be found, a
              MissingParameter exception is raised.

            - If any dihedral or improper parameters cannot be found, I will try
              inserting wildcards (at either end for dihedrals and as the two
              central atoms in impropers) and see if that matches.  Wild-cards
              will apply ONLY if specific parameters cannot be found.

            - This method will expand the dihedral_list attribute by adding a
              separate Dihedral object for each term for types that have a
              multi-term expansion
        """
        # First load the atom types
        types_are_int = False
        for atom in self.atom_list:
            try:
                if isinstance(atom.attype, int):
                    atype = parmset.atom_types_int[atom.attype]
                    types_are_int = True # if we have to change back
                else:
                    atype = parmset.atom_types_str[atom.attype]
            except KeyError:
                raise MissingParameter('Could not find atom type for %s' %
                                       atom.attype)
            atom.type = atype
            # Change to string attype to look up the rest of the parameters
            atom.type_to_str()

        # Next load all of the bonds
        for bond in self.bond_list:
            # Construct the key
            key = (min(bond.atom1.attype, bond.atom2.attype),
                   max(bond.atom1.attype, bond.atom2.attype))
            try:
                bond.bond_type = parmset.bond_types[key]
            except KeyError:
                raise MissingParameter('Missing bond type for %r' % bond)
        # Next load all of the angles. If a Urey-Bradley term is defined for
        # this angle, also build the urey_bradley and urey_bradley_type lists
        self.urey_bradley_list = TrackedList()
        for ang in self.angle_list:
            # Construct the key
            key = (min(ang.atom1.attype, ang.atom3.attype), ang.atom2.attype,
                   max(ang.atom1.attype, ang.atom3.attype))
            try:
                ang.angle_type = parmset.angle_types[key]
                ubt = parmset.urey_bradley_types[key]
                if ubt is not NoUreyBradley:
                    ub = UreyBradley(ang.atom1, ang.atom3, ubt)
                    self.urey_bradley_list.append(ub)
            except KeyError:
                raise MissingParameter('Missing angle type for %r' % ang)
        # Next load all of the dihedrals. This is a little trickier since we
        # need to back up the existing dihedral list and replace it with a
        # longer one that has only one Fourier term per Dihedral instance.
        dihedral_list = self.dihedral_list
        self.dihedral_list = TrackedList()
        active_dih_list = set()
        for dih in dihedral_list:
            # Store the atoms
            a1, a2, a3, a4 = dih.atom1, dih.atom2, dih.atom3, dih.atom4
            at1, at2, at3, at4 = a1.attype, a2.attype, a3.attype, a4.attype
            # First see if the exact dihedral is specified
            key = min((at1,at2,at3,at4), (at4,at3,at2,at1))
            if not key in parmset.dihedral_types:
                # Check for wild-cards
                key = min(('X',at2,at3,'X'), ('X',at3,at2,'X'))
                if not key in parmset.dihedral_types:
                    raise MissingParameter('No dihedral parameters found for '
                                           '%r' % dih)
            dtlist = parmset.dihedral_types[key]
            for i, dt in enumerate(dtlist):
                self.dihedral_list.append(Dihedral(a1, a2, a3, a4, dt))
                # See if we include the end-group interactions for this
                # dihedral. We do IFF it is the last or only dihedral term and
                # it is NOT in the angle/bond partners
                pair = tuple(sorted([a1.idx, a4.idx]))
                if i != len(dtlist) - 1:
                    self.dihedral_list[-1].end_groups_active = False
                elif a1 in a4.bond_partners or a1 in a4.angle_partners:
                    self.dihedral_list[-1].end_groups_active = False
                elif pair in active_dih_list:
                    self.dihedral_list[-1].end_groups_active = False
                else:
                    active_dih_list.add(pair)
        # Now do the impropers
        for imp in self.improper_list:
            # Store the atoms
            a1, a2, a3, a4 = imp.atom1, imp.atom2, imp.atom3, imp.atom4
            at1, at2, at3, at4 = a1.attype, a2.attype, a3.attype, a4.attype
            key = tuple(sorted([at1, at2, at3, at4]))
            if not key in parmset.improper_types:
                # Check for wild-cards
                for anchor in (at2, at3, at4):
                    key = tuple(sorted([at1, anchor, 'X', 'X']))
                    if key in parmset.improper_types:
                        break # This is the right key
            try:
                imp.improper_type = parmset.improper_types[key]
            except KeyError:
                raise MissingParameter('No improper parameters found for %r' %
                                       imp)
        # Now do the cmaps. These will not have wild-cards
        for cmap in self.cmap_list:
            # Store the atoms for easy reference
            if cmap.consecutive:
                a1, a2, a3, a4 = cmap.atom1, cmap.atom2, cmap.atom3, cmap.atom4
                a5, a6, a7, a8 = cmap.atom2, cmap.atom3, cmap.atom4, cmap.atom5
            else:
                a1, a2, a3, a4 = cmap.atom1, cmap.atom2, cmap.atom3, cmap.atom4
                a5, a6, a7, a8 = cmap.atom5, cmap.atom6, cmap.atom7, cmap.atom8
            at1, at2, at3, at4 = a1.attype, a2.attype, a3.attype, a4.attype
            at5, at6, at7, at8 = a5.attype, a6.attype, a7.attype, a8.attype
            # Construct the keys
            k1 = list(min((at1,at2,at3,at4), (at4,at3,at2,at1)))
            k2 = list(min((at5,at6,at7,at8), (at8,at7,at6,at5)))
            key = tuple(k1 + k2)
            try:
                cmap.cmap_type = parmset.cmap_types[key]
            except KeyError:
                raise MissingParameter('No CMAP parameters found for %r' % cmap)
        # If the types started out as integers, change them back
        if types_are_int:
            for atom in self.atom_list: atom.type_to_int()

    def set_coordinates(self, positions, velocities=None):
        """
        This method loads the coordinates and velocity information from an
        external object or passed data.

        Parameters:
            - positions (list of floats) : A 3-N length iterable with all of the
                coordinates in the order [x1, y1, z1, x2, y2, z2, ...].
            - velocities (list of floats) : If not None, is the velocity
                equivalent of the positions
        """
        if len(positions) / 3 != len(self.atom_list):
            raise ValueError('Coordinates given for %s atoms, but %d atoms '
                             'exist in this structure.' %
                             (len(positions)/3, len(self.atom_list)))
        # Now assign all of the atoms positions
        for i, atom in enumerate(self.atom_list):
            atom.xx = positions[3*i  ]
            atom.xy = positions[3*i+1]
            atom.xz = positions[3*i+2]

        # Do velocities if given
        if velocities is not None:
            if len(velocities) / 3 != len(self.atom_list):
                raise ValueError('Velocities given for %s atoms, but %d atoms '
                                 'exist in this structure.' %
                                 (len(velocities)/3, len(self.atom_list)))
            for i, atom in enumerate(self.atom_list):
                atom.vx = velocities[3*i  ]
                atom.vy = velocities[3*i+1]
                atom.vz = velocities[3*i+2]
            self.velocities = velocities

        self.positions = positions

    def set_box(self, box):
        """
        Sets the periodic box boundary conditions.

        Parameters:
            - box (list of 6 floats) : A list of 6 numbers representing a, b, c,
                alpha, beta, and gamma, respectively (box lengths and angles).
                If None, the system is assumed to be aperiodic
        Notes:
            The box can alternatively be a list of 3 numbers representing the
            box lengths. In this case, the angles are assumed to be 90 degrees

            The box here is copied via slicing, so changing the box that was
            passed in will have no effect after set_box is called.
        """
        if box is None:
            self.box = None
        elif len(box) == 6:
            self.box = list(box[:])
        elif len(box) == 3:
            self.box = list(box[:]) + [90.0, 90.0, 90.0]
        else:
            raise ValueError('set_box requires 3 box lengths, 3 box lengths '
                             'and 3 angles, or None for no box')

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def set_molecules(atom_list):
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
    setrecursionlimit(max(len(atom_list), getrecursionlimit()))

    # Unmark all atoms so we can track which molecule each goes into
    atom_list.unmark()

    # The molecule "ownership" list
    owner = []
    # The way I do this is via a recursive algorithm, in which
    # the "set_owner" method is called for each bonded partner an atom
    # has, which in turn calls set_owner for each of its partners and 
    # so on until everything has been assigned.
    molecule_number = 1 # which molecule number we are on
    for i in range(len(atom_list)):
        # If this atom has not yet been "owned", make it the next molecule
        # However, we only increment which molecule number we're on if 
        # we actually assigned a new molecule (obviously)
        if not atom_list[i].marked:
            tmp = [i]
            _set_owner(atom_list, tmp, i, molecule_number)
            # Make sure the atom indexes are sorted
            tmp.sort()
            owner.append(tmp)
            molecule_number += 1
    return owner

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _set_owner(atom_list, owner_array, atm, mol_id):
    """ Recursively sets ownership of given atom and all bonded partners """
    atom_list[atm].marked = mol_id
    for partner in atom_list[atm].bond_partners:
        if not partner.marked:
            owner_array.append(partner.idx)
            _set_owner(atom_list, owner_array, partner.idx, mol_id)
        elif partner.marked != mol_id:
            raise MoleculeError('Atom %d in multiple molecules' % 
                                partner.idx)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == '__main__':
    import doctest
    doctest.testmod()
