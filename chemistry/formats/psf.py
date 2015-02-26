"""
This module contains classes for reading and writing CHARMM- and XPLOR-style PSF
files
"""
from chemistry.charmm import CharmmPsfFile
# TODO -- move this functionality to a more centralized location
from chemistry.charmm.psf import set_molecules
from chemistry.formats.io import genopen, TextToBinaryFile
from chemistry.formats.registry import FileFormatType

class PSFFile(object):
    """
    CHARMM- or XPLOR-style PSF file parser and writer. This class is
    specifically a holder for the writing functionality and a vessel for
    automatic file type detection. If you wish to instantiate a PSF file
    directly, use :class:`chemistry.charmm.CharmmPsfFile` or the
    :func:`chemistry.formats.load_file` function instead.
    """
    __metaclass__ = FileFormatType

    #===================================================

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as a CHARMM PSF file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is a CHARMM or Xplor-style PSF file
        """
        f = genopen(filename, 'r')
        line = f.readline().decode()
        f.close()
        return line.strip().startswith('PSF')

    #===================================================

    @staticmethod
    def parse(filename):
        """ Read a CHARMM- or XPLOR-style PSF file

        Parameters
        ----------
        filename : str
            Name of the file to parse

        Returns
        -------
        psf_file : :class:`CharmmPsfFile`
            The PSF file instance with all information loaded
        """
        return CharmmPsfFile(filename)

    #===================================================

    @staticmethod
    def write(struct, dest, vmd=False):
        """
        Writes a PSF file from the stored molecule

        Parameters
        ----------
        struct : :class:`Structure`
            The Structure instance from which the PSF should be written
        dest : str or file-like
            The place to write the output PSF file.  If it has a "write"
            attribute, it will be used to print the PSF file. Otherwise, it will
            be treated like a string and a file will be opened, printed, then
            closed
        vmd : bool 
            If True, it will write out a PSF in the format that VMD prints it in
            (i.e., no NUMLP/NUMLPH or MOLNT sections)

        Examples
        --------
        >>> cs = CharmmPsfFile('testfiles/test.psf')
        >>> cs.write_psf('testfiles/test2.psf')
        """
        # See if this is an extended format
        try:
            ext = 'EXT' in struct.flags
        except AttributeError:
            ext = True
        own_handle = False
        # Index the atoms and residues
        if not hasattr(dest, 'write'):
            own_handle = True
            dest = TextToBinaryFile(genopen(dest, 'w'))

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
        dest.write('PSF ')
        if hasattr(struct, 'flags'):
            dest.write(' '.join(struct.flags))
        else:
            dest.write('EXT') # EXT is always active
        dest.write('\n\n')
        dest.write(intfmt % len(struct.title) + ' !NTITLE\n')
        dest.write('\n'.join(struct.title) + '\n\n')
        # Now time for the atoms
        dest.write(intfmt % len(struct.atoms) + ' !NATOM\n')
        # atmfmt1 is for CHARMM format (i.e., atom types are integers)
        # atmfmt is for XPLOR format (i.e., atom types are strings)
        for i, atom in enumerate(struct.atoms):
            typ = atom.type
            if isinstance(atom.type, str):
                fmt = atmfmt2
                if not atom.type: typ = atom.name
            else:
                fmt = atmfmt1
            if hasattr(atom, 'segid'):
                segid = atom.segid
            else:
                segid = 'SYS'
            atmstr = fmt % (i+1, segid, atom.residue.number,
                            atom.residue.name, atom.name, typ,
                            atom.charge, atom.mass)
            if hasattr(atom, 'props'):
                dest.write(atmstr + '   '.join(atom.props) + '\n')
            else:
                dest.write('\n')
        dest.write('\n')
        # Bonds
        dest.write(intfmt % len(struct.bonds) + ' !NBOND: bonds\n')
        for i, bond in enumerate(struct.bonds):
            dest.write((intfmt*2) % (bond.atom1.idx+1, bond.atom2.idx+1))
            if i % 4 == 3: # Write 4 bonds per line
                dest.write('\n')
        # See if we need to terminate
        if len(struct.bonds) % 4 != 0 or len(struct.bonds) == 0:
            dest.write('\n')
        dest.write('\n')
        # Angles
        dest.write(intfmt % len(struct.angles) + ' !NTHETA: angles\n')
        for i, angle in enumerate(struct.angles):
            dest.write((intfmt*3) % (angle.atom1.idx+1, angle.atom2.idx+1,
                                     angle.atom3.idx+1)
            )
            if i % 3 == 2: # Write 3 angles per line
                dest.write('\n')
        # See if we need to terminate
        if len(struct.angles) % 3 != 0 or len(struct.angles) == 0:
            dest.write('\n')
        dest.write('\n')
        # Dihedrals
        dest.write(intfmt % len(struct.dihedrals) + ' !NPHI: dihedrals\n')
        for i, dih in enumerate(struct.dihedrals):
            dest.write((intfmt*4) % (dih.atom1.idx+1, dih.atom2.idx+1,
                                     dih.atom3.idx+1, dih.atom4.idx+1)
            )
            if i % 2 == 1: # Write 2 dihedrals per line
                dest.write('\n')
        # See if we need to terminate
        if len(struct.dihedrals) % 2 != 0 or len(struct.dihedrals) == 0:
            dest.write('\n')
        dest.write('\n')
        # Impropers
        dest.write(intfmt % len(struct.impropers) + ' !NIMPHI: impropers\n')
        for i, imp in enumerate(struct.impropers):
            dest.write((intfmt*4) % (imp.atom1.idx+1, imp.atom2.idx+1,
                                     imp.atom3.idx+1, imp.atom4.idx+1)
            )
            if i % 2 == 1: # Write 2 dihedrals per line
                dest.write('\n')
        # See if we need to terminate
        if len(struct.impropers) % 2 != 0 or len(struct.impropers) == 0:
            dest.write('\n')
        dest.write('\n')
        # Donor section
        dest.write(intfmt % len(struct.donors) + ' !NDON: donors\n')
        for i, don in enumerate(struct.donors):
            dest.write((intfmt*2) % (don.atom1.idx+1, don.atom2.idx+1))
            if i % 4 == 3: # 4 donors per line
                dest.write('\n')
        if len(struct.donors) % 4 != 0 or len(struct.donors) == 0:
            dest.write('\n')
        dest.write('\n')
        # Acceptor section
        dest.write(intfmt % len(struct.acceptors) + ' !NACC: acceptors\n')
        for i, acc in enumerate(struct.acceptors):
            dest.write((intfmt*2) % (acc.atom1.idx+1, acc.atom2.idx+1))
            if i % 4 == 3: # 4 donors per line
                dest.write('\n')
        if len(struct.acceptors) % 4 != 0 or len(struct.acceptors) == 0:
            dest.write('\n')
        dest.write('\n')
        # NNB section ??
        dest.write(intfmt % 0 + ' !NNB\n\n')
        for i in xrange(len(struct.atoms)):
            dest.write(intfmt % 0)
            if i % 8 == 7: # Write 8 0's per line
                dest.write('\n')
        if len(struct.atoms) % 8 != 0: dest.write('\n')
        dest.write('\n')
        # Group section
        dest.write((intfmt*2) % (len(struct.groups), struct.groups.nst2))
        dest.write(' !NGRP NST2\n')
        for i, gp in enumerate(struct.groups):
            dest.write((intfmt*3) % (gp.bs, gp.type, gp.move))
            if i % 3 == 2: dest.write('\n')
        if len(struct.groups) % 3 != 0 or len(struct.groups) == 0:
            dest.write('\n')
        dest.write('\n')
        # The next two sections are never found in VMD prmtops...
        if not vmd:
            # Molecule section; first set molecularity
            set_molecules(struct.atoms)
            mollist = [a.marked for a in struct.atoms]
            dest.write(intfmt % max(mollist) + ' !MOLNT\n')
            for i, atom in enumerate(struct.atoms):
                dest.write(intfmt % atom.marked)
                if i % 8 == 7: dest.write('\n')
            if len(struct.atoms) % 8 != 0: dest.write('\n')
            dest.write('\n')
            # NUMLP/NUMLPH section
            dest.write((intfmt*2) % (0, 0) + ' !NUMLP NUMLPH\n')
            dest.write('\n')
        # CMAP section
        dest.write(intfmt % len(struct.cmaps) + ' !NCRTERM: cross-terms\n')
        for i, cmap in enumerate(struct.cmaps):
            dest.write((intfmt*8) % (cmap.atom1.idx+1, cmap.atom2.idx+1,
                                     cmap.atom3.idx+1, cmap.atom4.idx+1,
                                     cmap.atom2.idx+1, cmap.atom3.idx+1,
                                     cmap.atom4.idx+1, cmap.atom5.idx+1)
            )
            dest.write('\n')
        # Done!
        # If we opened our own handle, close it
        if own_handle:
            dest.close()

    #===================================================
