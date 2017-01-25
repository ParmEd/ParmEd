"""
This package contains classes responsible for reading and writing both PDB and
PDBx/mmCIF files.
"""
from __future__ import division, print_function, absolute_import

from contextlib import closing
import io
import ftplib
import numpy as np
from parmed.exceptions import PDBError, PDBWarning
from parmed.formats.pdbx import PdbxReader, PdbxWriter, containers
from parmed.formats.registry import FileFormatType
from parmed.periodic_table import AtomicNum, Mass, Element, element_by_name
from parmed.residue import AminoAcidResidue, RNAResidue, DNAResidue, WATER_NAMES
from parmed.modeller import StandardBiomolecularResidues
from parmed.structure import Structure
from parmed.topologyobjects import Atom, ExtraPoint, Bond
from parmed.symmetry import Symmetry
from parmed.utils.io import genopen
from parmed.utils.six import iteritems, string_types, add_metaclass, PY3
from parmed.utils.six.moves import range
import re
import warnings

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _compare_atoms(old_atom, new_atom, resname, resid, chain, segid, inscode):
    """
    Compares two atom instances, along with the residue name, number, and chain
    identifier, to determine if two atoms are actually the *same* atom, but
    simply different conformations

    Parameters
    ----------
    old_atom : :class:`Atom`
        The original atom that has been added to the structure already
    new_atom : :class:`Atom`
        The new atom that we want to see if it is the same as the old atom
    resname : ``str``
        The name of the residue that the new atom would belong to
    resid : ``int``
        The number of the residue that the new atom would belong to
    chain : ``str``
        The chain identifier that the new atom would belong to
    segid : ``str``
        The segment identifier for the molecule
    inscode : ``str``
        The insertion code for the residue

    Returns
    -------
    True if they are the same atom, False otherwise
    """
    if old_atom.name != new_atom.name: return False
    if old_atom.residue.name != resname: return False
    if old_atom.residue.number != resid: return False
    if old_atom.residue.chain != chain.strip(): return False
    if old_atom.residue.segid != segid.strip(): return False
    if old_atom.residue.insertion_code != inscode.strip(): return False
    return True

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _standardize_resname(resname):
    """ Looks up a standardized residue name for the given resname """
    try:
        return AminoAcidResidue.get(resname, abbronly=True).abbr, False
    except KeyError:
        try:
            return RNAResidue.get(resname).abbr, False
        except KeyError:
            try:
                return DNAResidue.get(resname).abbr, False
            except KeyError:
                if resname.strip() in WATER_NAMES:
                    return 'HOH', True
                else:
                    return resname, True

def _is_hetatm(resname):
    """ Sees if residue name is "standard", otherwise, we need to use HETATM to
    print atom records instead of ATOM
    """
    if len(resname) != 3:
        return not (RNAResidue.has(resname) or DNAResidue.has(resname))
    return not (AminoAcidResidue.has(resname) or RNAResidue.has(resname)
            or DNAResidue.has(resname))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

@add_metaclass(FileFormatType)
class PDBFile(object):
    """ Standard PDB file format parser and writer """
    #===================================================

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as a PDB file

        Parameters
        ----------
        filename : str or file object
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is a PDB file
        """
        if isinstance(filename, string_types):
            own_handle = True
            fileobject = genopen(filename, 'r')
        elif hasattr(filename, 'read'):
            own_handle = False
            fileobject = filename

        try:
            for line in fileobject:
                if line[:6] in ('CRYST1', 'END   ', 'END', 'HEADER', 'NUMMDL',
                        'MASTER', 'AUTHOR', 'CAVEAT', 'COMPND', 'EXPDTA',
                        'MDLTYP', 'KEYWDS', 'OBSLTE', 'SOURCE', 'SPLIT ',
                        'SPRSDE', 'TITLE ', 'ANISOU', 'CISPEP', 'CONECT',
                        'DBREF ', 'HELIX ', 'HET   ', 'LINK  ', 'MODRES',
                        'REVDAT', 'SEQADV', 'SHEET ', 'SSBOND', 'FORMUL',
                        'HETNAM', 'HETSYN', 'SEQRES', 'SITE  ', 'ENDMDL',
                        'MODEL ', 'TER   ', 'JRNL  ', 'REMARK', 'TER'):
                    continue
                # Hack to support reduce-added flags
                elif line[:6] == 'USER  ' and line[6:9] == 'MOD':
                    continue
                elif line[:5] in ('ORIGX', 'SCALE', 'MTRIX'):
                    if line[5] not in '123':
                        return False
                elif line[:6] in ('ATOM  ', 'HETATM'):
                    atnum, atname = line[6:11], line[12:16]
                    resname, resid = line[17:21], line[22:26]
                    x, y, z = line[30:38], line[38:46], line[46:54]
                    occupancy, bfactor = line[54:60], line[60:66]
                    elem = line[76:78]
                    # Check for various attributes. This is the first atom, so
                    # we can assume we haven't gotten into the regime of "weird"
                    # yet, like hexadecimal atom/residue indices.
                    if not atnum.strip().lstrip('-').isdigit(): return False
                    if atname.strip().isdigit(): return False
                    if not resname.strip(): return False
                    if not resid.strip().lstrip('-').isdigit(): return False
                    try:
                        float(x), float(y), float(z)
                    except ValueError:
                        return False
                    if occupancy.strip():
                        try:
                            float(occupancy)
                        except ValueError:
                            return False
                    if bfactor.strip():
                        try:
                            float(bfactor)
                        except ValueError:
                            return False
                    if elem.strip():
                        if any(x.isdigit() for x in elem):
                            return False
                    return True
                else:
                    return False
            return False
        finally:
            if own_handle:
                fileobject.close()

    #===================================================

    @staticmethod
    def download(pdb_id, timeout=10, saveto=None):
        """
        Goes to the wwPDB website and downloads the requested PDB, loading it
        as a :class:`Structure` instance

        Parameters
        ----------
        pdb_id : str
            The 4-letter PDB ID to try and download from the RCSB PDB database
        timeout : float, optional
            The number of seconds to wait before raising a timeout error.
            Default is 10 seconds
        saveto : str, optional
            If provided, this will be treated as a file name to which the PDB
            file will be saved. If None (default), no PDB file will be written.
            This will be a verbatim copy of the downloaded PDB file, unlike the
            somewhat-stripped version you would get by using
            :meth:`Structure.write_pdb <parmed.structure.Structure.write_pdb>`

        Returns
        -------
        struct : :class:`Structure <parmed.structure.Structure>`
            Structure instance populated by the requested PDB

        Raises
        ------
        socket.timeout if the connection times out while trying to contact the
        FTP server

        IOError if there is a problem retrieving the requested PDB or writing a
        requested ``saveto`` file

        ImportError if the gzip module is not available

        TypeError if pdb_id is not a 4-character string
        """
        import gzip
        if not isinstance(pdb_id, string_types) or len(pdb_id) != 4:
            raise ValueError('pdb_id must be the 4-letter PDB code')

        pdb_id = pdb_id.lower()
        ftp = ftplib.FTP('ftp.wwpdb.org', timeout=timeout)
        ftp.login()
        fileobj = io.BytesIO()
        try:
            ftp.retrbinary('RETR /pub/pdb/data/structures/divided/pdb/'
                           '%s/pdb%s.ent.gz' % (pdb_id[1:3], pdb_id),
                           fileobj.write)
        except ftplib.all_errors as err:
            raise IOError('Could not retrieve PDB ID %s; %s' % (pdb_id, err))
        finally:
            ftp.close()
        # Rewind, wrap it in a GzipFile and send it to parse
        fileobj.seek(0)
        if PY3:
            fileobj = io.TextIOWrapper(gzip.GzipFile(fileobj=fileobj, mode='r'))
        else:
            fileobj = gzip.GzipFile(fileobj=fileobj, mode='r')
        if saveto is not None:
            with closing(genopen(saveto, 'w')) as f:
                f.write(fileobj.read())
            fileobj.seek(0)
        return PDBFile.parse(fileobj)

    #===================================================

    _relatere = re.compile(r'RELATED ID: *(\w+) *RELATED DB: *(\w+)', re.I)

    @staticmethod
    def parse(filename, skip_bonds=False):
        """ Read a PDB file and return a populated `Structure` class

        Parameters
        ----------
        filename : str or file-like
            Name of the PDB file to read, or a file-like object that can iterate
            over the lines of a PDB. Compressed file names can be specified and
            are determined by file-name extension (e.g., file.pdb.gz,
            file.pdb.bz2)
        skip_bonds : bool, optional
            If True, skip trying to assign bonds. This can save substantial time
            when parsing large files with non-standard residue names. However,
            no bonds are assigned. This is OK if, for instance, the PDB file is
            being parsed simply for its coordinates. This may also reduce
            element assignment if element information is not present in the PDB
            file already. Default is False.

        Metadata
        --------
        The PDB parser also adds metadata to the returned Structure object that
        may be present in the PDB file

        experimental : ``str``
            EXPDTA record
        journal : ``str``
            JRNL record
        authors : ``str``
            AUTHOR records
        keywords : ``str``
            KEYWDS records
        doi : ``str``
            DOI from the JRNL record
        pmid : ``str``
            PMID from the JRNL record
        journal_authors : ``str``
            Author info from the JRNL record
        volume : ``str``
            Volume of the published article from the JRNL record
        page : ``str``
            Page of the published article from the JRNL record
        title : ``str``
            TITL section of the JRNL record
        year : ``int``
            Year that the article was published, from the JRNL record
        resolution : ``float``
            The X-RAY resolution in Angstroms, or None if not found
        related_entries : ``list of (str, str)``
            List of entries in other databases

        Returns
        -------
        structure : :class:`Structure`
            The Structure object initialized with all of the information from
            the PDB file.  No bonds or other topological features are added by
            default.
        """
        if isinstance(filename, string_types):
            own_handle = True
            fileobj = genopen(filename, 'r')
        else:
            own_handle = False
            fileobj = filename

        struct = Structure()
        # Add metadata fields
        struct.experimental = struct.journal = struct.authors = ''
        struct.keywords = struct.doi = struct.pmid = ''
        struct.journal_authors = struct.volume_page = struct.title = ''
        struct.year = struct.resolution = None
        struct.related_entries = []
        modelno = 1 # For PDB files with multiple MODELs
        atomno = 0
        coordinates = []
        all_coordinates = []

        # Support hexadecimal numbering like that printed by VMD
        last_atom = Atom()
        last_atom_added = None
        last_resid = 1
        resend = 26
        res_hex = False
        atom_hex = False
        atom_overflow = False
        ZEROSET = set('0')
        altloc_ids = set()
        _symmetry_lines = []

        try:
            for line in fileobj:
                if 'REMARK 290   SMTRY' in line:
                    _symmetry_lines.append(line)
                rec = line[:6]
                if rec == 'ATOM  ' or rec == 'HETATM':
                    atomno += 1
                    atnum, atname, altloc = line[6:11], line[12:16], line[16]
                    resname, chain = line[17:21], line[21]
                    resid, inscode = line[22:resend], line[26]
                    x, y, z = line[30:38], line[38:46], line[46:54]
                    occupancy, bfactor = line[54:60], line[60:66]
                    elem, chg = line[76:78], line[78:80]
                    segid = line[72:76].strip() # CHARMM-specific
                    atname = atname.strip()
                    altloc = altloc.strip()
                    resname = resname.strip()
                    chain = chain.strip()
                    inscode = inscode.strip()

                    elem = '%-2s' % elem # Make sure we have at least 2 chars
                    if elem[0] == ' ': elem = elem[1] + ' '
                    try:
                        atsym = (elem[0] + elem[1].lower()).strip()
                        atomic_number = AtomicNum[atsym]
                        mass = Mass[atsym]
                    except KeyError:
                        # Now try based on the atom name... but don't try too
                        # hard (e.g., don't try to differentiate b/w Ca and C)
                        elem = element_by_name(atname)
                        atomic_number = AtomicNum[elem]
                        mass = Mass[elem]
                    try:
                        bfactor = float(bfactor)
                    except ValueError:
                        bfactor = 0.0
                    try:
                        occupancy = float(occupancy)
                    except ValueError:
                        occupancy = 0.0
                    # Figure out what my residue number is and see if the PDB is
                    # outputting residue numbers in hexadecimal (e.g., VMD)
                    if last_resid >= 9999 and resend == 26:
                        if not res_hex and resid == '9999':
                            resid = 9999
                        elif not res_hex:
                            res_hex = int(resid, 16) == 10000
                        # So now we know if we use hexadecimal or not. If we do,
                        # convert. Otherwise, stay put
                        if res_hex:
                            try:
                                resid = int(resid, 16)
                            except ValueError:
                                if resid == '****':
                                    resid = None # Figure out by unique atoms
                                else:
                                    raise
                        elif resid == '1000' and line[26] == '0':
                            resend += 1
                            resid = 10000
                        else:
                            resid = int(resid)
                    elif resend > 26:
                        # VMD extends the field now... ugh.
                        if resid[0] == '1' and set(resid[1:]) == ZEROSET:
                            if line[resend] == '0':
                                resid = int(resid) * 10
                                resend += 1
                            else:
                                resid = int(resid)
                        else:
                            resid = int(resid)
                    else:
                        resid = int(resid)
                    # If the number has cycled, it too may be hexadecimal
                    if atom_hex:
                        try:
                            atnum = int(atnum, 16)
                        except ValueError:
                            if set(atnum) == set('*'):
                                atom_overflow = True
                                atnum = last_atom_added.number + 1
                            else:
                                raise
                    elif atom_overflow:
                        atnum = last_atom_added.number + 1
                    else:
                        try:
                            atnum = int(atnum)
                        except ValueError:
                            if set(atnum) == set('*'):
                                atom_overflow = True
                                atnum = last_atom_added.number + 1
                            else:
                                atnum = int(atnum, 16)
                                atom_hex = True
                    # It's possible that the residue number has cycled so much
                    # that it is now filled with ****'s. In that case, start a
                    # new residue if the current residue repeats the same atom
                    # name as # the 'last' residue.
                    if resid is None:
                        # If the last residue is number 0xffff, then this is the
                        # first residue that has overridden, so make it a new
                        # residue
                        if struct.residues[-1].number == 0xffff:
                            resid = struct.residues[-1].number + 1
                        else:
                            for atom in struct.residues[-1]:
                                if atom.name == atname:
                                    resid = last_resid + 1
                                    break
                    if resid is None:
                        # Still part of the last residue
                        resid = last_resid
                    last_resid = resid
                    try:
                        chg = float(chg)
                    except ValueError:
                        chg = 0.0
                    if atname in ('EP', 'LP'): # lone pair
                        atom = ExtraPoint(atomic_number=atomic_number,
                                name=atname, charge=chg, mass=mass,
                                occupancy=occupancy, bfactor=bfactor,
                                altloc=altloc, number=atnum)
                    else:
                        atom = Atom(atomic_number=atomic_number, name=atname,
                                charge=chg, mass=mass, occupancy=occupancy,
                                bfactor=bfactor, altloc=altloc, number=atnum)
                    atom.xx, atom.xy, atom.xz = float(x), float(y), float(z)
                    if (_compare_atoms(last_atom, atom, resname, resid, chain,
                                       segid, inscode) and altloc):
                        atom.residue = last_atom.residue
                        last_atom.other_locations[altloc] = atom
                        altloc_ids.add(atom.number)
                        last_atom_added = atom
                        continue
                    last_atom = last_atom_added = atom
                    if modelno == 1:
                        struct.add_atom(atom, resname, resid, chain,
                                        inscode, segid)
                    else:
                        try:
                            orig_atom = struct.atoms[atomno-1]
                        except IndexError:
                            raise PDBError('Extra atom in MODEL %d' % modelno)
                        if (orig_atom.residue.name != resname.strip()
                                or orig_atom.name != atname.strip()):
                            raise PDBError('Atom %d differs in MODEL %d [%s %s '
                                           'vs. %s %s]' % (atomno, modelno,
                                           orig_atom.residue.name,
                                           orig_atom.name, resname, atname))
                    coordinates.extend([atom.xx, atom.xy, atom.xz])
                elif rec == 'ANISOU':
                    try:
                        atnum = int(line[6:11])
                    except ValueError:
                        warnings.warn('Problem parsing atom number from ANISOU '
                                      'record', PDBWarning)
                        continue # Skip the rest of this record
                    aname = line[12:16].strip()
                    altloc = line[16].strip()
                    rname = line[17:21].strip()
                    chain = line[21].strip()
                    try:
                        resid = int(line[22:26])
                    except ValueError:
                        warnings.warn('Problem parsing residue number from '
                                      'ANISOU record', PDBWarning)
                        continue # Skip the rest of this record
                    icode = line[26].strip()
                    try:
                        u11 = int(line[28:35])
                        u22 = int(line[35:42])
                        u33 = int(line[42:49])
                        u12 = int(line[49:56])
                        u13 = int(line[56:63])
                        u23 = int(line[63:70])
                    except ValueError:
                        warnings.warn('Problem parsing anisotropic factors '
                                      'from ANISOU record', PDBWarning)
                        continue
                    if last_atom_added is None:
                        warnings.warn('Orphaned ANISOU record. Poorly '
                                      'formatted PDB file', PDBWarning)
                        continue
                    la = last_atom_added
                    if (la.name != aname or la.number != atnum or
                            la.altloc != altloc or la.residue.name != rname or
                            la.residue.chain != chain or
                            la.residue.insertion_code != icode):
                        warnings.warn('ANISOU record does not match previous '
                                      'atom', PDBWarning)
                        continue
                    la.anisou = np.array([u11/1e4, u22/1e4, u33/1e4,
                                          u12/1e4, u13/1e4, u23/1e4])
                elif rec.strip() == 'TER':
                    if modelno == 1: last_atom.residue.ter = True
                elif rec == 'ENDMDL':
                    # End the current model
                    if len(struct.atoms) == 0:
                        raise PDBError('MODEL ended before any atoms read in')
                    modelno += 1
                    if len(struct.atoms)*3 != len(coordinates):
                        raise PDBError(
                                'Inconsistent atom numbers in some PDB models')
                    all_coordinates.append(coordinates)
                    atomno = 0
                    coordinates = []
                    resend = 26
                    atom_overflow = False
                elif rec == 'MODEL ':
                    if modelno == 1 and len(struct.atoms) == 0: continue
                    if len(coordinates) > 0:
                        if len(struct.atoms)*3 != len(coordinates):
                            raise PDBError('Inconsistent atom numbers in '
                                           'some PDB models')
                        warnings.warn('MODEL not explicitly ended', PDBWarning)
                        all_coordinates.append(coordinates)
                        coordinates = []
                    modelno += 1
                    atomno = 0
                    resend = 26
                    atom_overflow = False
                elif rec == 'CRYST1':
                    a = float(line[6:15])
                    b = float(line[15:24])
                    c = float(line[24:33])
                    try:
                        A = float(line[33:40])
                        B = float(line[40:47])
                        C = float(line[47:54])
                    except (IndexError, ValueError):
                        A = B = C = 90.0
                    struct.box = [a, b, c, A, B, C]
                    struct.space_group = line[55:66].strip()
                elif rec == 'EXPDTA':
                    struct.experimental = line[6:].strip()
                elif rec == 'AUTHOR':
                    struct.authors += line[10:].strip()
                elif rec == 'JRNL  ':
                    part = line[12:16]
                    if part == 'AUTH':
                        struct.journal_authors += line[19:].strip()
                    elif part == 'TITL':
                        struct.title += ' %s' % line[19:].strip()
                    elif part == 'REF ':
                        struct.journal += ' %s' % line[19:47].strip()
                        if not line[16:18].strip():
                            struct.volume = line[51:55].strip()
                            struct.page = line[56:61].strip()
                            try:
                                struct.year = int(line[62:66])
                            except ValueError:
                                # Shouldn't happen, but don't throw a fit
                                pass
                    elif part == 'PMID':
                        struct.pmid = line[19:].strip()
                    elif part == 'DOI ':
                        struct.doi = line[19:].strip()
                elif rec == 'KEYWDS':
                    struct.keywords += '%s,' % line[10:]
                elif rec == 'REMARK':
                    if line[6:10] == ' 900':
                        # Related entries
                        rematch = PDBFile._relatere.match(line[11:])
                        if rematch:
                            struct.related_entries.append(rematch.groups())
                    elif line[6:10] == '   2':
                        # Resolution
                        if not line[11:22].strip(): continue
                        if struct.resolution is not None:
                            # Skip over comments
                            continue
                        if line[11:22] !=  'RESOLUTION.':
                            warnings.warn('Unrecognized RESOLUTION record in '
                                          'PDB file: %s' % line.strip())
                            continue
                        if line[23:38] == 'NOT APPLICABLE.':
                            # Not a diffraction experiment
                            continue
                        try:
                            struct.resolution = float(line[23:30])
                        except ValueError:
                            warnings.warn('Trouble converting resolution (%s) '
                                          'to float' % line[23:30])
                elif rec == 'CONECT' and not skip_bonds:
                    b = int(line[6:11])
                    try:
                        i = int(line[11:16])
                    except ValueError:
                        warnings.warn('Corrupt CONECT record', PDBWarning)
                        continue
                    # last 3 integers are optional and may not exist
                    j = line[16:21].strip()
                    k = line[21:26].strip()
                    l = line[26:31].strip()
                    if b in altloc_ids:
                        continue # Do not handle altloc bonds yet.
                    origin = _find_atom_index(struct, b)
                    if origin is None:
                        warnings.warn('CONECT record references non-existent '
                                      'origin atom %d' % b, PDBWarning)
                        continue # pragma: no cover
                    if i not in altloc_ids:
                        partner = _find_atom_index(struct, i)
                        if partner is None:
                            warnings.warn('CONECT record references non-existent '
                                          'destination atom %d' % i, PDBWarning)
                        elif partner not in origin.bond_partners:
                            struct.bonds.append(Bond(origin, partner))
                    # Other atoms are optional, so loop through the
                    # possibilities and bond them if they're set
                    for i in (j, k, l):
                        if not i: continue
                        i = int(i)
                        if i in altloc_ids: continue
                        partner = _find_atom_index(struct, i)
                        if partner is None:
                            warnings.warn('CONECT record references non-'
                                          'existent destination atom %d ' % i,
                                          PDBWarning)
                        elif partner not in origin.bond_partners:
                            struct.bonds.append(Bond(origin, partner))
        finally:
            # Make sure our file is closed if we opened it
            if own_handle: fileobj.close()

        # Assign bonds based on standard templates and simple distances
        if not skip_bonds:
            struct.assign_bonds()

        # Post-process some of the metadata to make it more reader-friendly
        struct.keywords = [s.strip() for s in struct.keywords.split(',')
                                            if s.strip()]
        struct.journal = struct.journal.strip()
        struct.title = struct.title.strip()

        struct.unchange()
        if coordinates:
            if len(coordinates) != 3*len(struct.atoms):
                raise PDBError('bad number of atoms in some PDB models')
            all_coordinates.append(coordinates)
        struct._coordinates = np.array(all_coordinates).reshape(
                        (-1, len(struct.atoms), 3))
        # process symmetry lines
        if _symmetry_lines:
            data = []
            for line in _symmetry_lines:
                if line.strip().startswith('REMARK 290   SMTRY'):
                    data.append(line.split()[4:])
            tensor = np.asarray(data, dtype='f8')
            struct.symmetry = Symmetry(tensor)
        return struct

    #===================================================

    @staticmethod
    def write(struct, dest, renumber=True, coordinates=None, altlocs='all',
              write_anisou=False, charmm=False, standard_resnames=False):
        """ Write a PDB file from a Structure instance

        Parameters
        ----------
        struct : :class:`Structure`
            The structure from which to write the PDB file
        dest : str or file-like
            Either a file name or a file-like object containing a `write`
            method to which to write the PDB file. If it is a filename that
            ends with .gz or .bz2, a compressed version will be written using
            either gzip or bzip2, respectively.
        renumber : bool, optional, default True
            If True, renumber the atoms and residues sequentially as they are
            stored in the structure.  If False, use the original numbering if
            it was assigned previously.
        coordinates : array-like of float, optional
            If provided, these coordinates will be written to the PDB file
            instead of the coordinates stored in the structure. These
            coordinates should line up with the atom order in the structure
            (not necessarily the order of the "original" PDB file if they
            differ)
        altlocs : str, optional, default 'all'
            Keyword controlling which alternate locations are printed to the
            resulting PDB file. Allowable options are:

                - 'all' : print all alternate locations
                - 'first' : print only the first alternate locations
                - 'occupancy' : print the one with the largest occupancy. If two
                  conformers have the same occupancy, the first one to occur is
                  printed

            Input is case-insensitive, and partial strings are permitted as long
            as it is a substring of one of the above options that uniquely
            identifies the choice.
        write_anisou : bool, optional, default False
            If True, an ANISOU record is written for every atom that has one. If
            False, ANISOU records are not written.
        charmm : bool, optional, default False
            If True, SEGID will be written in columns 73 to 76 of the PDB file
            in the typical CHARMM-style PDB output. This will be omitted for any
            atom that does not contain a SEGID identifier.
        standard_resnames : bool, optional, default False
            If True, common aliases for various amino and nucleic acid residues
            will be converted into the PDB-standard values.

        Notes
        -----
        If multiple coordinate frames are present, these will be written as
        separate models (but only the unit cell from the first model will be
        written, as the PDB standard dictates that only one set of unit cells
        shall be present).
        """
        if altlocs.lower() == 'all'[:len(altlocs)]:
            altlocs = 'all'
        elif altlocs.lower() == 'first'[:len(altlocs)]:
            altlocs = 'first'
        elif altlocs.lower() == 'occupancy'[:len(altlocs)]:
            altlocs = 'occupancy'
        else:
            raise ValueError("Illegal value of occupancy [%s]; expected 'all', "
                             "'first', or 'occupancy'" % altlocs)
        own_handle = False
        if not hasattr(dest, 'write'):
            dest = genopen(dest, 'w')
            own_handle = True
        if charmm:
            atomrec = ('ATOM  %5d %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f'
                       '%6.2f      %-4s%2s%-2s\n')
            hetatomrec = atomrec.replace('ATOM  ', 'HETATM')
            anisourec = ('ANISOU%5d %-4s%1s%-4s%1s%4d%1s %7d%7d%7d%7d%7d%7d'
                         '      %2s%-2s\n')
            terrec = ('TER   %5d      %-4s%1s%4d\n')
            reslen = 4
        else:
            atomrec = ('ATOM  %5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f'
                       '%6.2f      %-4s%2s%-2s\n')
            anisourec = ('ANISOU%5d %-4s%1s%-3s %1s%4d%1s %7d%7d%7d%7d%7d%7d'
                         '      %2s%-2s\n')
            terrec = ('TER   %5d      %-3s %1s%4d\n')
            reslen = 3
        hetatomrec = atomrec.replace('ATOM  ', 'HETATM')
        if struct.box is not None:
            dest.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4s\n' % (
                    struct.box[0], struct.box[1], struct.box[2], struct.box[3],
                    struct.box[4], struct.box[5], struct.space_group, ''))
        if struct.symmetry is not None:
            fmt = '%d%4d%10.6f%10.6f%10.6f%15.5f\n'
            for index, arr in enumerate(struct.symmetry.data):
                arr_list = [1 + index % 3, 1 + index//3] + arr.tolist()
                symm_line = "REMARK 290   SMTRY" + fmt % tuple(arr_list)
                dest.write(symm_line)
        if coordinates is not None:
            coords = np.array(coordinates, copy=False, subok=True)
            try:
                coords = coords.reshape((-1, len(struct.atoms), 3))
            except ValueError:
                raise TypeError("Coordinates has unexpected shape")
        else:
            coords = struct.get_coordinates('all')
        # Create a function to process each atom and return which one we want
        # to print, based on our alternate location choice
        if altlocs == 'all':
            def print_atoms(atom, coords):
                return atom, atom.other_locations, coords[atom.idx]
        elif altlocs == 'first':
            def print_atoms(atom, coords):
                return atom, dict(), coords[atom.idx]
        elif altlocs == 'occupancy':
            def print_atoms(atom, coords):
                occ = atom.occupancy
                a = atom
                for key, item in iteritems(atom.other_locations):
                    if item.occupancy > occ:
                        occ = item.occupancy
                        a = item
                return a, dict(), [a.xx, a.xy, a.xz]
        else:
            assert False, 'Should not be here'
        if standard_resnames:
            standardize = lambda x: _standardize_resname(x)[:reslen]
        else:
            standardize = lambda x: (x[:reslen], _is_hetatm(x))
        nmore = 0 # how many *extra* atoms have been added?
        last_number = 0
        last_rnumber = 0
        for model, coord in enumerate(coords):
            if coords.shape[0] > 1:
                dest.write('MODEL      %5d\n' % (model+1))
            for res in struct.residues:
                if renumber:
                    atoms = res.atoms
                else:
                    atoms = sorted(res.atoms, key=lambda atom: atom.number)
                if charmm:
                    segid = (res.segid or res.chain)[:4]
                else:
                    segid = ''
                for atom in atoms:
                    pa, others, (x, y, z) = print_atoms(atom, coord)
                    # Figure out the serial numbers we want to print
                    if renumber:
                        anum = (atom.idx + 1 + nmore)
                        rnum = (res.idx + 1)
                    else:
                        anum = (pa.number or last_number + 1)
                        rnum = (atom.residue.number or last_rnumber + 1)
                    anum = anum - anum // 100000 * 100000
                    rnum = rnum - rnum // 10000 * 10000
                    last_number = anum
                    last_rnumber = rnum
                    # Do any necessary name munging to respect the PDB spec
                    if len(pa.name) < 4 and len(Element[pa.atomic_number]) != 2:
                        aname = ' %-3s' % pa.name
                    else:
                        aname = pa.name[:4]
                    resname, hetatom = standardize(res.name)
                    if hetatom:
                        rec = hetatomrec
                    else:
                        rec = atomrec
                    dest.write(rec % (anum , aname, pa.altloc, resname,
                               res.chain[:1], rnum, res.insertion_code[:1],
                               x, y, z, pa.occupancy, pa.bfactor, segid,
                               Element[pa.atomic_number].upper(), ''))
                    if write_anisou and pa.anisou is not None:
                        anisou = [int(ani*1e4) for ani in pa.anisou]
                        dest.write(anisourec % (anum, aname, pa.altloc,
                                   resname, res.chain[:1], rnum,
                                   res.insertion_code[:1], anisou[0], anisou[1],
                                   anisou[2], anisou[3], anisou[4], anisou[5],
                                   Element[pa.atomic_number].upper(), ''))
                    for key in sorted(others.keys()):
                        oatom = others[key]
                        x, y, z = oatom.xx, oatom.xy, oatom.xz
                        if renumber:
                            nmore += 1
                            anum = (pa.idx + 1 + nmore)
                        else:
                            anum = oatom.number or last_number + 1
                        anum = anum - anum // 100000 * 100000
                        last_number = anum
                        if (len(oatom.name) < 4 and
                                len(Element[oatom.atomic_number]) != 2):
                            aname = ' %-3s' % oatom.name
                        else:
                            aname = oatom.name[:4]
                        dest.write(rec % (anum, aname, key, resname,
                                   res.chain[:1], rnum, res.insertion_code[:1],
                                   x, y, z, oatom.occupancy, oatom.bfactor, segid,
                                   Element[oatom.atomic_number].upper(), ''))
                        if write_anisou and oatom.anisou is not None:
                            anisou = [int(ani*1e4) for ani in oatom.anisou]
                            el = Element[oatom.atomic_number].upper()
                            dest.write(anisourec % (anum, aname,
                                oatom.altloc[:1], resname, res.chain[:1],
                                rnum, res.insertion_code[:1], anisou[0],
                                anisou[1], anisou[2], anisou[3],
                                anisou[4], anisou[5], el, ''))
                if res.ter or (len(struct.bonds) > 0 and _needs_ter_card(res)):
                    dest.write(terrec % (anum+1, resname, res.chain, rnum))
                    if renumber:
                        nmore += 1
                    else:
                        last_number += 1
            if coords.shape[0] > 1:
                dest.write('ENDMDL\n')

        dest.write("%-80s\n" % "END")
        if own_handle:
            dest.close()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

@add_metaclass(FileFormatType)
class CIFFile(object):
    """ Standard PDBx/mmCIF file format parser and writer """
    #===================================================

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as a PDBx/mmCIF file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is a PDBx/mmCIF file
        """
        f = genopen(filename)
        try:
            for line in f:
                if line.startswith('#'): continue
                if line[:5] == 'data_' and len(line.split()) == 1:
                    return True
                else:
                    return False
            return False
        finally:
            f.close()

    #===================================================

    @staticmethod
    def download(pdb_id, timeout=10, saveto=None):
        """
        Goes to the wwPDB website and downloads the requested PDBx/mmCIF,
        loading it as a :class:`Structure` instance

        Parameters
        ----------
        pdb_id : str
            The 4-letter PDB ID to try and download from the RCSB PDB database
        timeout : float, optional
            The number of seconds to wait before raising a timeout error.
            Default is 10 seconds
        saveto : str, optional
            If provided, this will be treated as a file name to which the PDB
            file will be saved. If None (default), no CIF file will be written.
            This will be a verbatim copy of the downloaded CIF file, unlike the
            somewhat-stripped version you would get by using
            :meth:`Structure.write_cif <parmed.structure.Structure.write_cif>`

        Returns
        -------
        struct : :class:`Structure <parmed.structure.Structure>`
            Structure instance populated by the requested PDBx/mmCIF

        Raises
        ------
        socket.timeout if the connection times out while trying to contact the
        FTP server

        IOError if there is a problem retrieving the requested PDB

        ImportError if the gzip module is not available

        TypeError if pdb_id is not a 4-character string
        """
        import gzip
        if not isinstance(pdb_id, string_types) or len(pdb_id) != 4:
            raise ValueError('pdb_id must be the 4-letter PDB code')

        pdb_id = pdb_id.lower()
        ftp = ftplib.FTP('ftp.wwpdb.org', timeout=timeout)
        ftp.login()
        fileobj = io.BytesIO()
        try:
            ftp.retrbinary('RETR /pub/pdb/data/structures/divided/mmCIF/'
                           '%s/%s.cif.gz' % (pdb_id[1:3], pdb_id),
                           fileobj.write)
        except ftplib.all_errors as err:
            raise IOError('Could not retrieve PDB ID %s; %s' % (pdb_id, err))
        finally:
            ftp.close()
        fileobj.seek(0)
        if PY3:
            fileobj = io.TextIOWrapper(gzip.GzipFile(fileobj=fileobj, mode='r'))
        else:
            fileobj = gzip.GzipFile(fileobj=fileobj, mode='r')
        if saveto is not None:
            with closing(genopen(saveto, 'w')) as f:
                f.write(fileobj.read())
            fileobj.seek(0)
        return CIFFile.parse(fileobj)

    #===================================================

    @staticmethod
    def parse(filename, skip_bonds=False):
        """
        Read a PDBx or mmCIF file and return a populated `Structure` class

        Parameters
        ----------
        filename : ``str or file-like``
            Name of PDB file to read, or a file-like object that can iterate
            over the lines of a PDB. Compressed file names can be specified and
            are determined by file-name extension (e.g., file.pdb.gz,
            file.pdb.bz2)
        skip_bonds : bool, optional
            If True, skip trying to assign bonds. This can save substantial time
            when parsing large files with non-standard residue names. However,
            no bonds are assigned. This is OK if, for instance, the CIF file is
            being parsed simply for its coordinates. Default is False.

        Metadata
        --------
        The PDB parser also adds metadata to the returned Structure object that
        may be present in the PDB file

        experimental : ``str``
            EXPDTA record
        journal : ``str``
            JRNL record
        authors : ``str``
            AUTHOR records
        keywords : ``str``
            KEYWDS records
        doi : ``str``
            DOI from the JRNL record
        pmid : ``str``
            PMID from the JRNL record
        journal_authors : ``str``
            Author info from the JRNL record
        volume : ``str``
            Volume of the published article from the JRNL record
        page : ``str``
            Page of the published article from the JRNL record
        title : ``str``
            TITL section of the JRNL record
        year : ``str``
            Year that the article was published, from the JRNL record
        resolution : ``float``
            The X-RAY resolution in Angstroms, or None if not found
        related_entries : ``list of (str, str)``
            List of entries in other databases

        Returns
        -------
        structure1 [, structure2 [, structure3 [, ...] ] ]

        structure# : :class:`Structure`
            The Structure object initialized with all of the information from
            the PDBx/mmCIF file.  No bonds or other topological features are
            added by default. If multiple structures are defined in the CIF
            file, multiple Structure instances will be returned as a tuple.

        Raises
        ------
        ValueError if the file severely violates the PDB format specification.
        If this occurs, check the formatting on each line and make sure it
        matches the others.
        """
        if isinstance(filename, string_types):
            own_handle = True
            fileobj = genopen(filename, 'r')
        else:
            own_handle = False
            fileobj = filename

        try:
            cifobj = PdbxReader(fileobj)
            data = []
            cifobj.read(data)
        finally:
            if own_handle: fileobj.close()

        structures = []
        for cont in data:
            struct = Structure()
            structures.append(struct)
            # Add metadata fields
            struct.experimental = struct.journal = struct.authors = ''
            struct.keywords = struct.doi = struct.pmid = ''
            struct.journal_authors = struct.volume_page = struct.title = ''
            struct.year = struct.resolution = None
            struct.related_entries = []

            # Now we have the data. First get the metadata if it exists
            exptl = cont.getObj('exptl')
            if exptl is not None:
                struct.experimental = exptl.getValue('method')
            auth = cont.getObj('audit_author')
            if auth is not None:
                nameidx = auth.getAttributeIndex('name')
                if nameidx != -1:
                    struct.authors = ', '.join([t[nameidx] for t in
                                               auth.getRowList()])
            reflns = cont.getObj('reflns')
            if reflns is not None:
                res = reflns.getValue('d_resolution_high')
                if res != '?':
                    try:
                        struct.resolution = float(res)
                    except ValueError:
                        warnings.warn('Could not convert resolution (%s) to '
                                      'float' % res)
            cite = cont.getObj('citation_author')
            if cite is not None:
                nameidx = cite.getAttributeIndex('name')
                if nameidx != -1:
                    journal_authors = []
                    for i in range(cite.getRowCount()):
                        a = cite.getRow(i)[nameidx]
                        if a not in journal_authors:
                            journal_authors.append(a)
                    struct.journal_authors = ', '.join(journal_authors)
            cite = cont.getObj('citation')
            if cite is not None:
                doiid = cite.getAttributeIndex('pdbx_database_id_DOI')
                pmiid = cite.getAttributeIndex('pdbx_database_id_PubMed')
                titlid = cite.getAttributeIndex('title')
                yearid = cite.getAttributeIndex('year')
                pageid = cite.getAttributeIndex('page_first')
                jrnlid = cite.getAttributeIndex('journal_abbrev')
                volid = cite.getAttributeIndex('journal_volume')
                rows = cite.getRowList()
                if doiid != -1:
                    struct.doi = ', '.join([row[doiid] for row in rows
                                                if row[doiid] != '?'])
                if pmiid != -1:
                    struct.pmid = ', '.join([row[pmiid] for row in rows
                                                if row[pmiid] != '?'])
                if titlid != -1:
                    struct.title = '; '.join([row[titlid] for row in rows])
                if yearid != -1:
                    struct.year = ', '.join([row[yearid] for row in rows])
                if pageid != -1:
                    struct.page = ', '.join([row[pageid] for row in rows])
                if jrnlid != -1:
                    struct.journal = '; '.join([row[jrnlid] for row in rows])
                if volid != -1:
                    struct.volume = ', '.join([row[volid] for row in rows])
            keywds = cont.getObj('struct_keywords')
            if keywds is not None:
                textid = keywds.getAttributeIndex('text')
                if textid != -1:
                    rows = keywds.getRowList()
                    struct.keywords = ', '.join([row[textid] for row in rows])
                    struct.keywords = [key.strip() for key in
                            struct.keywords.split(',') if key.strip()]
            dbase = cont.getObj('pdbx_database_related')
            if dbase is not None:
                dbid = dbase.getAttributeIndex('db_id')
                nameid = dbase.getAttributeIndex('db_name')
                if dbid != -1 and nameid != -1:
                    rows = dbase.getRowList()
                    struct.related_entries = [(r[dbid],r[nameid]) for r in rows]
            # Now go through all of the atoms. Any items that do *not* exist are
            # given an index of -1. So we append an empty string on the end of
            # each row so that the default value for any un-specified value is
            # the empty string. This avoids needing any conditionals inside the
            # loop
            atoms = cont.getObj('atom_site')
            atnumid = atoms.getAttributeIndex('id')
            elemid = atoms.getAttributeIndex('type_symbol')
            atnameid = atoms.getAttributeIndex('auth_atom_id')
            altlocid = atoms.getAttributeIndex('label_alt_id')
            resnameid = atoms.getAttributeIndex('auth_comp_id')
            chainid = atoms.getAttributeIndex('auth_asym_id')
            resnumid = atoms.getAttributeIndex('auth_seq_id')
            inscodeid = atoms.getAttributeIndex('pdbx_PDB_ins_code')
            xid = atoms.getAttributeIndex('Cartn_x')
            yid = atoms.getAttributeIndex('Cartn_y')
            zid = atoms.getAttributeIndex('Cartn_z')
            occupid = atoms.getAttributeIndex('occupancy')
            bfactorid = atoms.getAttributeIndex('B_iso_or_equiv')
            modelid = atoms.getAttributeIndex('pdbx_PDB_model_num')
            origmodel = None
            lastmodel = None
            all_coords = []
            xyz = []
            atommap = dict()
            last_atom = Atom()
            for i in range(atoms.getRowCount()):
                row = atoms.getRow(i) + ['']
                atnum = int(row[atnumid])
                elem = row[elemid]
                atname = row[atnameid]
                altloc = row[altlocid]
                if altloc == '.': altloc = ''
                resname = row[resnameid]
                chain = row[chainid]
                resnum = int(row[resnumid])
                inscode = row[inscodeid]
                if inscode in '?.': inscode = ''
                model = int(row[modelid])
                if origmodel is None:
                    origmodel = lastmodel = model
                x, y, z = float(row[xid]), float(row[yid]), float(row[zid])
                occup = float(row[occupid])
                bfactor = float(row[bfactorid])
                # Try to figure out the element
                elem = '%-2s' % elem # Make sure we have at least 2 characters
                if elem[0] == ' ': elem = elem[1] + ' '
                try:
                    atsym = (elem[0] + elem[1].lower()).strip()
                    atomic_number = AtomicNum[atsym]
                    mass = Mass[atsym]
                except KeyError:
                    # Now try based on the atom name... but don't try too hard
                    # (e.g., don't try to differentiate b/w Ca and C)
                    try:
                        atomic_number = AtomicNum[atname.strip()[0].upper()]
                        mass = Mass[atname.strip()[0].upper()]
                    except KeyError:
                        try:
                            sym = atname.strip()[:2]
                            sym = '%s%s' % (sym[0].upper(), sym[1].lower())
                            atomic_number = AtomicNum[sym]
                            mass = Mass[sym]
                        except KeyError:
                            atomic_number = 0 # give up
                            mass = 0.0
                if atname.startswith('EP') or atname.startswith('LP'):
                    atom = ExtraPoint(atomic_number=atomic_number, name=atname,
                                mass=mass, occupancy=occup, bfactor=bfactor,
                                altloc=altloc, number=atnum)
                else:
                    atom = Atom(atomic_number=atomic_number, name=atname,
                                mass=mass, occupancy=occup, bfactor=bfactor,
                                altloc=altloc, number=atnum)
                atom.xx, atom.xy, atom.xz = x, y, z
                if (_compare_atoms(last_atom, atom, resname, resnum,
                                   chain, '', inscode)
                        and altloc):
                    atom.residue = last_atom.residue
                    last_atom.other_locations[altloc] = atom
                else:
                    if model == origmodel:
                        # Only add the atoms once
                        struct.add_atom(atom, resname, resnum, chain, inscode)
                    last_atom = atom
                    if model == lastmodel:
                        xyz.extend([x, y, z])
                    else:
                        if all_coords and len(xyz) != len(all_coords[-1]):
                            raise ValueError('All frames must have same number '
                                             'of atoms')
                        all_coords.append(xyz)
                        xyz = [x, y, z]
                        lastmodel = model
                # Keep a mapping in case we need to go back and add attributes,
                # like anisotropic b-factors
                if model == origmodel:
                    key = (resnum,resname,inscode,chain,atnum,altloc,atname)
                    atommap[key] = atom
            # Check for unit cell parameters
            cell = cont.getObj('cell')
            if cell is not None:
                aid = cell.getAttributeIndex('length_a')
                bid = cell.getAttributeIndex('length_b')
                cid = cell.getAttributeIndex('length_c')
                alphaid = cell.getAttributeIndex('angle_alpha')
                betaid = cell.getAttributeIndex('angle_beta')
                gammaid = cell.getAttributeIndex('angle_gamma')
                row = cell.getRow(0)
                struct.box = np.array(
                        [float(row[aid]), float(row[bid]), float(row[cid]),
                         float(row[alphaid]), float(row[betaid]),
                         float(row[gammaid])]
                )
            symmetry = cont.getObj('symmetry')
            if symmetry is not None:
                spaceid = symmetry.getAttributeIndex('space_group_name_H-M')
                row = symmetry.getRow(0)
                if spaceid != -1:
                    struct.space_group = row[spaceid]
            # Check for anisotropic B-factors
            anisou = cont.getObj('atom_site_anisotrop')
            if anisou is not None:
                atnumid = anisou.getAttributeIndex('id')
                atnameid = anisou.getAttributeIndex('pdbx_auth_atom_id')
                altlocid = anisou.getAttributeIndex('pdbx_label_alt_id')
                resnameid = anisou.getAttributeIndex('pdbx_auth_comp_id')
                chainid = anisou.getAttributeIndex('pdbx_auth_asym_id')
                resnumid = anisou.getAttributeIndex('pdbx_auth_seq_id')
                inscodeid = anisou.getAttributeIndex('pdbx_PDB_ins_code')
                u11id = anisou.getAttributeIndex('U[1][1]')
                u22id = anisou.getAttributeIndex('U[2][2]')
                u33id = anisou.getAttributeIndex('U[3][3]')
                u12id = anisou.getAttributeIndex('U[1][2]')
                u13id = anisou.getAttributeIndex('U[1][3]')
                u23id = anisou.getAttributeIndex('U[2][3]')
                if -1 in (atnumid, atnameid, altlocid, resnameid, chainid,
                          resnumid, u11id, u22id, u33id, u12id, u13id, u23id):
                    warnings.warn('Incomplete anisotropic B-factor CIF '
                                  'section. Skipping', PDBWarning)
                else:
                    try:
                        for i in range(anisou.getRowCount()):
                            row = anisou.getRow(i) + ['']
                            atnum = int(row[atnumid])
                            atname = row[atnameid]
                            altloc = row[altlocid]
                            resname = row[resnameid]
                            chain = row[chainid]
                            resnum = int(row[resnumid])
                            inscode = row[inscodeid]
                            u11 = float(row[u11id])
                            u22 = float(row[u22id])
                            u33 = float(row[u33id])
                            u12 = float(row[u12id])
                            u13 = float(row[u13id])
                            u23 = float(row[u23id])
                            if altloc == '.': altloc = ''
                            if inscode in '?.': inscode = ''
                            key = (resnum, resname, inscode, chain, atnum,
                                   altloc, atname)
                            atommap[key].anisou = np.array(
                                    [u11, u22, u33, u12, u13, u23]
                            )
                    except (ValueError, KeyError):
                        # If at least one went wrong, set them all to None
                        for key, atom in iteritems(atommap):
                            atom.anisou = None
                        warnings.warn('Problem processing anisotropic '
                                      'B-factors. Skipping', PDBWarning)
            if xyz:
                if len(xyz) != len(struct.atoms) * 3:
                    raise ValueError('Corrupt CIF; all models must have the '
                                     'same atoms')
                all_coords.append(xyz)
            if all_coords:
                struct._coordinates = np.array(all_coords).reshape(
                            (-1, len(struct.atoms), 3))

        # Make sure we assign bonds for all of the structures we parsed
        if not skip_bonds:
            for struct in structures:
                struct.assign_bonds()
        # Build the return value
        if len(structures) == 1:
            return structures[0]
        return tuple(structures)

    #===================================================

    @staticmethod
    def write(struct, dest, renumber=True, coordinates=None,
              altlocs='all', write_anisou=False, standard_resnames=False):
        """
        Write a PDB file from the current Structure instance

        Parameters
        ----------
        struct : :class:`Structure`
            The structure from which to write the PDBx/mmCIF file
        dest : ``str or file-like``
            Either a file name or a file-like object containing a `write`
            method to which to write the PDB file. If it is a filename that
            ends with .gz or .bz2, a compressed version will be written using
            either gzip or bzip2, respectively.
        renumber : ``bool``
            If True, renumber the atoms and residues sequentially as they are
            stored in the structure.  If False, use the original numbering if
            it was assigned previously
        coordinates : ``array-like of float``
            If provided, these coordinates will be written to the PDB file
            instead of the coordinates stored in the structure. These
            coordinates should line up with the atom order in the structure
            (not necessarily the order of the "original" PDB file if they
            differ)
        altlocs : ``str``
            Keyword controlling which alternate locations are printed to the
            resulting PDB file. Allowable options are:

                - 'all' : (default) print all alternate locations
                - 'first' : print only the first alternate locations
                - 'occupancy' : print the one with the largest occupancy. If two
                  conformers have the same occupancy, the first one to occur is
                  printed

            Input is case-insensitive, and partial strings are permitted as long
            as it is a substring of one of the above options that uniquely
            identifies the choice.
        write_anisou : ``bool``
            If True, an ANISOU record is written for every atom that has one. If
            False, ANISOU records are not written
        standard_resnames : bool, optional
            If True, common aliases for various amino and nucleic acid residues
            will be converted into the PDB-standard values. Default is False

        Notes
        -----
        If multiple coordinate frames are present, these will be written as
        separate models (but only the unit cell from the first model will be
        written, as the PDBx standard dictates that only one set of unit cells
        shall be present).
        """
        if altlocs.lower() == 'all'[:len(altlocs)]:
            altlocs = 'all'
        elif altlocs.lower() == 'first'[:len(altlocs)]:
            altlocs = 'first'
        elif altlocs.lower() == 'occupancy'[:len(altlocs)]:
            altlocs = 'occupancy'
        else:
            raise ValueError("Illegal value of occupancy [%s]; expected 'all', "
                             "'first', or 'occupancy'" % altlocs)
        own_handle = False
        if not hasattr(dest, 'write'):
            dest = genopen(dest, 'w')
            own_handle = True
        # Make the main container
        cont = containers.DataContainer('cell')
        # Add cell info if applicable
        if struct.box is not None:
            cell = containers.DataCategory('cell')
            cell.appendAttribute('length_a')
            cell.appendAttribute('length_b')
            cell.appendAttribute('length_c')
            cell.appendAttribute('angle_alpha')
            cell.appendAttribute('angle_beta')
            cell.appendAttribute('angle_gamma')
            cell.append(struct.box[:])
            cont.append(cell)
        # symmetry
        sym = containers.DataCategory('symmetry')
        sym.appendAttribute('space_group_name_H-M')
        sym.append([struct.space_group])
        cont.append(sym)
        if coordinates is not None:
            coords = np.array(coordinates, copy=False, subok=True)
            try:
                coords = coords.reshape((-1, len(struct.atoms), 3))
            except ValueError:
                raise TypeError("Coordinates has unexpected shape")
        else:
            coords = struct.get_coordinates('all')
        # Create a function to process each atom and return which one we want
        # to print, based on our alternate location choice
        if altlocs == 'all':
            def print_atoms(atom, coords):
                return atom, atom.other_locations, coords[atom.idx]
        elif altlocs == 'first':
            def print_atoms(atom, coords):
                return atom, dict(), coords[atom.idx]
        elif altlocs == 'occupancy':
            def print_atoms(atom, coords):
                occ = atom.occupancy
                a = atom
                for key, item in iteritems(atom.other_locations):
                    if item.occupancy > occ:
                        occ = item.occupancy
                        a = item
                return a, dict(), [a.xx, a.xy, a.xz]
        else:
            assert False, 'Should not be here'
        if standard_resnames:
            standardize = lambda x: _standardize_resname(x)
        else:
            standardize = lambda x: (x, _is_hetatm(x))
        # Now add the atom section. Include all names that the CIF standard
        # usually includes, but put '?' in sections that contain data we don't
        # store in the Structure, Residue, or Atom classes
        cifatoms = containers.DataCategory('atom_site')
        cont.append(cifatoms)
        cifatoms.setAttributeNameList(
                ['group_PDB', 'id', 'type_symbol', 'label_atom_id',
                 'label_alt_id', 'label_comp_id', 'label_asym_id',
                 'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code',
                 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv',
                 'Cartn_x_esd', 'Cartn_y_esd', 'Cartn_z_esd', 'occupancy_esd',
                 'B_iso_or_equiv_esd', 'pdbx_formal_charge', 'auth_seq_id',
                 'auth_comp_id', 'auth_asym_id', 'auth_atom_id',
                 'pdbx_PDB_model_num']
        )
        write_anisou = write_anisou and any(atom.anisou is not None
                                            for atom in struct.atoms)
        if write_anisou:
            cifanisou = containers.DataCategory('atom_site_anisotrop')
            cont.append(cifanisou)
            cifanisou.setAttributeNameList(
                    ['id', 'type_symbol', 'pdbx_label_atom_id',
                    'pdbx_label_alt_id', 'pdbx_label_comp_id',
                    'pdbx_label_asym_id', 'pdbx_label_seq_id', 'U[1][1]',
                    'U[2][2]', 'U[3][3]', 'U[1][2]', 'U[1][3]', 'U[2][3]',
                    'U[1][1]_esd', 'U[2][2]_esd', 'U[3][3]_esd', 'U[1][2]_esd',
                    'U[1][3]_esd', 'U[2][3]_esd', 'pdbx_auth_seq_id',
                    'pdbx_auth_comp_id', 'pdbx_auth_asym_id',
                    'pdbx_auth_atom_id']
            )
        nmore = 0 # how many *extra* atoms have been added?
        last_number = 0
        last_rnumber = 0
        for model, coord in enumerate(coords):
            for res in struct.residues:
                if renumber:
                    atoms = res.atoms
                else:
                    atoms = sorted(res.atoms, key=lambda atom: atom.number)
                resname, hetatom = standardize(res.name)
                if hetatom:
                    atomrec = 'HETATM'
                else:
                    atomrec = 'ATOM  '
                for atom in atoms:
                    pa, others, (x, y, z) = print_atoms(atom, coord)
                    # Figure out the serial numbers we want to print
                    if renumber:
                        anum = (atom.idx + 1 + nmore)
                        rnum = (res.idx + 1)
                    else:
                        anum = (pa.number or last_number + 1)
                        rnum = (atom.residue.number or last_rnumber + 1)
                    last_number = anum
                    last_rnumber = rnum
                    cifatoms.append(
                            [atomrec, anum, Element[pa.atomic_number].upper(),
                             pa.name, pa.altloc, resname, res.chain, '?', rnum,
                             res.insertion_code, x, y, z, pa.occupancy,
                             pa.bfactor, '?', '?', '?', '?', '?', '', rnum,
                             resname, res.chain, pa.name, str(model+1)]
                    )
                    if write_anisou and pa.anisou is not None:
                        cifanisou.append(
                                [anum, Element[pa.atomic_number].upper(),
                                 pa.name, pa.altloc, resname, res.chain, rnum,
                                 pa.anisou[0], pa.anisou[1], pa.anisou[2],
                                 pa.anisou[3], pa.anisou[4], pa.anisou[5], '?',
                                 '?', '?', '?', '?', '?', rnum, resname,
                                 res.chain, pa.name]
                        )
                    for key in sorted(others.keys()):
                        oatom = others[key]
                        x, y, z = oatom.xx, oatom.xy, oatom.xz
                        if renumber:
                            nmore += 1
                            anum = (pa.idx + 1 + nmore)
                        else:
                            anum = oatom.number or last_number + 1
                        last_number = anum
                        el = Element[oatom.atomic_number].upper()
                        cifatoms.append(
                                [atomrec, anum, el, oatom.name, oatom.altloc,
                                 resname, res.chain, '?', rnum,
                                 res.insertion_code, x, y, z, oatom.occupancy,
                                 oatom.bfactor, '?', '?', '?', '?', '?', '',
                                 rnum, resname, res.chain, oatom.name, '1']
                        )
                        if write_anisou and oatom.anisou is not None:
                            cifanisou.append(
                                    [anum, Element[oatom.atomic_number].upper(),
                                     oatom.name, oatom.altloc, resname,
                                     res.chain, rnum, oatom.anisou[0],
                                     oatom.anisou[1], oatom.anisou[2],
                                     oatom.anisou[3], oatom.anisou[4],
                                     oatom.anisou[5], '?', '?', '?', '?', '?',
                                     '?', rnum, resname, res.chain, oatom.name]
                            )
        # Now write the PDBx file
        writer = PdbxWriter(dest)
        writer.write([cont])
        if own_handle:
            dest.close()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _find_atom_index(struct, idx):
    """
    Returns the atom with the given index in the structure. This is required
    because atom indices may not start from 1 and may contain gaps in a PDB
    file. This tries to find atoms quickly, assuming that indices *do* start
    from 1 and have no gaps. It then looks up or down, depending on whether we
    hit an index too high or too low. So it *assumes* that the sequence is
    monotonically increasing. If the atom can't be found, None is returned
    """
    idx0 = min(max(idx - 1, 0), len(struct.atoms)-1)
    if struct[idx0].number == idx:
        return struct[idx0]
    if struct[idx0].number < idx:
        idx0 += 1
        while idx0 < len(struct.atoms):
            if struct[idx0].number == idx:
                return struct[idx0]
            idx0 += 1
        return None # not found
    else:
        idx0 -= 1
        while idx0 > 0:
            if struct[idx0].number == idx:
                return struct[idx0]
            idx0 -= 1
        return None # not found

def _needs_ter_card(res):
    """ Determines if a TER card is needed by seeing if the residue is a
    polymeric residue that is *not* bonded to the next residue
    """
    # First see if it's in the list of standard biomolecular residues. If so,
    # and it has no tail, no TER is needed
    std_resname = _standardize_resname(res.name)[0]
    if std_resname in StandardBiomolecularResidues:
        is_std_res = True
        if StandardBiomolecularResidues[std_resname].tail is None:
            return False
    else:
        is_std_res = False
    my_res_idx = res.idx
    residxs = set()
    for atom in res.atoms:
        for bond in atom.bonds:
            residxs |= {bond.atom1.residue.idx, bond.atom2.residue.idx}
    if my_res_idx + 1 in residxs:
        return False # It's connected to next residue
    elif is_std_res:
        return True
    else:
        # Heuristic -- add a TER if it's bonded to the previous residue, which
        # indicates it's polymeric. Otherwise don't.
        return my_res_idx - 1 in residxs
