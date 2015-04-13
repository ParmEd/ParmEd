"""
This package contains classes responsible for reading and writing both PDB and
PDBx/mmCIF files.
"""
from chemistry.exceptions import PDBError, AnisouWarning, PDBWarning
from chemistry.formats.pdbx import PdbxReader, PdbxWriter, containers
from chemistry.formats.io import genopen, TextToBinaryFile
from chemistry.formats.registry import FileFormatType
from chemistry.periodic_table import AtomicNum, Mass, Element
from chemistry.structure import Structure
from chemistry.topologyobjects import Atom, ExtraPoint
import itertools
try:
    import numpy as np
    create_array = lambda x: np.array(x, dtype=np.float64)
except ImportError:
    create_array = lambda x: [float(v) for v in x]
import re
import warnings

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _compare_atoms(old_atom, new_atom, resname, resid, chain):
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

    Returns
    -------
    True if they are the same atom, False otherwise
    """
    if old_atom.name != new_atom.name: return False
    if old_atom.residue.name != resname: return False
    if old_atom.residue.number != resid: return False
    if old_atom.residue.chain != chain.strip(): return False
    return True

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PDBFile(object):
    """ Standard PDB file format parser and writer """
    __metaclass__ = FileFormatType

    #===================================================

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as a PDB file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is a PDB file
        """
        f = genopen(filename, 'r')
        lines = [f.readline().decode() for i in xrange(3)]
        f.close()

        for line in lines:
            if line[:6] in ('CRYST1', 'END   ', 'END', 'HEADER', 'NUMMDL',
                    'MASTER', 'ORIGXn', 'SCALEn', 'AUTHOR', 'CAVEAT', 'COMPND',
                    'EXPDTA', 'MDLTYP', 'KEYWDS', 'OBSLTE', 'SOURCE', 'SPLIT ',
                    'SPRSDE', 'TITLE ', 'ANISOU', 'ATOM  ', 'CISPEP', 'CONECT',
                    'DBREF ', 'HELIX ', 'HET   ', 'HETATM', 'LINK  ', 'MODRES',
                    'MTRIXn', 'REVDAT', 'SEQADV', 'SHEET ', 'SSBOND', 'FORMUL',
                    'HETNAM', 'HETSYN', 'SEQRES', 'SITE  ', 'ENDMDL', 'MODEL ',
                    'TER   ', 'TER', 'JRNL  ', 'REMARK'):
                continue
            return False
        return True

    #===================================================

    _relatere = re.compile(r'RELATED ID: *(\w+) *RELATED DB: *(\w+)', re.I)

    @staticmethod
    def parse(filename):
        """ Read a PDB file and return a populated `Structure` class

        Parameters
        ----------
        filename : str or file-like
            Name of the PDB file to read, or a file-like object that can iterate
            over the lines of a PDB. Compressed file names can be specified and
            are determined by file-name extension (e.g., file.pdb.gz,
            file.pdb.bz2)

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
        related_entries : ``list of (str, str)``
            List of entries in other databases
    
        Returns
        -------
        structure
        
        structure : :class:`Structure`
            The Structure object initialized with all of the information from
            the PDB file.  No bonds or other topological features are added by
            default.
    
        Notes
        -----
        The returned structure has an extra attribute, pdbxyz, that contains all
        of the coordinates for all of the frames in the PDB file as a list of
        NATOM*3 lists.
        """
        if isinstance(filename, basestring):
            own_handle = True
            fileobj = genopen(filename, 'r')
        else:
            own_handle = False
            fileobj = filename

        struct = Structure()
        # Add metadata fields
        struct.experimental = struct.journal = struct.authors = struct.keywords = ''
        struct.doi = struct.pmid = struct.journal_authors = struct.volume_page = ''
        struct.title = ''
        struct.year = None
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

        try:
            for line in fileobj:
                try:
                    line = line.decode('ascii')
                except AttributeError:
                    # Assume this is a string in Py3 which doesn't have 'decode'
                    pass
                rec = line[:6]
                if rec == 'ATOM  ' or rec == 'HETATM':
                    atomno += 1
                    atnum, atname, altloc = line[6:11], line[12:16], line[16]
                    resname, chain = line[17:21], line[21]
                    resid, inscode = line[22:resend], line[26]
                    x, y, z = line[30:38], line[38:46], line[47:54]
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
                        try:
                            atomic_number = AtomicNum[atname.strip()[0].upper()]
                            mass = Mass[atname.strip()[0].upper()]
                        except KeyError:
                            try:
                                sym = atname.strip()[:2]
                                sym = '%s%s' % (sym[0].upper(), sym[0].lower())
                                atomic_number = AtomicNum[sym]
                                mass = Mass[sym]
                            except KeyError:
                                atomic_number = 0 # give up
                                mass = 0.0
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
                            except ValueError, e:
                                if resid == '****':
                                    resid = None # Figure out by unique atoms
                                else:
                                    raise e
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
                                raise ValueError('Could not convert %s to int' %
                                                 atnum)
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
                    # name as # the 'last' residue. Do not worry about atom
                    # numbers going to *****'s, since that is >1M atoms.
                    if resid is None:
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
                    if segid: atom.segid = segid
                    if _compare_atoms(last_atom, atom, resname, resid, chain):
                        atom.residue = last_atom.residue
                        last_atom.other_locations[altloc] = atom
                        last_atom_added = atom
                        continue
                    last_atom = last_atom_added = atom
                    if modelno == 1:
                        struct.add_atom(atom, resname, resid, chain, inscode)
                    else:
                        try:
                            orig_atom = struct.atoms[atomno-1]
                        except IndexError:
                            raise PDBError('Atom %d differs in MODEL %d [%s %s '
                                           'vs. %s %s]' % (atomno, modelno,
                                           atom.residue.name, atom.name,
                                           resname, atname))
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
                                      'record', AnisouWarning)
                        continue # Skip the rest of this record
                    aname = line[12:16].strip()
                    altloc = line[16].strip()
                    rname = line[17:21].strip()
                    chain = line[21].strip()
                    try:
                        resid = int(line[22:26])
                    except ValueError:
                        warnings.warn('Problem parsing residue number from '
                                      'ANISOU record', AnisouWarning)
                        continue # Skip the rest of this record
                    icode = line[27].strip()
                    try:
                        u11 = int(line[28:35])
                        u22 = int(line[35:42])
                        u33 = int(line[42:49])
                        u12 = int(line[49:56])
                        u13 = int(line[56:63])
                        u23 = int(line[63:70])
                    except ValueError:
                        warnings.warn('Problem parsing anisotropic factors '
                                      'from ANISOU record', AnisouWarning)
                        continue
                    if last_atom_added is None:
                        warnings.warn('Orphaned ANISOU record. Poorly '
                                      'formatted PDB file', AnisouWarning)
                        continue
                    la = last_atom_added
                    if (la.name != aname or la.number != atnum or
                            la.altloc != altloc or la.residue.name != rname or
                            la.residue.chain != chain or
                            la.residue.insertion_code != icode):
                        warnings.warn('ANISOU record does not match previous '
                                      'atom', AnisouWarning)
                        continue
                    la.anisou = create_array([u11/1e4, u22/1e4, u33/1e4,
                                              u12/1e4, u13/1e4, u23/1e4])
                elif rec.strip() == 'TER':
                    if modelno == 1: last_atom.residue.ter = True
                elif rec == 'ENDMDL':
                    # End the current model
                    if len(struct.atoms) == 0:
                        raise PDBError('MODEL ended before any atoms read in')
                    modelno += 1
                    if len(struct.atoms)*3 != len(coordinates):
                        raise ValueError(
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
                            raise ValueError('Inconsistent atom numbers in '
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
                    try:
                        struct.space_group = line[55:66].strip()
                    except IndexError:
                        pass
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
                                pass
                    elif part == 'PMID':
                        struct.pmid = line[19:].strip()
                    elif part == 'DOI ':
                        struct.doi = line[19:].strip()
                elif rec == 'KEYWDS':
                    struct.keywords += '%s,' % line[10:]
                elif rec == 'REMARK' and line[6:10] == ' 900':
                    # Related entries
                    rematch = PDBFile._relatere.match(line[11:])
                    if rematch:
                        struct.related_entries.append(rematch.groups())
        finally:
            # Make sure our file is closed if we opened it
            if own_handle: fileobj.close()

        # Post-process some of the metadata to make it more reader-friendly
        struct.keywords = [s.strip() for s in struct.keywords.split(',')
                                            if s.strip()]
        struct.journal = struct.journal.strip()
        struct.title = struct.title.strip()

        struct.unchange()
        if coordinates:
            if len(coordinates) != 3*len(struct.atoms):
                raise ValueError('bad number of atoms in some PDB models')
            all_coordinates.append(coordinates)
        struct.pdbxyz = all_coordinates
        return struct

    #===================================================

    @staticmethod
    def write(struct, dest, renumber=True, coordinates=None, altlocs='all',
              write_anisou=False, charmm=False):
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
        renumber : bool, optional
            If True, renumber the atoms and residues sequentially as they are
            stored in the structure.  If False, use the original numbering if
            it was assigned previously. Default is True
        coordinates : array-like of float, optional
            If provided, these coordinates will be written to the PDB file
            instead of the coordinates stored in the structure. These
            coordinates should line up with the atom order in the structure
            (not necessarily the order of the "original" PDB file if they
            differ)
        altlocs : str, optional
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
        write_anisou : bool, optional
            If True, an ANISOU record is written for every atom that has one. If
            False, ANISOU records are not written. Default is False
        charmm : bool, optional
            If True, SEGID will be written in columns 73 to 76 of the PDB file
            in the typical CHARMM-style PDB output. This will be omitted for any
            atom that does not contain a SEGID identifier. Default is False
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
            dest = TextToBinaryFile(genopen(dest, 'w'))
            own_handle = True
        if charmm:
            atomrec = ('ATOM  %5d %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f'
                       '%6.2f      %-4s%2s%-2s\n')
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
        if struct.box is not None:
            dest.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4s\n' % (
                    struct.box[0], struct.box[1], struct.box[2], struct.box[3],
                    struct.box[4], struct.box[5], struct.space_group, ''))
        if coordinates is not None:
            try:
                crdsize = len(coordinates)
            except TypeError:
                raise TypeError("Cannot find length of coordinates")
            if crdsize == len(struct.atoms):
                try:
                    coords = coordinates.flatten()
                except AttributeError:
                    try:
                        coords = list(itertools.chain(*coordinates))
                    except TypeError:
                        raise TypeError("Unsupported coordinate dimensionality")
                if len(coords) != len(struct.atoms) * 3:
                    raise TypeError("Unsupported coordinate shape")
            elif crdsize == len(struct.atoms) * 3:
                coords = coordinates
            else:
                raise TypeError("Coordinates has unexpected shape")
        else:
            coords = [[a.xx, a.xy, a.xz] for a in struct.atoms]
            coords = list(itertools.chain(*coords))
        # Create a function to process each atom and return which one we want
        # to print, based on our alternate location choice
        if altlocs == 'all':
            def print_atoms(atom, coords):
                i3 = atom.idx * 3
                return atom, atom.other_locations, coords[i3:i3+3]
        elif altlocs == 'first':
            def print_atoms(atom, coords):
                i3 = atom.idx * 3
                return atom, dict(), coords[i3:i3+3]
        elif altlocs == 'occupancy':
            def print_atoms(atom, coords):
                occ = atom.occupancy
                a = atom
                for key, item in atom.other_locations.iteritems():
                    if item.occupancy > occ:
                        occ = item.occupancy
                        a = item
                return a, dict(), [a.xx, a.xy, a.xz]
        else:
            raise Exception("Should not be here!")
        nmore = 0 # how many *extra* atoms have been added?
        last_number = 0
        last_rnumber = 0
        for res in struct.residues:
            if renumber:
                atoms = res.atoms
            else:
                atoms = sorted(res.atoms, key=lambda atom: atom.number)
            for atom in atoms:
                pa, others, (x, y, z) = print_atoms(atom, coords)
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
                if charmm and hasattr(pa, 'segid'):
                    segid = pa.segid[:4]
                else:
                    segid = ''
                dest.write(atomrec % (anum , aname, pa.altloc,
                           res.name[:reslen], res.chain[:1], rnum,
                           res.insertion_code[:1], x, y, z, pa.occupancy,
                           pa.bfactor, segid,
                           Element[pa.atomic_number].upper(), ''))
                if write_anisou and pa.anisou is not None:
                    anisou = [int(ani*1e4) for ani in pa.anisou]
                    dest.write(anisourec % (anum, aname, pa.altloc,
                               res.name[:reslen], res.chain[:1], rnum,
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
                    if charmm and hasattr(oatom, 'segid'):
                        segid = oatom.segid[:4]
                    else:
                        segid = ''
                    dest.write(atomrec % (anum, aname, key, res.name[:reslen],
                               res.chain[:1], rnum, res.insertion_code[:1],
                               x, y, z, oatom.occupancy, oatom.bfactor, segid,
                               Element[oatom.atomic_number].upper(), ''))
                    if write_anisou and oatom.anisou is not None:
                        anisou = [int(ani*1e4) for ani in oatom.anisou]
                        dest.write(anisourec % (anum, aname,
                            oatom.altloc[:1], res.name[:reslen], res.chain[:1],
                            rnum, res.insertion_code[:1], anisou[0], anisou[1],
                            anisou[2], anisou[3], anisou[4], anisou[5],
                            Element[oatom.atomic_number].upper(), ''))
            if res.ter:
                dest.write(terrec % (anum+1, res.name[:reslen],
                                     res.chain, rnum))
                if renumber:
                    nmore += 1
                else:
                    last_number += 1

        dest.write("%-80s" % "END")
        if own_handle:
            dest.close()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CIFFile(object):
    """ Standard PDBx/mmCIF file format parser and writer """
    __metaclass__ = FileFormatType

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
                line = line.decode()
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
    def parse(filename):
        """
        Read a PDBx or mmCIF file and return a populated `Structure` class

        Parameters
        ----------
        filename : ``str or file-like``
            Name of PDB file to read, or a file-like object that can iterate
            over the lines of a PDB. Compressed file names can be specified and
            are determined by file-name extension (e.g., file.pdb.gz,
            file.pdb.bz2)

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

        Notes
        -----
        The returned structure has an extra attribute, pdbxyz, that contains all
        of the coordinates for all of the frames in the PDB file as a list of
        NATOM*3 lists.
        """
        if isinstance(filename, basestring):
            own_handle = True
            fileobj = TextToBinaryFile(genopen(filename, 'r'))
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
            cite = cont.getObj('citation_author')
            if cite is not None:
                nameidx = cite.getAttributeIndex('name')
                if nameidx != -1:
                    journal_authors = []
                    for i in xrange(cite.getRowCount()):
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
            chargeid = atoms.getAttributeIndex('pdbx_formal_charge')
            modelid = atoms.getAttributeIndex('pdbx_PDB_model_num')
            origmodel = None
            lastmodel = None
            xyz = []
            atommap = dict()
            last_atom = Atom()
            for i in xrange(atoms.getRowCount()):
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
                try:
                    model = int(row[modelid])
                except ValueError:
                    model = 0
                if origmodel is None:
                    origmodel = lastmodel = model
                x, y, z = float(row[xid]), float(row[yid]), float(row[zid])
                try:
                    occup = float(row[occupid])
                except ValueError:
                    occup = 0.0
                try:
                    bfactor = float(row[bfactorid])
                except ValueError:
                    bfactor = 0.0
                charge = row[chargeid]
                if not charge.strip() or charge.strip() in ('.', '?'):
                    charge = 0
                else:
                    try:
                        charge = float(charge)
                    except TypeError:
                        charge = 0.0
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
                            sym = '%s%s' % (sym[0].upper(), sym[0].lower())
                            atomic_number = AtomicNum[sym]
                            mass = Mass[sym]
                        except KeyError:
                            atomic_number = 0 # give up
                            mass = 0.0
                if atname.startswith('EP') or atname.startswith('LP'):
                    atom = ExtraPoint(atomic_number=atomic_number, name=atname,
                                charge=charge, mass=mass, occupancy=occup,
                                bfactor=bfactor, altloc=altloc, number=atnum)
                else:
                    atom = Atom(atomic_number=atomic_number, name=atname,
                                charge=charge, mass=mass, occupancy=occup,
                                bfactor=bfactor, altloc=altloc, number=atnum)
                atom.xx, atom.xy, atom.xz = x, y, z
                if _compare_atoms(last_atom, atom, resname, resnum, chain):
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
                        if lastmodel == origmodel:
                            struct.pdbxyz = [xyz]
                        else:
                            struct.pdbxyz.append(xyz)
                        xyz = []
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
                spaceid = cell.getAttributeIndex('space_group_name_H-M')
                row = cell.getRow(0)
                struct.box = create_array(
                        [float(row[aid]), float(row[bid]), float(row[cid]),
                         float(row[alphaid]), float(row[betaid]),
                         float(row[gammaid])]
                )
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
                                  'section. Skipping')
                else:
                    try:
                        for i in xrange(anisou.getRowCount()):
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
                            atommap[key].anisou = create_array(
                                    [u11, u22, u33, u12, u13, u23]
                            )
                    except (ValueError, KeyError):
                        # If at least one went wrong, set them all to None
                        for key, atom in atommap.iteritems():
                            atom.anisou = None
                        warnings.warn('Problem processing anisotropic '
                                      'B-factors. Skipping')
            if xyz:
                if len(xyz) != len(struct.atoms) * 3:
                    print len(xyz), len(struct.atoms)
                    raise ValueError('Corrupt CIF; all models must have the '
                                     'same atoms')
                try:
                    struct.pdbxyz.append(xyz)
                except AttributeError:
                    # Hasn't been assigned yet
                    struct.pdbxyz = xyz

        # Build the return value
        if len(structures) == 1:
            return structures[0]
        return tuple(structures)

    #===================================================

    @staticmethod
    def write(struct, dest, renumber=True, coordinates=None,
              altlocs='all', write_anisou=False):
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
            dest = TextToBinaryFile(genopen(dest, 'w'))
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
        if coordinates is not None:
            try:
                crdsize = len(coordinates)
            except TypeError:
                raise TypeError("Cannot find length of coordinates")
            if crdsize == len(struct.atoms):
                try:
                    coords = coordinates.flatten()
                except AttributeError:
                    try:
                        coords = list(itertools.chain(*coordinates))
                    except TypeError:
                        raise TypeError("Unsupported coordinate dimensionality")
                if len(coords) != len(struct.atoms) * 3:
                    raise TypeError("Unsupported coordinate shape")
            elif crdsize == len(struct.atoms) * 3:
                coords = coordinates
            else:
                raise TypeError("Coordinates has unexpected shape")
        else:
            coords = [[a.xx, a.xy, a.xz] for a in struct.atoms]
            coords = list(itertools.chain(*coords))
        # Create a function to process each atom and return which one we want
        # to print, based on our alternate location choice
        if altlocs == 'all':
            def print_atoms(atom, coords):
                i3 = atom.idx * 3
                return atom, atom.other_locations, coords[i3:i3+3]
        elif altlocs == 'first':
            def print_atoms(atom, coords):
                i3 = atom.idx * 3
                return atom, dict(), coords[i3:i3+3]
        elif altlocs == 'occupancy':
            def print_atoms(atom, coords):
                occ = atom.occupancy
                a = atom
                for key, item in atom.other_locations.iteritems():
                    if item.occupancy > occ:
                        occ = item.occupancy
                        a = item
                return a, dict(), [a.xx, a.xy, a.xz]
        else:
            raise Exception("Should not be here!")
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
        # Generator expression here instead of list comp, since the generator
        # need only execute until the first True (no point in wasting the time
        # or space creating the entire list just to see if one is True...)
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
        for res in struct.residues:
            if renumber:
                atoms = res.atoms
            else:
                atoms = sorted(res.atoms, key=lambda atom: atom.number)
            for atom in atoms:
                pa, others, (x, y, z) = print_atoms(atom, coords)
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
                        ['ATOM', anum, Element[pa.atomic_number].upper(),
                         pa.name, pa.altloc, res.name, res.chain, '?', rnum,
                         res.insertion_code, x, y, z, pa.occupancy, pa.bfactor,
                         '?', '?', '?', '?', '?', '', rnum, res.name, res.chain,
                         pa.name, '1']
                )
                if write_anisou and pa.anisou is not None:
                    cifanisou.append(
                            [anum, Element[pa.atomic_number].upper(), pa.name,
                             pa.altloc, res.name, res.chain, rnum, pa.anisou[0],
                             pa.anisou[1], pa.anisou[2], pa.anisou[3],
                             pa.anisou[4], pa.anisou[5], '?', '?', '?', '?',
                             '?', '?', rnum, res.name, res.chain, pa.name]
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
                    cifatoms.append(
                            ['ATOM', anum, Element[oatom.atomic_number].upper(),
                             oatom.name, oatom.altloc, res.name, res.chain, '?',
                             rnum, res.insertion_code, x, y, z, oatom.occupancy,
                             oatom.bfactor, '?', '?', '?', '?', '?', '',
                             rnum, res.name, res.chain, oatom.name, '1']
                    )
                    if write_anisou and oatom.anisou is not None:
                        cifanisou.append(
                                [anum, Element[oatom.atomic_number].upper(),
                                 oatom.name, oatom.altloc, res.name, res.chain,
                                 rnum, oatom.anisou[0], oatom.anisou[1],
                                 oatom.anisou[2], oatom.anisou[3],
                                 oatom.anisou[4], oatom.anisou[5], '?', '?',
                                 '?', '?', '?', '?', rnum, res.name, res.chain,
                                 oatom.name]
                        )
        # Now write the PDBx file
        writer = PdbxWriter(dest)
        writer.write([cont])
        if own_handle:
            dest.close()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Structure.write_pdb = PDBFile.write
Structure.write_cif = CIFFile.write
