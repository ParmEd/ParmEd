"""
This module contains the core base class for all of the chemical structures with
various topological and force field features.

Author: Jason Swails
Date: November 10, 2014
"""
try:
    import bz2
except ImportError:
    bz2 = None
from chemistry.exceptions import PDBError, PDBWarning
from chemistry.periodic_table import AtomicNum, Mass
from chemistry.topologyobjects import TrackedList, AtomList, ResidueList, Atom
try:
    import gzip
except ImportError:
    gzip = None
import re
import warnings

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

relatere = re.compile(r'RELATED ID: *(\w+) *RELATED DB: *(\w+)', re.I)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Structure(object):
    """
    A chemical structure composed of atoms, bonds, angles, torsions, and other
    topological features

    Attributes
    ----------
    atoms : AtomList
        List of all atoms in the structure
    residues : ResidueList
        List of all residues in the structure
    bonds : TrackedList(Bond)
        List of all bonds in the structure
    angles : TrackedList(Angle)
        List of all angles in the structure
    dihedrals : TrackedList(Dihedral)
        List of all dihedrals in the structure -- only one term per dihedral, so
        multi-term dihedral parameters will have the same 4 atoms appear
        multiple times in the list
    urey_bradleys : TrackedList(UreyBradley)
        List of all Urey-Bradley angle bends in the structure
    impropers : TrackedList(Improper)
        List of all CHARMM-style improper torsions in the structure
    cmaps : TrackedList(Cmap)
        List of all CMAP objects in the structure
    trigonal_angles : TrackedList(TrigonalAngle)
        List of all AMOEBA-style trigonal angles in the structure
    out_of_plane_bends : TrackedList(OutOfPlaneBends)
        List of all AMOEBA-style out-of-plane bending angles
    pi_torsions : TrackedList(PiTorsion)
        List of all AMOEBA-style pi-torsion angles
    stretch_bends : TrackedList(StretchBend)
        List of all AMOEBA-style stretch-bend compound bond/angle terms
    torsion_torsions : TrackedList(TorsionTorsion)
        List of all AMOEBA-style coupled torsion-torsion terms
    chiral_frames : TrackedList(ChiralFrame)
        List of all AMOEBA-style chiral frames defined in the structure
    multipole_frames : TrackedList(MultipoleFrame)
        List of all AMOEBA-style multipole frames defined in the structure
    box : list of 6 floats
        Box dimensions (a, b, c, alpha, beta, gamma) for the unit cell. If no
        box is defined, `box` is set to `None`

    This class also has a handful of type lists for each of the attributes above
    (excluding `atoms`, `residues`, `chiral_frames`, and `multipole_frames`).
    They are all TrackedList instances that are designed to hold the relevant
    parameter type. The list is:
        bond_types, angle_types, dihedral_types, urey_bradley_types,
        improper_types, cmap_types, trigonal_angle_types,
        out_of_plane_bend_types, pi_torsion_types, stretch_bend_types,
        torsion_torsion_types
    """

    def __init__(self):

        # Topological object lists
        self.atoms = AtomList()
        self.residues = ResidueList()
        self.bonds = TrackedList()
        self.angles = TrackedList()
        self.dihedrals = TrackedList()
        self.urey_bradleys = TrackedList()
        self.impropers = TrackedList()
        self.cmaps = TrackedList()
        self.trigonal_angles = TrackedList()
        self.out_of_plane_bends = TrackedList()
        self.pi_torsions = TrackedList()
        self.stretch_bends = TrackedList()
        self.torsion_torsions = TrackedList()
        self.chiral_frames = TrackedList()
        self.multipole_frames = TrackedList()

        # Parameter type lists
        self.bond_types = TrackedList()
        self.angle_types = TrackedList()
        self.dihedral_types = TrackedList()
        self.urey_bradley_types = TrackedList()
        self.improper_types = TrackedList()
        self.cmap_types = TrackedList()
        self.trigonal_angle_types = TrackedList()
        self.out_of_plane_bend_types = TrackedList()
        self.pi_torsion_types = TrackedList()
        self.stretch_bend_types = TrackedList()
        self.torsion_torsion_types = TrackedList()

        self.box = None

    def is_changed(self):
        """ Determines if any of the topology has changed for this structure """
        for attr in dir(self):
            if hasattr(getattr(self, attr), 'changed'):
                if getattr(getattr(self, attr), 'changed'):
                    return True
        return False

    def unchange(self):
        """ Toggles all lists so that they do not indicate any changes """
        for attr in dir(self):
            at = getattr(self, attr)
            if hasattr(at, 'changed'):
                setattr(at, 'changed', False)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def read_PDB(filename):
    """
    Read a PDB file and return a populated `Structure` class

    Parameters
    ----------
    filename : str or file-like
        Name of PDB file to read, or a file-like object that can iterate over
        the lines of a PDB. Compressed file names can be specified and are
        determined by file-name extension.

    Metadata
    --------
    The PDB parser also adds metadata to the returned Structure object that may
    be present in the PDB file

    experimental : str
        EXPDTA record
    journal : str
        JRNL record
    authors : str
        AUTHOR records
    keywords : str
        KEYWDS records
    doi : str
        DOI from the JRNL record
    pmid : str
        PMID from the JRNL record
    journal_authors : str
        Author info from the JRNL record
    volume_page : str
        Volume page from the JRNL record
    title : str
        TITL section of the JRNL record
    year : int=None
        Year that the article was published, from the JRNL record
    related_entries : list of (str, str)
        List of entries in other databases 

    Returns
    -------
    structure
    
    structure : Structure
        The Structure object initialized with all of the information from the
        PDB file.  No bonds or other topological features are added by default.

    Notes
    -----
    The returned structure has an extra attribute, pdbxyz, that contains all of
    the coordinates for all of the frames in the PDB file as a list of NATOM*3
    lists.
    """
    global relatere
    if isinstance(filename, basestring):
        own_handle = True
        if filename.endswith('.gz'):
            if gzip is None:
                raise ImportError('gzip is not available for compressed PDB')
            fileobj = gzip.open(filename, 'r')
        elif filename.endswith('.bz2'):
            if bz2 is None:
                raise ImportError('bz2 is not available for compressed PDB')
            fileobj = bz2.BZ2File(filename, 'r')
        else:
            fileobj = open(filename, 'r')
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
    last_resid = 1
    res_hex = False
    atom_hex = False

    try:
        for line in fileobj:
            try:
                line = line.decode('ascii')
            except AttributeError:
                # ssume this is a string in Py3 which doesn't have 'decode'
                pass
            rec = line[:6]
            if rec == 'ATOM  ' or rec == 'HETATM':
                atomno += 1
                atnum, atname, altloc = line[6:11], line[12:16], line[16]
                resname, chain, resid = line[17:20], line[21], line[22:26]
                inscode = line[26]
                x, y, z = line[30:38], line[38:46], line[47:54]
                occupancy, bfactor = line[54:60], line[60:66]
                elem, chg = line[76:78], line[78:80]
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
                if last_resid >= 9999:
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
                    else:
                        resid = int(resid)
                else:
                    resid = int(resid)
                # If the number has cycled, it too may be hexadecimal
                if atom_hex:
                    atnum = int(atnum, 16)
                else:
                    try:
                        atnum = int(atnum)
                    except ValueError:
                        atnum = int(atnum, 16)
                        atom_hex = True
                # It's possible that the residue number has cycled so much that
                # it is now filled with ****'s. In that case, start a new
                # residue if the current residue repeats the same atom name as
                # the 'last' residue. Do not worry about atom numbers going to
                # *****'s, since that is >1M atoms.
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
                    chg = 0
                if modelno == 1:
                    atom = Atom(atomic_number=atomic_number, name=atname,
                                charge=chg, mass=mass, occupancy=occupancy,
                                bfactor=bfactor, altloc=altloc)
                    atom.xx, atom.xy, atom.xz = float(x), float(y), float(z)
                    struct.residues.add_atom(atom, resname, resid,
                                             chain, inscode)
                    struct.atoms.append(atom)
                else:
                    try:
                        atom = struct.atoms[atomno-1]
                    except IndexError:
                        raise PDBError('Atom %d differs in MODEL %d [%s %s vs. '
                                       '%s %s]' % (atomno, modelno,
                                       atom.residue.name, atom.name, resname,
                                       atname))
                    if atom.residue.name != resname or atom.name != atname:
                        raise PDBError('Atom %d differs in MODEL %d [%s %s vs. '
                                       '%s %s]' % (atomno, modelno,
                                       atom.residue.name, atom.name, resname,
                                       atname))
                coordinates.extend([x, y, z])
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
            elif rec == 'MODEL ':
                if modelno == 1 and len(struct.atoms) == 0: continue
                if len(coordinates) > 0:
                    if len(struct.atoms)*3 != len(coordinates):
                        raise ValueError(
                                'Inconsistent atom numbers in some PDB models')
                    warnings.warn('MODEL not explicitly ended', PDBWarning)
                    all_coordinates.append(coordinates)
                    coordinates = []
                modelno += 1
                atomno = 0
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
                        struct.volume_page = line[50:61].strip()
                        try:
                            struct.year = int(line[62:66])
                        except ValueError:
                            pass
                elif part == 'PMID':
                    struct.pmid = line[19:].strip()
                elif part == 'DOI ':
                    struct.doi = line[19:].strip()
            elif rec == 'KEYWDS':
                struct.keywords += '%s ' % line[11:]
            elif rec == 'REMARK' and line[6:10] == ' 900':
                # Related entries
                rematch = relatere.match(line[11:])
                if rematch:
                    struct.related_entries.append(rematch.groups())
    finally:
        # Make sure our file is closed if we opened it
        if own_handle: fileobj.close()

    # Make the keywords into a list
    struct.keywords = [x.strip() for x in struct.keywords.split(',')
                                        if x.strip()]
    struct.unchange()
    if coordinates:
        if len(coordinates) != 3*len(struct.atoms):
            raise ValueError('bad number of atoms in some PDB models')
        all_coordinates.append(coordinates)
    struct.pdbxyz = all_coordinates
    return struct

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

