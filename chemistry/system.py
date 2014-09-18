"""
This is a series of classes representing a particle or collection of particles
for general use in the chemistry package

Author: Jason M. Swails
Contributors:
Date: May 27, 2014
"""

try:
    import bz2
except ImportError:
    bz2 = None
from chemistry.periodic_table import Element as _Element
from chemistry.periodic_table import AtomicNum as _AtomicNum
try:
    import gzip
except ImportError:
    gzip = None
import re

__all__ = ['Atom', 'Residue']

class _BaseParticle(object):
   
    _ordered_args = []
    _types = dict()
    _unordered_args = dict()

    def __init__(self, *args, **kwargs):
        """ Set all of the attributes from our arguments """
        # Go through our allowed ordered arguments
        for i, arg in enumerate(args):
            try:
                if not arg.strip(): continue
            except AttributeError:
                pass
            setattr(self, self._ordered_args[i],
                    self._types[self._ordered_args[i]](arg))
        for key in kwargs:
            # Skip over blank kwargs
            try:
                if not key.strip(): continue
            except AttributeError:
                pass
            if key in self._ordered_args:
                if not kwargs[key].strip():
                    setattr(self, key, self._types[key](kwargs[key]))
            elif key in self._unordered_args:
                setattr(self, key, self._unordered_args[key](kwargs[key]))
            else:
                setattr(self, key, kwargs[key])

class Atom(_BaseParticle):
    " An atom class that contains atomic properties "

    _ordered_args = ['number', 'name', 'atomic_number', 'charge']
    _types = dict(number=int, atomic_number=int, name=str, charge=str)
    _unordered_args = dict(x=float, y=float, z=float, occupancy=float,
                          bfactor=float, altloc=str)
   
    def __init__(self, *args, **kwargs):
        self.name = ''
        self.atomic_number = -1
        self.number = 0
        self.x, self.y, self.z = (0.0, 0.0, 0.0)
        self.occupancy = 0.0
        self.bfactor = 0.0
        self.charge = '0'  # as from PDB
        self.altloc = ''
        # Now go through the keyword arguments
        super(Atom, self).__init__(*args, **kwargs)

    def assign_to_residue(self, residue):
        """ Assigns this atom to a particular residue """
        self.residue = residue
   
    @property
    def element(self):
        if self.atomic_number == -1:
            return '????'
        return _Element[self.atomic_number]

    @property
    def coords(self):
        return (self.x, self.y, self.z)

    def __eq__(self, other):
        """ Two atoms are equivalent if... """
        if self.name != other.name: return False
        if self.atomic_number != other.atomic_number: return False
        if self.number != other.number: return False
        if self.occupancy != other.occupancy: return False
        if self.coords != other.coords: return False
        if self.altloc != other.altloc: return False
        if self.charge != other.charge: return False
        if hasattr(self, 'residue') and hasattr(other, 'residue'):
            if self.residue is not other.residue: return False
        return True

    def __repr__(self):
        return '<Atom %d [%s]; Elem %s; Bfac %.2f; Occ %.2f; Loc %s>' % (
                    self.number, self.name, self.element, self.bfactor,
                    self.occupancy, self.altloc
        )

    def __str__(self):
        return ('Atom %d [%s]:\n\tAtomic Number: %d\n\tTemp Factor:   %.2f\n\t'
                'Occupancy:     %.2f\n\tPosition:      [%.3f %.3f %.3f]\n\t'
                'Location:       %s' % (self.number, self.name,
                self.atomic_number, self.bfactor, self.occupancy, self.x,
                self.y, self.z, self.altloc)
      )

class Residue(_BaseParticle, list):
    """ A collection of atoms in a single 'residue' """

    _ordered_args = ['name', 'number', 'insertion_code', 'chain', 'model']
    _types = dict(name=str, number=int, insertion_code=str,
                  chain=str, model=int)

    def __init__(self, *args, **kwargs):
        self.name = ''
        self.number = ''
        self.insertion_code = ''
        self.chain = ''
        super(Residue, self).__init__(*args, **kwargs)

    def add_atom(self, atom):
        for atm in self:
            if atom == atm: return
        atom.assign_to_residue(self)
        super(Residue, self).append(atom)

    def __eq__(self, other):
        """
        Two residues are the same if their name, number, insertion code, and
        chain are the same
        """
        if self.name != other.name: return False
        if self.number != other.number: return False
        if self.insertion_code != other.insertion_code: return False
        if self.chain != other.chain: return False
        return True

    def __str__(self):
        return ('Residue:\n\tName:   %s\n\tNumber: %s\n\tChain:  %s\n\t'
                'InsCod: %s' % (self.name, self.number, self.chain,
                                self.insertion_code)
        )

    def __repr__(self):
        return ('<Residue %d (%s); Chain %s; InsCode %s>' %
                (self.number, self.name, self.chain, self.insertion_code)
        )
   
    def first_atoms(self):
        """
        Generator yielding each unique atom, and only the first location in each
        atom
        """
        last_atom = _DUMMY_ATOM
        for atom in self:
            # Require consective atoms to be sequential (as they always are in
            # the PDB)
            if atom.name == last_atom.name: continue
            last_atom = atom
            yield atom
      
class ChemicalSystem(list):
    """ A list of residues (this is a full system) """

    relatere = re.compile(r'RELATED ID: *(\w+) *RELATED DB: *(\w+)', re.I)

    def __init__(self, *args, **kwargs):
        self._models = 0
        self.respermodel = None
        super(ChemicalSystem, self).__init__(*args, **kwargs)
        self.box = None # [a, b, c, alpha, beta, gamma] or None

    @classmethod
    def load_from_pdb(cls, pdb):
        """
        Load a chemical system from a PDB file (filename provided). gzip or
        bzip2 compression are detected based on filename extension (.gz implies
        gzip and .bz2 implies bzip2). If either gzip or bz2 are unavailable,
        ImportError is raised.

        Parameters
        ----------
        pdb : str
            Name of the PDB file to parse

        Returns
        -------
        ChemicalSystem instance loaded from the PDB file
        """
        if pdb.endswith('.gz'):
            if gzip is None:
                raise ImportError('gzip not available for compressed PDB')
            return cls.load_from_open_pdb(gzip.open(pdb, 'r'))
        if pdb.endswith('.bz2'):
            if bz2 is None:
                raise ImportError('bzip is not available for compressed PDB')
            return cls.load_from_open_pdb(bz2.BZ2File(pdb, 'r'))
        return cls.load_from_open_pdb(open(pdb, 'r'))

    @classmethod
    def load_from_open_pdb(cls, pdb):
        inst = cls()
        last_resid = 1
        uses_hexadecimal = False
        atom_hexadecimal = False
        # Empty some information
        inst.experimental = inst.journal = inst.authors = inst.keywords = ''
        inst.doi = inst.pmid = inst.journal_authors = inst.volume_page = ''
        inst.title = ''
        inst.year = None
        # related_entries is a list of tuples of the form (ID, Database) where
        # both ID and Database are strings. Common values for Database are PDB,
        # BMRB, and EMDB
        inst.related_entries = []
        for line in pdb:
            try:
                line = line.decode('ascii') # In case we are parsing bytes (Py3)
            except AttributeError:
                # Assume this is a string in Py3 which doesn't have 'decode'
                pass
            rec = line[:6]
            if rec == 'ATOM  ' or rec == 'HETATM':
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
                    atomic_number = _AtomicNum[atsym]
                except KeyError:
                    atomic_number = -1
                # Figure out what my residue number is and see if the PDB is
                # outputting residue numbers in hexadecimal (VMD does this).
                if last_resid >= 9999:
                    if not uses_hexadecimal and resid == '9999':
                        resid = 9999
                    elif not uses_hexadecimal:
                        uses_hexadecimal = int(resid, 16) == 10000
                    # So now we know if we use hexadecimal or not. If we do,
                    # convert. Otherwise, stay put
                    if uses_hexadecimal:
                        try:
                            resid = int(resid, 16)
                        except ValueError, e:
                            if resid == '****':
                                resid = last_resid
                            else:
                                raise e
                    else:
                        resid = int(resid)
                else:
                    resid = int(resid)
                # If the atom number has cycled, it too may be hexadecimal
                if atom_hexadecimal:
                    atnum = int(atnum, 16)
                else:
                    try:
                        atnum = int(atnum)
                    except ValueError:
                        atnum = int(atnum, 16)
                        atom_hexadecimal = True
                last_resid = resid
                inst.add_atom(
                        Atom(atnum, atname, atomic_number, chg, x=x, y=y, z=z,
                             occupancy=occupancy, bfactor=bfactor,
                             altloc=altloc),
                        Residue(resname, resid, inscode, chain,
                                inst._models or 1),
                )
            elif rec == 'MODEL ':
                if inst.respermodel is None and len(inst) > 0:
                    inst.respermodel = len(inst)
                inst._models += 1
            elif rec == 'CRYST1':
                a = float(line[6:15])
                b = float(line[15:24])
                c = float(line[24:33])
                try:
                    A = float(line[33:40])
                    B = float(line[40:47])
                    C = float(line[47:54])
                except (IndexError, ValueError):
                    # Default to orthorhombic box
                    A = B = C = 90.0
                inst.box = [a, b, c, A, B, C]
            elif rec == 'EXPDTA':
                inst.experimental = line[6:].strip()
            elif rec == 'AUTHOR':
                inst.authors += line[10:].strip()
            elif rec == 'JRNL  ':
                part = line[12:16]
                if part == 'AUTH':
                    inst.journal_authors += line[19:].strip()
                elif part == 'TITL':
                    inst.title += ' %s' % line[19:].strip()
                elif part == 'REF ':
                    inst.journal += ' %s' % line[19:47].strip()
                    if not line[16:18].strip():
                        inst.volume_page = line[50:61].strip()
                        try:
                            inst.year = int(line[62:66])
                        except ValueError:
                            pass
                elif part == 'PMID':
                    inst.pmid = line[19:].strip()
                elif part == 'DOI ':
                    inst.doi = line[19:].strip()
            elif rec == 'KEYWDS':
                inst.keywords += '%s ' % line[11:].replace(',', ' ')
            elif rec == 'REMARK' and line[6:10] == ' 900':
                # Related entries
                rematch = ChemicalSystem.relatere.match(line[11:])
                if rematch:
                    inst.related_entries.append(rematch.groups())
        # Make the keywords a list
        inst.keywords = [x.strip() for x in inst.keywords.split(',')
                                        if x.strip()]
        # Strip off trailing and leading whitespace for some attributes
        inst.title = inst.title.strip()
        return inst
   
    @classmethod
    def load_from_pqr(cls, pqr):
        return cls.load_from_open_pqr(open(pqr, 'r'))

    @classmethod
    def load_from_open_pqr(cls, pqr):
        inst = cls()
        for line in pqr:
            # PQR files are whitespace-delimited
            words = line.split()
            if words[0] in ('ATOM', 'HETATM'):
                atnum, atname, resname, resid, x, y, z, chg = words[1:9]
                inst.add_atom(
                        Atom(atnum, atname, -1, chg, x=x, y=y, z=z),
                        Residue(resname, resid)
                )
        return inst

    def add_atom(self, atom, res):
        """ Adds an atom inside a residue """
        # Make sure our residue is new and not a continuation of the last one
        if len(self) != 0:
            if self[-1] == res:
                res = self[-1]
            else:
                self.append(res)
        else:
            self.append(res)
        res.add_atom(atom)

    @property
    def atoms(self):
        """ Generator that iterates through all of the atoms """
        for res in self:
            for atom in res:
                yield atom

    def positions(self, model=1):
        """
        Returns a list of the coordinates in the form [x1, y1, z1, x2, ...] for
        the requested model number. Model numbers start from 1
        """
        crds = []
        for res in self:
            if res.model != model: continue
            for atom in res:
                crds.extend([atom.x, atom.y, atom.z])
        if not crds:
            raise IndexError('Model %d out of range' % model)
        return crds

    @property
    def unique_atoms(self):
        """ Generator that iterates through only "unique" atoms """
        for res in self:
            for atom in res.first_atoms():
                yield atom

    def __len__(self):
        """
        Adjust the length to show number of residues -- correct for multiple
        models
        """
        if self.respermodel is None or self.respermodel <= 0:
            return super(ChemicalSystem, self).__len__()
        return self.respermodel

    @property
    def models(self):
        """
        Hack to get PDBs with only 1 structure to return the correct number of
        models (1)
        """
        return self._models or 1

_DUMMY_ATOM = Atom()
