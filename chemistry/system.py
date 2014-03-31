"""
This is a series of classes representing a particle or collection of particles
for general use in the chemistry package
"""

from chemistry.periodic_table import Element as _Element
from chemistry.periodic_table import AtomicNum as _AtomicNum

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
                if not arg.strip(): continue
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

    _ordered_args = ['name', 'number', 'insertion_code', 'chain', 'models']
    _types = dict(name=str, number=int, insertion_code=str,
                  chain=str, models=int)

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

    def __init__(self, *args, **kwargs):
        self.models = 0
        self.respermodel = None
        super(ChemicalSystem, self).__init__(*args, **kwargs)

    @classmethod
    def load_from_pdb(cls, pdb):
        return cls.load_from_open_pdb(open(pdb, 'r'))

    @classmethod
    def load_from_open_pdb(cls, pdb):
        inst = cls()
        for line in pdb:
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
                inst.add_atom(
                        Atom(atnum, atname, atomic_number, chg, x=x, y=y, z=z,
                             occupancy=occupancy, bfactor=bfactor,
                             altloc=altloc),
                        Residue(resname, resid, inscode, chain,
                                inst.models or 1),
                )
            if rec == 'MODEL ':
                if inst.respermodel is None and len(inst) > 0:
                    inst.respermodel = len(inst)
                inst.models += 1
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

_DUMMY_ATOM = Atom()
