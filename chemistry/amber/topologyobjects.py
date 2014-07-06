"""
This module contains objects that deal with system topology and parameters, such
as atoms, residues, bonds, angles, etc.

These are intended for use with the AmberParm and ChamberParm classes and are
used to recompute the parameter topology file when necessary (i.e., it provides
high-level access to Amber prmtop data)

Designed so it can be used like:
   from chemistry.amber.topologyobjects import *

by Jason Swails
"""
from __future__ import division

from chemistry.exceptions import (BondError, DihedralError, AmberParmError,
                                  CmapError)
from chemistry.amber.constants import NATOM, TINY
from chemistry.periodic_table import AtomicNum, Mass, Element as _Element
from compat24 import all, property

__all__ = ['Atom', 'Bond', 'BondType', 'Angle', 'AngleType', 'Dihedral',
           'DihedralType', 'Residue', 'ResidueList', 'AtomList', 'BondTypeList',
           'AngleTypeList', 'DihedralTypeList', 'TrackedList']

class Atom(object):
    """ 
    An atom. Only use these as elements in AtomList instances, since AtomList
    will keep track of when indexes and other stuff needs to be updated
    """
    #===================================================

    def __init__(self, parm, starting_index):
        self.atomic_number = 0
        self._bond_partners = set()
        self._angle_partners = set()
        self._dihedral_partners = set()
        self._exclusion_partners = set() # For arbitrary exclusions
        self.parm = parm
        self.idx = -1
        self.starting_index = starting_index
        self.load_from_parm()
        self.residue = None
        self.marked = 0 # For setting molecules
        self.bonds, self.angles, self.dihedrals = [], [], []
        # Chamber properties
        self.urey_bradleys, self.impropers, self.cmaps = [], [], []
        self.deleted = False
        self._has_loaded_exclusions = False
   
    #===================================================

    @property
    def bond_partners(self):
        """ Go through all bonded partners """
        return sorted(list(self._bond_partners))

    @property
    def angle_partners(self):
        """ List of all angle partners that are NOT bond partners """
        return sorted(list(self._angle_partners - self._bond_partners))

    @property
    def dihedral_partners(self):
        " List of all dihedral partners that are NOT angle or bond partners "
        return sorted(list(self._dihedral_partners - self._angle_partners -
                           self._bond_partners))

    @property
    def exclusion_partners(self):
        """
        List of all exclusions not otherwise excluded by bonds/angles/torsions
        """
        return sorted(list(self._exclusion_partners - self._dihedral_partners -
                           self._angle_partners - self._bond_partners))

    #===================================================

    def add_data(self):
        """ 
        Writes this atom's data to the AmberParm object. Don't pitch a fit if
        we're missing some useless (unused) flags.
        """
        # Determine how many excluded atoms we have. The only ones that count
        # are those with a smaller index (to avoid double-counting)
        numex = 0
        for atm in self.bond_partners:
            if atm.idx > self.idx: numex += 1
        for atm in self.angle_partners:
            if atm.idx > self.idx: numex += 1
        for atm in self.dihedral_partners:
            if atm.idx > self.idx: numex += 1
        for atm in self.exclusion_partners:
            if atm.idx > self.idx: numex += 1
        # For some reason, existing topology files follow the convention that
        # atoms with no exclusions (because all bonded partners have atom #s
        # lower than theirs) have num_excluded = 1, with a 0 placeholder in
        # EXCLUDED_ATOMS_LIST... Weird.
        if numex == 0: numex = 1
        # Make sure we're indexing from 0
        parm_data = self.parm.parm_data
        parm_data['ATOM_NAME'][self.idx] = self.atname[:4]
        parm_data['CHARGE'][self.idx] = self.charge
        parm_data['MASS'][self.idx] = self.mass
        parm_data['ATOM_TYPE_INDEX'][self.idx] = self.nb_idx
        parm_data['NUMBER_EXCLUDED_ATOMS'][self.idx] = numex
        parm_data['AMBER_ATOM_TYPE'][self.idx] = self.attype[:4]
        parm_data['JOIN_ARRAY'][self.idx] = 0
        parm_data['TREE_CHAIN_CLASSIFICATION'][self.idx] = self.tree[:4]
        parm_data['IROTAT'][self.idx] = 0
        parm_data['RADII'][self.idx] = self.radii
        parm_data['SCREEN'][self.idx] = self.screen
        try:
            parm_data['ATOMIC_NUMBER'][self.idx] = self.atomic_number
        except KeyError:
            pass

    #===================================================

    # Make 'element' an alias for 'atomic_number'

    @property
    def element(self):
        return self.atomic_number
    @element.setter
    def element(self, value):
        self.atomic_number = value

    #===================================================

    def load_from_parm(self):
        """ Load data from the AmberParm class """
        parm_data = self.parm.parm_data
        si = self.starting_index
        self.atname = parm_data['ATOM_NAME'][si]
        self.charge = parm_data['CHARGE'][si]
        self.mass = parm_data['MASS'][si]
        self.nb_idx = parm_data['ATOM_TYPE_INDEX'][si]
        self.attype = parm_data['AMBER_ATOM_TYPE'][si]
        self.tree = parm_data['TREE_CHAIN_CLASSIFICATION'][si]
        # Put in some 
        if 'RADII' in parm_data:
            self.radii = parm_data['RADII'][si]
        else:
            self.radii = 0.0 # dummy number
        if 'SCREEN' in parm_data:
            self.screen = parm_data['SCREEN'][si]
        else:
            self.screen = 0.0 # dummy number
        if 'ATOMIC_NUMBER' in parm_data:
            self.atomic_number = parm_data['ATOMIC_NUMBER'][si]
        else:
            # Determine from mass (unless already done)
            if self.atomic_number <= 0:
                self.atomic_number = AtomicNum[Element(self.mass)]

        # Load the positions and velocities if the amberParm object them
        if hasattr(self.parm, 'coords'):
            self.xx = self.parm.coords[si*3  ]
            self.xy = self.parm.coords[si*3+1]
            self.xz = self.parm.coords[si*3+2]
            if self.parm.hasvels:
                self.vx = self.parm.vels[si*3  ]
                self.vy = self.parm.vels[si*3+1]
                self.vz = self.parm.vels[si*3+2]

    #===================================================
      
    def load_exclusions(self):
        """
        Looks at the NUMBER_EXCLUDED_ATOMS and EXCLUDED_ATOMS_LIST to determine
        if the topology file defines more atom exclusions than what is defined
        simply by the bonds, angles, and dihedrals. It then adds these atoms to
        the exclusion_partners array. Only allow this to occur once, though,
        since exclusions are remembered for the life of the object
        """
        if self._has_loaded_exclusions:
            return
        excset = set()
        exclat = self.parm.parm_data['NUMBER_EXCLUDED_ATOMS']
        exclist = self.parm.parm_data['EXCLUDED_ATOMS_LIST']
        first_excl = sum(exclat[:self.starting_index])
        nexcl = exclat[self.starting_index]
        for i in range(nexcl):
            idx = exclist[first_excl+i] - 1
            # Skip over placeholders
            if idx < 0: continue
            excset.add(self.parm.atom_list[idx])
        # Now subtract off all of the bonds, angles, and dihedrals
        excset = (excset - self._bond_partners - self._angle_partners -
                  self._dihedral_partners)
        for atm in excset:
            self.exclude(atm)
        self._has_loaded_exclusions = True

    #===================================================

    def bond_to(self, other):
        """ Log this atom as bonded to another atom.  """
        if self is other:
            raise BondError("Cannot bond atom to itself!")
        self._bond_partners.add(other)

    #===================================================
      
    def angle_to(self, other):
        """ Log this atom as angled to another atom.  """
        if self is other:
            raise BondError("Cannot angle an atom with itself!")
        self._angle_partners.add(other)
   
    #===================================================

    def dihedral_to(self, other):
        """ Log this atom as dihedral-ed to another atom.  """
        if self is other:
            raise BondError("Cannot dihedral an atom with itself!")
        self._dihedral_partners.add(other)
      
    #===================================================

    def exclude(self, other):
        """ Add one atom to my arbitrary exclusion list """
        if self is other:
            raise BondError("Cannot exclude an atom from itself")
        self._exclusion_partners.add(other)
        # If he is excluded from me, then I am excluded from him
        other._exclusion_partners.add(self)

    #===================================================

    def reset_topology(self):
        """
        Deletes all of the bond, angle, and dihedral partners so they can be set
        up again with updated data. Keep the arbitrary exclusions, though
        """
        self._bond_partners = set()
        self._angle_partners = set()
        self._dihedral_partners = set()

    #===================================================

    # Comparisons are done by comparing the starting indexes

    def __eq__(self, other):
        return self.starting_index == other.starting_index
      
    def __ne__(self, other):
        return not Atom.__eq__(self, other)

    def __gt__(self, other):
        return self.starting_index > other.starting_index

    def __lt__(self, other):
        return self.starting_index < other.starting_index

    def __ge__(self, other):
        return not Atom.__lt__(self, other)

    def __le__(self, other):
        return not Atom.__gt__(self, other)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Bond(object):
    """ Bond class. Stores 2 atoms involved and force constant/equil value """

    #===================================================

    def __init__(self, atom1, atom2, bond_type):
        """ Bond constructor """
        # Make sure we're not bonding me to myself
        if atom1 is atom2:
            raise BondError('Cannot bond atom to itself!')
        # Order the atoms so the lowest atom # is first
        self.atom1 = atom1
        self.atom2 = atom2
        # Load this struct into the atoms
        self.atom1.bonds.append(self)
        self.atom2.bonds.append(self)
        # Load the force constant and equilibrium distance
        self.bond_type = bond_type
        self.register()

    #===================================================

    def write_info(self, parm, key, idx):
        """ Writes the bond info to the topology file. idx starts at 0 """
        parm.parm_data[key][3*idx  ] = 3*(self.atom1.idx)
        parm.parm_data[key][3*idx+1] = 3*(self.atom2.idx)
        parm.parm_data[key][3*idx+2] = self.bond_type.idx + 1

    #===================================================
   
    def register(self):
        """ Register each atom as bonded to the other """
        self.atom1.bond_to(self.atom2)
        self.atom2.bond_to(self.atom1)
   
    #===================================================
   
    def __contains__(self, thing):
        """ Quick and easy way to see if an Atom is in this Bond """
        return thing is self.atom1 or thing is self.atom2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondType(object):
    """ A bond type """

    #===================================================

    def __init__(self, k, req, idx):
        """BondType constructor. idx must start from 0!!! """
        self.idx = idx
        self.k = k
        self.req = req

    #===================================================

    def write_info(self, parm):
        """ Writes the bond parameters in the parameter file """
        # If our index is -1 (we're not being used), just return
        if self.idx == -1: return
        parm.parm_data['BOND_FORCE_CONSTANT'][self.idx] = self.k
        parm.parm_data['BOND_EQUIL_VALUE'][self.idx] = self.req

    #===================================================

    def __eq__(self, other):
        return self.k == other.k and self.req == other.req

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Angle(object):
    """ Angle class. Stores 3 atoms involved and force constant/equil value """
      
    #===================================================

    def __init__(self, atom1, atom2, atom3, angle_type):
        """ Angle constructor """
        # Make sure we're not angling me to myself
        if atom1 is atom2 or atom1 is atom3 or atom2 is atom3:
            raise BondError('Cannot angle atom to itself!')
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        # Log these angles in each atom
        atom1.angles.append(self)
        atom2.angles.append(self)
        atom3.angles.append(self)
        # Load the force constant and equilibrium angle
        self.angle_type = angle_type
        self.register()

    #===================================================

    def write_info(self, parm, key, idx):
        """ Write the info to the topology file """
        parm.parm_data[key][4*idx  ] = 3*(self.atom1.idx)
        parm.parm_data[key][4*idx+1] = 3*(self.atom2.idx)
        parm.parm_data[key][4*idx+2] = 3*(self.atom3.idx)
        parm.parm_data[key][4*idx+3] = self.angle_type.idx + 1

    #===================================================

    def register(self):
        """ Register each atom as angled to each other atom """
        self.atom1.angle_to(self.atom2)
        self.atom1.angle_to(self.atom3)
        self.atom2.angle_to(self.atom1)
        self.atom2.angle_to(self.atom3)
        self.atom3.angle_to(self.atom1)
        self.atom3.angle_to(self.atom2)
   
    #===================================================

    def __contains__(self, thing):
        """ Quick and easy way to see if an Atom or a Bond is in this Angle """
        if isinstance(thing, Bond):
            return ((self.atom1 in thing and self.atom2 in thing) or 
                    (self.atom2 in thing and self.atom3 in thing))
        # Assume it's an atom
        return thing is self.atom1 or thing is self.atom2 or thing is self.atom3

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleType(object):
    """ An angle type """
    #===================================================

    def __init__(self, k, theteq, idx):
        """ AngleType constructor. idx must start from 0!!! """
        self.k = k
        self.theteq = theteq
        self.idx = idx

    #===================================================

    def write_info(self, parm):
        """ Writes the bond parameters in the parameter file """
        # If we're not being used (idx == -1) just return
        if self.idx == -1: return
        parm.parm_data['ANGLE_FORCE_CONSTANT'][self.idx] = self.k
        parm.parm_data['ANGLE_EQUIL_VALUE'][self.idx] = self.theteq

    #===================================================

    def __eq__(self, other):
        return self.k == other.k and self.theteq == other.theteq

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Dihedral(object):
    """
    Dihedral class with 4 atoms involved and force constant/periodicity/phase
    """
      
    #===================================================

    def __init__(self, atom1, atom2, atom3, atom4, dihed_type, signs):
        """ Dihedral constructor. idx must start from 0!!! """
        # Make sure we're not dihedraling me to myself
        atmlist = [atom1, atom2, atom3, atom4]
        for i in range(len(atmlist)):
            for j in range(i+1, len(atmlist)):
                if atmlist[i] is atmlist[j]:
                    raise BondError('Cannot dihedral atom to itself!')
        # Set up instances
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        # Log these dihedrals in each atom
        atom1.dihedrals.append(self)
        atom2.dihedrals.append(self)
        atom3.dihedrals.append(self)
        atom4.dihedrals.append(self)
        # Load the force constant and equilibrium angle
        self.dihed_type = dihed_type
        self.signs = signs # is our 3rd or 4th term negative?
        self.register()

    #===================================================

    def write_info(self, parm, key, idx):
        """ Write the info to the topology file """
        # In order to preserve multi-term and improper dihedral definitions, we
        # have to make sure that atom index 0 is _never_ in either of the last
        # 2 positions (since you can't have '-0'. Therefore, if atoms 3 or 4 are
        # index 0 and their sign is negative, swap the dihedral
        if (self.atom3.idx == 0 and self.signs[0] == -1) or \
            (self.atom4.idx == 0 and self.signs[1] == -1):
            parm.parm_data[key][5*idx  ] = 3*(self.atom4.idx)
            parm.parm_data[key][5*idx+1] = 3*(self.atom3.idx)
            parm.parm_data[key][5*idx+2] = 3*(self.atom2.idx) * self.signs[0]
            parm.parm_data[key][5*idx+3] = 3*(self.atom1.idx) * self.signs[1]
        else:
            parm.parm_data[key][5*idx  ] = 3*(self.atom1.idx)
            parm.parm_data[key][5*idx+1] = 3*(self.atom2.idx)
            parm.parm_data[key][5*idx+2] = 3*(self.atom3.idx) * self.signs[0]
            parm.parm_data[key][5*idx+3] = 3*(self.atom4.idx) * self.signs[1]
        parm.parm_data[key][5*idx+4] = self.dihed_type.idx + 1

    #===================================================

    def register(self):
        """ Register each atom as dihedral-ed to each other atom """
        self.atom1.dihedral_to(self.atom2)
        self.atom1.dihedral_to(self.atom3)
        self.atom1.dihedral_to(self.atom4)
        self.atom2.dihedral_to(self.atom1)
        self.atom2.dihedral_to(self.atom3)
        self.atom2.dihedral_to(self.atom4)
        self.atom3.dihedral_to(self.atom2)
        self.atom3.dihedral_to(self.atom1)
        self.atom3.dihedral_to(self.atom4)
        self.atom4.dihedral_to(self.atom2)
        self.atom4.dihedral_to(self.atom3)
        self.atom4.dihedral_to(self.atom1)

    #===================================================

    def __contains__(self, thing):
        """
        Quick and easy way to find out if an Atom or Bond is in this Dihedral
        """
        if isinstance(thing, Bond):
            # A dihedral is made up of 3 bonds
            if self.signs[1] == -1:
                # An improper is different... Atom 3 is the central atom
                return ((self.atom1 in thing and self.atom3 in thing) or
                        (self.atom2 in thing and self.atom3 in thing) or
                        (self.atom4 in thing and self.atom3 in thing))
            return ((self.atom1 in thing and self.atom2 in thing) or
                    (self.atom2 in thing and self.atom3 in thing) or
                    (self.atom3 in thing and self.atom4 in thing))
        # Otherwise assume thing is an Atom
        return (thing is self.atom1 or thing is self.atom2 or
                thing is self.atom3 or thing is self.atom4)

    #===================================================

    def __eq__(self, thing):
        """
        A dihedral is equivalent if the 4 atoms are the same (or reverse) in
        order Allow comparison with another type of dihedral or with a list of
        4 integer indexes (or tuple instead of a list)
        """
        if isinstance(thing, Dihedral):
            # I'm comparing with another Dihedral here
            return ((self.atom1 is thing.atom1 and self.atom2 is thing.atom2 and
                     self.atom3 is thing.atom3 and self.atom4 is thing.atom4) or
                    (self.atom1 is thing.atom4 and self.atom2 is thing.atom3 and
                     self.atom4 is thing.atom1)
            )
        if isinstance(thing, list) or isinstance(thing, tuple):
            # Here, atoms are expected to index from 0 (Python standard) if we
            # are comparing with a list or tuple
            if len(thing) != 4:
                raise DihedralError('comparative %s has %d elements! Expect 4.'
                                    % (type(thing).__name__, len(thing)))
            # Compare starting_index, since we may not have an index right now
            return ( (self.atom1.starting_index == thing[0] and 
                    self.atom2.starting_index == thing[1] and
                    self.atom3.starting_index == thing[2] and
                    self.atom4.starting_index == thing[3]) or
                    (self.atom1.starting_index == thing[3] and 
                    self.atom2.starting_index == thing[2] and
                    self.atom3.starting_index == thing[1] and
                    self.atom4.starting_index == thing[0]) )

        raise TypeError('Cannot compare Dihedral with %s' %
                        type(thing).__name__)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralType(object):
    """ A type of dihedral """

    #===================================================
   
    def __init__(self, phi_k, per, phase, scee, scnb, idx):
        """ DihedralType constructor """
        self.phi_k = phi_k
        self.per = per
        self.phase = phase
        self.scee = scee
        self.scnb = scnb
        self.idx = idx

   #===================================================

    def write_info(self, parm):
        """ Write out the dihedral parameters """
        # If our idx == -1, we're not being used, so just return
        if self.idx == -1: return
        parm.parm_data['DIHEDRAL_FORCE_CONSTANT'][self.idx] = self.phi_k
        parm.parm_data['DIHEDRAL_PERIODICITY'][self.idx] = self.per
        parm.parm_data['DIHEDRAL_PHASE'][self.idx] = self.phase
        try:
            parm.parm_data['SCEE_SCALE_FACTOR'][self.idx] = self.scee
        except KeyError:
            pass
        try:
            parm.parm_data['SCNB_SCALE_FACTOR'][self.idx] = self.scnb
        except KeyError:
            pass
   
    #===================================================

    def __eq__(self, other):
        return (self.phi_k == other.phi_k and self.per == other.per and 
                self.phase == other.phase and self.scee == other.scee and
                self.scnb == other.scnb)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradley(object):
    " Urey-Bradley class. Stores 2 atoms involved, force constant/equil value "

   #===================================================

    def __init__(self, atom1, atom2, ub_type):
        """ Bond constructor """
        # Make sure we're not bonding me to myself
        if atom1 is atom2:
            raise BondError('Cannot angle atom to itself!')
        # Order the atoms so the lowest atom # is first
        self.atom1 = atom1
        self.atom2 = atom2
        # Log this urey-bradley in the atoms
        atom1.urey_bradleys.append(self)
        atom2.urey_bradleys.append(self)
        # Load the force constant and equilibrium distance
        self.ub_type = ub_type

    #===================================================

    def write_info(self, parm, key, idx):
        """ Writes the bond info to the topology file. idx starts at 0 """
        parm.parm_data[key][3*idx  ] = self.atom1.idx + 1
        parm.parm_data[key][3*idx+1] = self.atom2.idx + 1
        parm.parm_data[key][3*idx+2] = self.ub_type.idx + 1

    #===================================================
   
    def __contains__(self, thing):
        " Quick and easy way to see if an Atom or Bond is in this Urey-Bradley "
        # If this is a bond things are a bit more complicated since we don't
        # know the central atom of this Urey-Bradley. We need to make sure that
        # one of the atoms of the bond is either atom1 or atom2 and that the
        # OTHER atom in Bond "thing" has the OTHER atom in this Urey-Bradley in
        # its list of bonded partners.
        if isinstance(thing, Bond):
            if not thing.atom1 in self:
                if not thing.atom2 in self:
                    # Neither atom is in this Urey-Bradley...
                    return False
                # If we are here, thing.atom2 is in self
                end1 = thing.atom2
                cent = thing.atom1
            else:
                # If we are here, thing.atom1 is in self
                end1 = thing.atom1
                cent = thing.atom2
            # If we are here, end1 and cent are set. Look through the bond
            # partners of cent(er) and see if any of them is in this
            # Urey-Bradley (but ONLY if that atom is not the original end1)
            for atm in cent.bond_partners:
                if atm is end1: continue
                if atm in self:
                # If we got here, we found both atoms in this Urey-Bradley
                # separated by 2 bonds, so "thing" IS in this Urey-Bradley
                    return True
            # If we got here, we could not find the other atom connected by 2
            # bonds
            return False
        # If we are here, "thing" must not be a Bond and we therefore assume
        # it's an Atom
        return thing is self.atom1 or thing is self.atom2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradleyType(object):
    """ A Urey-Bradley type """

    #===================================================

    def __init__(self, k, req, idx):
        """UreyBradleyType constructor. idx must start from 0!!! """
        self.idx = idx
        self.k = k
        self.req = req

    #===================================================

    def write_info(self, parm):
        """ Writes the bond parameters in the parameter file """
        # If our index is -1 (we're not being used), just return
        if self.idx == -1: return
        parm.parm_data['CHARMM_UREY_BRADLEY_FORCE_CONSTANT'][self.idx] = self.k
        parm.parm_data['CHARMM_UREY_BRADLEY_EQUIL_VALUE'][self.idx] = self.req

    #===================================================

    def __eq__(self, other):
        return self.k == other.k and self.req == other.req

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Improper(object):
    """ Improper class with 4 atoms involved and force constant/phase """
      
    #===================================================

    def __init__(self, atom1, atom2, atom3, atom4, improp_type):
        """ Dihedral constructor. idx must start from 0!!! """
        # Make sure we're not dihedraling me to myself
        atmlist = [atom1, atom2, atom3, atom4]
        for i in range(len(atmlist)):
            for j in range(i+1, len(atmlist)):
                if atmlist[i] is atmlist[j]:
                    raise BondError('Cannot improper atom to itself!')
        # Set up instances
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        # Log these impropers in each atom
        atom1.impropers.append(self)
        atom2.impropers.append(self)
        atom3.impropers.append(self)
        atom4.impropers.append(self)
        # Load the force constant and equilibrium angle
        self.improp_type = improp_type

    #===================================================

    def write_info(self, parm, key, idx):
        """ Write the info to the topology file """
        parm.parm_data[key][5*idx  ] = self.atom1.idx + 1
        parm.parm_data[key][5*idx+1] = self.atom2.idx + 1
        parm.parm_data[key][5*idx+2] = self.atom3.idx + 1
        parm.parm_data[key][5*idx+3] = self.atom4.idx + 1
        parm.parm_data[key][5*idx+4] = self.improp_type.idx + 1

   #===================================================

    def __contains__(self, thing):
        """
        Quick and easy way to find out if an Atom or Bond is in this Improper
        """
        # An improper is defined as shown below
        #             A3
        #             |
        #             |
        #    A4 ----- A1 ----- A2
        #
        # So the bonds will either be between atom1 and any other atom
        if isinstance(thing, Bond):
            return ((self.atom1 in thing and self.atom2 in thing) or
                    (self.atom1 in thing and self.atom3 in thing) or
                    (self.atom1 in thing and self.atom4 in thing))
        # Here we assume that we just have an atom
        return (thing is self.atom1 or thing is self.atom2 or
                thing is self.atom3 or thing is self.atom4)

    #===================================================

    def __eq__(self, thing):
        """
        A dihedral is equivalent if the 4 atoms are the same (or reverse) in
        order Allow comparison with another type of dihedral or with a list of 4
        atoms (or tuple)
        """
        if isinstance(thing, Improper):
            # I'm comparing with another Improper here. Central atom must be
            # the same. Others can be in any order
            if self.atom3 != thing.atom3:
                return False
            # Make a set with the remaining atoms. If they are equal, the
            # impropers are equivalent
            selfset = set([self.atom1, self.atom2, self.atom4])
            otherset = set([thing.atom1, thing.atom2, thing.atom3])
            return selfset == otherset
        if isinstance(thing, list) or isinstance(thing, tuple):
            # Here, atoms are expected to index from 0 (Python standard) if we
            # are comparing with a list or tuple
            if len(thing) != 4:
                raise DihedralError('comparative %s has %d elements! Expect 4.'
                                    % (type(thing).__name__, len(thing)))
            if self.atom3.starting_index != thing[2]:
                return False
            selfset = set([self.atom1.starting_index, self.atom2.starting_index,
                            self.atom4.starting_index])
            otherset = set([thing[0], thing[1], thing[3]])
            return selfset == otherset

        raise TypeError('Cannot compare Improper with %s' %
                        type(thing).__name__)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ImproperType(object):
    """ A type of improper torsion """

    #===================================================
   
    def __init__(self, psi_k, psi_eq, idx):
        """ DihedralType constructor """
        self.psi_k = psi_k
        self.psi_eq = psi_eq
        self.idx = idx

    #===================================================

    def write_info(self, parm):
        """ Write out the dihedral parameters """
        # If our idx == -1, we're not being used, so just return
        if self.idx == -1: return
        parm.parm_data['CHARMM_IMPROPER_FORCE_CONSTANT'][self.idx] = self.psi_k
        parm.parm_data['CHARMM_IMPROPER_PHASE'][self.idx] = self.psi_eq
   
    #===================================================

    def __eq__(self, other):
        return self.psi_k == other.psi_k and self.psi_eq == other.psi_eq

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Cmap(object):
    """ A coupled-torsion correction map term """

    #===================================================

    def __init__(self, atom1, atom2, atom3, atom4, atom5, cmap_type):
        """ Takes a coupled-torsion (5 atoms) """
        # Make sure we're not CMAPping me to myself
        atmlist = [atom1, atom2, atom3, atom4, atom5]
        for i in range(len(atmlist)):
            for j in range(i+1, len(atmlist)):
                if atmlist[i] is atmlist[j]:
                    raise BondError('Cannot cmap atom to itself!')
        # Set up instances
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.atom5 = atom5
        # Log these cmaps in each atom
        atom1.cmaps.append(self)
        atom2.cmaps.append(self)
        atom3.cmaps.append(self)
        atom4.cmaps.append(self)
        atom4.cmaps.append(self)
        # Load the CMAP interpolation table
        self.cmap_type = cmap_type

    #===================================================

    def write_info(self, parm, key, idx):
        """ Write the info to the topology file """
        parm.parm_data[key][6*idx  ] = self.atom1.idx + 1
        parm.parm_data[key][6*idx+1] = self.atom2.idx + 1
        parm.parm_data[key][6*idx+2] = self.atom3.idx + 1
        parm.parm_data[key][6*idx+3] = self.atom4.idx + 1
        parm.parm_data[key][6*idx+4] = self.atom5.idx + 1
        parm.parm_data[key][6*idx+5] = self.cmap_type.idx + 1

    #===================================================

    def __contains__(self, thing):
        """
        Quick and easy way to find out an atom or bond is in this
        coupled-torsion
        """
        # If we have a bond, see if it is 1 of the 4 sequential bonds that make
        # up this coupled-torsion
        if isinstance(thing, Bond):
            return ((self.atom1 in thing and self.atom2 in thing) or
                    (self.atom2 in thing and self.atom3 in thing) or
                    (self.atom3 in thing and self.atom4 in thing) or
                    (self.atom4 in thing and self.atom5 in thing))
        # Otherwise assume we have an atom
        return (thing is self.atom1 or thing is self.atom2 or
                thing is self.atom3 or thing is self.atom4 or
                thing is self.atom5)

    #===================================================

    def __eq__(self, thing):
        """
        A coupled-torsion is equivalent if the 5 atoms are in the same or
        reverse order Allow comparison with another type of cmap or with a
        sequence of 5 indexes
        """
        if isinstance(thing, Cmap):
            # I'm comparing with another Improper here. Central atom must be the
            # same. Others can be in any order
            return ((self.atom1 is thing.atom1 and self.atom2 is thing.atom2 and
                    self.atom3 is thing.atom3 and self.atom4 is thing.atom4 and
                    self.atom5 is thing.atom5) or
                    (self.atom1 is thing.atom5 and self.atom2 is thing.atom4 and
                    self.atom3 is thing.atom3 and self.atom4 is thing.atom2 and
                    self.atom5 is thing.atom1))
        if isinstance(thing, list) or isinstance(thing, tuple):
            # Here, atoms are expected to index from 0 (Python standard) if we
            # are comparing with a list or tuple
            if len(thing) != 5:
                raise DihedralError('comparative %s has %d elements! Expect 4.'
                                    % (type(thing).__name__, len(thing)))
            return ((self.atom1.starting_index == thing[0] and
                     self.atom2.starting_index == thing[1] and
                     self.atom3.starting_index == thing[2] and
                     self.atom4.starting_index == thing[3] and
                     self.atom5.starting_index == thing[4]) or
                    (self.atom1.starting_index == thing[4] and
                     self.atom2.starting_index == thing[3] and
                     self.atom3.starting_index == thing[2] and
                     self.atom4.starting_index == thing[1] and
                     self.atom5.starting_index == thing[0]))
        raise TypeError('Cannot compare Improper with %s' %
                        type(thing).__name__)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CmapType(object):
    """ Contains a correction map interpolation grid """

    #===================================================

    def __init__(self, resolution, grid, cmts, idx):
        self.resolution = resolution
        self.grid = _CmapGrid(resolution, grid)
        self.comments = cmts # To avoid losing comments in the chamber prmtop
        if len(grid) != self.resolution * self.resolution:
            raise CmapError('CMAP grid does not match expected resolution')
        self.idx = idx

    #===================================================

    def write_info(self, parm):
        """ Write out the dihedral parameters """
        # If our idx == -1, we're not being used, so just return
        i = self.idx
        if i == -1: return
        parm.parm_data['CHARMM_CMAP_RESOLUTION'][i] = self.resolution
        parm.parm_data['CHARMM_CMAP_PARAMETER_%02d' % (i+1)] = self.grid[:]

    #===================================================

    def __eq__(self, other):
        return (self.resolution == other.resolution and
                all([abs(i - j) < TINY for i, j in zip(self.grid, other.grid)]))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _CmapGrid(object):
    """
    A grid object for storing Correction map data. Data can be accessed in one
    of two ways; either with 1 or 2 indexes. If 2 indexes are given, the index
    into the flattened array is i*resolution+j. Indexing starts from 0.

    The _CmapGrid usually has ranges for the two angles from -180 to 180. Some
    places will expect the range to be 0-360 degrees (e.g., OpenMM). The
    switch_range method returns a _CmapGrid with this range. This method will
    also work to go backwards.

    Example:
    >>> g = _CmapGrid(2, [0, 1, 2, 3])
    >>> print g[0], g[0,0]
    0 0
    >>> print g[1], g[0,1]
    1 1
    >>> print g[1,0]
    2
    >>> g[1,1] = 10
    >>> print g[3]
    10
    >>> print g.switch_range()
    [10.0000, 2.0000
     1.0000, 0.0000]
    """

    def __init__(self, resolution, data=None):
        self.resolution = resolution
        if data is None:
            self._data = [0 for i in range(self.resolution*self.resolution)]
        else:
            self._data = data

    @property
    def transpose(self):
        """ Returns the transpose of the grid """
        try:
            return self._transpose
        except AttributeError:
            pass
        _transpose = []
        size = len(self._data)
        for i in range(self.resolution):
            piece = [self[j] for j in range(i, size, self.resolution)]
            _transpose += piece
        self._transpose = _CmapGrid(self.resolution, _transpose)
        return self._transpose

    T = transpose

    def __getitem__(self, idx):
        if isinstance(idx, tuple):
            return self._data[self.resolution*idx[0]+idx[1]]
        return self._data[idx]

    def __setitem__(self, idx, val):
        if isinstance(idx, tuple):
            self._data[self.resolution*idx[0]+idx[1]] = val
        else:
            self._data[idx] = val

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __repr__(self):
        return '<_CmapGrid: %dx%d>' % (self.resolution, self.resolution)

    def __str__(self):
        retstr = '[%.4f,' % self._data[0]
        fmt = ' %.4f'
        for i, val in enumerate(self):
            if i == 0: continue
            retstr += fmt % val
            if (i+1) % self.resolution == 0 and i != len(self._data) - 1:
                retstr += '\n'
            elif i != len(self) - 1:
                retstr += ','
        return retstr + ']'

    def __eq__(self, other):
        try:
            if self.resolution != other.resolution:
                return False
            for x, y in zip(self, other):
                if abs(x - y) > TINY:
                    return False
            return True
        except AttributeError:
            return TypeError('Bad type comparison with _CmapGrid')

    def switch_range(self):
        """
        Returns a grid object whose range is 0 to 360 degrees in both dimensions
        instead of -180 to 180 degrees (or -180 to 180 degrees if the range is
        already 0 to 360 degrees)
        """
        res = self.resolution
        mid = res // 2
        newgrid = _CmapGrid(res)
        for i in range(res):
            for j in range(res):
                # Start from the middle
                newgrid[i, j] = self[(i+mid)%res, (j+mid)%res]
        return newgrid

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Residue(object):
    """ Residue class """

    def __init__(self, resname, idx):
        self.resname = resname
        self.idx = idx
        self.atoms = []

    def add_atom(self, atom):
        atom.residue = self
        self.atoms.append(atom)

    def delete_atom(self, atom):
        """
        If an atom is present in this residue, delete it from the list of
        atoms
        """
        self.atoms = [a for a in self.atoms if a is not atom]

    # Implement some container methods over the list of atoms
    def __contains__(self, thing):
        """ True if an atom is present in this residue """
        return thing in self.atoms

    def __len__(self):
        return len(self.atoms)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ResidueList(list):
    """ Array of Residues. """

    #===================================================

    def __init__(self, parm):
        list.__init__(self, [Residue(parm.parm_data['RESIDUE_LABEL'][i], i+1)
                             for i in range(parm.ptr('nres'))])
        for i, val in enumerate(parm.parm_data['RESIDUE_POINTER']):
            start = val - 1
            try:
                end = parm.parm_data['RESIDUE_POINTER'][i+1] - 1
            except IndexError:
                end = parm.parm_data['POINTERS'][NATOM]
            for j in range(start, end):
                self[i].add_atom(parm.atom_list[j])
        self.parm = parm
   
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AtomList(list):
    """ Array of Atoms """
    #===================================================

    def __init__(self, parm, fill_from=None):
        self.parm = parm
        if fill_from is None:
            list.__init__(self, [Atom(self.parm, i) for i in
                                range(self.parm.ptr('natom'))])
        else:
            list.__init__(self, [0 for i in range(self.parm.ptr('natom'))])
            for i, atm in enumerate(fill_from): self[i] = atm
        self.changed = False

    #===================================================

    def __delitem__(self, idx):
        """ Deletes this atom then re-indexes everybody else """
        self[idx].idx = -1
        # Delete this atom from its residue as well
        self[idx].residue.delete_atom(self[idx])
        list.__delitem__(self, idx)
        self.changed = True

    #===================================================
   
    def unmark(self):
        """ Unmark all atoms in this list """
        for atm in self: atm.marked = 0

    #===================================================
   
    def _index_us(self):
        """ We have deleted an atom, so now we have to re-index everybody """
        for i, atom in enumerate(self): atom.idx = atom.starting_index = i

    #===================================================

    def append(self, item):
        """ Don't allow this! """
        raise AmberParmError("Cannot add to an AtomList!")

    #===================================================

    def find_extra_exclusions(self):
        " Load all extra exclusions that may be stored in the topology file "
        for atom in self: atom.load_exclusions()

    #===================================================

    def write_to_parm(self):
        """ Writes all of the atom data to the topology file """
        # Write all of the arrays here
        self.parm.parm_data['POINTERS'][NATOM] = len(self)
        # Array slices are faster than copy() and creating new arrays
        # each time
        zeros = [0 for i in range(self.parm.parm_data['POINTERS'][NATOM])]
        self.parm.parm_data['ATOM_NAME'] = zeros[:]
        self.parm.parm_data['CHARGE'] = zeros[:]
        self.parm.parm_data['MASS'] = zeros[:]
        self.parm.parm_data['ATOM_TYPE_INDEX'] = zeros[:]
        self.parm.parm_data['NUMBER_EXCLUDED_ATOMS'] = zeros[:]
        self.parm.parm_data['AMBER_ATOM_TYPE'] = zeros[:]
        self.parm.parm_data['JOIN_ARRAY'] = zeros[:]
        self.parm.parm_data['TREE_CHAIN_CLASSIFICATION'] = zeros[:]
        self.parm.parm_data['IROTAT'] = zeros[:]
        self.parm.parm_data['RADII'] = zeros[:]
        self.parm.parm_data['SCREEN'] = zeros[:]
        self._index_us()
        self._determine_exclusions()
        for atm in self: 
            atm.add_data()
            atm.starting_index = atm.idx # arrays are updated...

    #===================================================

    def _determine_exclusions(self):
        """
        Figures out the EXCLUDED_ATOMS_LIST. Only do this right before you write
        the topology file, since it's expensive
        """
        self.parm.parm_data['EXCLUDED_ATOMS_LIST'] = []
        # We have to do something different for extra points. See the top of
        # extra_pts.f in the sander src/ directory. Effectively, the EP is
        # excluded from every atom that the atom it's attached to is excluded
        # from. So here we go through and exclude every extra point from every
        # other atom my bonded pair is excluded from:
        for atm in self:
            if not atm.attype[:2] in ['EP', 'LP']: continue
            partner = atm.bond_partners[0]
            # Now add all bond partners
            for patm in partner.bond_partners:
                # Don't add myself
                if patm is atm: continue
                atm.exclude(patm)
            # Now add all angle partners
            for patm in partner.angle_partners: atm.exclude(patm)
            # Now add all dihedral partners
            for patm in partner.dihedral_partners: atm.exclude(patm)
            # Now add all other arbitrary exclusions
            for patm in partner.exclusion_partners:
                if patm is atm: continue
                atm.exclude(patm)

        for atm in self:
            vals_to_add = []
            for member in atm.bond_partners:
                if member.idx > atm.idx: vals_to_add.append(member.idx+1)
            for member in atm.angle_partners:
                if member.idx > atm.idx: vals_to_add.append(member.idx+1)
            for member in atm.dihedral_partners:
                if member.idx > atm.idx: vals_to_add.append(member.idx+1)
            # Enable additional (arbitrary) exclusions
            for member in atm.exclusion_partners:
                if member.idx > atm.idx: vals_to_add.append(member.idx+1)
            vals_to_add.sort()
            # See comment above about numex = 0 --> numex = 1
            if not vals_to_add: vals_to_add = [0]
            self.parm.parm_data['EXCLUDED_ATOMS_LIST'].extend(vals_to_add)

    #===================================================

    def refresh_data(self):
        """
        Re-loads all data in parm.parm_data for each atom in case we changed
        any of it
        """
        for atm in self:
            atm.load_from_parm()

    #===================================================

    def __setitem__(self, idx, thing):
        """
        This means we changed things...
        """
        self.changed = True
        list.__setitem__(self, idx, thing)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _tracking(fcn):
    """ Decorator to indicate the list has changed """
    def new_fcn(self, *args):
        self.changed = True
        return fcn(self, *args)
    return new_fcn

class TrackedList(list):
    """
    This creates a list type that allows you to see if anything has changed
    """
    def __init__(self, arg=[]):
        self.changed = False
        list.__init__(self, arg)

    __delitem__ = _tracking(list.__delitem__)
    append = _tracking(list.append)
    extend = _tracking(list.extend)
    __setitem__ = _tracking(list.__setitem__)

# Python 3 does not have __delslice__, but make sure we override it for Python 2
if hasattr(TrackedList, '__delslice__'):
    TrackedList.__delslice__ = _tracking(TrackedList.__delslice__)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _TypeList(TrackedList):
    """ Base class for all type lists """

    #===================================================

    def __init__(self, parm):
        """ Constructs a list of bond types from the topology file """
        self.parm = parm
        self._make_array()
        self.changed = False

    #===================================================

    def _make_array(self):
        """ 
        This method fills self with whichever element we need. This MUST be
        overwritten, so I force it here
        """
        raise NotImplemented('Subclasses must implement _make_array')

    #===================================================

    def write_to_parm(self):
        """ Writes the data here to the parm data """
        for item in self: item.write_info(self.parm)

    #===================================================

    def reset(self):
        """ Reset indexes to -1 to allow multiple remake_parm calls """
        for thing in self: thing.idx = -1

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondTypeList(_TypeList):
    """ Bond type list """

    #===================================================

    def _make_array(self):
        kl = self.parm.parm_data['BOND_FORCE_CONSTANT']
        eql = self.parm.parm_data['BOND_EQUIL_VALUE']
        list.__init__(self, [BondType(k, eq, -1) for k, eq in zip(kl, eql)])
      
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleTypeList(_TypeList):
    """ Angle type list """

    #===================================================

    def _make_array(self):
        kl = self.parm.parm_data['ANGLE_FORCE_CONSTANT']
        eql = self.parm.parm_data['ANGLE_EQUIL_VALUE']
        list.__init__(self, [AngleType(k, eq, -1) for k, eq in zip(kl, eql)])
      
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralTypeList(_TypeList):
    """ Dihedral type list """

    #===================================================

    def _make_array(self):
        kl = self.parm.parm_data['DIHEDRAL_FORCE_CONSTANT']
        pel = self.parm.parm_data['DIHEDRAL_PERIODICITY']
        phl = self.parm.parm_data['DIHEDRAL_PHASE']
        if (not 'SCEE_SCALE_FACTOR' in self.parm.parm_data.keys() or
            not 'SCNB_SCALE_FACTOR' in self.parm.parm_data.keys()):
            list.__init__(self,
                    [DihedralType(k, pe, ph, 1.2, 2.0, -1)
                     for k, pe, ph in zip(kl, pel, phl)]
            )
        else:
            scel = self.parm.parm_data['SCEE_SCALE_FACTOR']
            scnl = self.parm.parm_data['SCNB_SCALE_FACTOR']
            list.__init__(self,
                    [DihedralType(k, pe, ph, sce, scn, -1)
                     for k, pe, ph, sce, scn in zip(kl, pel, phl, scel, scnl)]
            )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradleyTypeList(_TypeList):
    """ Urey-Bradley type list """

    #===================================================

    def _make_array(self):
        kl = self.parm.parm_data['CHARMM_UREY_BRADLEY_FORCE_CONSTANT']
        eql = self.parm.parm_data['CHARMM_UREY_BRADLEY_EQUIL_VALUE']
        list.__init__(self,
                [UreyBradleyType(k, eq, -1) for k, eq in zip(kl, eql)]
        )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ImproperTypeList(_TypeList):
    """ CHARMM Improper torsion type list """

    #===================================================

    def _make_array(self):
        kl = self.parm.parm_data['CHARMM_IMPROPER_FORCE_CONSTANT']
        pl = self.parm.parm_data['CHARMM_IMPROPER_PHASE']
        list.__init__(self, [ImproperType(k, p, -1) for k, p in zip(kl, pl)])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CmapTypeList(_TypeList):
    """ CHARMM correction-map type list """

    #===================================================

    def _make_array(self):
        # Need to get all of the CMAP types
        ncmaps = self.parm.parm_data['CHARMM_CMAP_COUNT'][1]
        list.__init__(self)
        for i in range(ncmaps):
            res = self.parm.parm_data['CHARMM_CMAP_RESOLUTION'][i]
            grid = self.parm.parm_data['CHARMM_CMAP_PARAMETER_%02d' % (i+1)]
            cmts = self.parm.parm_comments['CHARMM_CMAP_PARAMETER_%02d' % (i+1)]
            list.append(self, CmapType(res, grid, cmts, -1))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def Element(mass):
    """ Determines what element the given atom is based on its mass """

    diff = mass
    best_guess = 'EP'

    for element in _Element:
        if abs(Mass[element] - mass) < diff:
            best_guess = element
            diff = abs(Mass[element] - mass)

    return best_guess

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == '__main__':
    import doctest
    doctest.testmod()
