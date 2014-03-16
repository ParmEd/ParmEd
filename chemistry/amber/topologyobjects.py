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
      self.bond_partners = set()
      self.angle_partners = set()
      self.dihedral_partners = set()
      self.exclusion_partners = set() # For arbitrary exclusions
      self.parm = parm
      self.idx = -1
      self.starting_index = starting_index
      self.load_from_parm()
      self.residue = None
      self.marked = 0 # For setting molecules
      self.deleted = False
   
   #===================================================
   
   def bonds(self):
      """ Go through all bonded partners """
      return sorted(list(self.bond_partners))

   def angles(self):
      """ List of all angle partners that are NOT bond partners """
      return sorted(list(self.angle_partners - self.bond_partners))

   def dihedrals(self):
      """ List of all dihedral partners that are NOT angle or bond partners """
      return sorted(list(self.dihedral_partners - self.angle_partners -
                         self.bond_partners))

   def exclusions(self):
      " List of all exclusions not otherwise excluded by bonds/angles/torsions "
      return sorted(list(self.exclusion_partners - self.dihedral_partners -
                         self.angle_partners - self.bond_partners))

   #===================================================

   def add_data(self):
      """ 
      Writes this atom's data to the AmberParm object. Don't pitch a fit if
      we're missing some useless (unused) flags.
      """
      # Determine how many excluded atoms we have. The only ones that count are
      # those with a smaller index (to avoid double-counting)
      numex = 0
      for atm in self.bonds():
         if atm.idx > self.idx: numex += 1
      for atm in self.angles():
         if atm.idx > self.idx: numex += 1
      for atm in self.dihedrals():
         if atm.idx > self.idx: numex += 1
      for atm in self.exclusions():
         if atm.idx > self.idx: numex += 1
      # For some reason, existing topology files follow the convention that
      # atoms with no exclusions (because all bonded partners have atom #s lower
      # than theirs) have num_excluded = 1, with a 0 placeholder in
      # EXCLUDED_ATOMS_LIST... Weird.
      if numex == 0: numex = 1
      # Make sure we're indexing from 0
      self.parm.parm_data['ATOM_NAME'][self.idx] = self.atname[:4]
      self.parm.parm_data['CHARGE'][self.idx] = self.charge
      self.parm.parm_data['MASS'][self.idx] = self.mass
      self.parm.parm_data['ATOM_TYPE_INDEX'][self.idx] = self.nb_idx
      self.parm.parm_data['NUMBER_EXCLUDED_ATOMS'][self.idx] = numex
      self.parm.parm_data['AMBER_ATOM_TYPE'][self.idx] = self.attype[:4]
      self.parm.parm_data['JOIN_ARRAY'][self.idx] = 0
      self.parm.parm_data['TREE_CHAIN_CLASSIFICATION'][self.idx] = self.tree[:4]
      self.parm.parm_data['IROTAT'][self.idx] = 0
      self.parm.parm_data['RADII'][self.idx] = self.radii
      self.parm.parm_data['SCREEN'][self.idx] = self.screen
      try:
         self.parm.parm_data['ATOMIC_NUMBER'][self.idx] = self.atomic_number
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
      self.atname = self.parm.parm_data['ATOM_NAME'][self.starting_index]
      self.charge = self.parm.parm_data['CHARGE'][self.starting_index]
      self.mass = self.parm.parm_data['MASS'][self.starting_index]
      self.nb_idx = self.parm.parm_data['ATOM_TYPE_INDEX'][self.starting_index]
      self.attype = self.parm.parm_data['AMBER_ATOM_TYPE'][self.starting_index]
      self.tree = self.parm.parm_data['TREE_CHAIN_CLASSIFICATION'] \
                                                          [self.starting_index]
      # Put in some 
      if 'RADII' in self.parm.parm_data:
         self.radii = self.parm.parm_data['RADII'][self.starting_index]
      else:
         self.radii = 0.0 # dummy number
      if 'SCREEN' in self.parm.parm_data:
         self.screen = self.parm.parm_data['SCREEN'][self.starting_index]
      else:
         self.radii = 0.0 # dummy number
      if 'ATOMIC_NUMBER' in self.parm.parm_data:
         self.atomic_number = self.parm.parm_data['ATOMIC_NUMBER'] \
                                                          [self.starting_index]
      else:
         # Determine from mass (unless already done)
         if self.atomic_number <= 0:
            self.atomic_number = AtomicNum[Element(self.mass)]

      # Load the positions and velocities if the amberParm object them
      if hasattr(self.parm, 'coords'):
         self.xx = self.parm.coords[self.starting_index*3  ]
         self.xy = self.parm.coords[self.starting_index*3+1]
         self.xz = self.parm.coords[self.starting_index*3+2]
         if self.parm.hasvels:
            self.vx = self.parm.vels[self.starting_index*3  ]
            self.vy = self.parm.vels[self.starting_index*3+1]
            self.vz = self.parm.vels[self.starting_index*3+2]

   #===================================================
      
   def bond_to(self, other):
      """ 
      Log this atom as bonded to another atom. Check if this has already been
      added to the angle list. If so, remove it from there.
      """
      if self is other:
         raise BondError("Cannot bond atom to itself!")
      self.bond_partners.add(other)

   #===================================================
      
   def angle_to(self, other):
      """
      Log this atom as angled to another atom. Check if this has already been
      added to the bond list. If so, do nothing
      """
      if self is other:
         raise BondError("Cannot angle an atom with itself!")
      self.angle_partners.add(other)
   
   #===================================================

   def dihedral_to(self, other):
      """
      Log this atom as dihedral-ed to another atom. Check if this has already been
      added to the bond or angle list. If so, do nothing
      """
      if self is other:
         raise BondError("Cannot dihedral an atom with itself!")
      self.dihedral_partners.add(other)
      
   #===================================================

   def exclude(self, other):
      """
      Add one atom to my exclusion list, even if it's not bonded to, angled to, or
      dihedraled to it
      """
      if self is other:
         raise BondError("Cannot exclude an atom from itself")
      self.exclusion_partners.add(other)
      # If he is excluded from me, then I am excluded from him
      other.exclusion_partners.add(self)

   #===================================================

   def reset_topology(self):
      """
      Deletes all of the bond, angle, and dihedral partners so they can be set
      up again with updated data. Keep the arbitrary exclusions, though
      """
      self.bond_partners = set()
      self.angle_partners = set()
      self.dihedral_partners = set()

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
      # Assume it's a bond
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
   " Dihedral class with 4 atoms involved and force constant/periodicity/phase "
      
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
         return ((self.atom1 in thing and self.atom2 in thing) or
                 (self.atom2 in thing and self.atom3 in thing) or
                 (self.atom3 in thing and self.atom4 in thing))
      # Otherwise assume thing is an Atom
      return (thing is self.atom1 or thing is self.atom2 or
              thing is self.atom3 or thing is self.atom4)

   #===================================================

   def __eq__(self, thing):
      """
      A dihedral is equivalent if the 4 atoms are the same (or reverse) in order

      Allow comparison with another type of dihedral or with a list of 4 atoms
      (or tuple)
      """
      if isinstance(thing, Dihedral):
         # I'm comparing with another Dihedral here
         return ( (self.atom1 is thing.atom1 and self.atom2 is thing.atom2 and
                   self.atom3 is thing.atom3 and self.atom4 is thing.atom4) or
                  (self.atom1 is thing.atom4 and self.atom2 is thing.atom3 and
                   self.atom4 is thing.atom1) )
      if isinstance(thing, list) or isinstance(thing, tuple):
         # Here, atoms are expected to index from 0 (Python standard) if we
         # are comparing with a list or tuple
         if len(thing) != 4:
            raise DihedralError('comparative %s has %d elements! Expect 4.' % 
                                (type(thing).__name__, len(thing)))
         # Compare starting_index, since we may not have an index right now...
         return ( (self.atom1.starting_index == thing[0] and 
                   self.atom2.starting_index == thing[1] and
                   self.atom3.starting_index == thing[2] and
                   self.atom4.starting_index == thing[3]) or
                  (self.atom1.starting_index == thing[3] and 
                   self.atom2.starting_index == thing[2] and
                   self.atom3.starting_index == thing[1] and
                   self.atom4.starting_index == thing[0]) )

      raise TypeError('Cannot compare Dihedral with %s' % type(thing).__name__)

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
      # If this is a bond things are a bit more complicated since we don't know
      # the central atom if this Urey-Bradley. We need to make sure that one of
      # the atoms of the bond is either atom1 or atom2 and that the OTHER atom
      # in Bond "thing" has the OTHER atom in this Urey-Bradley in its list of
      # bonded partners.
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
         # If we are here, end1 and cent are set. Look through the bond partners
         # of cent(er) and see if any of them is in this Urey-Bradley (but ONLY
         # if that atom is not the original end1)
         for atm in cent.bonds():
            if atm is end1: continue
            if atm in self:
               # If we got here, we found both atoms in this Urey-Bradley
               # separated by 2 bonds, so "thing" IS in this Urey-Bradley
               return True
         # If we got here, we could not find the other atom connected by 2 bonds
         return False
      # If we are here, "thing" must not be a Bond and we therefore assume it's
      # an Atom
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
   " Improper class with 4 atoms involved and force constant/phase "
      
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
      " Quick and easy way to find out if an Atom or Bond is in this Improper "
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
      A dihedral is equivalent if the 4 atoms are the same (or reverse) in order

      Allow comparison with another type of dihedral or with a list of 4 atoms
      (or tuple)
      """
      if isinstance(thing, Improper):
         # I'm comparing with another Improper here. Central atom must be the
         # same. Others can be in any order
         if self.atom3 != thing.atom3:
            return False
         # Make a set with the remaining atoms. If they are equal, the impropers
         # are equivalent
         selfset = set([self.atom1, self.atom2, self.atom4])
         otherset = set([thing.atom1, thing.atom2, thing.atom3])
         return selfset == otherset
      if isinstance(thing, list) or isinstance(thing, tuple):
         # Here, atoms are expected to index from 0 (Python standard) if we
         # are comparing with a list or tuple
         if len(thing) != 4:
            raise DihedralError('comparative %s has %d elements! Expect 4.' % 
                                (type(thing).__name__, len(thing)))
         if self.atom3.starting_index != thing[2]:
            return False
         selfset = set([self.atom1.starting_index, self.atom2.starting_index,
                        self.atom4.starting_index])
         otherset = set([thing[0], thing[1], thing[3]])
         return selfset == otherset

      raise TypeError('Cannot compare Improper with %s' % type(thing).__name__)

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
      Quick and easy way to find out an atom or bond is in this coupled-torsion
      """
      # If we have a bond, see if it is 1 of the 4 sequential bonds that make up
      # this coupled-torsion
      if isinstance(thing, Bond):
         return ((self.atom1 in thing and self.atom2 in thing) or
                 (self.atom2 in thing and self.atom3 in thing) or
                 (self.atom3 in thing and self.atom4 in thing) or
                 (self.atom4 in thing and self.atom5 in thing))
      # Otherwise assume we have an atom
      return (thing is self.atom1 or thing is self.atom2 or
              thing is self.atom3 or thing is self.atom4 or thing is self.atom5)

   #===================================================

   def __eq__(self, thing):
      """
      A coupled-torsion is equivalent if the 5 atoms are in the same or reverse
      order

      Allow comparison with another type of cmap or with a sequence of 5 indexes
      """
      if isinstance(thing, Cmap):
         # I'm comparing with another Improper here. Central atom must be the
         # same. Others can be in any order
         return ( (self.atom1 is thing.atom1 and self.atom2 is thing.atom2 and
                   self.atom3 is thing.atom3 and self.atom4 is thing.atom4 and
                   self.atom5 is thing.atom5) or
                  (self.atom1 is thing.atom5 and self.atom2 is thing.atom4 and
                   self.atom3 is thing.atom3 and self.atom4 is thing.atom2 and
                   self.atom5 is thing.atom1))
      if isinstance(thing, list) or isinstance(thing, tuple):
         # Here, atoms are expected to index from 0 (Python standard) if we
         # are comparing with a list or tuple
         if len(thing) != 5:
            raise DihedralError('comparative %s has %d elements! Expect 4.' % 
                                (type(thing).__name__, len(thing)))
         return ( (self.atom1.starting_index == thing[0] and
                   self.atom2.starting_index == thing[1] and
                   self.atom3.starting_index == thing[2] and
                   self.atom4.starting_index == thing[3] and
                   self.atom5.starting_index == thing[4]) or
                  (self.atom1.starting_index == thing[4] and
                   self.atom2.starting_index == thing[3] and
                   self.atom3.starting_index == thing[2] and
                   self.atom4.starting_index == thing[1] and
                   self.atom5.starting_index == thing[0]))
      raise TypeError('Cannot compare Improper with %s' % type(thing).__name__)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CmapType(object):
   """ Contains a correction map interpolation grid """

   #===================================================

   def __init__(self, resolution, grid, cmts, idx):
      self.resolution = resolution
      self.grid = grid
      self.comments = cmts # To avoid losing any comments in the chamber prmtop
      if len(grid) != self.resolution * self.resolution:
         raise CmapError('CMAP grid does not match expected resolution')
      self.idx = idx

   #===================================================

   def write_info(self, parm):
      """ Write out the dihedral parameters """
      # If our idx == -1, we're not being used, so just return
      if self.idx == -1: return
      parm.parm_data['CHARMM_CMAP_RESOLUTION'][self.idx] = self.resolution
      parm.parm_data['CHARMM_CMAP_PARAMETER_%02d' % (self.idx+1)] = self.grid[:]

   #===================================================
   
   def __eq__(self, other):
      return (self.resolution == other.resolution and
              all([abs(i - j) < TINY for i, j in zip(self.grid, other.grid)]))

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
      list.__delitem__(self, idx)
      self.changed = True

   #===================================================
   
   def unmark(self):
      """ Unmark all atoms in this list """
      for atm in self: atm.marked = 0

   #===================================================
   
   def _index_us(self):
      """ We have deleted an atom, so now we have to re-index everybody """
      for i in range(len(self)): self[i].idx = self[i].starting_index = i

   #===================================================

   def append(self, item):
      """ Don't allow this! """
      raise AmberParmError("Cannot add to an AtomList!")

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
         partner = atm.bonds()[0]
         # Now add all bond partners
         for patm in partner.bonds():
            # Don't add myself
            if patm is atm: continue
            atm.exclude(patm)
         # Now add all angle partners
         for patm in partner.angles(): atm.exclude(patm)
         # Now add all dihedral partners
         for patm in partner.dihedrals(): atm.exclude(patm)
         # Now add all other arbitrary exclusions
         for patm in partner.exclusions():
            if patm is atm: continue
            atm.exclude(patm)

      for atm in self:
         vals_to_add = []
         for member in atm.bonds():
            if member.idx > atm.idx: vals_to_add.append(member.idx+1)
         for member in atm.angles():
            if member.idx > atm.idx: vals_to_add.append(member.idx+1)
         for member in atm.dihedrals():
            if member.idx > atm.idx: vals_to_add.append(member.idx+1)
         # Enable additional (arbitrary) exclusions
         for member in atm.exclusions():
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
   extend = _tracking(list.append)
   __delslice__ = _tracking(list.__delslice__)
   __setitem__ = _tracking(list.__setitem__)

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

   def deindex(self):
      """ Resets all of the type indexes to -1 """
      for item in self: item.idx = -1

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondTypeList(_TypeList):
   """ Bond type list """

   #===================================================

   def _make_array(self):
      list.__init__(self,
                    [BondType(self.parm.parm_data['BOND_FORCE_CONSTANT'][i],
                               self.parm.parm_data['BOND_EQUIL_VALUE'][i], -1)
                               for i in range(self.parm.ptr('numbnd'))]
                   )
      
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleTypeList(_TypeList):
   """ Angle type list """

   #===================================================

   def _make_array(self):
      list.__init__(self,
                    [AngleType(self.parm.parm_data['ANGLE_FORCE_CONSTANT'][i],
                                self.parm.parm_data['ANGLE_EQUIL_VALUE'][i], -1)
                                for i in range(self.parm.ptr('numang')) ])
      
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralTypeList(_TypeList):
   """ Dihedral type list """

   #===================================================

   def _make_array(self):
      if not 'SCEE_SCALE_FACTOR' in self.parm.parm_data.keys() or \
         not 'SCNB_SCALE_FACTOR' in self.parm.parm_data.keys():
         list.__init__(self,
              [DihedralType(self.parm.parm_data['DIHEDRAL_FORCE_CONSTANT'][i],
               self.parm.parm_data['DIHEDRAL_PERIODICITY'][i], 
               self.parm.parm_data['DIHEDRAL_PHASE'][i], 1.2, 2.0, -1)
               for i in range(self.parm.ptr('nptra')) ])
      else:
         list.__init__(self,
              [DihedralType(self.parm.parm_data['DIHEDRAL_FORCE_CONSTANT'][i],
               self.parm.parm_data['DIHEDRAL_PERIODICITY'][i], 
               self.parm.parm_data['DIHEDRAL_PHASE'][i], 
               self.parm.parm_data['SCEE_SCALE_FACTOR'][i],
               self.parm.parm_data['SCNB_SCALE_FACTOR'][i], -1)
               for i in range(self.parm.ptr('nptra')) ])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradleyTypeList(_TypeList):
   """ Urey-Bradley type list """

   #===================================================

   def _make_array(self):
      klist = self.parm.parm_data['CHARMM_UREY_BRADLEY_FORCE_CONSTANT']
      eqlist = self.parm.parm_data['CHARMM_UREY_BRADLEY_EQUIL_VALUE']
      list.__init__(self, [UreyBradleyType(k, eq, -1)
                     for k, eq in zip(klist, eqlist)])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ImproperTypeList(_TypeList):
   """ CHARMM Improper torsion type list """

   #===================================================

   def _make_array(self):
      klist = self.parm.parm_data['CHARMM_IMPROPER_FORCE_CONSTANT']
      plist = self.parm.parm_data['CHARMM_IMPROPER_PHASE']
      list.__init__(self, [ImproperType(k, eq, -1)
                     for k, eq in zip(klist, plist)])

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
