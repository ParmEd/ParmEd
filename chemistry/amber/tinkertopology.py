"""
This module contains objects that deal with system topology and parameters
specifically for the Amoeba force field.

These are intended for use with the TinkerParm class and is used to recompute
the parameter topology file when necessary (i.e., it provides high-level access
to topology file data)

Designed so it can be used like:
   from chemistry.amber.tinkertopology import *

by Jason Swails
"""

from compat24 import property
from chemistry.amber.constants import NATOM
from chemistry.exceptions import (BondError, AmberParmError, AmoebaError)
from chemistry.amber.topologyobjects import (_TypeList, TrackedList,
            ResidueList as _ResidueList,
            Bond as _Bond,
            UreyBradley as _UreyBradley,
            Angle as _Angle,
            UreyBradleyType as _UreyBradleyType,
            AngleType as _AngleType,
            Cmap as _Cmap)

class Atom(object):
   """
   An atom. Only use these as elements in AtomList instances, since AtomList
   will keep track of when indexes and other stuff needs to be updated.
   """

   #===================================================

   def __init__(self, parm, si):
      self.starting_index = si
      self.parm = parm
      self.bond_partners = set()
      self.angle_partners = set()
      self.dihedral_partners = set()
      self.tortor_partners = set()
      self.exclusion_partners = set()
      self.idx = -1
      self.load_from_parm()

   #===================================================

   def load_from_parm(self):
      parm = self.parm
      si = self.starting_index
      self.atomic_number = parm.parm_data['AMOEBA_ATOMIC_NUMBER'][si]
      self.atname = parm.parm_data['ATOM_NAME'][si]
      self.mass = parm.parm_data['MASS'][si]
      self.tree = parm.parm_data['TREE_CHAIN_CLASSIFICATION'][si]
      self.attype = parm.parm_data['AMBER_ATOM_TYPE'][si]
      self.typeidx = parm.parm_data['AMOEBA_ATOM_TYPE_INDEX'][si]
      self.classidx = parm.parm_data['AMOEBA_ATOM_CLASS_INDEX'][si]
      self.nb_idx = parm.parm_data['AMOEBA_VDW_ATOM_TYPES_LIST'][si]
      self._vdw_par = parm.parm_data['AMOEBA_VDW_ATOM_PARENT_LIST'][si]
      self.vdw_par_wt=parm.parm_data['AMOEBA_VDW_PARENT_COORD_WEIGHT_LIST'][si]
      self.multipoles = parm.parm_data['AMOEBA_LOCAL_FRAME_MULTIPOLES_LIST'] \
                                             [10*si:10*si+10]
      self.polar = parm.parm_data['AMOEBA_POLARIZABILITY_LIST'][si]
      if hasattr(self.parm, 'coords'):
         self.xx = self.parm.coords[si*3  ]
         self.xy = self.parm.coords[si*3+1]
         self.xz = self.parm.coords[si*3+2]
         if self.parm.hasvels:
            self.vx = self.parm.vels[si*3  ]
            self.vy = self.parm.vels[si*3+1]
            self.vz = self.parm.vels[si*3+2]
   
   #===================================================

   # Make 'element' an alias for 'atomic_number'

   @property
   def element(self):
      return self.atomic_number
   @element.setter
   def element(self, value):
      self.atomic_number = value

   #===================================================
   
   def __repr__(self):
      """ For debugging purposes -- print out where I am """
      return '<Atom %d [%s]; In residue %d %s>' % (self.starting_index + 1,
                  self.atname, self.residue.idx, self.residue.resname)

   #===================================================

   def add_data(self):
      """ Writes this atom's data to the TinkerParm object """
      # Determine how many excluded atoms we have. The only ones that count are
      # those with a smaller index (to avoid double-counting)
      numex = 0
      for atm in self.bonds():
         if atm.idx > self.idx: numex += 1
      for atm in self.angles():
         if atm.idx > self.idx: numex += 1
      for atm in self.dihedrals():
         if atm.idx > self.idx: numex += 1
      for atm in self.tortors():
         if atm.idx > self.idx: numex += 1
      for atm in self.exclusions():
         if atm.idx > self.idx: numex += 1
      # Make sure we are indexing from 0
      p = self.parm
      p.parm_data['ATOM_NAME'][self.idx] = self.atname[:4]
      p.parm_data['NUMBER_EXCLUDED_ATOMS'][self.idx] = numex
      p.parm_data['CHARGE'][self.idx] = 0.0
      p.parm_data['MASS'][self.idx] = self.mass
      p.parm_data['TREE_CHAIN_CLASSIFICATION'][self.idx] = self.tree
      p.parm_data['ATOM_TYPE_INDEX'][self.idx] = 1 # hard-coded
      p.parm_data['AMBER_ATOM_TYPE'][self.idx] = self.attype
      p.parm_data['AMOEBA_ATOMIC_NUMBER'][self.idx] = self.atomic_number
      p.parm_data['AMOEBA_ATOM_TYPE_INDEX'][self.idx] = self.typeidx
      p.parm_data['AMOEBA_ATOM_CLASS_INDEX'][self.idx] = self.classidx
      p.parm_data['AMOEBA_VDW_ATOM_TYPES_LIST'][self.idx] = self.nb_idx
      p.parm_data['AMOEBA_VDW_ATOM_PARENT_LIST'][self.idx] = self.vdw_par.idx+1
      p.parm_data['AMOEBA_VDW_PARENT_COORD_WEIGHT_LIST'][self.idx] = \
                                                               self.vdw_par_wt
      key = 'AMOEBA_LOCAL_FRAME_MULTIPOLES_LIST'
      for i in range(10):
         idx = 10 * self.idx + i
         p.parm_data[key][idx] = self.multipoles[i]
      p.parm_data['AMOEBA_POLARIZABILITY_LIST'][self.idx] = self.polar

   #===================================================

   def reset_topology(self):
      """ Resets all of the bonded partners and such """
      self.bond_partners = set()
      self.angle_partners = set()
      self.dihedral_partners = set()
      self.tortor_partners = set()

   #===================================================
   
   def bond_to(self, other):
      """ Log these atoms as bonded to each other """
      if self is other:
         raise BondError('Cannot bond atom to itself.')
      self.bond_partners.add(other)
      other.bond_partners.add(self)

   def angle_to(self, other):
      """ Log the other atom as angled to each other """
      if self is other:
         raise BondError('Cannot angle atom to itself.')
      self.angle_partners.add(other)
      other.angle_partners.add(self)

   def dihedral_to(self, other):
      """ Log the other atom as dihedraled to each other """
      if self is other:
         raise BondError('Cannot dihedral atom to itself.')
      self.dihedral_partners.add(other)
      other.dihedral_partners.add(self)

   def tortor_to(self, other):
      """ Log atoms as coupled-torsion partners """
      if self is other:
         raise BondError('Cannot coupled-dihedral atom to itself.')
      self.tortor_partners.add(other)
      other.tortor_partners.add(self)

   def exclude(self, other):
      """ Exclude atoms from each other """
      if self is other:
         raise BondError('Cannot exclude atom from itself')
      self.exclusion_partners.add(other)
      other.exclusion_partners.add(self)

   #===================================================
   
   def determine_all_exclusion_partners(self):
      """
      Amoeba topologies have complex exclusions due to the existence of
      polarization groups and the exclusions of 1-2, 1-3, 1-4, and 1-5 pairs.
      This method looks through the exclusion list and adds any atoms that are
      not naturally accounted for by the bond, angle, torsion, and
      coupled-torsion parameters
      """
      self.determine_exclusions_from_bonds()
      excset = set()
      exclat = self.parm.parm_data['NUMBER_EXCLUDED_ATOMS']
      exclist = self.parm.parm_data['EXCLUDED_ATOMS_LIST']
      first_excl = sum(exclat[:self.starting_index])
      nexcl = exclat[self.starting_index]
      for i in range(nexcl):
         idx = exclist[first_excl+i] - 1
         excset.add(self.parm.atom_list[idx])
      # Now subtract off all of the bonds, angles, torsions, and tortors
      excset = (excset - self.bond_partners - self.angle_partners -
                self.dihedral_partners - self.tortor_partners)
      for atm in excset:
         self.exclude(atm)

   #===================================================
   
   def determine_exclusions_from_bonds(self):
      """
      Add all atoms that are connected by 2, 3, and 4 bonds to the relevant
      exclusion list
      """
      for a1 in self.bonds():
         # Now loop through every bond of every connected atom, and make sure it
         # is part of the angle partners
         for a2 in a1.bonds():
            # All of these atoms are angle partners
            if a2 is not self:
               self.angle_partners.add(a2)
            for a3 in a2.bonds():
               # All of these atoms are dihedral partners
               if a3 is not self:
                  self.dihedral_partners.add(a3)
               for a4 in a3.bonds():
                  # All of these atoms are tortor partners
                  if a4 is not self:
                     self.tortor_partners.add(a4)
   
   #===================================================

   def register_vdw_parent(self, atom_list):
      """ Sets the vdw_par from the atom list based on the parent index """
      self.vdw_par = atom_list[self._vdw_par-1]

   #===================================================
   
   # Iterators

   def bonds(self):
      return sorted(list(self.bond_partners))

   def angles(self):
      return sorted(list(self.angle_partners - self.bond_partners))

   def dihedrals(self):
      return sorted(list(self.dihedral_partners - self.angle_partners -
                         self.bond_partners))

   def tortors(self):
      return sorted(list(self.tortor_partners - self.dihedral_partners -
                         self.angle_partners - self.bond_partners))

   def exclusions(self):
      return sorted(list(self.exclusion_partners - self.tortor_partners -
             self.dihedral_partners - self.angle_partners - self.bond_partners))

   #===================================================

   def __eq__(self, other):
      return id(self) == id(other)
      
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

class AtomList(list):
   """ Array of Atoms """
   #===================================================

   def __init__(self, parm, fill_from=None):
      self.parm = parm
      natom = self.parm.parm_data['POINTERS'][NATOM]
      if fill_from is None:
         list.__init__(self, [Atom(self.parm, i) for i in range(natom)])
      else:
         list.__init__(self, [0 for i in range(natom)])
         for i, atm in enumerate(fill_from): self[i] = atm
      self.changed = False
      # Register the vdW parent atom here, since we have to make sure that every
      # atom exists before we try to assign the parent atoms (since we can't be
      # assured that the necessary atom will already exist if we assign this at
      # list construction time)
      for atm in self: atm.register_vdw_parent(self)

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
      for i in range(len(self)): self[i].idx = i

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
      # Amoeba sections now
      self.parm.parm_data['AMEOBA_ATOM_TYPE_INDEX'] = zeros[:]
      self.parm.parm_data['AMOEBA_ATOM_CLASS_INDEX'] = zeros[:]
      self.parm.parm_data['AMOEBA_VDW_ATOM_TYPES_LIST'] = zeros[:]
      self.parm.parm_data['AMOEBA_VDW_ATOM_PARENT_LIST'] = zeros[:]
      self.parm.parm_data['AMOEBA_POLARIZABILITY_LIST'] = zeros[:]
      self.parm.parm_data['AMOEBA_LOCAL_FRAME_MULTIPOLES_LIST'] = zeros * 10
      self.parm.parm_data['AMOEBA_POLARIZABILITY_NUM_LIST'] = [len(self)]
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
      # The Amoeba force field does not have lone pairs. No need to determine EP
      # exclusions here like we did in topologyobjects.py
      for atm in self:
         vals_to_add = []
         for member in atm.bonds():
            if member.idx > atm.idx: vals_to_add.append(member.idx+1)
         for member in atm.angles():
            if member.idx > atm.idx: vals_to_add.append(member.idx+1)
         for member in atm.dihedrals():
            if member.idx > atm.idx: vals_to_add.append(member.idx+1)
         for member in atm.tortors():
            if member.idx > atm.idx: vals_to_add.append(member.idx+1)
         for member in atm.exclusions():
            if member.idx > atm.idx: vals_to_add.append(member.idx+1)
         # Nothing to do? Return
         if not vals_to_add:
            continue
         vals_to_add.sort()
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

class Bond(_Bond):
   """ No real differences from the standard Bond class """

   def write_info(self, parm, idx):
      """ Writes the bond info to the topology file. idx starts at 0 """
      key = 'AMOEBA_REGULAR_BOND_LIST'
      # Atom indexes are not multiplied by 3 here.
      parm.parm_data[key][3*idx  ] = self.atom1.idx + 1
      parm.parm_data[key][3*idx+1] = self.atom2.idx + 1
      parm.parm_data[key][3*idx+2] = self.bond_type.idx + 1

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondType(object):
   """ An AMOEBA bond type """

   def __init__(self, k, req, idx):
      self.k = k
      self.req = req
      self.idx = int(idx)

   def write_info(self, parm):
      """ Writes the bond parameters to the parameter file """
      parm.parm_data['AMOEBA_REGULAR_BOND_FORCE_CONSTANT'][self.idx] = self.k
      parm.parm_data['AMOEBA_REGULAR_BOND_EQUIL_VALUE'][self.idx] = self.req

   def __eq__(self, other):
      return self.k == other.k and self.req == other.req

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondTypeList(_TypeList):

   def __init__(self, parm):
      _TypeList.__init__(self, parm)
      self.degree = parm.parm_data['AMOEBA_REGULAR_BOND_FTAB_DEGREE'][0]
      self.coeffs = parm.parm_data['AMOEBA_REGULAR_BOND_FTAB_COEFFS'][:]
      if len(self.coeffs) != self.degree + 1:
         raise AmoebaError('Bond degree (%d) does not make sense with %d '
                           'coefficients.' % (self.degree, len(self.coeffs)))

   def _make_array(self):
      klist = self.parm.parm_data['AMOEBA_REGULAR_BOND_FORCE_CONSTANT']
      eqlist = self.parm.parm_data['AMOEBA_REGULAR_BOND_EQUIL_VALUE']
      list.__init__(self, [BondType(k, eq, -1) for k, eq in zip(klist, eqlist)])

   def write_to_parm(self):
      self.parm.parm_data['AMOEBA_REGULAR_BOND_FTAB_DEGREE'][0] = self.degree
      self.parm.parm_data['AMOEBA_REGULAR_BOND_FTAB_COEFFS'] = self.coeffs[:]
      super(BondTypeList, self).write_to_parm()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradley(_UreyBradley):
   """ An AMOEBA Urey-Bradley term """

   def write_info(self, parm, idx):
      """ Writes the Urey-Bradley parameters to the paramter file """
      _UreyBradley.write_info(self, parm, 'AMOEBA_UREY_BRADLEY_BOND_LIST', idx)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradleyType(_UreyBradleyType):
   
   def write_info(self, parm):
      """ Writes the Urey-Bradley parameters to the parameter file """
      # If our idx is -1 (not being used), just return
      if self.idx == -1: return
      parm.parm_data['AMOEBA_UREY_BRADLEY_BOND_FORCE_CONSTANT'][self.idx] = \
            self.k
      parm.parm_data['AMOEBA_UREY_BRADLEY_BOND_EQUIL_VALUE'][self.idx] = \
            self.req

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradleyTypeList(_TypeList):

   def __init__(self, parm):
      _TypeList.__init__(self, parm)
      self.degree = parm.parm_data['AMOEBA_UREY_BRADLEY_BOND_FTAB_DEGREE'][0]
      self.coeffs = parm.parm_data['AMOEBA_UREY_BRADLEY_BOND_FTAB_COEFFS'][:]
      if len(self.coeffs) != self.degree + 1:
         raise AmoebaError('Urey-Bradley degree (%d) does not make sense with '
                           '%d coefficients.' % (self.degree, len(self.coeffs)))

   def _make_array(self):
      klist = self.parm.parm_data['AMOEBA_UREY_BRADLEY_BOND_FORCE_CONSTANT']
      eqlist = self.parm.parm_data['AMOEBA_UREY_BRADLEY_BOND_EQUIL_VALUE']
      list.__init__(self, [UreyBradleyType(k, eq, -1)
                           for k, eq in zip(klist, eqlist)])

   def write_to_parm(self):
      self.parm.parm_data['AMOEBA_UREY_BRADLEY_BOND_FTAB_DEGREE'][0] = \
                           self.degree
      self.parm.parm_data['AMOEBA_UREY_BRADLEY_BOND_FTAB_COEFFS'] = \
                           self.coeffs[:]
      super(UreyBradleyTypeList, self).write_to_parm()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Angle(_Angle):
   """ A valence angle term """

   def write_info(self, parm, idx):
      key = 'AMOEBA_REGULAR_ANGLE_LIST'
      # Atom indexes are not multiplied by 3 here.
      parm.parm_data[key][4*idx  ] = self.atom1.idx + 1
      parm.parm_data[key][4*idx+1] = self.atom2.idx + 1
      parm.parm_data[key][4*idx+2] = self.atom3.idx + 1
      parm.parm_data[key][4*idx+3] = self.angle_type.idx + 1

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleType(_AngleType):
   """ A valence angle type """

   def write_info(self, parm):
      """ Writes the angle parameters in the file """
      # Do nothing if it's not being used (i.e., idx==-1)
      if self.idx == -1: return
      parm.parm_data['AMOEBA_REGULAR_ANGLE_FORCE_CONSTANT'][self.idx] = self.k
      parm.parm_data['AMOEBA_REGULAR_ANGLE_EQUIL_VALUE'][self.idx] = self.theteq

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleTypeList(_TypeList):

   def __init__(self, parm):
      _TypeList.__init__(self, parm)
      self.degree = parm.parm_data['AMOEBA_REGULAR_ANGLE_FTAB_DEGREE'][0]
      self.coeffs = parm.parm_data['AMOEBA_REGULAR_ANGLE_FTAB_COEFFS'][:]
      if len(self.coeffs) != self.degree + 1:
         raise AmoebaError('Angle degree (%d) does not make sense with %d '
                           'coefficients.' % (self.degree, len(self.coeffs)))

   def _make_array(self):
      klist = self.parm.parm_data['AMOEBA_REGULAR_ANGLE_FORCE_CONSTANT']
      eqlist = self.parm.parm_data['AMOEBA_REGULAR_ANGLE_EQUIL_VALUE']
      list.__init__(self, [AngleType(k, eq, -1)
                           for k, eq in zip(klist, eqlist)])

   def write_to_parm(self):
      self.parm.parm_data['AMOEBA_REGULAR_ANGLE_FTAB_DEGREE'][0] = self.degree
      self.parm.parm_data['AMOEBA_REGULAR_ANGLE_FTAB_COEFFS'] = self.coeffs[:]
      super(AngleTypeList, self).write_to_parm()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _FourAtomTerm(object):
   """ A base class for a parameter that spans 4 atoms """

   def __init__(self, atom1, atom2, atom3, atom4):
      self.atom1 = atom1
      self.atom2 = atom2
      self.atom3 = atom3
      self.atom4 = atom4

   def write_info(self, parm, key, idx):
      if idx == -1: return
      parm.parm_data[key][5*idx  ] = self.atom1.idx + 1
      parm.parm_data[key][5*idx+1] = self.atom2.idx + 1
      parm.parm_data[key][5*idx+2] = self.atom3.idx + 1
      parm.parm_data[key][5*idx+3] = self.atom4.idx + 1

   def __contains__(self, thing):
      return (self.atom1 is thing or self.atom2 is thing or
              self.atom3 is thing or self.atom4 is thing)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TrigonalAngle(_FourAtomTerm):
   """ In-plane angle term """

   def __init__(self, atom1, atom2, atom3, atom4, trigang_type):
      _FourAtomTerm.__init__(self, atom1, atom2, atom3, atom4)
      self.trigang_type = trigang_type

   def write_info(self, parm, idx):
      # If unused, skip this parameter
      if idx == -1: return
      key = 'AMOEBA_TRIGONAL_ANGLE_LIST'
      _FourAtomTerm.write_info(self, parm, key, idx)
      parm.parm_data[key][5*idx+4] = self.trigang_type.idx + 1

   def __eq__(self, other):
      if (self.atom1 is other.atom1 and self.atom2 is other.atom2 and
          self.atom3 is other.atom3 and self.atom4 is other.atom4):
         return self.trigang_type == other.trigang_type
      if (self.atom3 is other.atom1 and self.atom2 is other.atom2 and
          self.atom3 is other.atom1 and self.atom4 is other.atom4):
         return self.trigang_type == other.trigang_type
      return False

   def __contains__(self, thing):
      if isinstance(thing, Bond):
         return ((self.atom1 in thing and self.atom2 in thing) or
                 (self.atom2 in thing and self.atom3 in thing) or
                 (self.atom2 in thing and self.atom4 in thing))
      return super(TrigonalAngle, self).__contains__(thing)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TrigonalAngleType(_AngleType):

   def write_info(self, parm):
      """ Writes the trigonal angle parameters in the file """
      # Do nothing if it's not being used
      if self.idx == -1: return
      parm.parm_data['AMOEBA_TRIGONAL_ANGLE_FORCE_CONSTANT'][self.idx] = self.k
      parm.parm_data['AMOEBA_TRIGONAL_ANGLE_EQUIL_VALUE'][self.idx]=self.theteq

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TrigonalAngleTypeList(_TypeList):

   def __init__(self, parm):
      _TypeList.__init__(self, parm)
      self.degree = parm.parm_data['AMOEBA_TRIGONAL_ANGLE_FTAB_DEGREE'][0]
      self.coeffs = parm.parm_data['AMOEBA_TRIGONAL_ANGLE_FTAB_COEFFS'][:]
      if len(self.coeffs) != self.degree + 1:
         raise AmoebaError('Trigonal angle degree (%d) does not make sense with'
                          ' %d coefficients.' % (self.degree, len(self.coeffs)))
   
   def _make_array(self):
      klist = self.parm.parm_data['AMOEBA_TRIGONAL_ANGLE_FORCE_CONSTANT']
      eqlist = self.parm.parm_data['AMOEBA_TRIGONAL_ANGLE_EQUIL_VALUE']
      list.__init__(self, [TrigonalAngleType(k, eq, -1)
                           for k, eq in zip(klist, eqlist)])

   def write_to_parm(self):
      self.parm.parm_data['AMOEBA_TRIGONAL_ANGLE_FTAB_DEGREE'][0] = self.degree
      self.parm.parm_data['AMOEBA_TRIGONAL_ANGLE_FTAB_COEFFS'] = self.coeffs[:]
      super(TrigonalAngleTypeList, self).write_to_parm()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class OutOfPlaneBend(_FourAtomTerm):
   """ Out-of-plane bending term """

   def __init__(self, atom1, atom2, atom3, atom4, oopbend_type):
      _FourAtomTerm.__init__(self, atom1, atom2, atom3, atom4)
      self.oopbend_type = oopbend_type

   def write_info(self, parm, idx):
      # If unused, skip this parameter
      if idx == -1: return
      key = 'AMOEBA_OPBEND_ANGLE_LIST'
      _FourAtomTerm.write_info(self, parm, key, idx)
      parm.parm_data[key][5*idx+4] = self.oopbend_type.idx + 1

   def __eq__(self, other):
      if (self.atom1 is other.atom1 and self.atom2 is other.atom2 and
          self.atom3 is other.atom3 and self.atom4 is other.atom4):
         return self.trigang_type == other.trigang_type
      if (self.atom3 is other.atom1 and self.atom2 is other.atom2 and
          self.atom3 is other.atom1 and self.atom4 is other.atom4):
         return self.trigang_type == other.trigang_type
      return False

   def __contains__(self, thing):
      if isinstance(thing, Bond):
         return ((self.atom1 in thing and self.atom2 in thing) or
                 (self.atom2 in thing and self.atom3 in thing) or
                 (self.atom2 in thing and self.atom4 in thing))
      return super(OutOfPlaneBend, self).__contains__(thing)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class OutOfPlaneBendType(object):

   def __init__(self, k, idx):
      self.k = k
      self.idx = idx

   def write_info(self, parm):
      if self.idx == -1: return
      parm.parm_data['AMOEBA_OPBEND_ANGLE_FORCE_CONSTANT'][self.idx] = self.k

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class OutOfPlaneBendTypeList(_TypeList):

   def __init__(self, parm):
      _TypeList.__init__(self, parm)
      self.degree = parm.parm_data['AMOEBA_OPBEND_ANGLE_FTAB_DEGREE'][0]
      self.coeffs = parm.parm_data['AMOEBA_OPBEND_ANGLE_FTAB_COEFFS'][:]
      if len(self.coeffs) != self.degree + 1:
         raise AmoebaError('OOP-bend angle degree (%d) does not make sense with'
                          ' %d coefficients.' % (self.degree, len(self.coeffs)))
   
   def _make_array(self):
      klist = self.parm.parm_data['AMOEBA_OPBEND_ANGLE_FORCE_CONSTANT']
      list.__init__(self, [OutOfPlaneBendType(k, -1) for k in klist])

   def write_to_parm(self):
      self.parm.parm_data['AMOEBA_OPBEND_ANGLE_FTAB_DEGREE'][0] = self.degree
      self.parm.parm_data['AMOEBA_OPBEND_ANGLE_FTAB_COEFFS'] = self.coeffs[:]
      super(OutOfPlaneBendTypeList, self).write_to_parm()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Dihedral(_FourAtomTerm):
   """ Normal torsion term """

   def __init__(self, atom1, atom2, atom3, atom4, dihedral_type):
      _FourAtomTerm.__init__(self, atom1, atom2, atom3, atom4)
      self.dihedral_type = dihedral_type
      self.register()

   def write_info(self, parm, idx):
      # If unused, skip this parameter
      if idx == -1: return
      key = 'AMOEBA_TORSION_LIST'
      _FourAtomTerm.write_info(self, parm, key, idx)
      parm.parm_data[key][5*idx+4] = self.dihedral_type.idx + 1

   def __eq__(self, other):
      if (self.atom1 is other.atom1 and self.atom2 is other.atom2 and
          self.atom3 is other.atom3 and self.atom4 is other.atom4):
         return self.trigang_type == other.trigang_type
      if (self.atom3 is other.atom1 and self.atom2 is other.atom2 and
          self.atom3 is other.atom1 and self.atom4 is other.atom4):
         return self.trigang_type == other.trigang_type
      return False

   def register(self):
      self.atom1.dihedral_to(self.atom2)
      self.atom1.dihedral_to(self.atom3)
      self.atom1.dihedral_to(self.atom4)
      self.atom2.dihedral_to(self.atom3)
      self.atom2.dihedral_to(self.atom4)
      self.atom3.dihedral_to(self.atom4)

   def __contains__(self, thing):
      if isinstance(thing, Bond):
         return ((self.atom1 in thing and self.atom2 in thing) or
                 (self.atom2 in thing and self.atom3 in thing) or
                 (self.atom3 in thing and self.atom4 in thing))
      return super(Dihedral, self).__contains__(thing)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralType(object):

   def __init__(self, phi_k, per, phase, idx):
      self.phi_k = phi_k
      self.per = per
      self.phase = phase
      self.idx = idx

   def write_info(self, parm):
      if self.idx == -1: return
      parm.parm_data['AMOEBA_TORSION_FORCE_CONSTANT'][self.idx] = self.phi_k
      parm.parm_data['AMOEBA_TORSION_PERIODICITY'][self.idx] = self.per
      parm.parm_data['AMOEBA_TORSION_PHASE'][self.idx] = self.phase

   def __eq__(self, other):
      return (self.phi_k == other.phi_k and
              self.per == other.per and
              self.phase == other.phase)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralTypeList(_TypeList):

   def _make_array(self):
      klist = self.parm.parm_data['AMOEBA_TORSION_FORCE_CONSTANT']
      perlist = self.parm.parm_data['AMOEBA_TORSION_PERIODICITY']
      phaselist = self.parm.parm_data['AMOEBA_TORSION_PHASE']
      list.__init__(self, [DihedralType(k, per, phase, -1)
                           for k, per, phase in zip(klist, perlist, phaselist)])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PiTorsion(object):

   def __init__(self, atom1, atom2, atom3, atom4, atom5, atom6, pitor_type):
      self.atom1 = atom1
      self.atom2 = atom2
      self.atom3 = atom3
      self.atom4 = atom4
      self.atom5 = atom5
      self.atom6 = atom6
      self.pitor_type = pitor_type

   def write_info(self, parm, idx):
      if idx == -1: return
      key = 'AMOEBA_PI_TORSION_LIST'
      parm.parm_data[key][7*idx  ] = self.atom1.idx + 1
      parm.parm_data[key][7*idx+1] = self.atom2.idx + 1
      parm.parm_data[key][7*idx+2] = self.atom3.idx + 1
      parm.parm_data[key][7*idx+3] = self.atom4.idx + 1
      parm.parm_data[key][7*idx+4] = self.atom5.idx + 1
      parm.parm_data[key][7*idx+5] = self.atom6.idx + 1
      parm.parm_data[key][7*idx+6] = self.pitor_type.idx + 1

   def __contains__(self, thing):
      if isinstance(thing, Bond):
         return ((self.atom2 in thing and self.atom3 in thing) or
                 (self.atom3 in thing and self.atom4 in thing) or
                 (self.atom4 in thing and self.atom5 in thing) or
                 (self.atom4 in thing and self.atom6 in thing))
      return (thing is self.atom1 or thing is self.atom2 or thing is self.atom3
           or thing is self.atom4 or thing is self.atom5 or thing is self.atom6)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PiTorsionType(DihedralType):

   def write_info(self, parm):
      if self.idx == -1: return
      parm.parm_data['AMOEBA_PI_TORSION_FORCE_CONSTANT'][self.idx] = self.phi_k
      parm.parm_data['AMOEBA_PI_TORSION_PERIODICITY'][self.idx] = self.per
      parm.parm_data['AMOEBA_PI_TORSION_PHASE'][self.idx] = self.phase

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PiTorsionTypeList(_TypeList):

   def _make_array(self):
      klist = self.parm.parm_data['AMOEBA_PI_TORSION_FORCE_CONSTANT']
      perlist = self.parm.parm_data['AMOEBA_PI_TORSION_PERIODICITY']
      phaselist = self.parm.parm_data['AMOEBA_PI_TORSION_PHASE']
      list.__init__(self, [PiTorsionType(k, per, phase, -1)
                           for k, per, phase in zip(klist, perlist, phaselist)])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class StretchBend(object):

   def __init__(self, atom1, atom2, atom3, strbnd_type):
      self.atom1 = atom1
      self.atom2 = atom2
      self.atom3 = atom3
      self.strbnd_type = strbnd_type

   def write_info(self, parm, idx):
      if idx == -1: return
      key = 'AMOEBA_STRETCH_BEND_LIST'
      parm.parm_data[key][4*idx  ] = self.atom1.idx + 1
      parm.parm_data[key][4*idx+1] = self.atom2.idx + 1
      parm.parm_data[key][4*idx+2] = self.atom3.idx + 1
      parm.parm_data[key][4*idx+3] = self.strbnd_type.idx + 1

   def __contains__(self, thing):
      if isinstance(thing, Bond):
         return ((self.atom1 in thing and self.atom2 in thing) or
                 (self.atom2 in thing and self.atom3 in thing))
      return thing is self.atom1 or thing is self.atom2 or thing is self.atom3

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class StretchBendType(object):

   def __init__(self, k, theteq, bond1eq, bond2eq, idx):
      self.k = k
      self.theteq = theteq
      self.bond1eq = bond1eq
      self.bond2eq = bond2eq
      self.idx = idx

   def write_info(self, parm):
      if self.idx == -1: return
      parm.parm_data['AMOEBA_STRETCH_BEND_FORCE_CONSTANT'][self.idx] = self.k
      parm.parm_data['AMOEBA_STRETCH_BEND_ANGLE_EQUIL_VALUE'][self.idx] = \
            self.theteq
      parm.parm_data['AMOEBA_STRETCH_BEND_BOND1_EQUIL_VALUE'][self.idx] = \
            self.bond1eq
      parm.parm_data['AMOEBA_STRETCH_BEND_BOND2_EQUIL_VALUE'][self.idx] = \
            self.bond2eq

   def __eq__(self, other):
      return (self.k == other.k and self.theteq == other.theteq and
              self.bond1eq == other.bond1eq and self.bond2eq == other.bond2eq)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class StretchBendTypeList(_TypeList):

   def _make_array(self):
      klist = self.parm.parm_data['AMOEBA_STRETCH_BEND_FORCE_CONSTANT']
      eqlist = self.parm.parm_data['AMOEBA_STRETCH_BEND_ANGLE_EQUIL_VALUE']
      eq1list = self.parm.parm_data['AMOEBA_STRETCH_BEND_BOND1_EQUIL_VALUE']
      eq2list = self.parm.parm_data['AMOEBA_STRETCH_BEND_BOND2_EQUIL_VALUE']
      list.__init__(self, [StretchBendType(k, eq, eq1, eq2, -1)
                   for k, eq, eq1, eq2 in zip(klist, eqlist, eq1list, eq2list)])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TorsionTorsion(_Cmap):
   
   def __init__(self, *args, **kwargs):
      _Cmap.__init__(self, *args, **kwargs)
      self.register()

   def write_info(self, parm, idx):
      if idx == -1: return
      _Cmap.write_info(self, parm, 'AMOEBA_TORSION_TORSION_LIST', idx)

   def register(self):
      self.atom1.tortor_to(self.atom2)
      self.atom1.tortor_to(self.atom3)
      self.atom1.tortor_to(self.atom4)
      self.atom1.tortor_to(self.atom5)
      self.atom2.tortor_to(self.atom3)
      self.atom2.tortor_to(self.atom4)
      self.atom2.tortor_to(self.atom5)
      self.atom3.tortor_to(self.atom4)
      self.atom3.tortor_to(self.atom5)
      self.atom4.tortor_to(self.atom5)

   @property
   def tortor_type(self):
      return self.cmap_type

   @tortor_type.setter
   def tortor_type(self, thing):
      self.cmap_type = thing

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _TorTorTable(object):
   """
   A helper function for a coupled-torsion type: this is the interpolation
   table
   """

   def __init__(self, ang1, ang2, data):
      """
      Table for coupled-torsion interpolation.

      Parameters:
         ang1 (int) : Number of points in the first angle
         ang2 (int) : Number of points in the second angle
         data (list of floats) : Value of the potential grid at each point
                                 (ang2 varies fastest)
      """
      if len(data) != len(ang1) * len(ang2):
         raise AmoebaError('Coupled torsion parameter size mismatch. %dx%d grid'
                           ' expects %d elements (got %d)' % (len(ang1),
                           len(ang2), len(ang1)*len(ang2), len(data)))
      self.data = data
      self._indexes = dict()
      i = 0
      for a1 in ang1:
         for a2 in ang2:
            self._indexes[(a1, a2)] = i
            i += 1

   def __getitem__(self, idx, idx2=None):
      """
      Allow users to access table data from a tuple of 2 indexes or from two
      individual indexes
      """
      if idx2 is None:
         return self.data[self._indexes[idx]]
      return self.data[self._indexes[(idx, idx2)]]

   def __setitem__(self, idx, second, third=None):
      """ Sets a grid point """
      if third is not None:
         idx = self._indexes[(idx, second)]
         value = third
      else:
         idx = self._indexes[idx]
         value = second

      self.data[idx] = value

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TorsionTorsionType(object):
   """ Correction-map type """

   def __init__(self, dims, ang1, ang2, f, dfda1, dfda2, d2fda1da2, idx):
      """ Stores a coupled torsion with 4 interpolation tables passed in

      Parameters:
         dims (tuple of 2 ints): Table dimensions
         ang1 (list of floats): Angles in 1st dimension of interpolation table
         ang2 (list of floats): Angles in 2nd dimension of interpolation table
         f (list of floats): Interpolation table for the energy
         dfda1 (list of floats): Gradient of f w.r.t. angle 1
         dfda2 (list of floats): Gradient of f w.r.t. angle 2
         d2fda1da2 (list of floats): Second derivative of f w.r.t. both angles
         idx (int): Index of this type
      """

      self.dims = tuple(dims)
      if self.dims[0] != len(ang1):
         raise AmoebaError('Mismatched torsion-torsion dimensions [%d; got %d]'
               % (self.dims[0], len(ang1)))
      if self.dims[1] != len(ang2):
         raise AmoebaError('Mismatched torsion-torsion dimensions [%d; got %d]'
               % (self.dims[1], len(ang2)))

      self.ang1 = ang1
      self.ang2 = ang2
      self.f = _TorTorTable(ang1, ang2, f)
      self.dfda1 = _TorTorTable(ang1, ang2, dfda1)
      self.dfda2 = _TorTorTable(ang1, ang2, dfda2)
      self.d2fda1da2 = _TorTorTable(ang1, ang2, d2fda1da2)
      self.idx = idx

   def write_info(self, parm, after):
      if self.idx == -1: return
      # We need to get rid of the existing sections in the prmtop and overwrite
      # it with the data here
      i = self.idx + 1
      prefix = 'AMOEBA_TORSION_TORSION_TORTOR_TABLE_%02d'

      # Add these flags back with the correct data
      tblcmnt = ['dimension = (%d,%d)' % (self.dims)]
      parm.addFlag((prefix + '_DIMS') % i, '2I8', data=list(self.dims),
            comments=['dimension = (2)'], after=after)
      parm.addFlag((prefix + '_ANGLE1') % i, '5E16.8', data=self.ang1[:],
            comments=['dimension = (%d)' % self.dims[0]],
            after=(prefix + '_DIMS') % i)
      parm.addFlag((prefix + '_ANGLE2') % i, '5E16.8', data=self.ang2[:],
            comments=['dimension = (%d)' % self.dims[0]],
            after=(prefix + '_ANGLE1') % i)
      parm.addFlag((prefix + '_FUNC') % i, '5E16.8', data=self.f.data[:],
             comments=tblcmnt, after=(prefix + '_ANGLE2') % i)
      parm.addFlag((prefix + '_DFUNC_DANGLE1') % i, '5E16.8',
            data=self.dfda1.data[:], comments=tblcmnt,
            after=(prefix + '_FUNC') % i)
      parm.addFlag((prefix + '_DFUNC_DANGLE2') % i, '5E16.8',
            data=self.dfda2.data[:], comments=tblcmnt,
            after=(prefix + '_DFUNC_DANGLE1') % i)
      parm.addFlag((prefix + '_D2FUNC_DANGLE1_DANGLE2') % i, '5E16.8',
            data=self.d2fda1da2.data[:], comments=tblcmnt,
            after=(prefix + '_DFUNC_DANGLE2') % i)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TorsionTorsionTypeList(_TypeList):

   def _make_array(self):
      ntab = self.parm.parm_data['AMOEBA_TORSION_TORSION_NUM_PARAMS'][0]
      list.__init__(self)
      for i in range(ntab):
         prefix = 'AMOEBA_TORSION_TORSION_TORTOR_TABLE_%02d_' % (i + 1)
         dims = tuple(self.parm.parm_data[prefix + 'DIMS'])
         ang1 = self.parm.parm_data[prefix + 'ANGLE1']
         ang2 = self.parm.parm_data[prefix + 'ANGLE2']
         f = self.parm.parm_data[prefix + 'FUNC']
         dfda1 = self.parm.parm_data[prefix + 'DFUNC_DANGLE1']
         dfda2 = self.parm.parm_data[prefix + 'DFUNC_DANGLE2']
         d2fda1da2 = self.parm.parm_data[prefix + 'D2FUNC_DANGLE1_DANGLE2']
         self.append(TorsionTorsionType(dims, ang1, ang2, f, dfda1, dfda2,
                                        d2fda1da2, -1))

   def write_to_parm(self):
      after = 'AMOEBA_TORSION_TORSION_NUM_PARAMS'
      for i, typ in enumerate(self):
         typ.write_info(self.parm, after)
         after = ('AMOEBA_TORSION_TORSION_TORTOR_TABLE_%02d_D2FUNC_DANGLE1_'
                  'DANGLE2' % (i+1))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ChiralFrame(object):

   def __init__(self, atom1, atom2, chirality):
      self.atom1 = atom1
      self.atom2 = atom2
      self.chirality = chirality

   def write_info(self, parm, idx):
      if idx == -1: return
      key = 'AMOEBA_CHIRAL_FRAME_LIST'
      parm.parm_data[key][3*idx  ] = self.atom1.idx + 1
      parm.parm_data[key][3*idx+1] = self.atom2.idx + 1
      parm.parm_data[key][3*idx+2] = self.chirality # 1 or -1

   def __contains__(self, thing):
      return thing is self.atom1 or thing is self.atom2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class MultipoleFrame(object):
   
   def __init__(self, atom, frame_pt_num, vectail, vechead, nvec):
      self.atom = atom
      self.frame_pt_num = frame_pt_num
      self.vectail = vectail
      self.vechead = vechead
      self.nvec = nvec

   def write_info(self, parm, idx):
      if idx == -1: return
      key = 'AMOEBA_FRAME_DEF_LIST'
      parm.parm_data[key][5*idx  ] = self.atom.idx + 1
      parm.parm_data[key][5*idx+1] = self.frame_pt_num
      parm.parm_data[key][5*idx+2] = self.vectail
      parm.parm_data[key][5*idx+3] = self.vechead
      parm.parm_data[key][5*idx+4] = self.nvec

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ExclusionAssignment(object):
   
   def __init__(self, atom1, atom2, weight):
      self.atom1 = atom1
      self.atom2 = atom2
      self.weight = weight

   def write_info(self, parm, idx):
      if idx == -1: return
      data = parm.parm_data['AMOEBA_ADJUST_LIST']
      data[3*idx  ] = self.atom1.idx + 1
      data[3*idx+1] = self.atom2.idx + 1
      data[3*idx+2] = self.weight

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ExclusionAssignmentWeights(object):
   
   def __init__(self, parm):
      self.vdw = parm.parm_data['AMOEBA_ADJUST_VDW_WEIGHTS_LIST'][:]
      self.mpole = parm.parm_data['AMOEBA_ADJUST_MPOLE_WEIGHTS_LIST'][:]
      self.direct = parm.parm_data['AMOEBA_ADJUST_DIRECT_WEIGHTS_LIST'][:]
      self.polar = parm.parm_data['AMOEBA_ADJUST_POLAR_WEIGHTS_LIST'][:]
      self.mutual = parm.parm_data['AMOEBA_ADJUST_MUTUAL_WEIGHTS_LIST'][:]
      self.parm = parm

   def write_to_parm(self):
      data = self.parm.parm_data
      data['AMOEBA_ADJUST_VDW_WEIGHTS_LIST'] = self.vdw[:]
      data['AMOEBA_ADJUST_MPOLE_WEIGHTS_LIST'] = self.mpole[:]
      data['AMOEBA_ADJUST_DIRECT_WEIGHTS_LIST'] = self.direct[:]
      data['AMOEBA_ADJUST_POLAR_WEIGHTS_LIST'] = self.polar[:]
      data['AMOEBA_ADJUST_MUTUAL_WEIGHTS_LIST'] = self.mutual[:]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ResidueList = _ResidueList
