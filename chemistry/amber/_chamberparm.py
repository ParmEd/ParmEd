"""
This module contains a chamber prmtop class that will read in all
parameters and allow users to manipulate that data and write a new
prmtop object.

Copyright (C) 2010 - 2014  Jason Swails

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""
from __future__ import division

from chemistry.amber.topologyobjects import (TrackedList, UreyBradley, Improper,
            Cmap, UreyBradleyTypeList, ImproperTypeList, CmapTypeList)
from chemistry.amber._amberparm import AmberParm, _zeros

class ChamberParm(AmberParm):
   """
   Chamber Topology (parm7 format) class. Gives low, and some high, level access
   to topology data.
   """

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def initialize_topology(self, rst7_name=None):
      """
      Initializes topology data structures, like the list of atoms, bonds, etc.,
      after the topology file has been read. The following methods are called:
      """

      self.LJ_14_radius = []   # Radii for 1-4 atom pairs
      self.LJ_14_depth = []    # Depth for 1-4 atom pairs

      AmberParm.initialize_topology(self, rst7_name)
      
      # We now have the following instance arrays: All arrays are dynamic such
      # that removing an item propagates the indices if applicable. bond has
      # angle/dihed analogs. All non-dynamic lists have a check on any
      # modification function to track if they have been changed or not so we
      # know whether we have to reload the data before writing.
      #
      # atom_list          a dynamic list of all Atom objects
      # residue_list       a dynamic list of all Residue objects
      # bond_type_list     a dynamic list of all BondType objects
      # bonds_inc_h        list of all bonds including hydrogen
      # bonds_without_h    list of all bonds without hydrogen
      # angle_type_list
      # angles_inc_h
      # angles_without_h
      # dihedral_type_list
      # dihedrals_inc_h
      # dihedrals_without_h
      # urey_bradley
      # urey_bradley_type_list
      # improper
      # improper_type_list
      #
      # MAYBE the following (if CMAPs are present)
      # cmap
      # cmap_type_list

      if rst7_name is not None:
         self.LoadRst7(rst7_name)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def __copy__(self):
      """ Needs to copy a few additional data structures """
      other = AmberParm.__copy__(self)
      if other.valid:
         other.LJ_14_radius = self.LJ_14_radius[:]
         other.LJ_14_depth = self.LJ_14_depth[:]
      else:
         other.LJ_14_radius = []
         other.LJ_14_depth = []
      return other

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   def LoadPointers(self):
      """
      Loads the data in POINTERS section into a pointers dictionary with each
      key being the pointer name according to http://ambermd.org/formats.html
      """
      AmberParm.LoadPointers(self)
      # Other pointers
      self.pointers['NUB'] = self.parm_data['CHARMM_UREY_BRADLEY_COUNT'][0]
      self.pointers['NUBTYPES'] = self.parm_data['CHARMM_UREY_BRADLEY_COUNT'][1]
      self.pointers['NIMPHI'] = self.parm_data['CHARMM_NUM_IMPROPERS'][0]
      self.pointers['NIMPRTYPES']=self.parm_data['CHARMM_NUM_IMPR_TYPES'][0]
      # If CMAP is not present, don't load the pointers
      try:
         self.pointers['CMAP'] = self.parm_data['CHARMM_CMAP_COUNT'][0]
         self.pointers['CMAP_TYPES'] = self.parm_data['CHARMM_CMAP_COUNT'][1]
      except KeyError:
         # CMAP does not exist in this parm
         pass

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def _load_structure(self):
      """ 
      Loads all of the topology instance variables. This is necessary if we
      actually want to modify the topological layout of our system
      (like deleting atoms)
      """
      AmberParm._load_structure(self)
      # The AmberParm _load_structure routine loads everything except the
      # Urey-Bradley, Improper, and CMAP lists. The latter are only present in
      # chamber prmtops created with a force field that has correction map
      # potentials, so support topologies that do not have this extra potential
      # as well.

      # Do the Urey-Bradley terms now
      self.urey_bradley = TrackedList()
      self.urey_bradley_type_list = UreyBradleyTypeList(self)
      for i in range(self.ptr('NUB')):
         a1 = self.parm_data['CHARMM_UREY_BRADLEY'][3*i  ] - 1
         a2 = self.parm_data['CHARMM_UREY_BRADLEY'][3*i+1] - 1
         ty = self.parm_data['CHARMM_UREY_BRADLEY'][3*i+2] - 1
         self.urey_bradley.append(
            UreyBradley(self.atom_list[a1], self.atom_list[a2],
                        self.urey_bradley_type_list[ty])
         )

      # Now do the improper torsion terms
      self.improper = TrackedList()
      self.improper_type_list = ImproperTypeList(self)
      for i in range(self.ptr('NIMPHI')):
         a1 = self.parm_data['CHARMM_IMPROPERS'][5*i  ] - 1
         a2 = self.parm_data['CHARMM_IMPROPERS'][5*i+1] - 1
         a3 = self.parm_data['CHARMM_IMPROPERS'][5*i+2] - 1
         a4 = self.parm_data['CHARMM_IMPROPERS'][5*i+3] - 1
         ty = self.parm_data['CHARMM_IMPROPERS'][5*i+4] - 1
         self.improper.append(
               Improper(self.atom_list[a1], self.atom_list[a2],
                        self.atom_list[a3], self.atom_list[a4],
                        self.improper_type_list[ty])
         )

      # Mark all new lists unchanged
      self.urey_bradley.changed = False
      self.urey_bradley_type_list.changed = False
      self.improper.changed = False
      self.improper_type_list.changed = False

      # Now if we have CMAP, do those
      if 'CHARMM_CMAP_COUNT' in self.flag_list:
         self.cmap = TrackedList()
         self.cmap_type_list = CmapTypeList(self)
         for i in range(self.ptr('CMAP')):
            a1 = self.parm_data['CHARMM_CMAP_INDEX'][6*i  ] - 1
            a2 = self.parm_data['CHARMM_CMAP_INDEX'][6*i+1] - 1
            a3 = self.parm_data['CHARMM_CMAP_INDEX'][6*i+2] - 1
            a4 = self.parm_data['CHARMM_CMAP_INDEX'][6*i+3] - 1
            a5 = self.parm_data['CHARMM_CMAP_INDEX'][6*i+4] - 1
            ty = self.parm_data['CHARMM_CMAP_INDEX'][6*i+5] - 1
            self.cmap.append(
                  Cmap(self.atom_list[a1], self.atom_list[a2],
                       self.atom_list[a3], self.atom_list[a4],
                       self.atom_list[a5], self.cmap_type_list[ty])
            )
         # Mark the cmap lists unchanged
         self.cmap.changed = False
         self.cmap_type_list.changed = False

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def __str__(self):
      """ Returns the name of the topology file as its string representation """
      return self.prm_name

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def ptr(self,pointer):
      """
      Returns the value of the given pointer, and converts to upper-case so it's
      case-insensitive. A non-existent pointer meets with a KeyError
      """
      return self.pointers[pointer.upper()]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   flush_data_changes = _load_structure

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def remake_parm(self):
      """
      Re-fills the topology file arrays if we have changed the underlying
      structure
      """
      # Handle all parts that are common with Amber via AmberParm
      AmberParm.remake_parm(self)

      # Now rebuild the remaining parts -- Urey-Bradley, Impropers, and CMAP

      # Urey-Bradley
      ub_num = ub_type_num = 0
      self.urey_bradley_type_list.deindex()
      self.parm_data['CHARMM_UREY_BRADLEY'] = _zeros(len(self.urey_bradley)*3)
      for i, ub in enumerate(self.urey_bradley):
         if -1 in (ub.atom1.idx, ub.atom2.idx):
            continue
         if ub.ub_type.idx == -1:
            ub.ub_type.idx = ub_type_num
            ub_type_num += 1
         ub.write_info(self, 'CHARMM_UREY_BRADLEY', ub_num)
         ub_num += 1
      # Truncate our list to only include those Urey-Bradleys that remain
      self.parm_data['CHARMM_UREY_BRADLEY_COUNT'] = [ub_num, ub_type_num]
      self._truncate_array('CHARMM_UREY_BRADLEY', 3*ub_num)
      # type parameters
      for key in ('CHARMM_UREY_BRADLEY_FORCE_CONSTANT',
                  'CHARMM_UREY_BRADLEY_EQUIL_VALUE'):
         self.parm_data[key] = _zeros(ub_type_num)
      self.urey_bradley_type_list.write_to_parm()

      # Impropers
      imp_num = imp_type_num = 0
      self.improper_type_list.deindex()
      self.parm_data['CHARMM_IMPROPERS'] = _zeros(len(self.improper) * 5)
      for i, imp in enumerate(self.improper):
         if -1 in (imp.atom1.idx, imp.atom2.idx, imp.atom3.idx, imp.atom4.idx):
            continue
         if imp.improp_type.idx == -1:
            imp.improp_type.idx = imp_type_num
            imp_type_num += 1
         imp.write_info(self, 'CHARMM_IMPROPERS', imp_num)
         imp_num += 1
      # Truncate our list to only include those impropers that remain
      self.parm_data['CHARMM_NUM_IMPROPERS'] = [imp_num]
      self._truncate_array('CHARMM_IMPROPERS', 5*imp_num)
      # type parameters
      self.parm_data['CHARMM_NUM_IMPR_TYPES'] = [imp_type_num]
      for key in ('CHARMM_IMPROPER_FORCE_CONSTANT', 'CHARMM_IMPROPER_PHASE'):
         self.parm_data[key] = _zeros(imp_type_num)
      self.improper_type_list.write_to_parm()

      # These arrays should no longer appear changed
      self.urey_bradley.changed = False
      self.urey_bradley_type_list.changed = False
      self.improper.changed = False
      self.improper_type_list.changed = False

      # If we have no CMAP, we are done. Otherwise, press on
      if not 'CHARMM_CMAP_COUNT' in self.flag_list:
         self.LoadPointers() # update CHARMM pointers
         return

      # If we are here, then we have CMAP terms to do
      cmap_num = 0
      cmap_types = []
      self.parm_data['CHARMM_CMAP_INDEX'] = _zeros(len(self.cmap)*6)
      for i, cm in enumerate(self.cmap):
         if -1 in (cm.atom1.idx, cm.atom2.idx, cm.atom3.idx, cm.atom4.idx,
                   cm.atom5.idx):
            continue
         if cm.cmap_type.idx == -1:
            cm.cmap_type.idx = len(cmap_types)
            cmap_types.append(cm.cmap_type)
         cm.write_info(self, 'CHARMM_CMAP_INDEX', cmap_num)
         cmap_num += 1
      if cmap_num == 0:
         # We have deleted all cmaps. Get rid of them from the parm file and
         # bail out. This is probably pretty unlikely, though...
         self.deleteFlag('CHARMM_CMAP_COUNT')
         self.deleteFlag('CHARMM_CMAP_RESOLUTION')
         for flag in self.flag_list:
            if flag.startswith('CHARMM_CMAP_PARAMETER'):
               self.deleteFlag(flag)
         self.LoadPointers() # update CHARMM pointers
         return
      # Truncate our list to only include those cmaps that remain
      self._truncate_array('CHARMM_CMAP_INDEX', 6*cmap_num)
      self.parm_data['CHARMM_CMAP_COUNT'] = [cmap_num, len(cmap_types)]

      # Now comes the tricky part. We need to delete all of the
      # CHARMM_CMAP_PARAMETER_XX sections and then recreate them with the
      # correct size and comments. The comments have been stored in the cmap
      # types themselves to prevent them from being lost. We will also assume
      # that the Fortran format we're going to use is the same for all CMAP
      # types, so just pull it from CHARMM_CMAP_PARAMETER_01 (or fall back to
      # 8(F9.5))
      try:
         fmt = str(self.formats['CHARMM_CMAP_PARAMETER_01'])
      except KeyError:
         fmt = '8(F9.5)'
      for flag in self.flag_list:
         if flag.startswith('CHARMM_CMAP_PARAMETER'):
            self.deleteFlag(flag)
      # Now the parameters are gone, so go through and add those flags back
      # using the cmap_types we tagged along earlier.
      for i, ct in enumerate(cmap_types):
         self.addFlag('CHARMM_CMAP_PARAMETER_%02d' % (i+1), fmt, data=ct.grid,
                      comments=ct.comments)
      self.LoadPointers() # update CHARMM pointers

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   def _topology_changed(self):
      """ 
      Determines if any of the topological arrays have changed since the
      last upload
      """
      tc = (AmberParm._topology_changed(self) or self.urey_bradley.changed or
            self.urey_bradley_type_list.changed or self.improper.changed or
            self.improper_type_list.changed)
      if 'CHARMM_CMAP_COUNT' in self.flag_list:
         return tc or self.cmap.changed or self.cmap_type_list.changed
      return tc

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def writeOFF(self, off_file='off.lib'):
      """ Writes an OFF file from all of the residues found in a prmtop """
      raise NotImplemented('Cannot write OFF file from chamber topology')

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def fill_LJ(self):
      """
      Fills the LJ_radius, LJ_depth arrays and LJ_types dictionary with data
      from LENNARD_JONES_ACOEF and LENNARD_JONES_BCOEF sections of the prmtop
      files, by undoing the canonical combining rules.
      """
      AmberParm.fill_LJ(self)
      self.LJ_14_radius = []  # empty LJ_radii so it can be re-filled
      self.LJ_14_depth = []   # empty LJ_depths so it can be re-filled
      one_sixth = 1.0 / 6.0 # we need to raise some numbers to the 1/6th power

      for i in range(self.pointers["NTYPES"]):
         lj_index = self.parm_data["NONBONDED_PARM_INDEX"][
                     self.pointers["NTYPES"] * i + i] - 1
         if self.parm_data["LENNARD_JONES_14_ACOEF"][lj_index] < 1.0e-6:
            self.LJ_14_radius.append(0)
            self.LJ_14_depth.append(0)
         else:
            factor = (2 * self.parm_data["LENNARD_JONES_14_ACOEF"][lj_index] /
                          self.parm_data["LENNARD_JONES_14_BCOEF"][lj_index])
            self.LJ_14_radius.append(pow(factor, one_sixth) * 0.5)
            self.LJ_14_depth.append(
                self.parm_data["LENNARD_JONES_14_BCOEF"][lj_index] / 2 / factor)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def recalculate_LJ(self):
      """
      Takes the values of the LJ_radius and LJ_depth arrays and recalculates the
      LENNARD_JONES_A/BCOEF topology sections from the canonical combining
      rules.
      """
      AmberParm.recalculate_LJ(self)
      for i in range(self.pointers["NTYPES"]):
         for j in range(i,self.pointers["NTYPES"]):
            index = self.parm_data['NONBONDED_PARM_INDEX'][
                                                   self.ptr('ntypes')*i + j
                                                          ] - 1
            rij = self.combine_rmin(self.LJ_14_radius[i],self.LJ_14_radius[j])
            wdij = self.combine_epsilon(self.LJ_14_depth[i],self.LJ_14_depth[j])
            self.parm_data["LENNARD_JONES_14_ACOEF"][index] = wdij * rij ** 12
            self.parm_data["LENNARD_JONES_14_BCOEF"][index] = 2 * wdij * rij**6

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   # For compatibility with the old AmberParm class that used to instantiate
   # both Amber and Chamber topologies, the `chamber' attribute indicates
   # whether or not a topology file is a chamber topology. For ChamberParm, the
   # answer is always 'yes'
   @property
   def chamber(self):
      return True
   
   @property
   def amoeba(self):
      return False

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def ToMolecule(self):
      """ Translates an amber system into a molecule format """
      raise NotImplemented('ToMolecule disabled for ChamberParm')
