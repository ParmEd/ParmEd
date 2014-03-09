"""
This module contains an amber prmtop class that will read in all
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

from chemistry.amber._amberparm import AmberParm
from chemistry.amber.constants import NRES
from chemistry.amber.tinkertopology import (AtomList, ResidueList, TrackedList,
      Bond, BondTypeList, PiTorsionTypeList, PiTorsion, UreyBradleyTypeList,
      UreyBradley, AngleTypeList, Angle, TrigonalAngle, TrigonalAngleTypeList,
      OutOfPlaneBend, OutOfPlaneBendTypeList, Dihedral, DihedralTypeList,
      StretchBendTypeList, StretchBend, TorsionTorsion, TorsionTorsionTypeList,
      ChiralFrame, MultipoleFrame)
from chemistry.exceptions import (FormatError, AmberParmWarning)
from warnings import warn

class AmoebaParm(AmberParm):
   """
   Tinker Topology (parm7 format) class. Gives low, and some high, level access
   to topology data.
   """

   solvent_residues = ['WAT', 'HOH']

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def initialize_topology(self, rst7_name=None):
      """
      Initializes topology data structures, like the list of atoms, bonds, etc.,
      after the topology file has been read. The following methods are called:
      """
      try:
         if self.parm_data['AMOEBA_FORCEFIELD'][0] != 1:
            raise FormatError('Bad AMOEBA-format topology')
      except KeyError:
         raise FormatError('Bad AMOEBA-format topology')

      # Fill the residue container array
      self._fill_res_container()

      # We need to handle RESIDUE_ICODE properly since it may have picked up
      # some extra values
      if 'RESIDUE_ICODE' in self.flag_list:
         self.parm_data['RESIDUE_ICODE'] = \
               self.parm_data['RESIDUE_ICODE'][:self.parm_data['POINTERS'][NRES]]

      # instance variables other than those in AmberFormat
      self.pointers = {}       # list of all the pointers in the prmtop
      self.LJ_types = {}       # dict pairing atom name with its LJ atom type #
      self.LJ_radius = []      # ordered array of L-J radii in Ang -- indices
                               # are elements in LJ_types-1
      self.LJ_depth = []       # similarly ordered array of L-J depths

      # If we were given a prmtop, read it in
      if self.valid:
         try:
            self.LoadPointers()
            self.valid = True
         except KeyError:
            warn('POINTERS flag not found! Likely a bad AMBER topology file.',
                 AmberParmWarning)
            self.valid = False
         except IndexError:
            if (len(self.parm_data['POINTERS'])) < 30:
               warn('Fewer integers in POINTERS section than expected! Likely '
                    'a bad AMBER topology file.', AmberParmWarning)
               self.valid = False

         # Load the structure arrays
#        try:
         self._load_structure()
#        except (KeyError, IndexError, AttributeError):
#           warn('Error loading molecule topology. Cannot use delete_mask',
#                AmberParmWarning)
#           self.valid = False
      
      if rst7_name is not None:
         self.LoadRst7(rst7_name)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def _load_structure(self):
      """ Responsible for setting up the parameter and parameter type arrays """

      ##### First create our atoms #####
      self.atom_list = AtomList(self)
      ##### Next, load our residues #####
      self.residue_list = ResidueList(self)
      self.atom_list.changed = False
      self.residue_list.changed = False
      for i in range(self.ptr('natom')):
         residx = self.residue_container[i] - 1
         self.atom_list[i].residue = self.residue_list[residx]
      ##### Next create our list of bonds #####
      self.bond_type_list = BondTypeList(self)
      self.bond_list = TrackedList()
      # Add all of our bonds
      for i in range(self.parm_data['AMOEBA_REGULAR_BOND_NUM_LIST'][0]):
         id1 = self.parm_data['AMOEBA_REGULAR_BOND_LIST'][3*i  ] - 1
         id2 = self.parm_data['AMOEBA_REGULAR_BOND_LIST'][3*i+1] - 1
         typ = self.parm_data['AMOEBA_REGULAR_BOND_LIST'][3*i+2] - 1
         self.bond_list.append(
               Bond(self.atom_list[id1], self.atom_list[id2],
                    self.bond_type_list[typ])
         )
      self.bond_list.changed = False
      self.bond_type_list.changed = False
      ##### Next create our list of urey-bradleys #####
      self.urey_bradley_type_list = UreyBradleyTypeList(self)
      self.urey_bradley_list = TrackedList()
      # Add all of our urey-bradley terms
      for i in range(self.parm_data['AMOEBA_UREY_BRADLEY_BOND_NUM_LIST'][0]):
         id1 = self.parm_data['AMOEBA_UREY_BRADLEY_BOND_LIST'][3*i  ] - 1
         id2 = self.parm_data['AMOEBA_UREY_BRADLEY_BOND_LIST'][3*i+1] - 1
         typ = self.parm_data['AMOEBA_UREY_BRADLEY_BOND_LIST'][3*i+2] - 1
         self.urey_bradley_list.append(
               UreyBradley(self.atom_list[id1], self.atom_list[id2],
                           self.urey_bradley_type_list[typ])
         )
      self.urey_bradley_type_list.changed = False
      self.urey_bradley_list.changed = False
      ##### Next create our list of angles #####
      self.angle_type_list = AngleTypeList(self)
      self.angle_list = TrackedList()
      # Add all of our angles
      for i in range(self.parm_data['AMOEBA_REGULAR_ANGLE_NUM_LIST'][0]):
         id1 = self.parm_data['AMOEBA_REGULAR_ANGLE_LIST'][4*i  ] - 1
         id2 = self.parm_data['AMOEBA_REGULAR_ANGLE_LIST'][4*i+1] - 1
         id3 = self.parm_data['AMOEBA_REGULAR_ANGLE_LIST'][4*i+2] - 1
         typ = self.parm_data['AMOEBA_REGULAR_ANGLE_LIST'][4*i+3] - 1
         self.angle_list.append(
               Angle(self.atom_list[id1], self.atom_list[id2],
                     self.atom_list[id3], self.angle_type_list[typ])
         )
      self.angle_type_list.changed = False
      self.angle_list.changed = False
      ##### Next create our list of trigonal angles (in-plane)
      self.trigonal_angle_type_list = TrigonalAngleTypeList(self)
      self.trigonal_angle_list = TrackedList()
      # Add all trigonal angles
      for i in range(self.parm_data['AMOEBA_TRIGONAL_ANGLE_NUM_LIST'][0]):
         id1 = self.parm_data['AMOEBA_TRIGONAL_ANGLE_LIST'][5*i  ] - 1
         id2 = self.parm_data['AMOEBA_TRIGONAL_ANGLE_LIST'][5*i+1] - 1
         id3 = self.parm_data['AMOEBA_TRIGONAL_ANGLE_LIST'][5*i+2] - 1
         id4 = self.parm_data['AMOEBA_TRIGONAL_ANGLE_LIST'][5*i+3] - 1
         typ = self.parm_data['AMOEBA_TRIGONAL_ANGLE_LIST'][5*i+4] - 1
         self.trigonal_angle_list.append(
               TrigonalAngle(self.atom_list[id1], self.atom_list[id2],
                             self.atom_list[id3], self.atom_list[id4],
                             self.trigonal_angle_type_list[typ])
         )
      self.trigonal_angle_type_list.changed = False
      self.trigonal_angle_list.changed = False
      ##### Next create our list of out-of-plane bending terms #####
      self.oopbend_type_list = OutOfPlaneBendTypeList(self)
      self.oopbend_list = TrackedList()
      # Add the out-of-plane bending terms
      for i in range(self.parm_data['AMOEBA_OPBEND_ANGLE_NUM_LIST'][0]):
         id1 = self.parm_data['AMOEBA_OPBEND_ANGLE_LIST'][5*i  ] - 1
         id2 = self.parm_data['AMOEBA_OPBEND_ANGLE_LIST'][5*i+1] - 1
         id3 = self.parm_data['AMOEBA_OPBEND_ANGLE_LIST'][5*i+2] - 1
         id4 = self.parm_data['AMOEBA_OPBEND_ANGLE_LIST'][5*i+3] - 1
         typ = self.parm_data['AMOEBA_OPBEND_ANGLE_LIST'][5*i+4] - 1
         self.oopbend_list.append(
               OutOfPlaneBend(self.atom_list[id1], self.atom_list[id2],
                             self.atom_list[id3], self.atom_list[id4],
                             self.oopbend_type_list[typ])
         )
      self.oopbend_type_list.changed = False
      self.oopbend_list.changed = False
      ##### Next create our list of normal dihedrals #####
      self.dihedral_type_list = DihedralTypeList(self)
      self.dihedral_list = TrackedList()
      # Add the dihedrals
      for i in range(self.parm_data['AMOEBA_TORSION_NUM_LIST'][0]):
         id1 = self.parm_data['AMOEBA_TORSION_LIST'][5*i  ] - 1
         id2 = self.parm_data['AMOEBA_TORSION_LIST'][5*i+1] - 1
         id3 = self.parm_data['AMOEBA_TORSION_LIST'][5*i+2] - 1
         id4 = self.parm_data['AMOEBA_TORSION_LIST'][5*i+3] - 1
         typ = self.parm_data['AMOEBA_TORSION_LIST'][5*i+4] - 1
         self.dihedral_list.append(
               Dihedral(self.atom_list[id1], self.atom_list[id2],
                        self.atom_list[id3], self.atom_list[id4],
                        self.dihedral_type_list[typ])
         )
      self.dihedral_type_list.changed = False
      self.dihedral_list.changed = False
      ##### Next create our list of pi-torsions #####
      self.pitorsion_type_list = PiTorsionTypeList(self)
      self.pitorsion_list = TrackedList()
      # Add the pi-torsions
      for i in range(self.parm_data['AMOEBA_PI_TORSION_NUM_LIST'][0]):
         id1 = self.parm_data['AMOEBA_PI_TORSION_LIST'][7*i  ] - 1
         id2 = self.parm_data['AMOEBA_PI_TORSION_LIST'][7*i+1] - 1
         id3 = self.parm_data['AMOEBA_PI_TORSION_LIST'][7*i+2] - 1
         id4 = self.parm_data['AMOEBA_PI_TORSION_LIST'][7*i+3] - 1
         id5 = self.parm_data['AMOEBA_PI_TORSION_LIST'][7*i+4] - 1
         id6 = self.parm_data['AMOEBA_PI_TORSION_LIST'][7*i+5] - 1
         typ = self.parm_data['AMOEBA_PI_TORSION_LIST'][7*i+6] - 1
         self.pitorsion_list.append(
               PiTorsion(self.atom_list[id1], self.atom_list[id2],
                         self.atom_list[id3], self.atom_list[id4],
                         self.atom_list[id5], self.atom_list[id6],
                         self.pitorsion_type_list[typ])
         )
      self.pitorsion_type_list.changed = False
      self.pitorsion_list.changed = False
      ##### Next create stretch-bend terms #####
      self.stretch_bend_type_list = StretchBendTypeList(self)
      self.stretch_bend_list = TrackedList()
      # Add the stretch-bends
      for i in range(self.parm_data['AMOEBA_STRETCH_BEND_NUM_LIST'][0]):
         id1 = self.parm_data['AMOEBA_STRETCH_BEND_LIST'][4*i  ] - 1
         id2 = self.parm_data['AMOEBA_STRETCH_BEND_LIST'][4*i+1] - 1
         id3 = self.parm_data['AMOEBA_STRETCH_BEND_LIST'][4*i+2] - 1
         typ = self.parm_data['AMOEBA_STRETCH_BEND_LIST'][4*i+3] - 1
         self.stretch_bend_list.append(
               StretchBend(self.atom_list[id1], self.atom_list[id2],
                           self.atom_list[id3],
                           self.stretch_bend_type_list[typ])
         )
      self.stretch_bend_type_list.changed = False
      self.stretch_bend_list.changed = False
      ##### Next create the torsion-torsion parameters #####
      self.torsion_torsion_type_list = TorsionTorsionTypeList(self)
      self.torsion_torsion_list = TrackedList()
      # Add the torsion-torsions
      for i in range(self.parm_data['AMOEBA_TORSION_TORSION_NUM_LIST'][0]):
         id1 = self.parm_data['AMOEBA_TORSION_TORSION_LIST'][6*i  ] - 1
         id2 = self.parm_data['AMOEBA_TORSION_TORSION_LIST'][6*i+1] - 1
         id3 = self.parm_data['AMOEBA_TORSION_TORSION_LIST'][6*i+2] - 1
         id4 = self.parm_data['AMOEBA_TORSION_TORSION_LIST'][6*i+3] - 1
         id5 = self.parm_data['AMOEBA_TORSION_TORSION_LIST'][6*i+4] - 1
         typ = self.parm_data['AMOEBA_TORSION_TORSION_LIST'][6*i+5] - 1
         self.torsion_torsion_list.append(
               TorsionTorsion(self.atom_list[id1], self.atom_list[id2],
                              self.atom_list[id3], self.atom_list[id4],
                              self.atom_list[id5],
                              self.torsion_torsion_type_list[typ])
         )
      self.torsion_torsion_type_list.changed = False
      self.torsion_torsion_list.changed = False
      ##### Next create the chiral frame list #####
      self.chiral_frame_list = TrackedList()
      for i in range(self.parm_data['AMOEBA_CHIRAL_FRAME_NUM_LIST'][0]):
         id1 = self.parm_data['AMOEBA_CHIRAL_FRAME_LIST'][3*i  ] - 1
         id2 = self.parm_data['AMOEBA_CHIRAL_FRAME_LIST'][3*i+1] - 1
         chi = self.parm_data['AMOEBA_CHIRAL_FRAME_LIST'][3*i+2]
         self.chiral_frame_list.append(
               ChiralFrame(self.atom_list[id1], self.atom_list[id2], chi)
         )
      self.chiral_frame_list.changed = False
      ##### Finally create the multipole frame list #####
      self.multipole_frame_list = TrackedList()
      for i in range(self.parm_data['AMOEBA_FRAME_DEF_NUM_LIST'][0]):
         id1 = self.parm_data['AMOEBA_FRAME_DEF_LIST'][5*i  ] - 1
         fpn = self.parm_data['AMOEBA_FRAME_DEF_LIST'][5*i+1]
         vct = self.parm_data['AMOEBA_FRAME_DEF_LIST'][5*i+2]
         vch = self.parm_data['AMOEBA_FRAME_DEF_LIST'][5*i+3]
         nvc = self.parm_data['AMOEBA_FRAME_DEF_LIST'][5*i+4]
         self.multipole_frame_list.append(
               MultipoleFrame(self.atom_list[id1], fpn, vct, vch, nvc)
         )
      self.multipole_frame_list.changed = False
      ##### Finally we can determine the polar group exclusions #####
      for atom in self.atom_list: atom.determine_polar_partners()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def remake_parm(self):
      """ Recomputes the topology file parameters """
      raise NotImplemented('Amoeba topologies cannot currently be recomputed')

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   @property
   def chamber(self):
      return False
   
   @property
   def amoeba(self):
      return True
