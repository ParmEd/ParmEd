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

from chemistry.amber._amberparm import AmberParm, _zeros
from chemistry.amber.constants import NRES, NNB, NMXRS
from chemistry.amber.tinkertopology import (AtomList, ResidueList, TrackedList,
        Bond, BondTypeList, PiTorsionTypeList, PiTorsion, UreyBradleyTypeList,
        UreyBradley, AngleTypeList, Angle, TrigonalAngle, TrigonalAngleTypeList,
        OutOfPlaneBend, OutOfPlaneBendTypeList, Dihedral, DihedralTypeList,
        StretchBendTypeList, StretchBend, TorsionTorsion, ChiralFrame,
        TorsionTorsionTypeList, MultipoleFrame, ExclusionAssignment,
        ExclusionAssignmentWeights)
from chemistry.exceptions import (FormatError, AmberParmWarning)
from warnings import warn

class AmoebaParm(AmberParm):
    """
    Tinker Topology (parm7 format) class. Gives low, and some high, level
    access to topology data.
    """

    solvent_residues = ['WAT', 'HOH']

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def initialize_topology(self, rst7_name=None):
        """
        Initializes topology data structures, like the list of atoms, bonds,
        etc., after the topology file has been read. The following methods are
        called:
        """
        try:
            if self.parm_data['AMOEBA_FORCEFIELD'][0] != 1:
                raise FormatError('Bad AMOEBA-format topology')
        except KeyError:
            raise FormatError('Bad AMOEBA-format topology')

        # We need to handle RESIDUE_ICODE properly since it may have picked up
        # some extra values
        if 'RESIDUE_ICODE' in self.flag_list:
            self._truncate_array('RESIDUE_ICODE',
                                 self.parm_data['POINTERS'][NRES])

        # instance variables other than those in AmberFormat
        self.pointers = {}  # list of all the pointers in the prmtop
        self.LJ_types = {}  # dict pairing atom name with its LJ atom type #
        self.LJ_radius = [] # ordered array of L-J radii in Ang -- indices
                            # are elements in LJ_types-1
        self.LJ_depth = []  # similarly ordered array of L-J depths

        # If we were given a prmtop, read it in
        if self.valid:
            try:
                self.LoadPointers()
                self.valid = True
            except KeyError:
                warn('POINTERS flag not found! Likely a bad AMBER '
                     'topology file.', AmberParmWarning)
                self.valid = False
            except IndexError:
                if (len(self.parm_data['POINTERS'])) < 30:
                    warn('Fewer integers in POINTERS section than expected! '
                         'Likely a bad AMBER topology file.', AmberParmWarning)
                    self.valid = False

            # Load the structure arrays
            self._load_structure()
      
        if rst7_name is not None:
            self.LoadRst7(rst7_name)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def _load_structure(self):
        " Responsible for setting up the parameter and parameter type arrays "

        ##### First create our atoms #####
        self.atom_list = AtomList(self)
        ##### Next, load our residues #####
        self.residue_list = ResidueList(self)
        self.atom_list.changed = False
        self.residue_list.changed = False
        ##### Next create our list of bonds #####
        self.bond_type_list = BondTypeList(self)
        self.bond_list = TrackedList()
        # Add all of our bonds
        try:
            for i in range(self.parm_data['AMOEBA_REGULAR_BOND_NUM_LIST'][0]):
                id1 = self.parm_data['AMOEBA_REGULAR_BOND_LIST'][3*i  ] - 1
                id2 = self.parm_data['AMOEBA_REGULAR_BOND_LIST'][3*i+1] - 1
                typ = self.parm_data['AMOEBA_REGULAR_BOND_LIST'][3*i+2] - 1
                self.bond_list.append(
                        Bond(self.atom_list[id1], self.atom_list[id2],
                             self.bond_type_list[typ])
                )
        except KeyError:
            pass
        self.bond_list.changed = False
        self.bond_type_list.changed = False
        ##### Next create our list of urey-bradleys #####
        self.urey_bradley_type_list = UreyBradleyTypeList(self)
        self.urey_bradley_list = TrackedList()
        # Add all of our urey-bradley terms
        try:
            nub = self.parm_data['AMOEBA_UREY_BRADLEY_BOND_NUM_LIST'][0]
            for i in range(nub):
                id1 = self.parm_data['AMOEBA_UREY_BRADLEY_BOND_LIST'][3*i  ] - 1
                id2 = self.parm_data['AMOEBA_UREY_BRADLEY_BOND_LIST'][3*i+1] - 1
                typ = self.parm_data['AMOEBA_UREY_BRADLEY_BOND_LIST'][3*i+2] - 1
                self.urey_bradley_list.append(
                        UreyBradley(self.atom_list[id1], self.atom_list[id2],
                                    self.urey_bradley_type_list[typ])
                )
        except KeyError:
            pass
        self.urey_bradley_type_list.changed = False
        self.urey_bradley_list.changed = False
        ##### Next create our list of angles #####
        self.angle_type_list = AngleTypeList(self)
        self.angle_list = TrackedList()
        # Add all of our angles
        try:
            for i in range(self.parm_data['AMOEBA_REGULAR_ANGLE_NUM_LIST'][0]):
                id1 = self.parm_data['AMOEBA_REGULAR_ANGLE_LIST'][4*i  ] - 1
                id2 = self.parm_data['AMOEBA_REGULAR_ANGLE_LIST'][4*i+1] - 1
                id3 = self.parm_data['AMOEBA_REGULAR_ANGLE_LIST'][4*i+2] - 1
                typ = self.parm_data['AMOEBA_REGULAR_ANGLE_LIST'][4*i+3] - 1
                self.angle_list.append(
                        Angle(self.atom_list[id1], self.atom_list[id2],
                              self.atom_list[id3], self.angle_type_list[typ])
                )
        except KeyError:
            pass
        self.angle_type_list.changed = False
        self.angle_list.changed = False
        ##### Next create our list of trigonal angles (in-plane)
        self.trigonal_angle_type_list = TrigonalAngleTypeList(self)
        self.trigonal_angle_list = TrackedList()
        # Add all trigonal angles
        try:
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
        except KeyError:
            pass
        self.trigonal_angle_type_list.changed = False
        self.trigonal_angle_list.changed = False
        ##### Next create our list of out-of-plane bending terms #####
        self.oopbend_type_list = OutOfPlaneBendTypeList(self)
        self.oopbend_list = TrackedList()
        # Add the out-of-plane bending terms
        try:
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
        except KeyError:
            pass
        self.oopbend_type_list.changed = False
        self.oopbend_list.changed = False
        ##### Next create our list of normal dihedrals #####
        self.dihedral_type_list = DihedralTypeList(self)
        self.dihedral_list = TrackedList()
        # Add the dihedrals
        try:
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
        except KeyError:
            pass
        self.dihedral_type_list.changed = False
        self.dihedral_list.changed = False
        ##### Next create our list of pi-torsions #####
        self.pitorsion_type_list = PiTorsionTypeList(self)
        self.pitorsion_list = TrackedList()
        # Add the pi-torsions
        try:
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
        except KeyError:
            pass
        self.pitorsion_type_list.changed = False
        self.pitorsion_list.changed = False
        ##### Next create stretch-bend terms #####
        self.stretch_bend_type_list = StretchBendTypeList(self)
        self.stretch_bend_list = TrackedList()
        # Add the stretch-bends
        try:
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
        except KeyError:
            pass
        self.stretch_bend_type_list.changed = False
        self.stretch_bend_list.changed = False
        ##### Next create the torsion-torsion parameters #####
        self.torsion_torsion_type_list = TorsionTorsionTypeList(self)
        self.torsion_torsion_list = TrackedList()
        # Add the torsion-torsions
        try:
            ntt = self.parm_data['AMOEBA_TORSION_TORSION_NUM_LIST'][0]
            for i in range(ntt):
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
        except KeyError:
            pass
        self.torsion_torsion_type_list.changed = False
        self.torsion_torsion_list.changed = False
        ##### Next create the chiral frame list #####
        self.chiral_frame_list = TrackedList()
        try:
            for i in range(self.parm_data['AMOEBA_CHIRAL_FRAME_NUM_LIST'][0]):
                id1 = self.parm_data['AMOEBA_CHIRAL_FRAME_LIST'][3*i  ] - 1
                id2 = self.parm_data['AMOEBA_CHIRAL_FRAME_LIST'][3*i+1] - 1
                chi = self.parm_data['AMOEBA_CHIRAL_FRAME_LIST'][3*i+2]
                self.chiral_frame_list.append(
                        ChiralFrame(self.atom_list[id1],
                                    self.atom_list[id2], chi)
                )
        except KeyError:
            pass
        self.chiral_frame_list.changed = False
        ##### Create the multipole frame list #####
        self.multipole_frame_list = TrackedList()
        try:
            for i in range(self.parm_data['AMOEBA_FRAME_DEF_NUM_LIST'][0]):
                id1 = self.parm_data['AMOEBA_FRAME_DEF_LIST'][5*i  ] - 1
                fpn = self.parm_data['AMOEBA_FRAME_DEF_LIST'][5*i+1]
                vct = self.parm_data['AMOEBA_FRAME_DEF_LIST'][5*i+2]
                vch = self.parm_data['AMOEBA_FRAME_DEF_LIST'][5*i+3]
                nvc = self.parm_data['AMOEBA_FRAME_DEF_LIST'][5*i+4]
                self.multipole_frame_list.append(
                        MultipoleFrame(self.atom_list[id1], fpn, vct, vch, nvc)
                )
        except KeyError:
            pass
        self.multipole_frame_list.changed = False
        ##### Create the "adjust" (detailed exclusion) list #####
        self.adjust_list = TrackedList()
        for i in range(self.parm_data['AMOEBA_ADJUST_NUM_LIST'][0]):
            id1 = self.parm_data['AMOEBA_ADJUST_LIST'][3*i  ] - 1
            id2 = self.parm_data['AMOEBA_ADJUST_LIST'][3*i+1] - 1
            wt = self.parm_data['AMOEBA_ADJUST_LIST'][3*i+2]
            self.adjust_list.append(
                    ExclusionAssignment(self.atom_list[id1],
                                        self.atom_list[id2], wt)
            )
        self.adjust_weights = ExclusionAssignmentWeights(self)
        self.adjust_weights.changed = False
        self.adjust_list.changed = False
        ##### Finally we can determine the full exclusion list #####
        for atom in self.atom_list: atom.determine_all_exclusion_partners()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def remake_parm(self):
        """ Recomputes the topology file parameters """
        # First thing we have to do is load any of our old atom parameters into
        # our atom_list to preserve any changes we've made directly to the data
        self.atom_list.refresh_data()
        # Now delete all of the bond/angle/dihedral partner information and
        # refresh it to make sure we get the exclusions right
        for atm in self.atom_list: atm.reset_topology()
        for bnd in self.bond_list: bnd.register()
        # Reset all of the type lists
        self.bond_type_list.reset()
        self.urey_bradley_type_list.reset()
        self.angle_type_list.reset()
        self.trigonal_angle_type_list.reset()
        self.oopbend_type_list.reset()
        self.dihedral_type_list.reset()
        self.pitorsion_type_list.reset()
        self.stretch_bend_type_list.reset()
        self.torsion_torsion_type_list.reset()
        # Shortcut here -- for each atom now just determine the angle/dihedral
        # partners directly from the bonded network
        for atm in self.atom_list: atm.determine_exclusions_from_bonds()
        # Now fill up the arrays of atomic properties. This will also adjust
        # NATOM for us if we've deleted atoms
        self.atom_list.write_to_parm()
        nnb = len(self.parm_data['EXCLUDED_ATOMS_LIST'])
        self.parm_data['POINTERS'][NNB] = nnb
        # Reset the residue indexes
        for res in self.residue_list: res.idx = -1
        # Write the residue arrays
        num_res = 0
        for i, atm in enumerate(self.atom_list):
            if atm.residue.idx == -1:
                self.parm_data['RESIDUE_LABEL'][num_res] = atm.residue.resname
                self.parm_data['RESIDUE_POINTER'][num_res] = i + 1
                atm.residue.idx = num_res
                num_res += 1
        self.parm_data['POINTERS'][NRES] = num_res
        self._truncate_array('RESIDUE_LABEL', num_res)
        self._truncate_array('RESIDUE_POINTER', num_res)
        nmxrs = max([len(r) for r in self.residue_list])
        self.parm_data['POINTERS'][NMXRS] = nmxrs
    
        # Now write all of the bond arrays. We will loop through all of the
        # bonds to make sure that all of their atoms still exist (atm.idx > -1).
        # At the same time, we will start applying indexes to the bond_types so
        # we only print out the bond types that will be used. To do this, we
        # need a copule counters. Different bond types will have an index of -1
        # until we find out that they are needed. Then we assign them an index
        # and write out that bond info. We also have to make sure that every
        # array is at least large enough, so give it enough elements to cover
        # every bond in the list, which will be reduced in size of not every
        # bond is actually added
        if self.bond_list:
            bond_num = typenum = 0
            self.parm_data['AMOEBA_REGULAR_BOND_LIST'] = \
                                            _zeros(len(self.bond_list)*3)
            for i, bnd in enumerate(self.bond_list):
                if -1 in (bnd.atom1.idx, bnd.atom2.idx): continue
                if bnd.bond_type.idx == -1:
                    bnd.bond_type.idx = typenum
                    typenum += 1
                bnd.write_info(self, bond_num)
                bond_num += 1
            # At this point, bond_num is one past the last index used, but since
            # Python indexes from 0 bond_num is actually now equal to the total
            # number of bonds that were written
            self.parm_data['AMOEBA_REGULAR_BOND_NUM_LIST'] = [bond_num]
            self._truncate_array('AMOEBA_REGULAR_BOND_LIST', bond_num*3)
            # Now write the types
            self.parm_data['AMOEBA_REGULAR_BOND_NUM_PARAMS'] = [typenum]
            self.parm_data['AMOEBA_REGULAR_BOND_FORCE_CONSTANT'] = \
                                                _zeros(typenum)
            self.parm_data['AMOEBA_REGULAR_BOND_EQUIL_VALUE'] = _zeros(typenum)
            self.bond_type_list.write_to_parm()
        else:
            self.deleteFlag('AMOEBA_REGULAR_BOND_LIST')
            self.deleteFlag('AMOEBA_REGULAR_BOND_NUM_LIST')
            self.deleteFlag('AMOEBA_REGULAR_BOND_NUM_PARAMS')
            self.deleteFlag('AMOEBA_REGULAR_BOND_FORCE_CONSTANT')
            self.deleteFlag('AMOEBA_REGULAR_BOND_EQUIL_VALUE')

        # Now time for urey-bradleys
        if self.urey_bradley_list:
            ub_num = typenum = 0
            self.parm_data['AMOEBA_UREY_BRADLEY_BOND_LIST'] = \
                                        _zeros(len(self.urey_bradley_list)*3)
            for i, ub in enumerate(self.urey_bradley_list):
                if -1 in (ub.atom1, ub.atom2): continue
                if ub.ub_type.idx == -1:
                    ub.ub_type.idx = typenum
                    typenum += 1
                ub.write_info(self, ub_num)
                ub_num += 1
            self.parm_data['AMOEBA_UREY_BRADLEY_BOND_NUM_LIST'] = [ub_num]
            self._truncate_array('AMOEBA_UREY_BRADLEY_BOND_LIST', ub_num*3)
            # Now the types
            self.parm_data['AMOEBA_UREY_BRADLEY_NUM_PARAMS'] = [typenum]
            self.parm_data['AMOEBA_UREY_BRADLEY_BOND_FORCE_CONSTANT'] = \
                                            _zeros(typenum)
            self.parm_data['AMOEBA_UREY_BRADLEY_BOND_EQUIL_VALUE'] = \
                                            _zeros(typenum)
            self.urey_bradley_type_list.write_to_parm()
        else:
            self.deleteFlag('AMOEBA_UREY_BRADLEY_BOND_NUM_LIST')
            self.deleteFlag('AMOEBA_UREY_BRADLEY_BOND_LIST')
            self.deleteFlag('AMOEBA_UREY_BRADLEY_NUM_PARAMS')
            self.deleteFlag('AMOEBA_UREY_BRADLEY_BOND_FORCE_CONSTANT')
            self.deleteFlag('AMOEBA_UREY_BRADLEY_BOND_EQUIL_VALUE')

        # Now time for angles
        if self.angle_list:
            ang_num = typenum = 0
            self.parm_data['AMOEBA_REGULAR_ANGLE_LIST'] = \
                                            _zeros(len(self.angle_list)*4)
            for i, ang in enumerate(self.angle_list):
                if -1 in (ang.atom1.idx, ang.atom2.idx, ang.atom3.idx): continue
                if ang.angle_type.idx == -1:
                    ang.angle_type.idx = typenum
                    typenum += 1
                ang.write_info(self, ang_num)
                ang_num += 1
            self.parm_data['AMOEBA_REGULAR_ANGLE_NUM_LIST'] = [ang_num]
            self._truncate_array('AMOEBA_REGULAR_ANGLE_LIST', ang_num*4)
            # Now the types
            self.parm_data['AMOEBA_REGULAR_ANGLE_NUM_PARAMS'] = [typenum]
            self.parm_data['AMOEBA_REGULAR_ANGLE_FORCE_CONSTANT'] = \
                                                _zeros(typenum)
            self.parm_data['AMOEBA_REGULAR_ANGLE_EQUIL_VALUE'] = _zeros(typenum)
            self.angle_type_list.write_to_parm()
        else:
            self.deleteFlag('AMOEBA_REGULAR_ANGLE_NUM_LIST')
            self.deleteFlag('AMOEBA_REGULAR_ANGLE_LIST')
            self.deleteFlag('AMOEBA_REGULAR_ANGLE_NUM_PARAMS')
            self.deleteFlag('AMOEBA_REGULAR_ANGLE_FORCE_CONSTANT')
            self.deleteFlag('AMOEBA_REGULAR_ANGLE_EQUIL_VALUE')

        # Now time for the trigonal angles
        if self.trigonal_angle_list:
            ang_num = typenum = 0
            self.parm_data['AMOEBA_TRIGONAL_ANGLE_LIST'] = \
                                        _zeros(len(self.trigonal_angle_list)*5)
            for i, ang in enumerate(self.trigonal_angle_list):
                if -1 in (ang.atom1.idx, ang.atom2.idx,
                          ang.atom3.idx, ang.atom4.idx):
                    continue
                if ang.trigang_type.idx == -1:
                    ang.trigang_type.idx = typenum
                    typenum += 1
                ang.write_info(self, ang_num)
                ang_num += 1
            self.parm_data['AMOEBA_TRIGONAL_ANGLE_NUM_LIST'] = [ang_num]
            self._truncate_array('AMOEBA_TRIGONAL_ANGLE_LIST', ang_num*5)
            # Now the types
            self.parm_data['AMOEBA_TRIGONAL_ANGLE_NUM_PARAMS'] = [typenum]
            self.parm_data['AMOEBA_TRIGONAL_ANGLE_FORCE_CONSTANT'] = \
                                        _zeros(typenum)
            self.parm_data['AMOEBA_TRIGONAL_ANGLE_EQUIL_VALUE'] = \
                                        _zeros(typenum)
            self.trigonal_angle_type_list.write_to_parm()
        else:
            self.deleteFlag('AMOEBA_TRIGONAL_ANGLE_NUM_LIST')
            self.deleteFlag('AMOEBA_TRIGONAL_ANGLE_LIST')
            self.deleteFlag('AMOEBA_TRIGONAL_ANGLE_NUM_PARAMS')
            self.deleteFlag('AMOEBA_TRIGONAL_ANGLE_FORCE_CONSTANT')
            self.deleteFlag('AMOEBA_TRIGONAL_ANGLE_EQUIL_VALUE')

        # Now time for the out-of-plane bending terms
        if self.oopbend_list:
            oop_num = typenum = 0
            self.parm_data['AMOEBA_OPBEND_ANGLE_LIST'] = \
                                            _zeros(len(self.oopbend_list)*5)
            for i, oop in enumerate(self.oopbend_list):
                if -1 in (oop.atom1.idx, oop.atom2.idx,
                          oop.atom3.idx, oop.atom4.idx):
                    continue
                if oop.oopbend_type.idx == -1:
                    oop.oopbend_type.idx = typenum
                    typenum += 1
                oop.write_info(self, oop_num)
                oop_num += 1
            self.parm_data['AMOEBA_OPBEND_ANGLE_NUM_LIST'] = [oop_num]
            self._truncate_array('AMOEBA_OPBEND_ANGLE_LIST', oop_num*5)
            # Now the types
            self.parm_data['AMOEBA_OPBEND_ANGLE_NUM_PARAMS'] = [typenum]
            self.parm_data['AMOEBA_OPBEND_ANGLE_FORCE_CONSTANT'] = \
                                            _zeros(typenum)
            self.parm_data['AMOEBA_OPBEND_ANGLE_EQUIL_VALUE'] = _zeros(typenum)
            self.oopbend_type_list.write_to_parm()
        else:
            self.deleteFlag('AMOEBA_OPBEND_ANGLE_NUM_LIST')
            self.deleteFlag('AMOEBA_OPBEND_ANGLE_LIST')
            self.deleteFlag('AMOEBA_OPBEND_ANGLE_NUM_PARAMS')
            self.deleteFlag('AMOEBA_OPBEND_ANGLE_FORCE_CONSTANT')
            self.deleteFlag('AMOEBA_OPBEND_ANGLE_EQUIL_VALUE')

        # Now time for torsions
        if self.dihedral_list:
            dih_num = typenum = 0
            self.parm_data['AMOEBA_TORSION_LIST'] = \
                                            _zeros(len(self.dihedral_list)*5)
            for i, dih in enumerate(self.dihedral_list):
                if -1 in (dih.atom1.idx, dih.atom2.idx,
                          dih.atom3.idx, dih.atom4.idx):
                    continue
                if dih.dihedral_type.idx == -1:
                    dih.dihedral_type.idx = typenum
                    typenum += 1
                dih.write_info(self, dih_num)
                dih_num += 1
            self.parm_data['AMOEBA_TORSION_NUM_LIST'] = [dih_num]
            self._truncate_array('AMOEBA_TORSION_LIST', dih_num*5)
            # Now the types
            self.parm_data['AMOEBA_TORSION_NUM_PARAMS'] = [typenum]
            self.parm_data['AMOEBA_TORSION_FORCE_CONSTANT'] = _zeros(typenum)
            self.parm_data['AMOEBA_TORSION_PERIODICITY'] = _zeros(typenum)
            self.parm_data['AMOEBA_TORSION_PHASE'] = _zeros(typenum)
            self.dihedral_type_list.write_to_parm()
        else:
            self.deleteFlag('AMOEBA_TORSION_NUM_LIST')
            self.deleteFlag('AMOEBA_TORSION_LIST')
            self.deleteFlag('AMOEBA_TORSION_NUM_PARAMS')
            self.deleteFlag('AMOEBA_TORSION_FORCE_CONSTANT')
            self.deleteFlag('AMOEBA_TORSION_PERIODICITY')
            self.deleteFlag('AMOEBA_TORSION_PHASE')

        # Now time for pi-torsions
        if self.pitorsion_list:
            tor_num = typenum = 0
            self.parm_data['AMOEBA_PI_TORSION_LIST'] = \
                                            _zeros(len(self.pitorsion_list)*7)
            for i, tor in enumerate(self.pitorsion_list):
                if -1 in (tor.atom1.idx, tor.atom2.idx, tor.atom3.idx,
                        tor.atom4.idx, tor.atom5.idx, tor.atom6.idx):
                    continue
                if tor.pitor_type.idx == -1:
                    tor.pitor_type.idx = typenum
                    typenum += 1
                tor.write_info(self, tor_num)
                tor_num += 1
            self.parm_data['AMOEBA_PI_TORSION_NUM_LIST'] = [tor_num]
            self._truncate_array('AMOEBA_PI_TORSION_LIST', tor_num*7)
            # Now the types
            self.parm_data['AMOEBA_PI_TORSION_NUM_PARAMS'] = [typenum]
            self.parm_data['AMOEBA_PI_TORSION_FORCE_CONSTANT'] = _zeros(typenum)
            self.parm_data['AMOEBA_PI_TORSION_PERIODICITY'] = _zeros(typenum)
            self.parm_data['AMOEBA_PI_TORSION_PHASE'] = _zeros(typenum)
            self.pitorsion_type_list.write_to_parm()
        else:
            self.deleteFlag('AMOEBA_PI_TORSION_NUM_LIST')
            self.deleteFlag('AMOEBA_PI_TORSION_LIST')
            self.deleteFlag('AMOEBA_PI_TORSION_NUM_PARAMS')
            self.deleteFlag('AMOEBA_PI_TORSION_FORCE_CONSTANT')
            self.deleteFlag('AMOEBA_PI_TORSION_PERIODICITY')
            self.deleteFlag('AMOEBA_PI_TORSION_PHASE')

        # Now time for the stretch-bends
        if self.stretch_bend_list:
            strb_num = typenum = 0
            self.parm_data['AMOEBA_STRETCH_BEND_LIST'] = \
                                        _zeros(len(self.stretch_bend_list)*4)
            for i, strb in enumerate(self.stretch_bend_list):
                if -1 in (strb.atom1.idx, strb.atom2.idx, strb.atom3.idx):
                    continue
                if strb.strbnd_type.idx == -1:
                    strb.strbnd_type.idx = typenum
                    typenum += 1
                strb.write_info(self, strb_num)
                strb_num += 1
            self.parm_data['AMOEBA_STRETCH_BEND_NUM_LIST'] = [strb_num]
            self._truncate_array('AMOEBA_STRETCH_BEND_LIST', strb_num*4)
            # Now the types
            self.parm_data['AMOEBA_STRETCH_BEND_NUM_PARAMS'] = [typenum]
            self.parm_data['AMOEBA_STRETCH_BEND_FORCE_CONSTANT'] = \
                                        _zeros(typenum)
            self.parm_data['AMOEBA_STRETCH_BEND_EQUIL_VALUE'] = _zeros(typenum)
            self.parm_data['AMOEBA_STRETCH_BEND_BOND1_EQUIL_VALUE'] = \
                                        _zeros(typenum)
            self.parm_data['AMOEBA_STRETCH_BEND_BOND2_EQUIL_VALUE'] = \
                                        _zeros(typenum)
            self.stretch_bend_type_list.write_to_parm()
        else:
            self.deleteFlag('AMOEBA_STRETCH_BEND_NUM_LIST')
            self.deleteFlag('AMOEBA_STRETCH_BEND_LIST')
            self.deleteFlag('AMOEBA_STRETCH_BEND_NUM_PARAMS')
            self.deleteFlag('AMOEBA_STRETCH_BEND_FORCE_CONSTANT')
            self.deleteFlag('AMOEBA_STRETCH_BEND_EQUIL_VALUE')
            self.deleteFlag('AMOEBA_STRETCH_BEND_BOND1_EQUIL_VALUE')
            self.deleteFlag('AMOEBA_STRETCH_BEND_BOND2_EQUIL_VALUE')

        # Now time for the coupled torsions
        if self.torsion_torsion_list:
            tor_num = 0
            typelist = []
            self.parm_data['AMOEBA_TORSION_TORSION_LIST'] = \
                                    _zeros(len(self.torsion_torsion_list)*6)
            for i, tor in enumerate(self.torsion_torsion_list):
                if -1 in (tor.atom1.idx, tor.atom2.idx, tor.atom3.idx,
                          tor.atom4.idx, tor.atom5.idx):
                    continue
                if tor.tortor_type.idx == -1:
                    tor.tortor_type.idx = len(typelist)
                    typelist.append(tor.tortor_type)
                tor.write_info(self, tor_num)
                tor_num += 1
            self.parm_data['AMOEBA_TORSION_TORSION_NUM_LIST'] = [tor_num]
            self._truncate_array('AMOEBA_TORSION_TORSION_LIST', tor_num*6)
            # Now the types
            self.parm_data['AMOEBA_TORSION_TORSION_NUM_PARAMS'] = \
                                    [len(typelist)]
            # We have to delete all of the 'old' types
            for flag in self.flag_list:
                if flag.startswith('AMOEBA_TORSION_TORSION_TORTOR_TABLE'):
                    self.deleteFlag(flag)
            # Now add them back for all of the types
            after = 'AMOEBA_TORSION_TORSION_NUM_PARAMS'
            for tortype in typelist:
                tortype.write_info(self, after)
                after = ('AMOEBA_TORSION_TORSION_TORTOR_TABLE_%02d_'
                         'D2FUNC_DANGLE1_DANGLE2' % (tortype.idx+1))
        else:
            self.deleteFlag('AMOEBA_TORSION_TORSION_NUM_LIST')
            self.deleteFlag('AMOEBA_TORSION_TORSION_LIST')
            self.deleteFlag('AMOEBA_TORSION_TORSION_NUM_PARAMS')
            for flag in self.flag_list:
                if flag.startswith('AMOEBA_TORSION_TORSION_TORTOR_TABLE'):
                    self.deleteFlag(flag)

        # Now time for the chiral frames
        if self.chiral_frame_list:
            chi_num = 0
            self.parm_data['AMOEBA_CHIRAL_FRAME_LIST'] = \
                                        _zeros(len(self.chiral_frame_list)*3)
            for i, chi in enumerate(self.chiral_frame_list):
                if -1 in (chi.atom1, chi.atom2): continue
                chi.write_info(self, chi_num)
                chi_num += 1
            self.parm_data['AMOEBA_CHIRAL_FRAME_NUM_LIST'] = [chi_num]
            self._truncate_array('AMOEBA_CHIRAL_FRAME_LIST', chi_num*3)
        else:
            self.deleteFlag('AMOEBA_CHIRAL_FRAME_NUM_LIST')
            self.deleteFlag('AMOEBA_CHIRAL_FRAME_LIST')

        # Now time for the multipole frames
        if self.multipole_frame_list:
            mul_num = 0
            self.parm_data['AMOEBA_FRAME_DEF_LIST'] = \
                                    _zeros(len(self.multipole_frame_list)*5)
            for i, mul in enumerate(self.multipole_frame_list):
                if mul.atom.idx == -1: continue
                mul.write_info(self, mul_num)
                mul_num += 1
            self.parm_data['AMOEBA_FRAME_DEF_NUM_LIST'] = [mul_num]
            self._truncate_array('AMOEBA_FRAME_DEF_NUM_LIST', mul_num*5)
        else:
            self.deleteFlag('AMOEBA_FRAME_DEF_NUM_LIST')
            self.deleteFlag('AMOEBA_FRAME_DEF_LIST')

        # Now time for the adjust array
        if self.adjust_list:
            adj_num = 0
            self.parm_data['AMOEBA_ADJUST_LIST'] = \
                                    _zeros(len(self.adjust_list)*3)
            for i, adj in enumerate(self.adjust_list):
                if -1 in (adj.atom1.idx, adj.atom2.idx): continue
                adj.write_info(self, adj_num)
                adj_num += 1
            self.parm_data['AMOEBA_ADJUST_NUM_LIST'] = [adj_num]
            self._truncate_array('AMOEBA_ADJUST_LIST', adj_num*3)
            # Now for the weights
            self.adjust_weights.write_to_parm()
        else:
            self.deleteFlag('AMOEBA_ADJUST_NUM_LIST')
            self.deleteFlag('AMOEBA_ADJUST_LIST')

        # Mark all lists as *not* changed
        self.bond_type_list.changed = False
        self.bond_list.changed = False
        self.urey_bradley_type_list.changed = False
        self.urey_bradley_list.changed = False
        self.angle_type_list.changed = False
        self.angle_list.changed = False
        self.trigonal_angle_type_list.changed = False
        self.trigonal_angle_list.changed = False
        self.oopbend_type_list.changed = False
        self.oopbend_list.changed = False
        self.dihedral_type_list.changed = False
        self.dihedral_list.changed = False
        self.pitorsion_type_list.changed = False
        self.pitorsion_list.changed = False
        self.stretch_bend_type_list.changed = False
        self.stretch_bend_list.changed = False
        self.torsion_torsion_type_list.changed = False
        self.torsion_torsion_list.changed = False
        self.chiral_frame_list.changed = False
        self.multipole_frame_list.changed = False
        self.adjust_weights.changed = False
        self.adjust_list.changed = False

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def _topology_changed(self):
        return (self.bond_type_list.changed or
                self.bond_list.changed or
                self.urey_bradley_type_list.changed or
                self.urey_bradley_list.changed or
                self.angle_type_list.changed or
                self.angle_list.changed or
                self.trigonal_angle_list.changed or
                self.trigonal_angle_type_list.changed or
                self.oopbend_type_list.changed or
                self.oopbend_list.changed or
                self.dihedral_type_list.changed or
                self.dihedral_list.changed or
                self.pitorsion_type_list.changed or
                self.pitorsion_list.changed or
                self.stretch_bend_type_list.changed or
                self.stretch_bend_list.changed or
                self.torsion_torsion_type_list.changed or
                self.torsion_torsion_list.changed or
                self.chiral_frame_list.changed or
                self.adjust_weights.changed or
                self.adjust_list.changed)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def mdin_skeleton(self):
        """
        Returns the skeleton of an mdin file with the &amoeba namelist set up
        correctly for the potential terms that are present in this topology
        file.

        Returns
        -------
        str
            A skeleton MDIN file with all of the do_* variables in the &amoeba
            section set correctly. It is commented for easy editing
        """
        return ('Input file for AMOEBA simulations.\n'
                ' &cntrl\n'
                '     ! Add whatever variables you need here\n'
                '     ntb=1, ntt=1, ntp=0, ! PBC, thermostat, barostat\n'
                '     irest=0, ntx=1,      ! restart flags\n'
                ' /\n'
                ' &amoeba\n'
                '     ! Some basic potential parameters. For better\n'
                '     ! energy conservation you need to adjust these\n'
                '     ! defaults\n'
                '     beeman_integrator=1,   ! Use Beeman integrator\n'
                '     dipole_scf_tol=0.01,   ! 10e-6 gives good NVE\n'
                '\n'
                '     ! You should not generally modify these variables:\n'
                '     do_valence=1, do_bond=%d, do_ureyb=%d,\n'
                '     do_reg_angle=%d, do_trig_angle=%d, do_opbend=%d,\n'
                '     do_torsion=%d, do_pi_torsion=%d, do_strbend=%d,\n'
                '     do_torsion_torsion=%d,\n'
                ' /\n' % (bool(self.bond_list), bool(self.urey_bradley_list),
                bool(self.angle_list), bool(self.trigonal_angle_list),
                bool(self.oopbend_list), bool(self.dihedral_list),
                bool(self.pitorsion_list), bool(self.stretch_bend_list),
                bool(self.torsion_torsion_list))
        )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    @property
    def bonds_inc_h(self):
        """
        For compatibility with Amber and Chamber-style prmtops. This is a
        generator, so it cannot be used everywhere that a bonds_inc_h array can
        be used
        """
        for bnd in self.bond_list:
            if 1 in (bnd.atom1.element, bnd.atom2.element):
                yield bnd
   
    @bonds_inc_h.setter
    def bonds_inc_h(self, thing):
        raise NotImplemented('bonds_inc_h is a generator property '
                             'in AmoebaParm')

    @property
    def bonds_without_h(self):
        for bnd in self.bond_list:
            if not 1 in (bnd.atom1.element, bnd.atom2.element):
                yield bnd

    @bonds_without_h.setter
    def bonds_without_h(self, thing):
        raise NotImplemented('bonds_without_h is a generator property in '
                             'AmoebaParm')

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    @property
    def chamber(self):
        return False
   
    @property
    def amoeba(self):
        return True
