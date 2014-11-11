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

from chemistry import periodic_table
from chemistry.amber.topologyobjects import (Bond, Angle,
            Dihedral, ResidueList, AtomList, BondTypeList,
            AngleTypeList, DihedralTypeList, TrackedList)
from chemistry.amber.constants import (NATOM, NTYPES, NBONH, MBONA, NTHETH,
            MTHETA, NPHIH, MPHIA, NHPARM, NPARM, NEXT, NRES, NBONA, NTHETA,
            NPHIA, NUMBND, NUMANG, NPTRA, NATYP, NPHB, IFPERT, NBPER, NGPER,
            NDPER, MBPER, MGPER, MDPER, IFBOX, NMXRS, IFCAP, NUMEXTRA, NCOPY,
            NNB)
from chemistry.amber.amberformat import AmberFormat
from chemistry.exceptions import (AmberParmWarning, AmberParmError, ReadError,
                                  MoleculeError, MoleculeWarning)
from warnings import warn
from math import sqrt

class AmberParm(AmberFormat):
    """
    Amber Topology (parm7 format) class. Gives low, and some high, level access
    to topology data.
    """

    solvent_residues = ['WAT', 'HOH']

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def __init__(self, prm_name=None, rst7_name=None):
        """
        Instantiates an AmberParm object from data in prm_name and establishes
        validity based on presence of POINTERS and CHARGE sections. In general,
        you should use LoadParm from the readparm module instead. LoadParm will
        correctly dispatch the object to the 'correct' flavor of AmberParm
        """
        AmberFormat.__init__(self, prm_name)
        self.hasvels = self.hasbox = False
        if prm_name is not None:
            self.initialize_topology(rst7_name)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def initialize_topology(self, rst7_name=None):
        """
        Initializes topology data structures, like the list of atoms, bonds,
        etc., after the topology file has been read.
        """
        # We need to handle RESIDUE_ICODE properly since it may have picked up
        # some extra values
        if 'RESIDUE_ICODE' in self.flag_list:
            self._truncate_array('RESIDUE_ICODE',
                                 self.parm_data['POINTERS'][NRES])

        # instance variables other than those in AmberFormat
        self.pointers = {}   # list of all the pointers in the prmtop
        self.LJ_types = {}   # dict pairing atom name with its LJ atom type #
        self.LJ_radius = []  # ordered array of L-J radii in Ang -- indices
                             # are elements in LJ_types-1
        self.LJ_depth = []   # similarly ordered array of L-J depths

        # If we were given a prmtop, read it in
        self.LoadPointers()
        self.fill_LJ()

        # Load the structure arrays
        self._load_structure()
        # Find any extra exclusion rules that may be defined
        self.atom_list.find_extra_exclusions()

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

        # If we have coordinates or velocities, load them into the atom list
        if hasattr(self, 'coords'):
            for i, atom in enumerate(self.atom_list):
                i3 = i * 3
                atom.xx, atom.xy, atom.xz = self.coords[i3:i3+3]
        if hasattr(self, 'vels'):
            for i, atom in enumerate(self.atom_list):
                i3 = i * 3
                atom.vx, atom.vy, atom.vz = self.coords[i3:i3+3]

        if rst7_name is not None:
            self.LoadRst7(rst7_name)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    @classmethod
    def load_from_rawdata(cls, rawdata):
        """
        Take the raw data from a AmberFormat object and initialize an AmberParm
        from that data.

        Parameters:
            - rawdata (AmberFormat): Already has a parsed file

        Returns:
            Populated AmberParm instance
        """
        inst = cls()
        inst.prm_name = rawdata.prm_name
        inst.version = rawdata.version
        inst.formats = rawdata.formats
        inst.parm_data = rawdata.parm_data
        inst.parm_comments = rawdata.parm_comments
        inst.flag_list = rawdata.flag_list
        inst.initialize_topology()
        # Convert charges if necessary due to differences in electrostatic
        # scaling factors
        chgscale = rawdata.CHARGE_SCALE / cls.CHARGE_SCALE
        for i in xrange(len(inst.parm_data['CHARGE'])):
            inst.parm_data['CHARGE'][i] *= chgscale
        # See if the rawdata has any kind of structural attributes, like rst7
        # (coordinates) and an atom list with positions and/or velocities
        if hasattr(rawdata, 'rst7'):
            inst.rst7 = rawdata.rst7
        if hasattr(rawdata, 'coords'):
            inst.load_coordinates(rawdata.coords)
        if hasattr(rawdata, 'vels'):
            inst.load_velocities(rawdata.vels)
        if hasattr(rawdata, 'box'):
            inst.box = rawdata.box
        if hasattr(rawdata, 'hasbox'):
            inst.hasbox = rawdata.hasbox
        if hasattr(rawdata, 'hasvels'):
            inst.hasvels = rawdata.hasvels
        return inst
   
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def __copy__(self):
        """ Needs to copy a few additional data structures """
        other = super(AmberParm, self).__copy__()
        other.pointers = {}
        other.LJ_types = {}
        other.LJ_radius = self.LJ_radius[:]
        other.LJ_depth = self.LJ_depth[:]
        other.hasvels = self.hasvels
        other.hasbox = self.hasbox

        # Now fill the LJ and other data structures
        for p in self.pointers: other.pointers[p] = self.pointers[p]
        for typ in self.LJ_types: other.LJ_types[typ] = self.LJ_types[typ]
        try:
            other._load_structure()
        except (KeyError, IndexError, AttributeError):
            raise AmberParmError('Could not set up topology for parm copy')
        # See if we have a restart file
        if hasattr(self, 'rst7'):
            other.rst7 = Rst7.copy_from(self.rst7)
        # Now we should have a full copy
        return other

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
    def LoadPointers(self):
        """
        Loads the data in POINTERS section into a pointers dictionary with each
        key being the pointer name according to http://ambermd.org/formats.html
        """
        self.pointers["NATOM"] = self.parm_data["POINTERS"][NATOM]
        self.pointers["NTYPES"] = self.parm_data["POINTERS"][NTYPES]
        self.pointers["NBONH"] = self.parm_data["POINTERS"][NBONH]
        self.pointers["MBONA"] = self.parm_data["POINTERS"][MBONA]
        self.pointers["NTHETH"] = self.parm_data["POINTERS"][NTHETH]
        self.pointers["MTHETA"] = self.parm_data["POINTERS"][MTHETA]
        self.pointers["NPHIH"] = self.parm_data["POINTERS"][NPHIH]
        self.pointers["MPHIA"] = self.parm_data["POINTERS"][MPHIA]
        self.pointers["NHPARM"] = self.parm_data["POINTERS"][NHPARM]
        self.pointers["NPARM"] = self.parm_data["POINTERS"][NPARM]
        self.pointers["NEXT"] = self.parm_data["POINTERS"][NEXT]
        self.pointers["NNB"] = self.parm_data["POINTERS"][NNB] # alias for above
        self.pointers["NRES"] = self.parm_data["POINTERS"][NRES]
        self.pointers["NBONA"] = self.parm_data["POINTERS"][NBONA]
        self.pointers["NTHETA"] = self.parm_data["POINTERS"][NTHETA]
        self.pointers["NPHIA"] = self.parm_data["POINTERS"][NPHIA]
        self.pointers["NUMBND"] = self.parm_data["POINTERS"][NUMBND]
        self.pointers["NUMANG"] = self.parm_data["POINTERS"][NUMANG]
        self.pointers["NPTRA"] = self.parm_data["POINTERS"][NPTRA]
        self.pointers["NATYP"] = self.parm_data["POINTERS"][NATYP]
        self.pointers["NPHB"] = self.parm_data["POINTERS"][NPHB]
        self.pointers["IFPERT"] = self.parm_data["POINTERS"][IFPERT]
        self.pointers["NBPER"] = self.parm_data["POINTERS"][NBPER]
        self.pointers["NGPER"] = self.parm_data["POINTERS"][NGPER]
        self.pointers["NDPER"] = self.parm_data["POINTERS"][NDPER]
        self.pointers["MBPER"] = self.parm_data["POINTERS"][MBPER]
        self.pointers["MGPER"] = self.parm_data["POINTERS"][MGPER]
        self.pointers["MDPER"] = self.parm_data["POINTERS"][MDPER]
        self.pointers["IFBOX"] = self.parm_data["POINTERS"][IFBOX]
        self.pointers["NMXRS"] = self.parm_data["POINTERS"][NMXRS]
        self.pointers["IFCAP"] = self.parm_data["POINTERS"][IFCAP]
        self.pointers["NUMEXTRA"] = self.parm_data["POINTERS"][NUMEXTRA]
        if self.parm_data['POINTERS'][IFBOX] > 0:
            self.pointers['IPTRES'] = self.parm_data['SOLVENT_POINTERS'][0]
            self.pointers['NSPM'] = self.parm_data['SOLVENT_POINTERS'][1]
            self.pointers['NSPSOL'] = self.parm_data['SOLVENT_POINTERS'][2]
        # The next is probably only there for LES-prmtops
        try:
            self.pointers["NCOPY"] = self.parm_data["POINTERS"][NCOPY]
        except:
            pass

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def _load_structure(self):
        """ 
        Loads all of the topology instance variables. This is necessary if we
        actually want to modify the topological layout of our system
        (like deleting atoms)
        """
        ##### First create our atoms #####
        self.atom_list = AtomList(self)
        ##### Next, load our residues #####
        self.residue_list = ResidueList(self)
        ##### Next create our list of bonds #####
        self.bond_type_list = BondTypeList(self)
        self.bonds_inc_h, self.bonds_without_h = TrackedList(), TrackedList()
        # Array of bonds with hydrogen
        for i in xrange(self.ptr('nbonh')):
            blist = self.parm_data['BONDS_INC_HYDROGEN']
            self.bonds_inc_h.append(
                    Bond(self.atom_list[blist[3*i  ]//3],
                         self.atom_list[blist[3*i+1]//3],
                         self.bond_type_list[blist[3*i+2]-1])
            )
        # Array of bonds without hydrogen
        for i in xrange(self.ptr('mbona')):
            blist = self.parm_data['BONDS_WITHOUT_HYDROGEN']
            self.bonds_without_h.append(
                    Bond(self.atom_list[blist[3*i  ]//3],
                         self.atom_list[blist[3*i+1]//3],
                         self.bond_type_list[blist[3*i+2]-1])
            )
        # We haven't changed yet...
        self.bonds_inc_h.changed = self.bonds_without_h.changed = False
        ##### Next create our list of angles #####
        self.angle_type_list = AngleTypeList(self)
        self.angles_inc_h, self.angles_without_h = TrackedList(), TrackedList()
        # Array of angles with hydrogen
        for i in xrange(self.ptr('ntheth')):
            alist = self.parm_data['ANGLES_INC_HYDROGEN']
            self.angles_inc_h.append(
                    Angle(self.atom_list[alist[4*i  ]//3],
                          self.atom_list[alist[4*i+1]//3],
                          self.atom_list[alist[4*i+2]//3],
                          self.angle_type_list[alist[4*i+3]-1])
            )
        # Array of angles without hydrogen
        for i in xrange(self.ptr('mtheta')):
            alist = self.parm_data['ANGLES_WITHOUT_HYDROGEN']
            self.angles_without_h.append(
                    Angle(self.atom_list[alist[4*i  ]//3],
                          self.atom_list[alist[4*i+1]//3],
                          self.atom_list[alist[4*i+2]//3],
                          self.angle_type_list[alist[4*i+3]-1])
            )
        # We haven't changed yet
        self.angles_inc_h.changed = self.angles_without_h.changed = False
        ##### Next create our list of dihedrals #####
        self.dihedral_type_list = DihedralTypeList(self)
        self.dihedrals_inc_h = TrackedList()
        self.dihedrals_without_h = TrackedList()
        # Array of dihedrals with hydrogen
        for i in xrange(self.ptr('nphih')):
            dlist = self.parm_data['DIHEDRALS_INC_HYDROGEN']
            signs = [1,1]
            if dlist[5*i+2] < 0: signs[0] = -1
            if dlist[5*i+3] < 0: signs[1] = -1
            self.dihedrals_inc_h.append(
                    Dihedral(self.atom_list[dlist[5*i  ]//3],
                             self.atom_list[dlist[5*i+1]//3],
                             self.atom_list[abs(dlist[5*i+2]//3)],
                             self.atom_list[abs(dlist[5*i+3]//3)],
                             self.dihedral_type_list[dlist[5*i+4]-1],
                             signs)
            )
        # Array of dihedrals without hydrogen
        for i in xrange(self.ptr('mphia')):
            dlist = self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN']
            signs = [1,1]
            if dlist[5*i+2] < 0: signs[0] = -1
            if dlist[5*i+3] < 0: signs[1] = -1
            self.dihedrals_without_h.append(
                    Dihedral(self.atom_list[dlist[5*i  ]//3],
                             self.atom_list[dlist[5*i+1]//3],
                             self.atom_list[abs(dlist[5*i+2]//3)],
                             self.atom_list[abs(dlist[5*i+3]//3)],
                             self.dihedral_type_list[dlist[5*i+4]-1],
                             signs)
            )
        # We haven't changed yet
        self.dihedrals_inc_h.changed = self.dihedrals_without_h.changed = False

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def __str__(self):
        " Returns the name of the topology file as its string representation "
        if self.prm_name is not None:
            return self.prm_name
        return repr(self)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def ptr(self,pointer):
        """
        Returns the value of the given pointer, and converts to upper-case so
        it's case-insensitive. A non-existent pointer meets with a KeyError
        """
        return self.pointers[pointer.upper()]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def writeRst7(self, name, netcdf=None):
        """
        Writes a restart file with the current coordinates and velocities and
        box info if it's present
        """
        # By default, determine file type by extension (.ncrst is NetCDF)
        netcdf = netcdf or (netcdf is None and name.endswith('.ncrst'))

        # Check that we have a rst7 loaded, then overwrite it with a new one if
        # necessary
        if not hasattr(self, 'rst7'):
#           raise AmberParmError('No coordinates loaded. Cannot write restart')
            self.rst7 = Rst7(hasbox=self.hasbox)
            if self.hasbox:
                self.rst7.box = self.box

        # Now fill in the rst7 coordinates
        self.rst7.natom = len(self.atom_list)
        self.rst7.coordinates = [0.0 for i in xrange(len(self.atom_list)*3)]
        if self.rst7.hasvels:
            self.rst7.velocities = [0.0 for i in xrange(len(self.atom_list)*3)]

        for i, at in enumerate(self.atom_list):
            self.rst7.coordinates[3*i  ] = at.xx
            self.rst7.coordinates[3*i+1] = at.xy
            self.rst7.coordinates[3*i+2] = at.xz
            if self.rst7.hasvels:
                self.rst7.velocities[3*i  ] = at.vx
                self.rst7.velocities[3*i+1] = at.vy
                self.rst7.velocities[3*i+2] = at.vz

        # Now write the restart file
        self.rst7.write(name, netcdf)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    flush_data_changes = _load_structure

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def writeParm(self, name):
        """
        Writes the current data in parm_data into a new topology file with a
        given name.
        """
        if hasattr(self, 'atom_list') and self._topology_changed(): 
            self.remake_parm()
            # Reload the structure now that we've recalculated it to flush all
            # data structures to what they *should* be.
            self._load_structure()
            # Now we have to redo the ATOMS_PER_MOLECULE/SOLVENT_POINTERS
            # sections
            if self.ptr('ifbox'): self.rediscover_molecules()

        AmberFormat.writeParm(self, name)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def remake_parm(self):
        """
        Re-fills the topology file arrays if we have changed the underlying
        structure
        """
        # First thing we have to do is load any of our old atom parameters into
        # our atom_list to preserve any changes we've made directly to the data
        self.atom_list.refresh_data()
        # Now delete all of the bond/angle/dihedral partner information and
        # refresh it to make sure we get the exclusions right
        for atm in self.atom_list: atm.reset_topology()
        for bnd in self.bonds_inc_h: bnd.register()
        for bnd in self.bonds_without_h: bnd.register()
        for ang in self.angles_inc_h: ang.register()
        for ang in self.angles_without_h: ang.register()
        for dih in self.dihedrals_inc_h: dih.register()
        for dih in self.dihedrals_without_h: dih.register()
        # Reset all type lists
        self.bond_type_list.reset()
        self.angle_type_list.reset()
        self.dihedral_type_list.reset()
        # Fill up the atom arrays. This will also adjust NATOM for us if 
        # we've deleted atoms
        self.atom_list.write_to_parm()

        # Recount number of extra points
        nextra = 0
        for atm in self.atom_list:
            if atm.attype[:2] in ['EP', 'LP']: nextra += 1
        self.parm_data['POINTERS'][NUMEXTRA] = nextra

        nnb = self.parm_data['EXCLUDED_ATOMS_LIST']
        self.parm_data['POINTERS'][NNB] = len(nnb)
        # Reset the residue indexes
        for res in self.residue_list: res.idx = -1
        # Write the residue arrays
        num_res = 0
        for i, atm in enumerate(self.atom_list):
            if atm.residue.idx == -1:
                self.parm_data['RESIDUE_LABEL'][num_res] = atm.residue.resname
                self.parm_data['RESIDUE_POINTER'][num_res] = i+1
                num_res += 1
                atm.residue.idx = num_res # index from 1
        self.parm_data['POINTERS'][NRES] = num_res
        self._truncate_array('RESIDUE_LABEL', num_res)
        self._truncate_array('RESIDUE_POINTER', num_res)

        # Adjust NMXRS (number of atoms in largest residue) in case that changed
        bigres = max([len(r) for r in self.residue_list])
        self.parm_data['POINTERS'][NMXRS] = bigres

        # Now write all of the bond arrays. We will loop through all of the
        # bonds to make sure that all of their atoms still exist (atm.idx > -1).
        # At the same time, we will start applying indexes to the bond_types so
        # we only print out the bond types that will be used. To do this, we
        # need a couple counters. Different bond types will have an index of -1
        # until we find out they are needed. Then we assign them an index and
        # write out that bond info. We also have to make sure that every array
        # is at least large enough, so give it enough elements to cover every
        # bond in the list which will be reduced in size if not every bond is
        # actually added
        bond_num = 0
        bond_type_num = 0
        self.parm_data['BONDS_INC_HYDROGEN'] = _zeros(len(self.bonds_inc_h)*3)
        for i, bnd in enumerate(self.bonds_inc_h):
            if -1 in (bnd.atom1.idx, bnd.atom2.idx): continue
            if bnd.bond_type.idx == -1:
                bnd.bond_type.idx = bond_type_num
                bond_type_num += 1
            bnd.write_info(self, 'BONDS_INC_HYDROGEN', bond_num)
            bond_num += 1
        self.parm_data['POINTERS'][NBONH] = bond_num
        self._truncate_array('BONDS_INC_HYDROGEN', 3*bond_num)
        # Now we know how many bonds with hydrogen we have. Note that bond_num
        # is +1 past the last index used, but that last index is -1 from total
        # bond number due to indexing from 0, so it's just right now. So is the
        # bond_type index, but that is applicable for the bonds_without_h
        # as well.
        bond_num = 0
        self.parm_data['BONDS_WITHOUT_HYDROGEN'] = \
                        _zeros(len(self.bonds_without_h)*3)
        for i, bnd in enumerate(self.bonds_without_h):
            if -1 in (bnd.atom1.idx, bnd.atom2.idx): continue
            if bnd.bond_type.idx == -1:
                bnd.bond_type.idx = bond_type_num
                bond_type_num += 1
            bnd.write_info(self, 'BONDS_WITHOUT_HYDROGEN', bond_num)
            bond_num += 1
        # Make sure BOND_FORCE_CONSTANT and BOND_EQUIL_VALUE arrays are big
        # enough
        self.parm_data['BOND_FORCE_CONSTANT'] = _zeros(bond_type_num)
        self.parm_data['BOND_EQUIL_VALUE'] = _zeros(bond_type_num)
        # Now we can write all of the bond types out
        self.bond_type_list.write_to_parm()
        # Now we know how many bonds without H we have and our # of bond types
        self.parm_data['POINTERS'][MBONA] = bond_num
        self.parm_data['POINTERS'][NBONA] = bond_num
        self.parm_data['POINTERS'][NUMBND] = bond_type_num
        self._truncate_array('BONDS_WITHOUT_HYDROGEN', 3*bond_num)

        # Now do all of the angle arrays
        angle_num = 0
        angle_type_num = 0
        # Make sure we have enough ANGLES_INC_HYDROGEN
        self.parm_data['ANGLES_INC_HYDROGEN'] = _zeros(len(self.angles_inc_h)*4)
        for i, ang in enumerate(self.angles_inc_h):
            if -1 in (ang.atom1.idx, ang.atom2.idx, ang.atom3.idx):
                continue
            if ang.angle_type.idx == -1:
                ang.angle_type.idx = angle_type_num
                angle_type_num += 1
            ang.write_info(self, 'ANGLES_INC_HYDROGEN', angle_num)
            angle_num += 1
        self.parm_data['POINTERS'][NTHETH] = angle_num
        self._truncate_array('ANGLES_INC_HYDROGEN', 4*angle_num)
        # Time for Angles without H
        angle_num = 0
        self.parm_data['ANGLES_WITHOUT_HYDROGEN'] = \
                                          _zeros(len(self.angles_without_h)*4)
        for i, ang in enumerate(self.angles_without_h):
            if -1 in (ang.atom1.idx, ang.atom2.idx, ang.atom3.idx):
                continue
            if ang.angle_type.idx == -1:
                ang.angle_type.idx = angle_type_num
                angle_type_num += 1
            ang.write_info(self, 'ANGLES_WITHOUT_HYDROGEN', angle_num)
            angle_num += 1
        # Make sure BOND_FORCE_CONSTANT and BOND_EQUIL_VALUE arrays are big
        # enough
        self.parm_data['ANGLE_FORCE_CONSTANT'] = _zeros(angle_type_num)
        self.parm_data['ANGLE_EQUIL_VALUE'] = _zeros(angle_type_num)
        # Write angle type info to parm
        self.angle_type_list.write_to_parm()
        self.parm_data['POINTERS'][NTHETA] = angle_num
        self.parm_data['POINTERS'][MTHETA] = angle_num
        self.parm_data['POINTERS'][NUMANG] = angle_type_num
        self._truncate_array('ANGLES_WITHOUT_HYDROGEN', 4*angle_num)

        # Now do all of the dihedral arrays
        dihedral_num = 0
        dihedral_type_num = 0
        self.parm_data['DIHEDRALS_INC_HYDROGEN'] = \
                                          _zeros(len(self.dihedrals_inc_h)*5)
        for i, dih in enumerate(self.dihedrals_inc_h):
            if -1 in (dih.atom1.idx, dih.atom2.idx,
                      dih.atom3.idx, dih.atom4.idx):
                continue
            if dih.dihed_type.idx == -1:
                dih.dihed_type.idx = dihedral_type_num
                dihedral_type_num += 1
            dih.write_info(self, 'DIHEDRALS_INC_HYDROGEN', dihedral_num)
            dihedral_num += 1
        self.parm_data['POINTERS'][NPHIH] = dihedral_num
        self._truncate_array('DIHEDRALS_INC_HYDROGEN', 5*dihedral_num)
        # Time for dihedrals without H
        dihedral_num = 0
        self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'] = \
                                       _zeros(len(self.dihedrals_without_h)*5)
        for i, dih in enumerate(self.dihedrals_without_h):
            if -1 in (dih.atom1.idx, dih.atom2.idx,
                      dih.atom3.idx, dih.atom4.idx):
                continue
            if dih.dihed_type.idx == -1:
                dih.dihed_type.idx = dihedral_type_num
                dihedral_type_num += 1
            dih.write_info(self, 'DIHEDRALS_WITHOUT_HYDROGEN', dihedral_num)
            dihedral_num += 1
        self.parm_data['POINTERS'][NPHIA] = dihedral_num
        self.parm_data['POINTERS'][MPHIA] = dihedral_num
        self.parm_data['POINTERS'][NPTRA] = dihedral_type_num
        self._truncate_array('DIHEDRALS_WITHOUT_HYDROGEN', 5*dihedral_num)

        # Adjust lengths of the dihedral arrays to make sure they're long enough
        for key in ('DIHEDRAL_FORCE_CONSTANT', 'DIHEDRAL_PERIODICITY',
                    'DIHEDRAL_PHASE', 'SCEE_SCALE_FACTOR', 'SCNB_SCALE_FACTOR'):
            if not key in self.parm_data.keys(): continue
            self.parm_data[key] = _zeros(dihedral_type_num)
      
        self.dihedral_type_list.write_to_parm()
      
        # Load the pointers now
        self.LoadPointers()
        # Mark atom list as unchanged
        self.atom_list.changed = False
        self.bond_type_list.changed = False
        self.bonds_inc_h.changed = False
        self.bonds_without_h.changed = False
        self.angle_type_list.changed = False
        self.angles_inc_h.changed = False
        self.angles_without_h.changed = False
        self.dihedral_type_list.changed = False
        self.dihedrals_inc_h.changed = False
        self.dihedrals_without_h.changed = False

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
    def _topology_changed(self):
        """ 
        Determines if any of the topological arrays have changed since the
        last upload
        """
        topology_changed = (self.atom_list.changed or
                            self.bond_type_list.changed or
                            self.bonds_inc_h.changed or
                            self.bonds_without_h.changed or
                            self.angle_type_list.changed or
                            self.angles_inc_h.changed or
                            self.angles_without_h.changed or
                            self.dihedral_type_list.changed or
                            self.dihedrals_inc_h.changed or
                            self.dihedrals_without_h.changed)
        if topology_changed and hasattr(self, '_topology'):
            del self._topology
        return topology_changed

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def delete_mask(self, mask):
        """ Deletes all of the atoms corresponding to an entire mask """
        from chemistry.amber.mask import AmberMask
        # Determine if we were given an AmberMask object or a string. If the
        # latter, turn it into an AmberMask and get the selection
        if isinstance(mask, AmberMask):
            # Make sure the AmberMask's parm is this one!
            if id(self) != id(mask.parm):
                raise AmberParmError('Mask belongs to different prmtop!')
            selection = mask.Selection()
        else:
            selection = AmberMask(self, mask).Selection()

        # Delete all of the atoms
        for i in reversed(range(len(selection))):
            if not selection[i]: continue
            self.atom_list[i].deleted = True
            del self.atom_list[i]

        # Remake the topology file and re-set the molecules if we have periodic
        # boxes (or delete the Molecule info if we removed all solvent)
        self.remake_parm()

        # Reconstruct the coordinates and velocities from the remaining atoms
        if hasattr(self, 'coords'):
            self.coords = []
            if self.hasvels: self.vels = []
            for atm in self.atom_list:
                self.coords.extend([atm.xx, atm.xy, atm.xz])
                if self.hasvels: self.vels.extend([atm.vx, atm.vy, atm.vz])

        self._load_structure()
        if self.ptr('ifbox'): self.rediscover_molecules()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def rediscover_molecules(self, solute_ions=True, fix_broken=True):
        """
        This determines the molecularity and sets the ATOMS_PER_MOLECULE and
        SOLVENT_POINTERS sections of the prmtops. Returns the new atom sequence
        in terms of the 'old' atom indexes if re-ordering was necessary to fix
        the tleap bug. Returns None otherwise.
        """
        # Bail out of we are not doing a solvated prmtop
        if not self.ptr('ifbox'): return None

        owner = set_molecules(self)
        ions = ['Br-','Cl-','Cs+','F-','I-','K+','Li+','Mg+','Na+','Rb+','IB',
                'CIO','MG2']
        indices = []
        for res in self.solvent_residues:
            try:
                indices.append(self.parm_data['RESIDUE_LABEL'].index(res))
            except ValueError:
                pass
        # Add ions to list of solvent if necessary
        if not solute_ions:
            for ion in ions:
                if ion in self.parm_data['RESIDUE_LABEL']:
                    indices.append(self.parm_data['RESIDUE_LABEL'].index(ion))
        # If we have no water, we do not have a molecules section!
        if not indices:
            self.parm_data['POINTERS'][IFBOX] = 0
            self.LoadPointers()
            self.deleteFlag('SOLVENT_POINTERS')
            self.deleteFlag('ATOMS_PER_MOLECULE')
            self.deleteFlag('BOX_DIMENSIONS')
            self.hasbox = False
            try: 
                self.rst7.hasbox = False
                del self.box
                del self.rst7.box
            except AttributeError:
                # So we don't have box information... doesn't matter :)
                pass
            return None
        # Now remake our SOLVENT_POINTERS and ATOMS_PER_MOLECULE section
        self.parm_data['SOLVENT_POINTERS'] = [min(indices), len(owner), 0]
        first_solvent = self.parm_data['RESIDUE_POINTER'][min(indices)]
        # Find the first solvent molecule
        for i, mol in enumerate(owner):
            if first_solvent-1 == mol[0]:
                self.parm_data['SOLVENT_POINTERS'][2] = i + 1
                break
        else: # this else belongs to 'for', not 'if'
            raise MoleculeError('Could not find first solvent atom!')

        # Now set up ATOMS_PER_MOLECULE and catch any errors
        self.parm_data['ATOMS_PER_MOLECULE'] = [len(mol) for mol in owner]

        # Check that all of our molecules are contiguous, because we have to
        # re-order atoms if they're not
        try:
            for mol in owner:
                for i in xrange(1, len(mol)):
                    if mol[i] != mol[i-1] + 1:
                        raise StopIteration()
        except StopIteration:
            if not fix_broken:
                raise MoleculeError('Molecule atoms are not contiguous!')
            # Non-contiguous molecules detected... time to fix (ugh!)
            warn('Molecule atoms are not contiguous! I am attempting to fix '
                 'this, but it may take a while.', MoleculeWarning)
            new_atoms = AtomList(self, fill_from=self.atom_list)
            i = 0
            for mol in owner:
                for atm in mol:
                    new_atoms[i] = self.atom_list[atm]
                    i += 1
            self.atom_list = new_atoms
            return owner

        return None

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def writeOFF(self, off_file='off.lib'):
        """ Writes an OFF file from all of the residues found in a prmtop """
        from chemistry.amber.residue import ToResidue
   
        off_file = open(off_file, 'w')
   
        # keep track of all the residues we have to print to the OFF file
        residues = []
   
        # First create a Molecule object from the prmtop
        mol = self.ToMolecule()
   
        # Now loop through all of the residues in the Molecule object and add
        # unique ones to the list of residues to print
        for i in xrange(len(mol.residues)):
            res = ToResidue(mol, i)
            present = False
            for compres in residues:
                if res == compres:
                    present = True
   
            if not present:
                residues.append(res)
      
        # Now that we have all of the residues that we need to add, put their
        # names in the header of the OFF file
        off_file.write('!!index array str\n')
        for res in residues:
            off_file.write(' "%s"\n' % res.name)
   
        # Now write the OFF strings to the file
        for res in residues:
            off_file.write(res.OFF())

        off_file.close()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def fill_LJ(self):
        """
        Fills the LJ_radius, LJ_depth arrays and LJ_types dictionary with data
        from LENNARD_JONES_ACOEF and LENNARD_JONES_BCOEF sections of the prmtop
        files, by undoing the canonical combining rules.
        """
        self.LJ_radius = []  # empty LJ_radii so it can be re-filled
        self.LJ_depth = []   # empty LJ_depths so it can be re-filled
        self.LJ_types = {}   # empty LJ_types so it can be re-filled
        one_sixth = 1 / 6    # we need to raise some numbers to the 1/6th power

        pd = self.parm_data
        acoef = pd['LENNARD_JONES_ACOEF']
        bcoef = pd['LENNARD_JONES_BCOEF']
        natom = self.pointers['NATOM']
        ntypes = self.pointers['NTYPES']
        for i in xrange(natom): # fill the LJ_types array
            self.LJ_types[pd["AMBER_ATOM_TYPE"][i]] = pd["ATOM_TYPE_INDEX"][i]
         
        for i in xrange(ntypes):
            lj_index = pd["NONBONDED_PARM_INDEX"][ntypes*i+i] - 1
            if pd["LENNARD_JONES_ACOEF"][lj_index] < 1.0e-10:
                self.LJ_radius.append(0)
                self.LJ_depth.append(0)
            else:
                factor = 2 * acoef[lj_index] / bcoef[lj_index]
                self.LJ_radius.append(pow(factor, one_sixth) * 0.5)
                self.LJ_depth.append(bcoef[lj_index] / 2 / factor)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def fill_14_LJ(self):
        """
        Fills the LJ_14_radius, LJ_14_depth arrays with data (LJ_types is
        identical) from LENNARD_JONES_14_ACOEF and LENNARD_JONES_14_BCOEF
        sections of the prmtop files, by undoing the canonical combining rules.
        """
        if not self.chamber:
            raise TypeError('fill_14_LJ() only valid on a chamber prmtop!')

        pd = self.parm_data
        acoef = pd['LENNARD_JONES_14_ACOEF']
        bcoef = pd['LENNARD_JONES_14_BCOEF']
        ntypes = self.pointers['NTYPES']

        self.LJ_14_radius = []  # empty LJ_radii so it can be re-filled
        self.LJ_14_depth = []   # empty LJ_depths so it can be re-filled
        one_sixth = 1.0 / 6.0 # we need to raise some numbers to the 1/6th power

        for i in xrange(ntypes):
            lj_index = pd["NONBONDED_PARM_INDEX"][ntypes*i+i] - 1
            if acoef[lj_index] < 1.0e-6:
                self.LJ_14_radius.append(0)
                self.LJ_14_depth.append(0)
            else:
                factor = 2 * acoef[lj_index] / bcoef[lj_index]
                self.LJ_14_radius.append(pow(factor, one_sixth) * 0.5)
                self.LJ_14_depth.append(bcoef[lj_index] / 2 / factor)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def recalculate_LJ(self):
        """
        Takes the values of the LJ_radius and LJ_depth arrays and recalculates
        the LENNARD_JONES_A/BCOEF topology sections from the canonical combining
        rules.
        """
        pd = self.parm_data
        ntypes = self.pointers['NYTPES']
        for i in xrange(ntypes):
            for j in xrange(i, ntypes):
                index = pd['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                rij = self.combine_rmin(self.LJ_radius[i], self.LJ_radius[j])
                wdij = self.combine_epsilon(self.LJ_depth[i], self.LJ_depth[j])
                pd["LENNARD_JONES_ACOEF"][index] = wdij * rij**12
                pd["LENNARD_JONES_BCOEF"][index] = 2 * wdij * rij**6

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def recalculate_14_LJ(self):
        """
        Takes the values of the LJ_radius and LJ_depth arrays and recalculates
        the LENNARD_JONES_A/BCOEF topology sections from the canonical combining
        rules for the 1-4 LJ interactions (CHAMBER only)
        """
        if not self.chamber:
            raise TypeError('recalculate_14_LJ() requires a CHAMBER prmtop!')

        pd = self.parm_data
        ntypes = self.pointers['NYTPES']
        for i in xrange(ntypes):
            for j in xrange(i, ntypes):
                index = pd['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                rij = self.combine_rmin(
                        self.LJ_14_radius[i],self.LJ_14_radius[j]
                )
                wdij = self.combine_epsilon(
                        self.LJ_14_depth[i],self.LJ_14_depth[j]
                )
                pd["LENNARD_JONES_14_ACOEF"][index] = wdij * rij**12
                pd["LENNARD_JONES_14_BCOEF"][index] = 2 * wdij * rij**6

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def combine_rmin(self, rmin1, rmin2):
        """ Define the combining rule for Rmin """
        return rmin1 + rmin2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def combine_epsilon(self, eps1, eps2):
        """ Define the combining rule for Epsilon """
        return sqrt(eps1 * eps2)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def LoadRst7(self, rst7):
        """ Loads coordinates into the AmberParm class """
        if isinstance(rst7, Rst7):
            self.rst7 = rst7
        elif isinstance(rst7, basestring):
            self.rst7 = Rst7.open(rst7)
        self.load_coordinates(self.rst7.coordinates)
        self.hasvels = self.rst7.hasvels
        self.hasbox = self.rst7.hasbox
        if self.hasbox:
            self.box = self.rst7.box
        if self.hasvels:
            self.load_velocities(self.rst7.velocities)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def load_coordinates(self, coords):
        """ Loads the coordinates into the atom list """
        self.coords = coords
        for i, atom in enumerate(self.atom_list):
            i3 = 3 * i
            atom.xx = coords[i3  ]
            atom.xy = coords[i3+1]
            atom.xz = coords[i3+2]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def load_velocities(self, vels):
        """ Loads the coordinates into the atom list """
        self.hasvels = True
        self.vels = vels
        for i, atom in enumerate(self.atom_list):
            i3 = 3 * i
            atom.vx = vels[i3  ]
            atom.vy = vels[i3+1]
            atom.vz = vels[i3+2]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # For backwards-compatibility. AmberParm used to be a single container for
    # chamber and amber topologies, with a property 'chamber' that indicated
    # whether or not it was a chamber-style topology file or not. Now, this is
    # never a chamber topology

    @property
    def chamber(self):
        return False

    @property
    def amoeba(self):
        return False

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def _truncate_array(self, section, length):
        """ Truncates an array to get the given length """
        self.parm_data[section] = self.parm_data[section][:length]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def ToMolecule(self):
        """ Translates an amber system into a molecule format """
        from chemistry.molecule import Molecule

        # Remake the topology file if it's changed
        if self._topology_changed():
            self.remake_parm()
            if self.ptr('ifbox'): self.rediscover_molecules()
            self._load_structure()

        all_bonds = []        # bond array in Molecule format
        residue_pointers = [] # residue pointers adjusted for indexing from 0
        radii = []

        # Set up initial, blank, bond array
        for i in xrange(self.pointers['NATOM']):
            all_bonds.append([])
      
        # Fill up bond arrays with bond partners excluding H atoms
        for i in xrange(self.pointers['MBONA']):
            atom1 = self.parm_data['BONDS_WITHOUT_HYDROGEN'][3*i  ]//3
            atom2 = self.parm_data['BONDS_WITHOUT_HYDROGEN'][3*i+1]//3
            all_bonds[atom1].append(atom2)
            all_bonds[atom2].append(atom1)

        # Fill up bond arrays with bond partners including H atoms
        for i in xrange(self.pointers['NBONH']):
            atom1 = self.parm_data['BONDS_INC_HYDROGEN'][3*i  ]//3
            atom2 = self.parm_data['BONDS_INC_HYDROGEN'][3*i+1]//3
            all_bonds[atom1].append(atom2)
            all_bonds[atom2].append(atom1)

        # Sort bond arrays
        for i in xrange(len(all_bonds)):
            all_bonds[i].sort()

        # Adjust RESIDUE_POINTER for indexing from 0
        for i in xrange(len(self.parm_data['RESIDUE_POINTER'])):
            residue_pointers.append(self.parm_data['RESIDUE_POINTER'][i]-1)

        # Determine which element each atom is
        elements = [periodic_table.Element[atm.atomic_number]
                    for atm in self.atom_list]

        # Put together the title
        title = ''
        try:
            for i in xrange(len(self.parm_data['TITLE'])):
                title += self.parm_data['TITLE'][i]
        except KeyError:
            for i in xrange(len(self.parm_data['CTITLE'])):
                title += self.parm_data['CTITLE'][i]

        # Fill the VDW radii array
        self.fill_LJ()
        for atm in self.atom_list:
            radii.append(self.LJ_radius[self.LJ_types[atm.attype]-1])
        try:
            return Molecule(atoms=self.parm_data['ATOM_NAME'][:],
                            atom_types=self.parm_data['AMBER_ATOM_TYPE'][:],
                            charges=self.parm_data['CHARGE'][:],
                            residues=self.parm_data['RESIDUE_LABEL'][:],
                            bonds=all_bonds,
                            residue_pointers=residue_pointers,
                            coords=self.coords[:],
                            elements=elements,
                            title=title,
                            radii=radii
            )
        except AttributeError: # use dummy list if no coords are loaded
            return Molecule(atoms=self.parm_data['ATOM_NAME'][:],
                            atom_types=self.parm_data['AMBER_ATOM_TYPE'][:],
                            charges=self.parm_data['CHARGE'][:],
                            residues=self.parm_data['RESIDUE_LABEL'][:], 
                            bonds=all_bonds,
                            residue_pointers=residue_pointers,
                            coords=list(xrange(self.pointers['NATOM']*3)),
                            elements=elements,
                            title=title,
                            radii=radii
            )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Rst7(object):
    """
    Amber input coordinate (or restart coordinate) file. Front-end for the
    readers and writers, supports both NetCDF and ASCII restarts.
    """

    def __init__(self, filename=None, natom=None, title='', hasvels=False,
                 hasbox=False, time=0.0):
        """
        Optionally takes a filename to read. This is deprecated, though, as the
        alternative constructor "open" should be used instead
        """
        self.coordinates = []
        self.velocities = []
        self.box = []
        self.hasvels = hasvels
        self.hasbox = hasbox
        self.natom = natom
        self.title = title
        self.time = 0
        if filename is not None:
            self.filename = filename
            warn('Use Rst7.open() constructor instead of default constructor '
                 'to parse restart files.', DeprecationWarning)
            self._read(filename)

    @classmethod
    def open(cls, filename):
        """
        Constructor that opens and parses an input coordinate file
        """
        inst = cls()
        inst.filename = filename
        inst._read(filename)
        return inst

    def _read(self, filename):
        """
        Open and parse an input coordinate file in either ASCII or NetCDF format
        """
        from chemistry.amber.asciicrd import AmberAsciiRestart
        from chemistry.amber.netcdffiles import NetCDFRestart
        try:
            f = AmberAsciiRestart(filename, 'r')
            self.natom = f.natom
        except ValueError:
            # Maybe it's a NetCDF file?
            try:
                f = NetCDFRestart.open_old(filename)
                self.natom = f.atom
            except ImportError:
                raise ReadError('Could not parse %s as an ASCII restart and '
                                'could not find any NetCDF-Python packages to '
                                'attempt to parse as a NetCDF Restart.'
                                % filename)
            except RuntimeError:
                raise ReadError('Could not parse restart file %s' % filename)

        self.coordinates = f.coordinates
        self.hasvels = f.hasvels
        self.hasbox = f.hasbox
        if f.hasvels:
            self.velocities = f.velocities
        if f.hasbox:
            self.box = f.box
        self.title = f.title
        self.time = f.time

    @property
    def coords(self):
        """ Deprecated for coordinates now """
        warn('coords attribute of Rst7 is deprecated. Use coordinates instead',
             DeprecationWarning)
        return self.coordinates
   
    @property
    def vels(self):
        """ Deprecated for velocities now """
        warn('vels attribute of Rst7 is deprecated. Use velocities instead',
             DeprecationWarning)
        return self.velocities

    @classmethod
    def copy_from(cls, thing):
        """
        Copies the coordinates, velocities, and box information from another
        instance
        """
        inst = cls()
        inst.natom = thing.natom
        inst.title = thing.title
        inst.coordinates = thing.coordinates[:]
        inst.hasvels = thing.hasvels
        if hasattr(thing, 'velocities'): inst.velocities = thing.velocities[:]
        inst.hasbox = thing.hasbox
        if hasattr(thing, 'box'): inst.box = thing.box[:]
        inst.time = thing.time

        return inst

    def __copy__(self):
        """ Copy constructor """
        return type(self).copy_from(self)

    def write(self, fname, netcdf=False):
        """ Writes the coordinates and/or velocities to a restart file """
        from chemistry.amber.asciicrd import AmberAsciiRestart
        from chemistry.amber.netcdffiles import NetCDFRestart
        if netcdf:
            if self.natom is None:
                raise RuntimeError('Number of atoms must be set for NetCDF '
                                   'Restart files before write time')
            f = NetCDFRestart.open_new(fname, self.natom, self.hasbox,
                                       self.hasvels, self.title)
        else:
            f = AmberAsciiRestart(fname, 'w', natom=self.natom,
                                  title=self.title)

        f.time = self.time
        # Now write the coordinates
        f.coordinates = self.coordinates
        if self.hasvels:
            f.velocities = self.velocities
        if self.hasbox:
            f.box = self.box
        f.close()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def set_molecules(parm):
    """
    Correctly sets the ATOMS_PER_MOLECULE and SOLVENT_POINTERS sections of the
    topology file.
    """
    from sys import setrecursionlimit, getrecursionlimit
    # Since we use a recursive function here, we make sure that the recursion
    # limit is large enough to handle the maximum possible recursion depth we'll
    # need (NATOM). We don't want to shrink it, though, since we use list
    # comprehensions in list constructors in some places that have an implicit
    # (shallow) recursion, therefore, reducing the recursion limit too much here
    # could raise a recursion depth exceeded exception during a _Type/Atom/XList
    # creation. Therefore, set the recursion limit to the greater of the current
    # limit or the number of atoms
    setrecursionlimit(max(parm.ptr('natom'), getrecursionlimit()))

    # Unmark all atoms so we can track which molecule each goes into
    parm.atom_list.unmark()

    if not parm.ptr('ifbox'):
        raise MoleculeError('Only periodic prmtops can have '
                            'Molecule definitions')
    # The molecule "ownership" list
    owner = []
    # The way I do this is via a recursive algorithm, in which
    # the "set_owner" method is called for each bonded partner an atom
    # has, which in turn calls set_owner for each of its partners and 
    # so on until everything has been assigned.
    molecule_number = 1 # which molecule number we are on
    for i in xrange(parm.ptr('natom')):
        # If this atom has not yet been "owned", make it the next molecule
        # However, we only increment which molecule number we're on if 
        # we actually assigned a new molecule (obviously)
        if not parm.atom_list[i].marked:
            tmp = [i]
            _set_owner(parm, tmp, i, molecule_number)
            # Make sure the atom indexes are sorted
            tmp.sort()
            owner.append(tmp)
            molecule_number += 1
    return owner

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _set_owner(parm, owner_array, atm, mol_id):
    """ Recursively sets ownership of given atom and all bonded partners """
    parm.atom_list[atm].marked = mol_id
    for partner in parm.atom_list[atm].bond_partners:
        if not partner.marked:
            owner_array.append(partner.starting_index)
            _set_owner(parm, owner_array, partner.starting_index, mol_id)
        elif partner.marked != mol_id:
            raise MoleculeError('Atom %d in multiple molecules' % 
                                partner.starting_index)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _zeros(length):
    """ Returns an array of zeros of the given length """
    return [0 for i in xrange(length)]
