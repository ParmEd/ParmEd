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

from chemistry.amber._amberparm import AmberParm, Rst7
from chemistry.amber.constants import (NATOM, NTYPES, NBONH, MBONA, NTHETH,
                MTHETA, NPHIH, MPHIA, NNB, NRES, NBONA, NTHETA, NPHIA, NUMBND,
                NUMANG, NPTRA, NATYP, IFBOX, NMXRS, CHARMM_ELECTROSTATIC)
from chemistry.topologyobjects import (UreyBradley, Improper, Cmap, BondType,
                                       ImproperType, CmapType)
from chemistry.exceptions import CharmmPSFError
from math import pi, sqrt
import warnings

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ChamberParm(AmberParm):
    """
    Chamber Topology (parm7 format) class. Gives low, and some high, level
    access to topology data or interact with some of the high-level classes
    comprising the system topology and parameters.  The ChamberParm class uses
    the same attributes that the AmberParm class uses, and only the ones unique
    to ChamberParm will be shown below.

    Parameters
    ----------
    prm_name : str=None
        If provided, this file is parsed and the data structures will be loaded
        from the data in this file
    rst7_name : str=None
        If provided, the coordinates and unit cell dimensions from the provided
        Amber inpcrd/restart file will be loaded into the molecule

    Attributes
    ----------
    In addition to the attributes listed for `AmberParm`, the following
    attributes are available:

    LJ_14_radius : list(float)
        The same as LJ_radius, except specific for 1-4 nonbonded parameters,
        which may differ in the CHARMM force field
    LJ_14_depth : list(float)
        The same as LJ_depth, except specific for 1-4 nonbonded parameters,
        which may differ in the CHARMM force field
    urey_bradleys : TrackedList(UreyBradley)
        List of Urey-Bradley terms between two atoms in a valence angle
    impropers : TrackedList(Improper)
        List of CHARMM-style improper torsions
    cmaps : TrackedList(Cmap)
        List of coupled-torsion correction map parameters
    urey_bradley_types : TrackedList(UreyBradleyType)
        List of parameters defining the Urey-Bradley terms
    improper_types : TrackedList(Improper)
        List of parameters defining the Improper terms
    cmap_types : TrackedList(CmapType)
        List of parameters defining the CMAP terms
    chamber : bool=True
        On ChamberParm instances, this is always True to indicate that it is a
        CHAMBER-style topology file
    amoeba : bool=False
        On ChamberParm instances, this is always False to indicate that it is
        not an AMOEBA-style topology file
    has_cmap : bool
        True if CMAP parameters are present in this system; False otherwise
    """

    CHARGE_SCALE = CHARMM_ELECTROSTATIC

    #===================================================

    def initialize_topology(self, rst7_name=None):
        """
        Initializes topology data structures, like the list of atoms, bonds,
        etc., after the topology file has been read. The following methods are
        called:
        """
        self.LJ_14_radius = []
        self.LJ_14_depth = []
        AmberParm.initialize_topology(self, rst7_name)
      
    #===================================================

    def __copy__(self):
        """ Needs to copy a few additional data structures """
        other = AmberParm.__copy__(self)
        other.LJ_14_radius = self.LJ_14_radius[:]
        other.LJ_14_depth = self.LJ_14_depth[:]
        return other

    #===================================================
   
    def load_pointers(self):
        """
        Loads the data in POINTERS section into a pointers dictionary with each
        key being the pointer name according to http://ambermd.org/formats.html
        """
        AmberParm.load_pointers(self)
        # Other pointers
        nub, nubtypes = self.parm_data['CHARMM_UREY_BRADLEY_COUNT'][:2]
        self.pointers['NUB'] = nub
        self.pointers['NUBTYPES'] = nubtypes
        self.pointers['NIMPHI'] = self.parm_data['CHARMM_NUM_IMPROPERS'][0]
        self.pointers['NIMPRTYPES']=self.parm_data['CHARMM_NUM_IMPR_TYPES'][0]
        # If CMAP is not present, don't load the pointers
        if self.has_cmap:
            self.pointers['CMAP'] = self.parm_data['CHARMM_CMAP_COUNT'][0]
            self.pointers['CMAP_TYPES'] = self.parm_data['CHARMM_CMAP_COUNT'][1]

    #===================================================

    def load_structure(self):
        """ 
        Loads all of the topology instance variables. This is necessary if we
        actually want to modify the topological layout of our system
        (like deleting atoms)
        """
        super(ChamberParm, self).load_structure()
        self._load_urey_brad_info()
        self._load_improper_info()
        self._load_cmap_info()
        super(ChamberParm, self).unchange()

    #===================================================

    def _load_urey_brad_info(self):
        """ Loads the Urey-Bradley types and array """
        del self.urey_bradleys[:]
        del self.urey_bradley_types[:]
        for k, req in zip(self.parm_data['CHARMM_UREY_BRADLEY_FORCE_CONSTANT'],
                          self.parm_data['CHARMM_UREY_BRADLEY_EQUIL_VALUE']):
            self.urey_bradley_types.append(
                    BondType(k, req, self.urey_bradley_types)
            )
        ulist = self.parm_data['CHARMM_UREY_BRADLEY']
        for i in xrange(0, 3*self.pointers['NUB'], 3):
            self.urey_bradleys.append(
                    UreyBradley(self.atoms[ulist[i  ]-1],
                                self.atoms[ulist[i+1]-1],
                                self.urey_bradley_types[ulist[i+2]-1])
            )

    #===================================================

    def _load_improper_info(self):
        """ Loads the CHARMM Improper types and array """
        del self.impropers[:]
        del self.improper_types[:]
        for k, eq in zip(self.parm_data['CHARMM_IMPROPER_FORCE_CONSTANT'],
                         self.parm_data['CHARMM_IMPROPER_PHASE']):
            self.improper_types.append(
                    ImproperType(k, eq, self.improper_types)
            )
        ilist = self.parm_data['CHARMM_IMPROPERS']
        for i in xrange(0, 5*self.pointers['NIMPHI'], 5):
            self.impropers.append(
                    Improper(self.atoms[ilist[i  ]-1],
                             self.atoms[ilist[i+1]-1],
                             self.atoms[ilist[i+2]-1],
                             self.atoms[ilist[i+3]-1],
                             self.improper_types[ilist[i+4]-1])
            )

    #===================================================

    def _load_cmap_info(self):
        """ Loads the CHARMM CMAP types and array """
        if not self.has_cmap: return
        del self.cmaps[:]
        del self.cmap_types[:]
        for i in xrange(self.pointers['CMAP_TYPES']):
            resolution = self.parm_data['CHARMM_CMAP_RESOLUTION'][i]
            grid = self.parm_data['CHARMM_CMAP_PARAMETER_%02d' % (i+1)]
            cmts = self.parm_comments['CHARMM_CMAP_PARAMETER_%02d' % (i+1)]
            self.cmap_types.append(
                    CmapType(resolution, grid, cmts, list=self.cmap_types)
            )
        clist = self.parm_data['CHARMM_CMAP_INDEX']
        for i in xrange(0, 6*self.pointers['CMAP'], 6):
            self.cmaps.append(
                    Cmap(self.atoms[clist[i  ]-1], self.atoms[clist[i+1]-1],
                         self.atoms[clist[i+2]-1], self.atoms[clist[i+3]-1],
                         self.atoms[clist[i+4]-1],
                         self.cmap_types[clist[i+5]-1])
            )

    #===================================================

    def remake_parm(self):
        """
        Re-fills the topology file arrays if we have changed the underlying
        structure
        """
        super(ChamberParm, self).remake_parm()
        self._xfer_urey_bradley_properties()
        self._xfer_improper_properties()
        self._xfer_cmap_properties()
        super(ChamberParm, self).unchange()

    #===================================================

    def writeOFF(self, off_file='off.lib'):
        """ Writes an OFF file from all of the residues found in a prmtop """
        warnings.warn('OFF files from chamber topologies should NOT be used '
                      'in LEaP unless you KNOW what you are doing...')
        return super(ChamberParm, self).writeOFF(off_file)

    #===================================================

    def fill_LJ(self):
        """
        Fills the LJ_radius, LJ_depth arrays and LJ_types dictionary with data
        from LENNARD_JONES_ACOEF and LENNARD_JONES_BCOEF sections of the prmtop
        files, by undoing the canonical combining rules.
        """
        AmberParm.fill_LJ(self)

        acoef = self.parm_data['LENNARD_JONES_14_ACOEF']
        bcoef = self.parm_data['LENNARD_JONES_14_BCOEF']
        ntypes = self.pointers['NTYPES']

        self.LJ_14_radius = []  # empty LJ_radii so it can be re-filled
        self.LJ_14_depth = []   # empty LJ_depths so it can be re-filled
        one_sixth = 1.0 / 6.0 # we need to raise some numbers to the 1/6th power

        for i in xrange(ntypes):
            lj_index = self.parm_data["NONBONDED_PARM_INDEX"][ntypes*i+i] - 1
            if acoef[lj_index] < 1.0e-6:
                self.LJ_14_radius.append(0)
                self.LJ_14_depth.append(0)
            else:
                factor = 2 * acoef[lj_index] / bcoef[lj_index]
                self.LJ_14_radius.append(pow(factor, one_sixth) * 0.5)
                self.LJ_14_depth.append(bcoef[lj_index] / 2 / factor)

    #===================================================

    def recalculate_LJ(self):
        """
        Takes the values of the LJ_radius and LJ_depth arrays and recalculates
        the LENNARD_JONES_A/BCOEF topology sections from the canonical
        combining rules.
        """
        AmberParm.recalculate_LJ(self)
        ntypes = self.pointers['NTYPES']
        acoef = self.parm_data['LENNARD_JONES_14_ACOEF']
        bcoef = self.parm_data['LENNARD_JONES_14_BCOEF']
        for i in xrange(ntypes):
            for j in xrange(i, ntypes):
                index = self.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                rij = self.combine_rmin(self.LJ_14_radius[i],
                                        self.LJ_14_radius[j])
                wdij = self.combine_epsilon(self.LJ_14_depth[i],
                                            self.LJ_14_depth[j])
                acoef[index] = wdij * rij ** 12
                bcoef[index] = 2 * wdij * rij**6

    #===================================================

    @property
    def chamber(self):
        return True
   
    @property
    def amoeba(self):
        return False

    @property
    def has_cmap(self):
        return 'CHARMM_CMAP_COUNT' in self.flag_list

    #===========  PRIVATE INSTANCE METHODS  ============

    def _xfer_urey_bradley_properties(self):
        """
        Sets the various topology file section data from the Urey-Bradley arrays
        """
        data = self.parm_data
        for urey_type in self.urey_bradley_types:
            urey_type.used = False
        for urey in self.urey_bradleys:
            urey.type.used = True
        self.urey_bradley_types.prune_unused()
        data['CHARMM_UREY_BRADLEY_FORCE_CONSTANT'] = \
                [type.k for type in self.urey_bradley_types]
        data['CHARMM_UREY_BRADLEY_EQUIL_VALUE'] = \
                [type.req for type in self.urey_bradley_types]
        data['CHARMM_UREY_BRADLEY'] = bond_array = []
        for i, urey in enumerate(self.urey_bradleys):
            bond_array.extend([urey.atom1.idx+1, urey.atom2.idx+1,
                               urey.type.idx+1])
        data['CHARMM_UREY_BRADLEY_COUNT'] = [i+1, len(self.urey_bradley_types)]
        self.pointers['NUB'] = i + 1
        self.pointers['NUBTYPES'] = len(self.urey_bradley_types)

    #===================================================

    def _xfer_improper_properties(self):
        """ Sets the topology file section data from the improper arrays """
        data = self.parm_data
        for improper_type in self.improper_types:
            improper_type.used = False
        for improper in self.impropers:
            improper.type.used = True
        self.improper_types.prune_unused()
        data['CHARMM_IMPROPER_FORCE_CONSTANT'] = \
                [type.psi_k for type in self.improper_types]
        data['CHARMM_IMPROPER_PHASE'] = \
                [type.psi_eq for type in self.improper_types]
        data['CHARMM_IMPROPERS'] = improper_array = []
        for i, imp in enumerate(self.impropers):
            improper_array.extend([imp.atom1.idx+1, imp.atom2.idx+1,
                                   imp.atom3.idx+1, imp.atom4.idx+1,
                                   imp.type.idx+1])
        data['CHARMM_NUM_IMPROPERS'] = [i+1]
        data['CHARMM_NUM_IMPR_TYPES'] = [len(self.improper_types)]
        self.pointers['NIMPHI'] = i + 1
        self.pointers['NIMPRTYPES'] = len(self.improper_types)

    #===================================================

    def _xfer_cmap_properties(self):
        """ Sets the topology file section data from the cmap arrays """
        if not self.has_cmap: return
        data = self.parm_data
        for ct in self.cmap_types:
            ct.used = False
        for cmap in self.cmaps:
            cmap.type.used = True
        self.cmap_types.prune_unused()
        # If we have no cmaps, delete all remnants of the CMAP terms in the
        # prmtop and bail out
        if len(self.cmaps) == 0:
            # We have deleted all cmaps. Get rid of them from the parm file and
            # bail out. This is probably pretty unlikely, though...
            self.delete_flag('CHARMM_CMAP_COUNT')
            self.delete_flag('CHARMM_CMAP_RESOLUTION')
            for flag in self.flag_list:
                if flag.startswith('CHARMM_CMAP_PARAMETER'):
                    self.delete_flag(flag)
            del self.pointers['CMAP']
            del self.pointers['CMAP_TYPES']
            return
        # All of our CMAP types are in different topology file sections. We need
        # to delete all of the CHARMM_CMAP_PARAMETER_XX sections and then
        # recreate them with the correct size and comments.  The comments have
        # been stored in the CMAP types themselves to prevent them from being
        # lost. We will also assume that the Fortran format we're going to use
        # is the same for all CMAP types, so just pull it from
        # CHARMM_CMAP_PARAMETER_01 (or fall back to 8(F9.5))
        try:
            fmt = str(self.formats['CHARMM_CMAP_PARAMETER_01'])
        except KeyError:
            fmt = '8(F9.5)'
        for flag in self.flag_list:
            if flag.startswith('CHARMM_CMAP_PARAMETER'):
                self.delete_flag(flag)
        # Now add them back
        for i, ct in enumerate(self.cmap_types):
            self.add_flag('CHARMM_CMAP_PARAMETER_%02d' % (i+1), fmt,
                          data=ct.grid, comments=ct.comments)
        # Now do the CMAP_INDEX section
        data['CHARMM_CMAP_INDEX'] = cmap_array = []
        for i, cm in enumerate(self.cmaps):
            cmap_array.extend([cmap.atom1.idx+1, cmap.atom2.idx+1,
                               cmap.atom3.idx+1, cmap.atom4.idx+1,
                               cmap.atom5.idx+1, cmap.type.idx+1])
        data['CHARMM_CMAP_COUNT'] = [i+1, len(self.cmap_types)]
        self.pointers['CMAP'] = i+1
        self.pointers['CMAP_TYPES'] = len(self.cmap_types)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

charmm_solvent = ('WAT', 'TIP3', 'HOH', 'TIP4', 'TIP5', 'SPCE', 'SPC')

def ConvertFromPSF(struct, frcfield, vmd=False, title=''):
    """
    This function instantiates a ChamberParm instance from a data structure
    instantiated by a CHARMM PSF.

    Parameters:
        - struct (ProteinStructure) : The structure loaded from a PSF file.
        - frcfield (str) : The name of the CHARMM force field being used
        - vmd (bool) : Print a VMD-compatible topology file
        - title (str) : The title of the topology file

    Returns:
        ChamberParm with all parameters loaded
    """
    # Make sure that frcfield is the right type
    if isinstance(frcfield, list):
        for i, x in enumerate(frcfield):
            if i % 2 == 0 and not isinstance(x, int):
                raise TypeError('frcfield must be [int, str, [int, str, ...]]')
            elif i % 2 == 1 and not isinstance(x, str):
                raise TypeError('frcfield must be [int, str, [int, str, ...]]')
    # Make sure the structure has had its box info loaded
    if not hasattr(struct, 'box'):
        raise CharmmPSFError('Box and coordinate info must be set in the PSF '
                             'prior to converting to a topology.')
    struct.atom_list.assign_indexes()
    # First define the pointers
    # Split parameters b/w those with hydrogen and those without
    bonds_inc_h, bonds_without_h = [], []
    for bond in struct.bond_list:
        if (bond.atom1.type.atomic_number == 1 or
            bond.atom2.type.atomic_number == 1):
            bonds_inc_h.append(bond)
        else:
            bonds_without_h.append(bond)
    angles_inc_h, angles_without_h = [], []
    for angle in struct.angle_list:
        if (angle.atom1.type.atomic_number == 1 or
            angle.atom2.type.atomic_number == 1 or
            angle.atom3.type.atomic_number == 1):
            angles_inc_h.append(angle)
        else:
            angles_without_h.append(angle)
    dihedrals_inc_h, dihedrals_without_h = [], []
    for dihed in struct.dihedral_parameter_list:
        if (dihed.atom1.type.atomic_number == 1 or
            dihed.atom2.type.atomic_number == 1 or
            dihed.atom3.type.atomic_number == 1 or
            dihed.atom4.type.atomic_number == 1):
            dihedrals_inc_h.append(dihed)
        else:
            dihedrals_without_h.append(dihed)
    pointers = [0 for i in xrange(32)]
    # Define the exclusion list from the 1-2, 1-3, and 1-4 partners
    excluded_atoms_list = []
    num_excluded_atoms = [0 for i in xrange(len(struct.atom_list))]
    for i, atom in enumerate(struct.atom_list):
        vals_to_add = []
        for partner in atom.bond_partners:
            if partner.idx > atom.idx: vals_to_add.append(partner.idx + 1)
        for partner in atom.angle_partners:
            if partner.idx > atom.idx: vals_to_add.append(partner.idx + 1)
        for partner in atom.dihedral_partners:
            if partner.idx > atom.idx: vals_to_add.append(partner.idx + 1)
        if not vals_to_add: vals_to_add = [0]
        num_excluded_atoms[i] = len(vals_to_add)
        excluded_atoms_list.extend(sorted(vals_to_add))
    # Determine the number of LJ atom types by condensing them into as few atom
    # types as possible
    lj_idx_list = [0 for atom in struct.atom_list]
    num_lj_types = 0
    lj_radii, lj_depths = [], []
    lj_radii14, lj_depths14 = [], []
    lj_type_list = []
    for i, atom in enumerate(struct.atom_list):
        atom = atom.type
        # If we have already assigned the LJ type, skip this one
        if lj_idx_list[i] > 0: continue
        # This is a NEW nonbonded LJ type
        num_lj_types += 1
        lj_idx_list[i] = num_lj_types
        ljtype = (atom.rmin, atom.rmin_14, atom.epsilon, atom.epsilon_14)
        lj_type_list.append(atom)
        lj_radii.append(atom.rmin)
        lj_radii14.append(atom.rmin_14)
        lj_depths.append(atom.epsilon)
        lj_depths14.append(atom.epsilon_14)
        # Look through the rest of the atoms and assign any equivalent types as
        # the same
        for j in xrange(i+1, len(struct.atom_list)):
            atom2 = struct.atom_list[j].type
            if lj_idx_list[j] > 0: continue
            elif atom2 is atom: # same exact type!
                lj_idx_list[j] = num_lj_types
            elif not atom.nbfix:
                # Only non-NBFIXed atom types can be compressed
                ljtype2 = (atom2.rmin, atom2.rmin_14,
                           atom2.epsilon, atom2.epsilon_14)
                if ljtype == ljtype2:
                    lj_idx_list[j] = num_lj_types
    # Now we should be ready to set some of our pointers (not # of types yet)
    pointers[NATOM] = len(struct.atom_list)
    pointers[NTYPES] = pointers[NATYP] = num_lj_types
    pointers[NBONH] = len(bonds_inc_h)
    pointers[MBONA] = pointers[NBONA] = len(bonds_without_h)
    pointers[NTHETH] = len(angles_inc_h)
    pointers[MTHETA] = pointers[NTHETA] = len(angles_without_h)
    pointers[NPHIH] = len(dihedrals_inc_h)
    pointers[MPHIA] = pointers[NPHIA] = len(dihedrals_without_h)
    pointers[NNB] = len(excluded_atoms_list)
    pointers[NRES] = len(struct.residue_list)
    if struct.box is None:
        pointers[IFBOX] = 0
    elif struct.box[3:] == [90, 90, 90]:
        pointers[IFBOX] = 1
    else:
        # triclinic
        pointers[IFBOX] = 2
    pointers[NMXRS] = max([len(r) for r in struct.residue_list])
    # First add all of the flags
    parm = ChamberParm()
    if vmd:
        parm.add_flag('TITLE', '20a4', data=[title])
    else:
        parm.add_flag('CTITLE', 'a80', data=[title])
    parm.add_flag('POINTERS', '10I8', data=pointers)
    parm.add_flag('FORCE_FIELD_TYPE', 'i2,a78', data=frcfield)
    parm.add_flag('ATOM_NAME','20a4',
                  data=[a.name[:4] for a in struct.atom_list])
    chgdat = [a.charge for a in struct.atom_list]
    if vmd:
        parm.add_flag('CHARGE', '5E16.8', data=chgdat)
    else:
        parm.add_flag('CHARGE', '3E24.16', data=chgdat,
                     comments=['Atomic charge multiplied by sqrt(332.0716D0) '
                               '(CCELEC)']
        )
    parm.add_flag('MASS', '5E16.8', data=[a.mass for a in struct.atom_list])
    parm.add_flag('ATOM_TYPE_INDEX', '10I8', data=lj_idx_list)
    parm.add_flag('NUMBER_EXCLUDED_ATOMS', '10I8', data=num_excluded_atoms)
    parm.add_flag('EXCLUDED_ATOMS_LIST', '10I8', data=excluded_atoms_list)
    # Assign the nonbonded parm index table
    holder = [0 for i in xrange(num_lj_types*num_lj_types)]
    idx = 0
    for i in xrange(num_lj_types):
        for j in xrange(i+1):
            idx += 1
            holder[num_lj_types*i+j] = idx
            holder[num_lj_types*j+i] = idx
    holder2 = [0 for i in xrange(num_lj_types*num_lj_types)]
    idx = 0
    for i in xrange(num_lj_types):
        for j in xrange(num_lj_types):
            holder2[idx] = holder[i*num_lj_types+j]
            idx += 1
    del holder
    parm.add_flag('NONBONDED_PARM_INDEX', '10I8', data=holder2)
    parm.add_flag('RESIDUE_LABEL', '20a4',
                  data=[res.resname[:4] for res in struct.residue_list])
    resptr = [0 for i in xrange(len(struct.residue_list))]
    n = 1
    for i, res in enumerate(struct.residue_list):
        resptr[i] = n
        n += len(res)
    parm.add_flag('RESIDUE_POINTER', '10I8', data=resptr)
    # Assign the bond types and bond type indexes
    numbnd = 0
    bond_frc_cnst, bond_equil = [], []
    for bond in struct.bond_list: bond.bond_type.idx = -1
    for bond in struct.bond_list:
        if bond.bond_type.idx == -1:
            bond_frc_cnst.append(bond.bond_type.k)
            bond_equil.append(bond.bond_type.req)
            bond.bond_type.idx = numbnd
            numbnd += 1
    parm.parm_data['POINTERS'][NUMBND] = numbnd
    parm.add_flag('BOND_FORCE_CONSTANT', '5E16.8', data=bond_frc_cnst)
    parm.add_flag('BOND_EQUIL_VALUE', '5E16.8', data=bond_equil)
    # Assign the angle types and angle type indexes
    numang = 0
    angle_frc_cnst, angle_equil = [], []
    for angle in struct.angle_list: angle.angle_type.idx = -1
    for angle in struct.angle_list:
        if angle.angle_type.idx == -1:
            angle_frc_cnst.append(angle.angle_type.k)
            angle_equil.append(angle.angle_type.theteq * pi/180) # to radians
            angle.angle_type.idx = numang
            numang += 1
    parm.parm_data['POINTERS'][NUMANG] = numang
    parm.add_flag('ANGLE_FORCE_CONSTANT', '5E16.8', data=angle_frc_cnst)
    parm.add_flag('ANGLE_EQUIL_VALUE', '3E25.17', data=angle_equil)
    # Assign the Urey-Bradley terms
    nub = 0
    ub_frc_cnst, ub_equil = [], []
    for ureybrad in struct.urey_bradley_list: ureybrad.ub_type.idx = -1
    for ureybrad in struct.urey_bradley_list:
        if ureybrad.ub_type.idx == -1:
            ub_frc_cnst.append(ureybrad.ub_type.k)
            ub_equil.append(ureybrad.ub_type.req)
            ureybrad.ub_type.idx = nub
            nub += 1
    ubcnt = [len(struct.urey_bradley_list), nub]
    parm.add_flag('CHARMM_UREY_BRADLEY_COUNT', '2I8', data=ubcnt,
                  comments=['V(ub) = K_ub(r_ik - R_ub)**2',
                            'Number of Urey Bradley terms and types']
    )
    parm.add_flag('CHARMM_UREY_BRADLEY', '10I8', num_items=ubcnt[0]*3,
                  comments=['List of the two atoms and its parameter index',
                            'in each UB term: i,k,index']
    )
    parm.add_flag('CHARMM_UREY_BRADLEY_FORCE_CONSTANT', '5E16.8',
                  data=ub_frc_cnst, comments=['K_ub: kcal/mol/A**2'])
    parm.add_flag('CHARMM_UREY_BRADLEY_EQUIL_VALUE', '5E16.8', data=ub_equil,
                  comments='r_ub: A')
    for i, ureybrad in enumerate(struct.urey_bradley_list):
        parm.parm_data['CHARMM_UREY_BRADLEY'][3*i  ] = ureybrad.atom1.idx + 1
        parm.parm_data['CHARMM_UREY_BRADLEY'][3*i+1] = ureybrad.atom2.idx + 1
        parm.parm_data['CHARMM_UREY_BRADLEY'][3*i+2] = ureybrad.ub_type.idx + 1
    # Assign the dihedral constants
    nphi = 0
    dih_frc_cnst, dih_per, dih_phase = [], [], []
    for dihed in struct.dihedral_parameter_list: dihed.dihedral_type.idx = -1
    for dihed in struct.dihedral_parameter_list:
        if dihed.dihedral_type.idx == -1:
            dih_frc_cnst.append(dihed.dihedral_type.phi_k)
            dih_per.append(dihed.dihedral_type.per)
            dih_phase.append(dihed.dihedral_type.phase * pi/180)
            dihed.dihedral_type.idx = nphi
            nphi += 1
    parm.parm_data['POINTERS'][NPTRA] = nphi
    parm.add_flag('DIHEDRAL_FORCE_CONSTANT', '5E16.8', data=dih_frc_cnst)
    parm.add_flag('DIHEDRAL_PERIODICITY', '5E16.8', data=dih_per)
    parm.add_flag('DIHEDRAL_PHASE', '5E16.8', data=dih_phase)
    # No 1-4 scaling 
    parm.add_flag('SCEE_SCALE_FACTOR', '5E16.8', data=[1 for i in xrange(nphi)])
    parm.add_flag('SCNB_SCALE_FACTOR', '5E16.8', data=[1 for i in xrange(nphi)])
    # Assign impropers
    nimpt = 0
    imp_frc_cnst, imp_equil = [], []
    for imp in struct.improper_list: imp.improper_type.idx = -1
    for imp in struct.improper_list:
        if imp.improper_type.idx == -1:
            imp_frc_cnst.append(imp.improper_type.k)
            imp_equil.append(imp.improper_type.phieq)
            imp.improper_type.idx = nimpt
            nimpt += 1
    nimp = len(struct.improper_list)
    parm.add_flag('CHARMM_NUM_IMPROPERS', '10I8',
                  data=[nimp], comments=['Number of terms contributing to the',
                  'quadratic four atom improper energy term:',
                  'V(improper) = K_psi(psi - psi_0)**2']
    )
    parm.add_flag('CHARMM_IMPROPERS', '10I8', num_items=nimp*5,
                  comments=['List of the four atoms in each improper term',
                            'i,j,k,l,index  i,j,k,l,index',
                            'where index is into the following two lists:',
                            'CHARMM_IMPROPER_{FORCE_CONSTANT,IMPROPER_PHASE}']
    )
    parm.add_flag('CHARMM_NUM_IMPR_TYPES', '1I8', data=[nimpt],
                  comments=['Number of unique parameters contributing to the',
                            'quadratic four atom improper energy term']
    )
    parm.add_flag('CHARMM_IMPROPER_FORCE_CONSTANT', '5E16.8', data=imp_frc_cnst,
                  comments=['K_psi: kcal/mole/rad**2'])
    parm.add_flag('CHARMM_IMPROPER_PHASE', '5E16.8', data=imp_equil,
                  comments=['psi: degrees'])
    # Add the impropers
    for i, imp in enumerate(struct.improper_list):
        parm.parm_data['CHARMM_IMPROPERS'][5*i  ] = imp.atom1.idx + 1
        parm.parm_data['CHARMM_IMPROPERS'][5*i+1] = imp.atom2.idx + 1
        parm.parm_data['CHARMM_IMPROPERS'][5*i+2] = imp.atom3.idx + 1
        parm.parm_data['CHARMM_IMPROPERS'][5*i+3] = imp.atom4.idx + 1
        parm.parm_data['CHARMM_IMPROPERS'][5*i+4] = imp.improper_type.idx + 1

    parm.add_flag('SOLTY','5E16.8',num_items=parm.parm_data['POINTERS'][NATYP])
    if vmd:
        fmt = '5E16.8'
    else:
        fmt = '3E24.16'
    nttyp = num_lj_types * (num_lj_types + 1) // 2
    parm.add_flag('LENNARD_JONES_ACOEF', fmt, num_items=nttyp)
    parm.add_flag('LENNARD_JONES_BCOEF', fmt, num_items=nttyp)
    parm.add_flag('LENNARD_JONES_14_ACOEF', fmt, num_items=nttyp)
    parm.add_flag('LENNARD_JONES_14_BCOEF', fmt, num_items=nttyp)
    for i in xrange(num_lj_types):
        for j in xrange(i, num_lj_types):
            index = parm.parm_data['NONBONDED_PARM_INDEX'][num_lj_types*i+j] - 1
            typi = lj_type_list[i]
            typj = lj_type_list[j]
            # Get any NBFIXes we may have
            try:
                rij, wdij, wdij14, rij14 = typi.nbfix[typj.name]
            except KeyError:
                rij = lj_radii[i] + lj_radii[j]
                wdij = sqrt(lj_depths[i] * lj_depths[j])
                rij14 = lj_radii14[i] + lj_radii14[j]
                wdij14 = sqrt(lj_depths14[i] * lj_depths14[j])
            a = wdij * rij**12
            a14 = wdij14 * rij14**12
            b = 2 * wdij * rij**6
            b14 = 2 * wdij14 * rij14**6
            parm.parm_data['LENNARD_JONES_ACOEF'][index] = a
            parm.parm_data['LENNARD_JONES_BCOEF'][index] = b
            parm.parm_data['LENNARD_JONES_14_ACOEF'][index] = a14
            parm.parm_data['LENNARD_JONES_14_BCOEF'][index] = b14
    # Bonds with and without hydrogen now
    parm.add_flag('BONDS_INC_HYDROGEN', '10I8', num_items=len(bonds_inc_h)*3)
    parm.add_flag('BONDS_WITHOUT_HYDROGEN', '10I8',
                  num_items=len(bonds_without_h)*3)
    pd = parm.parm_data
    for i, bond in enumerate(bonds_inc_h):
        pd['BONDS_INC_HYDROGEN'][3*i  ] = bond.atom1.idx * 3
        pd['BONDS_INC_HYDROGEN'][3*i+1] = bond.atom2.idx * 3
        pd['BONDS_INC_HYDROGEN'][3*i+2] = bond.bond_type.idx + 1
    for i, bond in enumerate(bonds_without_h):
        pd['BONDS_WITHOUT_HYDROGEN'][3*i  ] = bond.atom1.idx * 3
        pd['BONDS_WITHOUT_HYDROGEN'][3*i+1] = bond.atom2.idx * 3
        pd['BONDS_WITHOUT_HYDROGEN'][3*i+2] = bond.bond_type.idx + 1
    # Angles with and without hydrogen now
    parm.add_flag('ANGLES_INC_HYDROGEN', '10I8', num_items=len(angles_inc_h)*4)
    parm.add_flag('ANGLES_WITHOUT_HYDROGEN', '10I8',
                  num_items=len(angles_without_h)*4)
    for i, ang in enumerate(angles_inc_h):
        pd['ANGLES_INC_HYDROGEN'][4*i  ] = ang.atom1.idx * 3
        pd['ANGLES_INC_HYDROGEN'][4*i+1] = ang.atom2.idx * 3
        pd['ANGLES_INC_HYDROGEN'][4*i+2] = ang.atom3.idx * 3
        pd['ANGLES_INC_HYDROGEN'][4*i+3] = ang.angle_type.idx + 1
    for i, ang in enumerate(angles_without_h):
        pd['ANGLES_WITHOUT_HYDROGEN'][4*i  ] = ang.atom1.idx * 3
        pd['ANGLES_WITHOUT_HYDROGEN'][4*i+1] = ang.atom2.idx * 3
        pd['ANGLES_WITHOUT_HYDROGEN'][4*i+2] = ang.atom3.idx * 3
        pd['ANGLES_WITHOUT_HYDROGEN'][4*i+3] = ang.angle_type.idx + 1
    # Dihedrals with and without hydrogen now
    parm.add_flag('DIHEDRALS_INC_HYDROGEN', '10I8',
                 num_items=len(dihedrals_inc_h)*5)
    parm.add_flag('DIHEDRALS_WITHOUT_HYDROGEN', '10I8',
                 num_items=len(dihedrals_without_h)*5)
    for i, dih in enumerate(dihedrals_inc_h):
        if dih.atom3.idx == 0 or dih.atom4.idx == 0:
            a1, a2, a3, a4 = dih.atom4, dih.atom3, dih.atom2, dih.atom1
        else:
            a1, a2, a3, a4 = dih.atom1, dih.atom2, dih.atom3, dih.atom4
        pd['DIHEDRALS_INC_HYDROGEN'][5*i  ] = a1.idx * 3
        pd['DIHEDRALS_INC_HYDROGEN'][5*i+1] = a2.idx * 3
        if not dih.end_groups_active:
            pd['DIHEDRALS_INC_HYDROGEN'][5*i+2] = -a3.idx * 3
        else:
            pd['DIHEDRALS_INC_HYDROGEN'][5*i+2] = a3.idx * 3
        pd['DIHEDRALS_INC_HYDROGEN'][5*i+3] = a4.idx * 3
        pd['DIHEDRALS_INC_HYDROGEN'][5*i+4] = dih.dihedral_type.idx + 1
    for i, dih in enumerate(dihedrals_without_h):
        if dih.atom3.idx == 0 or dih.atom4.idx == 0:
            a1, a2, a3, a4 = dih.atom4, dih.atom3, dih.atom2, dih.atom1
        else:
            a1, a2, a3, a4 = dih.atom1, dih.atom2, dih.atom3, dih.atom4
        pd['DIHEDRALS_WITHOUT_HYDROGEN'][5*i  ] = a1.idx * 3
        pd['DIHEDRALS_WITHOUT_HYDROGEN'][5*i+1] = a2.idx * 3
        if not dih.end_groups_active:
            pd['DIHEDRALS_WITHOUT_HYDROGEN'][5*i+2] = -a3.idx * 3
        else:
            pd['DIHEDRALS_WITHOUT_HYDROGEN'][5*i+2] = a3.idx * 3
        pd['DIHEDRALS_WITHOUT_HYDROGEN'][5*i+3] = a4.idx * 3
        pd['DIHEDRALS_WITHOUT_HYDROGEN'][5*i+4] = dih.dihedral_type.idx + 1
    parm.add_flag('HBOND_ACOEF', '5E16.8', num_items=0)
    parm.add_flag('HBOND_BCOEF', '5E16.8', num_items=0)
    parm.add_flag('HBCUT', '5E16.8', num_items=0)
    parm.add_flag('AMBER_ATOM_TYPE', '20a4',
                 data=[a.type.name[:4] for a in struct.atom_list])
    parm.add_flag('TREE_CHAIN_CLASSIFICATION', '20a4',
                 data=['BLA' for a in struct.atom_list])
    parm.add_flag('JOIN_ARRAY', '10I8', data=[0 for a in struct.atom_list])
    parm.add_flag('IROTAT', '10I8', data=[0 for a in struct.atom_list])
    if hasattr(struct, 'cmap_list'):
        # Do the CMAP terms if we have cmaps
        ncmapt = 0
        ncmap = len(struct.cmap_list)
        cmap_types = []
        for cmap in struct.cmap_list: cmap.cmap_type.idx = -1
        for cmap in struct.cmap_list:
            if cmap.cmap_type.idx == -1:
                cmap_types.append(cmap.cmap_type)
                cmap.cmap_type.idx = ncmapt
                ncmapt += 1
        parm.add_flag('CHARMM_CMAP_COUNT', '2I8', data=[ncmap, ncmapt],
                     comments=['Number of CMAP terms, number of unique CMAP '
                               'parameters']
        )
        parm.add_flag('CHARMM_CMAP_RESOLUTION', '20I4',
                     data=[cm.cmap_type.resolution for cm in struct.cmap_list],
                     comments=['Number of steps along each phi/psi CMAP axis',
                               'for each CMAP_PARAMETER grid']
        )
        # Now add each of the cmap flags
        for i, ct in enumerate(cmap_types):
            parm.add_flag('CHARMM_CMAP_PARAMETER_%02d' % (i+1), '8(F9.5)',
                         data=ct.grid[:])
        # Now add the cmap data
        parm.add_flag('CHARMM_CMAP_INDEX', '6I8', num_items=ncmap*6,
                     comments=['Atom index i,j,k,l,m of the cross term',
                               'and then pointer to CHARMM_CMAP_PARAMETER_n']
        )
        for i, cmap in enumerate(struct.cmap_list):
            if not cmap.consecutive:
                raise CharmmPSFError('Amber-CHARMM implementation only '
                                     'supports consecutive coupled torsions.')
            parm.parm_data['CHARMM_CMAP_INDEX'][6*i  ] = cmap.atom1.idx + 1
            parm.parm_data['CHARMM_CMAP_INDEX'][6*i+1] = cmap.atom2.idx + 1
            parm.parm_data['CHARMM_CMAP_INDEX'][6*i+2] = cmap.atom3.idx + 1
            parm.parm_data['CHARMM_CMAP_INDEX'][6*i+3] = cmap.atom4.idx + 1
            parm.parm_data['CHARMM_CMAP_INDEX'][6*i+4] = cmap.atom5.idx + 1
            parm.parm_data['CHARMM_CMAP_INDEX'][6*i+5] = cmap.cmap_type.idx + 1
    if struct.box is not None:
        nspm = max([atom.marked for atom in struct.atom_list])
        iptres = len(struct.residue_list)
        for i, res in enumerate(struct.residue_list):
            if res.resname in charmm_solvent:
                iptres = i + 1
                break
        nspsol = nspm
        for atom in struct.atom_list:
            if atom.residue.resname in charmm_solvent:
                nspsol = atom.marked
                break
        parm.add_flag('SOLVENT_POINTERS', '3I8', data=[iptres, nspm, nspsol])
        parm.add_flag('ATOMS_PER_MOLECULE', '10i8', num_items=nspm)
        for atom in struct.atom_list:
            parm.parm_data['ATOMS_PER_MOLECULE'][atom.marked-1] += 1
        box_info = [struct.box[3]] + struct.box[:3]
        parm.add_flag('BOX_DIMENSIONS', '5E16.8', data=box_info)
        parm.initialize_topology()
        parm.hasbox = True
        parm.box = struct.box[:]
    else:
        parm.initialize_topology()
        parm.hasbox = False

    # See if we need to set a rst7 attribute (if coordinates are set)
    if hasattr(struct.atom_list[0], 'xx'):
        # Assume if the first atom has an X-coordinate, all coordinates are
        # present
        parm.rst7 = Rst7(natom=len(struct.atom_list))
        crds = []
        for atom in struct.atom_list:
            crds += [atom.xx, atom.xy, atom.xz]
        parm.coords = crds
        parm.hasvels = False
    return parm
