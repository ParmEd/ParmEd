"""
This module contains classes that can serve as a drop-in replacement for the
Amber file reading classes provided by OpenMM for creating an OpenMM System
for performing simulations.

It also pulls the box information from the restart file instead of the topology
file.
"""
from __future__ import division

from chemistry.amber.constants import TINY
from chemistry.amber.readparm import AmberParm, ChamberParm, Rst7
from chemistry.exceptions import APIError, OpenMMError
import chemistry.periodic_table as pt
from math import asin, cos, sin, sqrt, pi
import simtk.openmm as mm
from simtk.openmm.vec3 import Vec3
import simtk.unit as u
from simtk.openmm.app import (forcefield as ff, Topology, element)
from simtk.openmm.app.amberprmtopfile import HCT, OBC1, OBC2, GBn, GBn2
from simtk.openmm.app.internal.customgbforces import (GBSAHCTForce,
                GBSAOBC1Force, GBSAOBC2Force, GBSAGBnForce, GBSAGBn2Force)

WATNAMES = ('WAT', 'HOH')
EPNAMES = ('EP', 'LP')

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class OpenMMAmberParm(AmberParm):
    """
    OpenMM-compatible subclass of AmberParm. This object should still work with
    the ParmEd API while also being compatible with OpenMM's environment
    """
   
    # Define default force groups for all of the bonded terms. This allows them
    # to be turned on and off selectively. This is a way to implement per-term
    # energy decomposition to compare individual components

    BOND_FORCE_GROUP = 0
    ANGLE_FORCE_GROUP = 1
    DIHEDRAL_FORCE_GROUP = 2
    NONBONDED_FORCE_GROUP = 3
    GB_FORCE_GROUP = 3

    def openmm_LJ(self):
        """
        Same as fill_LJ, except it uses 0.5 for the LJ radius for H-atoms with
        no vdW parameters (per OpenMM's standard)

        Returns:
            list, list : The 1st list is the list of all Rmin/2 terms. The
                         2nd is the list of all epsilon (or well depth) terms.
        """
        LJ_radius = []  # empty LJ_radii so it can be re-filled
        LJ_depth = []   # empty LJ_depths so it can be re-filled
        one_sixth = 1 / 6    # we need to raise some numbers to the 1/6th power

        ntypes = self.pointers['NTYPES']
        acoef = self.parm_data['LENNARD_JONES_ACOEF']
        bcoef = self.parm_data['LENNARD_JONES_BCOEF']

        for i in range(ntypes):
            lj_index = self.parm_data["NONBONDED_PARM_INDEX"][ntypes*i+i] - 1
            if acoef[lj_index] < 1.0e-10:
                LJ_radius.append(0.5)
                LJ_depth.append(0)
            else:
                factor = (2 * acoef[lj_index] / bcoef[lj_index])
                LJ_radius.append(pow(factor, one_sixth) * 0.5)
                LJ_depth.append(bcoef[lj_index] / 2 / factor)
      
        # Now check that we haven't modified any off-diagonals, since that will
        # not work with OpenMM
        for i in range(ntypes):
            for j in range(ntypes):
                idx = self.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                rij = LJ_radius[i] + LJ_radius[j]
                wdij = sqrt(LJ_depth[i] * LJ_depth[j])
                a = acoef[idx]
                b = bcoef[idx]
                if a == 0 or b == 0:
                    if a != 0 or b != 0 or (wdij != 0 and rij != 0):
                        raise OpenMMError('Off-diagonal LJ modifications '
                                          'detected. These are incompatible '
                                          'with the OpenMM API')
                elif (abs((a - (wdij * rij**12)) / a) > 1e-6 or
                      abs((b - (2 * wdij * rij**6)) / b) > 1e-6):
                    raise OpenMMError(
                            'Off-diagonal LJ modifications detected. These are '
                            'incompatible with the OpenMM API. Acoef=%s; '
                            'computed=%s. Bcoef=%s; computed=%s' %
                            (acoef, wdij*rij**12, bcoef, 2*wdij*rij**6)
                    )

        return LJ_radius, LJ_depth

    def openmm_14_LJ(self):
        """
        Returns the radii and depths for the LJ interactions between 1-4 pairs.
        For Amber topology files this is the same as the normal LJ parameters,
        but is done here so that OpenMMChamberParm can inherit and override this
        behavior without having to duplicate all of the system creation code.
        """
        return self.openmm_LJ()

    @property
    def topology(self):
        """
        The OpenMM Topology object. Cached when possible, but any changes to the
        topology object lists results in the topology being deleted and rebuilt
        """
        # If anything changed, rebuild the topology
        if not self._topology_changed():
            try:
                return self._topology
            except AttributeError:
                pass
        else:
            self.remake_parm()
      
        self._topology = Topology()

        # Add all of the atoms to the topology file in the same chain
        chain = self._topology.addChain()
        last_residue = None
        for i, atm in enumerate(self.atom_list):
            resnum = atm.residue.idx
            if last_residue != resnum:
                last_residue = resnum
                resname = atm.residue.resname
                res = self._topology.addResidue(resname, chain)
            elem = element.get_by_symbol(pt.Element[atm.element])
            self._topology.addAtom(atm.atname, elem, res)

        # Add bonds to the topology (both with and without hydrogen)
        atoms = list(self._topology.atoms())
        for bnd in self.bonds_inc_h + self.bonds_without_h:
            self._topology.addBond(atoms[bnd.atom1.starting_index],
                                   atoms[bnd.atom2.starting_index])
      
        # Set the box dimensions
        if self.ptr('ifbox'):
            if hasattr(self, 'rst7'):
                self._topology.setUnitCellDimensions(
                        self.rst7.box[:3]*u.angstrom
                )
            else:
                self._topology.setUnitCellDimensions(
                        self.parm_data['BOX_DIMENSIONS'][1:4]*u.angstrom
                )

        return self._topology
   
    def _get_gb_params(self, gb_model=HCT):
        """ Gets the GB parameters. Need this method to special-case GB neck """
        if gb_model is GBn:
            screen = [0.5 for atom in self.atom_list]
            for i, atom in enumerate(self.atom_list):
                if atom.element == 6:
                    screen[i] = 0.48435382330
                elif atom.element == 1:
                    screen[i] = 1.09085413633
                elif atom.element == 7:
                    screen[i] = 0.700147318409
                elif atom.element == 8:
                    screen[i] = 1.06557401132
                elif atom.element == 16:
                    screen[i] = 0.602256336067
        elif gb_model is GBn2:
            # Add non-optimized values as defaults
            alpha = [1.0 for i in self.atom_list]
            beta = [0.8 for i in self.atom_list]
            gamma = [4.85 for i in self.atom_list]
            screen = [0.5 for i in self.atom_list]
            for i, atom in enumerate(self.atom_list):
                if atom.element == 6:
                    screen[i] = 1.058554
                    alpha[i] = 0.733756
                    beta[i] = 0.506378
                    gamma[i] = 0.205844
                elif atom.element == 1:
                    screen[i] = 1.425952
                    alpha[i] = 0.788440
                    beta[i] = 0.798699
                    gamma[i] = 0.437334
                elif atom.element == 7:
                    screen[i] = 0.733599
                    alpha[i] = 0.503364
                    beta[i] = 0.316828
                    gamma[i] = 0.192915
                elif atom.element == 8:
                    screen[i] = 1.061039
                    alpha[i] = 0.867814
                    beta[i] = 0.876635
                    gamma[i] = 0.387882
                elif atom.element == 16:
                    screen[i] = -0.703469
                    alpha[i] = 0.867814
                    beta[i] = 0.876635
                    gamma[i] = 0.387882
        else:
            screen = self.parm_data['SCREEN']

        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        radii = [rad * length_conv for rad in self.parm_data['RADII']]

        if gb_model is GBn2:
            return zip(radii, screen, alpha, beta, gamma)
        return zip(radii, screen)

    def createSystem(self, nonbondedMethod=ff.NoCutoff,
                     nonbondedCutoff=1.0*u.nanometer,
                     constraints=None,
                     rigidWater=True,
                     implicitSolvent=None,
                     implicitSolventKappa=None,
                     implicitSolventSaltConc=0.0*u.moles/u.liter,
                     temperature=298.15*u.kelvin,
                     soluteDielectric=1.0,
                     solventDielectric=78.5,
                     removeCMMotion=True,
                     hydrogenMass=None,
                     ewaldErrorTolerance=0.0005,
                     flexibleConstraints=True,
                     verbose=False,
                     forceNBFIX=False):
        """
        Construct an OpenMM System representing the topology described by the
        prmtop file.

        Parameters:
         -  nonbondedMethod (object=NoCutoff) The method to use for nonbonded
               interactions. Allowed values are NoCutoff, CutoffNonPeriodic,
               CutoffPeriodic, Ewald, or PME.
         -  nonbondedCutoff (distance=1*nanometer) The cutoff distance to use
               for nonbonded interactions.
         -  constraints (object=None) Specifies which bonds or angles should be
               implemented with constraints. Allowed values are None, HBonds,
               AllBonds, or HAngles.
         -  rigidWater (boolean=True) If true, water molecules will be fully
               rigid regardless of the value passed for the constraints argument
         -  implicitSolvent (object=None) If not None, the implicit solvent
               model to use. Allowed values are HCT, OBC1, OBC2, or GBn
         -  implicitSolventKappa (float=None): Debye screening parameter to
               model salt concentrations in GB solvent.
         -  implicitSolventSaltConc (float=0.0*u.moles/u.liter): Salt
               concentration for GB simulations. Converted to Debye length
               `kappa'
         -  temperature (float=298.15*u.kelvin): Temperature used in the salt
               concentration-to-kappa conversion for GB salt concentration term
         -  soluteDielectric (float=1.0) The solute dielectric constant to use
               in the implicit solvent model.
         -  solventDielectric (float=78.5) The solvent dielectric constant to
               use in the implicit solvent model.
         -  removeCMMotion (boolean=True) If true, a CMMotionRemover will be
               added to the System.
         -  hydrogenMass (mass=None) The mass to use for hydrogen atoms bound to
               heavy atoms. Any mass added to a hydrogen is subtracted from the
               heavy atom to keep their total mass the same.
         -  ewaldErrorTolerance (float=0.0005) The error tolerance to use if the
               nonbonded method is Ewald or PME.
         -  flexibleConstraints (bool=True) Are our constraints flexible or not?
         -  verbose (bool=False) Optionally prints out a running progress report
         -  forceNBFIX (bool=False) Debugging option used to force all systems
               to follow the NBFIX code path. When NBFIX is not needed, the
               resulting system will run significantly slower with this option
        """
        # Rebuild the topology file if necessary, and flush the atom property
        # data to the atom list
        if self._topology_changed():
            self.remake_parm()
        else:
            self.atom_list.refresh_data()
        try:
            LJ_radius, LJ_depth = self.openmm_LJ() # Get our LJ parameters
            LJ_14_radius, LJ_14_depth = self.openmm_14_LJ()
            has_nbfix = False
        except OpenMMError:
            has_nbfix = True

        # Set the cutoff distance in nanometers
        cutoff = None
        if nonbondedMethod is not ff.NoCutoff:
            cutoff = nonbondedCutoff
            # Remove units from cutoff
            if u.is_quantity(cutoff):
                cutoff = cutoff.value_in_unit(u.nanometers)

        if nonbondedMethod not in (ff.NoCutoff, ff.CutoffNonPeriodic,
                                   ff.CutoffPeriodic, ff.Ewald, ff.PME):
            raise ValueError('Illegal value for nonbonded method')
        if self.ptr('ifbox') == 0 and nonbondedMethod in (ff.CutoffPeriodic,
                                                          ff.Ewald, ff.PME):
            raise ValueError('Illegal nonbonded method for a '
                             'non-periodic system')
        if implicitSolvent not in (HCT, OBC1, OBC2, GBn, GBn2, None):
            raise ValueError('Illegal implicit solvent model choice.')
        if not constraints in (None, ff.HAngles, ff.HBonds, ff.AllBonds):
            raise ValueError('Illegal constraints choice')
      
        # Define conversion factors
        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        _ambfrc = u.kilocalorie_per_mole/(u.angstrom*u.angstrom)
        _openmmfrc = u.kilojoule_per_mole/(u.nanometer*u.nanometer)
        bond_frc_conv = _ambfrc.conversion_factor_to(_openmmfrc)
        _ambfrc = u.kilocalorie_per_mole/(u.radians*u.radians)
        _openmmfrc = u.kilojoule_per_mole/(u.radians*u.radians)
        angle_frc_conv = _ambfrc.conversion_factor_to(_openmmfrc)
        dihe_frc_conv = u.kilocalorie_per_mole.conversion_factor_to(
                            u.kilojoule_per_mole)
        ene_conv = dihe_frc_conv
      
        # Create the system
        system = mm.System()
        if verbose: print('Adding particles...')
        for mass in self.parm_data['MASS']:
            system.addParticle(mass)
        # Set up the constraints
        if verbose and (constraints is not None and not rigidWater):
            print('Adding constraints...')
        if constraints in (ff.HBonds, ff.AllBonds, ff.HAngles):
            for bond in self.bonds_inc_h:
                system.addConstraint(bond.atom1.starting_index,
                                     bond.atom2.starting_index,
                                     bond.bond_type.req*length_conv)
        if constraints in (ff.AllBonds, ff.HAngles):
            for bond in self.bonds_without_h:
                system.addConstraint(bond.atom1.starting_index,
                                     bond.atom2.starting_index,
                                     bond.bond_type.req*length_conv)
        if rigidWater and constraints is None:
            for bond in self.bonds_inc_h:
                if (bond.atom1.residue.resname in WATNAMES and
                    bond.atom2.residue.resname in WATNAMES):
                    system.addConstraint(bond.atom1.starting_index,
                                         bond.atom2.starting_index,
                                         bond.bond_type.req*length_conv)
        # Add Bond forces
        if verbose: print('Adding bonds...')
        force = mm.HarmonicBondForce()
        force.setForceGroup(self.BOND_FORCE_GROUP)
        if flexibleConstraints or (constraints not in (ff.HBonds, ff.AllBonds,
                                                       ff.HAngles)):
            for bond in self.bonds_inc_h:
                force.addBond(bond.atom1.starting_index,
                              bond.atom2.starting_index,
                              bond.bond_type.req*length_conv,
                              2*bond.bond_type.k*bond_frc_conv)
        if flexibleConstraints or (constraints not in (ff.AllBonds,ff.HAngles)):
            for bond in self.bonds_without_h:
                force.addBond(bond.atom1.starting_index,
                              bond.atom2.starting_index,
                              bond.bond_type.req*length_conv,
                              2*bond.bond_type.k*bond_frc_conv)
        system.addForce(force)
        # Add Angle forces
        if verbose: print('Adding angles...')
        force = mm.HarmonicAngleForce()
        force.setForceGroup(self.ANGLE_FORCE_GROUP)
        if constraints is ff.HAngles:
            num_constrained_bonds = system.getNumConstraints()
            atom_constraints = [[]] * system.getNumParticles()
            for i in range(num_constrained_bonds):
                c = system.getConstraintParameters(i)
                dist = c[2].value_in_unit(u.nanometer)
                atom_constraints[c[0]].append((c[1], dist))
                atom_constraints[c[1]].append((c[0], dist))
        for angle in self.angles_inc_h:
            if constraints is ff.HAngles:
                a1 = angle.atom1.element
                a2 = angle.atom2.element
                a3 = angle.atom3.element
                nh = int(a1==1) + int(a2==1) + int(a3==1)
                constrained = (nh >= 2 or (nh == 1 and a2 == 8))
            else:
                constrained = False # no constraints
            if constrained:
                l1 = l2 = None
                for bond in angle.atom2.bonds:
                    if bond.atom1 is angle.atom1 or bond.atom2 is angle.atom1:
                        l1 = bond.bond_type.req * length_conv
                    elif bond.atom1 is angle.atom3 or bond.atom2 is angle.atom3:
                        l2 = bond.bond_type.req * length_conv
                # Compute the distance between the atoms and add a constraint
                length = sqrt(l1*l1 + l2*l2 - 2*l1*l2*
                              cos(angle.angle_type.theteq))
                system.addConstraint(bond.atom1.starting_index,
                                     bond.atom2.starting_index, length)
            if flexibleConstraints or not constrained:
                force.addAngle(angle.atom1.starting_index,
                               angle.atom2.starting_index,
                               angle.atom3.starting_index,
                               angle.angle_type.theteq,
                               2*angle.angle_type.k*angle_frc_conv)
        for angle in self.angles_without_h:
            force.addAngle(angle.atom1.starting_index,
                           angle.atom2.starting_index,
                           angle.atom3.starting_index,
                           angle.angle_type.theteq,
                           2*angle.angle_type.k*angle_frc_conv)
        system.addForce(force)
        # Add dihedral forces
        if verbose: print('Adding torsions...')
        force = mm.PeriodicTorsionForce()
        force.setForceGroup(self.DIHEDRAL_FORCE_GROUP)
        for tor in self.dihedrals_inc_h + self.dihedrals_without_h:
            force.addTorsion(tor.atom1.starting_index,
                             tor.atom2.starting_index,
                             tor.atom3.starting_index,
                             tor.atom4.starting_index,
                             int(tor.dihed_type.per),
                             tor.dihed_type.phase,
                             tor.dihed_type.phi_k*dihe_frc_conv)
        system.addForce(force)

        # Add nonbonded terms now
        if verbose: print('Adding nonbonded interactions...')
        force = mm.NonbondedForce()
        force.setForceGroup(self.NONBONDED_FORCE_GROUP)
        if self.ptr('ifbox') == 0: # non-periodic
            if nonbondedMethod is ff.NoCutoff:
                force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
            elif nonbondedMethod is ff.CutoffNonPeriodic:
                if cutoff is None:
                    raise ValueError('No cutoff value specified')
                force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
                force.setCutoffDistance(cutoff)
            else:
                raise ValueError('Illegal nonbonded method for non-periodic '
                                 'system')
        else: # periodic
            # Set up box vectors (from inpcrd if available, or fall back to
            # prmtop definitions
            system.setDefaultPeriodicBoxVectors(*self.box_vectors)

            # Set cutoff
            if cutoff is None:
                # Compute cutoff automatically
                box = self.box_lengths
                min_box_width = min((box[0]/u.nanometers,
                                     box[1]/u.nanometers,
                                     box[2]/u.nanometers))
                CLEARANCE_FACTOR = 0.97
                cutoff = u.Quantity((min_box_width*CLEARANCE_FACTOR)/2.0,
                                    u.nanometers)
            if nonbondedMethod is not ff.NoCutoff:
                force.setCutoffDistance(cutoff)

            # Set nonbonded method.
            if nonbondedMethod is ff.NoCutoff:
                force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
            elif nonbondedMethod is ff.CutoffNonPeriodic:
                force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
            elif nonbondedMethod is ff.CutoffPeriodic:
                force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
            elif nonbondedMethod is ff.Ewald:
                force.setNonbondedMethod(mm.NonbondedForce.Ewald)
            elif nonbondedMethod is ff.PME:
                force.setNonbondedMethod(mm.NonbondedForce.PME)
            else:
                raise ValueError('Cutoff method is not understood')

            if ewaldErrorTolerance is not None:
                force.setEwaldErrorTolerance(ewaldErrorTolerance)

        # Add per-particle nonbonded parameters (LJ params)
        sigma_scale = 2**(-1/6) * 2 * length_conv
        if not (has_nbfix or forceNBFIX):
            for i, atm in enumerate(self.atom_list):
                force.addParticle(atm.charge,
                                  sigma_scale*LJ_radius[atm.nb_idx-1],
                                  LJ_depth[atm.nb_idx-1]*ene_conv)
        else:
            for i, atm in enumerate(self.atom_list):
                force.addParticle(atm.charge, 1.0, 0.0)

        excluded_atom_pairs = set() # save these pairs so we don't zero them out
        sigma_scale = 2**(-1/6) * length_conv
        if not (has_nbfix or forceNBFIX):
            # Add 1-4 interactions
            for tor in self.dihedrals_inc_h + self.dihedrals_without_h:
                if min(tor.signs) < 0: continue # multi-terms and impropers
                charge_prod = (tor.atom1.charge * tor.atom4.charge /
                               tor.dihed_type.scee)
                epsilon = (sqrt(LJ_14_depth[tor.atom1.nb_idx-1] *
                                LJ_14_depth[tor.atom4.nb_idx-1]) * ene_conv /
                                tor.dihed_type.scnb)
                sigma = (LJ_14_radius[tor.atom1.nb_idx-1] +
                         LJ_14_radius[tor.atom4.nb_idx-1]) * sigma_scale
                force.addException(tor.atom1.starting_index,
                                   tor.atom4.starting_index,
                                   charge_prod, sigma, epsilon)
                excluded_atom_pairs.add(
                        min( (tor.atom1.starting_index, tor.atom4.starting_index),
                             (tor.atom4.starting_index, tor.atom1.starting_index) )
                )
        else:
            # Add 1-4 interactions to the NonbondedForce object. Support both
            # chamber and regular parms here (this way the ChamberParm
            # createSystem call can still call this function)
            try:
                parm_acoef = self.parm_data['LENNARD_JONES_14_ACOEF']
                parm_bcoef = self.parm_data['LENNARD_JONES_14_BCOEF']
            except KeyError:
                parm_acoef = self.parm_data['LENNARD_JONES_ACOEF']
                parm_bcoef = self.parm_data['LENNARD_JONES_BCOEF']
            nbidx = self.parm_data['NONBONDED_PARM_INDEX']
            ntypes = self.ptr('ntypes')
            for tor in self.dihedrals_inc_h + self.dihedrals_without_h:
                if min(tor.signs) < 0: continue
                charge_prod = (tor.atom1.charge * tor.atom4.charge /
                               tor.dihed_type.scee)
                typ1 = tor.atom1.nb_idx - 1
                typ2 = tor.atom4.nb_idx - 1
                idx = nbidx[ntypes*typ1+typ2] - 1
                b = parm_bcoef[idx]
                a = parm_acoef[idx]
                try:
                    epsilon = b * b / (4 * a) * ene_conv / tor.dihed_type.scnb
                    sigma = (2*a/b)**(1/6) * sigma_scale
                except ZeroDivisionError:
                    if a != 0 or b != 0:
                        raise RuntimeError('Cannot have only one of '
                                    'A-coefficient or B-coefficient be 0.')
                    epsilon = sigma = 0
                force.addException(tor.atom1.starting_index,
                                   tor.atom4.starting_index,
                                   charge_prod, sigma, epsilon)
                excluded_atom_pairs.add(
                        min( (tor.atom1.starting_index, tor.atom4.starting_index),
                             (tor.atom4.starting_index, tor.atom1.starting_index) )
                )

        # Add excluded atoms
        for atom in self.atom_list:
            # Exclude all bonds and angles
            for atom2 in atom.bond_partners:
                if atom2.starting_index > atom.starting_index:
                    force.addException(atom.starting_index,
                                       atom2.starting_index, 0.0, 0.1, 0.0)
            for atom2 in atom.angle_partners:
                if atom2.starting_index > atom.starting_index:
                    force.addException(atom.starting_index,
                                       atom2.starting_index, 0.0, 0.1, 0.0)
            for atom2 in atom.exclusion_partners:
                if atom2.starting_index > atom.starting_index:
                    force.addException(atom.starting_index,
                                       atom2.starting_index, 0.0, 0.1, 0.0)
            for atom2 in atom.dihedral_partners:
                if atom2.starting_index <= atom.starting_index: continue
                if ((atom.starting_index, atom2.starting_index) in
                    excluded_atom_pairs):
                    continue
                force.addException(atom.starting_index,
                                   atom2.starting_index, 0.0, 0.1, 0.0)
        system.addForce(force)

        if has_nbfix or forceNBFIX:
            # Now we need to add a CustomNonbondedForce to handle the vdW
            # potential
            parm_acoef = self.parm_data['LENNARD_JONES_ACOEF']
            parm_bcoef = self.parm_data['LENNARD_JONES_BCOEF']
            acoef = [0 for i in range(ntypes*ntypes)]
            bcoef = acoef[:]
            # Take sqrt of A coefficient to reduce taxing single precision
            # limits (since length_conv is 10, we don't want to needlessly
            # multiply by 10^12 when 10^6 will do)
            afac = sqrt(ene_conv) * length_conv**6
            bfac = ene_conv * length_conv**6
            for i in range(ntypes):
                for j in range(ntypes):
                    idx = nbidx[ntypes*i+j] - 1
                    acoef[i*ntypes+j] = sqrt(parm_acoef[idx]) * afac
                    bcoef[i*ntypes+j] = parm_bcoef[idx] * bfac
            cforce = mm.CustomNonbondedForce('(a/r6)^2-b/r6; r6=r^6;'
                                             'a=acoef(type1, type2);'
                                             'b=bcoef(type1, type2);')
            cforce.addTabulatedFunction('acoef',
                    mm.Discrete2DFunction(ntypes, ntypes, acoef))
            cforce.addTabulatedFunction('bcoef',
                    mm.Discrete2DFunction(ntypes, ntypes, bcoef))
            cforce.addPerParticleParameter('type')
            cforce.setForceGroup(self.NONBONDED_FORCE_GROUP)
            for atom in self.atom_list:
                cforce.addParticle((atom.nb_idx - 1,)) # index from 0
            # Now add the exclusions
            for i in range(force.getNumExceptions()):
                ii, jj, q, eps, sig = force.getExceptionParameters(i)
                cforce.addExclusion(ii, jj)
            if (nonbondedMethod is ff.PME or nonbondedMethod is ff.Ewald or
                        nonbondedMethod is ff.CutoffPeriodic):
                cforce.setNonbondedMethod(cforce.CutoffPeriodic)
                cforce.setCutoffDistance(nonbondedCutoff)
                cforce.setUseLongRangeCorrection(True)
            elif nonbondedMethod is ff.NoCutoff:
                cforce.setNonbondedMethod(cforce.NoCutoff)
            elif nonbondedMethod is ff.CutoffNonPeriodic:
                cforce.setNonbondedMethod(cforce.CutoffNonPeriodic)
                cforce.setCutoffDistance(nonbondedCutoff)
            else:
                raise ValueError('Illegal nonbonded method')
            system.addForce(cforce)

        # Add virtual sites for water
        # First tag the residues that have an extra point in them
        for res in self.residue_list: res.has_ep = False
        ep = [atom for atom in self.atom_list if atom.atname in EPNAMES]
        for atom in ep: atom.residue.has_ep = True
        if len(ep) > 0:
            numRes = ep[-1].residue.idx + 1
            waterO = [[] for i in range(numRes)]
            waterH = [[] for i in range(numRes)]
            waterEP = [[] for i in range(numRes)]
            for atom in self.atom_list:
                if atom.residue.has_ep:
                    if atom.element == 8:
                        waterO[res].append(atom)
                    elif atom.element == 1:
                        waterH[res].append(atom)
                    elif atom.element == 0:
                        waterEP[res].append(atom)
            # Record bond lengths for faster access
            distOH = [None] * numRes
            distHH = [None] * numRes
            distOE = [None] * numRes
            for bond in self.bonds_inc_h + self.bonds_without_h:
                a1 = bond.atom1
                a2 = bond.atom2
                if a1.residue.has_ep:
                    res = a1.residue.idx
                    if a1.element == 1 or a2.element == 1:
                        if a1.element == 1 and a2.element == 1:
                            distHH[res] = bond.bond_type.req * u.angstroms
                        if a1.element == 8 or a2.element == 8:
                            distOH[res] = bond.bond_type.req * u.angstroms
                    elif ((a1.element == 8 or a2.element == 8) and
                          (a1.element == 0 or a2.element == 0)):
                        distOE[res] = bond.bond_type.req * u.angstroms
            # Loop over residues and add the virtual points
            out_of_plane_angle = 54.735 * u.degrees
            cosOOP = u.cos(out_of_plane_angle)
            sinOOP = u.sin(out_of_plane_angle)
            for residue in self.residue_list:
                if not residue.has_ep: continue
                res = residue.idx
                if len(waterO[res]) == 1 and len(waterH[res]) == 2:
                    if len(waterEP[res]) == 1:
                        # 4-point water
                        weightH = distOE[res] / sqrt(distOH[res] * distOH[res] -
                                               0.25 * distHH[res] * distHH[res])
                        system.setVirtualSite(
                                waterEP[res][0],
                                mm.ThreeParticleAverageSite(waterO[res][0],
                                waterH[res][0], waterH[res][1],
                                1-weightH, weightH/2, weightH/2)
                        )
                elif len(waterEP[res]) == 2:
                    # 5-point water
                    weightH = (cosOOP * distOE[res] /
                               sqrt(distOH[res] * distOH[res] -
                                 0.25 * distHH[res] * distHH[res])
                    )
                    angleHOH = 2 * asin(0.5 * distHH[res] / distOH[res])
                    lenCross = distOH[res] * distOH[res] * sin(angleHOH)
                    weightCross = sinOOP * distOE[res] / lenCross
                    site1 = mm.OutOfPlaneSite(waterO[res][0], waterH[res][0],
                            waterH[res][1], weightH/2, weightH/2, weightCross)
                    site2 = mm.OutOfPlaneSite(waterO[res][0], waterH[res][0],
                            waterH[res][1], weightH/2, weightH/2, -weightCross)
                    system.setVirtualSite(waterEP[res][0], site1)
                    system.setVirtualSite(waterEP[res][1], site2)

        # Add GB model if we're doing one
        if implicitSolvent is not None:
            if verbose: print('Adding GB parameters...')
            gb_parms = self._get_gb_params(implicitSolvent)

            # If implicitSolventKappa is None, compute it from salt
            # concentration
            if implicitSolventKappa is None:
                if u.is_quantity(implicitSolventSaltConc):
                    sc = implicitSolventSaltConc.value_in_unit(u.moles/u.liter)
                    implicitSolventSaltConc = sc
                if u.is_quantity(temperature):
                    temperature = temperature.value_in_unit(u.kelvin)
                # The constant is 1 / sqrt( epsilon_0 * kB / (2 * NA * q^2 *
                # 1000) ) where NA is avogadro's number, epsilon_0 is the
                # permittivity of free space, q is the elementary charge (this
                # number matches sander/pmemd's kappa conversion factor)
                implicitSolventKappa = 50.33355 * sqrt(implicitSolventSaltConc /
                                                solventDielectric / temperature)
                # Multiply by 0.73 to account for ion exclusions, and multiply
                # by 10 to convert to 1/nm from 1/angstroms
                implicitSolventKappa *= 7.3
            elif implicitSolvent is None:
                implicitSolventKappa = 0.0

            if u.is_quantity(implicitSolventKappa):
                implicitSolventKappa = implicitSolventKappa.value_in_unit(
                                            (1.0/u.nanometer).unit)
            if implicitSolvent is HCT:
                gb = GBSAHCTForce(solventDielectric, soluteDielectric, None,
                                  cutoff, kappa=implicitSolventKappa)
            elif implicitSolvent is OBC1:
                gb = GBSAOBC1Force(solventDielectric, soluteDielectric, None,
                                   cutoff, kappa=implicitSolventKappa)
            elif implicitSolvent is OBC2:
                gb = GBSAOBC2Force(solventDielectric, soluteDielectric, None,
                                   cutoff, kappa=implicitSolventKappa)
            elif implicitSolvent is GBn:
                gb = GBSAGBnForce(solventDielectric, soluteDielectric, None,
                                  cutoff, kappa=implicitSolventKappa)
            elif implicitSolvent is GBn2:
                gb = GBSAGBn2Force(solventDielectric, soluteDielectric, None,
                                   cutoff, kappa=implicitSolventKappa)
            for i, atom in enumerate(self.atom_list):
                gb.addParticle([atom.charge] + list(gb_parms[i]))
            # Set cutoff method
            if nonbondedMethod is ff.NoCutoff:
                gb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
            elif nonbondedMethod is ff.CutoffNonPeriodic:
                gb.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
                gb.setCutoffDistance(cutoff)
            elif nonbondedMethod is ff.CutoffPeriodic:
                gb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
                gb.setCutoffDistance(cutoff)
            else:
                raise ValueError('Illegal nonbonded method for use with GBSA')
            gb.setForceGroup(self.GB_FORCE_GROUP)
            system.addForce(gb)
            force.setReactionFieldDielectric(1.0) # applies to NonbondedForce

        # See if we repartition the hydrogen masses
        if hydrogenMass is not None:
            if not u.is_quantity(hydrogenMass):
                hydrogenMass *= u.daltons
            for bond in self.bonds_inc_h:
                atom1, atom2 = bond.atom1, bond.atom2
                if atom1.element == 1:
                    atom1, atom2 = atom2, atom1 # now atom2 is hydrogen for sure
                if atom1.element != 1:
                    transfer_mass = hydrogenMass - atom2.mass * u.daltons
                    new_mass1 = (system.getParticleMass(atom1.starting_index) -
                                 transfer_mass)
                    system.setParticleMass(atom2.starting_index, hydrogenMass)
                    system.setParticleMass(atom1.starting_index, new_mass1)
        # See if we want to remove COM motion
        if removeCMMotion:
            system.addForce(mm.CMMotionRemover())

        # Cache our system for easy access
        self._system = system

        return system

    @property
    def system(self):
        """
        Return the cached system class -- it needs to be initialized via
        "createSystem" first!
        """
        try:
            return self._system
        except AttributeError:
            raise APIError('You must initialize the system with createSystem '
                           'before accessing the cached object.')

    @property
    def positions(self):
        """
        Return the cached positions or create new ones from the atoms
        """
        try:
            if len(self._positions) == len(self.atom_list):
                return self._positions
        except AttributeError:
            pass

        self._positions = tuple([Vec3(a.xx, a.xy, a.xz)
                               for a in self.atom_list]) * u.angstroms
        return self._positions

    @positions.setter
    def positions(self, stuff):
        """
        Update the cached positions and assign the coordinates to the atoms
        """
        self._positions = stuff
        for i, pos in enumerate(stuff.value_in_unit(u.angstroms)):
            i3 = i * 3
            atom = self.atom_list[i]
            atom.xx, atom.xy, atom.xz = pos
            self.coords[i3], self.coords[i3+1], self.coords[i3+2] = pos

    @property
    def velocities(self):
        """ Same as for positions, but for velocities """
        try:
            if len(self._velocities) == len(self.atom_list):
                return self._velocities
        except AttributeError:
            pass

        self._velocities = tuple([Vec3(a.vx, a.vy, a.vz)
                    for a in self.atom_list]) * (u.angstroms/u.picosecond) 
        return self._velocities

    @velocities.setter
    def velocities(self, stuff):
        self._velocities = stuff
        for atom, vel in zip(self.atom_list, stuff):
            atom.vx, atom.vy, atom.vz = vel.value_in_unit(
                    u.angstroms/u.picoseconds)

    @property
    def box_vectors(self):
        """ Return tuple of box vectors """
        if hasattr(self, 'rst7'):
            box = [x*u.angstrom for x in self.rst7.box[:3]]
            ang = [self.rst7.box[3], self.rst7.box[4], self.rst7.box[5]]
            return _box_vectors_from_lengths_angles(box[0], box[1], box[2],
                                                    ang[0], ang[1], ang[2])
        else:
            box = [x*u.angstrom for x in self.parm_data['BOX_DIMENSIONS'][1:]]
            ang = [self.parm_data['BOX_DIMENSIONS'][0]] * 3
            return _box_vectors_from_lengths_angles(box[0], box[1], box[2],
                                                    ang[0], ang[1], ang[2])

    @property
    def box_lengths(self):
        """ Return tuple of 3 units """
        if hasattr(self, 'rst7'):
            box = [x*u.angstrom for x in self.rst7.box[:3]]
        else:
            box = [x*u.angstrom for x in self.parm_data['BOX_DIMENSIONS'][1:]]
        return tuple(box)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class OpenMMChamberParm(ChamberParm, OpenMMAmberParm):
    """
    OpenMM-compatible subclass of AmberParm. This object should still work with
    the ParmEd API while also being compatible with OpenMM's environment
    """
   
    # Define default force groups for all of the bonded terms. This allows them
    # to be turned on and off selectively. This is a way to implement per-term
    # energy decomposition to compare individual components

    BOND_FORCE_GROUP = 0
    ANGLE_FORCE_GROUP = 1
    DIHEDRAL_FORCE_GROUP = 2
    UREY_BRADLEY_FORCE_GROUP = 3
    IMPROPER_FORCE_GROUP = 4
    CMAP_FORCE_GROUP = 5
    NONBONDED_FORCE_GROUP = 6
    GB_FORCE_GROUP = 6

    def openmm_14_LJ(self):
        """
        Same as fill_14_LJ, except it uses 0.5 for the LJ radius for H-atoms
        with no vdW parameters (per OpenMM's standard)

        Returns:
            list, list : The 1st list is the list of all Rmin/2 terms. The
                         2nd is the list of all epsilon (or well depth) terms.
        """
        LJ_radius = []  # empty LJ_radii so it can be re-filled
        LJ_depth = []   # empty LJ_depths so it can be re-filled
        one_sixth = 1 / 6    # we need to raise some numbers to the 1/6th power

        ntypes = self.pointers['NTYPES']
        acoef = self.parm_data['LENNARD_JONES_14_ACOEF']
        bcoef = self.parm_data['LENNARD_JONES_14_BCOEF']
        for i in range(ntypes):
            lj_index = self.parm_data["NONBONDED_PARM_INDEX"][ntypes*i+i] - 1
            if acoef[lj_index] < 1.0e-10:
                LJ_radius.append(0.5)
                LJ_depth.append(0)
            else:
                factor = 2 * acoef[lj_index] / bcoef[lj_index]
                LJ_radius.append(pow(factor, one_sixth) * 0.5)
                LJ_depth.append(bcoef[lj_index] / 2 / factor)

        # Now check that we haven't modified any off-diagonals, since that will
        # not work with OpenMM
        for i in range(ntypes):
            for j in range(ntypes):
                idx = self.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                rij = LJ_radius[i] + LJ_radius[j]
                wdij = sqrt(LJ_depth[i] * LJ_depth[j])
                a = acoef[idx]
                b = bcoef[idx]
                if a == 0 or b == 0:
                    if a != 0 or b != 0 or (wdij != 0 and rij != 0):
                        raise OpenMMError('Off-diagonal LJ modifications '
                                          'detected. These are incompatible '
                                          'with the OpenMM API')
                elif (abs((a - (wdij * rij**12)) / a) > 1e-6 or
                      abs((b - (2 * wdij * rij**6)) / b) > 1e-6):
                    raise OpenMMError(
                            'Off-diagonal LJ modifications detected. These are'
                            ' incompatible with the OpenMM API. Acoef=%s; '
                            'computed=%s. Bcoef=%s; computed=%s' %
                            (acoef, wdij*rij**12, bcoef, 2*wdij*rij**6)
                    )

        return LJ_radius, LJ_depth

    def createSystem(self, nonbondedMethod=ff.NoCutoff,
                     nonbondedCutoff=1.0*u.nanometer,
                     switchDistance=0.0*u.nanometer,
                     constraints=None,
                     rigidWater=True,
                     implicitSolvent=None,
                     implicitSolventKappa=None,
                     implicitSolventSaltConc=0.0*u.moles/u.liter,
                     temperature=298.15*u.kelvin,
                     soluteDielectric=1.0,
                     solventDielectric=78.5,
                     removeCMMotion=True,
                     hydrogenMass=None,
                     ewaldErrorTolerance=0.0005,
                     flexibleConstraints=True,
                     verbose=False,
                     forceNBFIX=False):
        """
        Construct an OpenMM System representing the topology described by the
        prmtop file.

        Parameters:
         -  nonbondedMethod (object=NoCutoff) The method to use for nonbonded
               interactions. Allowed values are NoCutoff, CutoffNonPeriodic,
               CutoffPeriodic, Ewald, or PME.
         -  nonbondedCutoff (distance=1*nanometer) The cutoff distance to use
               for nonbonded interactions.
         -  switchDistance (distance=0*nanometer) The distance at which the
               switching function is active for van der Waals interactions. If
               the switchDistance evaluates to boolean False (e.g., if it is 0),
               no switching function will be used. Illegal values raise a
               ValueError
         -  constraints (object=None) Specifies which bonds or angles should be
               implemented with constraints. Allowed values are None, HBonds,
               AllBonds, or HAngles.
         -  rigidWater (boolean=True) If true, water molecules will be fully
               rigid regardless of the value passed for the constraints argument
         -  implicitSolvent (object=None) If not None, the implicit solvent
               model to use. Allowed values are HCT, OBC1, OBC2, or GBn
         -  implicitSolventKappa (float=None): Debye screening parameter to
               model salt concentrations in GB solvent.
         -  implicitSolventSaltConc (float=0.0*u.moles/u.liter): Salt
               concentration for GB simulations. Converted to Debye length
               `kappa'
         -  temperature (float=298.15*u.kelvin): Temperature used in the salt
               concentration-to-kappa conversion for GB salt concentration term
         -  soluteDielectric (float=1.0) The solute dielectric constant to use
               in the implicit solvent model.
         -  solventDielectric (float=78.5) The solvent dielectric constant to
               use in the implicit solvent model.
         -  removeCMMotion (boolean=True) If true, a CMMotionRemover will be
               added to the System.
         -  hydrogenMass (mass=None) The mass to use for hydrogen atoms bound to
               heavy atoms. Any mass added to a hydrogen is subtracted from the
               heavy atom to keep their total mass the same.
         -  ewaldErrorTolerance (float=0.0005) The error tolerance to use if the
               nonbonded method is Ewald or PME.
         -  flexibleConstraints (bool=True) Are our constraints flexible or not?
         -  verbose (bool=False) Optionally prints out a running progress report
         -  forceNBFIX (bool=False) Debugging option used to force all systems
               to follow the NBFIX code path. When NBFIX is not needed, the
               resulting system will run significantly slower with this option
        """
        # Start from the system created by the AmberParm class
        system = super(OpenMMChamberParm, self).createSystem(
                        nonbondedMethod, nonbondedCutoff, constraints,
                        rigidWater, implicitSolvent, implicitSolventKappa,
                        implicitSolventSaltConc, temperature, soluteDielectric,
                        solventDielectric, removeCMMotion, hydrogenMass,
                        ewaldErrorTolerance, flexibleConstraints, verbose,
                        forceNBFIX)

        # Define conversion factors
        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        _ambfrc = u.kilocalorie_per_mole/(u.angstrom*u.angstrom)
        _openmmfrc = u.kilojoule_per_mole/(u.nanometer*u.nanometer)
        bond_frc_conv = _ambfrc.conversion_factor_to(_openmmfrc)
        _ambfrc = u.kilocalorie_per_mole/(u.radians*u.radians)
        _openmmfrc = u.kilojoule_per_mole/(u.radians*u.radians)
        dihe_frc_conv = u.kilocalorie_per_mole.conversion_factor_to(
                                    u.kilojoule_per_mole)
        ene_conv = dihe_frc_conv

        # Add the urey-bradley terms
        if verbose: print('Adding Urey-Bradley terms')
        force = mm.HarmonicBondForce()
        force.setForceGroup(self.UREY_BRADLEY_FORCE_GROUP)
        for ub in self.urey_bradley:
            force.addBond(ub.atom1.starting_index,
                          ub.atom2.starting_index,
                          ub.ub_type.req*length_conv,
                          2*ub.ub_type.k*bond_frc_conv)
        system.addForce(force)

        if verbose: print('Adding impropers...')
        # Ick. OpenMM does not have an improper torsion class. Need to
        # construct one from CustomTorsionForce
        force = mm.CustomTorsionForce('k*(theta-theta0)^2')
        force.addPerTorsionParameter('k')
        force.addPerTorsionParameter('theta0')
        force.setForceGroup(self.IMPROPER_FORCE_GROUP)
        for imp in self.improper:
            force.addTorsion(imp.atom1.starting_index,
                             imp.atom2.starting_index,
                             imp.atom3.starting_index,
                             imp.atom4.starting_index,
                             (imp.improp_type.psi_k*dihe_frc_conv,
                              imp.improp_type.psi_eq*pi/180)
            )
        system.addForce(force)

        if self.has_cmap:
            if verbose: print('Adding CMAP coupled torsions...')
            force = mm.CMAPTorsionForce()
            force.setForceGroup(self.CMAP_FORCE_GROUP)
            # First get the list of cmap maps we're going to use. Just store the
            # IDs so we have simple integer comparisons to do later
            cmap_type_list = []
            cmap_map = dict()
            for cmap in self.cmap:
                if not id(cmap.cmap_type) in cmap_type_list:
                    ct = cmap.cmap_type
                    cmap_type_list.append(id(ct))
                    # Our torsion correction maps need to go from 0 to 360
                    # degrees
                    grid = ct.grid.switch_range().T
                    m = force.addMap(ct.resolution, [x*ene_conv for x in grid])
                    cmap_map[id(ct)] = m
            # Now add in all of the cmaps
            for cmap in self.cmap:
                force.addTorsion(cmap_map[id(cmap.cmap_type)],
                                 cmap.atom1.starting_index,
                                 cmap.atom2.starting_index,
                                 cmap.atom3.starting_index,
                                 cmap.atom4.starting_index, # end of 1st torsion
                                 cmap.atom2.starting_index,
                                 cmap.atom3.starting_index,
                                 cmap.atom4.starting_index,
                                 cmap.atom5.starting_index  # end of 2nd torsion
                )
            system.addForce(force)

        # Now see if we need to toggle the switching function
        if switchDistance and nonbondedMethod is not ff.NoCutoff:
            for force in system.getForces():
                if isinstance(force, mm.NonbondedForce): break
            if switchDistance >= nonbondedCutoff:
                raise ValueError('switchDistance must be smaller than the '
                                 'cutoff!')
            if abs(switchDistance) != switchDistance:
                # Identifies negative floats and Quantity's
                raise ValueError('switchDistance must be non-negative!')
            force.setUseSwitchingFunction(True)
            force.setSwitchingDistance(switchDistance)

        return system

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class OpenMMRst7(Rst7):
    """
    Contains positions and maybe velocities and box information from a restart
    file using OpenMM data structures
    """

    @property
    def positions(self):
        """ Positions as a sequence of Vec3 objects with angstrom units """
        # Return the cached copy of positions
        try:
            return self._positions
        except AttributeError:
            pass

        # The first time it's called, cache the data structure
        self._positions = tuple([Vec3(self.coordinates[3*i],
                self.coordinates[3*i+1], self.coordinates[3*i+2])
                for i in xrange(self.natom)]) * u.angstroms

        return self._positions

    @property
    def velocities(self):
        """ Same as positions, but for velocities """
        try:
            return self._velocities
        except AttributeError:
            pass
      
        self._velocities = tuple([Vec3(self._ambvels[3*i],
            self._ambvels[3*i+1], self._ambvels[3*i+2])
            for i in xrange(self.natom)]) * u.angstroms / u.picoseconds 

        return self._velocities

    @velocities.setter
    def velocities(self, stuff):
        """
        This is a hack to work around the new Amber Rst7 class using the name
        'velocities'. This forces the Amber object to store the velocities in a
        special '_ambvels' attribute instead
        """
        self._ambvels = stuff

    @property
    def box_vectors(self):
        """ Return tuple of box vectors """
        box = [x*u.angstrom for x in self.box[:3]]
        ang = [self.box[3], self.box[4], self.box[5]]
        return _box_vectors_from_lengths_angles(box[0], box[1], box[2],
                                                ang[0], ang[1], ang[2])

    @property
    def box_lengths(self):
        """ Return tuple of floats """
        box = [x*u.angstrom for x in self.rst7.box[:3]]
        return tuple(box)

def _box_vectors_from_lengths_angles(a, b, c, alpha, beta, gamma):
    """
    This method takes the lengths and angles from a unit cell and creates unit
    cell vectors.

    Parameters:
        - a (unit, dimension length): Length of the first vector
        - b (unit, dimension length): Length of the second vector
        - c (unit, dimension length): Length of the third vector
        - alpha (float): Angle between b and c in degrees
        - beta (float): Angle between a and c in degrees
        - gamma (float): Angle between a and b in degrees

    Returns:
        Tuple of box vectors (as Vec3 instances)
    """
    if not (u.is_quantity(a) and u.is_quantity(b) and u.is_quantity(c)):
        raise TypeError('a, b, and c must be units of dimension length')
    if u.is_quantity(alpha): alpha = alpha.value_in_unit(u.degree)
    if u.is_quantity(beta): beta = beta.value_in_unit(u.degree)
    if u.is_quantity(gamma): gamma = gamma.value_in_unit(u.degree)
    a = a.value_in_unit(u.angstrom)
    b = b.value_in_unit(u.angstrom)
    c = c.value_in_unit(u.angstrom)

    if alpha <= 2 * pi and beta <= 2 * pi and gamma <= 2 * pi:
        raise ValueError('box angles must be given in degrees')

    alpha *= pi / 180
    beta *= pi / 180
    gamma *= pi / 180

    av = Vec3(a, 0.0, 0.0) * u.angstrom
    bx = b * cos(gamma)
    by = b * sin(gamma)
    bz = 0.0
    cx = c * cos(beta)
    cy = c * (cos(alpha) - cos(beta) * cos(gamma))
    cz = sqrt(c * c - cx * cx - cy * cy)
   
    # Make sure any components that are close to zero are set to zero exactly
    if abs(bx) < TINY: bx = 0.0
    if abs(by) < TINY: by = 0.0
    if abs(cx) < TINY: cx = 0.0
    if abs(cy) < TINY: cy = 0.0
    if abs(cz) < TINY: cz = 0.0

    bv = Vec3(bx, by, bz) * u.angstrom
    cv = Vec3(cx, cy, cz) * u.angstrom

    return (av, bv, cv)
