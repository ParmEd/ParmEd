"""
This module contains OpenMM-enabled components of the CHARMM classes

Author: Jason M. Swails
Contributors:
Date: April 20, 2014
"""
from __future__ import division

from compat24 import property
from chemistry.charmm.charmmcrds import CharmmCrdFile, CharmmRstFile
from chemistry.charmm.parameters import element_by_mass
from chemistry.charmm.psf import CharmmPsfFile
from chemistry import TrackedList
from chemistry.exceptions import APIError
from math import sqrt, cos, pi, sin, acos
import re
import simtk.openmm as mm
from simtk.openmm.vec3 import Vec3
import simtk.unit as u
from simtk.openmm.app import (forcefield as ff, Topology, element)
from simtk.openmm.app.amberprmtopfile import HCT, OBC1, OBC2, GBn, GBn2
from simtk.openmm.app.internal.customgbforces import (GBSAHCTForce,
                GBSAOBC1Force, GBSAOBC2Force, GBSAGBnForce, GBSAGBn2Force,
                convertParameters)

TINY = 1e-8
WATNAMES = ('WAT', 'HOH', 'TIP3', 'TIP4', 'TIP5', 'SPCE', 'SPC')

class OpenMMCharmmPsfFile(CharmmPsfFile):
    """ Same as CharmmPsfFile but implements OpenMM functionality on top """

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
    
    @property
    def topology(self):
        """ Create an OpenMM Topology object from the stored bonded network """
        try:
            return self._topology
        except AttributeError:
            # If none exists, we need to create it
            pass
        # Cache the topology for easy returning later
        self._topology = topology = Topology()
        
        last_chain = None
        last_residue = None
        # Add each chain (separate 'system's) and residue
        for atom in self.atoms:
            if atom.residue.chain != last_chain:
                chain = topology.addChain()
                last_chain = atom.residue.chain
            if atom.residue.idx != last_residue:
                last_residue = atom.residue.idx
                residue = topology.addResidue(atom.residue.name, chain)
            element_name = element_by_mass(atom.mass)
            elem = element.get_by_symbol(element_name)
            topology.addAtom(atom.name, elem, residue)

        # Add all of the bonds
        atoms = list(topology.atoms())
        # Assign atom indexes to make sure they're current
        for bond in self.bonds:
            topology.addBond(atoms[bond.atom1.idx], atoms[bond.atom2.idx])

        # Add the periodic box if there is one
        if hasattr(self, 'box') and self.box is not None:
            topology.setUnitCellDimensions(self.box[:3] * u.angstroms)
        else:
            self.box_vectors = self.box = None

        return topology

    def _get_gb_params(self, gb_model=HCT):
        """ Gets the GB parameters. Need this method to special-case GB neck """
        screen = [0 for atom in self.atoms]
        if gb_model is GBn:
            radii = _bondi_radii(self.atoms)
            screen = [0.5 for atom in self.atoms]
            for i, atom in enumerate(self.atoms):
                if atom.atomic_number == 6:
                    screen[i] = 0.48435382330
                elif atom.atomic_number == 1:
                    screen[i] = 1.09085413633
                elif atom.atomic_number == 7:
                    screen[i] = 0.700147318409
                elif atom.atomic_number == 8:
                    screen[i] = 1.06557401132
                elif atom.atomic_number == 16:
                    screen[i] = 0.602256336067
        elif gb_model is GBn2:
            radii = _mbondi3_radii(self.atoms)
            # Add non-optimized values as defaults
            alpha = [1.0 for i in self.atoms]
            beta = [0.8 for i in self.atoms]
            gamma = [4.85 for i in self.atoms]
            screen = [0.5 for i in self.atoms]
            for i, atom in enumerate(self.atoms):
                if atom.atomic_number == 6:
                    screen[i] = 1.058554
                    alpha[i] = 0.733756
                    beta[i] = 0.506378
                    gamma[i] = 0.205844
                elif atom.atomic_number == 1:
                    screen[i] = 1.425952
                    alpha[i] = 0.788440
                    beta[i] = 0.798699
                    gamma[i] = 0.437334
                elif atom.atomic_number == 7:
                    screen[i] = 0.733599
                    alpha[i] = 0.503364
                    beta[i] = 0.316828
                    gamma[i] = 0.192915
                elif atom.atomic_number == 8:
                    screen[i] = 1.061039
                    alpha[i] = 0.867814
                    beta[i] = 0.876635
                    gamma[i] = 0.387882
                elif atom.atomic_number == 16:
                    screen[i] = -0.703469
                    alpha[i] = 0.867814
                    beta[i] = 0.876635
                    gamma[i] = 0.387882
        else:
            # Set the default screening parameters
            for i, atom in enumerate(self.atoms):
                if atom.atomic_number == 1:
                    screen[i] = 0.85
                elif atom.atomic_number == 6:
                    screen[i] = 0.72
                elif atom.atomic_number == 7:
                    screen[i] = 0.79
                elif atom.atomic_number == 8:
                    screen[i] = 0.85
                elif atom.atomic_number == 9:
                    screen[i] = 0.88
                elif atom.atomic_number == 15:
                    screen[i] = 0.86
                elif atom.atomic_number == 16:
                    screen[i] = 0.96
                else:
                    screen[i] = 0.8
            # Determine which radii set we need
            if gb_model is OBC1 or gb_model is OBC2:
                radii = _mbondi2_radii(self.atoms)
            elif gb_model is HCT:
                radii = _mbondi_radii(self.atoms)

        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        radii = [x * length_conv for x in radii]

        if gb_model is GBn2:
            return zip(radii, screen, alpha, beta, gamma)
        return zip(radii, screen)

    def createSystem(self, params, nonbondedMethod=ff.NoCutoff,
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
                     verbose=False):
        """
        Construct an OpenMM System representing the topology described by the
        prmtop file.

        Parameters:
         -  params (CharmmParameterSet) All of the parameters for the system
         -  nonbondedMethod (object=NoCutoff) The method to use for nonbonded
               interactions. Allowed values are NoCutoff, CutoffNonPeriodic,
               CutoffPeriodic, Ewald, or PME.
         -  nonbondedCutoff (distance=1*nanometer) The cutoff distance to use
               for nonbonded interactions.
         -  switchDistance (distance=0*nanometer) The distance at which the
               switching function is active for nonbonded interactions. If the
               switchDistance evaluates to boolean False (if it is 0), no
               switching function will be used. Illegal values will raise a
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
        """
        # Load the parameter set
        self.load_parameters(params.condense())

        hasbox = self.topology.getUnitCellDimensions() is not None
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
        if not hasbox and nonbondedMethod in (ff.CutoffPeriodic,
                                              ff.Ewald, ff.PME):
            raise ValueError('Illegal nonbonded method for a '
                             'non-periodic system')
        if nonbondedMethod is not ff.NoCutoff:
            if cutoff is None or cutoff <= 0:
                raise ValueError('Illegal cutoff value [%s] when cutoff is '
                                 'employed' % (nonbondedMethod))
        if implicitSolvent not in (HCT, OBC1, OBC2, GBn, GBn2, None):
            raise ValueError('Illegal implicit solvent model choice.')
        if not constraints in (None, ff.HAngles, ff.HBonds, ff.AllBonds):
            raise ValueError('Illegal constraints choice')
      
        # Define conversion factors
        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        _chmfrc = u.kilocalorie_per_mole/(u.angstrom*u.angstrom)
        _openmmfrc = u.kilojoule_per_mole/(u.nanometer*u.nanometer)
        bond_frc_conv = _chmfrc.conversion_factor_to(_openmmfrc)
        _chmfrc = u.kilocalorie_per_mole/(u.radians*u.radians)
        _openmmfrc = u.kilojoule_per_mole/(u.radians*u.radians)
        angle_frc_conv = _chmfrc.conversion_factor_to(_openmmfrc)
        dihe_frc_conv = u.kilocalorie_per_mole.conversion_factor_to(
                            u.kilojoule_per_mole)
        ene_conv = dihe_frc_conv
      
        # Create the system and determine if any of our atoms have NBFIX (and
        # therefore requires a CustomNonbondedForce instead)
        typenames = set()
        system = mm.System()
        if verbose: print('Adding particles...')
        for atom in self.atoms:
            typenames.add(str(atom.atom_type))
            system.addParticle(atom.mass)
        has_nbfix_terms = False
        typenames = list(typenames)
        try:
            for i, typename in enumerate(typenames):
                typ = params.atom_types_str[typename]
                for j in xrange(i, len(typenames)):
                    if typenames[j] in typ.nbfix:
                        has_nbfix_terms = True
                        raise StopIteration
        except StopIteration:
            pass
        # Set up the constraints
        if verbose and (constraints is not None and not rigidWater):
            print('Adding constraints...')
        if constraints in (ff.HBonds, ff.AllBonds, ff.HAngles):
            for bond in self.bonds:
                if (bond.atom1.atomic_number != 1 and
                    bond.atom2.atomic_number != 1):
                    continue
                system.addConstraint(bond.atom1.idx, bond.atom2.idx,
                                     bond.type.req*length_conv)
        if constraints in (ff.AllBonds, ff.HAngles):
            for bond in self.bonds:
                if (bond.atom1.atomic_number == 1 or
                    bond.atom2.atomic_number == 1):
                    continue
                system.addConstraint(bond.atom1.idx, bond.atom2.idx,
                                     bond.type.req*length_conv)
        if rigidWater and constraints is None:
            for bond in self.bonds:
                if (bond.atom1.atomic_number != 1 and
                    bond.atom2.atomic_number != 1):
                    continue
                if (bond.atom1.residue.name in WATNAMES and
                    bond.atom2.residue.name in WATNAMES):
                    system.addConstraint(bond.atom1.idx, bond.atom2.idx,
                                         bond.type.req*length_conv)
        # Add Bond forces
        if verbose: print('Adding bonds...')
        force = mm.HarmonicBondForce()
        force.setForceGroup(self.BOND_FORCE_GROUP)
        # See which, if any, energy terms we omit
        omitall = not flexibleConstraints and constraints is ff.AllBonds
        omith = omitall or (flexibleConstraints and constraints in
                            (ff.HBonds, ff.AllBonds, ff.HAngles))
        for bond in self.bonds:
            if omitall: continue
            if omith and (bond.atom1.atomic_number == 1 or
                          bond.atom2.atomic_number == 1):
                continue
            force.addBond(bond.atom1.idx, bond.atom2.idx,
                          bond.type.req*length_conv,
                          2*bond.type.k*bond_frc_conv)
        system.addForce(force)
        # Add Angle forces
        if verbose: print('Adding angles...')
        force = mm.HarmonicAngleForce()
        force.setForceGroup(self.ANGLE_FORCE_GROUP)
        if constraints is ff.HAngles:
            num_constrained_bonds = system.getNumConstraints()
            atom_constraints = [[]] * system.getNumParticles()
            for i in xrange(num_constrained_bonds):
                c = system.getConstraintParameters(i)
                dist = c[2].value_in_unit(u.nanometer)
                atom_constraints[c[0]].append((c[1], dist))
                atom_constraints[c[1]].append((c[0], dist))
        for angle in self.angles:
            # Only constrain angles including hydrogen here
            if (angle.atom1.atomic_number != 1 and
                angle.atom2.atomic_number != 1 and
                angle.atom3.atomic_number != 1):
                continue
            if constraints is ff.HAngles:
                a1 = angle.atom1.atomic_number
                a2 = angle.atom2.atomic_number
                a3 = angle.atom3.atomic_number
                nh = int(a1==1) + int(a2==1) + int(a3==1)
                constrained = (nh >= 2 or (nh == 1 and a2 == 8))
            else:
                constrained = False # no constraints
            if constrained:
                l1 = l2 = None
                for bond in angle.atom2.bonds:
                    if bond.atom1 is angle.atom1 or bond.atom2 is angle.atom1:
                        l1 = bond.type.req * length_conv
                    elif bond.atom1 is angle.atom3 or bond.atom2 is angle.atom3:
                        l2 = bond.type.req * length_conv
                # Compute the distance between the atoms and add a constraint
                length = sqrt(l1*l1 + l2*l2 - 2*l1*l2*
                              cos(angle.type.theteq))
                system.addConstraint(bond.atom1.idx, bond.atom2.idx, length)
            if flexibleConstraints or not constrained:
                force.addAngle(angle.atom1.idx, angle.atom2.idx,
                               angle.atom3.idx, angle.type.theteq,
                               2*angle.type.k*angle_frc_conv)
        for angle in self.angles:
            # Already did the angles with hydrogen above. So skip those here
            if (angle.atom1.atomic_number == 1 or
                angle.atom2.atomic_number == 1 or
                angle.atom3.atomic_number == 1):
                continue
            force.addAngle(angle.atom1.idx, angle.atom2.idx,
                           angle.atom3.idx, angle.type.theteq,
                           2*angle.type.k*angle_frc_conv)
        system.addForce(force)

        # Add the urey-bradley terms
        if verbose: print('Adding Urey-Bradley terms')
        force = mm.HarmonicBondForce()
        force.setForceGroup(self.UREY_BRADLEY_FORCE_GROUP)
        for ub in self.urey_bradleys:
            force.addBond(ub.atom1.idx, ub.atom2.idx,
                          ub.type.req*length_conv,
                          2*ub.type.k*bond_frc_conv)
        system.addForce(force)

        # Add dihedral forces
        if verbose: print('Adding torsions...')
        force = mm.PeriodicTorsionForce()
        force.setForceGroup(self.DIHEDRAL_FORCE_GROUP)
        for tor in self.dihedrals:
            for typ in tor.type:
                force.addTorsion(tor.atom1.idx, tor.atom2.idx, tor.atom3.idx,
                                 tor.atom4.idx, int(typ.per), typ.phase,
                                 typ.phi_k*dihe_frc_conv)
        system.addForce(force)

        if verbose: print('Adding impropers...')
        # Ick. OpenMM does not have an improper torsion class. Need to
        # construct one from CustomTorsionForce
        force = mm.CustomTorsionForce('k*(theta-theta0)^2')
        force.addPerTorsionParameter('k')
        force.addPerTorsionParameter('theta0')
        force.setForceGroup(self.IMPROPER_FORCE_GROUP)
        for imp in self.impropers:
            force.addTorsion(imp.atom1.idx, imp.atom2.idx,
                             imp.atom3.idx, imp.atom4.idx,
                             (imp.type.psi_k*dihe_frc_conv,
                              imp.type.psi_eq)
            )
        system.addForce(force)

        if hasattr(self, 'cmaps'):
            if verbose: print('Adding CMAP coupled torsions...')
            force = mm.CMAPTorsionForce()
            force.setForceGroup(self.CMAP_FORCE_GROUP)
            # First get the list of cmap maps we're going to use. Just store the
            # IDs so we have simple integer comparisons to do later
            cmap_types = []
            cmap_map = dict()
            for cmap in self.cmaps:
                if not id(cmap.type) in cmap_types:
                    ct = cmap.type
                    cmap_types.append(id(ct))
                    # Our torsion correction maps need to go from 0 to 360
                    # degrees
                    grid = ct.grid.switch_range().T
                    m = force.addMap(ct.resolution, [x*ene_conv for x in grid])
                    cmap_map[id(ct)] = m
            # Now add in all of the cmaps
            for cmap in self.cmaps:
                id1, id2 = cmap.atom1.idx, cmap.atom2.idx
                id3, id4 = cmap.atom3.idx, cmap.atom4.idx
                id5, id6 = cmap.atom2.idx, cmap.atom3.idx
                id7, id8 = cmap.atom4.idx, cmap.atom5.idx
                force.addTorsion(cmap_map[id(cmap.type)],
                                 id1, id2, id3, id4, id5, id6, id7, id8)
            system.addForce(force)
        # Add nonbonded terms now
        if verbose: print('Adding nonbonded interactions...')
        force = mm.NonbondedForce()
        force.setForceGroup(self.NONBONDED_FORCE_GROUP)
        if not hasbox: # non-periodic
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

            # See if we need to use a switching function
            if switchDistance and nonbondedMethod is not ff.NoCutoff:
                # make sure it's legal
                if switchDistance >= nonbondedCutoff:
                    raise ValueError('switchDistance is too large compared '
                                     'to the cutoff!')
                    force.setUseSwitchingFunction(True)
                    force.setSwitchingDistance(switchDistance)

        else: # periodic
            # Set up box vectors (from inpcrd if available, or fall back to
            # prmtop definitions
            system.setDefaultPeriodicBoxVectors(*self.boxVectors)

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

            # See if we need to use a switching function
            if switchDistance and nonbondedMethod is not ff.NoCutoff:
                # make sure it's legal
                if switchDistance >= nonbondedCutoff:
                    raise ValueError('switchDistance is too large compared '
                                     'to the cutoff!')
                    force.setUseSwitchingFunction(True)
                    force.setSwitchingDistance(switchDistance)

            if ewaldErrorTolerance is not None:
                force.setEwaldErrorTolerance(ewaldErrorTolerance)

        # Add per-particle nonbonded parameters (LJ params)
        sigma_scale = 2**(-1/6) * 2
        if not has_nbfix_terms:
            for atm in self.atoms:
                force.addParticle(atm.charge,
                                  sigma_scale*atm.atom_type.rmin*length_conv,
                                  abs(atm.atom_type.epsilon*ene_conv))
        else:
            for atm in self.atoms:
                force.addParticle(atm.charge, 1.0, 0.0)
            # Now add the custom nonbonded force that implements NBFIX. First
            # thing we need to do is condense our number of types
            lj_idx_list = [0 for atom in self.atoms]
            lj_radii, lj_depths = [], []
            num_lj_types = 0
            lj_type_list = []
            for i, atom in enumerate(self.atoms):
                atom = atom.type
                if lj_idx_list[i]: continue # already assigned
                num_lj_types += 1
                lj_idx_list[i] = num_lj_types
                ljtype = (atom.rmin, atom.epsilon)
                lj_type_list.append(atom)
                lj_radii.append(atom.rmin)
                lj_depths.append(atom.epsilon)
                for j in xrange(i+1, len(self.atoms)):
                    atom2 = self.atoms[j].type
                    if lj_idx_list[j] > 0: continue # already assigned
                    if atom2 is atom:
                        lj_idx_list[j] = num_lj_types
                    elif not atom.nbfix:
                        # Only non-NBFIXed atom types can be compressed
                        ljtype2 = (atom2.rmin, atom2.epsilon)
                        if ljtype == ljtype2:
                            lj_idx_list[j] = num_lj_types
            # Now everything is assigned. Create the A-coefficient and
            # B-coefficient arrays
            acoef = [0 for i in xrange(num_lj_types*num_lj_types)]
            bcoef = acoef[:]
            for i in xrange(num_lj_types):
                for j in xrange(num_lj_types):
                    namej = lj_type_list[j].name
                    try:
                        rij, wdij, rij14, wdij14 = lj_type_list[i].nbfix[namej]
                    except KeyError:
                        rij = (lj_radii[i] + lj_radii[j]) * length_conv
                        wdij = sqrt(lj_depths[i] * lj_depths[j]) * ene_conv
                    else:
                        rij *= length_conv
                        wdij *= ene_conv
                    acoef[i+num_lj_types*j] = sqrt(wdij) * rij**6
                    bcoef[i+num_lj_types*j] = 2 * wdij * rij**6
            cforce = mm.CustomNonbondedForce('(a/r6)^2-b/r6; r6=r^6;'
                                             'a=acoef(type1, type2);'
                                             'b=bcoef(type1, type2)')
            cforce.addTabulatedFunction('acoef',
                    mm.Discrete2DFunction(num_lj_types, num_lj_types, acoef))
            cforce.addTabulatedFunction('bcoef',
                    mm.Discrete2DFunction(num_lj_types, num_lj_types, bcoef))
            cforce.addPerParticleParameter('type')
            cforce.setForceGroup(self.NONBONDED_FORCE_GROUP)
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
                raise ValueError('Unrecognized nonbonded method')
            if switchDistance and nonbondedMethod is not ff.NoCutoff:
                # make sure it's legal
                if switchDistance >= nonbondedCutoff:
                    raise ValueError('switchDistance is too large compared '
                                     'to the cutoff!')
                    cforce.setUseSwitchingFunction(True)
                    cforce.setSwitchingDistance(switchDistance)
            for i in lj_idx_list:
                cforce.addParticle((i - 1,)) # adjust for indexing from 0

        # Add 1-4 interactions
        excluded_atom_pairs = set() # save these pairs so we don't zero them out
        sigma_scale = 2**(-1/6)
        for tor in self.dihedrals:
            # First check to see if atoms 1 and 4 are already excluded because
            # they are 1-2 or 1-3 pairs (would happen in 6-member rings or
            # fewer). Then check that they're not already added as exclusions
            if tor.atom1 in tor.atom4.bond_partners: continue
            if tor.atom1 in tor.atom4.angle_partners: continue
            key = min((tor.atom1.idx, tor.atom4.idx),
                      (tor.atom4.idx, tor.atom1.idx))
            if key in excluded_atom_pairs: continue # multiterm...
            charge_prod = (tor.atom1.charge * tor.atom4.charge)
            epsilon = (sqrt(abs(tor.atom1.atom_type.epsilon_14) * ene_conv *
                            abs(tor.atom4.atom_type.epsilon_14) * ene_conv))
            sigma = (tor.atom1.atom_type.rmin_14+tor.atom4.atom_type.rmin_14)*(
                     length_conv * sigma_scale)
            force.addException(tor.atom1.idx, tor.atom4.idx,
                               charge_prod, sigma, epsilon)
            excluded_atom_pairs.add(
                    min((tor.atom1.idx, tor.atom4.idx),
                        (tor.atom4.idx, tor.atom1.idx))
            )

        # Add excluded atoms
        for atom in self.atoms:
            # Exclude all bonds and angles
            for atom2 in atom.bond_partners:
                if atom2.idx > atom.idx:
                    force.addException(atom.idx, atom2.idx, 0.0, 0.1, 0.0)
            for atom2 in atom.angle_partners:
                if atom2.idx > atom.idx:
                    force.addException(atom.idx, atom2.idx, 0.0, 0.1, 0.0)
            for atom2 in atom.dihedral_partners:
                if atom2.idx <= atom.idx: continue
                if ((atom.idx, atom2.idx) in excluded_atom_pairs):
                    continue
                force.addException(atom.idx, atom2.idx, 0.0, 0.1, 0.0)
        system.addForce(force)
        # If we needed a CustomNonbondedForce, map all of the exceptions from
        # the NonbondedForce to the CustomNonbondedForce
        if has_nbfix_terms:
            for i in xrange(force.getNumExceptions()):
                ii, jj, q, eps, sig = force.getExceptionParameters(i)
                cforce.addExclusion(ii, jj)
            system.addForce(cforce)

        # Add GB model if we're doing one
        if implicitSolvent is not None:
            if verbose: print('Adding GB parameters...')
            gb_parms = self._get_gb_params(implicitSolvent)
            gb_parms = convertParameters(gb_parms, str(implicitSolvent))

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
            for i, atom in enumerate(self.atoms):
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
            for bond in self.bonds:
                # Only take the ones with at least one hydrogen
                if (bond.atom1.atomic_number != 1 and
                    bond.atom2.atomic_number != 1):
                    continue
                atom1, atom2 = bond.atom1, bond.atom2
                if atom1.atomic_number == 1:
                    atom1, atom2 = atom2, atom1 # now atom2 is hydrogen for sure
                if atom1.atomic_number != 1:
                    transfer_mass = hydrogenMass - atom2.mass
                    new_mass1 = (system.getParticleMass(atom1.idx) -
                                 transfer_mass)
                    system.setParticleMass(atom2.idx, hydrogenMass)
                    system.setParticleMass(atom1.idx, new_mass1)
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
            return None

    @property
    def positions(self):
        """
        Return the cached positions or create new ones from the atoms
        """
        try:
            if len(self._positions) == len(self.atoms):
                return self._positions
        except AttributeError:
            pass

        try:
            self._positions = tuple([Vec3(a.xx, a.xy, a.xz)
                                     for a in self.atoms]) * u.angstroms
        except AttributeError:
            return None

        return self._positions

    @positions.setter
    def positions(self, stuff):
        """
        Replace the cached positions and the positions of each atom. If no units
        are applied to "stuff", it is assumed to be Angstroms.
        """
        if not u.is_quantity(stuff):
            # Assume this is Angstroms
            stuff *= u.angstroms

        # If we got a 1-D array, reshape it into an natom list of Vec3's
        if len(stuff) == len(self.atoms) * 3:
            stuff = [Vec3(stuff[i*3], stuff[i*3+1], stuff[i*3+2])
                     for i in xrange(len(self.atoms))]
        self._positions = stuff
        for atom, pos in zip(self.atoms, stuff):
            atom.xx, atom.xy, atom.xz = pos.value_in_unit(u.angstrom)

    @property
    def velocities(self):
        """ Same as for positions, but for velocities """
        try:
            if len(self._velocities) == len(self.atoms):
                return self._velocities
        except AttributeError:
            pass

        try:
            self._velocities = tuple([Vec3(a.vx, a.vy, a.vz)
                        for a in self.atoms]) * (u.angstroms/u.picosecond) 
        except AttributeError:
            return None

        return self._velocities

    @property
    def box_lengths(self):
        """ Return tuple of 3 units """
        if self.box_vectors is None:
            return None
        if u.is_quantity(self.box_vectors):
            vecs = self.box_vectors.value_in_unit(u.nanometers)
        elif u.is_quantity(self.box_vectors[0]):
            vecs = [x.value_in_unit(u.nanometers) for x in self.box_vectors]
        else:
            vecs = self.box_vectors
        a = sqrt(vecs[0][0]*vecs[0][0] + vecs[0][1]*vecs[0][1] +
                 vecs[0][2]*vecs[0][2])
        b = sqrt(vecs[1][0]*vecs[1][0] + vecs[1][1]*vecs[1][1] +
                 vecs[1][2]*vecs[1][2])
        c = sqrt(vecs[2][0]*vecs[2][0] + vecs[2][1]*vecs[2][1] +
                 vecs[2][2]*vecs[2][2])
        return (a, b, c) * u.nanometers

    @property
    def box_vectors(self):
        try:
            return self._box_vectors
        except AttributeError:
            return None

    @box_vectors.setter
    def box_vectors(self, thing):
        self._box_vectors = thing

    def setBox(self, a, b, c, alpha=90.0*u.degrees, beta=90.0*u.degrees,
               gamma=90.0*u.degrees):
        """
        Sets the periodic box boundary conditions.

        Parameters:
            - a, b, c (float) : Lengths of the unit cell
            - alpha, beta, gamma (float) : Angles between the unit cell vectors
        """
        self._box_vectors = _box_vectors_from_lengths_angles(a, b, c,
                                                            alpha, beta, gamma)
        # Now call "set_box" for the non-OpenMM object
        self.set_box(a.value_in_unit(u.angstroms),
                     b.value_in_unit(u.angstroms),
                     c.value_in_unit(u.angstroms),
                     alpha.value_in_unit(u.degrees),
                     beta.value_in_unit(u.degrees),
                     gamma.value_in_unit(u.degrees))
        # If we already have a _topology instance, then we have possibly changed
        # the existence of box information (whether or not this is a periodic
        # system), so delete any cached reference to a topology so it's
        # regenerated with updated information
        if hasattr(self, '_topology'):
            del self._topology

    @property
    def boxVectors(self):
        return self.box_vectors

    @boxVectors.setter
    def boxVectors(self, stuff):
        self.box_vectors = stuff
        # Now set the box
        lengths, angles = _box_vectors_to_lengths_angles(*stuff)
        self.box = [lengths[0], lengths[1], lengths[2],
                    angles[0], angles[1], angles[2]]

    # For consistency with OpenMM's API
    boxLengths = box_lengths
    def loadParameters(self, parmset):
        super(OpenMMCharmmPsfFile, self).load_parameters(parmset)
    def deleteCmap(self):
        """ Deletes the CMAP terms """
        del self.cmaps[:]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Routines that set the necessary radii lists based on a list of atoms with
# their connectivities

def _bondi_radii(atoms):
    """ Sets the bondi radii """
    radii = [0.0 for atom in atoms]
    for i, atom in enumerate(atoms):
        if atom.atomic_number == 6:
            radii[i] = 1.7
        elif atom.atomic_number == 1:
            radii[i] = 1.2
        elif atom.atomic_number == 7:
            radii[i] = 1.55
        elif atom.atomic_number == 8:
            radii[i] = 1.5
        elif atom.atomic_number == 9:
            radii[i] = 1.5
        elif atom.atomic_number == 14:
            radii[i] = 2.1
        elif atom.atomic_number == 15:
            radii[i] = 1.85
        elif atom.atomic_number == 16:
            radii[i] = 1.8
        elif atom.atomic_number == 17:
            radii[i] = 1.5
        else:
            radii[i] = 1.5
    return radii  # converted to nanometers above

def _mbondi_radii(atoms):
    """ Sets the mbondi radii """
    radii = [0.0 for atom in atoms]
    for i, atom in enumerate(atoms):
        # Radius of H atom depends on element it is bonded to
        if atom.atomic_number == 1:
            bondeds = list(atom.bond_partners)
            if bondeds[0].atomic_number in (6, 7): # C or N
                radii[i] = 1.3
            elif bondeds[0].atomic_number in (8, 16): # O or S
                radii[i] = 0.8
            else:
                radii[i] = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            radii[i] = 1.7
        # All other elements have fixed radii for all types/partners
        elif atom.atomic_number == 7:
            radii[i] = 1.55
        elif atom.atomic_number == 8:
            radii[i] = 1.5
        elif atom.atomic_number == 9:
            radii[i] = 1.5
        elif atom.atomic_number == 14:
            radii[i] = 2.1
        elif atom.atomic_number == 15:
            radii[i] = 1.85
        elif atom.atomic_number == 16:
            radii[i] = 1.8
        elif atom.atomic_number == 17:
            radii[i] = 1.5
        else:
            radii[i] = 1.5
    return radii  # converted to nanometers above

def _mbondi2_radii(atoms):
    """ Sets the mbondi2 radii """
    radii = [0.0 for atom in atoms]
    for i, atom in enumerate(atoms):
        # Radius of H atom depends on element it is bonded to
        if atom.atomic_number == 1:
            if atom.bond_partners[0].atomic_number == 7:
                radii[i] = 1.3
            else:
                radii[i] = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            radii[i] = 1.7
        # All other elements have fixed radii for all types/partners
        elif atom.atomic_number == 7:
            radii[i] = 1.55
        elif atom.atomic_number == 8:
            radii[i] = 1.5
        elif atom.atomic_number == 9:
            radii[i] = 1.5
        elif atom.atomic_number == 14:
            radii[i] = 2.1
        elif atom.atomic_number == 15:
            radii[i] = 1.85
        elif atom.atomic_number == 16:
            radii[i] = 1.8
        elif atom.atomic_number == 17:
            radii[i] = 1.5
        else:
            radii[i] = 1.5
    return radii  # Converted to nanometers above

def _mbondi3_radii(atoms):
    """ Sets the mbondi3 radii """
    radii = _mbondi2_radii(atoms)
    for i, atom in enumerate(atoms):
        # Adjust OE (GLU), OD (ASP) and HH/HE (ARG)
        if atom.residue.name in ('GLU', 'ASP'):
            if atom.name.startswith('OE') or atom.name.startswith('OD'):
                radii[i] = 1.4
        elif atom.residue.name == 'ARG':
            if atom.name.startswith('HH') or atom.name.startswith('HE'):
                radii[i] = 1.17
        # Adjust carboxylate O radii on C-termini
        if atom.name == 'OXT':
            radii[i] = 1.4
    return radii  # Converted to nanometers above

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _box_vectors_from_lengths_angles(a, b, c, alpha, beta, gamma):
    """
    This method takes the lengths and angles from a unit cell and creates unit
    cell vectors.

    Parameters:
        - a (unit, dimension length): Length of the first vector
        - b (unit, dimension length): Length of the second vector
        - c (unit, dimension length): Length of the third vector
        - alpha (float): Angle between b and c
        - beta (float): Angle between a and c
        - gamma (float): Angle between a and b

    Returns:
        Tuple of box vectors (as Vec3 instances)

    Notes:
        If all angles given are less than 6.28 (2*pi), we assume the angles are
        given in radians since a unit cell with ALL 3 angles < 2 pi is HIGHLY
        unusual.  If, however, you want a box whose angles are all smaller than
        2 pi degrees, first convert to radians.
   """
    if not (u.is_quantity(a) and u.is_quantity(b) and u.is_quantity(c)):
        raise TypeError('a, b, and c must be units of dimension length')
    if u.is_quantity(alpha): alpha = alpha.value_in_unit(u.degree)
    if u.is_quantity(beta): beta = beta.value_in_unit(u.degree)
    if u.is_quantity(gamma): gamma = gamma.value_in_unit(u.degree)
    a = a.value_in_unit(u.angstrom)
    b = b.value_in_unit(u.angstrom)
    c = c.value_in_unit(u.angstrom)

    # Convert to radians if it's not already there
    if not (alpha <= 2 * pi and beta <= 2 * pi and gamma <= 2 * pi):
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class OpenMMCharmmCrdFile(CharmmCrdFile):
    """ An OpenMM-enabled version of CharmmCrdFile """

    @property
    def positions(self):
        try:
            return self._positions
        except AttributeError:
            pass
        self._positions = []
        for i in xrange(self.natom):
            i3 = i * 3
            crd = Vec3(self.coords[i3], self.coords[i3+1], self.coords[i3+2])
            self._positions.append(crd)
        self._positions = self._positions * u.angstroms
        return self._positions
    
    @positions.setter
    def positions(self, stuff):
        self._positions = stuff
        for i, crd in enumerate(stuff.value_in_unit(u.angstroms)):
            i3 = i * 3
            self.coords[i3  ] = crd[0]
            self.coords[i3+1] = crd[1]
            self.coords[i3+2] = crd[2]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class OpenMMCharmmRstFile(CharmmRstFile):
    """ An OpenMM-enabled version of CharmmRstFile """

    @property
    def positions(self):
        try:
            return self._positions
        except AttributeError:
            pass
        self._positions = []
        for i in xrange(self.natom):
            i3 = i * 3
            crd = Vec3(self.coords[i3], self.coords[i3+1], self.coords[i3+2])
            self._positions.append(crd)
        self._positions = self._positions * u.angstroms
        return self._positions

    @positions.setter
    def positions(self, stuff):
        self._positions = stuff
        for i, crd in enumerate(stuff.value_in_unit(u.angstroms)):
            i3 = i * 3
            self.coords[i3  ] = crd[0]
            self.coords[i3+1] = crd[1]
            self.coords[i3+2] = crd[2]

    @property
    def velocities(self):
        try:
            return self._velocities
        except AttributeError:
            pass
        self._velocities = []
        for i in xrange(self.natom):
            i3 = i * 3
            vel = Vec3(self.vels[i3], self.vels[i3+1], self.vels[i3+2])
            self._velocities.append(vel)
        self._velocities = self._velocities * u.angstroms / u.picoseconds
        return self._velocities

    @velocities.setter
    def velocities(self, stuff):
        self._velocities = stuff
        for i, vel in enumerate(stuff.value_in_unit(u.angstroms/u.picoseconds)):
            i3 = i * 3
            self.vels[i3  ] = vel[0]
            self.vels[i3+1] = vel[1]
            self.vels[i3+2] = vel[2]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _box_vectors_to_lengths_angles(a, b, c):
    """
    This converts 3 periodic box vectors, converting to units of angstroms,
    and returns a 2-element tuple of box lengths and box angles. Both will be
    3-element tuples themselves
    """
    RADDEG = 180.0 / pi
    a = a.value_in_unit(u.angstrom)
    b = b.value_in_unit(u.angstrom)
    c = c.value_in_unit(u.angstrom)
    # Get the lengths of the unit vectors (now in angstroms)
    la = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
    lb = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2])
    lc = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2])
    # Angles
    alpha = acos((b[0]*c[0] + b[1]*c[1] + b[2]*c[2]) / (lb * lc)) * RADDEG
    beta = acos((a[0]*c[0] + a[1]*c[1] + a[2]*c[2]) / (la * lc)) * RADDEG
    gamma = acos((a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / (la * lb)) * RADDEG

    return ( (la, lb, lc), (alpha, beta, gamma) )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _OpenMMRestartFile(object):
    """
    This is a private class that is basically a container for a standardized
    "Restart" file for OpenMM simulations.  Basically, it stores the
    coordinates, velocities, and box vectors if present and loads it into this
    data structure if they're present

    Attributes
    ----------
        - boxVectors (tuple of 3 Vec3's) : The periodic box vectors
        - positions (natom-length list of Vec3's) : The coordinates
        - velocities (natom-length list of Vec3's) : The velocities

    The format of the file is:

    - BEGIN -
    OpenMM Restart File
    natom=# vels=[True|False] box=[True|False]
    <x1> <y1> <z1>
    <x2> <y2> <z2>
    ...
    <vx1> <vy1> <vz1>
    <vx2> <vy2> <vz2>
    ...
    <boxx1> <boxy1> <boxz1>
    <boxx2> <boxy2> <boxz2>
    <boxx3> <boxy3> <boxz3>
    - END -

    The velocities and box are only present if vels=True and box=True is found
    on the second line of the file. The units are nanometers for distance and
    nanometers/picosecond for velocities
    """
    infolinere = re.compile(r'natom=(\d+) vels=(True|False) box=(True|False)')
    tagline = 'OpenMM Restart File\n'

    @classmethod
    def read(cls, fname):
        """ Read a restart file """
        positions = []
        velocities = []
        boxVectors = None
        with open(fname, 'r') as f:
            if f.readline() != cls.tagline:
                raise ValueError('%s is not a recognized OpenMM Restart' %
                                 fname)
            rematch = cls.infolinere.match(f.readline())
            if not rematch:
                raise RuntimeError('Corrupt or incompatible OpenMM Restart '
                                   'file. Cannot parse header.')
            natom, hasvels, hasbox = rematch.groups()
            # Because the infolinere is strict about matching, the following is
            # safe
            natom = int(natom)
            hasvels = eval(hasvels)
            hasbox = eval(hasbox)

            i = 0
            while i < natom:
                i += 1
                line = f.readline()
                if not line:
                    raise RuntimeError('Unexpected EOF in restart file!')
                x, y, z = line.split()
                positions.append(Vec3(float(x), float(y), float(z)) *
                                      u.nanometers)
            if hasvels:
                i = 0
                while i < natom:
                    i += 1
                    line = f.readline()
                    if not line:
                        raise RuntimeError('Unexpected EOF in restart file!')
                    x, y, z = line.split()
                    velocities.append(Vec3(float(x), float(y), float(z)) *
                                           u.nanometers / u.picoseconds)
            if hasbox:
                box1 = f.readline().split()
                box2 = f.readline().split()
                box3 = f.readline().split()
                boxVectors = (Vec3(float(box1[0]), float(box1[1]),
                                   float(box1[2])) * u.nanometers,
                              Vec3(float(box2[0]), float(box2[1]),
                                   float(box2[2])) * u.nanometers,
                              Vec3(float(box3[0]), float(box3[1]),
                                   float(box3[2])) * u.nanometers
                )
            return cls(positions, velocities, boxVectors)

    def __init__(self, positions, velocities, boxVectors):
        self.natom = len(positions)
        self.positions = positions
        if velocities and len(velocities) != self.natom:
            raise ValueError('Bad number of velocities')
        self.velocities = velocities
        if boxVectors is not None and len(boxVectors) != 3:
            raise ValueError('Bad number of box vectors')
        self.boxVectors = boxVectors

    def write(self, fname):
        """ Writes the information to a new restart file """
        with open(fname, 'w') as f:
            f.write(self.tagline)
            f.write('natom=%d ' % self.natom)
            if self.velocities:
                f.write('vels=True ')
            else:
                f.write('vels=False ')
            if self.boxVectors:
                f.write('box=True\n')
            else:
                f.write('box=False\n')
            for crd in self.positions:
                crd = crd.value_in_unit(u.nanometers)
                f.write('%s %s %s\n' % tuple(crd))
            if self.velocities:
                for vel in self.velocities:
                    vel = vel.value_in_unit(u.nanometers/u.picoseconds)
                    f.write('%s %s %s\n' % tuple(vel))
            if self.boxVectors:
                for box in self.boxVectors:
                    box = box.value_in_unit(u.nanometers)
                    f.write('%s %s %s\n' % tuple(box))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def read_restart(fname):
    """ Reads a restart file.

    Returns
    -------
        restart object with the following properties:
            - positions
            - velocities (empty list if no velocities are present)
            - boxVectors (None if box is not present)
    """
    return _OpenMMRestartFile.read(fname)

def write_restart(fname, sim):
    """ Writes an OpenMM Restart file

    Parameters
    ----------
        - fname (str) : File name of restart file to write
        - sim (Simulation) : OpenMM Simulation object
    """
    hasbox = sim.topology.getUnitCellDimensions() is not None
    state = sim.context.getState(getPositions=True, getVelocities=True,
                                 enforcePeriodicBox=hasbox)
    _OpenMMRestartFile(state.getPositions(), state.getVelocities(),
                       state.getPeriodicBoxVectors()).write(fname)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
