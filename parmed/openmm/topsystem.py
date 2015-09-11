"""
Convert an OpenMM Topology into a Structure instance, optionally filling in
parameters from a System
"""
from __future__ import division, print_function, absolute_import

__all__ = ['load_topology']

from parmed.exceptions import OpenMMWarning
from parmed.geometry import box_vectors_to_lengths_and_angles
from parmed.periodic_table import Element
from parmed.structure import Structure
from parmed.topologyobjects import (Atom, Bond, BondType, Angle, AngleType,
        Dihedral, DihedralType, Improper, ImproperType, AtomType, ExtraPoint,
        UreyBradley, NonbondedExceptionType, NonbondedException, Cmap, CmapType,
        RBTorsionType)
from parmed import unit as u
from parmed.utils.decorators import needs_openmm
from parmed.utils.six import iteritems, string_types
from parmed.utils.six.moves import range
from collections import defaultdict
try:
    import numpy as np
    def create_array(array):
        return np.asarray(array)
except ImportError:
    np = None
    def create_array(array):
        return array
try:
    import simtk.openmm as mm
except ImportError:
    pass
import warnings

@needs_openmm
def load_topology(topology, system=None):
    """
    Creates a :class:`parmed.structure.Structure` instance from an OpenMM
    Topology, optionally filling in parameters from a System

    Parameters
    ----------
    topology : :class:`simtk.openmm.app.Topology`
        The Topology instance with the list of atoms and bonds for this system
    system : :class:`simtk.openmm.System` or str, optional
        If provided, parameters from this System will be applied to the
        Structure. If a string is given, it will be interpreted as the file name
        of an XML-serialized System, and it will be deserialized into a System
        before used to supply parameters

    Returns
    -------
    struct : :class:`Structure <parmed.structure.Structure>`
        The structure from the provided topology

    Raises
    ------
    OpenMMWarning if parameters are found that cannot be interpreted or
    processed by ParmEd

    TypeError if there are any mismatches between the provided topology and
    system (e.g., they have different numbers of atoms)

    IOError if system is a string and it is not an existing file

    Notes
    -----
    Due to its flexibility with CustomForces, it is entirely possible that the
    functional form of the potential will be unknown to ParmEd. This function
    will try to use the energy expression to identify supported potential types
    that are implemented as CustomForce objects. In particular, quadratic
    improper torsions, when recognized, will be extracted.

    Other CustomForces, including the CustomNonbondedForce used to implement
    NBFIX (off-diagonal L-J modifications) and the 12-6-4 potential, will not be
    processed and will result in an unknown functional form
    """
    struct = Structure()
    atommap = dict()
    for c in topology.chains():
        chain = c.id
        for r in c.residues():
            residue = r.name
            resid = r.index
            for a in r.atoms():
                if a.element is None:
                    atom = ExtraPoint(name=a.name)
                else:
                    atom = Atom(atomic_number=a.element.atomic_number,
                                name=a.name, mass=a.element.mass)
                struct.add_atom(atom, residue, resid, chain)
                atommap[a] = atom
    for a1, a2 in topology.bonds():
        struct.bonds.append(Bond(atommap[a1], atommap[a2]))

    vectors = topology.getPeriodicBoxVectors()
    if vectors is not None:
        leng, ang = box_vectors_to_lengths_and_angles(*vectors)
        leng = leng.value_in_unit(u.angstroms)
        ang = ang.value_in_unit(u.degrees)
        struct.box = [leng[0], leng[1], leng[2], ang[0], ang[1], ang[2]]

    if struct.box is not None:
        struct.box = create_array(struct.box)

    if system is None:
        return struct

    if isinstance(system, string_types):
        with open(system, 'r') as f:
            system = mm.XmlSerializer.deserialize(f.read())

    # We have a system, try to extract parameters from it
    if len(struct.atoms) != system.getNumParticles():
        raise TypeError('Topology and System have different numbers of atoms '
                '(%d vs. %d)' % (len(struct.atoms), system.getNumParticles()))

    processed_forces = set()
    ignored_forces = (mm.CMMotionRemover, mm.AndersenThermostat,
                      mm.MonteCarloBarostat, mm.MonteCarloAnisotropicBarostat,
                      mm.MonteCarloMembraneBarostat, mm.CustomExternalForce,
                      mm.GBSAOBCForce, mm.CustomGBForce)

    if system.usesPeriodicBoundaryConditions():
        vectors = system.getDefaultPeriodicBoxVectors()
        leng, ang = box_vectors_to_lengths_and_angles(*vectors)
        leng = leng.value_in_unit(u.angstroms)
        ang = ang.value_in_unit(u.degrees)
        struct.box = create_array(
                [leng[0], leng[1], leng[2], ang[0], ang[1], ang[2]]
        )
    else:
        struct.box = None

    for force in system.getForces():
        if isinstance(force, mm.HarmonicBondForce):
            if mm.HarmonicBondForce in processed_forces:
                # Try to process this HarmonicBondForce as a Urey-Bradley term
                _process_urey_bradley(struct, force)
            else:
                _process_bond(struct, force)
        elif isinstance(force, mm.HarmonicAngleForce):
            _process_angle(struct, force)
        elif isinstance(force, mm.PeriodicTorsionForce):
            _process_dihedral(struct, force)
        elif isinstance(force, mm.RBTorsionForce):
            _process_rbtorsion(struct, force)
        elif isinstance(force, mm.CustomTorsionForce):
            if not _process_improper(struct, force):
                struct.unknown_functional = True
                warnings.warn('Unknown functional form of CustomTorsionForce',
                              OpenMMWarning)
        elif isinstance(force, mm.CMAPTorsionForce):
            _process_cmap(struct, force)
        elif isinstance(force, mm.NonbondedForce):
            _process_nonbonded(struct, force)
        elif isinstance(force, ignored_forces):
            continue
        else:
            struct.unknown_functional = True
            warnings.warn('Unsupported Force type %s' % type(force).__name__,
                          OpenMMWarning)
        processed_forces.add(type(force))

    return struct

def _process_bond(struct, force):
    """ Adds bond parameters to the structure """
    typemap = dict()
    for ii in range(force.getNumBonds()):
        i, j, req, k = force.getBondParameters(ii)
        ai, aj = struct.atoms[i], struct.atoms[j]
        key = (req._value, k._value)
        if key in typemap:
            bond_type = typemap[key]
        else:
            bond_type = BondType(k*0.5, req)
            typemap[key] = bond_type
            struct.bond_types.append(bond_type)
        if aj in ai.bond_partners:
            for bond in ai.bonds:
                if aj in bond:
                    break
            else:
                raise RuntimeError('aj in ai.bond_partners, but couldn\'t find '
                                   'that bond!')
            bond.type = bond_type
        else:
            struct.bonds.append(Bond(ai, aj, type=bond_type))
    struct.bond_types.claim()

def _process_angle(struct, force):
    """ Adds angle parameters to the structure """
    typemap = dict()
    for ii in range(force.getNumAngles()):
        i, j, k, theteq, frc_k = force.getAngleParameters(ii)
        key = (theteq._value, frc_k._value)
        ai, aj, ak = struct.atoms[i], struct.atoms[j], struct.atoms[k]
        if key in typemap:
            angle_type = typemap[key]
        else:
            angle_type = AngleType(frc_k*0.5, theteq)
            typemap[key] = angle_type
            struct.angle_types.append(angle_type)
        struct.angles.append(Angle(ai, aj, ak, type=angle_type))
    struct.angle_types.claim()

def _process_urey_bradley(struct, force):
    """ Adds Urey-Bradley parameters to the structure """
    if not struct.angles:
        warnings.warn('Adding what seems to be Urey-Bradley terms before '
                      'Angles. This is unexpected, but the parameters will '
                      'all be present in one form or another.', OpenMMWarning)
    typemap = dict()
    for ii in range(force.getNumBonds()):
        i, j, req, k = force.getBondParameters(ii)
        ai, aj = struct.atoms[i], struct.atoms[j]
        key = (req._value, k._value)
        if struct.angles and ai not in aj.angle_partners:
            warnings.warn('Adding what seems to be Urey-Bradley terms, but '
                          'atoms %d and %d do not appear to be angled to each '
                          'other. Parameters will all be present, but may not '
                          'be in expected places.' % (ai.idx, aj.idx),
                          OpenMMWarning)
        if key in typemap:
            urey_type = typemap[key]
        else:
            urey_type = BondType(k*0.5, req)
            typemap[key] = urey_type
            struct.urey_bradley_types.append(urey_type)
        struct.urey_bradleys.append(UreyBradley(ai, aj, type=urey_type))
    struct.urey_bradley_types.claim()

def _process_dihedral(struct, force):
    """ Adds periodic torsions to the structure """
    typemap = dict()
    for ii in range(force.getNumTorsions()):
        i, j, k, l, per, phase, phi_k = force.getTorsionParameters(ii)
        ai, aj = struct.atoms[i], struct.atoms[j]
        ak, al = struct.atoms[k], struct.atoms[l]
        key = (per, phase._value, phi_k._value)
        if key in typemap:
            dihed_type = typemap[key]
        else:
            dihed_type = DihedralType(phi_k, per, phase)
            typemap[key] = dihed_type
            struct.dihedral_types.append(dihed_type)
        improper = (ai in ak.bond_partners and aj in ak.bond_partners and
                    al in ak.bond_partners)
        struct.dihedrals.append(Dihedral(ai, aj, ak, al, improper=improper,
                                         type=dihed_type))
    struct.dihedral_types.claim()

def _process_rbtorsion(struct, force):
    """ Adds Ryckaert-Bellemans torsions to the structure """
    typemap = dict()
    for ii in range(force.getNumTorsions()):
        i, j, k, l, c0, c1, c2, c3, c4, c5 = force.getTorsionParameters(ii)
        ai, aj = struct.atoms[i], struct.atoms[j]
        ak, al = struct.atoms[k], struct.atoms[l]
        # TODO -- Fix this when OpenMM is fixed
        try:
            key = (c0._value, c1._value, c2._value,
                   c3._value, c4._value, c5._value)
            f = 1
        except AttributeError:
            key = (c0, c1, c2, c3, c4, c5)
            f = u.kilojoules_per_mole
        if key in typemap:
            dihed_type = typemap[key]
        else:
            dihed_type = RBTorsionType(c0*f, c1*f, c2*f, c3*f, c4*f, c5*f)
            typemap[key] = dihed_type
            struct.rb_torsion_types.append(dihed_type)
        struct.rb_torsions.append(Dihedral(ai, aj, ak, al, type=dihed_type))
    struct.rb_torsion_types.claim()

def _process_improper(struct, force):
    """ Processes a CustomTorsionForce and looks at the energy expression to see
    if it's a quadratic improper torsion. Then adds the parameters if applicable

    Returns
    -------
    is_improper : bool
        Returns True if the energy expression is recognized as a quadratic
        improper, and False otherwise
    """
    eqn = force.getEnergyFunction().replace(' ', '')
    # Don't try to be fancy with regexes for fear of making a possible mistake.
    # ParmEd and OpenMM use only these two eqns for the improper torsions:
    #  k*(theta-theta0)^2 vs. 0.5*k*(theta-theta0)^2
    # So only recognize the above 2 forms
    if eqn not in ('0.5*k*(theta-theta0)^2', 'k*(theta-theta0)^2'):
        return False
    if eqn == '0.5*k*(theta-theta0)^2':
        fac = 0.5
    else:
        fac = 1
    typemap = dict()
    for ii in range(force.getNumTorsions()):
        # All formulations put k first and equilibrium angle second, in radians
        i, j, k, l, (psi_k, psi_eq) = force.getTorsionParameters(ii)
        ai, aj = struct.atoms[i], struct.atoms[j]
        ak, al = struct.atoms[k], struct.atoms[l]
        key = (psi_k, psi_eq)
        if key in typemap:
            imp_type = typemap[key]
        else:
            imp_type = ImproperType(psi_k*fac*u.kilojoule_per_mole/u.radian**2,
                                    psi_eq*u.radian)
            typemap[key] = imp_type
            struct.improper_types.append(imp_type)
        struct.impropers.append(Improper(ai, aj, ak, al, type=imp_type))
    struct.improper_types.claim()
    return True

def _process_cmap(struct, force):
    """ Adds CMAPs to the structure """
    # store the list of cmap types
    cmap_types = []
    for ii in range(force.getNumMaps()):
        size, grid = force.getMapParameters(ii)
        # Future-proof in case units start getting added to these maps
        if u.is_quantity(grid):
            typ = CmapType(size, grid)
        else:
            typ = CmapType(size, grid*u.kilojoules_per_mole)
        cmap_types.append(typ)
        typ.grid = typ.grid.T.switch_range()
        typ.used = False
    # Add all cmaps
    for ii in range(force.getNumTorsions()):
        mapidx, ii, ij, ik, il, ji, jj, jk, jl = force.getTorsionParameters(ii)
        if ij != ji or ik != jj or il != jk:
            warnings.warn('Non-continuous CMAP torsions detected. Not '
                          'supported.', OpenMMWarning)
            continue
        ai, aj, ak = struct.atoms[ii], struct.atoms[ij], struct.atoms[ik]
        al, am = struct.atoms[il], struct.atoms[jl]
        cmap_type = cmap_types[mapidx]
        cmap_type.used = True
        struct.cmaps.append(Cmap(ai, aj, ak, al, am, type=cmap_type))
    for cmap_type in cmap_types:
        if cmap_type.used:
            struct.cmap_types.append(cmap_type)
    struct.cmap_types.claim()

def _process_nonbonded(struct, force):
    """ Adds nonbonded parameters to the structure """
    typemap = dict()
    element_typemap = defaultdict(int)
    assert force.getNumParticles() == len(struct.atoms), "Atom # mismatch"
    for i in range(force.getNumParticles()):
        atom = struct.atoms[i]
        chg, sig, eps = force.getParticleParameters(i)
        atype_name = Element[atom.atomic_number]
        key = (atype_name, sig._value, eps._value)
        if key in typemap:
            atom_type = typemap[key]
        else:
            element_typemap[atype_name] += 1
            atype_name = '%s%d' % (atype_name, element_typemap[atype_name])
            atom_type = AtomType(atype_name, None, atom.mass,
                                 atom.atomic_number)
        atom.charge = chg.value_in_unit(u.elementary_charge)
        rmin = sig.value_in_unit(u.angstroms) * 2**(1/6) / 2 # to rmin/2
        eps = eps.value_in_unit(u.kilocalories_per_mole)
        atom_type.set_lj_params(eps, rmin)
        atom.atom_type = atom_type
        atom.type = atom_type.name

    explicit_exceptions = defaultdict(set)
    bond_graph_exceptions = defaultdict(set)
    for atom in struct.atoms:
        for a2 in atom.bond_partners:
            bond_graph_exceptions[atom].add(a2)
            for a3 in a2.bond_partners:
                bond_graph_exceptions[atom].add(a3)

    # TODO should we compress exception types?
    for ii in range(force.getNumExceptions()):
        i, j, q, sig, eps = force.getExceptionParameters(ii)
        q = q.value_in_unit(u.elementary_charge**2)
        sig = sig.value_in_unit(u.angstrom)
        eps = eps.value_in_unit(u.kilocalorie_per_mole)
        ai, aj = struct.atoms[i], struct.atoms[j]
        if q == 0 and (sig == 0 or eps == 0):
            explicit_exceptions[ai].add(aj)
            explicit_exceptions[aj].add(ai)
            continue
        try:
            chgscale = q / (ai.charge * aj.charge)
        except ZeroDivisionError:
            if q != 0:
                raise TypeError('Cannot scale charge product of 0 to match '
                                '%s' % q)
            chgscale = 1
        nbtype = NonbondedExceptionType(sig*2**(1/6), eps, chgscale)
        struct.adjusts.append(NonbondedException(ai, aj, type=nbtype))
        struct.adjust_types.append(nbtype)
    struct.adjust_types.claim()

    # Check that all of our exceptions are accounted for
    for ai, exceptions in iteritems(bond_graph_exceptions):
        if exceptions - explicit_exceptions[ai]:
            warnings.warn('Detected incomplete exceptions. Not supported.',
                          OpenMMWarning)
            struct.unknown_functional = True
