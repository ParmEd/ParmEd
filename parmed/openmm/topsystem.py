"""
Convert an OpenMM Topology into a Structure instance, optionally filling in
parameters from a System
"""
from __future__ import absolute_import, division, print_function

import warnings
from collections import defaultdict

import numpy as np

from .. import unit as u
from ..exceptions import OpenMMWarning
from ..formats import load_file
from ..geometry import box_vectors_to_lengths_and_angles
from ..periodic_table import Element
from ..structure import Structure
from ..topologyobjects import (Angle, AngleType, Atom, AtomType, Bond, BondType, Cmap, CmapType,
                               Dihedral, DihedralType, ExtraPoint, Improper, ImproperType,
                               NonbondedException, NonbondedExceptionType, RBTorsionType,
                               UreyBradley)
from ..utils.decorators import needs_openmm
from ..utils.six import integer_types, iteritems, string_types
from ..utils.six.moves import range

__all__ = ['load_topology']


@needs_openmm
def load_topology(topology, system=None, xyz=None, box=None):
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
    xyz : str or array of float
        Name of a file containing coordinate information or an array of
        coordinates. If file has unit cell information, it also uses that
        information unless ``box`` (below) is also specified
    box : array of 6 floats
        Unit cell dimensions

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
    processed and will result in an unknown functional form.

    If an OpenMM Atom.id attribute is populated by a non-integer, it will be
    used to name the corresponding ParmEd AtomType object.
    """
    import simtk.openmm as mm
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
                    atype = a.id if not isinstance(a.id, integer_types) else ''
                    atom = Atom(atomic_number=a.element.atomic_number,
                                name=a.name, mass=a.element.mass, type=atype)
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

    loaded_box = False

    if xyz is not None:
        if isinstance(xyz, string_types):
            xyz = load_file(xyz, skip_bonds=True)
            struct.coordinates = xyz.coordinates
            if struct.box is not None:
                if xyz.box is not None:
                    loaded_box = True
                    struct.box = xyz.box
        else:
            struct.coordinates = xyz

    if box is not None:
        loaded_box = True
        struct.box = box

    if struct.box is not None:
        struct.box = np.asarray(struct.box)

    if system is None:
        return struct

    if isinstance(system, string_types):
        system = load_file(system)

    if not isinstance(system, mm.System):
        raise TypeError('system must be an OpenMM System object or serialized '
                        'XML of an OpenMM System object')

    # We have a system, try to extract parameters from it
    if len(struct.atoms) != system.getNumParticles():
        raise TypeError('Topology and System have different numbers of atoms (%d vs. %d)' %
                        (len(struct.atoms), system.getNumParticles()))

    processed_forces = set()
    ignored_forces = (mm.CMMotionRemover, mm.AndersenThermostat, mm.MonteCarloBarostat,
                      mm.MonteCarloAnisotropicBarostat, mm.MonteCarloMembraneBarostat,
                      mm.CustomExternalForce, mm.GBSAOBCForce, mm.CustomGBForce)

    if system.usesPeriodicBoundaryConditions():
        if not loaded_box:
            vectors = system.getDefaultPeriodicBoxVectors()
            leng, ang = box_vectors_to_lengths_and_angles(*vectors)
            leng = leng.value_in_unit(u.angstroms)
            ang = ang.value_in_unit(u.degrees)
            struct.box = np.asarray([leng[0], leng[1], leng[2], ang[0], ang[1], ang[2]])
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
                warnings.warn('Unknown functional form of CustomTorsionForce', OpenMMWarning)
        elif isinstance(force, mm.CMAPTorsionForce):
            _process_cmap(struct, force)
        elif isinstance(force, mm.NonbondedForce):
            _process_nonbonded(struct, force)
        elif isinstance(force, ignored_forces):
            continue
        else:
            struct.unknown_functional = True
            warnings.warn('Unsupported Force type %s' % type(force).__name__, OpenMMWarning)
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
        warnings.warn('Adding what seems to be Urey-Bradley terms before ' # pragma: no cover
                      'Angles. This is unexpected, but the parameters will '
                      'all be present in one form or another.', OpenMMWarning)
    typemap = dict()
    for ii in range(force.getNumBonds()):
        i, j, req, k = force.getBondParameters(ii)
        ai, aj = struct.atoms[i], struct.atoms[j]
        key = (req._value, k._value)
        if struct.angles and ai not in aj.angle_partners:
            warnings.warn('Adding what seems to be Urey-Bradley terms, but ' # pragma: no cover
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
            key = (c0._value, c1._value, c2._value, c3._value, c4._value, c5._value)
            f = 1                          # pragma: no cover
        except AttributeError:             # pragma: no cover
            key = (c0, c1, c2, c3, c4, c5) # pragma: no cover
            f = u.kilojoules_per_mole      # pragma: no cover
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
    if ';' in eqn: # Just look at the first segment of the equation if there are multiple
        eqn = eqn[:eqn.index(';')]
    # Don't try to be fancy with regexes for fear of making a possible mistake.
    # ParmEd and OpenMM use only these two eqns for the improper torsions:
    #  k*(theta-theta0)^2 vs. 0.5*k*(theta-theta0)^2
    # So only recognize the above 2 forms
    if eqn not in ('0.5*k*(theta-theta0)^2', 'k*(theta-theta0)^2', 'k*dtheta_torus^2',
                   '0.5*k*dtheta_torus^2'):
        return False
    if eqn.startswith('0.5'):
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
            typ = CmapType(size, grid)                       # pragma: no cover
        else:
            typ = CmapType(size, grid*u.kilojoules_per_mole) # pragma: no cover
        cmap_types.append(typ)
        typ.grid = typ.grid.T.switch_range()
        typ.used = False
    # Add all cmaps
    for ii in range(force.getNumTorsions()):
        mapidx, ii, ij, ik, il, ji, jj, jk, jl = force.getTorsionParameters(ii)
        if ij != ji or ik != jj or il != jk:
            warnings.warn('Non-continuous CMAP torsions detected. Not ' # pragma: no cover
                          'supported.', OpenMMWarning)
            continue # pragma: no cover
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
        atype_name = (atom.type if atom.type != ''
                      else Element[atom.atomic_number])
        key = (atype_name, sig._value, eps._value)
        if key in typemap:
            atom_type = typemap[key]
        else:
            if atom.type == '':
                element_typemap[atype_name] += 1
                atype_name = '%s%d' % (atype_name, element_typemap[atype_name])
            typemap[key] = atom_type = AtomType(atype_name, None, atom.mass,
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
            if atom is not a2:
                bond_graph_exceptions[atom].add(a2)
            for a3 in a2.bond_partners:
                if a3 is atom: continue
                if atom is not a3:
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
                raise ValueError("Can't scale charge product 0 to match %s" % q)
            chgscale = None
        nbtype = NonbondedExceptionType(sig*2**(1/6), eps, chgscale)
        struct.adjusts.append(NonbondedException(ai, aj, type=nbtype))
        struct.adjust_types.append(nbtype)
    struct.adjust_types.claim()
    # Go through all adjust_types and replace any chgscale values that are None with the first
    # non-None value present. If all are None, set all to 1.0. This way we maximize the likelihood
    # that all generated systems have 1 scaling factor which makes it easier to translate to other
    # formats that may not support multiple scaling factors in electrostatic scaling (like GROMACS)
    first_scaling_factor = 1.0
    for adjust_type in struct.adjust_types:
        if adjust_type.chgscale is not None:
            first_scaling_factor = adjust_type.chgscale
            break
    # Now go through and set all Nones to first_scaling_factor
    for adjust_type in struct.adjust_types:
        if adjust_type.chgscale is None:
            adjust_type.chgscale = first_scaling_factor

    # Check that all of our exceptions are accounted for
    for ai, exceptions in iteritems(bond_graph_exceptions):
        if exceptions - explicit_exceptions[ai]:
            struct.unknown_functional = True
            warnings.warn('Detected incomplete exceptions. Not supported.',
                          OpenMMWarning)
            break
