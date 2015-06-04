"""
Convert an OpenMM Topology into a Structure instance, optionally filling in
parameters from a System
"""
from __future__ import division, print_function, absolute_import

__all__ = ['load_topology']

try:
    import simtk.openmm as mm
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False

from chemistry.structure import Structure
from chemistry.topologyobjects import (Atom, Bond, BondType, Angle, AngleType,
        Dihedral, DihedralType, Improper, ImproperType, Residue, AtomType,
        ExtraPoint)
from chemistry import unit as u

def load_topology(topology, system=None):
    """
    Creates a :class:`chemistry.structure.Structure` instance from an OpenMM
    Topology, optionally filling in parameters from a System

    Parameters
    ----------
    topology : :class:`simtk.openmm.app.Topology`
        The Topology instance with the list of atoms and bonds for this system
    system : :class:`simtk.openmm.System`, optional
        If provided, parameters from this System will be applied to the
        Structure

    Returns
    -------
    struct : :class:`Structure <chemistry.structure.Structure>`
        The structure from the provided topology

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

    return struct
