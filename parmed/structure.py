"""
This module contains the core base class for all of the chemical structures with
various topological and force field features.

Author: Jason Swails
Contributors: Pawel Janowski

Copyright (C) 2014 - 2015  Jason Swails

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
59 Temple Place - Suite 330
Boston, MA 02111-1307, USA.
"""
from __future__ import division, print_function

from parmed.constants import DEG_TO_RAD, SMALL
from parmed.exceptions import ParameterError
from parmed.geometry import (box_lengths_and_angles_to_vectors,
        box_vectors_to_lengths_and_angles)
from parmed.residue import WATER_NAMES
from parmed.topologyobjects import (AtomList, ResidueList, TrackedList,
        DihedralTypeList, Bond, Angle, Dihedral, UreyBradley, Improper, Cmap,
        TrigonalAngle, OutOfPlaneBend, PiTorsion, StretchBend, TorsionTorsion,
        NonbondedException, AcceptorDonor, Group, Atom, ExtraPoint,
        TwoParticleExtraPointFrame, ChiralFrame, MultipoleFrame,
        ThreeParticleExtraPointFrame, OutOfPlaneExtraPointFrame)
from parmed import unit as u
from parmed.utils import tag_molecules
from parmed.utils.decorators import needs_openmm
from parmed.utils.six import string_types, integer_types, iteritems
from parmed.utils.six.moves import zip, range
from parmed.vec3 import Vec3
from copy import copy
import math
import re
import numpy as np

try:
    import pandas as pd
except ImportError:
    pd = None

# Try to import the OpenMM modules
try:
    from simtk.openmm import app
    from simtk import openmm as mm
    from simtk.openmm.app.internal.unitcell import reducePeriodicBoxVectors
except ImportError:
    app = mm = None

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Private attributes and methods

relatere = re.compile(r'RELATED ID: *(\w+) *RELATED DB: *(\w+)', re.I)

def _compare_atoms(old_atom, new_atom, resname, resid, chain):
    """
    Compares two atom instances, along with the residue name, number, and chain
    identifier, to determine if two atoms are actually the *same* atom, but
    simply different conformations

    Parameters
    ----------
    old_atom : :class:`Atom`
        The original atom that has been added to the structure already
    new_atom : :class:`Atom`
        The new atom that we want to see if it is the same as the old atom
    resname : ``str``
        The name of the residue that the new atom would belong to
    resid : ``int``
        The number of the residue that the new atom would belong to
    chain : ``str``
        The chain identifier that the new atom would belong to

    Returns
    -------
    True if they are the same atom, False otherwise
    """
    if old_atom.name != new_atom.name: return False
    if old_atom.residue.name != resname: return False
    if old_atom.residue.number != resid: return False
    if old_atom.residue.chain != chain.strip(): return False
    return True

def _bondi(atom):
    if atom.atomic_number == 6: return 1.7
    if atom.atomic_number == 1: return 1.2
    if atom.atomic_number == 7: return 1.55
    if atom.atomic_number == 14: return 2.1
    if atom.atomic_number == 15: return 1.85
    if atom.atomic_number == 16: return 1.8
    return 1.5

def _mbondi(atom):
    if atom.atomic_number == 1:
        bondeds = atom.bond_partners
        if bondeds[0].atomic_number in (6, 7):
            return 1.3
        if bondeds[0].atomic_number in (8, 16):
            return 0.8
        return 1.2
    return _bondi(atom)

def _mbondi2(atom):
    if atom.atomic_number == 1:
        if atom.bond_partners[0].atomic_number == 7:
            return 1.3
        return 1.2
    return _bondi(atom)

def _mbondi3(atom):
    if atom.residue.name in ('GLU', 'ASP', 'GL4', 'AS4'):
        if atom.name.startswith('OE') or atom.name.startswith('OD'):
            return 1.4
    elif atom.residue.name == 'ARG':
        if atom.name.startswith('HH') or atom.name.startswith('HE'):
            return 1.17
    if atom.name == 'OXT':
        return 1.4
    return _mbondi2(atom)

@needs_openmm
def _gb_rad_screen(atom, model):
    """
    Gets the default GB parameters for a given atom according to a specific
    Generalized Born model

    Parameters
    ----------
    atom : :class:`Atom`
        The atom to get the default GB parameters for
    model : ``app.HCT, app.OBC1, or app.OBC2``
        The GB model to get the default parameters for (app.GBn and app.GBn2 are
        already handled in Structure._get_gb_parameters)

    Returns
    -------
    radius, screen [,alpha, beta, gamma] : ``float, float [,float, float, float]``
        The intrinsic radius of the atom and the screening factor of the atom.
        If the model is GBn2, alpha, beta, and gamma parameters are also
        returned
    """
    if model in (app.OBC1, app.OBC2):
        rad = _mbondi2(atom)
    else:
        rad = _mbondi(atom)
    if atom.atomic_number == 1: return rad, 0.85
    if atom.atomic_number == 6: return rad, 0.72
    if atom.atomic_number == 7: return rad, 0.79
    if atom.atomic_number == 8: return rad, 0.85
    if atom.atomic_number == 9: return rad, 0.88
    if atom.atomic_number == 15: return rad, 0.86
    if atom.atomic_number == 16: return rad, 0.96
    return rad, 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Structure(object):
    """
    A chemical structure composed of atoms, bonds, angles, torsions, and other
    topological features

    Attributes
    ----------
    atoms : :class:`AtomList`
        List of all atoms in the structure
    residues : :class:`ResidueList`
        List of all residues in the structure
    bonds : :class:`TrackedList` (:class:`Bond`)
        List of all bonds in the structure
    angles : :class:`TrackedList` (:class:`Angle`)
        List of all angles in the structure
    dihedrals : :class:`TrackedList` (:class:`Dihedral`)
        List of all dihedrals in the structure -- only one term per dihedral, so
        multi-term dihedral parameters will have the same 4 atoms appear
        multiple times in the list
    urey_bradleys : :class:`TrackedList` (:class:`UreyBradley`)
        List of all Urey-Bradley angle bends in the structure
    impropers : :class:`TrackedList` (:class:`Improper`)
        List of all CHARMM-style improper torsions in the structure
    cmaps : :class:`TrackedList` (:class:`Cmap`)
        List of all CMAP objects in the structure
    trigonal_angles : :class:`TrackedList` (:class:`TrigonalAngle`)
        List of all AMOEBA-style trigonal angles in the structure
    out_of_plane_bends : :class:`TrackedList` (:class:`OutOfPlaneBends`)
        List of all AMOEBA-style out-of-plane bending angles
    pi_torsions : :class:`TrackedList` (:class:`PiTorsion`)
        List of all AMOEBA-style pi-torsion angles
    stretch_bends : :class:`TrackedList` (:class:`StretchBend`)
        List of all AMOEBA-style stretch-bend compound bond/angle terms
    torsion_torsions : :class:`TrackedList` (:class:`TorsionTorsion`)
        List of all AMOEBA-style coupled torsion-torsion terms
    chiral_frames : :class:`TrackedList` (:class:`ChiralFrame`)
        List of all AMOEBA-style chiral frames defined in the structure
    multipole_frames : :class:`TrackedList` (:class:`MultipoleFrame`)
        List of all AMOEBA-style multipole frames defined in the structure
    adjusts : :class:`TrackedList` (:class:`NonbondedException`)
        List of all nonbonded pair-exception rules
    acceptors : :class:`TrackedList` (:class:`AcceptorDonor`)
        List of all H-bond acceptors, if that information is present
    donors : :class:`TrackedList` (:class:`AcceptorDonor`)
        List of all H-bond donors, if that information is present
    groups : :class:`TrackedList` (:class:`Group`)
        List of all CHARMM-style GROUP objects (whatever those are used for)
    box : ``list of 6 floats``
        Box dimensions (a, b, c, alpha, beta, gamma) for the unit cell. If no
        box is defined, `box` is set to `None`
    space_group : ``str``
        The space group of the structure (default is "P 1")
    nrexcl : ``int``
        The number of bonds away that an atom can be in order to be excluded
        from the direct nonbonded summation
    title : ``str``
        Cosmetic only, it is an arbitrary title assigned to the system. Default
        value is an empty string
    positions : u.Quantity(list(Vec3), u.angstroms)
        Unit-bearing atomic coordinates. If not all atoms have coordinates, this
        property is None
    coordinates : np.ndarray of shape (nframes, natom, 3)
        If no coordinates are set, this is set to None. The first frame will
        match the coordinates present on the atoms.

    Notes
    -----
    This class also has a handful of type lists for each of the attributes above
    (excluding `atoms`, `residues`, `chiral_frames`, and `multipole_frames`).
    They are all TrackedList instances that are designed to hold the relevant
    parameter type. The list is:
        bond_types, angle_types, dihedral_types, urey_bradley_types,
        improper_types, cmap_types, trigonal_angle_types,
        out_of_plane_bend_types, pi_torsion_types, stretch_bend_types,
        torsion_torsion_types, adjust_types

    dihedral_types _may_ be a list of :class:`DihedralType` instances, since
    torsion profiles are often represented by a Fourier series with multiple
    terms
    """
    # Force groups assigned to each type of force
    BOND_FORCE_GROUP = 0
    ANGLE_FORCE_GROUP = 1
    DIHEDRAL_FORCE_GROUP = 2
    UREY_BRADLEY_FORCE_GROUP = 3
    IMPROPER_FORCE_GROUP = 4
    CMAP_FORCE_GROUP = 5
    TRIGONAL_ANGLE_FORCE_GROUP = 6
    OUT_OF_PLANE_BEND_FORCE_GROUP = 7
    PI_TORSION_FORCE_GROUP = 8
    STRETCH_BEND_FORCE_GROUP = 9
    TORSION_TORSION_FORCE_GROUP = 10
    NONBONDED_FORCE_GROUP = 11
    RB_TORSION_FORCE_GROUP = 12

    #===================================================

    def __init__(self):

        # Topological object lists
        self.atoms = AtomList()
        self.residues = ResidueList()
        self.bonds = TrackedList()
        self.angles = TrackedList()
        self.dihedrals = TrackedList()
        self.rb_torsions = TrackedList()
        self.urey_bradleys = TrackedList()
        self.impropers = TrackedList()
        self.cmaps = TrackedList()
        self.trigonal_angles = TrackedList()
        self.out_of_plane_bends = TrackedList()
        self.pi_torsions = TrackedList()
        self.stretch_bends = TrackedList()
        self.torsion_torsions = TrackedList()
        self.chiral_frames = TrackedList()
        self.multipole_frames = TrackedList()
        self.adjusts = TrackedList()
        # Extraneous information stored in CHARMM PSF files... not used as far
        # as I can tell for much of anything
        self.acceptors = TrackedList()
        self.donors = TrackedList()
        self.groups = TrackedList()

        # Parameter type lists
        self.bond_types = TrackedList()
        self.angle_types = TrackedList()
        self.dihedral_types = TrackedList()
        self.urey_bradley_types = TrackedList()
        self.improper_types = TrackedList()
        self.rb_torsion_types = TrackedList()
        self.cmap_types = TrackedList()
        self.trigonal_angle_types = TrackedList()
        self.out_of_plane_bend_types = TrackedList()
        self.pi_torsion_types = TrackedList()
        self.stretch_bend_types = TrackedList()
        self.torsion_torsion_types = TrackedList()
        self.adjust_types = TrackedList()

        self._box = None
        self._coordinates = None
        self.space_group = "P 1"
        self.unknown_functional = False
        self.nrexcl = 3
        self.title = ''

    #===================================================

    def __repr__(self):
        natom = len(self.atoms)
        nres = len(self.residues)
        nextra = sum([isinstance(a, ExtraPoint) for a in self.atoms])
        retstr = ['<%s %d atoms' % (type(self).__name__, natom)]
        if nextra > 0:
            retstr.append(' [%d EPs]' % nextra)
        retstr.append('; %d residues' % nres)
        nbond = len(self.bonds)
        retstr.append('; %d bonds' % nbond)
        if self.box is not None:
            retstr.append('; PBC')
            if (abs(self.box[3]-90) > SMALL or abs(self.box[4]-90) > SMALL or
                    abs(self.box[5]-90) > SMALL):
                retstr.append(' (triclinic)')
            else:
                retstr.append(' (orthogonal)')
        # Just assume that if the first bond has a defined type, so does
        # everything else... we don't want __repr__ to be super expensive
        if len(self.bonds) > 0 and self.bonds[0].type is not None:
            retstr.append('; parametrized>')
        else:
            retstr.append('; NOT parametrized>')
        return ''.join(retstr)

    #===================================================

    def add_atom(self, atom, resname, resnum, chain='', inscode=''):
        """
        Adds a new atom to the Structure, adding a new residue to `residues` if
        it has a different name or number as the last residue added and adding
        it to the `atoms` list.

        Parameters
        ----------
        atom : :class:`Atom`
            The atom to add to this residue list
        resname : ``str``
            The name of the residue this atom belongs to
        resnum : ``int``
            The number of the residue this atom belongs to
        chain : ``str``
            The chain ID character for this residue
        inscode : ``str``
            The insertion code ID character for this residue (it is stripped)

        Notes
        -----
        If the residue name and number differ from the last residue in this
        list, a new residue is added and the atom is added to that residue
        """
        self.residues.add_atom(atom, resname, resnum, chain, inscode)
        self.atoms.append(atom)

    #===================================================

    def add_atom_to_residue(self, atom, residue):
        """
        Adds a new atom to the Structure at the end if the given residue

        Parameters
        ----------
        atom : :class:`Atom`
            The atom to add to the system
        residue : :class:`Residue`
            The residue to which to add this atom. It MUST be part of this
            Structure instance already or a ValueError is raised

        Notes
        -----
        This atom is added at the end of the residue and is inserted into the
        `atoms` list in such a way that all residues are composed of atoms
        contiguous in the atoms list. For large systems, this may be a
        relatively expensive operation
        """
        # Make sure residue belongs to this list
        if residue.list is not self.residues:
            raise ValueError('Residue is not part of the structure')
        last_atom = residue.atoms[-1]
        residue.add_atom(atom)
        # Special-case if this is going to be the last atom
        if not self.atoms or last_atom is self.atoms[-1]:
            self.atoms.append(atom)
        else:
            self.atoms.insert(last_atom.idx + 1, atom)

    #===================================================

    def __copy__(self):
        """ A deep copy of the Structure """
        return self.copy(type(self))

    #===================================================

    def copy(self, cls, split_dihedrals=False):
        """
        Makes a copy of the current structure as an instance of a specified
        subclass

        Parameters
        ----------
        cls : Structure subclass
            The returned object is a copy of this structure as a `cls` instance
        split_dihedrals : ``bool``
            If True, then the Dihedral entries will be split up so that each one
            is paired with a single DihedralType (rather than a
            DihedralTypeList)

        Returns
        -------
        *cls* instance
            The instance of the Structure subclass `cls` with a copy of the
            current Structure's topology information
        """
        c = cls()
        for atom in self.atoms:
            res = atom.residue
            a = copy(atom)
            c.add_atom(a, res.name, res.number, res.chain, res.insertion_code)
        # Now copy all of the types
        for bt in self.bond_types:
            c.bond_types.append(copy(bt))
        c.bond_types.claim()
        for at in self.angle_types:
            c.angle_types.append(copy(at))
        c.angle_types.claim()
        if split_dihedrals:
            ndt = 0
            mapdt = {}
            for idt, dt in enumerate(self.dihedral_types):
                if hasattr(dt, '__iter__'):
                    for t in dt:
                        c.dihedral_types.append(copy(t))
                        mapdt.setdefault(idt, []).append(ndt)
                        ndt += 1
                else:
                    mapdt.setdefault(idt, []).append(ndt)
                    ndt += 1
                    c.dihedral_types.append(copy(dt))
        else:
            for dt in self.dihedral_types:
                c.dihedral_types.append(copy(dt))
        c.dihedral_types.claim()
        for ut in self.urey_bradley_types:
            c.urey_bradley_types.append(copy(ut))
        c.urey_bradley_types.claim()
        for it in self.improper_types:
            c.improper_types.append(copy(it))
        c.improper_types.claim()
        for rt in self.rb_torsion_types:
            c.rb_torsion_types.append(copy(rt))
        c.rb_torsion_types.claim()
        for ct in self.cmap_types:
            c.cmap_types.append(copy(ct))
        c.cmap_types.claim()
        for ta in self.trigonal_angle_types:
            c.trigonal_angle_types.append(copy(ta))
        c.trigonal_angle_types.claim()
        for ot in self.out_of_plane_bend_types:
            c.out_of_plane_bend_types.append(copy(ot))
        c.out_of_plane_bend_types.claim()
        for pt in self.pi_torsion_types:
            c.pi_torsion_types.append(copy(pt))
        c.pi_torsion_types.claim()
        for st in self.stretch_bend_types:
            c.stretch_bend_types.append(copy(st))
        c.stretch_bend_types.claim()
        for tt in self.torsion_torsion_types:
            c.torsion_torsion_types.append(copy(tt))
        c.torsion_torsion_types.claim()
        for at in self.adjust_types:
            c.adjust_types.append(copy(at))
        c.adjust_types.claim()
        # Now create the topological objects
        atoms = c.atoms
        for b in self.bonds:
            c.bonds.append(
                    Bond(atoms[b.atom1.idx], atoms[b.atom2.idx],
                         b.type and c.bond_types[b.type.idx])
            )
            c.bonds[-1].funct = b.funct
        for a in self.angles:
            c.angles.append(
                    Angle(atoms[a.atom1.idx], atoms[a.atom2.idx],
                          atoms[a.atom3.idx],
                          a.type and c.angle_types[a.type.idx])
            )
            c.angles[-1].funct = a.funct
        if split_dihedrals:
            for d in self.dihedrals:
                if hasattr(d.type, '__iter__'):
                    for i in range(len(d.type)):
                        ie = d.ignore_end or i < len(d.type) - 1
                        ti = mapdt[d.type.idx][i]
                        c.dihedrals.append(
                                Dihedral(c.atoms[d.atom1.idx],
                                         c.atoms[d.atom2.idx],
                                         c.atoms[d.atom3.idx],
                                         c.atoms[d.atom4.idx],
                                         improper=d.improper, ignore_end=ie,
                                         type=c.dihedral_types[ti])
                        )
                        c.dihedrals[-1]._funct = d._funct
                else:
                    ti = mapdt[d.type.idx][0]
                    c.dihedrals.append(
                        Dihedral(c.atoms[d.atom1.idx], c.atoms[d.atom2.idx],
                                 c.atoms[d.atom3.idx], c.atoms[d.atom4.idx],
                                 improper=d.improper, ignore_end=d.ignore_end,
                                 type=d.type and c.dihedral_types[ti])
                    )
                    c.dihedrals[-1]._funct = d._funct
        else:
            for d in self.dihedrals:
                if d.type is None:
                    typ = None
                else:
                    typ = c.dihedral_types[d.type.idx]
                c.dihedrals.append(
                        Dihedral(atoms[d.atom1.idx], atoms[d.atom2.idx],
                                 atoms[d.atom3.idx], atoms[d.atom4.idx],
                                 improper=d.improper, ignore_end=d.ignore_end,
                                 type=typ)
                )
                c.dihedrals[-1]._funct = d._funct
        for ub in self.urey_bradleys:
            c.urey_bradleys.append(
                    UreyBradley(atoms[ub.atom1.idx], atoms[ub.atom2.idx],
                                ub.type and c.urey_bradley_types[ub.type.idx])
            )
        for i in self.impropers:
            c.impropers.append(
                    Improper(atoms[i.atom1.idx], atoms[i.atom2.idx],
                             atoms[i.atom3.idx], atoms[i.atom4.idx],
                             i.type and c.improper_types[i.type.idx])
            )
        for r in self.rb_torsions:
            c.rb_torsions.append(
                    Dihedral(atoms[r.atom1.idx], atoms[r.atom2.idx],
                             atoms[r.atom3.idx], atoms[r.atom4.idx],
                             type=r.type and c.rb_torsion_types[r.type.idx])
            )
            c.rb_torsions[-1]._funct = r._funct
        for cm in self.cmaps:
            c.cmaps.append(
                    Cmap(atoms[cm.atom1.idx], atoms[cm.atom2.idx],
                         atoms[cm.atom3.idx], atoms[cm.atom4.idx],
                         atoms[cm.atom5.idx],
                         cm.type and c.cmap_types[cm.type.idx])
            )
            c.cmaps[-1].funct = cm.funct
        for t in self.trigonal_angles:
            c.trigonal_angles.append(
                    TrigonalAngle(atoms[t.atom1.idx], atoms[t.atom2.idx],
                                  atoms[t.atom3.idx], atoms[t.atom4.idx],
                                  t.type and c.trigonal_angle_types[t.type.idx])
            )
        for o in self.out_of_plane_bends:
            c.out_of_plane_bends.append(
                    OutOfPlaneBend(atoms[o.atom1.idx], atoms[o.atom2.idx],
                                   atoms[o.atom3.idx], atoms[o.atom4.idx],
                                   o.type and c.out_of_plane_bend_types[o.type.idx])
            )
        for p in self.pi_torsions:
            c.pi_torsions.append(
                    PiTorsion(atoms[p.atom1.idx], atoms[p.atom2.idx],
                              atoms[p.atom3.idx], atoms[p.atom4.idx],
                              atoms[p.atom5.idx], atoms[p.atom6.idx],
                              p.type and c.pi_torsion_types[p.type.idx])
            )
        for s in self.stretch_bends:
            c.stretch_bends.append(
                    StretchBend(atoms[s.atom1.idx], atoms[s.atom2.idx],
                                atoms[s.atom3.idx],
                                s.type and c.stretch_bend_types[s.type.idx])
            )
        for t in self.torsion_torsions:
            c.torsion_torsions.append(
                    TorsionTorsion(atoms[t.atom1.idx], atoms[t.atom2.idx],
                                   atoms[t.atom3.idx], atoms[t.atom4.idx],
                                   atoms[t.atom5.idx],
                                   t.type and c.torsion_torsion_types[t.type.idx])
            )
        for ch in self.chiral_frames:
            c.chiral_frames.append(
                    ChiralFrame(atoms[ch.atom1.idx], atoms[ch.atom2.idx],
                                ch.chirality)
            )
        for m in self.multipole_frames:
            c.multipole_frames.append(
                    MultipoleFrame(atoms[m.atom.idx], m.frame_pt_num, m.vectail,
                                   m.vechead, m.nvec)
            )
        for a in self.adjusts:
            c.adjusts.append(
                    NonbondedException(atoms[a.atom1.idx], atoms[a.atom2.idx],
                                       a.type and c.adjust_types[a.type.idx])
            )
        for a in self.acceptors:
            c.acceptors.append(
                    AcceptorDonor(atoms[a.atom1.idx], atoms[a.atom2.idx])
            )
        for d in self.donors:
            c.donors.append(
                    AcceptorDonor(atoms[d.atom1.idx], atoms[d.atom2.idx])
            )
        for g in self.groups:
            c.groups.append(Group(g.bs, g.type, g.move))
        c.box = copy(self.box)
        c._coordinates = copy(self._coordinates)
        return c

    #===================================================

    def to_dataframe(self):
        """ Creates a pandas.DataFrame object from the current Structure

        Returns
        -------
        dataframe : :class:`pandas.DataFrame`

        Notes
        -----
        The DataFrame will be over all atoms. The columns will be the attributes
        of the atom (as well as its containing residue). Some columns will
        *always* exist. Others will only exist if those attributes have been set
        on the Atom instances (see the :class:`Atom` docs for possible
        attributes and their meaning). The columns that will always be present
        are:

            - number : int
            - name : str
            - type : str
            - atomic_number : int
            - charge : float
            - mass : float
            - nb_idx : int
            - radii : float
            - screen : float
            - occupancy : float
            - bfactor : float
            - altloc : str
            - tree : str
            - join : int
            - irotat : int
            - rmin : float
            - epsilon : float
            - rmin_14 : float
            - epsilon_14 : float
            - resname : str (name of the containing residue)
            - resid : int (Sequential index of the containing residue)
            - resnum : int (original residue number in the input structure)
            - chain : str (chain ID that the containing residue belongs to)

        The following attributes are optionally present if they were present in
        the original file defining the structure:

            - xx : float (x-coordinate position)
            - xy : float (y-coordinate position)
            - xz : float (z-coordinate position)
            - vx : float (x-coordinate velocity)
            - vy : float (y-coordinate velocity)
            - vz : float (z-coordinate velocity)
            - type_idx : int (integer type index for AMOEBA)
            - class_idx : int (integer class type index for AMOEBA)
            - multipole_111 : float (Monopole)
            - multipole_211 : float (1,1 Dipole component)
            - multipole_212 : float (1,2 Dipole component)
            - multipole_222 : float (2,2 Dipole component)
            - multipole_411 : float (1,1 Quadrupole component)
            - multipole_412 : float (1,2 Quadrupole component)
            - multipole_422 : float (2,2 Quadrupole component)
            - multipole_413 : float (1,3 Quadrupole component)
            - multipole_423 : float (2,3 Quadrupole component)
            - multipole_433 : float (3,3 Quadrupole component)
            - polarizability : float (dipole polarizability)
            - vdw_parent : int (index of the vdW parent atom of this atom)
            - segid : segment ID (similar to chain, but for CHARMM)
            - U11 : float (U[1][1] of anisotropic b-factor tensor)
            - U22 : float (U[2][2] of anisotropic b-factor tensor)
            - U33 : float (U[3][3] of anisotropic b-factor tensor)
            - U12 : float (U[1][2] of anisotropic b-factor tensor)
            - U13 : float (U[1][3] of anisotropic b-factor tensor)
            - U23 : float (U[2][3] of anisotropic b-factor tensor)
        """
        if pd is None:
            raise ImportError('pandas is not available; cannot create a pandas '
                              'DataFrame from this Structure')
        ret = pd.DataFrame()

        ret['number'] = [atom.number for atom in self.atoms]
        ret['name'] = [atom.name for atom in self.atoms]
        ret['type'] = [atom.type for atom in self.atoms]
        ret['atomic_number'] = [atom.atomic_number for atom in self.atoms]
        ret['charge'] = [atom.charge for atom in self.atoms]
        ret['mass'] = [atom.mass for atom in self.atoms]
        ret['nb_idx'] = [atom.nb_idx for atom in self.atoms]
        ret['radii'] = [atom.radii for atom in self.atoms]
        ret['screen'] = [atom.screen for atom in self.atoms]
        ret['occupancy'] = [atom.occupancy for atom in self.atoms]
        ret['bfactor'] = [atom.bfactor for atom in self.atoms]
        ret['altloc'] = [atom.altloc for atom in self.atoms]
        ret['tree'] = [atom.tree for atom in self.atoms]
        ret['join'] = [atom.join for atom in self.atoms]
        ret['irotat'] = [atom.irotat for atom in self.atoms]
        ret['rmin'] = [atom.rmin for atom in self.atoms]
        ret['epsilon'] = [atom.epsilon for atom in self.atoms]
        ret['rmin_14'] = [atom.rmin_14 for atom in self.atoms]
        ret['epsilon_14'] = [atom.epsilon_14 for atom in self.atoms]
        ret['resname'] = [atom.residue.name for atom in self.atoms]
        ret['resid'] = [atom.residue.idx for atom in self.atoms]
        ret['resnum'] = [atom.residue.number for atom in self.atoms]
        ret['chain'] = [atom.residue.chain for atom in self.atoms]

        # Now for optional attributes
        # Coordinates
        try:
            coords = pd.DataFrame(
                    [[atom.xx, atom.xy, atom.xz] for atom in self.atoms],
                    columns=['xx', 'xy', 'xz']
            )
        except AttributeError:
            pass
        else:
            ret = ret.join(coords)
        # Velocities
        try:
            vels = pd.DataFrame(
                    [[atom.vx, atom.vy, atom.vz] for atom in self.atoms],
                    columns=['vx', 'vy', 'vz']
            )
        except AttributeError:
            pass
        else:
            ret = ret.join(vels)
        # AMOEBA LJ type
        try:
            ret['type_idx'] = [atom.type_idx for atom in self.atoms]
        except AttributeError:
            pass
        # AMOEBA class type
        try:
            ret['class_idx'] = [atom.class_idx for atom in self.atoms]
        except AttributeError:
            pass
        # AMOEBA multipoles
        try:
            multipoles = pd.DataFrame(
                    [atom.multipoles for atom in self.atoms],
                    columns=['multipole_111', 'multipole_211', 'multipole_212',
                             'multipole_222', 'multipole_411', 'multipole_412',
                             'multipole_422', 'multipole_413', 'multipole_423',
                             'multipole_433']
            )
        except AttributeError:
            pass
        else:
            ret = ret.join(multipoles)
        # AMOEBA polarizabilities
        try:
            ret['polarizability'] = [atom.polarizability for atom in self.atoms]
        except AttributeError:
            pass
        # AMOEBA vdw parent atom
        try:
            ret['vdw_parent'] = [atom.vdw_parent.idx for atom in self.atoms]
        except AttributeError:
            pass
        # SEGID (CHARMM)
        try:
            ret['segid'] = [atom.segid for atom in self.atoms]
        except AttributeError:
            pass
        # anisotropic b-factors
        none6 = [None] * 6
        anisos = [atom.anisou for atom in self.atoms]
        for i, aniso in enumerate(anisos):
            if hasattr(aniso, 'tolist'):
                anisos[i] = aniso.tolist()
        all_nones = True
        for i, aniso in enumerate(anisos):
            if aniso is None:
                anisos[i] = none6
            elif all_nones:
                all_nones = False
        if not all_nones:
            ret = ret.join(
                    pd.DataFrame(anisos,
                        columns=['U11', 'U22', 'U33', 'U12', 'U13', 'U23'])
            )
        return ret

    #===================================================

    def is_changed(self):
        """ Determines if any of the topology has changed for this structure """
        return (self.atoms.changed or self.residues.changed or
                self.bonds.changed or self.trigonal_angles.changed or
                self.dihedrals.changed or self.urey_bradleys.changed or
                self.impropers.changed or self.cmaps.changed or
                self.angles.changed or self.out_of_plane_bends.changed or
                self.pi_torsions.changed or self.stretch_bends.changed or
                self.torsion_torsions.changed or self.chiral_frames.changed or
                self.multipole_frames.changed or self.adjusts.changed or
                self.acceptors.changed or self.donors.changed or
                self.groups.changed or self.bond_types.changed or
                self.angle_types.changed or self.dihedral_types.changed or
                self.urey_bradley_types.changed or self.cmap_types.changed or
                self.improper_types.changed or self.adjust_types.changed or
                self.trigonal_angle_types.changed or
                self.out_of_plane_bends.changed or
                self.stretch_bend_types.changed or
                self.torsion_torsion_types.changed or
                self.pi_torsion_types.changed or self.rb_torsions.changed or
                self.rb_torsion_types.changed)

    #===================================================

    def unchange(self):
        """ Toggles all lists so that they do not indicate any changes """
        self.atoms.changed = False
        self.residues.changed = False
        self.bonds.changed = False
        self.angles.changed = False
        self.dihedrals.changed = False
        self.urey_bradleys.changed = False
        self.impropers.changed = False
        self.cmaps.changed = False
        self.trigonal_angles.changed = False
        self.out_of_plane_bends.changed = False
        self.pi_torsions.changed = False
        self.stretch_bends.changed = False
        self.torsion_torsions.changed = False
        self.chiral_frames.changed = False
        self.multipole_frames.changed = False
        self.adjusts.changed = False
        self.acceptors.changed = False
        self.donors.changed = False
        self.groups.changed = False

        # Parameter type lists
        self.bond_types.changed = False
        self.angle_types.changed = False
        self.dihedral_types.changed = False
        self.urey_bradley_types.changed = False
        self.improper_types.changed = False
        self.cmap_types.changed = False
        self.trigonal_angle_types.changed = False
        self.out_of_plane_bend_types.changed = False
        self.pi_torsion_types.changed = False
        self.stretch_bend_types.changed = False
        self.torsion_torsion_types.changed = False
        self.adjust_types.changed = False

    #===================================================

    def prune_empty_terms(self):
        """
        Looks through all of the topological lists and gets rid of terms
        in which at least one of the atoms is None or has an `idx` attribute set
        to -1 (indicating that it has been removed from the `atoms` atom list)
        """
        self._prune_empty_bonds()
        self._prune_empty_angles()
        self._prune_empty_dihedrals()
        self._prune_empty_rb_torsions()
        self._prune_empty_ureys()
        self._prune_empty_impropers()
        self._prune_empty_cmaps()
        self._prune_empty_trigonal_angles()
        self._prune_empty_out_of_plane_bends()
        self._prune_empty_pi_torsions()
        self._prune_empty_stretch_bends()
        self._prune_empty_torsion_torsions()
        self._prune_empty_chiral_frames()
        self._prune_empty_multipole_frames()
        self._prune_empty_adjusts()

    #===================================================

    def update_dihedral_exclusions(self):
        """
        Nonbonded exclusions and exceptions have the following priority:

        bond -> angle -> dihedral

        Since bonds and angles are completely excluded, any ring systems in
        which two atoms are attached by a bond or angle as well as a dihedral
        should be completely excluded as the bond and angle exclusion rules take
        precedence.  If a Bond or Angle was _added_ to the structure between a
        pair of atoms previously connected only by a dihedral term, it's
        possible that those two atoms have both an exclusion *and* an exception
        defined. The result of this scenario is that sander and pmemd will
        happily compute an energy, _including_ the 1-4 nonbonded terms between
        atoms now connected by a bond or an Angle.  OpenMM, on the other hand,
        will complain about an exception specified multiple times. This method
        scans through all of the dihedrals in which `ignore_end` is `False` and
        turns it to `True` if the two end atoms are in the bond or angle
        partners arrays
        """
        set14 = set()
        for dihedral in self.dihedrals:
            if dihedral.ignore_end : continue
            if (dihedral.atom1 in dihedral.atom4.bond_partners or
                dihedral.atom1 in dihedral.atom4.angle_partners):
                dihedral.ignore_end = True
            elif (dihedral.atom1.idx, dihedral.atom4.idx) in set14:
                # Avoid double counting of 1-4 in a six-membered ring
                dihedral.ignore_end = True
            else:
                set14.add((dihedral.atom1.idx, dihedral.atom4.idx))
                set14.add((dihedral.atom4.idx, dihedral.atom1.idx))

    #===================================================

    def strip(self, selection):
        """
        Deletes a subset of the atoms corresponding to an atom-based selection.

        Parameters
        ----------
        selection : :class:`AmberMask`, ``str``, or ``iterable``
            This is the selection of atoms that will be deleted from this
            structure. If it is a string, it will be interpreted as an
            AmberMask. If it is an AmberMask, it will be converted to a
            selection of atoms. If it is an iterable, it must be the same length
            as the `atoms` list.
        """
        from parmed.amber import AmberMask
        if isinstance(selection, AmberMask):
            if selection.parm is not self:
                raise TypeError('passed mask does not belong to Structure')
            sel = selection.Selection()
        elif isinstance(selection, string_types):
            sel = AmberMask(self, selection).Selection()
        else:
            try:
                sel = list(selection)
            except TypeError:
                raise TypeError('Selection not a supported type [%s]' %
                                type(selection))
            if len(sel) != len(self.atoms):
                raise ValueError('Selection iterable wrong length')
        atomlist = sorted([i for i, s in enumerate(sel) if s])
        for i in reversed(atomlist):
            del self.atoms[i]
        self.prune_empty_terms()
        self.residues.prune()
        self.unchange()
        # Slice out coordinates if present
        if self._coordinates is not None:
            self._coordinates = self._coordinates[np.array(sel)==0]

    #===================================================

    def __getitem__(self, selection):
        """
        Allows extracting a single atom from the structure or a slice of atoms
        as a new Structure instance. The following syntaxes are allowed:

            - struct[str] : str is interpreted as an Amber selection mask
            - struct[[sel,[sel,]],sel] : sel can be a list of indices (or
                            strings for chain selection) an integer, or a slice

        Parameters
        ----------
        selection : str, or slice|iter|str, slice|iter|int, slice|iter|int
            Which atoms to select for the new structure

        Returns
        -------
        struct : :class:`Structure`
            If more than one atom is selected, the resulting return value is a
            new structure containing all selected atoms and all valence terms
            (and types) in which all involved atoms are present

        or

        atom : :class:`Atom`
            If only a single atom is selected

        or

        None
            If no atoms are selected

        Notes
        -----
        When selecting more than 1 atom, this is a costly operation. It is
        currently implemented by making a full copy of the object and then
        stripping the unused atoms.

        Raises
        ------
        ``TypeError`` : if selection (or any part of the selection) is not an
                        allowed type

        ``ValueError`` : if the selection is a boolean-like list and its length
                         is not the same as the number of atoms in the system
        """
        from parmed.amber import AmberMask
        if isinstance(selection, integer_types):
            return self.atoms[selection]

        # Now we select a subset of atoms. Convert "selection" into a natom list
        # with 0s and 1s, depending on what the input selection is
        if isinstance(selection, string_types):
            mask = AmberMask(self, selection)
            selection = mask.Selection()
        elif isinstance(selection, slice):
            sel = [0 for a in self.atoms]
            for idx in list(range(len(self.atoms)))[selection]:
                sel[idx] = 1
            selection = sel
        elif isinstance(selection, tuple) and len(selection) in (2, 3):
            # This is a selection of chains, and/or residues, and/or atoms
            if len(selection) == 2:
                # Select just residues and atoms
                ressel, atomsel = selection
                has_chain = False
            elif len(selection) == 3:
                chainsel, ressel, atomsel = selection
                chainmap = dict() # First residue in each chain
                for r in self.residues:
                    if r.chain not in chainmap:
                        chainmap[r.chain] = r.idx
                if isinstance(chainsel, string_types):
                    if chainsel in chainmap:
                        chainset = set([chainsel])
                    else:
                        # The selected chain is not in the system. Cannot
                        # possibly select *any* atoms. Bail now for speed
                        return None
                elif isinstance(chainsel, slice):
                    # Build an ordered set of chains
                    chains = [self.residues[0].chain]
                    for res in self.residues:
                        if res.chain != chains[-1]:
                            chains.append(res.chain)
                    chainset = set(chains[chainsel])
                else:
                    chainset = set(chainsel)
                for chain in chainset:
                    if chain in chainmap:
                        break
                else:
                    # No requested chain is present. Bail now for speed
                    return None
                has_chain = True
            # Residue selection can either be by name or index
            if isinstance(ressel, slice):
                resset = set(list(range(len(self.residues)))[ressel])
            elif isinstance(ressel, string_types) or isinstance(ressel,
                    integer_types):
                resset = set([ressel])
            else:
                resset = set(ressel)
            if isinstance(atomsel, slice):
                atomset = set(list(range(len(self.atoms)))[atomsel])
            elif isinstance(atomsel, string_types) or isinstance(atomsel,
                    integer_types):
                atomset = set([atomsel])
            else:
                atomset = set(atomsel)
            if has_chain:
                selection = [
                        (a.residue.chain in chainset) and
                        (a.residue.name in resset or
                         a.residue.idx-chainmap[a.residue.chain] in resset) and
                        (a.name in atomset or a.idx-a.residue[0].idx in atomset)
                            for a in self.atoms
                ]
            else:
                selection = [
                    (a.name in atomset or a.idx-a.residue[0].idx in atomset)
                    and (a.residue.name in resset or a.residue.idx in resset)
                            for a in self.atoms
                ]
        else:
            # Assume it is an iterable. If it is the same length as the atoms,
            # it is a boolean mask array. Otherwise, it is a list of atom
            # indices to select
            sel = [0 for atom in self.atoms]
            selection = list(selection)
            if len(selection) == len(self.atoms):
                for i, val in enumerate(selection):
                    if val: sel[i] = 1
            elif len(selection) > len(self.atoms):
                raise ValueError('Selection iterable is too long')
            else:
                try:
                    for val in selection:
                        sel[val] = 1
                except IndexError:
                    raise ValueError('Selected atom out of range')
            selection = sel
        # Make sure we have an integer array of 0s and 1s
        selection = [int(bool(x)) for x in selection]
        sumsel = sum(selection)
        if sumsel == 0:
            # No atoms selected. Return None
            return None
        if sumsel == 1:
            # 1 atom selected; return that atom
            return self.atoms[selection.index(1)]
        # The cumulative sum of selection will give our index + 1 of each
        # selected atom into the new structure
        scan = [selection[0]]
        for i in range(1, len(selection)):
            scan.append(scan[i-1] + selection[i])
        # Zero-out the unselected atoms
        scan = [x * y for x, y in zip(scan, selection)]
        # Copy all parameters
        struct = type(self)()
        for i, atom in enumerate(self.atoms):
            if not selection[i]: continue
            res = atom.residue
            if res.number == 0:
                num = res.idx
            else:
                num = res.number
            struct.add_atom(copy(atom), res.name, num, res.chain,
                            res.insertion_code)
        def copy_valence_terms(oval, otyp, sval, styp, attrlist):
            """ Copies the valence terms from one list to another;
            oval=Other VALence; otyp=Other TYPe; sval=Self VALence;
            styp=Self TYPe; attrlist=ATTRibute LIST (atom1, atom2, ...)
            """
            otypcp = [copy(typ) for typ in styp]
            used_types = [False for typ in otypcp]
            for val in sval:
                ats = [getattr(val, attr) for attr in attrlist]
                # Make sure all of our atoms in this valence term is "selected"
                indices = [scan[at.idx] for at in ats if isinstance(at, Atom)]
                if not all(indices):
                    continue
                # Add the type if applicable
                kws = dict()
                if otypcp and val.type is not None:
                    kws['type'] = otypcp[val.type.idx]
                    used_types[val.type.idx] = True
                for i, at in enumerate(ats):
                    if isinstance(at, Atom):
                        ats[i] = struct.atoms[scan[at.idx]-1]
                oval.append(type(val)(*ats, **kws))
                if hasattr(val, 'funct'):
                    oval[-1].funct = val.funct
            # Now tack on the "new" types copied from `other`
            for used, typ in zip(used_types, otypcp):
                if used: otyp.append(typ)
            if hasattr(otyp, 'claim'):
                otyp.claim()
        copy_valence_terms(struct.bonds, struct.bond_types, self.bonds,
                           self.bond_types, ['atom1', 'atom2'])
        copy_valence_terms(struct.angles, struct.angle_types, self.angles,
                           self.angle_types, ['atom1', 'atom2', 'atom3'])
        copy_valence_terms(struct.dihedrals, struct.dihedral_types,
                           self.dihedrals, self.dihedral_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'improper',
                           'ignore_end'])
        copy_valence_terms(struct.rb_torsions, struct.rb_torsion_types,
                           self.rb_torsions, self.rb_torsion_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'improper',
                           'ignore_end'])
        copy_valence_terms(struct.urey_bradleys, struct.urey_bradley_types,
                           self.urey_bradleys, self.urey_bradley_types,
                           ['atom1', 'atom2'])
        copy_valence_terms(struct.impropers, struct.improper_types,
                           self.impropers, self.improper_types,
                           ['atom1', 'atom2', 'atom3', 'atom4'])
        copy_valence_terms(struct.cmaps, struct.cmap_types,
                           self.cmaps, self.cmap_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
        copy_valence_terms(struct.trigonal_angles, struct.trigonal_angle_types,
                           self.trigonal_angles, self.trigonal_angle_types,
                           ['atom1', 'atom2', 'atom3', 'atom4'])
        copy_valence_terms(struct.out_of_plane_bends,
                           struct.out_of_plane_bend_types,
                           self.out_of_plane_bends,
                           self.out_of_plane_bend_types,
                           ['atom1', 'atom2', 'atom3', 'atom4'])
        copy_valence_terms(struct.pi_torsions, struct.pi_torsion_types,
                           self.pi_torsions, self.pi_torsion_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'atom5',
                            'atom6'])
        copy_valence_terms(struct.stretch_bends, struct.stretch_bend_types,
                           self.stretch_bends, self.stretch_bend_types,
                           ['atom1', 'atom2', 'atom3'])
        copy_valence_terms(struct.torsion_torsions, struct.torsion_torsion_types,
                           self.torsion_torsions, self.torsion_torsion_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
        copy_valence_terms(struct.chiral_frames, [], self.chiral_frames, [],
                           ['atom1', 'atom2', 'chirality'])
        copy_valence_terms(struct.multipole_frames, [], self.multipole_frames,
                           [], ['atom', 'frame_pt_num', 'vectail', 'vechead',
                           'nvec'])
        copy_valence_terms(struct.adjusts, struct.adjust_types, self.adjusts,
                           self.adjust_types, ['atom1', 'atom2'])
        copy_valence_terms(struct.donors, [], self.donors, [],
                           ['atom1', 'atom2'])
        copy_valence_terms(struct.acceptors, [], self.acceptors, [],
                           ['atom1', 'atom2'])
        copy_valence_terms(struct.groups, [], self.groups, [],
                           ['bs', 'type', 'move'])
        return struct

    #===================================================

    def split(self):
        """
        Split the current Structure into separate Structure instances for each
        unique molecule. A molecule is defined as all atoms connected by a graph
        of covalent bonds.

        Returns
        -------
        [structs, counts] : list of (:class:`Structure`, int) tuples
            List of all molecules in the order that they appear first in the
            parent structure accompanied by the number of times that molecule
            appears in the Structure
        """
        tag_molecules(self)
        mollist = [atom.marked for atom in self.atoms]
        nmol = max(mollist)
        structs = []
        counts = []
        res_molecules = dict()
        # Divvy up the atoms into their separate molecules
        molatoms = [[] for i in range(nmol)]
        for atom in self.atoms:
            molatoms[atom.marked-1].append(atom)
        for i in range(nmol):
            sel = molatoms[i]
            involved_residues = set(atom.residue.idx for atom in sel)
            # Shortcut -- keep names of single-residue molecules and the number
            # of atoms present in those residues in a dict and use that to see
            # if two molecules are unique -- this gives drastic speedup for
            # systems with many duplicated single-residue molecules
            # (i.e., just about every solvent out there)
            if len(involved_residues) == 1:
                res = sel[0].residue
                names = tuple(a.name for a in res)
                if (res.name, len(res), names) in res_molecules:
                    counts[res_molecules[(res.name, len(res), names)]] += 1
                    continue
                else:
                    res_molecules[(res.name, len(res), names)] = len(structs)
            is_duplicate = False
            for j, struct in enumerate(structs):
                if len(struct.atoms) == len(sel):
                    for a1, a2 in zip(struct.atoms, sel):
                        if a1.residue is None:
                            if a2.residue is not None:
                                break
                        elif a2.residue is None:
                            break
                        elif a1.residue.name != a2.residue.name:
                            break
                        if not a1.type and not a2.type:
                            if a1.name != a2.name: break
                        else:
                            if a1.type != a2.type: break
                    else:
                        counts[j] += 1
                        is_duplicate = True
                        break
            if not is_duplicate:
                mol = self[[atom.marked == i+1 for atom in self.atoms]]
                if isinstance(mol, Atom):
                    s = type(self)()
                    s.add_atom(copy(mol), mol.residue.name,
                               mol.residue.number, mol.residue.chain,
                               mol.residue.insertion_code)
                    mol = s
                structs.append(mol)
                counts.append(1)
        return list(zip(structs, counts))

    #===================================================

    @property
    @needs_openmm
    def topology(self):
        """
        The OpenMM Topology object. Cached when possible, but any changes to the
        Structure instance results in the topology being deleted and rebuilt

        Notes
        -----
        This function calls ``prune_empty_terms`` if any topology lists have
        changed.
        """
        if not self.is_changed():
            try:
                return self._topology
            except AttributeError:
                pass
        else:
            self.prune_empty_terms()
            self.unchange()

        self._topology = top = app.Topology()
        chain = top.addChain()
        try:
            last_chain = self.residues[0].chain
            last_residue = None
            last_omm_residue = None
        except IndexError:
            # Empty topology
            return self._topology
        # Add the atoms
        for i, atom in enumerate(self.atoms):
            # See if we need to add a new residue
            if atom.residue is not last_residue:
                # See if we need a new chain
                if last_chain != atom.residue.chain:
                    last_chain = atom.residue.chain
                    chain = top.addChain()
                last_residue = atom.residue
                last_omm_residue = top.addResidue(atom.residue.name, chain)
            try:
                elem = app.element.Element.getByAtomicNumber(atom.atomic_number)
            except KeyError:
                elem = None
            top.addAtom(atom.name, elem, last_omm_residue)
        # Add the bonds
        atoms = list(top.atoms())
        for bond in self.bonds:
            top.addBond(atoms[bond.atom1.idx], atoms[bond.atom2.idx])
        # Set the unit cell dimensions
        if self.box is not None:
            top.setPeriodicBoxVectors(
                    reducePeriodicBoxVectors(
                        box_lengths_and_angles_to_vectors(*self.box)
                    )
            )
        return top

    #===================================================

    @property
    def positions(self):
        """
        A list of 3-element Quantity tuples of dimension length representing the
        atomic positions for every atom in the system. If set with unitless
        numbers, those numbers are assumed to be in angstroms. If any atoms do
        not have coordinates, this is simply ``None``.
        """
        try:
            return [Vec3(a.xx,a.xy,a.xz) for a in self.atoms] * u.angstroms
        except AttributeError:
            return None

    @positions.setter
    def positions(self, value):
        """
        A list of 3-element Quantity tuples of dimension length representing the
        atomic positions for every atom in the system. If set with unitless
        numbers, those numbers are assumed to be in angstroms.

        Raises
        ------
        ValueError if the positions are not either a (natom, 3)-shape iterable
        of iterables or a (natom*3)-length iterable.
        """
        if u.is_quantity(value):
            value = value.value_in_unit(u.angstroms)
        # See if the array is flattened
        if len(value) == len(self.atoms):
            # It had better all be 3-length iterables
            for i, atom in enumerate(self.atoms):
                atom.xx, atom.xy, atom.xz = value[i]
        elif len(value) == 3 * len(self.atoms):
            for i, atom in enumerate(self.atoms):
                i3 = i * 3
                atom.xx, atom.xy, atom.xz = value[i3:i3+3]
        else:
            raise ValueError('Wrong shape for position array')

    @property
    def coordinates(self):
        if self._coordinates is None:
            try:
                coords = [[a.xx, a.xy, a.xz] for a in self.atoms]
            except AttributeError:
                return None
            else:
                return np.array(coords)
        elif self.is_changed():
            # Make sure our first frame matches our atomic coordinates. If not,
            # delete those coordinates
            coords = np.array([[a.xx, a.xy, a.xz] for a in self.atoms])
            if np.abs(self._coordinates[0] - coords).max() > SMALL:
                self._coordinates = None
                return coords
        assert len(self._coordinates.shape) == 3, \
                'Internal coordinate shape wrong'
        assert self._coordinates.shape[1] == len(self.atoms), \
                'Coordinate shape different from number of atoms'
        return self._coordinates[0]

    @coordinates.setter
    def coordinates(self, value):
        """ Setting coordinates will also set xx, xy, and xz on the atoms """
        if value is None:
            # Wipe out coordinates
            self._coordinates = None
            for atom in self.atoms:
                try:
                    del atom.xx, atom.xy, atom.xz
                except AttributeError:
                    pass
        else:
            if u.is_quantity(value):
                value = value.value_in_unit(u.angstroms)
            value = list(value)
            coords = np.array(value, dtype=np.float64, copy=False, subok=True)
            coords = coords.reshape((-1, len(self.atoms), 3))
            for a, xyz in zip(self.atoms, coords[0]):
                a.xx, a.xy, a.xz = xyz
            self._coordinates = coords

    def get_coordinates(self, frame='all'):
        """
        In some cases, multiple conformations may be stored in the Structure.
        This function retrieves a particular frame's coordinates

        Parameters
        ----------
        frame : int or 'all', optional
            The frame number whose coordinates should be retrieved. Default is
            'all'

        Returns
        -------
        coords : np.ndarray, shape([#,] natom, 3) or None
            If frame is 'all', all coordinates are returned with shape
            (#, natom, 3). Otherwise the requested frame is returned with shape
            (natom, 3). If no coordinates exist and 'all' is requested, None is
            returned

        Raises
        ------
        IndexError if there are fewer than ``frame`` coordinates
        """
        if self.is_changed() and self._coordinates is not None:
            try:
                coords = np.array([[a.xx, a.xy, a.xz] for a in self.atoms])
            except AttributeError:
                self._coordinates = None
            else:
                if np.abs(self._coordinates[0] - coords).max() > SMALL:
                    self._coordinates = None
        if frame == 'all':
            if self._coordinates is not None:
                return self._coordinates
            try:
                return np.array([[a.xx, a.xy, a.xz]
                    for a in self.atoms]).reshape((1, len(self.atoms), 3))
            except AttributeError:
                return None
        elif self._coordinates is None:
            if frame == 0:
                try:
                    return np.array([[a.xx, a.xy, a.xz] for a in self.atoms])
                except AttributeError:
                    raise IndexError('No coordinate frames present')
            # We requested *not* the first frame
            raise IndexError('No coordinate frames present')
        return self._coordinates[frame]

    @property
    def box(self):
        return self._box

    @box.setter
    def box(self, value):
        if value is None:
            self._box = None
        else:
            if isinstance(value, np.ndarray):
                box = value
            else:
                box = list(value)
                if len(box) != 6:
                    raise ValueError('Box information must be 6 floats')
                if u.is_quantity(box[0]):
                    box[0] = box[0].value_in_unit(u.angstroms)
                if u.is_quantity(box[1]):
                    box[1] = box[1].value_in_unit(u.angstroms)
                if u.is_quantity(box[2]):
                    box[2] = box[2].value_in_unit(u.angstroms)
                if u.is_quantity(box[3]):
                    box[3] = box[3].value_in_unit(u.degrees)
                if u.is_quantity(box[4]):
                    box[4] = box[4].value_in_unit(u.degrees)
                if u.is_quantity(box[5]):
                    box[5] = box[5].value_in_unit(u.degrees)
            box = np.array(box, dtype=np.float64, copy=False, subok=True)
            if box.shape != (6,):
                raise ValueError('Box information must be 6 floats')
            self._box = box

    #===================================================

    @property
    def velocities(self):
        """
        A (natom, 3)-shape numpy array with atomic velocities for every atom in
        the system (in units of angstrom/picosecond), or None if there are no
        velocities
        """
        try:
            return np.array([[a.vx, a.vy, a.vz] for a in self.atoms])
        except AttributeError:
            return None

    @velocities.setter
    def velocities(self, value):
        """
        A list of 3-element Quantity tuples of dimension length representing the
        atomic velocities for every atom in the system
        """
        if u.is_quantity(value):
            value = value.value_in_unit(u.angstroms/u.picoseconds)
        if value is None:
            for atom in self.atoms:
                try:
                    del atom.vx, atom.vy, atom.vz
                except AttributeError:
                    pass
        else:
            value = np.array(value, copy=False).reshape((-1,len(self.atoms),3))
            for atom, xyz in zip(self.atoms, value[0]):
                atom.vx, atom.vy, atom.vz = xyz

    #===================================================

    def has_NBFIX(self):
        """
        Returns whether or not any pairs of atom types have their LJ
        interactions modified by an NBFIX definition

        Returns
        -------
        has_nbfix : bool
            If True, at least two atom types have NBFIX mod definitions
        """
        typemap = dict()
        for a in self.atoms:
            if a.atom_type is None: continue
            typemap[str(a.atom_type)] = a.atom_type
        # Now we have a map of all atom types that we have defined in our
        # system. Look through all of the atom types and see if any of their
        # NBFIX definitions are also keys in typemap
        for key, type in iteritems(typemap):
            for key in type.nbfix:
                if key in typemap:
                    return True
        return False

    #===================================================

    @property
    def box_vectors(self):
        """
        3, 3-element tuple of unit cell vectors that are Quantity objects of
        dimension length
        """
        if self._box is None: return None
        return box_lengths_and_angles_to_vectors(*self.box)

    @box_vectors.setter
    def box_vectors(self, value):
        """
        3, 3-element tuple of unit cell vectors that are Quantity objects of
        dimension length
        """
        (a, b, c), (A, B, G) = box_vectors_to_lengths_and_angles(*value)
        a = a.value_in_unit(u.angstroms)
        b = b.value_in_unit(u.angstroms)
        c = c.value_in_unit(u.angstroms)
        A = A.value_in_unit(u.degrees)
        B = B.value_in_unit(u.degrees)
        G = G.value_in_unit(u.degrees)
        self._box = np.array([a, b, c, A, B, G], dtype=np.float64)

    #===================================================

    @needs_openmm
    def createSystem(self, nonbondedMethod=None,
                     nonbondedCutoff=8.0*u.angstroms,
                     switchDistance=0.0*u.angstroms,
                     constraints=None,
                     rigidWater=True,
                     implicitSolvent=None,
                     implicitSolventKappa=None,
                     implicitSolventSaltConc=0.0*u.moles/u.liters,
                     temperature=298.15*u.kelvin,
                     soluteDielectric=1.0,
                     solventDielectric=78.5,
                     useSASA=False,
                     removeCMMotion=True,
                     hydrogenMass=None,
                     ewaldErrorTolerance=0.0005,
                     flexibleConstraints=True,
                     verbose=False,
                     forceNBFIX=False):
        """
        Construct an OpenMM System representing the topology described by the
        prmtop file.

        Parameters
        ----------
        nonbondedMethod : cutoff method
            This is the cutoff method. It can be either the NoCutoff,
            CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald objects from the
            simtk.openmm.app namespace
        nonbondedCutoff : float or distance Quantity
            The nonbonded cutoff must be either a floating point number
            (interpreted as nanometers) or a Quantity with attached units. This
            is ignored if nonbondedMethod is NoCutoff.
        switchDistance : float or distance Quantity
            The distance at which the switching function is turned on for van
            der Waals interactions. This is ignored when no cutoff is used, and
            no switch is used if switchDistance is 0, negative, or greater than
            the cutoff
        constraints : None, app.HBonds, app.HAngles, or app.AllBonds
            Which type of constraints to add to the system (e.g., SHAKE). None
            means no bonds are constrained. HBonds means bonds with hydrogen are
            constrained
        rigidWater : bool=True
            If True, water is kept rigid regardless of the value of constraints.
            A value of False is ignored if constraints is not None.
        implicitSolvent : None, app.HCT, app.OBC1, app.OBC2, app.GBn, app.GBn2
            The Generalized Born implicit solvent model to use.
        implicitSolventKappa : float or 1/distance Quantity = None
            This is the Debye kappa property related to modeling saltwater
            conditions in GB. It should have units of 1/distance (1/nanometers
            is assumed if no units present). A value of None means that kappa
            will be calculated from implicitSolventSaltConc (below)
        implicitSolventSaltConc : float or amount/volume Quantity=0 moles/liter
            If implicitSolventKappa is None, the kappa will be computed from the
            salt concentration. It should have units compatible with mol/L
        temperature : float or temperature Quantity = 298.15 kelvin
            This is only used to compute kappa from implicitSolventSaltConc
        soluteDielectric : float=1.0
            The dielectric constant of the protein interior used in GB
        solventDielectric : float=78.5
            The dielectric constant of the water used in GB
        useSASA : bool=False
            If True, use the ACE non-polar solvation model. Otherwise, use no
            SASA-based nonpolar solvation model.
        removeCMMotion : bool=True
            If True, the center-of-mass motion will be removed periodically
            during the simulation. If False, it will not.
        hydrogenMass : float or mass quantity = None
            If not None, hydrogen masses will be changed to this mass and the
            difference subtracted from the attached heavy atom (hydrogen mass
            repartitioning)
        ewaldErrorTolerance : float=0.0005
            When using PME or Ewald, the Ewald parameters will be calculated
            from this value
        flexibleConstraints : bool=True
            If False, the energies and forces from the constrained degrees of
            freedom will NOT be computed. If True, they will (but those degrees
            of freedom will *still* be constrained).
        verbose : bool=False
            If True, the progress of this subroutine will be printed to stdout

        Notes
        -----
        This function calls prune_empty_terms if any Topology lists have changed
        """
        if self.unknown_functional:
            raise ParameterError('Cannot createSystem: unknown functional')
        # Establish defaults
        if nonbondedMethod is None:
            nonbondedMethod = app.NoCutoff
        system = mm.System()
        # Make sure periodic simulations have a box
        if nonbondedMethod in (app.CutoffPeriodic, app.PME, app.Ewald):
            if self.box is None:
                raise ValueError('No periodic boundary conditions detected')
        # Do hydrogen mass repartitioning if necessary
        masses = [atom.mass for atom in self.atoms]
        if hydrogenMass is not None:
            if u.is_quantity(hydrogenMass):
                hydrogenMass = hydrogenMass.value_in_unit(u.dalton)
            if hydrogenMass <= 0:
                raise ValueError('Hydrogen mass must be positive')
            for atom in self.atoms:
                if atom.element == 1:
                    heavy_atom = None
                    for a2 in atom.bond_partners:
                        if a2.element != 1:
                            heavy_atom = a2
                            break
                    if heavy_atom is not None:
                        masses[atom.idx] = hydrogenMass
                        masses[heavy_atom.idx] -= hydrogenMass - atom.mass
        for mass in masses: system.addParticle(mass)
        self.omm_add_constraints(system, constraints, rigidWater)
        # Prune empty terms if we have changed
        if self.is_changed():
            self.prune_empty_terms()
            self.unchange()
        # Add the various types of forces
        if verbose: print('Adding bonds...')
        self._add_force_to_system(system,
                self.omm_bond_force(constraints, rigidWater,
                                    flexibleConstraints)
        )
        if verbose: print('Adding angles...')
        self._add_force_to_system(system,
                self.omm_angle_force(constraints, flexibleConstraints)
        )
        if verbose: print('Adding dihedrals...')
        self._add_force_to_system(system, self.omm_dihedral_force())
        if verbose: print('Adding Ryckaert-Bellemans torsions...')
        self._add_force_to_system(system, self.omm_rb_torsion_force())
        if verbose: print('Adding Urey-Bradleys...')
        self._add_force_to_system(system, self.omm_urey_bradley_force())
        if verbose: print('Adding improper torsions...')
        self._add_force_to_system(system, self.omm_improper_force())
        if verbose: print('Adding CMAP torsions...')
        self._add_force_to_system(system, self.omm_cmap_force())
        if verbose: print('Adding trigonal angle terms...')
        self._add_force_to_system(system, self.omm_trigonal_angle_force())
        if verbose: print('Adding out-of-plane bends...')
        self._add_force_to_system(system, self.omm_out_of_plane_bend_force())
        if verbose: print('Adding pi-torsions...')
        self._add_force_to_system(system, self.omm_pi_torsion_force())
        if verbose: print('Adding stretch-bends...')
        self._add_force_to_system(system, self.omm_stretch_bend_force())
        if verbose: print('Adding torsion-torsions...')
        self._add_force_to_system(system, self.omm_torsion_torsion_force())
        if verbose: print('Adding Nonbonded force...')
        if implicitSolvent is not None:
            rf_dielc = 1.0
        else:
            rf_dielc = 78.5
        self._add_force_to_system(system,
                self.omm_nonbonded_force(nonbondedMethod, nonbondedCutoff,
                                         switchDistance,ewaldErrorTolerance,
                                         rf_dielc)
        )
        if implicitSolvent is not None:
            if verbose: print('Adding GB force...')
            self._add_force_to_system(system,
                    self.omm_gbsa_force(implicitSolvent, nonbondedMethod,
                                        nonbondedCutoff, soluteDielectric,
                                        solventDielectric, implicitSolventKappa,
                                        implicitSolventSaltConc, temperature,
                                        useSASA)
            )
        if removeCMMotion:
            system.addForce(mm.CMMotionRemover())
        if self.box is not None:
            system.setDefaultPeriodicBoxVectors(
                    *reducePeriodicBoxVectors(self.box_vectors)
            )
        self.omm_set_virtual_sites(system)
        return system

    #===================================================

    @needs_openmm
    def omm_add_constraints(self, system, constraints, rigidWater):
        """ Adds constraints to a given system

        Parameters
        ----------
        system : mm.System
            The OpenMM system for which constraints should be added
        constraints : None, app.HBonds, app.AllBonds, or app.HAngles
            Which kind of constraints should be used
        rigidWater : bool
            If True, water bonds are constrained regardless of whether
            constrains is None
        """
        if constraints is None and not rigidWater: return
        if constraints not in (None, app.HBonds, app.AllBonds, app.HAngles):
            raise ValueError("Unrecognized constraints option (%s)" %
                             constraints)
        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        # Rigid water only
        if constraints is None:
            for bond in self.bonds:
                # Skip all extra points... don't constrain those
                if isinstance(bond.atom1, ExtraPoint): continue
                if isinstance(bond.atom2, ExtraPoint): continue
                if (bond.atom1.residue.name in WATER_NAMES or
                        bond.atom2.residue.name in WATER_NAMES):
                    system.addConstraint(bond.atom1.idx, bond.atom2.idx,
                                         bond.type.req*length_conv)
            return
        # Other types of constraints
        for bond in self.bonds:
            if constraints is not app.HBonds or (bond.atom1.element == 1
                    or bond.atom2.element == 1):
                system.addConstraint(bond.atom1.idx, bond.atom2.idx,
                                     bond.type.req*length_conv)
        if constraints is app.HAngles:
            for angle in self.angles:
                num_h = (angle.atom1.element == 1) + (angle.atom3.element == 1)
                if num_h == 2 or (num_h == 1 and angle.atom2.element == 8):
                    # Constrain this angle
                    l1 = l2 = None
                    for bond in angle.atom2.bonds:
                        if bond in angle and angle.atom1 in bond:
                            l1 = bond.type.req * length_conv
                        elif bond in angle and angle.atom3 in bond:
                            l2 = bond.type.req * length_conv
                    # Law of cosines to find the constraint distance
                    if l1 is None or l2 is None: continue # no bonds found...
                    cost = math.cos(angle.type.theteq*DEG_TO_RAD)
                    length = math.sqrt(l1*l1 + l2*l2 - 2*l1*l2*cost)*length_conv
                    system.addConstraint(angle.atom1.idx,angle.atom3.idx,length)

    #===================================================

    @needs_openmm
    def omm_set_virtual_sites(self, system):
        """
        Sets the virtual sites in a given OpenMM `System` object from the extra
        points defined in this system

        Parameters
        ----------
        system : mm.System
            The system for which the virtual sites will be set. All particles
            must have already been added to this System before calling this
            method
        """
        if system.getNumParticles() != len(self.atoms):
            raise ValueError('OpenMM System does not correspond to Structure')
        for atom in self.atoms:
            if not isinstance(atom, ExtraPoint): continue
            # This is a virtual site... get its frame type
            typ = atom.frame_type
            weights = typ.get_weights()
            refatoms = typ.get_atoms()
            if isinstance(typ, TwoParticleExtraPointFrame):
                a1, a2 = refatoms
                w1, w2 = weights
                system.setVirtualSite(atom.idx,
                        mm.TwoParticleAverageSite(a1.idx, a2.idx, w1, w2)
                )
            elif isinstance(typ, ThreeParticleExtraPointFrame):
                a1, a2, a3 = refatoms
                w1, w2, w3 = weights
                system.setVirtualSite(atom.idx,
                        mm.ThreeParticleAverageSite(a1.idx, a2.idx, a3.idx,
                                                    w1, w2, w3)
                )
            elif isinstance(typ, OutOfPlaneExtraPointFrame):
                a1, a2, a3 = refatoms
                w1, w2, w3 = weights
                system.setVirtualSite(atom.idx,
                        mm.OutOfPlaneSite(a1.idx, a2.idx, a3.idx, w1, w2, w3)
                )

    #===================================================

    @needs_openmm
    def omm_bond_force(self, constraints=None, rigidWater=True,
                       flexibleConstraints=True):
        """
        Creates an OpenMM Bond Force object (or AmoebaBondForce if the bonds are
        for an Amoeba-parametrized system)

        Parameters
        ----------
        constraints : None, app.HBonds, app.AllBonds, or app.HAngles
            The types of constraints that are on the system. If
            flexibleConstraints is False, then the constrained bonds will not be
            added to the resulting Force
        rigidWater : bool=True
            Should water-H bonds be constrained regardless of `constraints`?
        flexibleConstraints : bool=True
            If True, all bonds are added to the force regardless of
            `constraints`

        Returns
        -------
        force
            HarmonicBondForce (or AmoebaBondForce if this is an Amoeba system),
            or None if there are no bonds to add
        """
        if not flexibleConstraints and constraints in (app.HAngles,
                app.AllBonds) or not self.bonds:
            return None # No bonds to add
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        _ambfrc = u.kilocalorie_per_mole/u.angstrom**2
        _ommfrc = u.kilojoule_per_mole/u.nanometer**2
        frc_conv = _ambfrc.conversion_factor_to(_ommfrc)
        # See if we need to add Amoeba bonds or regular bonds
        if (hasattr(self.bond_types, 'degree') and
                hasattr(self.bond_types, 'coeffs')):
            force = mm.AmoebaBondForce()
            force.setGlobalBondCubic(self.bond_types.coeffs[3]/length_conv)
            force.setGlobalBondQuartic(self.bond_types.coeffs[4]/length_conv**2)
        else:
            force = mm.HarmonicBondForce()
        force.setForceGroup(self.BOND_FORCE_GROUP)
        # Add the bonds
        for bond in self.bonds:
            if (bond.atom1.element == 1 or bond.atom2.element == 1) and (
                    not flexibleConstraints and constraints is app.HBonds):
                continue
            if bond.type is None:
                raise ParameterError('Cannot find necessary parameters')
            force.addBond(bond.atom1.idx, bond.atom2.idx,
                          bond.type.req*length_conv, 2*bond.type.k*frc_conv)
        # Done adding the force
        if force.getNumBonds() == 0:
            return None
        return force

    #===================================================

    @needs_openmm
    def omm_angle_force(self, constraints=None, flexibleConstraints=True):
        """
        Creates an OpenMM HarmonicAngleForce object (or AmoebaAngleForce if the
        angles are for an Amoeba-parametrized system)

        Parameters
        ----------
        constraints : None, app.HBonds, app.AllBonds, or app.HAngles
            The types of constraints that are on the system. If
            flexibleConstraints is False, then the constrained bonds will not be
            added to the resulting Force
        flexibleConstraints : bool=True
            If True, all bonds are added to the force regardless of
            `constraints`

        Returns
        -------
        force
            HarmonicAngleForce (or AmoebaAngleForce if this is an Amoeba
            system), or None if there are no angles to add
        """
        if not self.angles: return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        if (hasattr(self.angle_types, 'degree') and
                hasattr(self.angle_types, 'coeffs')):
            c = self.angle_types.coeffs
            force = mm.AmoebaAngleForce()
            force.setAmoebaGlobalAngleCubic(c[3])
            force.setAmoebaGlobalAngleQuartic(c[4])
            force.setAmoebaGlobalAnglePentic(c[5])
            force.setAmoebaGlobalAngleSextic(c[6])
        else:
            force = mm.HarmonicAngleForce()
        force.setForceGroup(self.ANGLE_FORCE_GROUP)
        for angle in self.angles:
            num_h = (angle.atom1.element == 1) + (angle.atom3.element == 1)
            if constraints is app.HAngles and (num_h == 2 or (num_h == 1 and
                    angle.atom2.element == 8) and not flexibleConstraints):
                continue
            if angle.type is None:
                raise ParameterError('Cannot find angle parameters')
            force.addAngle(angle.atom1.idx, angle.atom2.idx, angle.atom3.idx,
                           angle.type.theteq*DEG_TO_RAD, 2*angle.type.k*frc_conv)
        if force.getNumAngles() == 0:
            return None
        return force

    #===================================================

    @needs_openmm
    def omm_dihedral_force(self):
        """ Creates the OpenMM PeriodicTorsionForce modeling dihedrals

        Returns
        -------
        PeriodicTorsionForce
            Or None if no torsions are present in this system
        """
        if not self.dihedrals: return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        force = mm.PeriodicTorsionForce()
        force.setForceGroup(self.DIHEDRAL_FORCE_GROUP)
        for tor in self.dihedrals:
            if tor.type is None:
                raise ParameterError('Cannot find torsion parameters')
            if isinstance(tor.type, DihedralTypeList):
                for typ in tor.type:
                    force.addTorsion(tor.atom1.idx, tor.atom2.idx,
                                     tor.atom3.idx, tor.atom4.idx,
                                     int(typ.per), typ.phase*DEG_TO_RAD,
                                     typ.phi_k*frc_conv)
            else:
                force.addTorsion(tor.atom1.idx, tor.atom2.idx, tor.atom3.idx,
                                 tor.atom4.idx, int(tor.type.per),
                                 tor.type.phase*DEG_TO_RAD,
                                 tor.type.phi_k*frc_conv)
        return force

    #===================================================

    @needs_openmm
    def omm_rb_torsion_force(self):
        """ Creates the OpenMM RBTorsionForce for Ryckaert-Bellemans torsions

        Returns
        -------
        RBTorsionForce
            Or None if no torsions are present in this system
        """
        if not self.rb_torsions: return None
        conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        force = mm.RBTorsionForce()
        force.setForceGroup(self.RB_TORSION_FORCE_GROUP)
        for tor in self.rb_torsions:
            if tor.type is None:
                raise ParameterError('Cannot find R-B torsion parameters')
            force.addTorsion(tor.atom1.idx, tor.atom2.idx, tor.atom3.idx,
                             tor.atom4.idx, tor.type.c0*conv, tor.type.c1*conv,
                             tor.type.c2*conv, tor.type.c3*conv,
                             tor.type.c4*conv, tor.type.c5*conv)
        return force

    #===================================================

    @needs_openmm
    def omm_urey_bradley_force(self):
        """ Creates the OpenMM Urey-Bradley force

        Returns
        -------
        HarmonicBondForce
            Or None, if no urey-bradleys are present
        """
        if not self.urey_bradleys: return None
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        _ambfrc = u.kilocalorie_per_mole/u.angstrom**2
        _ommfrc = u.kilojoule_per_mole/u.nanometer**2
        frc_conv = _ambfrc.conversion_factor_to(_ommfrc)
        force = mm.HarmonicBondForce()
        force.setForceGroup(self.UREY_BRADLEY_FORCE_GROUP)
        for urey in self.urey_bradleys:
            if urey.type is None:
                raise ParameterError('Cannot find urey-bradley parameters')
            force.addBond(urey.atom1.idx, urey.atom2.idx,
                          urey.type.req*length_conv, 2*urey.type.k*frc_conv)
        return force

    #===================================================

    @needs_openmm
    def omm_improper_force(self):
        """ Creates the OpenMM improper torsion force (quadratic bias)

        Returns
        -------
        CustomTorsionForce
            With the formula k*(phi-phi0)^2, or None if there are no impropers
        """
        if not self.impropers: return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        force = mm.CustomTorsionForce("k*(theta-theta0)^2")
        force.addPerTorsionParameter('k')
        force.addPerTorsionParameter('theta0')
        force.setForceGroup(self.IMPROPER_FORCE_GROUP)
        for imp in self.impropers:
            if imp.type is None:
                raise ParameterError('Cannot find improper torsion '
                                       'parameters')
            force.addTorsion(imp.atom1.idx, imp.atom2.idx, imp.atom3.idx,
                             imp.atom4.idx, (imp.type.psi_k*frc_conv,
                             imp.type.psi_eq*DEG_TO_RAD))
        return force

    #===================================================

    @needs_openmm
    def omm_cmap_force(self):
        """ Creates the OpenMM CMAP torsion force

        Returns
        -------
        CMAPTorsionForce
            Or None, if no CMAP terms are present
        """
        if not self.cmaps: return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        force = mm.CMAPTorsionForce()
        force.setForceGroup(self.CMAP_FORCE_GROUP)
        # First get the list of cmap maps we're going to use. Just store the IDs
        # so we have simple integer comparisons to do later
        cmap_type_list = []
        cmap_map = dict()
        for cmap in self.cmaps:
            if cmap.type is None:
                raise ParameterError('Cannot find CMAP torsion parameters')
            if not id(cmap.type) in cmap_type_list:
                ct = cmap.type
                cmap_type_list.append(id(ct))
                # OpenMM takes the correction maps in the range 0 to 360
                # degrees, but we store it in -180 -- 180 (and in the transpose
                # of what OpenMM expects).
                grid = ct.grid.switch_range().T
                m = force.addMap(ct.resolution, [x*frc_conv for x in grid])
                cmap_map[id(ct)] = m
        # Now add all of the cmaps
        for cmap in self.cmaps:
            force.addTorsion(cmap_map[id(cmap.type)],
                             cmap.atom1.idx, cmap.atom2.idx, cmap.atom3.idx,
                             cmap.atom4.idx, cmap.atom2.idx, cmap.atom3.idx,
                             cmap.atom4.idx, cmap.atom5.idx)
        return force

    #===================================================

    @needs_openmm
    def omm_nonbonded_force(self, nonbondedMethod=None,
                            nonbondedCutoff=8*u.angstroms,
                            switchDistance=0*u.angstroms,
                            ewaldErrorTolerance=0.0005,
                            reactionFieldDielectric=78.5):
        """ Creates the OpenMM NonbondedForce instance

        Parameters
        ----------
        nonbondedMethod : cutoff method
            This is the cutoff method. It can be either the NoCutoff,
            CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald objects from the
            simtk.openmm.app namespace
        nonbondedCutoff : float or distance Quantity
            The nonbonded cutoff must be either a floating point number
            (interpreted as nanometers) or a Quantity with attached units. This
            is ignored if nonbondedMethod is NoCutoff.
        switchDistance : float or distance Quantity
            The distance at which the switching function is turned on for van
            der Waals interactions. This is ignored when no cutoff is used, and
            no switch is used if switchDistance is 0, negative, or greater than
            the cutoff
        ewaldErrorTolerance : float=0.0005
            When using PME or Ewald, the Ewald parameters will be calculated
            from this value
        reactionFieldDielectric : float=78.5
            If the nonbondedMethod is CutoffPeriodic or CutoffNonPeriodic, the
            region beyond the cutoff is treated using a reaction field method
            with this dielectric constant. It should be set to 1 if another
            implicit solvent model is being used (e.g., GB)

        Returns
        -------
        NonbondedForce
            This just implements the very basic NonbondedForce with the typical
            charge-charge and 12-6 Lennard-Jones interactions with the
            Lorentz-Berthelot combining rules.

        Notes
        -----
        This method calls update_dihedral_exclusions, which might alter the
        ``ignore_end`` attribute of some
        :class:`parmed.topologyobjects.Dihedral` instances based on any
        changes that might have been made to the bonds and angles since the last
        time it was called.

        Subclasses of Structure for which this nonbonded treatment is inadequate
        should override this method to implement what is needed.

        If nrexcl is set to 3 and no exception parameters are stored in the
        adjusts list, the 1-4 interactions are determined from the list of
        dihedrals
        """
        if not self.atoms: return None
        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        ene_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        force = mm.NonbondedForce()
        force.setForceGroup(self.NONBONDED_FORCE_GROUP)
        if u.is_quantity(nonbondedCutoff):
            nonbondedCutoff = nonbondedCutoff.value_in_unit(u.nanometers)
        if nonbondedMethod is None or nonbondedMethod is app.NoCutoff:
            force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        elif nonbondedMethod is app.CutoffNonPeriodic:
            force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
            force.setCutoffDistance(nonbondedCutoff)
        elif nonbondedMethod is app.CutoffPeriodic:
            force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
            force.setCutoffDistance(nonbondedCutoff)
        elif nonbondedMethod is app.PME:
            force.setNonbondedMethod(mm.NonbondedForce.PME)
            force.setCutoffDistance(nonbondedCutoff)
            force.setEwaldErrorTolerance(ewaldErrorTolerance)
        elif nonbondedMethod is app.Ewald:
            force.setNonbondedMethod(mm.NonbondedForce.Ewald)
            force.setCutoffDistance(nonbondedCutoff)
            force.setEwaldErrorTolerance(ewaldErrorTolerance)
        else:
            raise ValueError('Unrecognized nonbondedMethod (%s)' %
                             nonbondedMethod)
        force.setReactionFieldDielectric(reactionFieldDielectric)
        # Now add the particles
        sigma_scale = length_conv * 2 * 2**(-1/6)
        for atom in self.atoms:
            force.addParticle(atom.charge, atom.rmin*sigma_scale,
                              abs(atom.epsilon*ene_conv))
        # Add exclusions from the bond graph out to nrexcl-1 bonds away (atoms
        # nrexcl bonds away will be exceptions defined later)
        def exclude_to(origin, atom, level, end):
            if level >= end: return
            for partner in atom.bond_partners:
                if isinstance(partner, ExtraPoint): continue # Handled later
                if partner is origin: continue
                force.addException(origin.idx, atom.idx, 0.0, 0.5, 0.0, True)
                # Exclude EP children, too
                for child in origin.children:
                    force.addException(partner.idx, child.idx, 0.0, 0.5, 0.0,
                                       True)
                for child in partner.children:
                    force.addException(child.idx, origin.idx, 0.0, 0.5, 0.0,
                                       True)
                for child in origin.children:
                    for child2 in partner.children:
                        force.addException(child.idx, child2.idx, 0.0, 0.5,
                                           0.0, True)
                exclude_to(origin, partner, level+1, end)
        for atom in self.atoms:
            if isinstance(atom, ExtraPoint): continue # Handled separately
            exclude_to(atom, atom, 0, self.nrexcl)
        # Add the exceptions from the dihedral list IFF no explicit exceptions
        # (or *adjusts*) have been specified. If dihedral.ignore_end is False, a
        # 1-4 with the appropriate scaling factor is used as the exception.
        sigma_scale = 2**(-1/6) * length_conv
        if not self.adjusts:
            for dih in self.dihedrals:
                if dih.ignore_end: continue
                if isinstance(dih.type, DihedralTypeList):
                    scee = scnb = 0
                    i = 0
                    while (scee == 0 or scnb == 0) and i < len(dih.type):
                        scee = dih.type[i].scee
                        scnb = dih.type[i].scnb
                        i += 1
                    if scee == 0 or scnb == 0:
                        raise ValueError('Detected scaling constants of 0 for '
                                         'dihedral containing 1-4 info!')
                else:
                    scee = dih.type.scee
                    scnb = dih.type.scnb
                try:
                    rij, wdij, rij14, wdij14 = dih.atom1.atom_type.nbfix[
                                                    str(dih.atom4.atom_type)]
                except KeyError:
                    epsprod = abs(dih.atom1.epsilon_14 * dih.atom4.epsilon_14)
                    epsprod = math.sqrt(epsprod) * ene_conv / scnb
                    sigprod = (dih.atom1.rmin_14+dih.atom4.rmin_14)*sigma_scale
                else:
                    epsprod = wdij14 * ene_conv / scnb
                    sigprod = rij * length_conv * sigma_scale
                chgprod = dih.atom1.charge * dih.atom4.charge / scee
                force.addException(dih.atom1.idx, dih.atom4.idx, chgprod,
                                   sigprod, epsprod, True)
                for child in dih.atom1.children:
                    epsprod = abs(child.epsilon_14 * dih.atom4.epsilon_14)
                    epsprod = math.sqrt(epsprod) * ene_conv / scnb
                    sigprod = (child.rmin_14 + dih.atom4.rmin_14) * sigma_scale
                    chgprod = (child.charge * dih.atom4.charge) / scee
                    force.addException(child.idx, dih.atom4.idx, chgprod,
                                       sigprod, epsprod, True)
                for child in dih.atom4.children:
                    epsprod = abs(child.epsilon_14 * dih.atom1.epsilon_14)
                    epsprod = math.sqrt(epsprod) * ene_conv / scnb
                    sigprod = (child.rmin_14 + dih.atom1.rmin_14) * sigma_scale
                    chgprod = child.charge * dih.atom1.charge / scee
                    force.addException(child.idx, dih.atom1.idx, chgprod,
                                       sigprod, epsprod, True)
                for c1 in dih.atom1.children:
                    for c2 in dih.atom2.children:
                        epsprod = abs(c1.epsilon_14 * c2.epsilon_14)
                        epsprod = math.sqrt(epsprod) * ene_conv / scnb
                        sigprod = (c1.rmin_14 + c2.rmin_14) * sigma_scale
                        chgprod = c1.charge * c2.charge / scee
                        force.addException(c1.idx, c2.idx, chgprod, sigprod,
                                           epsprod, True)
        # Allow our specific exceptions (in adjusts) to override anything that
        # came before
        for pair in self.adjusts:
            chgprod = pair.atom1.charge * pair.atom2.charge * pair.type.chgscale
            force.addException(pair.atom1.idx, pair.atom2.idx, chgprod,
                               pair.type.rmin*sigma_scale,
                               pair.type.epsilon*ene_conv, True)

        # Any exclusion partners we already added will zero-out existing
        # exclusions/exceptions
        for a2 in atom.exclusion_partners:
            force.addException(atom.idx, a2.idx, 0.0, 0.5, 0.0, True)
            for c1 in atom.children:
                force.addException(c1.idx, a2.idx, 0.0, 0.5, 0.0, True)
            for c2 in a2.children:
                force.addException(atom.idx, c2.idx, 0.0, 0.5, 0.0, True)
            for c1 in atom.children:
                for c2 in a2.children:
                    force.addException(c1.idx, c2.idx, 0.0, 0.5, 0.0, True)

        if switchDistance and nonbondedMethod is not app.NoCutoff:
            if u.is_quantity(switchDistance):
                switchDistance = switchDistance.value_in_unit(u.nanometers)
            if 0 < switchDistance < nonbondedCutoff:
                force.setUseSwitchingFunction(True)
                force.setSwitchingDistance(switchDistance)

        # Now see if we have any NBFIXes that we need to implement via a
        # CustomNonbondedForce
        if self.has_NBFIX():
            return force, self._omm_nbfixed_force(force, nonbondedMethod)

        return force

    #===================================================

    def _omm_nbfixed_force(self, nonbfrc, nonbondedMethod):
        """ Private method for creating a CustomNonbondedForce with a lookup
        table. This should not be called by users -- you have been warned.

        Parameters
        ----------
        nonbfrc : NonbondedForce
            NonbondedForce for the "standard" nonbonded interactions. This will
            be modified (specifically, L-J ixns will be zeroed)
        nonbondedMethod : Nonbonded Method (e.g., NoCutoff, PME, etc.)
            The nonbonded method to apply here. Ewald and PME will be
            interpreted as CutoffPeriodic for the CustomNonbondedForce.

        Returns
        -------
        force : CustomNonbondedForce
            The L-J force with NBFIX implemented as a lookup table
        """
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        ene_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        # We need a CustomNonbondedForce to implement the NBFIX functionality.
        # First derive the type lookup tables
        lj_idx_list = [0 for atom in self.atoms]
        lj_radii, lj_depths = [], []
        num_lj_types = 0
        lj_type_list = []
        for i, atom in enumerate(self.atoms):
            atype = atom.atom_type
            if lj_idx_list[i]: continue # already assigned
            num_lj_types += 1
            lj_idx_list[i] = num_lj_types
            ljtype = (atype.rmin, abs(atype.epsilon))
            lj_type_list.append(atype)
            lj_radii.append(atype.rmin)
            lj_depths.append(abs(atype.epsilon))
            for j in range(i+1, len(self.atoms)):
                atype2 = self.atoms[j].atom_type
                if lj_idx_list[j] > 0: continue # already assigned
                if atype2 is atype:
                    lj_idx_list[j] = num_lj_types
                elif not atype.nbfix:
                    # Only non-NBFIXed atom types can be compressed
                    ljtype2 = (atype2.rmin, abs(atype2.epsilon))
                    if ljtype == ljtype2:
                        lj_idx_list[j] = num_lj_types
        # Now everything is assigned. Create the A- and B-coefficient arrays
        acoef = [0 for i in range(num_lj_types*num_lj_types)]
        bcoef = acoef[:]
        for i in range(num_lj_types):
            for j in range(num_lj_types):
                namej = lj_type_list[j].name
                try:
                    rij, wdij, rij14, wdij14 = lj_type_list[i].nbfix[namej]
                except KeyError:
                    rij = (lj_radii[i] + lj_radii[j]) * length_conv
                    wdij = math.sqrt(lj_depths[i] * lj_depths[j]) * ene_conv
                else:
                    rij *= length_conv
                    wdij *= ene_conv
                rij6 = rij**6
                acoef[i+num_lj_types*j] = math.sqrt(wdij) * rij6
                bcoef[i+num_lj_types*j] = 2 * wdij * rij6
        force = mm.CustomNonbondedForce('(a/r6)^2-b/r6; r6=r2*r2*r2; r2=r^2; '
                                        'a=acoef(type1, type2); '
                                        'b=bcoef(type1, type2)')
        force.addTabulatedFunction('acoef',
                mm.Discrete2DFunction(num_lj_types, num_lj_types, acoef))
        force.addTabulatedFunction('bcoef',
                mm.Discrete2DFunction(num_lj_types, num_lj_types, bcoef))
        force.addPerParticleParameter('type')
        force.setForceGroup(self.NONBONDED_FORCE_GROUP)
        if (nonbondedMethod is app.PME or nonbondedMethod is app.Ewald or
                nonbondedMethod is app.CutoffPeriodic):
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        elif nonbondedMethod is app.NoCutoff:
            force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
        elif nonbondedMethod is app.CutoffNonPeriodic:
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        else:
            raise ValueError('Unrecognized nonbonded method [%s]' %
                             nonbondedMethod)
        # Add the particles
        for i in lj_idx_list:
            force.addParticle((i-1,))
        # Now wipe out the L-J parameters in the nonbonded force
        for i in range(nonbfrc.getNumParticles()):
            chg, sig, eps = nonbfrc.getParticleParameters(i)
            nonbfrc.setParticleParameters(i, chg, 0.5, 0.0)
        # Now transfer the exclusions
        for ii in range(nonbfrc.getNumExceptions()):
            i, j, qq, ss, ee = nonbfrc.getExceptionParameters(ii)
            force.addExclusion(i, j)
        # Now transfer the other properties (cutoff, switching function, etc.)
        force.setUseLongRangeCorrection(True)
        if nonbondedMethod is app.NoCutoff:
            force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
        elif nonbondedMethod is app.CutoffNonPeriodic:
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        elif nonbondedMethod in (app.PME, app.Ewald, app.CutoffPeriodic):
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        else:
            raise ValueError('Unsupported nonbonded method %s' %
                             nonbondedMethod)
        force.setCutoffDistance(nonbfrc.getCutoffDistance())
        if nonbfrc.getUseSwitchingFunction():
            force.setUseSwitchingFunction(True)
            force.setSwitchingDistance(nonbfrc.getSwitchingDistance())

        return force

    #===================================================

    @needs_openmm
    def omm_gbsa_force(self, implicitSolvent,
                       nonbondedMethod=None,
                       nonbondedCutoff=30.0*u.angstroms,
                       soluteDielectric=1.0,
                       solventDielectric=78.5,
                       implicitSolventKappa=None,
                       implicitSolventSaltConc=0.0*u.moles/u.liter,
                       temperature=298.15*u.kelvin,
                       useSASA=True):
        """
        Creates a Generalized Born force for running implicit solvent
        calculations

        Parameters
        ----------
        implicitSolvent : app.HCT, app.OBC1, app.OBC2, app.GBn, app.GBn2
            The Generalized Born implicit solvent model to use.
        nonbondedMethod : cutoff method
            This is the cutoff method. It can be either the NoCutoff,
            CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald objects from the
            simtk.openmm.app namespace. Default is NoCutoff
        nonbondedCutoff : float or distance Quantity
            The nonbonded cutoff must be either a floating opint number
            (interpreted as nanometers) or a Quantity with attached units. This
            is ignored if nonbondedMethod is NoCutoff
        implicitSolventKappa : float or 1/distance Quantity = None
            This is the Debye kappa property related to modeling saltwater
            conditions in GB. It should have units of 1/distance (1/nanometers
            is assumed if no units present). A value of None means that kappa
            will be calculated from implicitSolventSaltConc (below)
        implicitSolventSaltConc : float or amount/volume Quantity=0 moles/liter
            If implicitSolventKappa is None, the kappa will be computed from the
            salt concentration. It should have units compatible with mol/L
        temperature : float or temperature Quantity = 298.15 kelvin
            This is only used to compute kappa from implicitSolventSaltConc
        soluteDielectric : float=1.0
            The dielectric constant of the protein interior used in GB
        solventDielectric : float=78.5
            The dielectric constant of the water used in GB
        """
        from simtk.openmm.app.internal.customgbforces import (GBSAHCTForce,
                GBSAOBC1Force, GBSAOBC2Force, GBSAGBnForce, GBSAGBn2Force,
                convertParameters)
        if implicitSolvent is None: return None
        if useSASA:
            sasa = 'ACE'
        else:
            sasa = None
        if nonbondedMethod is None:
            nonbondedMethod = app.NoCutoff
        if implicitSolvent not in (app.HCT, app.OBC1, app.OBC2, app.GBn,
                app.GBn2):
            raise ValueError('Unrecognized implicit solvent model')
        gb_parms = convertParameters(self._get_gb_parameters(implicitSolvent),
                                     str(implicitSolvent))
        if implicitSolventKappa is None:
            if u.is_quantity(implicitSolventSaltConc):
                sc = implicitSolventSaltConc.value_in_unit(u.moles/u.liter)
                implicitSolventSaltConc = sc
            if u.is_quantity(temperature):
                temperature = temperature.value_in_unit(u.kelvin)
            # The constant is 1 / sqrt(eps_0 * kB / (2*NA*q^2*1000)) where NA is
            # Avogadro's number, eps_0 is the permittivity of free space, q is
            # the charge (this # matches Amber's conversion factor)
            implicitSolventKappa = 50.33355 * math.sqrt(implicitSolventSaltConc
                                          / solventDielectric / temperature)
            # Multiply by 0.73 to account for ion exclusions, and multiply by 10
            # to convert to 1/nm from 1/angstroms
            implicitSolventKappa *= 7.3
        elif u.is_quantity(implicitSolventKappa):
            implicitSolventKappa = implicitSolventKappa.value_in_unit(
                    u.nanometer**-1)

        if nonbondedMethod is app.NoCutoff:
            cutoff = None
        elif u.is_quantity(nonbondedCutoff):
            cutoff = nonbondedCutoff.value_in_unit(u.nanometers)
        else:
            cutoff = nonbondedCutoff
        if implicitSolvent is app.HCT:
            force = GBSAHCTForce(solventDielectric, soluteDielectric, sasa,
                                 cutoff, kappa=implicitSolventKappa)
        elif implicitSolvent is app.OBC1:
            force = GBSAOBC1Force(solventDielectric, soluteDielectric, sasa,
                                  cutoff, kappa=implicitSolventKappa)
        elif implicitSolvent is app.OBC2:
            force = GBSAOBC2Force(solventDielectric, soluteDielectric, sasa,
                                  cutoff, kappa=implicitSolventKappa)
        elif implicitSolvent is app.GBn:
            force = GBSAGBnForce(solventDielectric, soluteDielectric, sasa,
                                 cutoff, kappa=implicitSolventKappa)
        elif implicitSolvent is app.GBn2:
            force = GBSAGBn2Force(solventDielectric, soluteDielectric, sasa,
                                  cutoff, kappa=implicitSolventKappa)
        else:
            raise ValueError('Unexpected implicit solvent model... '
                             'should not be here')
        for i, atom in enumerate(self.atoms):
            force.addParticle([atom.charge] + list(gb_parms[i]))
        # Set cutoff method
        if nonbondedMethod is app.NoCutoff:
            force.setNonbondedMethod(mm.CustomGBForce.NoCutoff)
        elif nonbondedMethod is app.CutoffNonPeriodic:
            force.setNonbondedMethod(mm.CustomGBForce.CutoffNonPeriodic)
            force.setCutoffDistance(cutoff)
        else: # cutoff periodic (PME, CutoffPeriodic, Ewald)
            force.setNonbondedMethod(mm.CustomGBForce.CutoffPeriodic)
            force.setCutoffDistance(cutoff)
        force.setForceGroup(self.NONBONDED_FORCE_GROUP)

        return force

    #===================================================

    # Amoeba-specific forces below

    @needs_openmm
    def omm_trigonal_angle_force(self):
        """ Creates the Amoeba trigonal-angle force

        Returns
        -------
        AmoebaInPlaneAngleForce
            The trigonal in-plane Angle force
        """
        if not self.trigonal_angles: return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        if (not hasattr(self.trigonal_angle_types, 'degree') or not
                hasattr(self.trigonal_angle_types, 'coeffs')):
            raise ParameterError('Do not have the trigonal angle force '
                                   'table parameters')
        force = mm.AmoebaInPlaneAngleForce()
        c = self.trigonal_angle_types.coeffs
        force.setAmoebaGlobalInPlaneAngleCubic(c[3])
        force.setAmoebaGlobalInPlaneAngleQuartic(c[4])
        force.setAmoebaGlobalInPlaneAngleQuintic(c[5])
        force.setAmoebaGlobalInPlaneAngleSextic(c[6])
        force.setForceGroup(self.TRIGONAL_ANGLE_FORCE_GROUP)
        for ang in self.trigonal_angles:
            if ang.type is None:
                raise ParameterError('Missing trigonal angle parameters')
            force.addAngle(ang.atom1.idx, ang.atom2.idx, ang.atom3.idx,
                           ang.atom4.idx, ang.type.theteq,
                           ang.type.k*frc_conv)
        return force

    #===================================================

    @needs_openmm
    def omm_out_of_plane_bend_force(self):
        """ Creates the Amoeba out-of-plane bend force

        Returns
        -------
        AmoebaOutOfPlaneBendForce
            The out-of-plane bend Angle force
        """
        if not self.out_of_plane_bends: return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        if (not hasattr(self.out_of_plane_bend_types, 'degree') or not
                hasattr(self.out_of_plane_bend_types, 'coeffs')):
            raise ParameterError('Do not have the trigonal angle force '
                                   'table parameters')
        force = mm.AmoebaOutOfPlaneBendForce()
        c = self.out_of_plane_bend_types.coeffs
        force.setAmoebaGlobalOutOfPlaneBendCubic(c[3])
        force.setAmoebaGlobalOutOfPlaneBendQuartic(c[4])
        force.setAmoebaGlobalOutOfPlaneBendQuintic(c[5])
        force.setAmoebaGlobalOutOfPlaneBendSextic(c[6])
        force.setForceGroup(self.OUT_OF_PLANE_BEND_FORCE_GROUP)
        for ang in self.out_of_plane_bends:
            if ang.type is None:
                raise ParameterError('Missing out-of-plane bend parameters')
            force.addOutOfPlaneBend(ang.atom1.idx, ang.atom2.idx, ang.atom3.idx,
                                    ang.atom4.idx, 2*ang.type.k*frc_conv)
        return force

    #===================================================

    @needs_openmm
    def omm_pi_torsion_force(self):
        """ Creates the Amoeba pi-torsion force

        Returns
        -------
        AmoebaPiTorsionForce
            The pi-torsion force
        """
        if not self.pi_torsions: return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        force = mm.AmoebaPiTorsionForce()
        force.setForceGroup(self.PI_TORSION_FORCE_GROUP)
        for ang in self.pi_torsions:
            if ang.type is None:
                raise ParameterError('Missing pi-torsion parameters')
            force.addPiTorsion(ang.atom1.idx, ang.atom2.idx, ang.atom3.idx,
                               ang.atom4.idx, ang.atom5.idx, ang.atom6.idx,
                               ang.type.phi_k*frc_conv)
        return force

    #===================================================

    @needs_openmm
    def omm_stretch_bend_force(self):
        """ Create the OpenMM Amoeba stretch-bend force for this system

        Returns
        -------
        AmoebaStretchBendForce
            The stretch-bend force containing all terms in this system
        """
        if not self.stretch_bends: return None
        # Conversion factor taken from pyopenmm/processTinkerForceField.py
        frc_conv = math.pi / 180 * 41.84 # 4.184 * 10
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        force = mm.AmoebaStretchBendForce()
        force.setForceGroup(self.STRETCH_BEND_FORCE_GROUP)
        for strbnd in self.stretch_bends:
            if strbnd.type is None:
                raise ParameterError("Missing stretch-bend parameters")
            force.addStretchBend(strbnd.atom1.idx, strbnd.atom2.idx,
                                 strbnd.atom3.idx, strbnd.type.req1*length_conv,
                                 strbnd.type.req2*length_conv,
                                 strbnd.type.k*frc_conv)

    #===================================================

    @needs_openmm
    def omm_torsion_torsion_force(self):
        """ Create the OpenMM Amoeba coupled-torsion (CMAP) force

        Returns
        -------
        AmoebaTorsionTorsionForce
            The torsion-torsion (CMAP) force with all coupled-torsion parameters
            for this system
        """
        if not self.torsion_torsions: return None
        # Not implemented yet...
        raise NotImplementedError("Torsion-torsions not yet implemented")

    #===================================================

    def _prune_empty_bonds(self):
        """ Gets rid of any empty bonds """
        for i in reversed(range(len(self.bonds))):
            bond = self.bonds[i]
            if bond.atom1 is None and bond.atom2 is None:
                del self.bonds[i]
            elif bond.atom1.idx == -1 or bond.atom2.idx == -1:
                bond.delete()
                del self.bonds[i]

    #===================================================

    def _prune_empty_angles(self):
        """ Gets rid of any empty angles """
        for i in reversed(range(len(self.angles))):
            angle = self.angles[i]
            if (angle.atom1 is None and angle.atom2 is None and
                    angle.atom3 is None):
                del self.angles[i]
            elif (angle.atom1.idx == -1 or angle.atom2.idx == -1 or
                    angle.atom3.idx == -1):
                angle.delete()
                del self.angles[i]

    #===================================================

    def _prune_empty_dihedrals(self):
        """ Gets rid of any empty dihedrals """
        for i in reversed(range(len(self.dihedrals))):
            dihed = self.dihedrals[i]
            if (dihed.atom1 is None and dihed.atom2 is None and
                    dihed.atom3 is None and dihed.atom4 is None):
                del self.dihedrals[i]
            elif (dihed.atom1.idx == -1 or dihed.atom2.idx == -1 or
                    dihed.atom3.idx == -1 or dihed.atom4.idx == -1):
                dihed.delete()
                del self.dihedrals[i]

    #===================================================

    def _prune_empty_rb_torsions(self):
        """ Gets rid of any empty R-B torsions """
        for i in reversed(range(len(self.rb_torsions))):
            dihed = self.rb_torsions[i]
            if (dihed.atom1 is None and dihed.atom2 is None and
                    dihed.atom3 is None and dihed.atom4 is None):
                del self.rb_torsions[i]
            elif (dihed.atom1.idx == -1 or dihed.atom2.idx == -1 or
                    dihed.atom3.idx == -1 or dihed.atom4.idx == -1):
                dihed.delete()
                del self.rb_torsions[i]

    #===================================================

    def _prune_empty_ureys(self):
        """ Gets rid of any empty Urey-Bradley terms """
        for i in reversed(range(len(self.urey_bradleys))):
            ub = self.urey_bradleys[i]
            if ub.atom1 is None and ub.atom2 is None:
                del self.urey_bradleys[i]
            elif ub.atom1.idx == -1 or ub.atom2.idx == -1:
                ub.delete()
                del self.urey_bradleys[i]

    #===================================================

    def _prune_empty_impropers(self):
        """ Gets rid of any empty improper torsions """
        for i in reversed(range(len(self.impropers))):
            imp = self.impropers[i]
            if (imp.atom1 is None and imp.atom2 is None and imp.atom3 is None
                    and imp.atom4 is None):
                del self.impropers[i]
            elif (imp.atom1.idx == -1 or imp.atom2.idx == -1 or
                    imp.atom3.idx == -1 or imp.atom4.idx == -1):
                imp.delete()
                del self.impropers[i]

    #===================================================

    def _prune_empty_cmaps(self):
        """ Gets rid of any empty CMAP terms """
        for i in reversed(range(len(self.cmaps))):
            cmap = self.cmaps[i]
            if (cmap.atom1 is None and cmap.atom2 is None and cmap.atom3 is None
                    and cmap.atom4 is None and cmap.atom5 is None):
                del self.cmaps[i]
            elif (cmap.atom1.idx == -1 or cmap.atom2.idx == -1 or
                    cmap.atom3.idx == -1 or cmap.atom4.idx == -1 or
                    cmap.atom5.idx == -1):
                cmap.delete()
                del self.cmaps[i]

    #===================================================

    def _prune_empty_trigonal_angles(self):
        """ Gets rid of any empty trigonal angles """
        for i in reversed(range(len(self.trigonal_angles))):
            ta = self.trigonal_angles[i]
            if (ta.atom1 is None and ta.atom2 is None and ta.atom3 is None and
                    ta.atom4 is None):
                del self.trigonal_angles[i]
            elif (ta.atom1.idx == -1 or ta.atom2.idx == -1 or
                    ta.atom3.idx == -1 or ta.atom4.idx == -1):
                # Not stored anywhere, no need to call delete()
                del self.trigonal_angles[i]

    #===================================================

    def _prune_empty_out_of_plane_bends(self):
        """ Gets rid of any empty out-of-plane bends """
        for i in reversed(range(len(self.out_of_plane_bends))):
            oop = self.out_of_plane_bends[i]
            if (oop.atom1 is None and oop.atom2 is None and oop.atom3 is None
                    and oop.atom4 is None):
                del self.out_of_plane_bends[i]
            elif (oop.atom1.idx == -1 or oop.atom2.idx == -1 or
                    oop.atom3.idx == -1 or oop.atom4.idx == -1):
                # Not stored anywhere, no need to call delete()
                del self.out_of_plane_bends[i]

    #===================================================

    def _prune_empty_pi_torsions(self):
        """ Gets rid of any empty pi-torsions """
        for i in reversed(range(len(self.pi_torsions))):
            pit = self.pi_torsions[i]
            if (pit.atom1 is None and pit.atom2 is None and
                    pit.atom3 is None and pit.atom4 is None and
                    pit.atom5 is None and pit.atom6 is None):
                del self.pi_torsions[i]
            elif (pit.atom1.idx == -1 or pit.atom2.idx == -1 or
                    pit.atom3.idx == -1 or pit.atom4.idx == -1 or
                    pit.atom5.idx == -1 or pit.atom6.idx == -1):
                # Not stored anywhere, no need to call delete()
                del self.pi_torsions[i]

    #===================================================

    def _prune_empty_stretch_bends(self):
        """ Gets rid of any empty stretch-bend terms """
        for i in reversed(range(len(self.stretch_bends))):
            sb = self.stretch_bends[i]
            if sb.atom1 is None and sb.atom2 is None and sb.atom3 is None:
                del self.stretch_bends[i]
            elif (sb.atom1.idx == -1 or sb.atom2.idx == -1 or
                    sb.atom3.idx == -1):
                # Not stored anywhere, no need to call delete()
                del self.stretch_bends[i]

    #===================================================

    def _prune_empty_torsion_torsions(self):
        """ Gets rid of any empty torsion-torsion terms """
        for i in reversed(range(len(self.torsion_torsions))):
            tt = self.torsion_torsions[i]
            if (tt.atom1 is None and tt.atom2 is None and tt.atom3 is None
                    and tt.atom4 is None and tt.atom5 is None):
                del self.torsion_torsions[i]
            elif (tt.atom1.idx == -1 or tt.atom2.idx == -1 or
                    tt.atom3.idx == -1 or tt.atom4.idx == -1 or
                    tt.atom5.idx == -1):
                tt.delete()
                del self.torsion_torsions[i]

    #===================================================

    def _prune_empty_chiral_frames(self):
        """ Gets rid of any empty chiral frame terms """
        for i in reversed(range(len(self.chiral_frames))):
            cf = self.chiral_frames[i]
            if cf.atom1 is None or cf.atom2 is None:
                del self.chiral_frames[i]
            elif cf.atom1.idx == -1 or cf.atom2.idx == -1:
                del self.chiral_frames[i]

    #===================================================

    def _prune_empty_multipole_frames(self):
        """ Gets rid of any empty multipole frame terms """
        for i in reversed(range(len(self.multipole_frames))):
            mf = self.multipole_frames[i]
            if mf.atom is None or mf.atom.idx == -1:
                del self.multipole_frames[i]

    #===================================================

    def _prune_empty_adjusts(self):
        """ Gets rid of any empty nonbonded exception adjustments """
        for i in reversed(range(len(self.adjusts))):
            adj = self.adjusts[i]
            if adj.atom1 is None or adj.atom2 is None:
                del self.adjusts[i]
            elif adj.atom1.idx == -1 or adj.atom2.idx == -1:
                del self.adjusts[i]

    #===================================================

    @needs_openmm
    def _get_gb_parameters(self, implicitSolvent):
        """ Gets the GB parameters for the requested GB model used by OpenMM

        Parameters
        ----------
        implicitSolvent : app.HCT, app.OBC1, app.OBC2, app.GBn, or app.GBn2
            The object specifying a particular GB model in OpenMM

        Returns
        -------
        parameters : list of float
            List of parameters for the requested GB model
        """
        if implicitSolvent is app.GBn:
            screen = [0.5 for atom in self.atoms]
            radii = [atom.radii for atom in self.atoms]
            for i, atom in enumerate(self.atoms):
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
                if radii[i] == 0:
                    radii[i] = _bondi(atom)
        elif implicitSolvent is app.GBn2:
            # Add non-optimized values as defaults
            alpha = [1.0 for i in self.atoms]
            beta = [0.8 for i in self.atoms]
            gamma = [4.85 for i in self.atoms]
            screen = [0.5 for i in self.atoms]
            radii = [atom.radii for atom in self.atoms]
            for i, atom in enumerate(self.atoms):
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
                if not radii[i]:
                    radii[i] = _mbondi3(atom)
        else:
            radii = [atom.radii for atom in self.atoms]
            screen = [atom.screen for atom in self.atoms]
            for i, atom in enumerate(self.atoms):
                if not radii[i] or not screen[i]:
                    # Replace with defaults
                    radii[i], screen[i] = _gb_rad_screen(atom, implicitSolvent)

        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        radii = [x * length_conv for x in radii]

        if implicitSolvent is app.GBn2:
            return list(zip(radii, screen, alpha, beta, gamma))
        return list(zip(radii, screen))

    #===================================================

    def __str__(self):
        if hasattr(self, 'name') and self.name:
            return self.name
        return repr(self)

    #===================================================

    # Support concatenating and replicating the Structure

    def __add__(self, other):
        cp = copy(self)
        cp += other
        return cp

    def __iadd__(self, other):
        if not isinstance(other, Structure):
            raise TypeError('Must concatenate with Structure')
        # The basic approach taken here is to extend the atom list, then scan
        # through all of the valence terms in `other`, adding them to the
        # corresponding arrays of `self`, using an offset to look into the atom
        # array. The types are done a little differently. The types are all
        # copied into a "temporary" array so they have the same index as the
        # original Structure `other`, and as the valence terms are created the
        # reference to that type is added. Afterwards, those types are tacked on
        # to the end of the types list on `self`
        aoffset = len(self.atoms)
        try:
            roffset = self.residues[-1].number + 1
        except IndexError:
            # No residues... must be an empty structure
            roffset = 0
        for atom in other.atoms:
            res = atom.residue
            self.add_atom(copy(atom), res.name, res.idx+roffset, res.chain,
                          res.insertion_code)
        def copy_valence_terms(oval, otyp, sval, styp, attrlist):
            """ Copies the valence terms from one list to another;
            oval=Other VALence; otyp=Other TYPe; sval=Self VALence;
            styp=Self TYPe; attrlist=ATTRibute LIST (atom1, atom2, ...)
            """
            otypcp = [copy(typ) for typ in otyp]
            for val in oval:
                ats = [getattr(val, attr) for attr in attrlist]
                # Replace any atoms with the ones from self.atoms
                for i, at in enumerate(ats):
                    if isinstance(at, Atom):
                        ats[i] = self.atoms[at.idx+aoffset]
                kws = dict()
                if otypcp and val.type is not None:
                    kws['type'] = otypcp[val.type.idx]
                sval.append(type(val)(*ats, **kws))
                if hasattr(val, 'funct'):
                    sval[-1].funct = val.funct
            # Now tack on the "new" types copied from `other`
            styp.extend(otypcp)
            if hasattr(styp, 'claim'):
                styp.claim()
        copy_valence_terms(other.bonds, other.bond_types, self.bonds,
                           self.bond_types, ['atom1', 'atom2'])
        copy_valence_terms(other.angles, other.angle_types, self.angles,
                           self.angle_types, ['atom1', 'atom2', 'atom3'])
        copy_valence_terms(other.dihedrals, other.dihedral_types,
                           self.dihedrals, self.dihedral_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'improper',
                           'ignore_end'])
        copy_valence_terms(other.rb_torsions, other.rb_torsion_types,
                           self.rb_torsions, self.rb_torsion_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'improper',
                           'ignore_end'])
        copy_valence_terms(other.urey_bradleys, other.urey_bradley_types,
                           self.urey_bradleys, self.urey_bradley_types,
                           ['atom1', 'atom2'])
        copy_valence_terms(other.impropers, other.improper_types,
                           self.impropers, self.improper_types,
                           ['atom1', 'atom2', 'atom3', 'atom4'])
        copy_valence_terms(other.cmaps, other.cmap_types,
                           self.cmaps, self.cmap_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
        copy_valence_terms(other.trigonal_angles, other.trigonal_angle_types,
                           self.trigonal_angles, self.trigonal_angle_types,
                           ['atom1', 'atom2', 'atom3', 'atom4'])
        copy_valence_terms(other.out_of_plane_bends,
                           other.out_of_plane_bend_types,
                           self.out_of_plane_bends,
                           self.out_of_plane_bend_types,
                           ['atom1', 'atom2', 'atom3', 'atom4'])
        copy_valence_terms(other.pi_torsions, other.pi_torsion_types,
                           self.pi_torsions, self.pi_torsion_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'atom5',
                            'atom6'])
        copy_valence_terms(other.stretch_bends, other.stretch_bend_types,
                           self.stretch_bends, self.stretch_bend_types,
                           ['atom1', 'atom2', 'atom3'])
        copy_valence_terms(other.torsion_torsions, other.torsion_torsion_types,
                           self.torsion_torsions, self.torsion_torsion_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
        copy_valence_terms(other.chiral_frames, [], self.chiral_frames, [],
                           ['atom1', 'atom2', 'chirality'])
        copy_valence_terms(other.multipole_frames, [], self.multipole_frames,
                           [], ['atom', 'frame_pt_num', 'vectail', 'vechead',
                           'nvec'])
        copy_valence_terms(other.adjusts, other.adjust_types, self.adjusts,
                           self.adjust_types, ['atom1', 'atom2'])
        copy_valence_terms(other.donors, [], self.donors, [],
                           ['atom1', 'atom2'])
        copy_valence_terms(other.acceptors, [], self.acceptors, [],
                           ['atom1', 'atom2'])
        copy_valence_terms(other.groups, [], self.groups, [],
                           ['bs', 'type', 'move'])
        if self._coordinates is None or other._coordinates is None:
            self._coordinates = None
        elif self._coordinates.shape[0] != other._coordinates.shape[0]:
            self._coordinates = None
        else:
            np.concatenate((self._coordinates, other.coordinates))
        return self

    def __mul__(self, ncopies):
        """ Replicates the current Structure `ncopies` times """
        cp = copy(self)
        return cp.__imul__(ncopies, self)

    def __imul__(self, ncopies, other=None):
        if not isinstance(ncopies, integer_types):
            raise TypeError('Can only multiply a structure by an integer')
        # The basic approach here is similar to what we used in __iadd__, except
        # we don't have to extend the type arrays at all -- we just point to the
        # same one that the first copy pointed to.
        def copy_valence_terms(oval, aoffset, sval, styp, attrlist):
            """ Copies the valence terms from one list to another;
            oval=Other VALence; otyp=Other TYPe; sval=Self VALence;
            styp=Self TYPe; attrlist=ATTRibute LIST (atom1, atom2, ...)
            """
            for val in oval:
                ats = [getattr(val, attr) for attr in attrlist]
                # Replace any atoms with the ones from self.atoms
                for i, at in enumerate(ats):
                    if isinstance(at, Atom):
                        ats[i] = self.atoms[at.idx+aoffset]
                kws = dict()
                if styp and val.type is not None:
                    kws['type'] = styp[val.type.idx]
                sval.append(type(val)(*ats, **kws))
                if hasattr(val, 'funct'):
                    sval[-1].funct = val.funct
        if other is None: other = copy(self)
        for i in range(ncopies-1):
            aoffset = len(self.atoms)
            roffset = self.residues[-1].number + 1
            for atom in other.atoms:
                res = atom.residue
                self.add_atom(copy(atom), res.name, res.idx+roffset, res.chain,
                              res.insertion_code)
            copy_valence_terms(other.bonds, aoffset, self.bonds,
                               self.bond_types, ['atom1', 'atom2'])
            copy_valence_terms(other.angles, aoffset, self.angles,
                               self.angle_types, ['atom1', 'atom2', 'atom3'])
            copy_valence_terms(other.dihedrals, aoffset, self.dihedrals,
                               self.dihedral_types, ['atom1', 'atom2', 'atom3',
                               'atom4', 'improper', 'ignore_end'])
            copy_valence_terms(other.rb_torsions, aoffset,
                               self.rb_torsions, self.rb_torsion_types,
                               ['atom1', 'atom2', 'atom3', 'atom4', 'improper',
                                'ignore_end'])
            copy_valence_terms(other.urey_bradleys, aoffset,
                               self.urey_bradleys, self.urey_bradley_types,
                               ['atom1', 'atom2'])
            copy_valence_terms(other.impropers, aoffset,
                               self.impropers, self.improper_types,
                               ['atom1', 'atom2', 'atom3', 'atom4'])
            copy_valence_terms(other.cmaps, aoffset,
                               self.cmaps, self.cmap_types,
                               ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
            copy_valence_terms(other.trigonal_angles, aoffset,
                               self.trigonal_angles, self.trigonal_angle_types,
                               ['atom1', 'atom2', 'atom3', 'atom4'])
            copy_valence_terms(other.out_of_plane_bends,
                               aoffset, self.out_of_plane_bends,
                               self.out_of_plane_bend_types,
                               ['atom1', 'atom2', 'atom3', 'atom4'])
            copy_valence_terms(other.pi_torsions, aoffset,
                               self.pi_torsions, self.pi_torsion_types,
                               ['atom1', 'atom2', 'atom3', 'atom4', 'atom5',
                                'atom6'])
            copy_valence_terms(other.stretch_bends, aoffset,
                               self.stretch_bends, self.stretch_bend_types,
                               ['atom1', 'atom2', 'atom3'])
            copy_valence_terms(other.torsion_torsions, aoffset,
                               self.torsion_torsions, self.torsion_torsion_types,
                               ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
            copy_valence_terms(other.chiral_frames, aoffset, self.chiral_frames,
                               [], ['atom1', 'atom2', 'chirality'])
            copy_valence_terms(other.multipole_frames, aoffset,
                               self.multipole_frames, [],
                               ['atom', 'frame_pt_num', 'vectail', 'vechead',
                               'nvec'])
            copy_valence_terms(other.adjusts, aoffset, self.adjusts,
                               self.adjust_types, ['atom1', 'atom2'])
            copy_valence_terms(other.donors, aoffset, self.donors, [],
                               ['atom1', 'atom2'])
            copy_valence_terms(other.acceptors, aoffset, self.acceptors, [],
                               ['atom1', 'atom2'])
            copy_valence_terms(other.groups, aoffset, self.groups, [],
                               ['bs', 'type', 'move'])
        if self._coordinates is not None:
            coords = np.tile(self._coordinates.ravel(), ncopies).reshape(
                    (-1, len(self.atoms), 3))
            self._coordinates = coords
        return self

    #===================================================

    @staticmethod
    def _add_force_to_system(system, force):
        """ Adds an OpenMM force to a system IFF the force is not None """
        if force is None: return
        if isinstance(force, tuple) or isinstance(force, list):
            # It's possible we got multiple forces to add
            for f in force:
                system.addForce(f)
            return
        system.addForce(force)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
