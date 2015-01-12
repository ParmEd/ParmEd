"""
This module contains the core base class for all of the chemical structures with
various topological and force field features.

Author: Jason Swails
Date: November 10, 2014
"""
try:
    import bz2
except ImportError:
    bz2 = None
from chemistry.exceptions import PDBError, PDBWarning
from chemistry.periodic_table import AtomicNum, Mass, Element
from chemistry.topologyobjects import *
import copy
try:
    import gzip
except ImportError:
    gzip = None
import itertools
import re
import warnings

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
    old_atom : Atom
        The original atom that has been added to the structure already
    new_atom : Atom
        The new atom that we want to see if it is the same as the old atom
    resname : str
        The name of the residue that the new atom would belong to
    resid : int
        The number of the residue that the new atom would belong to
    chain : str
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Structure(object):
    """
    A chemical structure composed of atoms, bonds, angles, torsions, and other
    topological features

    Attributes
    ----------
    atoms : AtomList
        List of all atoms in the structure
    residues : ResidueList
        List of all residues in the structure
    bonds : TrackedList(Bond)
        List of all bonds in the structure
    angles : TrackedList(Angle)
        List of all angles in the structure
    dihedrals : TrackedList(Dihedral)
        List of all dihedrals in the structure -- only one term per dihedral, so
        multi-term dihedral parameters will have the same 4 atoms appear
        multiple times in the list
    urey_bradleys : TrackedList(UreyBradley)
        List of all Urey-Bradley angle bends in the structure
    impropers : TrackedList(Improper)
        List of all CHARMM-style improper torsions in the structure
    cmaps : TrackedList(Cmap)
        List of all CMAP objects in the structure
    trigonal_angles : TrackedList(TrigonalAngle)
        List of all AMOEBA-style trigonal angles in the structure
    out_of_plane_bends : TrackedList(OutOfPlaneBends)
        List of all AMOEBA-style out-of-plane bending angles
    pi_torsions : TrackedList(PiTorsion)
        List of all AMOEBA-style pi-torsion angles
    stretch_bends : TrackedList(StretchBend)
        List of all AMOEBA-style stretch-bend compound bond/angle terms
    torsion_torsions : TrackedList(TorsionTorsion)
        List of all AMOEBA-style coupled torsion-torsion terms
    chiral_frames : TrackedList(ChiralFrame)
        List of all AMOEBA-style chiral frames defined in the structure
    multipole_frames : TrackedList(MultipoleFrame)
        List of all AMOEBA-style multipole frames defined in the structure
    adjusts : TrackedList(NonbondedException)
        List of all AMOEBA-style nonbonded pair-exception rules
    acceptors : TrackedList(AcceptorDonor)
        List of all H-bond acceptors, if that information is present
    donors : TrackedList(AcceptorDonor)
        List of all H-bond donors, if that information is present
    groups : TrackedList(Group)
        List of all CHARMM-style GROUP objects (whatever those are used for)
    box : list of 6 floats
        Box dimensions (a, b, c, alpha, beta, gamma) for the unit cell. If no
        box is defined, `box` is set to `None`
    space_group : str
        The space group of the structure (default is "P 1")

    This class also has a handful of type lists for each of the attributes above
    (excluding `atoms`, `residues`, `chiral_frames`, and `multipole_frames`).
    They are all TrackedList instances that are designed to hold the relevant
    parameter type. The list is:
        bond_types, angle_types, dihedral_types, urey_bradley_types,
        improper_types, cmap_types, trigonal_angle_types,
        out_of_plane_bend_types, pi_torsion_types, stretch_bend_types,
        torsion_torsion_types, adjust_types

    Notes
    -----
    dihedral_types _may_ be a list of DihedralType instances, since torsion
    profiles are often represented by a Fourier series with multiple terms
    """
    #===================================================

    def __init__(self):

        # Topological object lists
        self.atoms = AtomList()
        self.residues = ResidueList()
        self.bonds = TrackedList()
        self.angles = TrackedList()
        self.dihedrals = TrackedList()
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
        self.cmap_types = TrackedList()
        self.trigonal_angle_types = TrackedList()
        self.out_of_plane_bend_types = TrackedList()
        self.pi_torsion_types = TrackedList()
        self.stretch_bend_types = TrackedList()
        self.torsion_torsion_types = TrackedList()
        self.adjust_types = TrackedList()

        self.box = None
        self.space_group = "P 1"

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
        split_dihedrals : bool=False
            If True, then the Dihedral entries will be split up so that each one
            is paired with a single DihedralType (rather than a
            DihedralTypeList)

        Returns
        -------
        cls instance
            The instance of the Structure subclass `cls` with a copy of the
            current Structure's topology information
        """
        c = cls()
        for atom in self.atoms:
            res = atom.residue
            a = copy.copy(atom)
            c.residues.add_atom(a, res.name, res.number, res.chain,
                                res.insertion_code)
            c.atoms.append(a)
        # Now copy all of the types
        for bt in self.bond_types:
            c.bond_types.append(BondType(bt.k, bt.req, c.bond_types))
        for at in self.angle_types:
            c.angle_types.append(AngleType(at.k, at.theteq, c.angle_types))
        if split_dihedrals:
            idx = 0
            for dt in self.dihedral_types:
                dt.starting_index = idx
                if hasattr(dt, '__iter__'):
                    for t in dt:
                        c.dihedral_types.append(
                                DihedralType(t.phi_k, t.per, t.phase,
                                             t.scee, t.scnb,
                                             list=c.dihedral_types)
                        )
                        idx += 1
                else:
                    c.dihedral_types.append(
                            DihedralType(t.phi_k, t.per, t.phase,
                                         t.scee, t.scnb,
                                         list=c.dihedral_types)
                    )
                    idx += 1
        else:
            for dt in self.dihedral_types:
                if hasattr(dt, '__iter__'):
                    dtl = DihedralTypeList()
                    for t in dt:
                        dtl.append(DihedralType(t.phi_k, t.per, t.phase, t.scee,
                                                t.scnb))
                else:
                    dtl = DihedralType(dt.phi_k, dt.per, dt.phase, dt.scee,
                                       dt.scnb)
                dtl.list = c.dihedral_types
                c.dihedral_types.append(dtl)
        for ut in self.urey_bradley_types:
            c.urey_bradley_types.append(BondType(ut.k, ut.req,
                                                 c.urey_bradley_types))
        for it in self.improper_types:
            c.improper_types.append(ImproperType(it.psi_k, it.psi_eq,
                                                 c.improper_types))
        for ct in self.cmap_types:
            c.cmap_types.append(CmapType(ct.resolution, list(ct.grid),
                                         list=c.cmap_types))
        for ta in self.trigonal_angle_types:
            c.trigonal_angle_types.append(
                    AngleType(ta.k, ta.theteq, c.trigonal_angle_types)
            )
        for ot in self.out_of_plane_bend_types:
            c.out_of_plane_bend_types.append(
                    OutOfPlaneBendType(ot.k, self.out_of_plane_bend_types)
            )
        for pt in self.pi_torsion_types:
            c.pi_torsion_types.append(
                    DihedralType(pt.phi_k, pt.per, pt.phase,
                                 list=c.pi_torsion_types)
            )
        for st in self.stretch_bend_types:
            c.stretch_bend_types.append(
                    StretchBendType(st.k, st.req1, st.req2, st.theteq,
                                    c.stretch_bend_types)
            )
        for tt in self.torsion_torsion_types:
            c.torsion_torsion_types.append(
                    TorsionTorsionType(tt.dims, tt.ang1[:], tt.ang2[:],
                                       tt.f.data[:], tt.dfda1.data[:],
                                       tt.dfda2.data[:], tt.d2fda1da2.data[:],
                                       c.torsion_torsion_types)
            )
        for at in self.adjust_types:
            c.adjust_types.append(
                    NonbondedExceptionType(at.vdw_weight, at.multipole_weight,
                                           at.direct_weight, at.polar_weight,
                                           at.mutual_weight, c.adjust_types)
            )
        # Now create the topological objects
        atoms = c.atoms
        for b in self.bonds:
            c.bonds.append(
                    Bond(atoms[b.atom1.idx], atoms[b.atom2.idx],
                         c.bond_types[b.type.idx])
            )
        for a in self.angles:
            c.angles.append(
                    Angle(atoms[a.atom1.idx], atoms[a.atom2.idx],
                          atoms[a.atom3.idx], c.angle_types[a.type.idx])
            )
        if split_dihedrals:
            for d in self.dihedrals:
                if hasattr(d.type, '__iter__'):
                    for i in xrange(len(d.type)):
                        ie = d.ignore_end or i < len(d.type) - 1
                        ti = d.type.starting_index + i
                        c.dihedrals.append(
                                Dihedral(c.atoms[d.atom1.idx],
                                         c.atoms[d.atom2.idx],
                                         c.atoms[d.atom3.idx],
                                         c.atoms[d.atom4.idx],
                                         improper=d.improper, ignore_end=ie,
                                         type=c.dihedral_types[ti])
                        )
                else:
                    c.dihedrals.append(
                        Dihedral(c.atoms[d.atom1.idx], c.atoms[d.atom2.idx],
                                 c.atoms[d.atom3.idx], c.atoms[d.atom4.idx],
                                 improper=d.improper, ignore_end=d.ignore_end,
                                 type=c.dihedral_types[d.type.starting_index])
                    )
        else:
            for d in self.dihedrals:
                c.dihedrals.append(
                        Dihedral(atoms[d.atom1.idx], atoms[d.atom2.idx],
                                 atoms[d.atom3.idx], atoms[d.atom4.idx],
                                 improper=d.improper, ignore_end=d.ignore_end,
                                 type=c.dihedral_types[d.type.idx])
                )
        for u in self.urey_bradleys:
            c.urey_bradleys.append(
                    UreyBradley(atoms[u.atom1.idx], atoms[u.atom2.idx],
                                c.urey_bradley_types[u.type.idx])
            )
        for i in self.impropers:
            c.impropers.append(
                    Improper(atoms[i.atom1.idx], atoms[i.atom2.idx],
                             atoms[i.atom3.idx], atoms[i.atom4.idx],
                             c.improper_types[i.type.idx])
            )
        for cm in self.cmaps:
            c.cmaps.append(
                    Cmap(atoms[cm.atom1.idx], atoms[cm.atom2.idx],
                         atoms[cm.atom3.idx], atoms[cm.atom4.idx],
                         atoms[cm.atom5.idx], c.cmap_types[cm.type.idx])
            )
        for t in self.trigonal_angles:
            c.trigonal_angles.append(
                    TrigonalAngle(atoms[t.atom1.idx], atoms[t.atom2.idx],
                                  atoms[t.atom3.idx], atoms[t.atom4.idx],
                                  c.trigonal_angle_types[t.type.idx])
            )
        for o in self.out_of_plane_bends:
            c.out_of_plane_bends.append(
                    OutOfPlaneBend(atoms[t.atom1.idx], atoms[t.atom2.idx],
                                   atoms[t.atom3.idx], atoms[t.atom4.idx],
                                   c.out_of_plane_bend_types[o.type.idx])
            )
        for p in self.pi_torsions:
            c.pi_torsions.append(
                    PiTorsion(atoms[p.atom1.idx], atoms[p.atom2.idx],
                              atoms[p.atom3.idx], atoms[p.atom4.idx],
                              atoms[p.atom5.idx], atoms[p.atom6.idx],
                              c.pi_torsion_types[p.type.idx])
            )
        for s in self.stretch_bends:
            c.stretch_bends.append(
                    StretchBend(atoms[s.atom1.idx], atoms[s.atom2.idx],
                                atoms[s.atom3.idx],
                                c.stretch_bend_types[s.type.idx])
            )
        for t in self.torsion_torsions:
            c.torsion_torsions.append(
                    TorsionTorsion(atoms[t.atom1.idx], atoms[t.atom2.idx],
                                   atoms[t.atom3.idx], atoms[t.atom4.idx],
                                   atoms[t.atom5.idx],
                                   c.torsion_torsion_types[t.type.idx])
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
                                       c.adjust_types[a.type.idx])
            )
        for a in self.acceptors:
            c.acceptors.append(
                    AcceptorDonor(atoms[a.atom1.idx], atoms[a.atom2.idx])
            )
        for d in self.donors:
            c.acceptors.append(
                    AcceptorDonor(atoms[d.atom1.idx], atoms[d.atom2.idx])
            )
        for g in self.groups:
            c.groups.append(Group(g.bs, g.type, g.move))
        c.box = copy.copy(self.box)
        return c

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
                self.pi_torsion_types.changed)

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
        for dihedral in self.dihedrals:
            if not dihedral.ignore_end: continue
            if (dihedral.atom1 in dihedral.atom4.bond_partners or
                dihedral.atom1 in dihedral.atom4.angle_partners):
                dihedral.ignore_end = True

    #===================================================

    def write_pdb(self, dest, renumber=True, coordinates=None,
                  altlocs='all'):
        """
        Write a PDB file from the current Structure instance

        Parameters
        ----------
        dest : str or file-like
            Either a file name or a file-like object containing a `write`
            method to which to write the PDB file. If it is a filename that
            ends with .gz or .bz2, a compressed version will be written using
            either gzip or bzip2, respectively.
        renumber : bool=True
            If True, renumber the atoms and residues sequentially as they are
            stored in the structure.  If False, use the original numbering if
            it was assigned previously
        coordinates : array-like of float=None
            If provided, these coordinates will be written to the PDB file
            instead of the coordinates stored in the structure. These
            coordinates should line up with the atom order in the structure
            (not necessarily the order of the "original" PDB file if they
            differ)
        altlocs : str='all'
            Keyword controlling which alternate locations are printed to the
            resulting PDB file. Allowable options are:
                - 'all' : (default) print all alternate locations
                - 'first' : print only the first alternate locations
                - 'occupancy' : print the one with the largest occupancy. If two
                  conformers have the same occupancy, the first one to occur is
                  printed
            Input is case-insensitive, and partial strings are permitted as long
            as it is a substring of one of the above options that uniquely
            identifies the choice.
        """
        if altlocs.lower() == 'all'[:len(altlocs)]:
            altlocs = 'all'
        elif altlocs.lower() == 'first'[:len(altlocs)]:
            altlocs = 'first'
        elif altlocs.lower() == 'occupancy'[:len(altlocs)]:
            altlocs = 'occupancy'
        else:
            raise ValueError("Illegal value of occupancy [%s]; expected 'all', "
                             "'first', or 'occupancy'" % altlocs)
        own_handle = False
        if not hasattr(dest, 'write'):
            if dest.endswith('.gz'):
                dest = gz.open(dest, 'w')
            elif dest.endswith('.bz2'):
                dest = bz2.BZ2File(dest, 'w')
            else:
                dest = open(dest, 'w')
            own_handle = True
        atomrec = ('ATOM  %5d %-4s%-1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f'
                   '%6.2f          %2s%-2s\n')
        terrec = ('TER   %5d      %-3s %1s%4d\n')
        nchains = len(set([res.chain for res in self.residues if res.chain]))
        if self.box is not None:
            a, b, c, alpha, beta, gamma = self.box
            dest.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n' % (
                       self.box[0], self.box[1], self.box[2], self.box[3],
                       self.box[4], self.box[5], self.space_group, nchains))
        if coordinates is not None:
            try:
                crdsize = len(coordinates)
            except TypeError:
                raise TypeError("Cannot find length of coordinates")
            if crdsize == len(self.atoms):
                try:
                    coords = coordinates.flatten()
                except AttributeError:
                    try:
                        coords = list(itertools.chain(*coordinates))
                    except TypeError:
                        raise TypeError("Unsupported coordinate dimensionality")
                if len(coords) != len(self.atoms) * 3:
                    raise TypeError("Unsupported coordinate shape")
            elif crdsize == len(self.atoms) * 3:
                coords = coordinates
            else:
                raise TypeError("Coordinates has unexpected shape")
        else:
            coords = [[a.xx, a.xy, a.xz] for a in self.atoms]
            coords = list(itertools.chain(*coords))
        # Create a function to process each atom and return which one we want
        # to print, based on our alternate location choice
        if altlocs == 'all':
            def print_atoms(atom, coords):
                i3 = atom.idx * 3
                return atom, atom.other_locations, coords[i3:i3+3]
        elif altlocs == 'first':
            def print_atoms(atom, coords):
                i3 = atom.idx * 3
                return atom, dict(), coords[i3:i3+3]
        elif altlocs == 'occupancy':
            def print_atoms(atom, coords):
                occ = atom.occupancy
                a = atom
                for key, item in atom.other_locations.iteritems():
                    if item.occupancy > occ:
                        occ = item.occupancy
                        a = item
                return a, dict(), [a.xx, a.xy, a.xz]
        else:
            raise Exception("Should not be here!")
        if renumber:
            nmore = 0 # how many *extra* atoms have been added?
            for res in self.residues:
                for atom in res.atoms:
                    i3 = atom.idx * 3
                    anum = (atom.idx + 1 + nmore) % 100000
                    rnum = (res.idx + 1) % 10000
                    pa, others, (x, y, z) = print_atoms(atom, coords)
                    dest.write(atomrec % (anum , pa.name, pa.altloc,
                               res.name, res.chain, rnum, res.insertion_code,
                               x, y, z, pa.occupancy, pa.bfactor,
                               Element[atom.atomic_number].upper(), ''))
                    for key in sorted(others.keys()):
                        oatom = others[key]
                        nmore += 1
                        anum = (pa.idx + 1 + nmore) % 100000
                        x, y, z = oatom.xx, oatom.xy, oatom.xz
                        dest.write(atomrec % (anum, oatom.name, key, res.name,
                                   res.chain, rnum, res.insertion_code, x, y,
                                   z, oatom.occupancy, oatom.bfactor,
                                   Element[oatom.atomic_number].upper(), ''))
                if res.ter:
                    dest.write(terrec % (anum+1, res.name, res.chain, rnum))
                    nmore += 1
        else:
            def acmp(x, y):
                xn = x.number or x.idx + 1
                yn = y.number or y.idx + 1
                return xn - yn
            last_number = 0
            last_rnumber = 0
            for res in self.residues:
                for atom in sorted(res.atoms, cmp=acmp):
                    i3 = atom.idx * 3
                    rnum = atom.residue.number or last_rnumber + 1
                    pa, others, (x, y, z) = print_atoms(atom, coords)
                    num = pa.number or last_number + 1
                    dest.write(atomrec % (num % 100000, pa.name, pa.altloc,
                               res.name, res.chain, rnum % 10000,
                               res.insertion_code, x, y, z,
                               pa.occupancy, pa.bfactor,
                               Element[atom.atomic_number].upper(), ''))
                    last_number = num
                    for key in sorted(others.keys()):
                        oatom = others[key]
                        anum = oatom.number or last_number + 1
                        x, y, z = oatom.xx, oatom.xy, oatom.xz
                        dest.write(atomrec % (anum % 100000, oatom.name, key,
                                   res.name, res.chain, rnum,
                                   res.insertion_code, x, y, z,
                                   oatom.occupancy, oatom.bfactor,
                                   Element[oatom.atomic_number].upper(), ''))
                        last_number = anum
                if res.ter:
                    dest.write(terrec % (last_number+1, res.name, res.chain,
                               rnum))
                    last_number += 1
                last_rnumber = rnum
        if own_handle:
            dest.close()

    #===================================================

    def _prune_empty_bonds(self):
        """ Gets rid of any empty bonds """
        for i in reversed(xrange(len(self.bonds))):
            bond = self.bonds[i]
            if bond.atom1 is None and bond.atom2 is None:
                del self.bonds[i]
            elif bond.atom1.idx == -1 or bond.atom2.idx == -1:
                bond.delete()
                del self.bonds[i]

    #===================================================

    def _prune_empty_angles(self):
        """ Gets rid of any empty angles """
        for i in reversed(xrange(len(self.angles))):
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
        for i in reversed(xrange(len(self.dihedrals))):
            dihed = self.dihedrals[i]
            if (dihed.atom1 is None and dihed.atom2 is None and
                    dihed.atom3 is None and dihed.atom4 is None):
                del self.dihedrals[i]
            elif (dihed.atom1.idx == -1 or dihed.atom2.idx == -1 or
                    dihed.atom3.idx == -1 or dihed.atom4.idx == -1):
                dihed.delete()
                del self.dihedrals[i]

    #===================================================

    def _prune_empty_ureys(self):
        """ Gets rid of any empty Urey-Bradley terms """
        for i in reversed(xrange(len(self.urey_bradleys))):
            ub = self.urey_bradleys[i]
            if ub.atom1 is None and ub.atom2 is None:
                del self.urey_bradleys[i]
            elif ub.atom1.idx == -1 or ub.atom2.idx == -1:
                ub.delete()
                del self.urey_bradleys[i]

    #===================================================

    def _prune_empty_impropers(self):
        """ Gets rid of any empty improper torsions """
        for i in reversed(xrange(len(self.impropers))):
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
        for i in reversed(xrange(len(self.cmaps))):
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
        for i in reversed(xrange(len(self.trigonal_angles))):
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
        for i in reversed(xrange(len(self.out_of_plane_bends))):
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
        for i in reversed(xrange(len(self.pi_torsions))):
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
        for i in reversed(xrange(len(self.stretch_bends))):
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
        for i in reversed(xrange(len(self.torsion_torsions))):
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
        for i in reversed(xrange(len(self.chiral_frames))):
            cf = self.chiral_frames[i]
            if cf.atom1 is None or cf.atom2 is None:
                del self.chiral_frames[i]
            elif cf.atom1.idx == -1 or cf.atom2.idx == -1:
                del self.chiral_frames[i]

    #===================================================

    def _prune_empty_multipole_frames(self):
        """ Gets rid of any empty multipole frame terms """
        for i in reversed(xrange(len(self.multipole_frames))):
            mf = self.multipole_frames[i]
            if mf.atom is None or mf.atom.idx == -1:
                del self.multipole_frames[i]

    #===================================================

    def _prune_empty_adjusts(self):
        """ Gets rid of any empty nonbonded exception adjustments """
        for i in reversed(xrange(len(self.adjusts))):
            adj = self.adjusts[i]
            if adj.atom1 is None or adj.atom2 is None:
                del self.adjusts[i]
            elif adj.atom1.idx == -1 or adj.atom2.idx == -1:
                del self.adjusts[i]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def read_PDB(filename):
    """
    Read a PDB file and return a populated `Structure` class

    Parameters
    ----------
    filename : str or file-like
        Name of PDB file to read, or a file-like object that can iterate over
        the lines of a PDB. Compressed file names can be specified and are
        determined by file-name extension (e.g., file.pdb.gz, file.pdb.bz2)

    Metadata
    --------
    The PDB parser also adds metadata to the returned Structure object that may
    be present in the PDB file

    experimental : str
        EXPDTA record
    journal : str
        JRNL record
    authors : str
        AUTHOR records
    keywords : str
        KEYWDS records
    doi : str
        DOI from the JRNL record
    pmid : str
        PMID from the JRNL record
    journal_authors : str
        Author info from the JRNL record
    volume : str
        Volume of the published article from the JRNL record
    page : str
        Page of the published article from the JRNL record
    title : str
        TITL section of the JRNL record
    year : int=None
        Year that the article was published, from the JRNL record
    related_entries : list of (str, str)
        List of entries in other databases 

    Returns
    -------
    structure
    
    structure : Structure
        The Structure object initialized with all of the information from the
        PDB file.  No bonds or other topological features are added by default.

    Notes
    -----
    The returned structure has an extra attribute, pdbxyz, that contains all of
    the coordinates for all of the frames in the PDB file as a list of NATOM*3
    lists.
    """
    global relatere
    if isinstance(filename, basestring):
        own_handle = True
        if filename.endswith('.gz'):
            if gzip is None:
                raise ImportError('gzip is not available for compressed PDB')
            fileobj = gzip.open(filename, 'r')
        elif filename.endswith('.bz2'):
            if bz2 is None:
                raise ImportError('bz2 is not available for compressed PDB')
            fileobj = bz2.BZ2File(filename, 'r')
        else:
            fileobj = open(filename, 'r')
    else:
        own_handle = False
        fileobj = filename

    struct = Structure()
    # Add metadata fields
    struct.experimental = struct.journal = struct.authors = struct.keywords = ''
    struct.doi = struct.pmid = struct.journal_authors = struct.volume_page = ''
    struct.title = ''
    struct.year = None
    struct.related_entries = []
    modelno = 1 # For PDB files with multiple MODELs
    atomno = 0
    coordinates = []
    all_coordinates = []

    # Support hexadecimal numbering like that printed by VMD
    last_atom = Atom()
    last_resid = 1
    resend = 26
    res_hex = False
    atom_hex = False
    ZEROSET = set('0')

    try:
        for line in fileobj:
            try:
                line = line.decode('ascii')
            except AttributeError:
                # Assume this is a string in Py3 which doesn't have 'decode'
                pass
            rec = line[:6]
            if rec == 'ATOM  ' or rec == 'HETATM':
                atomno += 1
                atnum, atname, altloc = line[6:11], line[12:16], line[16]
                resname, chain, resid = line[17:20], line[21], line[22:resend]
                inscode = line[26]
                x, y, z = line[30:38], line[38:46], line[47:54]
                occupancy, bfactor = line[54:60], line[60:66]
                elem, chg = line[76:78], line[78:80]
                elem = '%-2s' % elem # Make sure we have at least 2 characters
                if elem[0] == ' ': elem = elem[1] + ' '
                try:
                    atsym = (elem[0] + elem[1].lower()).strip()
                    atomic_number = AtomicNum[atsym]
                    mass = Mass[atsym]
                except KeyError:
                    # Now try based on the atom name... but don't try too hard
                    # (e.g., don't try to differentiate b/w Ca and C)
                    try:
                        atomic_number = AtomicNum[atname.strip()[0].upper()]
                        mass = Mass[atname.strip()[0].upper()]
                    except KeyError:
                        try:
                            sym = atname.strip()[:2]
                            sym = '%s%s' % (sym[0].upper(), sym[0].lower())
                            atomic_number = AtomicNum[sym]
                            mass = Mass[sym]
                        except KeyError:
                            atomic_number = 0 # give up
                            mass = 0.0
                try:
                    bfactor = float(bfactor)
                except ValueError:
                    bfactor = 0.0
                try:
                    occupancy = float(occupancy)
                except ValueError:
                    occupancy = 0.0
                # Figure out what my residue number is and see if the PDB is
                # outputting residue numbers in hexadecimal (e.g., VMD)
                if last_resid >= 9999 and resend == 26:
                    if not res_hex and resid == '9999':
                        resid = 9999
                    elif not res_hex:
                        res_hex = int(resid, 16) == 10000
                    # So now we know if we use hexadecimal or not. If we do,
                    # convert. Otherwise, stay put
                    if res_hex:
                        try:
                            resid = int(resid, 16)
                        except ValueError, e:
                            if resid == '****':
                                resid = None # Figure out by unique atoms
                            else:
                                raise e
                    elif resid == '1000' and line[26] == '0':
                        resend += 1
                        resid = 10000
                    else:
                        resid = int(resid)
                elif resend > 26:
                    # VMD extends the field now... ugh.
                    if resid[0] == '1' and set(resid[1:]) == ZEROSET:
                        if line[resend] == '0':
                            resid = int(resid) * 10
                            resend += 1
                        else:
                            resid = int(resid)
                    else:
                        resid = int(resid)
                else:
                    resid = int(resid)
                # If the number has cycled, it too may be hexadecimal
                if atom_hex:
                    atnum = int(atnum, 16)
                else:
                    try:
                        atnum = int(atnum)
                    except ValueError:
                        atnum = int(atnum, 16)
                        atom_hex = True
                # It's possible that the residue number has cycled so much that
                # it is now filled with ****'s. In that case, start a new
                # residue if the current residue repeats the same atom name as
                # the 'last' residue. Do not worry about atom numbers going to
                # *****'s, since that is >1M atoms.
                if resid is None:
                    for atom in struct.residues[-1]:
                        if atom.name == atname:
                            resid = last_resid + 1
                            break
                if resid is None:
                    # Still part of the last residue
                    resid = last_resid
                last_resid = resid
                try:
                    chg = float(chg)
                except ValueError:
                    chg = 0
                atom = Atom(atomic_number=atomic_number, name=atname,
                            charge=chg, mass=mass, occupancy=occupancy,
                            bfactor=bfactor, altloc=altloc, number=atnum)
                atom.xx, atom.xy, atom.xz = float(x), float(y), float(z)
                if _compare_atoms(last_atom, atom, resname, resid, chain):
                    atom.residue = last_atom.residue
                    last_atom.other_locations[altloc] = atom
                    continue
                last_atom = atom
                if modelno == 1:
                    struct.residues.add_atom(atom, resname, resid,
                                             chain, inscode)
                    struct.atoms.append(atom)
                else:
                    try:
                        orig_atom = struct.atoms[atomno-1]
                    except IndexError:
                        raise PDBError('Atom %d differs in MODEL %d [%s %s vs. '
                                       '%s %s]' % (atomno, modelno,
                                       atom.residue.name, atom.name, resname,
                                       atname))
                    if (orig_atom.residue.name != resname.strip()
                            or orig_atom.name != atname.strip()):
                        raise PDBError('Atom %d differs in MODEL %d [%s %s vs. '
                                       '%s %s]' % (atomno, modelno,
                                       orig_atom.residue.name, orig_atom.name,
                                       resname, atname))
                coordinates.extend([atom.xx, atom.xy, atom.xz])
            elif rec.strip() == 'TER':
                if modelno == 1: last_atom.residue.ter = True
            elif rec == 'ENDMDL':
                # End the current model
                if len(struct.atoms) == 0:
                    raise PDBError('MODEL ended before any atoms read in')
                modelno += 1
                if len(struct.atoms)*3 != len(coordinates):
                    raise ValueError(
                            'Inconsistent atom numbers in some PDB models')
                all_coordinates.append(coordinates)
                atomno = 0
                coordinates = []
                resend = 26
            elif rec == 'MODEL ':
                if modelno == 1 and len(struct.atoms) == 0: continue
                if len(coordinates) > 0:
                    if len(struct.atoms)*3 != len(coordinates):
                        raise ValueError(
                                'Inconsistent atom numbers in some PDB models')
                    warnings.warn('MODEL not explicitly ended', PDBWarning)
                    all_coordinates.append(coordinates)
                    coordinates = []
                modelno += 1
                atomno = 0
                resend = 26
            elif rec == 'CRYST1':
                a = float(line[6:15])
                b = float(line[15:24])
                c = float(line[24:33])
                try:
                    A = float(line[33:40])
                    B = float(line[40:47])
                    C = float(line[47:54])
                except (IndexError, ValueError):
                    A = B = C = 90.0
                struct.box = [a, b, c, A, B, C]
                try:
                    struct.space_group = line[55:66].strip()
                except IndexError:
                    pass
            elif rec == 'EXPDTA':
                struct.experimental = line[6:].strip()
            elif rec == 'AUTHOR':
                struct.authors += line[10:].strip()
            elif rec == 'JRNL  ':
                part = line[12:16]
                if part == 'AUTH':
                    struct.journal_authors += line[19:].strip()
                elif part == 'TITL':
                    struct.title += ' %s' % line[19:].strip()
                elif part == 'REF ':
                    struct.journal += ' %s' % line[19:47].strip()
                    if not line[16:18].strip():
                        struct.volume = line[51:55].strip()
                        struct.page = line[56:61].strip()
                        try:
                            struct.year = int(line[62:66])
                        except ValueError:
                            pass
                elif part == 'PMID':
                    struct.pmid = line[19:].strip()
                elif part == 'DOI ':
                    struct.doi = line[19:].strip()
            elif rec == 'KEYWDS':
                struct.keywords += '%s,' % line[10:]
            elif rec == 'REMARK' and line[6:10] == ' 900':
                # Related entries
                rematch = relatere.match(line[11:])
                if rematch:
                    struct.related_entries.append(rematch.groups())
    finally:
        # Make sure our file is closed if we opened it
        if own_handle: fileobj.close()

    # Post-process some of the metadata to make it more reader-friendly
    struct.keywords = [s.strip() for s in struct.keywords.split(',')
                                        if s.strip()]
    struct.journal = struct.journal.strip()
    struct.title = struct.title.strip()

    struct.unchange()
    if coordinates:
        if len(coordinates) != 3*len(struct.atoms):
            raise ValueError('bad number of atoms in some PDB models')
        all_coordinates.append(coordinates)
    struct.pdbxyz = all_coordinates
    return struct

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def write_PDB(struct, dest, renumber=True, coordinates=None, altlocs='all'):
    """
    Write a PDB file from a structure instance

    Parameters
    ----------
    struct : Structure
        A Structure instance from which to write the PDB file
    dest : str or file-like
        Either a file name or a file-like object containing a `write`
        method to which to write the PDB file
    renumber : bool=True
        If True, renumber the atoms and residues sequentially as they are
        stored in the structure.  If False, use the original numbering if
        it was assigned previously
    coordinates : array-like of float=None
        If provided, these coordinates will be written to the PDB file
        instead of the coordinates stored in the structure. These
        coordinates should line up with the atom order in the structure
        (not necessarily the order of the "original" PDB file if they
        differ)
    altlocs : str='all'
        Keyword controlling which alternate locations are printed to the
        resulting PDB file. Allowable options are:
            - 'all' : (default) print all alternate locations
            - 'first' : print only the first alternate locations
            - 'occupancy' : print the one with the largest occupancy. If two
              conformers have the same occupancy, the first one to occur is
              printed
        Input is case-insensitive, and partial strings are permitted as long
        as it is a substring of one of the above options that uniquely
        identifies the choice.
    """
    struct.write_pdb(dest, renumber, coordinates, altlocs)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
