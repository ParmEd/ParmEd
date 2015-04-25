"""
This package contains classes responsible for loading and dumping PyRosetta
Pose objects.
"""

from __future__ import print_function

from chemistry.exceptions import RosettaError, RosettaWarning
from chemistry.periodic_table import AtomicNum, Mass
from chemistry.structure import Structure
from chemistry.topologyobjects import Atom, Bond

try:
    from rosetta import Pose, AtomID
except ImportError:
    RosettaWarning('Cannot find the PyRosetta module.')


class RosettaPose(object):

    @staticmethod
    def load(pose):
        """ Load a Pose object and return a populated `Structure` class
        """

        if not isinstance(pose, Pose):
            raise RosettaError('Object is not a PyRosetta Pose object.')

        struct = Structure()
        struct.experimental = struct.journal = struct.authors = ''
        struct.keywords = struct.doi = struct.pmid = ''
        struct.journal_authors = struct.volume_page = ''
        struct.title = ''
        struct.year = None
        struct.related_entries = []

        # try:
        atnum = 1
        conf = pose.conformation()
        for resid in xrange(1, pose.total_residue()+1):
            res = pose.residue(resid)
            resname = res.name3().strip()
            chain = chr(res.chain()+ord('A')-1)
            for atno, at in enumerate(res.atoms(), start=1):
                try:
                    atname = res.atom_name(atno).strip()
                    atsym = res.atom_type(atno).element()
                    rmin = res.atom_type(atno).lj_radius()
                    epsilon = res.atom_type(atno).lj_wdepth()
                    atomic_number = AtomicNum[atsym]
                    mass = Mass[atsym]
                except KeyError:
                    raise RosettaError('')

                atom = Atom(atomic_number=atomic_number, name=atname,
                            charge=0.0, mass=mass, occupancy=0.0,
                            bfactor=0.0, altloc='', number=atnum,
                            rmin=rmin, epsilon=epsilon)
                atom.xx, atom.xy, atom.xz = tuple(at.xyz())
                struct.add_atom(atom, resname, resid, chain, '')
                for nbr in conf.bonded_neighbor_all_res(AtomID(atno,
                                                               resid)):
                    if nbr.rsd() <= resid and nbr.atomno() < atno:
                        Bond(struct.atoms[sum([pose.residue(i).natoms()
                                          for i in xrange(1, nbr.rsd())])
                                          + nbr.atomno() - 1],
                             atom)
                atnum += 1

        # except:
            # raise RosettaError('')

        struct.unchange()
        return struct
