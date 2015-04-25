"""
This package contains classes responsible for loading and dumping PyRosetta
Pose objects.
"""

from __future__ import print_function

from chemistry.exceptions import RosettaError
from chemistry.periodic_table import AtomicNum, Mass
from chemistry.structure import Structure
from chemistry.topologyobjects import Atom, ExtraPoint, Bond

try:
    from rosetta import Pose, AtomID
except ImportError:
    Pose = AtomID = None
    pass


class RosettaPose(object):

    @staticmethod
    def load(self, pose):
        """ Load a Pose object and return a populated `Structure` class
        """
        if not Pose or not AtomID:
            raise ImportError('Could not load the PyRosetta module.')
        if not isinstance(pose, Pose):
            raise TypeError('Object is not a PyRosetta Pose object.')

        self.pose = pose
        struct = Structure()

        try:
            atnum = 1
            conf = pose.conformation()
            for resid in xrange(1, pose.total_residue()+1):
                res = pose.residue(resid)
                resname = res.name3().strip()
                chain = chr(res.chain()+ord('A')-1)
                for atno, at in enumerate(res.atoms(), start=1):
                    try:
                        atinfo = res.atom_type(atno)
                        if atinfo.is_virtual():
                            atname = 'EP'
                        else:
                            atname = res.atom_name(atno).strip()
                        atsym = atinfo.element()
                        rmin = atinfo.lj_radius()
                        epsilon = atinfo.lj_wdepth()
                        atomic_number = AtomicNum[atsym]
                        mass = Mass[atsym]
                    except KeyError:
                        raise RosettaError('Could not recognize element: %s.'
                                           % atsym)
                    if atinfo.is_virtual():
                        atom = ExtraPoint(atomic_number=atomic_number,
                                          name=atname, charge=0.0, mass=mass,
                                          occupancy=0.0, bfactor=0.0,
                                          altloc='', number=atnum, rmin=rmin,
                                          epsilon=epsilon)
                    else:
                        atom = Atom(atomic_number=atomic_number, name=atname,
                                    charge=0.0, mass=mass, occupancy=0.0,
                                    bfactor=0.0, altloc='', number=atnum,
                                    rmin=rmin, epsilon=epsilon)
                    atom.xx, atom.xy, atom.xz = tuple(at.xyz())
                    struct.add_atom(atom, resname, resid, chain, '')
                    try:
                        for nbr in conf.bonded_neighbor_all_res(AtomID(atno,
                                                                       resid)):
                            if nbr.rsd() <= resid and nbr.atomno() < atno:
                                struct.bonds.append(
                                    Bond(struct.atoms[self._n_prior_(nbr)],
                                         atom))
                        atnum += 1
                    except:
                        raise RosettaError('Could not add bonds.')

        except:
            raise RosettaError('Could not load structure.')

        struct.unchange()
        return struct

        def _n_prior_(self, nbr):
            prior = 0
            for i in xrange(1, nbr.rsd()):
                prior += self.pose.residue(i).natoms()
            return prior + nbr.atomno() - 1
