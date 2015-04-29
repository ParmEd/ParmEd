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


def _n_prior(pose, nbr):
    prior = -1
    for i in xrange(1, nbr.rsd()):
        prior += pose.residue(i).natoms()
    return prior + nbr.atomno()


class RosettaPose(object):

    @staticmethod
    def load(pose):
        """ Load a :class:`Pose` object and return a populated :class:`Structure` class

            Parameters
            ----------
            pose : :class:`Pose`
                PyRosetta :class:`Pose` object to convert
        """
        if not Pose or not AtomID:
            raise ImportError('Could not load the PyRosetta module.')
        if not isinstance(pose, Pose):
            raise TypeError('Object is not a PyRosetta Pose object.')

        struct = Structure()

        atnum = 1
        conf = pose.conformation()
        for resid in xrange(1, pose.total_residue()+1):
            res = pose.residue(resid)
            resname = res.name3().strip()
            chain = chr(res.chain()+ord('A')-1)
            for atno, at in enumerate(res.atoms(), start=1):
                try:
                    atinfo = res.atom_type(atno)
                    atname = res.atom_name(atno).strip()
                    if atinfo.is_virtual():
                        atsym = 'EP'
                    else:
                        atsym = atinfo.element()
                    rmin = atinfo.lj_radius()
                    epsilon = atinfo.lj_wdepth()
                    atomic_number = AtomicNum[atsym]
                    mass = Mass[atsym]
                except KeyError:
                    raise RosettaError('Could not recognize element: %s.'
                                       % atsym)

                params = dict(atomic_number=atomic_number, name=atname,
                              charge=0.0, mass=mass, occupancy=0.0,
                              bfactor=0.0, altloc='', number=atnum,
                              rmin=rmin, epsilon=epsilon)

                if atinfo.is_virtual():
                    atom = ExtraPoint(**params)
                else:
                    atom = Atom(**params)

                atom.xx, atom.xy, atom.xz = tuple(at.xyz())

                struct.add_atom(atom, resname, resid, chain, '')
                atnum += 1
                try:
                    for nbr in conf.bonded_neighbor_all_res(AtomID(atno,
                                                                   resid)):
                        if nbr.rsd() < resid or (nbr.rsd() == resid
                                                 and nbr.atomno() < atno):
                            struct.bonds.append(
                                Bond(struct.atoms[_n_prior(pose, nbr)],
                                     atom))
                except:
                    raise RosettaError('Could not add bonds.')

        struct.unchange()
        return struct
