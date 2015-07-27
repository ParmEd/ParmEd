"""
Contains all of the class objects for the Tinker topology
"""
from __future__ import division
from parmed.utils.six.moves import range, zip
from collections import OrderedDict

class Atom(object):
    """ An atom in the system """

    def __init__(self, symbol, type_, class_, 
                 atomic_number, mass, valence, desc):
        self.symbol = str(symbol).strip()
        self.type = int(type_)
        self.class_ = int(class_)
        self.atomic_number = int(atomic_number)
        self.mass = float(mass)
        self.valence = int(valence)
        self.desc = str(desc).strip()

    def add_vdw(self, radius, epsilon, radius14, epsilon14, reduction):
        """ Adds van der Waals terms to the atom """
        self.radius = float(radius)
        self.epsilon = float(epsilon)
        # These values might not exist
        if not str(radius14).strip():
            self.radius14 = None
        else:
            self.radius14 = float(radius14)
        if not str(epsilon14).strip():
            self.epsilon14 = None
        else:
            self.epsilon14 = float(epsilon14)
        if not str(reduction).strip():
            self.reduction = None
        else:
            self.reduction = float(reduction)

    def __repr__(self):
        return '<Atom %s [#%d]; Mass=%.2f; "%s">' % (
                        self.symbol, self.atomic_number, self.mass, self.desc)

    def __str__(self):
        retstr =  ('Atom %s [#%s] %s\n'
                    '   Type:        %d\n'
                    '   Class:       %d\n'
                    '   Mass:        %.2f\n'
                    '   Valence:     %d' % 
                    (self.symbol, self.atomic_number, self.desc, self.type,
                     self.class_, self.mass, self.valence)
        )

        if hasattr(self, 'radius'):
            retstr += ('\n   vdW Rad:     %.2f\n'
                        '   vdW Eps:     %.2f' % (self.radius, self.epsilon))
            if self.radius14 is not None or self.epsilon14 is not None:
                retstr += ('\n   vdW 1-4 Rad: %.2f\n   vdW 1-4 Eps: %.2f' %
                           (self.radius14, self.epsilon14)
                )

        if self.reduction is not None:
            retstr += '\n   vdW Reduct:  %.2f' % (self.reduction)

        return retstr

class _ParamTypeList(list):
    """
    A list for all types of parameters. This should be subclassed for each
    parameter type that's necessary. It performs type checking as well as other
    tasks.
    """

    typeclass = None

    def __init__(self):
        super(_ParamTypeList, self).__init__()

    def add(self, *args, **kwargs):
        self.append(self.typeclass(*args, **kwargs))

    def append(self, thing):
        if not isinstance(thing, self.typeclass):
            raise TypeError('Can only append "%s" objects to %s' %
                            (self.typeclass.__name__, type(self).__name__))
        super(_ParamTypeList, self).append(thing)

    def extend(self, things):
        for thing in things:
            self.append(thing)

class AtomList(_ParamTypeList):
    """ List for Atom instances """
    typeclass = Atom

class BondStretch(object):
    """ A bond-stretching term """
    def __init__(self, atom1, atom2, k, req):
        self.atom1 = atom1
        self.atom2 = atom2
        self.k = float(k)
        self.req = float(req)

    def __repr__(self):
        return '<%s [%s-%s]; k=%.2f; Req=%.2f>' % (type(self).__name__,
                self.atom1.symbol, self.atom2.symbol, self.k, self.req)

    def __str__(self):
        return ('%s %r --- %r\n'
                '     K = %.4f\n'
                '   Req = %.4f' % (type(self).__name__, self.atom1, self.atom2,
                                   self.k, self.req)
        )

class BondList(_ParamTypeList):
    typeclass = BondStretch

class AngleBend(object):
    """ An angle-bending term """
    def __init__(self, atom1, atom2, atom3, k, theteq, fold, type_):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.k = float(k)
        self.theteq = float(theteq)
        if str(fold).strip():
            self.fold = float(fold)
        else:
            self.fold = None
        self.type = type_.strip() or 'Harmonic'

    def __repr__(self):
        retstr = '<AngleBend [%s-%s-%s]; k=%.2f; Theta_eq=%.2f' % (
                self.atom1.symbol, self.atom2.symbol, self.atom3.symbol,
                self.k, self.theteq
        )
        if self.fold is not None:
            retstr += '; fold=%.1f' % self.fold
        return retstr + '; %s>' % self.type

    def __str__(self):
        retstr = ('AngleBend %r --- %r --- %r\n'
                  '     K = %.4f\n'
                  ' Theta = %.4f\n' % (self.atom1, self.atom2, self.atom3,
                                       self.k, self.theteq)
        )
        if self.fold is not None:
            retstr += '  Fold = %.1f\n' % (self.fold)
        return retstr + '  Type = %s' % (self.type)

class AngleList(_ParamTypeList):
    typeclass = AngleBend

class StretchBend(object):
    """ A stretch-bending term """
    def __init__(self, atom1, atom2, atom3, k, theteq, r1eq, r2eq):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.k = float(k)
        self.theteq = float(theteq)
        self.r1eq = float(r1eq)
        self.r2eq = float(r2eq)

    def __repr__(self):
        return ('<StretchBend [%s-%s-%s]; k=%.2f; Theta_eq=%.2f; R1eq=%.2f; '
                'R2eq=%.2f>' % (self.atom1.symbol, self.atom2.symbol,
                                self.atom3.symbol, self.k, self.theteq,
                                self.r1eq, self.r2eq)
        )

    def __str__(self):
        return ('StretchBend %r --- %r --- %r\n'
                '     K = %.4f\n'
                ' Theta = %.4f\n'
                '  R1eq = %.4f\n'
                '  R2eq = %.4f' % (self.atom1, self.atom2, self.atom3, self.k,
                                   self.theteq, self.r1eq, self.r2eq))

class StretchBendList(_ParamTypeList):
    typeclass = StretchBend

class UreyBradley(BondStretch):
    """ A urey-bradley term -- functionally identical to a bond-stretch """

class UreyBradleyList(_ParamTypeList):
    typeclass = UreyBradley

class OutOfPlaneBend(object):
    """ An out-of-plane bending term """
    def __init__(self, atom1, atom2, atom3, atom4, k):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.k = float(k)

    def __repr__(self):
        return '<%s [%s-%s-%s-%s]; k=%.2f>' % (
                type(self).__name__, self.atom1.symbol, self.atom2.symbol,
                self.atom3.symbol, self.atom4.symbol, self.k
        )

    def __str__(self):
        return ('%s %r --- %r --- %r --- %r\n'
                '   K = %.2f' % (type(self).__name__, self.atom1, self.atom2,
                                 self.atom3, self.atom4, self.k)
        )

class OutOfPlaneBendList(_ParamTypeList):
    typeclass = OutOfPlaneBend

class OutOfPlaneDist(OutOfPlaneBend):
    """ An out-of-plane distance (functionally equivalent to OOP Bend """

class OutOfPlaneDistList(_ParamTypeList):
    typeclass = OutOfPlaneDist

class TorsionAngle(object):
    """ Torsional Angle parameter """
    def __init__(self, atom1, atom2, atom3, atom4, args):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        # We don't know how many terms to expect, but it will be a multiple of 3
        if len(args) % 3 != 0:
            raise TypeError('TorsionAngle expects an equal number of '
                            'Amplitudes, phases, and periodicities.')
        nterms = len(args) // 3
        self.amplitude = tuple([float(args[3*i]) for i in range(nterms)])
        self.phase = tuple([float(args[3*i+1]) for i in range(nterms)])
        self.periodicity = tuple([int(args[3*i+2]) for i in range(nterms)])
        if (len(self.amplitude) != len(self.phase) or
                len(self.amplitude) != len(self.periodicity)):
            raise RuntimeError('BUGBUG!! Inconsistent # of terms in torsion')

    def __repr__(self):
        if len(self.amplitude) == 0:
            return '<TorsionAngle [%s-%s-%s-%s]; No Terms>' % (
                        self.atom1.symbol, self.atom2.symbol,
                        self.atom3.symbol, self.atom4.symbol,
            )
        return '<TorsionAngle [%s-%s-%s-%s]; Ampl=%s; Phase=%s; Per=%s>' % (
                        self.atom1.symbol, self.atom2.symbol,
                        self.atom3.symbol, self.atom4.symbol, self.amplitude,
                        self.phase, self.periodicity
        )

    def __str__(self):
        if len(self.amplitude) == 0:
            return ('TorsionAngle %r --- %r --- %r --- %r\n'
                    '   No Terms' % (self.atom1, self.atom2, self.atom3,
                                     self.atom4)
            )
        retstr = 'TorsionAngle %r --- %r --- %r --- %r' % (self.atom1,
                self.atom2, self.atom3, self.atom4)
        seq = range(len(self.amplitude))
        for i, amp, phase, per in enumerate(zip(seq, self.amplitude, self.phase,
                                                self.periodicity)):
            retstr += ('\n   Term %d\n'
                       '   Amplitude = %.4f\n'
                       '       Phase = %.4f\n'
                       ' Periodicity = %.4f' % (i + 1, amp, phase, per))
        return retstr

class TorsionAngleList(_ParamTypeList):
    typeclass = TorsionAngle

class PiTorsion(object):
    """ A Pi-Orbital Torsion parameter """
    def __init__(self, atom1, atom2, amplitude):
        self.atom1, self.atom2 = atom1, atom2
        self.amplitude = float(amplitude)

    def __repr__(self):
        return '<PiTorsion [%s-%s]; Ampl=%.2f>' % (self.atom1.symbol,
                self.atom2.symbol, self.amplitude)

    def __str__(self):
        return ('PiTorsion %r --- %r\n   Amplitude = %.4f' %
                (self.atom1, self.atom2, self.amplitude))

class PiTorsionList(_ParamTypeList):
    typeclass = PiTorsion

class TorsionTorsion(object):
    """ A coupled-torsion parameter (like CMAP) """
    def __init__(self, atom1, atom2, atom3, atom4, atom5, spline1, spline2):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.atom5 = atom5
        self.spline1 = int(spline1)
        self.spline2 = int(spline2)

    def __repr__(self):
        return '<TorsionTorsion [%s-%s-%s-%s-%s]; Spline=(%d,%d)>' % (
                        self.atom1.symbol, self.atom2.symbol,
                        self.atom3.symbol, self.atom4.symbol,
                        self.atom5.symbol, self.spline1, self.spline2,
        )

    def __str__(self):
        return ('TorsionTorsion %r --- %r --- %r ---\n'
                '                      %r --- %r\n'
                '   Spline Grid = (%d x %d)' %
                    (self.atom1, self.atom2, self.atom3, self.atom4,
                     self.atom5, self.spline1, self.spline2)
        )

class TorsionTorsionList(_ParamTypeList):
    typeclass = TorsionTorsion

class TorsionTorsionGrid(OrderedDict):
    """
    The interpolation grid of a coupled-torsion correction map. Every unique
    grid is cached and if a duplicate grid is instantiated, a reference to the
    original grid is returned. As a result, all unique TorsionTorsionGrid
    instances are singletons and should be compared for equality with "is"
    """
    # Because these grids are space-consuming, we only hold onto unique grids
    _typelist = list()

    @classmethod
    def new(cls, data):
        inst = cls()
        for d in data:
            inst[tuple(d[:2])] = d[2]
        # Potentially expensive comparison of all grids.
        for ttg in TorsionTorsionGrid._typelist:
            if ttg == inst:
                return ttg
        TorsionTorsionGrid._typelist.append(inst)
        return inst
   
    def __eq__(self, other):
        if self.keys() != other.keys(): return False
        for key in self:
            if abs(self[key] - other[key]) > 1e-8: return False
        return True
   
    def __ne__(self, other):
        return not self == other

    # No comparisons are implemented
    def __gt__(self, other):
        raise NotImplemented('TorsionTorsionGrid instances are not '
                             'well-ordered')
    __lt__ = __ge__ = __le__ = __gt__

class AtomicMultipole(object):
    """ Atomic multipole parameter """
    def __init__(self, atom, frame, definition, moments):
        self.atom = atom
        self.frame = [int(i) for i in frame if str(i).strip()]
        self.definition = definition.strip()
        self.moment = [float(x) for x in moments]

    def __repr__(self):
        return '<AtomicMultipole [%s] %s; Frame=%s; Moments=%s>' % (
                        self.atom.symbol, self.definition,
                        self.frame, self.moment
        )

    def __str__(self):
        return ('AtomicMultipole %r "%s"\n'
                '     Frame = %s\n'
                '   Moments = %8.5f\n'
                '             %8.5f %8.5f %8.5f\n'
                '             %8.5f\n'
                '             %8.5f %8.5f\n'
                '             %8.5f %8.5f %8.5f' % (
                        self.atom, self.definition, self.frame,
                        self.moment[0], self.moment[1], self.moment[2],
                        self.moment[3], self.moment[4], self.moment[5],
                        self.moment[6], self.moment[7], self.moment[8],
                        self.moment[9])
        )

class AtomicMultipoleList(_ParamTypeList):
    typeclass = AtomicMultipole

class DipolePolarizability(object):
    """ A dipole polarizability parameter """
    def __init__(self, atom, alpha, group):
        self.atom = atom
        self.alpha = float(alpha)
        self.group = [int(g) for g in group]

    def __repr__(self):
        return '<DipolePolarizability [%s]; Alpha=%.2f; Group=%s>' % (
                self.atom.symbol, self.alpha, self.group)

    def __str__(self):
        return ('DipolePolarizability %r\n'
                '    Alpha = %.4f\n'
                '    Group = %s' % (self.atom, self.alpha, self.group))


class DipolePolarizabilityList(_ParamTypeList):
    typeclass = DipolePolarizability
