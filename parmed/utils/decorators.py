""" A list of helpful decorators for use throughout ParmEd """

__all__ = ['needs_openmm']

from parmed.utils.six import wraps
try:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
else:
    del mm, app
import warnings

def needs_openmm(fcn):
    global HAS_OPENMM
    @wraps(fcn)
    def new_fcn(*args, **kwargs):
        if not HAS_OPENMM:
            raise ImportError('Could not find or import OpenMM')
        return fcn(*args, **kwargs)

    return new_fcn

def deprecated(fcn):
    @wraps(fcn)
    def new_fcn(*args, **kwargs):
        warnings.warn('%s is deprecated and will be removed in the future' %
                      fcn.__name__, DeprecationWarning)
        return fcn(*args, **kwargs)
    return new_fcn
