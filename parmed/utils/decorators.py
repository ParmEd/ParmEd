""" A list of helpful decorators for use throughout ParmEd """
from functools import wraps

__all__ = ['needs_openmm']

try:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    HAS_OPENMM = True
    try:
        from simtk.openmm.app.internal import unitcell
    except ImportError:
        unitcell = None
        HAS_OPENMM = False
except ImportError:
    HAS_OPENMM = False
else:
    del mm, app, unitcell
import warnings

def needs_openmm(fcn):
    @wraps(fcn)
    def new_fcn(*args, **kwargs):
        if not HAS_OPENMM:
            raise ImportError('Could not find or import OpenMM version 6.3+')
        return fcn(*args, **kwargs)

    return new_fcn

def deprecated(fcn):
    @wraps(fcn)
    def new_fcn(*args, **kwargs):
        warnings.warn(
            f'{fcn.__name__} is deprecated and will be removed in the future', DeprecationWarning
        )
        return fcn(*args, **kwargs)
    return new_fcn
