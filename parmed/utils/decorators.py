""" A list of helpful decorators for use throughout ParmEd """

__all__ = ['needs_openmm']

from parmed.utils.six import wraps
try:
    import simtk.openmm as mm
    HAS_OPENMM = True
    del mm
except ImportError:
    HAS_OPENMM = False

def needs_openmm(fcn):
    global HAS_OPENMM
    @wraps(fcn)
    def new_fcn(*args, **kwargs):
        if not HAS_OPENMM:
            raise ImportError('Could not find or import OpenMM')
        return fcn(*args, **kwargs)

    return new_fcn

