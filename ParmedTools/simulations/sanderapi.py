"""
This module contains functionality needed to compute energies and forces with
the sander-Python bindings
"""
from __future__ import division
import sys
import warnings

try:
    import sander
except ImportError:
    sander = None

from ParmedTools.exceptions import (SimulationError, SimulationWarning,
               UnhandledArgumentWarning)

HAS_SANDER = sander is not None

def energy(parm, args, output=sys.stdout):
    """
    Compute a single-point energy using sander and print the result to the
    desired output
    """
    global HAS_SANDER
    if not HAS_SANDER:
        raise SimulationError('Could not import sander')

    cutoff = args.get_key_float('cutoff', None)
    igb = args.get_key_int('igb', 5)
    saltcon = args.get_key_float('saltcon', 0.0)
    do_ewald = args.has_key('Ewald')
    vdw_longrange = not args.has_key('nodisper')
    has_1264 = 'LENNARD_JONES_CCOEF' in parm.parm_data

    # Get any unmarked arguments
    unmarked_cmds = args.unmarked()
    if len(unmarked_cmds) > 0:
        warnings.warn("Un-handled arguments: " + ' '.join(unmarked_cmds),
                      UnhandledArgumentWarning)

    gbmeth, kappa = None, 0.0
    if parm.ptr('ifbox') == 0:
        if not igb in (0, 1, 2, 5, 6, 7, 8):
            raise SimulationError('Bad igb value. Must be 0, 1, 2, 5, '
                                  '6, 7, or 8')
        # Force vacuum electrostatics down the GB code path
        if igb == 0:
            igb = 6
        inp = sander.gas_input(igb)
        if cutoff is None:
            cutoff = 1000.0
        if cutoff <= 0:
            raise SimulationError('cutoff must be > 0')
        inp.cut = cutoff
        if saltcon < 0:
            raise SimulationError('salt concentration must be >= 0')
        inp.saltcon = saltcon
    elif parm.ptr('ifbox') > 0:
        inp = sander.pme_input()
        if cutoff is None:
            cutoff = 8.0
        elif cutoff <= 0:
            raise SimulationError('cutoff must be > 0')
        inp.cut = cutoff
        inp.ew_type = int(do_ewald)
        inp.vdwmeth = int(vdw_longrange)
        inp.lj1264 = int(has_1264)

    if not hasattr(parm, 'coords'):
        raise SimulationError('No coordinates are loaded')
    # Time to set up sander
    sander.setup(parm, parm.coords, parm.box, inp)

    e, f = sander.energy_forces()

    if parm.chamber:
        output.write('Bond          = %20.7f     Angle         = %20.7f\n'
                     'Dihedral      = %20.7f     Urey-Bradley  = %20.7f\n'
                     'Improper      = %20.7f     ' % (e.bond, e.angle,
                     e.dihedral, e.angle_ub, e.imp))
        if parm.has_cmap:
            output.write('CMAP         = %20.7f\n' % e.cmap)
        output.write('1-4 vdW       = %20.7f     1-4 Elec.     = %20.7f\n'
                     'Lennard-Jones = %20.7f     Electrostatic = %20.7f\n'
                     'TOTAL         = %20.7f\n' % (e.vdw_14, e.elec_14,
                     e.vdw, e.elec, e.tot))
    else:
        output.write('Bond     = %20.7f     Angle    = %20.7f\n'
                     'Dihedral = %20.7f     1-4 vdW  = %20.7f\n'
                     '1-4 Elec = %20.7f     vdWaals  = %20.7f\n'
                     'Elec.    = %20.7f' % (e.bond, e.angle, e.dihedral,
                      e.vdw_14, e.elec_14, e.vdw, e.elec))
        if igb != 0 and inp.ntb == 0:
            output.write('     Egb      = %20.7f' % e.gb)
        elif e.hbond != 0:
            output.write('     EHbond   = %20.7f' % e.hbond)
        output.write('\nTOTAL    = %20.7f\n' % e.tot)

    sander.cleanup()
