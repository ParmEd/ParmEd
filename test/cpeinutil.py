#!/usr/bin/env python
from __future__ import division

import os
import sys

__author__ = 'Vinicius Wilian D Cruzeiro: work adapted from Jason Swails'
__version__ = '16.0'

try:
    from argparse import ArgumentParser
    from parmed.amber import AmberParm
    from parmed.amber import titratable_residues as residues
    from parmed.exceptions import ParmedError, AmberError
    from parmed.tools.actions import changeRadii, change
except ImportError:
    if sys.version_info < (2, 7):
        raise ImportError('%s requires Python 2.7 or later' %
                          (os.path.split(sys.argv[0])[1]))
    amberhome = os.getenv('AMBERHOME') or '$AMBERHOME'
    raise ImportError('Could not import Amber Python modules. Please make sure '
                      'you have sourced %s/amber.sh (if you are using sh/ksh/'
                      'bash/zsh) or %s/amber.csh (if you are using csh/tcsh)' %
                      (amberhome, amberhome))
else:
    TitratableResidueList = residues.TitratableResidueList
    LineBuffer = residues._LineBuffer

parser = ArgumentParser(epilog='''This program will read a topology file and
                        generate a cpein file for constant pH and Redox Potential simulations with
                        sander or pmemd''', usage='%(prog)s [Options]')
parser.add_argument('-v', '--version', action='version', version='%s: %s' %
                    (parser.prog, __version__))
parser.add_argument('-d', '--debug', dest='debug', action='store_const',
                    help='Enable verbose tracebacks to debug this program',
                    const=True, default=False)
group = parser.add_argument_group('Output files')
group.add_argument('-o', '--output', dest='output', metavar='FILE',
                   help='Output file. Defaults to standard output')
group = parser.add_argument_group('Required Arguments')
group.add_argument('-p', dest='prmtop', metavar='FILE', required=False,
                   help='Topology file to be used in constant pH and Redox Potential simulation',
                   type=str, default='prmtop')
group = parser.add_argument_group('Simulation Options')
group.add_argument('-igb', dest='igb', metavar='IGB', required=False, type=int,
                   help='Generalized Born model which you intend to use to '
                   'evaluate dynamics (or protonation state swaps). Default '
                   'is 2.', default=2)
group.add_argument('-intdiel', dest='intdiel', metavar='DIEL', type=float,
                   default=1.0, help='''Internal dielectric constant to use in
                   the evaluation of the GB potential. Default %(default)s.''')
group = parser.add_argument_group('Residue Selection Options')
group.add_argument('-resnames', dest='resnames', metavar='RES', nargs='*',
                   help='Residue names to include in CPEIN file', default=None)
group.add_argument('-notresnames', dest='notresnames', metavar='RES', nargs='*',
                   help='Residue names to exclude from CPEIN file', default=None)
group.add_argument('-resnums', dest='resnums', metavar='NUM',
                   nargs='*', help='Residue numbers to include in CPEIN file',
                   default=None)
group.add_argument('-notresnums', dest='notresnums', nargs='*', metavar='NUM',
                 help='Residue numbers to exclude from CPEIN file', default=None)
group = parser.add_argument_group('System Information')
group.add_argument('-states', dest='resstates', metavar='NUM', nargs='*',
                 help='List of default states to assign to titratable residues')
group.add_argument('-system', dest='system', metavar='<system name>',
                   help='Name of system to titrate. No effect on simulation.',
                   default='Unknown')
group = parser.add_argument_group('Residue Information', '''If any options here
              are used, no CPEIN file will be written. These arguments take
              precedence and are mutually exclusive with each other.''')
group.add_argument('--describe', dest='descres', metavar='RESNAME', nargs='*',
                   help='Print out the details of given residues', default=None)
group.add_argument('-l', '--list', dest='list', default=False,
                   action='store_true', help='List all titratable residues')

def print_residues(resnames,mode):
    for resname in resnames:
        if not hasattr(residues, resname):
            print ('%s is not titratable\n' % resname)
            sys.exit(0)
        print (str(getattr(residues, resname)) + '\n')

def list_residues():
    """ Lists all titratable residues defined in residues.py """
    line = LineBuffer(sys.stdout)
    strarray = []
    for resname in residues.titratable_residues:
        strarray.append(resname)
    line.add_words(', '.join(strarray).split(),
                   space_delimited=True)
    line.flush()

def process_arglist(arglist, argtype):
   """
   This processes an argument list with an arbitrary number of arguments that
   may or may not be comma-delimited (with any number of spaces)
   """
   # If the argument list is not set, just return None
   if arglist is None:
      return None
   # Otherwise, process the arguments
   processed_args = []
   for arg in arglist:
      # Delete any spaces, split on commas, and add this to processed arg list
      for arg in arg.replace(' ', '').split(','):
         if not arg: continue
         try:
            processed_args.append(argtype(arg))
         except ValueError:
            raise AmberError('Expected type %s for argument. Got %s' %
                             (argtype.__name__, arg))

   return processed_args

def main(opt):
    # Check all of the arguments
    if not os.path.exists(opt.prmtop):
        raise AmberError('prmtop file (%s) is missing' % opt.prmtop)

    # Process the arguments that take multiple args
    resstates = process_arglist(opt.resstates, int)
    resnums = process_arglist(opt.resnums, int)
    notresnums = process_arglist(opt.notresnums, int)
    resnames = process_arglist(opt.resnames, str)
    notresnames = process_arglist(opt.notresnames, str)

    if not opt.igb in (1, 2, 5, 7, 8):
        raise AmberError('-igb must be 1, 2, 5, 7, or 8!')

    if resnums is not None and notresnums is not None:
        raise AmberError('Cannot specify -resnums and -notresnums together')

    if resnames is not None and notresnames is not None:
        raise AmberError('Cannot specify -resnames and -notresnames together')

    if opt.intdiel != 1 and opt.intdiel != 2:
        raise AmberError('-intdiel must be either 1 or 2 currently')

    # Set the list of residue names we will be willing to titrate
    titratable_residues = []
    if notresnames is not None:
        for resname in residues.titratable_residues:
            if resname in notresnames: continue
            titratable_residues.append(resname)
    elif resnames is not None:
        for resname in resnames:
            if not resname in residues.titratable_residues:
                raise AmberError('%s is not a titratable residue!' % resname)
            titratable_residues.append(resname)
    else:
        for resname in residues.titratable_residues:
            titratable_residues.append(resname)

    solvent_ions = ['WAT', 'Na+', 'Br-', 'Cl-', 'Cs+', 'F-', 'I-', 'K+', 'Li+',
                    'Mg+', 'Rb+', 'CIO', 'IB', 'MG2']

    # Make sure we still have a couple residues
    if len(titratable_residues) == 0:
        raise AmberError('No titratable residues fit your criteria!')

    # Load the topology file
    parm = AmberParm(opt.prmtop)

    # Replace an un-set notresnums with an empty list so we get __contains__()
    if notresnums is None:
        notresnums = []

    # If we have a list of residue numbers, check that they're all titratable
    if resnums is not None:
        for resnum in resnums:
            if resnum > parm.ptr('nres'):
                raise AmberError('%s only has %d residues. (%d chosen)' %
                                     (parm, parm.ptr('nres'), resnum))
            if resnum <= 0:
                raise AmberError('Cannot select negative residue numbers.')
            resname = parm.parm_data['RESIDUE_LABEL'][resnum-1]
            if not resname in titratable_residues:
                raise AmberError('Residue number %s [%s] is not titratable'
                                     % (resnum, resname))
    else:
        # Select every residue except those in notresnums
        resnums = []
        for resnum in range(1, parm.ptr('nres') + 1):
            if resnum in notresnums: continue
            resnums.append(resnum)

    solvated = False
    first_solvent = 0
    if 'WAT' in parm.parm_data['RESIDUE_LABEL']:
        solvated = True
        for i, res in enumerate(parm.parm_data['RESIDUE_LABEL']):
            if res in solvent_ions:
                first_solvent = parm.parm_data['RESIDUE_POINTER'][i]
                break
    main_reslist = TitratableResidueList(system_name=opt.system,
                        solvated=solvated, first_solvent=first_solvent)
    trescnt = 0
    for resnum in resnums:
        resname = parm.parm_data['RESIDUE_LABEL'][resnum-1]
        if not resname in titratable_residues: continue
        res = getattr(residues, resname)
        # Filter out termini (make sure the residue in the prmtop has as many
        # atoms as the titratable residue defined in residues.py)
        if resnum == parm.ptr('nres'):
            natoms = (parm.ptr('natom') + 1 -
                      parm.parm_data['RESIDUE_POINTER'][resnum-1])
        else:
            natoms = (parm.parm_data['RESIDUE_POINTER'][resnum] -
                      parm.parm_data['RESIDUE_POINTER'][resnum-1])
        if natoms != len(res.atom_list): continue
        # If we have gotten this far, add it to the list.
        main_reslist.add_residue(res, resnum,
                                 parm.parm_data['RESIDUE_POINTER'][resnum-1])
        trescnt += 1

    # Prints a warning if the number of titratable residues is larger than 50
    if trescnt > 50: sys.stderr.write('Warning: Your CPEIN file has more than 50 titratable residues! pmemd and sander have a\n'
                                      '         default limit of 50 titrable residues, thus this CPEIN file can only be used\n'
                                      '         if the definitions for this limit are modified at the top of\n'
                                      '         $AMBERHOME/src/pmemd/src/constante.F90 or $AMBERHOME/AmberTools/src/sander/constantphe.F90.\n'
                                      '         AMBER needs to be recompilied after these files are modified.\n')

    # Set the states if requested
    if resstates is not None:
        main_reslist.set_states(resstates)

    # Open the output file
    if opt.output is None:
        output = sys.stdout
    else:
        output = open(opt.output, 'w')

    main_reslist.write_cpin(output, opt.igb, opt.intdiel, False, "phredox")

    if opt.output is not None:
        output.close()

    sys.stderr.write('CPEIN generation complete!\n')

if __name__ == '__main__':
    opt = parser.parse_args()

    # List all residues
    if opt.list:
        list_residues()
        sys.exit(0)

    # Describe requested residues
    if opt.descres is not None:
        if len(opt.descres) == 0 or (len(opt.descres) == 1 and
                                     opt.descres[0].upper() == 'ALL'):
            print_residues(residues.titratable_residues,0)
        else:
            opt.descres = process_arglist(opt.descres, str)
            print_residues(opt.descres,1)
        sys.exit(0)

    # Go ahead and make the CPEIN file.
    try:
        main(opt)
    except ParmedError as e:
        sys.exit('%s: %s' % (type(e).__name__, e))
    sys.exit(0)
