""" 
All of the prmtop actions used in PARMED. Each class is a separate action.
"""
from __future__ import division

import sys
from ParmedTools.exceptions import (WriteOFFError, ParmError, ParmWarning,
              ParmedMoleculeError, ChangeStateError, ChangeLJPairError,
              ParmedChangeError, SetParamError, DeleteDihedralError, NoArgument,
              NonexistentParm, AmbiguousParmError, IncompatibleParmsError,
              ArgumentError, AddPDBError, AddPDBWarning, HMassRepartitionError,
              SimulationError, UnhandledArgumentWarning, SeriousParmWarning,
              FileExists, NonexistentParmWarning, LJ12_6_4Error, ChamberError,
              FileDoesNotExist, InputError,
#             CoarseGrainError,
                                   )
from ParmedTools.argumentlist import ArgumentList
from ParmedTools.parmlist import ParmList, ChamberParm
from chemistry.amber.mask import AmberMask
from chemistry.amber.readparm import AmberFormat
from chemistry.amber.topologyobjects import (Bond, BondType, Angle, AngleType,
                                             Dihedral, DihedralType)
from chemistry.exceptions import ChemError
from chemistry.periodic_table import Element as _Element
from compat24 import any
import math
import os
import warnings

# Add a help dictionary entry for each additional Action added to this class!
# Each help entry should be a list with 2 elements: [Usage, description]
Usages = {
              'help' : 'help [<action>]',
           'parmout' : 'parmout <prmtop_name> [<inpcrd_name>] [netcdf]',
      'setoverwrite' : 'setOverwrite [True|False]', 
       'writefrcmod' : 'writeFrcmod <frcmod_name>',
        'loadrestrt' : 'loadRestrt <restrt_filename>',
          'writeoff' : 'writeOFF <OFF_filename>',
       'changeradii' : 'changeRadii <radii_set>',
      'changeljpair' : 'changeLJPair <mask1> <mask2> <Rmin> <epsilon>',
    'changelj14pair' : 'changeLJ14Pair <mask1> <mask2> <Rmin> <epsilon>',
     'checkvalidity' : 'checkValidity',
            'change' : 'change <property> <mask> <new_value> [quiet]',
         'printinfo' : 'printInfo <flag>',
         'addljtype' : 'addLJType <mask> [radius <new_radius>] '
                       '[epsilon <new_epsilon>] [radius_14 <new_radius14>] '
                       '[epsilon_14 <new_epsilon14>]',
           'outparm' : 'outparm <prmtop_name> [<inpcrd_name>] [netcdf]',
      'printljtypes' : 'printLJTypes [<mask>|<type name>]',
              'scee' : 'scee <scee_value>',
              'scnb' : 'scnb <scnb_value>',
'changeljsingletype' : 'changeLJSingleType <mask> <radius> <depth>',
      'printdetails' : 'printDetails <mask>',
        'printflags' : 'printFlags',
     'printpointers' : 'printPointers',
        'printbonds' : 'printBonds <mask>',
       'printangles' : 'printAngles <mask>',
    'printdihedrals' : 'printDihedrals <mask>',
      'setmolecules' : 'setMolecules [solute_ions True|False]',
  'combinemolecules' : 'combineMolecules <mol_id1> [<mol_id2>]',
                'go' : 'go',
              'quit' : 'quit',
#   'addcoarsegrain' : 'addCoarseGrain <parameter_file>',
   'changeprotstate' : 'changeProtState <mask> <state #>',
         'netcharge' : 'netCharge [<mask>]',
             'strip' : 'strip <mask>',
     'definesolvent' : 'defineSolvent <residue list>',
     'addexclusions' : 'addExclusions <mask1> <mask2>',
       'adddihedral' : 'addDihedral <mask1> <mask2> <mask3> <mask4> <phi_k> '
                       '<per> <phase> [<scee>] [<scnb>] [type <type>]',
           'setbond' : 'setBond <mask1> <mask2> <k> <Req>',
          'setangle' : 'setAngle <mask1> <mask2> <mask3> <k> <THETeq>',
   'addatomicnumber' : 'addAtomicNumber',
    'deletedihedral' : 'deleteDihedral <mask1> <mask2> <mask3> <mask4>',
        'deletebond' : 'deleteBond <mask1> <mask2>',
     'printljmatrix' : 'printLJMatrix <mask>|<index>',
            'source' : 'source <file>',
              'parm' : 'parm <filename> [<filename> [<filename> ...]] || '
                       'parm copy <filename>|<index> || parm select '
                       '<filename>|<index>',
                'ls' : 'ls [Unix ls options]',
                'cd' : 'cd <directory>',
         'listparms' : 'listParms',
           'timerge' : 'tiMerge <mol1mask> <mol2mask> <scmask1> <scmask2>'
                            ' [<scmask1N>] [<scmask2N>] [tol <tol>]',
       'interpolate' : 'interpolate <nparm> [parm2 <other_parm>] [eleconly]'
                       ' [prefix <prefix>] [startnum <num>]',
           'summary' : 'summary',
             'scale' : 'scale <FLAG> <factor>',
              'lmod' : 'lmod',
            'addpdb' : 'addPDB <filename> [elem] [strict] [allicodes]',
         'deletepdb' : 'deletePDB',
         'add12_6_4' : 'add12_6_4 [<divalent ion mask>] '
                       '[c4file <C4 Param. File> | watermodel <water model>] '
                       '[polfile <Pol. Param File>] [tunfactor <tunfactor>]',
  'hmassrepartition' : 'HMassRepartition [<mass>] [dowater]',
            'openmm' : 'OpenMM -p <parm>|<parm index> [sander/pmemd options] '
                       '-platform <platform> -precision <precision model> '
                       '[dcd] [progress] [script <script_file.py>] [norun]',
            'energy' : 'energy [cutoff <cut>] [[igb <IGB>] [saltcon <conc>] | '
                       '[Ewald]] [nodisper] [applayer] [platform <platform>] '
                       '[precision <precision model>] [decompose]',
       'fixtopology' : 'fixTopology',
           'chamber' : 'chamber -top <CHARMM.top> -param <CHARMM.par> [-str '
                       '<CHARMM.str>] -psf <CHARMM.psf> [-crd <CHARMM.pdb>] '
                       '[-nocmap] [usechamber] [box a,b,c[,alpha,beta,gamma]]',
          'minimize' : 'minimize [cutoff <cut>] [[igb <IGB>] [saltcon <conc>]] '
                       '[[restrain <mask>] [weight <k>]] [norun] '
                       '[script <script_file.py>] [platform <platform>] '
                       '[precision <precision model>]',
#             'heat' : 'heat [cutoff <cut>] [[igb <IGB>] [saltcon <conc>]] '
#                      '[[restrain <mask>] [weight <k>]] [langevin | '
#                      'anderson] [nvt | npt] [anisotropic] [norun] [script '
#                      '<script_file.py>] [platform <platform>] '
#                      '[precision <precision model>]',
#               'md' : 'MD [cutoff <cut>] [[igb <IGB>] [saltcon <conc>] | '
#                      '[Ewald]] [[restrain <mask>] [weight <k>]] [npt | nve | '
#                      'nvt] [langevin | anderson] [anisotropic] [norun] '
#                      '[script <script_file.py>] [platform <platform>] '
#                      '[precision <precision model>]',
}

# Add help and go as a class here to basically hold docstrings for the help
# function

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

lawsuit = object

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Action(lawsuit):
    """ Base action class """
    stderr = sys.stderr
    # Does this action need a populated parm_list? If yes, bail out if it's
    # unpopulated
    needs_parm = True
    # Do we allow any action to overwrite existing files? Static variable.
    overwrite = True
    # Set a list of which classes of topology files is supported by this action.
    # Use class names rather than "isinstance" since ChamberParm and AmoebaParm
    # both inherit from AmberParm
    supported_classes = ('AmberParm', 'AmoebaParm', 'ChamberParm')
    def __init__(self, input_parm, arg_list=None, *args, **kwargs):
        """ Constructor """
        # Is this action 'valid' (i.e., can we run 'execute')?
        self.valid = False
        # Accept both an AmberParm or ParmList instance to simplify the API
        if isinstance(input_parm, ParmList):
            self.parm_list = input_parm
        elif isinstance(input_parm, AmberFormat):
            self.parm_list = ParmList()
            self.parm_list.add_parm(input_parm)
        # Set active parm
        self.parm = self.parm_list.parm # Current active parm
        if self.needs_parm and self.parm_list.empty():
            raise ParmError('Action requires a loaded topology file')
        # Fill the argument list
        if args or kwargs:
            if arg_list is None:
                arg_list = ''
            else:
                arg_list = '%s ' % arg_list
            arg_list += ' '.join([str(a) for a in args])
            for kw in kwargs:
                arg_list += ' %s %s ' % (kw, kwargs[kw])
        elif arg_list is None:
            arg_list = ArgumentList('')
        # If our arg_list is a string, convert it to an ArgumentList
        # (but coerce to string if it is string-like)
        try:
            arg_list = arg_list.decode()
            arg_list = ArgumentList(arg_list)
        except AttributeError:
            if isinstance(arg_list, str):
                arg_list = ArgumentList(arg_list)
        # Now that we have the argument list, see if we have requested a
        # specific parm. If so, set self.parm equal to that object, instead
        parm = arg_list.get_key_string('parm', None)
        # If it is an integer, see if it is a valid index. If it's not a valid
        # index, try it as a string. Otherwise just try it as a string.
        if parm is not None:
            try:
                if int(parm) >= 0 and int(parm) < len(self.parm_list):
                    print 'Using parm %s' % self.parm_list[int(parm)]
                    self.parm = self.parm_list[int(parm)]
                elif parm in self.parm_list:
                    print 'Using parm %s' % self.parm_list[parm]
                    self.parm = self.parm_list
                else:
                    warnings.warn('Cannot find parm %s. Skipping this action'
                                  % parm, SeriousParmWarning)
                    return
            except ValueError:
                if parm in self.parm_list:
                    print 'Using parm %s' % parm
                    self.parm = self.parm_list[parm]
                else:
                    warnings.warn('Cannot find parm %s. Skipping this action'
                                  % parm, SeriousParmWarning)
                    return

        # Make sure our topology file type is supported. OpenMM* classes should
        # be treated the same as their non-OpenMM counterparts
        typename = type(self.parm).__name__.replace('OpenMM', '')
        if self.needs_parm and typename not in self.supported_classes:
            raise ParmError('%s topologies are not supported by this action' %
                            type(self.parm).__name__)
        try:
            self.init(arg_list)
        except NoArgument:
            try:
                usage = Usages[type(self).__name__]
                cmdname = usage.split()[0]
            except KeyError:
                usage = ''
                cmdname = type(self).__name__
            Action.stderr.write("Bad command %s:\n\t%s\n" % (cmdname, usage))
            self.__str__ = Action.__str__
            self.execute = Action.execute
            return

        # Check any unmarked commands
        unmarked_cmds = arg_list.unmarked()
        if len(unmarked_cmds) > 0:
            warnings.warn(' '.join(unmarked_cmds), UnhandledArgumentWarning)
        self.valid = True

    def init(self, arg_list):
        """ This should be overridden if anything needs to be done """
        raise NotImplemented('init must be overridden by Action subclass')

    def execute(self):
        """ Commands involved in executing the action """
        pass

    def __str__(self):
        return ''
   
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class parmout(Action):
    """
    Final prmtop written after all actions are complete
    """
    def init(self, arg_list):
        self.filename = arg_list.get_next_string()
        self.rst_name = arg_list.get_next_string(optional=True)
        if arg_list.has_key('netcdf'):
            self.netcdf = True
        else:
            self.netcdf = None

    def __str__(self):
        if self.rst_name is not None:
            return 'Outputting Amber topology file %s and restart %s' % (
                        self.filename, self.rst_name)
        return 'Outputting Amber topology file %s' % self.filename

    def execute(self):
        if not Action.overwrite and os.path.exists(self.filename):
            raise FileExists('%s exists; not overwriting.' % self.filename)
        if self.rst_name is not None:
            if not Action.overwrite and os.path.exists(self.rst_name):
                raise FileExists('%s exists; not overwriting.' % self.rst_name)
        self.parm.writeParm(self.filename)
        if self.rst_name is not None:
            self.parm.writeRst7(self.rst_name, netcdf=self.netcdf)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class setoverwrite(Action):
    """
    Necessary to overwrite original topology file name.
    """
    def init(self, arg_list):
        # see if we got an argument
        arg = arg_list.get_next_string(optional=True)
        if arg is not None:
            if not arg.lower() in ('false', 'true'):
                warnings.warn("setOverwrite: unrecognized argument. "
                              "Assuming False", SeriousParmWarning)
        self._overwrite = (arg is None or arg.lower() == "true")

    def __str__(self):
        if self._overwrite:
            return 'Files are overwritable'
        else:
            return 'Files are NOT overwritable'

    def execute(self):
        Action.overwrite = self._overwrite

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class writefrcmod(Action):
    """
    Writes an frcmod file from all of the parameters in the topology file.
    """
    supported_classes = ('AmberParm',)
    def init(self, arg_list):
        self.frcmod_name = arg_list.get_next_string(optional=True)
        if self.frcmod_name is None: self.frcmod_name = 'frcmod'
        # Emit warnings about 10-12 prmtops if we detect any 10-12 parameters
        try:
            if (len(self.parm.parm_data['HBOND_ACOEF']) + 
                len(self.parm.parm_data['HBOND_BCOEF']) +
                len(self.parm.parm_data['HBCUT']) > 0):
                warnings.warn('Frcmod dumping does not work with 10-12 prmtops',
                              SeriousParmWarning)
        except KeyError:
            pass

    def __str__(self):
        return 'Dumping FRCMOD file %s with parameters from %s' % (
                self.frcmod_name, self.parm)

    def execute(self):
        """ Writes the frcmod file """
        from chemistry.amber.parameters import ParameterSet
        if not Action.overwrite and os.path.exists(self.frcmod_name):
            raise FileExists('%s exists; not overwriting' % self.frcmod_name)
        parmset = ParameterSet()
        parmset.load_from_parm(self.parm)
        frcmod = open(self.frcmod_name, 'w')
        frcmod.write('Force field parameters from %s\n' % self.parm)
        parmset.write(frcmod)
        frcmod.close()

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class loadrestrt(Action):
    """
    Loads a restart file so we have coordinates. Necessary for distance-based
    mask criteria and writeOFF
    """
    def init(self, arg_list):
        self.rst_name = arg_list.get_next_string()

    def __str__(self):
        return 'Loading restart file %s' % self.rst_name

    def execute(self):
        self.parm.LoadRst7(self.rst_name)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class writeoff(Action):
    """
    Writes an Amber OFF Library with all of the residues found in the topology
    """
    supported_classes = ('AmberParm',)
    def init(self, arg_list):
        self.off_file = arg_list.get_next_string()

    def __str__(self):
        return 'Writing Amber OFF file %s' % self.off_file

    def execute(self):
        if not Action.overwrite and os.path.exists(self.off_file):
            raise FileExists('%s exists; not overwriting' % self.off_file)
        try:
            self.parm.rst7
        except:
            raise WriteOFFError('You must load a restart for WriteOFF!')

        self.parm.writeOFF(self.off_file)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class changeradii(Action):
    """
    Changes intrinsic GB radii to the specified set: bondi, mbondi, amber6,
    mbondi2, or mbondi3
    """
    supported_classes = ('AmberParm', 'ChamberParm')
    def init(self, arg_list):
        self.radii = arg_list.get_next_string()

    def __str__(self):
        return 'Changing PB/GB radii to %s' % self.radii

    def execute(self):
        from ParmedTools.changeradii import ChRad
        # Add RADIUS_SET to prmtop if it's not there already, and a blank 
        # description, since it's about to be set here
        if not 'RADIUS_SET' in self.parm.flag_list:
            self.parm.addFlag('RADIUS_SET', '1a80', num_items=1)
        # Add RADII prmtop section if it doesn't exist already. Just make it a
        # zeroed array, since it's all about to be set here
        if not 'RADII' in self.parm.flag_list:
            self.parm.addFlag('RADII', '5E16.8',
                              num_items=self.parm.ptr('natom'))
        if not 'SCREEN' in self.parm.flag_list:
            self.parm.addFlag('SCREEN', '5E16.8',
                              num_items=self.parm.ptr('natom'))
        ChRad(self.parm, self.radii)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class changeljpair(Action):
    """
    Changes a particular Lennard Jones pair based on a given (pre-combined)
    epsilon/Rmin
    """
    supported_classes = ('AmberParm', 'ChamberParm')
    def init(self, arg_list):
        self.mask1 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask2 = AmberMask(self.parm, arg_list.get_next_mask())
        self.rmin = arg_list.get_next_float()
        self.eps = arg_list.get_next_float()

    def __str__(self):
        return ('Setting LJ %s-%s pairwise interaction to have ' +
                'Rmin = %16.5f and Epsilon = %16.5f') % (self.mask1, self.mask2,
                self.rmin, self.eps)

    def execute(self):
        from ParmedTools.changeljpair import ChLJPair
        selection1 = self.mask1.Selection()
        selection2 = self.mask2.Selection()
        if sum(selection1) == 0 or sum(selection2) == 0:
            Action.stderr.write('Skipping empty masks in changeLJPair\n')
            return 0
        # Make sure we only selected 1 atom type in each mask
        attype1 = None
        attype2 = None
        for i in range(self.parm.ptr('natom')):
            if selection1[i] == 1:
                if not attype1:
                    attype1 = self.parm.parm_data['ATOM_TYPE_INDEX'][i]
                else:
                    if attype1 != self.parm.parm_data['ATOM_TYPE_INDEX'][i]:
                        raise ChangeLJPairError(
                                      'First mask matches multiple atom types!')
            if selection2[i] == 1:
                if not attype2:
                    attype2 = self.parm.parm_data['ATOM_TYPE_INDEX'][i]
                else:
                    if attype2 != self.parm.parm_data['ATOM_TYPE_INDEX'][i]:
                        raise ChangeLJPairError(
                                     'Second mask matches multiple atom types!')
        ChLJPair(self.parm, attype1, attype2, self.rmin, self.eps)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class changelj14pair(Action):
    """
    Changes a particular 1-4 Lennard Jones pair based on a given (pre-combined)
    epsilon/Rmin. Only valid for CHAMBER prmtops
    """
    supported_classes = ('ChamberParm',)
    def init(self, arg_list):
        # Make sure this is a chamber prmtop
        if not self.parm.chamber:
            raise ChangeLJPairError('Changing 1-4 NB pairs makes no sense for '
                                    'non-chamber created prmtops!')
        # If not, initiate instance data
        self.mask1 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask2 = AmberMask(self.parm, arg_list.get_next_mask())
        self.rmin = arg_list.get_next_float()
        self.eps = arg_list.get_next_float()

    def __str__(self):
        if self.parm.chamber:
            return ('Setting LJ 1-4 %s-%s pairwise interaction to have '
                    '1-4 Rmin = %16.5f and 1-4 Epsilon = %16.5f' %
                    (self.mask1, self.mask2, self.rmin, self.eps))
        return 'Not a chamber topology. Nothing to do.'

    def execute(self):
        from ParmedTools.changeljpair import ChLJPair
        selection1 = self.mask1.Selection()
        selection2 = self.mask2.Selection()
        if sum(selection1) == 0 or sum(selection2) == 0:
            Action.stderr.write('Skipping empty masks in changeLJ14Pair')
            return None
        # Make sure we only selected 1 atom type, and figure out what it is
        attype1 = None
        attype2 = None
        for i in range(self.parm.ptr('natom')):
            if selection1[i] == 1:
                if not attype1:
                    attype1 = self.parm.parm_data['ATOM_TYPE_INDEX'][i]
                else:
                    if attype1 != self.parm.parm_data['ATOM_TYPE_INDEX'][i]:
                        raise ChangeLJPairError(
                                      'First mask matches multiple atom types!')
            if selection2[i] == 1:
                if not attype2:
                    attype2 = self.parm.parm_data['ATOM_TYPE_INDEX'][i]
                else:
                    if attype2 != self.parm.parm_data['ATOM_TYPE_INDEX'][i]:
                        raise ChangeLJPairError(
                                     'Second mask matches multiple atom types!')

        # Adjust 1-4 non-bonded terms as well if we're using a chamber-prmtop
        ChLJPair(self.parm, attype1, attype2, self.rmin, self.eps, True)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class checkvalidity(Action):
    """
    Basic checks for prmtop validity.
    """
    output = sys.stdout
    supported_classes = ('AmberParm', 'ChamberParm')

    def init(self, arg_list):
        pass

    def __str__(self):
        return 'Determining validity of prmtop'

    def execute(self):
        from ParmedTools.checkvalidity import check_validity
        from ParmedTools.exceptions import WarningList
        # Clear our warnings and start logging them, since check_validity
        # reports concerns about the prmtop through the warning system.
        warning_log = WarningList(empty_msg=('%s looks OK to me!' % self.parm))
        check_validity(self.parm, warning_log)
        warning_log.dump(self.output)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class change(Action):
    """
    Changes the property of given atoms to a new value. <property> can be
    CHARGE, MASS, RADII, SCREEN, ATOM_NAME, AMBER_ATOM_TYPE, ATOM_TYPE_INDEX,
    or ATOMIC_NUMBER (note, changing elements with this command will NOT
    change their assignment for SHAKE!).
   
    If given, the [quiet] keyword will prevent ParmEd from printing out a
    summary with every property changed for every atom that was changed
    useful for suppressing overwhelming output if you are zeroing every
    charge, for instance)
    """
    def init(self, arg_list):
        self.quiet = arg_list.has_key('quiet')
        self.mask = AmberMask(self.parm, arg_list.get_next_mask())
        self.prop = arg_list.get_next_string().upper()
        if self.prop in ('CHARGE', 'RADII', 'SCREEN', 'MASS'):
            self.new_val = arg_list.get_next_float()
            self.new_val_str = '%.4f' % self.new_val
        elif self.prop in ('ATOM_TYPE_INDEX', 'ATOMIC_NUMBER'):
            self.new_val = arg_list.get_next_int()
            self.new_val_str = '%4i' % self.new_val
        elif self.prop in ('ATOM_NAME', 'AMBER_ATOM_TYPE'):
            self.new_val = arg_list.get_next_string()
            if len(self.new_val) > 4:
                warnings.warn('Only 4 letters allowed for %s entries!'
                              'Truncating remaining letters.' % self.prop,
                              ParmWarning)
                self.new_val = self.new_val[:4]
            self.new_val_str = '%-4s' % self.new_val
        else:
            raise ParmedChangeError(
                        'You may only use "change" with CHARGE, MASS, RADII, '
                        'SCREEN, ATOM_NAME, AMBER_ATOM_TYPE, ATOM_TYPE_INDEX, '
                        'or ATOMIC_NUMBER!')

    def __str__(self):
        atnums = self.mask.Selection()
        if sum(atnums) == 0:
            return "change %s: Nothing to do" % self.prop
        if self.quiet:
            return "Changing %s of %s to %s" % (self.prop, self.mask,
                                                self.new_val)
        string = '\n'
        for i in range(self.parm.ptr('natom')):
            if atnums[i] == 1:
                string += "Changing %s of atom # %d (%s) from %s to %s\n" % (
                        self.prop, i+1, self.parm.parm_data['ATOM_NAME'][i],
                        self.parm.parm_data[self.prop][i], self.new_val_str)
        return string

    def execute(self):
        atnums = self.mask.Selection()
        if sum(atnums) == 0:
            warnings.warn('change %s: %s matches no atoms' %
                          (self.prop,self.mask), ParmWarning)
            return
        for i in range(len(atnums)):
            if atnums[i] == 1:
                self.parm.parm_data[self.prop][i] = self.new_val
        # Update the atom properties in the atom list
        try:
            self.parm.atom_list.refresh_data()
        except AttributeError:
            pass

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printinfo(Action):
    """
    Prints all prmtop data corresponding to the given %FLAG
    """
    outfile = sys.stdout
    def init(self, arg_list):
        self.flag = arg_list.get_next_string().upper()
        if not self.flag in self.parm.flag_list:
            warnings.warn('%%FLAG %s not found!' % self.flag,
                          SeriousParmWarning)
            self.found = False
        else:
            if self.parm.formats[self.flag].type is float:
                self.format = '%16.5f '
            else:
                self.format = '%-16s '

            self.found = True

    def __str__(self):
        ret_str = ''
        for i, item in enumerate(self.parm.parm_data[self.flag]):
            ret_str += self.format % item
            if i % 5 == 4:
                ret_str += '\n'

        return ret_str

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class addljtype(Action):
    """
    Turns given mask into a new LJ atom type. It uses the radius and Rmin from
    the first atom type in <mask> if new_radius or new_epsilon aren't provided
    """
    supported_classes = ('AmberParm', 'ChamberParm')
    def init(self, arg_list):
        self.new_radius_14 = arg_list.get_key_float('radius_14', None)
        self.new_epsilon_14 = arg_list.get_key_float('epsilon_14', None)
        self.new_radius = arg_list.get_key_float('radius', None)
        self.new_epsilon = arg_list.get_key_float('epsilon', None)
        self.mask = AmberMask(self.parm, arg_list.get_next_mask())

    def __str__(self):
        return 'Making atoms %s into a new LJ atom type' % self.mask

    def execute(self):
        from ParmedTools.addljtype import AddLJType
        # Find the first atom that's selected in this selection. We've
        # already made sure that at least one atom was selected
        sel_atms = self.mask.Selection()
        for i in range(self.parm.ptr('natom')):
            if sel_atms[i] == 1: 
                first_atm = i
                break
        # If either the radius or epsilon were not specified, then pull it
        # from the *first* atom in the selection
        if self.new_radius is None:
            self.new_radius = self.parm.LJ_radius[
                    self.parm.parm_data['ATOM_TYPE_INDEX'][first_atm]-1]
        else:
            self.new_radius = self.new_radius
        if self.new_epsilon is None:
            self.new_epsilon = self.parm.LJ_depth[
                    self.parm.parm_data['ATOM_TYPE_INDEX'][first_atm]-1]
        else:
            self.new_epsilon = self.new_epsilon
        # Now do the same for chamber prmtops
        if self.new_radius_14 is None and self.parm.chamber:
            self.new_radius_14 = self.parm.LJ_14_radius[
                    self.parm.parm_data['ATOM_TYPE_INDEX'][first_atm]-1]
        elif not self.parm.chamber:
            self.new_radius_14 = None
        if self.new_epsilon_14 is None and self.parm.chamber:
            self.new_epsilon_14 = self.parm.LJ_14_depth[
                    self.parm.parm_data['ATOM_TYPE_INDEX'][first_atm]-1]
        elif not self.parm.chamber:
            self.new_epsilon_14 = None
      
        AddLJType(self.parm, sel_atms, self.new_radius, self.new_epsilon, 
                  self.new_radius_14, self.new_epsilon_14)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class outparm(parmout):
    """
    Prints a new prmtop like parmout, but keeps its place in the action stack so
    several can be written out in 1 parmed session
    """

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printljtypes(Action):
    """
    Prints the Lennard Jones type index for the given atom mask or, if no value
    is given, the LJ type index for each atom.
    """
    supported_classes = ('AmberParm', 'ChamberParm')
    def init(self, arg_list):
        # Compile the list of indices
        try:
            self.mask = AmberMask(self.parm, arg_list.get_next_mask())
            self.type_list = None
            try_type = False
        except NoArgument:
            try_type = True
        if try_type:
            try:
                self.type_list = arg_list.get_next_string()
                self.mask = None
            except NoArgument:
                # Our default is print all indices...
                self.mask = AmberMask(self.parm, ':*')

        if self.mask is None:
            self.type_list = []
            type_fields = self.type_list.strip().split(',')
            for field in type_fields:
                if len(field.strip()) == 0: continue
                if '-' in field:
                    begin = int(field.split('-')[0])
                    end = min(int(field.split('-')[1]), self.parm.ptr('ntypes'))
                    if begin < 0 or end < begin: 
                        raise ParmError('printLJTypes: Bad atom type range')
                    self.type_list.extend([i for i in range(begin, end+1)])
                else:
                    self.type_list.append(int(field))

    def __str__(self):
        # Construct the atom selections and related atom types lists
        if self.mask:
            selection = self.mask.Selection()
        elif self.type_list:
            selection = [0 for i in range(self.parm.ptr('natom'))]
            for item in self.type_list:
                selection[item-1] = 1
        else:
            return 'Nothing to do for printLJTypes'

        self.idx = []

        for i in range(self.parm.ptr('natom')):
            if selection[i] == 1:
                if not self.parm.parm_data['ATOM_TYPE_INDEX'][i] in self.idx:
                    self.idx.append(self.parm.parm_data['ATOM_TYPE_INDEX'][i])

        string = '\n%15s %4s %4s\n' % ("  ATOM NUMBER  ", 'NAME', 'TYPE')
        string += '---------------------------------------------\n'
        for i in range(self.parm.ptr('natom')):
            if self.parm.parm_data['ATOM_TYPE_INDEX'][i] in self.idx:
                string += 'ATOM %-10d %-4s %-4s: Type index: %d\n' % (
                                      i+1, self.parm.parm_data['ATOM_NAME'][i],
                                      self.parm.parm_data['AMBER_ATOM_TYPE'][i],
                                      self.parm.parm_data['ATOM_TYPE_INDEX'][i])

        return string

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class scee(Action):
    """
    Sets the 1-4 EEL scaling factor in the prmtop
    """
    supported_classes = ('AmberParm', 'ChamberParm')
    def init(self, arg_list):
        self.scee_value = arg_list.get_next_float()

    def __str__(self):
        return ("Setting default SCEE electrostatic scaling value to %.4f" % 
                self.scee_value)

    def execute(self):
        if not 'SCEE_SCALE_FACTOR' in self.parm.flag_list:
            self.parm.addFlag('SCEE_SCALE_FACTOR', '5E16.8',
                              data=[self.scee_value
                                    for i in range(self.parm.ptr('nptra'))]
            )
        else:
            self.parm.parm_data['SCEE_SCALE_FACTOR'] = [self.scee_value 
                                for i in range(self.parm.ptr('nptra'))]
        # Now add it to each of the torsions
        for dih in self.parm.dihedrals_inc_h:
            dih.dihed_type.scee = self.scee_value
        for dih in self.parm.dihedrals_without_h:
            dih.dihed_type.scee = self.scee_value

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class scnb(Action):
    """
    Sets the 1-4 VDW scaling factor in the prmtop
    """
    supported_classes = ('AmberParm', 'ChamberParm')
    def init(self, arg_list):
        self.scnb_value = arg_list.get_next_float()

    def __str__(self):
        return ("Setting default SCNB van der Waals scaling value to %.4f" %
                 self.scnb_value)

    def execute(self):
      if not 'SCNB_SCALE_FACTOR' in self.parm.flag_list:
         self.parm.addFlag('SCNB_SCALE_FACTOR','5E16.8', data=[self.scnb_value 
                                        for i in range(self.parm.ptr('nptra'))])
      else:
         self.parm.parm_data['SCNB_SCALE_FACTOR'] = [self.scnb_value 
                                         for i in range(self.parm.ptr('nptra'))]
      # Now add it to each of the torsions
      for dih in self.parm.dihedrals_inc_h:
         dih.dihed_type.scnb = self.scnb_value
      for dih in self.parm.dihedrals_without_h:
         dih.dihed_type.scnb = self.scnb_value

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class changeljsingletype(Action):
    """
    Allows you to change the radius/well depth of a single LJ type specified by
    <mask>
    """
    supported_classes = ('AmberParm', 'ChamberParm')
    def init(self, arg_list):
        self.mask = AmberMask(self.parm, arg_list.get_next_mask())
        self.radius = arg_list.get_next_float()
        self.depth = arg_list.get_next_float()
        self.orig_radius = self.orig_depth = None
        if 1 in self.mask.Selection():
            first_loc = self.mask.Selection().index(1)
            attype = self.parm.parm_data['ATOM_TYPE_INDEX'][first_loc]
            self.orig_radius = self.parm.LJ_radius[attype-1]
            self.orig_depth = self.parm.LJ_depth[attype-1]

    def __str__(self):
        if sum(self.mask.Selection()) == 0:
            return "No atoms selected in %s. Nothing to do." % self.mask
        return ("Changing %s Lennard-Jones well depth from %.4f to %.4f "
                "(kal/mol) and radius from %.4f to %.4f (Angstroms)" %
                (self.mask, self.orig_depth, self.depth,
                 self.orig_radius, self.radius)
        )

    def execute(self):
        from math import sqrt
        from ParmedTools.exceptions import LJ_TypeError
        # If this is an empty mask do nothing
        if self.orig_radius is None: return
        # Make sure we've only selected a single atom type with our mask
        attype = None
        for i, sel in enumerate(self.mask.Selection()):
            if sel == 1:
                if attype is None:
                    attype = self.parm.parm_data['ATOM_TYPE_INDEX'][i]
                else:
                    if attype != self.parm.parm_data['ATOM_TYPE_INDEX'][i]:
                        raise LJ_TypeError('changeLJSingleType: Selection '
                                           'mask has multiple atom types!')
        # Fill the Lennard-Jones radius and depth arrays to make sure they're
        # up-to-date
        self.parm.fill_LJ()
        self.parm.LJ_radius[attype-1] = self.radius
        self.parm.LJ_depth[attype-1] = self.depth

        for i in range(self.parm.ptr('ntypes')):
            lj_index = self.parm.parm_data['NONBONDED_PARM_INDEX'][
                                    self.parm.ptr('ntypes')*i+attype-1] - 1
            rij = self.parm.LJ_radius[i] + self.radius
            wij = sqrt(self.parm.LJ_depth[i] * self.depth)
            acoef = wij * rij ** 12
            bcoef = 2 * wij * rij ** 6
            self.parm.parm_data['LENNARD_JONES_ACOEF'][lj_index] = acoef
            self.parm.parm_data['LENNARD_JONES_BCOEF'][lj_index] = bcoef

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printdetails(Action):
    """
    Returns information about all atoms in a given mask
    """
    def init(self, arg_list):
        self.mask = AmberMask(self.parm, arg_list.get_next_mask())

    def __str__(self):
        selection = self.mask.Selection()
        retstr = "\nThe mask %s matches %d atoms:\n\n" % (
                        self.mask, sum(selection))
        retstr += ("%7s%7s%9s%6s%6s%12s%12s%10s%10s%10s%10s\n" %
                   ('ATOM', 'RES', 'RESNAME', 'NAME', 'TYPE', 'LJ Radius',
                    'LJ Depth', 'Mass', 'Charge','GB Radius','GB Screen')
        )
        for i, atm in enumerate(self.parm.atom_list):
            if selection[i] == 1:
                retstr += (
                        "%7d%7d%9s%6s%6s%12.4f%12.4f%10.4f%10.4f%10.4f%10.4f\n"
                        % (i+1, atm.residue.idx, atm.residue.resname,
                           atm.atname, atm.attype,
                           self.parm.LJ_radius[atm.nb_idx-1], 
                           self.parm.LJ_depth[atm.nb_idx-1], atm.mass,
                           atm.charge, atm.radii, atm.screen)
                )
        return retstr

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printflags(Action):
    """
    Prints all %FLAGs found in the topology file
    """
    def init(self, arg_list):
        pass

    def __str__(self):
        string = '\n'
        for flag in self.parm.flag_list:
            string += '%%FLAG %s\n' % flag
        return string

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printpointers(Action):
    """
    Prints a list of all the POINTERS and their values
    """
    def init(self, arg_list):
        pass

    def __str__(self):
        ptrs = self.parm.parm_data['POINTERS']
        ret_str = """\nNATOM (number of atoms in system)................= %d
NTYPES (number of atom type names)...............= %d
NBONH (number of bonds containing H).............= %d
MBONA (number of bonds without H)................= %d
NTHETH (number of angles containing H)...........= %d
MTHETA (number of angles without H)..............= %d
NPHIH (number of dihedrals containing H).........= %d
MPHIA (number of dihedrals without H)............= %d
NHPARM (currently unused)........................= %d
NPARM (1 if made with addles, 0 if not)..........= %d
NNB (number of excluded atoms)...................= %d
NRES (number of residues in system)..............= %d
NBONA (MBONA + constraint bonds).................= %d
NTHETA (MTHETA + constraint angles)..............= %d
NPHIA (MPHIA + constraint dihedrals).............= %d
NUMBND (number of unique bond types).............= %d
NUMANG (number of unique angle types)............= %d
NPTRA (number of unique dihedral types)..........= %d
NATYP (number of nonbonded atom types)...........= %d
NPHB (number of distinct 10-12 H-bond pairs).....= %d
IFPERT (1 if prmtop is perturbed; not used)......= %d
NBPER (perturbed bonds; not used)................= %d
NGPER (perturbed angles; not used)...............= %d
NDPER (perturbed dihedrals; not used)............= %d
MBPER (bonds in perturbed group; not used).......= %d
MGPER (angles in perturbed group; not used)......= %d
MDPER (diheds in perturbed group; not used)......= %d
IFBOX (Type of box: 1=orthogonal, 2=not, 0=none).= %d
NMXRS (number of atoms in largest residue).......= %d
IFCAP (1 if solvent cap exists)..................= %d
NUMEXTRA (number of extra points in topology)....= %d
""" % (ptrs[0], ptrs[1], ptrs[2], ptrs[3], ptrs[4], ptrs[5], 
       ptrs[6], ptrs[7], ptrs[8], ptrs[9], ptrs[10], ptrs[11], 
       ptrs[12], ptrs[13], ptrs[14], ptrs[15], ptrs[16], 
       ptrs[17], ptrs[18], ptrs[19], ptrs[20], ptrs[21], 
       ptrs[22], ptrs[23], ptrs[24], ptrs[25], ptrs[26], 
       ptrs[27], ptrs[28], ptrs[29], ptrs[30])
        if len(ptrs) == 32: ret_str += \
            "NCOPY (number of PIMD slices/number of beads)....= %d\n" % ptrs[31]
        if self.parm.ptr('IFBOX'):
            ret_str += "\nSOLVENT POINTERS\n" + """
IPTRES (Final solute residue)....................= %d
NSPM (Total number of molecules).................= %d
NSPSOL (The first solvent "molecule")............= %d
""" % (self.parm.parm_data['SOLVENT_POINTERS'][0],
       self.parm.parm_data['SOLVENT_POINTERS'][1],
       self.parm.parm_data['SOLVENT_POINTERS'][2])
         
        return ret_str

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class setmolecules(Action):
    """
    Determines the molecularity of the system based on the bonding network and
    correctly determines the SOLVENT_POINTERS and ATOMS_PER_MOLECULE sections of
    the topology file. It will consider the ions to be part of the solute if
    True is passed or not if False is passed. Defaults to False.
    """
    def init(self, arg_list):
        # Allow solute_ions to be a keyword argument or the only argument.
        # Default is True
        solute_ions = arg_list.get_key_string('solute_ions', None)
        if solute_ions is None:
            solute_ions = arg_list.get_next_string(optional=True)
        if solute_ions is None:
            solute_ions = True
        elif solute_ions.lower() == 'true':
            self.solute_ions = True
        elif solute_ions.lower() == 'false':
            self.solute_ions = False
        else:
            warnings.warn("Value of solute_ions is unrecognized [%s]! "
                          "Assuming True" % solute_ions, SeriousParmWarning)
            self.solute_ions = True

    def __str__(self):
        return ("Setting MOLECULE properties of the prmtop (SOLVENT_POINTERS "
                "and ATOMS_PER_MOLECULE)")

    def execute(self):
        owner = self.parm.rediscover_molecules(self.solute_ions)
        if owner is not None:
            if not hasattr(self.parm, 'rst7'):
                warnings.warn(
                        'The atoms in %s were reordered to correct molecule '
                        'ordering. Any topology printed from now on will _not_ '
                        'work with the original inpcrd or trajectory files '
                        'created with this prmtop! Consider quitting and '
                        'loading a restart prior to using setMolecules' %
                        self.parm, ParmWarning
                )
   
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class combinemolecules(Action):
    """
    Combines the molecule <mol_id1> with the adjacent molecule <mol_id2> to
    ensure that those two molecules are imaged together. Most commonly used with
    protein/small ligand systems to keep them together during wrapping.  This
    will slightly affect the pressure calculation, but not by very much for
    small ligands. Note that <mol_id1> and <mol_id2> must be sequential if
    <mol_id2> is supplied
    """
    def init(self, arg_list):
        self.mol_id = arg_list.get_next_int()
        mol_id2 = arg_list.get_next_int(optional=True)
        if mol_id2 is not None:
            if self.mol_id + 1 != mol_id2:
                raise ParmedMoleculeError('Can only combine adjacent '
                                          'molecules!')

    def __str__(self):
        return "Combining molecules %d and %d into the same molecule" % (
                                    self.mol_id, self.mol_id + 1)

    def execute(self):
        from ParmedTools.mod_molecules import combineMolecules as cm
        cm(self.parm, self.mol_id)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
#class addcoarsegrain(Action):
#    """
#    Adds coarse graining information to the Amber topology file according to
#    a given coarse graining parameter file
#    """
#    def init(self, arg_list):
#        self.cg_param_file = arg_list.get_next_string()
#        if not os.path.exists(self.cg_param_file):
#            raise CoarseGrainError('Cannot find parameter file %s' % 
#                                    self.cg_param_file)
#        # Check to see if we've already made this a coarsegrained file...
#        if 'ANGLE_COEF_A' in self.parm.flag_list:
#            warnings.warn('Prmtop already has coarse grained sections',
#                        ParmWarning)
#
#    def __str__(self):
#        return ("Setting up coarse graining for topology file using parameter "
#                "file " + self.cg_param_file)
#
#    def execute(self):
#        from ParmedTools.coarsegrain import addCoarseGrain as addCG
#        if 'ANGLE_COEF_A' in self.parm.flag_list: return
#        addCG(self.parm, self.cg_param_file)
#
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class changeprotstate(Action):
    """
    Changes the protonation state of a given titratable residue that can be
    treated via constant pH MD in Amber. This corresponds to all residues found
    in $AMBERHOME/AmberTools/src/etc/cpin_data.py
    """
    supported_classes = ('AmberParm',)
    def init(self, arg_list):
        self.state = arg_list.get_next_int()
        self.mask = AmberMask(self.parm, arg_list.get_next_mask())

    def __str__(self):
        sel = self.mask.Selection()
        if sum(sel) == 0:
            return "No residues selected for state change"
        res = self.atom_list[self.index(1)].residue
        return 'Changing protonation state of residue %d (%s) to %d' % (res.idx,
                            res.resname, self.state)
   
    @staticmethod
    def _add_ash_glh(residues):
        """
        Adds ASH and GLH to the titratable residue list unless it's already
        there
        """
        if 'ASH' in residues.titratable_residues: return None
        ash = residues.TitratableResidue('ASH',
                    ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG',
                     'OD1', 'OD2', 'HD21', 'C', 'O'], pka=4.0)
        ash.add_state(protcnt=0, refene=None,
                      charges=[-0.4157, 0.2719, 0.0341, 0.0864, -0.1783,
                               -0.0122, -0.0122, 0.7994, -0.8014, -0.8014, 0.0,
                               0.5973, -0.5679]
        )
        ash.add_state(protcnt=1, refene=None,
                      charges=[-0.4157, 0.2719, 0.0341, 0.0864, -0.0316, 0.0488,
                               0.0488, 0.6462, -0.5554, -0.6376, 0.4747, 0.5973,
                               -0.5679]
        )

        glh = residues.TitratableResidue('GLH',
                    ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2',
                     'HG3', 'CD', 'OE1', 'OE2', 'HE21', 'C', 'O'], pka=4.4)
        glh.add_state(protcnt=0, refene=None,
                      charges=[-0.4157, 0.2719, 0.0145, 0.0779, -0.0398,
                               -0.0173, -0.0173, 0.0136, -0.0425, -0.0425,
                               0.8054, -0.8188, -0.8188, 0.0, 0.5973, -0.5679]
        )
        glh.add_state(protcnt=1, refene=None,
                      charges=[-0.4157, 0.2719, 0.0145, 0.0779, -0.0071, 0.0256,
                               0.0256, -0.0174, 0.0430, 0.0430, 0.6801, -0.5838,
                               -0.6511, 0.4641, 0.5973, -0.5679]
        )
        residues.ASH, residues.GLH = ash, glh
        residues.titratable_residues.extend(['ASH', 'GLH'])

    def execute(self):
        from cpinutils import residues
        changeprotstate._add_ash_glh(residues)
        sel = self.mask.Selection()
        # If we didn't select any residues, just return
        if sum(sel) == 0: return
        res = self.parm.atom_list[sel.index(1)].residue
        resnum = res.idx
        resname = res.resname
        # Get the charges from cpin_data. The first 2 elements are energy and 
        # proton count so the charges are chgs[2:]
        if not resname in residues.titratable_residues:
            raise ChangeStateError("Residue %s isn't defined as a titratable "
                                   "residue in cpin_data.py" % resname)

        res = getattr(residues, resname)

        if self.state >= len(res.states):
            raise ChangeStateError('Residue %s only has titratable states '
                                   '0--%d. You chose state %d' %
                                   (resname, len(res.states)-1, self.state))

        if sum(sel) != len(res.states[self.state].charges):
            raise ChangeStateError('You must select one and only one entire '
                                   'titratable residue')
      
        chgnum = 0
        for i in range(self.parm.parm_data['RESIDUE_POINTER'][resnum-1]-1,
                       self.parm.parm_data['RESIDUE_POINTER'][resnum]-1):
            if sel[i] != 1:
                raise ChangeStateError('You must select 1 and only 1 entire '
                                       'residue to change the protonation '
                                       'state of')
            # Actually make the change
            self.parm.parm_data['CHARGE'][i] = \
                        res.states[self.state].charges[chgnum]
            chgnum += 1

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class netcharge(Action):
    """
    Prints the total charge of all of the atoms given by the mask. Defaults to
    all atoms
    """
    supported_classes = ('AmberParm', 'ChamberParm')
    outfile = sys.stdout
    def init(self, arg_list):
        mask = arg_list.get_next_mask(optional=True)
        if mask is None: mask = ':*'
        self.mask = AmberMask(self.parm, mask)

    def execute(self):
        """ Calculates the charge of all atoms selected in mask """
        sel = self.mask.Selection()

        netchg = 0.0
        for i in range(len(sel)):
            if sel[i]: netchg += self.parm.parm_data['CHARGE'][i]

        self.outfile.write('  The net charge of %s is %.4f\n' % 
                           (self.mask, netchg))
        return netchg

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class strip(Action):
    """
    Deletes the atoms specified by <mask> from the topology file and rebuilds
    the topology file according to the parameters that remain.
    """
    def init(self, arg_list):
        self.mask = AmberMask(self.parm, arg_list.get_next_mask())
        self.num_atms = sum(self.mask.Selection())

    def __str__(self):
        return "Removing mask '%s' (%d atoms) from the topology file." % (
                                    self.mask, self.num_atms)

    def execute(self):
        self.parm.delete_mask(self.mask)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class definesolvent(Action):
    """
    Allows you to change what parmed.py will consider to be "solvent". 
    <residue list> must be a comma-separated set of residue names with no
    spaces between them.
    """
    def init(self, arg_list):
        res_list = arg_list.get_next_string()
        res_list.replace(' ', '')
        if res_list.endswith(','): self.res_list = res_list[:len(res_list)-1]
        else: self.res_list = res_list
        self.parm.solvent_residues = res_list.split(',')

    def __str__(self):
        return "Residues %s are now considered to be solvent" % self.res_list

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class addexclusions(Action):
    """ 
    Allows you to add arbitrary exclusions to the exclusion list. Every atom in
    <mask2> is added to the exclusion list for each atom in <mask1> so that
    non-bonded interactions between those atom pairs will not be computed. NOTE
    that this ONLY applies to direct-space (short-range) non-bonded potentials.
    For PME simulations, long-range electrostatics between these atom pairs are
    still computed (in different unit cells).
    """
    def init(self, arg_list):
        self.mask1 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask2 = AmberMask(self.parm, arg_list.get_next_mask())

    def __str__(self):
        return 'Adding atoms from %s to exclusion lists of atoms in %s' % (
            self.mask2, self.mask1)
   
    def execute(self):
        sel1 = self.mask1.Selection()
        sel2 = self.mask2.Selection()
        # Loop through both selections and add each selected atom in sel2 to
        # the exclusion list for selected atoms in sel1 (and vice-versa).
        for i in range(len(sel1)):
            if sel1[i]:
                for j in range(len(sel2)):
                    if sel2[j]:
                        # Make sure that this atom isn't already in the
                        # exclusion list by virtue of being a bonded partner.
                        atm1 = self.parm.atom_list[i]
                        atm2 = self.parm.atom_list[j]
                        # Skip over atm1 == atm2
                        if atm1 is atm2: continue
                        # Add each other to each other's exclusion lists.
                        atm1.exclude(atm2)
                        self.parm.atom_list.changed = True

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printbonds(Action):
    """
    Prints all of the bonds (with their details) for the given atoms in the
    mask
    """
    def init(self, arg_list):
        self.mask = AmberMask(self.parm, arg_list.get_next_mask())

    def __str__(self):
        retstr = '%-20s %-20s %-10s %-10s\n' % (
                                    'Atom 1', 'Atom 2', 'R eq', 'Frc Cnst')
        # Loop through all of the bonds without and inc hydrogen
        atomsel = self.mask.Selection()
        for bond in self.parm.bonds_without_h:
            idx1 = bond.atom1.starting_index
            idx2 = bond.atom2.starting_index
            if not (atomsel[idx1] or atomsel[idx2]): continue

            atm1 = self.parm.atom_list[idx1]
            atm2 = self.parm.atom_list[idx2]
            retstr += '%7d %4s (%4s) %7d %4s (%4s) %10.4f %10.4f\n' % (
                    idx1+1, atm1.atname, atm1.attype, idx2+1, atm2.atname,
                    atm2.attype, bond.bond_type.req, bond.bond_type.k)

        for bond in self.parm.bonds_inc_h:
            idx1 = bond.atom1.starting_index
            idx2 = bond.atom2.starting_index
            if not (atomsel[idx1] or atomsel[idx2]): continue

            atm1 = self.parm.atom_list[idx1]
            atm2 = self.parm.atom_list[idx2]
            retstr += '%7d %4s (%4s)  %7d %4s (%4s) %10.4f %10.4f\n' % (
                        idx1+1, atm1.atname, atm1.attype, idx2+1, atm2.atname,
                        atm2.attype, bond.bond_type.req, bond.bond_type.k)

        return retstr

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printangles(Action):
    """
    Prints all of the angles (with their details) for the given atoms in the
    mask
    """
    def init(self, arg_list):
        self.mask = AmberMask(self.parm, arg_list.get_next_mask())

    def __str__(self):
        retstr = '%-20s %-20s %-20s %-10s %-10s\n' % (
                        'Atom 1', 'Atom 2', 'Atom 3', 'Frc Cnst', 'Theta eq')
        # Loop through all of the bonds without and inc hydrogen
        atomsel = self.mask.Selection()
        for angle in self.parm.angles_without_h:
            idx1 = angle.atom1.starting_index
            idx2 = angle.atom2.starting_index
            idx3 = angle.atom3.starting_index
            if not (atomsel[idx1] or atomsel[idx2] or atomsel[idx3]): continue

            atm1 = self.parm.atom_list[idx1]
            atm2 = self.parm.atom_list[idx2]
            atm3 = self.parm.atom_list[idx3]
            retstr += ('%7d %4s (%4s)  %7d %4s (%4s)  %7d %4s (%4s) '
                       '%10.4f %10.4f\n' % (idx1+1, atm1.atname, atm1.attype,
                       idx2+1, atm2.atname, atm2.attype, idx3+1, atm3.atname,
                       atm3.attype, angle.angle_type.k,
                       angle.angle_type.theteq*180/math.pi)
            )

        for angle in self.parm.angles_inc_h:
            idx1 = angle.atom1.starting_index
            idx2 = angle.atom2.starting_index
            idx3 = angle.atom3.starting_index
            if not (atomsel[idx1] or atomsel[idx2] or atomsel[idx3]): continue

            atm1 = self.parm.atom_list[idx1]
            atm2 = self.parm.atom_list[idx2]
            atm3 = self.parm.atom_list[idx3]
            retstr += ('%7d %4s (%4s)  %7d %4s (%4s)  %7d %4s (%4s) '
                       '%10.4f %10.4f\n' % (idx1+1, atm1.atname, atm1.attype,
                       idx2+1, atm2.atname, atm2.attype, idx3+1, atm3.atname,
                       atm3.attype, angle.angle_type.k,
                       angle.angle_type.theteq*180/math.pi)
            )
      
        return retstr

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printdihedrals(Action):
    """
    Prints all of the dihedrals (with their details) for the given atoms in the
    mask
    """
    def init(self, arg_list):
        self.mask = AmberMask(self.parm, arg_list.get_next_mask())

    def __str__(self):
        retstr = '%-20s %-20s %-20s %-20s  %-10s %-10s %-10s %-10s %-10s\n' % (
                'Atom 1', 'Atom 2', 'Atom 3', 'Atom 4', 'Height', 'Periodic.',
                'Phase', 'EEL Scale', 'VDW Scale')
        # Loop through all of the bonds without and inc hydrogen
        atomsel = self.mask.Selection()
        for dihedral in self.parm.dihedrals_without_h:
            idx1 = dihedral.atom1.starting_index
            idx2 = dihedral.atom2.starting_index
            idx3 = dihedral.atom3.starting_index
            idx4 = dihedral.atom4.starting_index
            if not (atomsel[idx1] or atomsel[idx2] or atomsel[idx3] or 
                    atomsel[idx4]):
                continue

            atm1 = self.parm.atom_list[idx1]
            atm2 = self.parm.atom_list[idx2]
            atm3 = self.parm.atom_list[idx3]
            atm4 = self.parm.atom_list[idx4]
            # Determine if it's an Improper, Multiterm, or neither
            if dihedral.signs[1] < 0:
                char = 'I'
            elif dihedral.signs[0] < 0:
                char = 'M'
            else:
                char = ' '
            retstr += ('%1s %7d %4s (%4s)  %7d %4s (%4s)  %7d %4s (%4s)  '
                       '%7d %4s (%4s) %10.4f %10.4f %10.4f %10.4f %10.4f\n' %
                       (char, idx1+1, atm1.atname, atm1.attype, idx2+1,
                        atm2.atname, atm2.attype, idx3+1, atm3.atname,
                        atm3.attype, idx4+1, atm4.atname, atm4.attype,
                        dihedral.dihed_type.phi_k, dihedral.dihed_type.per,
                        dihedral.dihed_type.phase*180/math.pi,
                        dihedral.dihed_type.scee, dihedral.dihed_type.scnb)
            )

        for dihedral in self.parm.dihedrals_inc_h:
            idx1 = dihedral.atom1.starting_index
            idx2 = dihedral.atom2.starting_index
            idx3 = dihedral.atom3.starting_index
            idx4 = dihedral.atom4.starting_index
            if not (atomsel[idx1] or atomsel[idx2] or atomsel[idx3] or 
                    atomsel[idx4]):
                continue

            atm1 = self.parm.atom_list[idx1]
            atm2 = self.parm.atom_list[idx2]
            atm3 = self.parm.atom_list[idx3]
            atm4 = self.parm.atom_list[idx4]
            if dihedral.signs[1] < 0:
                char = 'I'
            elif dihedral.signs[0] < 0:
                char = 'M'
            else:
                char = ' '
            retstr += ('%1s %7d %4s (%4s)  %7d %4s (%4s)  %7d %4s (%4s)  '
                       '%7d %4s (%4s) %10.4f %10.4f %10.4f %10.4f %10.4f\n' %
                      (char, idx1+1, atm1.atname, atm1.attype, idx2+1,
                       atm2.atname, atm2.attype, idx3+1, atm3.atname,
                       atm3.attype, idx4+1, atm4.atname, atm4.attype,
                       dihedral.dihed_type.phi_k, dihedral.dihed_type.per,
                       dihedral.dihed_type.phase*180/math.pi,
                       dihedral.dihed_type.scee, dihedral.dihed_type.scnb)
            )

        return retstr

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class setbond(Action):
    """
    Changes (or adds a non-existent) bond in the topology file. Each mask must
    select the same number of atoms, and a bond will be placed between the
    atoms in mask1 and mask2 (one bond between atom1 from mask1 and atom1
    from mask2 and another bond between atom2 from mask1 and atom2 from mask2,
    etc.)
    """
    supported_classes = ('AmberParm', 'ChamberParm')

    def init(self, arg_list):
        self.k = arg_list.get_next_float()
        self.req = arg_list.get_next_float()
        self.mask1 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask2 = AmberMask(self.parm, arg_list.get_next_mask())
   
    def __str__(self):
        return ('Set a bond between %s and %s with k = %f kcal/(mol Angstrom'
                '**2) and Req = %f Angstroms') % (self.mask1, self.mask2,
                self.k, self.req)

    def execute(self):
        sel1 = self.mask1.Selection()
        sel2 = self.mask2.Selection()

        if sum(sel1) != sum(sel2):
            raise SetParamError('setBond: Each mask must select the same '
                                'number of atoms!')

        # If no atoms, nothing to do
        if sum(sel1) == 0: return

        # Create the new bond type
        new_bnd_typ = BondType(self.k, self.req, -1)
        # Does that bond type exist in the list already? If it does, re-bind
        # new_bnd to that bond type reference
        exists = False
        for bnd_typ in self.parm.bond_type_list:
            if new_bnd_typ == bnd_typ:
                new_bnd_typ = bnd_typ
                exists = True
                break
        # If the bond is new, add it to the type list
        if not exists:
            self.parm.bond_type_list.append(new_bnd_typ)

        atnum1, atnum2 = -1, -1
        # Loop through all of the selected atoms
        for it in range(sum(sel1)):
            # Collect the atoms involved
            atnum1 = sel1.index(1, atnum1+1)
            atnum2 = sel2.index(1, atnum2+1)
            atm1 = self.parm.atom_list[atnum1]
            atm2 = self.parm.atom_list[atnum2]

            # See if any atom is Hydrogen (allow for deuteriums)
            if atm1.element == 1 or atm2.element == 1:
                bond_list = self.parm.bonds_inc_h
            else:
                bond_list = self.parm.bonds_without_h
   
            # See if the bond exists in the first place, and if so, replace its
            # bond type with our new bond type (new_bnd)
            if atm2 in atm1.bond_partners or atm1 in atm2.bond_partners:
                for bnd in bond_list:
                    if atm1 in bnd and atm2 in bnd:
                        bnd.bond_type = new_bnd_typ
                        bond_list.changed = True
                        break
   
            # Otherwise, it doesn't exist, so we just create a new one
            else:
                bond_list.append(Bond(atm1, atm2, new_bnd_typ))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class setangle(Action):
    """
    Changes (or adds a non-existent) angle in the topology file. Each mask must
    select the same number of atoms, and an angle will be placed between the
    atoms in mask1, mask2, and mask3 (one angle between atom1 from mask1, atom1
    from mask2, and atom1 from mask3, another angle between atom2 from mask1,
    atom2 from mask2, and atom2 from mask3, etc.)
    """
    supported_classes = ('AmberParm', 'ChamberParm')

    def init(self, arg_list):
        self.k = arg_list.get_next_float()
        self.theteq = arg_list.get_next_float() * math.pi / 180.0
        self.mask1 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask2 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask3 = AmberMask(self.parm, arg_list.get_next_mask())
   
    def __str__(self):
        return ('Set an angle between %s, %s and %s with k = %f kcal/(mol '
                'rad**2) and THETeq = %f degrees') % (self.mask1, self.mask2,
                self.mask3, self.k, self.theteq * 180/math.pi)

    def execute(self):
        sel1 = self.mask1.Selection()
        sel2 = self.mask2.Selection()
        sel3 = self.mask3.Selection()

        if sum(sel1) != sum(sel2) or sum(sel1) != sum(sel3):
            raise SetParamError('Each mask in setAngle must select the same '
                                'number of atoms!')

        if sum(sel1) == 0: return

        # Create the new angle type
        new_ang_typ = AngleType(self.k, self.theteq, -1)
        # Does that angle type exist in the list already? If it does, re-bind
        # new_ang to that angle type reference
        exists = False
        for ang_typ in self.parm.angle_type_list:
            if new_ang_typ == ang_typ:
                new_ang_typ = ang_typ
                exists = True
                break
        # If the angle is new, add it to the type list
        if not exists:
            self.parm.angle_type_list.append(new_ang_typ)

        atnum1, atnum2, atnum3 = -1, -1, -1

        # Loop through all of the selections
        for it in range(sum(sel1)):
            # Collect the atoms involved
            atnum1 = sel1.index(1, atnum1+1)
            atnum2 = sel2.index(1, atnum2+1)
            atnum3 = sel3.index(1, atnum3+1)
            atm1 = self.parm.atom_list[atnum1]
            atm2 = self.parm.atom_list[atnum2]
            atm3 = self.parm.atom_list[atnum3]
            # See if any atom is Hydrogen (allow for deuteriums)
            if atm1.element == 1 or atm2.element == 1 or atm3.element == 1:
                angle_list = self.parm.angles_inc_h
            else:
                angle_list = self.parm.angles_without_h
   
            # See if the angle exists in the first place, and if so, replace its
            # angle type with our new angle type (new_ang)
            if ((atm1 in atm2.bond_partners and atm1 in atm3.angle_partners) and
                (atm2 in atm3.bond_partners)):
                for ang in angle_list:
                    if atm1 in ang and atm2 in ang and atm3 in ang:
                        ang.angle_type = new_ang_typ
                        angle_list.changed = True
                        break
   
            # Otherwise, it doesn't exist, so we just create a new one
            else:
                angle_list.append(Angle(atm1, atm2, atm3, new_ang_typ))
   
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class adddihedral(Action):
    """
    Adds a dihedral between mask1, mask2, mask3, and mask4. Each mask must
    specify the same number of atoms, and the dihedral is defined around the
    bond between atoms in mask 2 and 3. If each mask selects 2 atoms, for
    instance, a dihedral will be placed around atom1 in mask 1, atom1 in mask 2,
    atom1 in mask 3, and atom1 in mask 4.  A second dihedral will be placed
    around atom2 in mask 1, atom2 in mask 2, atom2 in mask 3, and atom2 in
    mask4. dihed_type can either be "normal", "multiterm", or "improper". Note
    for ring systems of 6 or fewer atoms, you'll need to use "multiterm" to
    avoid double-counting 1-4 interactions for some of the dihedrals.  <phi_k>
    is the barrier height of the dihedral term, <per> is the periodicity,
    <phase> is the phase offset, <scee> is the 1-4 EEL scaling factor for this
    dihedral (default for AMBER is 1.2, default for GLYCAM is 1.0), and <scnb>
    is the 1-4 VDW scaling factor for this dihedral (default for AMBER is 2.0,
    default for GLYCAM is 1.0)
    """
    supported_classes = ('AmberParm', 'ChamberParm')

    def init(self, arg_list):
        self.phi_k = arg_list.get_next_float()
        self.per = arg_list.get_next_float()
        self.phase = arg_list.get_next_float() * math.pi/180.0
        self.scee = arg_list.get_next_float(optional=True, default=1.2)
        self.scnb = arg_list.get_next_float(optional=True, default=2.0)
        self.mask1 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask2 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask3 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask4 = AmberMask(self.parm, arg_list.get_next_mask())
        dihed_type = arg_list.get_key_string('type', 'normal')
        if dihed_type.lower() == 'normal'[:len(dihed_type)]:
            self.improper = False
            self.multiterm = False
            self.type = 'a normal'
        elif dihed_type.lower() == 'improper'[:len(dihed_type)]:
            self.improper = True
            self.multiterm = True
            self.type = 'an improper'
        elif dihed_type.lower() == 'multiterm'[:len(dihed_type)]:
            self.improper = False
            self.multiterm = True
            self.type = 'a multiterm'
   
    def __str__(self):
        return ('Set %s dihedral between %s, %s, %s, and %s with phi_k = %f '
                'kcal/mol periodicity = %f phase = %f degrees scee = %f '
                'scnb = %f' %
               (self.type, self.mask1, self.mask2, self.mask3, self.mask4,
                self.phi_k, self.per, self.phase * 180/math.pi, self.scee,
                self.scnb)
        )

    def execute(self):
        sel1 = self.mask1.Selection()
        sel2 = self.mask2.Selection()
        sel3 = self.mask3.Selection()
        sel4 = self.mask4.Selection()

        if (sum(sel1) != sum(sel2) or sum(sel1) != sum(sel3) or
                                      sum(sel1) != sum(sel4)):
            raise SetParamError('addDihedral: Each mask must select the same '
                                'number of atoms!')
      
        # If we do nothing, just return
        if sum(sel1) == 0: return
   
        # Create the new dihedral type
        new_dih_typ = DihedralType(self.phi_k, self.per, self.phase, self.scee,
                                   self.scnb, -1)
        self.parm.dihedral_type_list.append(new_dih_typ)

        # Loop through all of the atoms
        atnum1, atnum2, atnum3, atnum4 = -1, -1, -1, -1

        for it in range(sum(sel1)):
            # Collect the atoms involved
            atnum1 = sel1.index(1, atnum1+1)
            atnum2 = sel2.index(1, atnum2+1)
            atnum3 = sel3.index(1, atnum3+1)
            atnum4 = sel4.index(1, atnum4+1)
            atm1 = self.parm.atom_list[atnum1]
            atm2 = self.parm.atom_list[atnum2]
            atm3 = self.parm.atom_list[atnum3]
            atm4 = self.parm.atom_list[atnum4]
            if (atm1 == atm2 or atm1 == atm3 or atm1 == atm4 or
                atm2 == atm3 or atm2 == atm4 or atm3 == atm4):
                raise SetParamError('addDihedral: Duplicate atoms found!')
   
            # Make sure atom 1 doesn't occur in the 3rd or 4th spot since we'll
            # have the -0 effect biting us...
            if atm3.starting_index == 0 or atm4.starting_index == 0:
                atm1, atm2, atm3, atm4 = atm4, atm3, atm2, atm1
   
            # See if any atom is Hydrogen (allow for deuterium)
            if (atm1.element == 1 or atm2.element == 1 or 
                atm3.element == 1 or atm4.element == 1):
                dihed_list = self.parm.dihedrals_inc_h
            else:
                dihed_list = self.parm.dihedrals_without_h
   
            # See what signs has to be.  signs is a 2-element list with a 1 or
            # -1 for the 3rd or 4th atom index. A -1 for the 3rd means the end
            # groups are not calculated (multiterm, 6-or-lower-membered rings,
            # impropers, etc.) and a -1 for the 4th means improper (-1 for 4th
            # always means -1 for 3rd)
            if self.improper:
                signs = [-1,-1]
            elif self.multiterm:
                signs = [-1,1]
            else:
                signs = [1,1]
            # Create our new dihedral!
            dihed_list.append(
                    Dihedral(atm1, atm2, atm3, atm4, new_dih_typ, signs)
            )

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class addatomicnumber(Action):
    """
    Adds the atomic number of each atom to a new section titled "ATOMIC_NUMBER"
    in the topology file. Elements are identified by the atomic mass found in
    the MASS section of the topology files.  Elements are matched by picking the
    element whose average atomic mass in the periodic table is closest to each
    atom, which should work appropriately for all isotopes of all atoms, except
    possibly Tritium
    """
    def init(self, arg_list):
        self.present = 'ATOMIC_NUMBER' in self.parm.flag_list

    def __str__(self):
        if self.present:
            return 'ATOMIC_NUMBER already in [%s] -- Doing nothing.' % self.parm
        return "Adding ATOMIC_NUMBER to [%s]" % self.parm

    def execute(self):
        if self.present: return
        self.parm.addFlag('ATOMIC_NUMBER', '10I8',
                          num_items=self.parm.ptr('natom'))
        for i, atm in enumerate(self.parm.atom_list):
            self.parm.parm_data['ATOMIC_NUMBER'][i] = atm.atomic_number

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class deletedihedral(Action):
    """
    Deletes the dihedral around <mask2> and <mask3> in which the end-groups are
    <mask1> and <mask4>. For multi-term dihedrals, it removes each term.
    """
    supported_classes = ('AmberParm', 'ChamberParm')
    def init(self, arg_list):
        self.mask1 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask2 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask3 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask4 = AmberMask(self.parm, arg_list.get_next_mask())
        if (sum(self.mask1.Selection()) != sum(self.mask2.Selection()) or
            sum(self.mask1.Selection()) != sum(self.mask3.Selection()) or
            sum(self.mask1.Selection()) != sum(self.mask4.Selection())):
            raise DeleteDihedralError('All masks must select the same number '
                                      'of atoms!. They selected %d, %d, %d, '
                                      'and %d, respectively' % (
                    sum(self.mask1.Selection()), sum(self.mask2.Selection()),
                    sum(self.mask3.Selection()), sum(self.mask4.Selection()))
            )

    def __str__(self):
        if sum(self.mask1.Selection()) == 0:
            return 'No specified dihedrals to delete'
        return ('Deleting dihedral terms involving [%s]-[%s]-[%s]-[%s]'
                ' (At most %d total, distinct, dihedrals)' %
               (self.mask1, self.mask2, self.mask3, self.mask4,
                sum(self.mask1.Selection()))
        )

    def execute(self):
        sel1, sel2 = self.mask1.Selection(), self.mask2.Selection()
        sel3, sel4 = self.mask3.Selection(), self.mask4.Selection()
        # Bail out if we're deleting nothing
        if sum(sel1) == 0: return

        # Keep track of the dihedrals we want to delete from each
        # dihedral list (dihedrals_inc_h, dihedrals_without_h)
        deleting_dihedrals = [[],[]]
        # We have already checked that they are the same number of atoms
        # Now, loop through the atoms and see if any dihedrals match that spec
        atnum1 = atnum2 = atnum3 = atnum4 = -1
        total_diheds = 0
        for i in range(sum(sel1)):
            # Collect the atoms involved
            atnum1 = sel1.index(1, atnum1+1)
            atnum2 = sel2.index(1, atnum2+1)
            atnum3 = sel3.index(1, atnum3+1)
            atnum4 = sel4.index(1, atnum4+1)
            atm1 = self.parm.atom_list[atnum1]
            atm2 = self.parm.atom_list[atnum2]
            atm3 = self.parm.atom_list[atnum3]
            atm4 = self.parm.atom_list[atnum4]
            # Make sure none of the indices are the same
            if (atm1 == atm2 or atm1 == atm3 or atm1 == atm4 or 
                atm2 == atm3 or atm2 == atm4 or atm3 == atm4):
                warnings.warn('Skipping %d-%d-%d-%d dihedral deletion -- '
                              'duplicate atoms!' %
                              (atnum1, atnum2, atnum3, atnum4),
                              SeriousParmWarning)
                continue
            # This helps us keep track of multi-term dihedrals so we don't
            # confuse users
            found_this_dihedral = False
            # Figure out if our dihedral would have hydrogen or not (limits what
            # dihedral list we have to search...)
            if (atm1.element == 1 or atm2.element == 1 or
                atm3.element == 1 or atm4.element == 1):
                dihed_list = self.parm.dihedrals_inc_h
                dihed_list_idx = 0
            else:
                dihed_list = self.parm.dihedrals_without_h
                dihed_list_idx = 1
            # Now search through our dihedral list to see which indexes (if any)
            # we have to remove. Keep tabs of them so we can pop them in reverse
            # order (so we don't have to re-figure indices) afterwards
            proposed_dihedral = (atnum1, atnum2, atnum3, atnum4)
            for j, dihed in enumerate(dihed_list):
                if dihed == proposed_dihedral:
                    if not found_this_dihedral:
                        print 'Matched dihedral number %d' % j
                        found_this_dihedral = True
                        total_diheds += 1
                    else:
                        print '  Matched multi-term dihedral number %d' % j
                    deleting_dihedrals[dihed_list_idx].append(j)

        if not deleting_dihedrals[0] and not deleting_dihedrals[1]:
            print 'No dihedrals matched -- not deleting any dihedrals'
            return

        print 'Deleting %d dihedrals' % (len(deleting_dihedrals[0]) + 
                                         len(deleting_dihedrals[1])),
        print ' (%d distinct dihedrals)' % total_diheds

        # At this point, we've collected all of our dihedrals, now sort them
        if deleting_dihedrals[0]: deleting_dihedrals[0].sort()
        if deleting_dihedrals[1]: deleting_dihedrals[1].sort()
        # deleting_dihedrals now contains all of our dihedral indexes
        if deleting_dihedrals[0]:
            while deleting_dihedrals[0]:
                del self.parm.dihedrals_inc_h[deleting_dihedrals[0].pop()]
        if deleting_dihedrals[1]:
            while deleting_dihedrals[1]:
                del self.parm.dihedrals_without_h[deleting_dihedrals[1].pop()]

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printljmatrix(Action):
    """ 
    This function prints out how every atom type interacts with the atom type(s)
    in <mask>. The atom types are printed as all type names that have at least
    one atom with the Lennard Jones index given in square brackets at the end.
    Alternatively, you can request a particular atom type index
    """
    supported_classes = ('AmberParm', 'ChamberParm')
    def init(self, arg_list):
        self.idx = arg_list.get_next_int(optional=True)
        self.mask = None
        if self.idx is None:
            self.mask = AmberMask(self.parm, arg_list.get_next_mask())
   
    def __str__(self):
        ntypes = self.parm.ptr('NTYPES')
        ret_str = ''
        if self.idx is not None:
            sel = [0 for i in self.parm.parm_data['ATOM_TYPE_INDEX']]
            for i, idx in enumerate(self.parm.parm_data['ATOM_TYPE_INDEX']):
                if idx == self.idx: sel[i] = 1
        else:
            sel = self.mask.Selection()
        # If we selected no atoms, bail out
        if sum(sel) == 0: return ''
        # Figure out which types correspond to which names
        typenames = [set() for i in range(self.parm.ptr('NTYPES'))]
        for i, ty in enumerate(self.parm.parm_data['ATOM_TYPE_INDEX']):
            typenames[ty-1].add(self.parm.parm_data['AMBER_ATOM_TYPE'][i])
        # Otherwise, collect our list of atom types that we selected
        sel_types = set()
        for i, val in enumerate(sel):
            if not val: continue
            attype = self.parm.parm_data['ATOM_TYPE_INDEX'][i]
            sel_types.add(attype)
        sel_types = sorted(list(sel_types)) # sort the atom types
        # Convert all of the typenames into strings, then find the longest one so
        # we can properly format the string
        maxlen = 0
        for i, names in enumerate(typenames):
            typenames[i] = ','.join(sorted(list(names))) + ' [%d]' % (i+1)
            maxlen = max(maxlen, len(typenames[i]))
        ret_str = '\n%%%ds %%%ds %%15s %%15s %%10s %%10s' % (maxlen, maxlen) % (
                    'Atom Type 1', 'Atom Type 2', 'A coefficient',
                    'B coefficient', 'R i,j', 'Eps i,j')
        ret_str += '\n' + '-'*len(ret_str) + '\n'
        for ty in sel_types:
            for ty2 in range(1,ntypes+1):
                type1, type2 = min(ty, ty2), max(ty, ty2)
                idx = self.parm.parm_data['NONBONDED_PARM_INDEX'][
                            ntypes*(type1-1)+type2-1]
                acoef = self.parm.parm_data['LENNARD_JONES_ACOEF'][idx-1]
                bcoef = self.parm.parm_data['LENNARD_JONES_BCOEF'][idx-1]
                if bcoef == 0 or acoef == 0:
                    rij = eij = 0.0
                else:
                    rij = (2 * acoef / bcoef) ** (1 / 6)
                    eij = (bcoef * bcoef / (4 * acoef))
                ret_str += ('%%%ds %%%ds %%15.6f %%15.6f %%10.6f %%10.6f\n' %
                            (maxlen, maxlen) %
                            (typenames[type1-1], typenames[type2-1],
                             acoef, bcoef, rij, eij)
                )

        return ret_str

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class timerge(Action):
    """
    Merges molecules removing redundant bonding terms.  Input amber masks
    corresponding to molecules 1/2 <mol1mask>/<mol2mask>, and the soft core
    atoms in each molecule as <scmask1>/<scmask2>. The input topology can be
    created using leap, with the two molecules to be merged adjacent to each
    other in residue number. This improves the efficiency for pmemd TI when only
    part of a molecule is being perturbed.
       
    <scmask1/2N> are for softcore molecules that are not going to be merged.
    These options will just add these atoms to the timask output, correcting for
    any changes in atom number.

    This can also be used for non-softcore simulations, where
    <scmask1>/<scmask2> represent the perturbed atoms. The output will give the
    scmask1/scmask2 flags, which can just be ignored.

    <tol> is the tolerence to use when matching coordinates (default 0.0001).
    This is used when the atoms in molecules 1/2 are not in the same order and
    for checking the input coordinates.
    """
    supported_classes = ('AmberParm',)

    def init(self, arg_list):
        self.tol = arg_list.get_key_float('tol', 0.0001)
        self.molmask1 = AmberMask(self.parm, arg_list.get_next_mask())
        self.molmask2 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask1 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask2 = AmberMask(self.parm, arg_list.get_next_mask())
        self.sc_mask1 = ''
        self.sc_mask2 = ''  

        molmask1N = arg_list.get_next_mask(optional=True)
        if molmask1N is not None:
            self.molmask1N = AmberMask(self.parm, molmask1N)
        else:
            self.molmask1N = None

        molmask2N = arg_list.get_next_mask(optional=True)
        if molmask2N is not None:
            self.molmask2N = AmberMask(self.parm, molmask2N)
        else:
            self.molmask2N = None

        self.err = None   

    def __str__(self):
        sel1 = self.mask1.Selection()
        sel2 = self.mask2.Selection()
        molsel1 = self.molmask1.Selection()
        molsel2 = self.molmask2.Selection()

        natom = self.parm.ptr('natom')

        if self.molmask1N is not None:
            molsel1N = self.molmask1N.Selection()
        else:
            molsel1N = [0 for i in range(natom)]

        if self.molmask2N is not None:
            molsel2N = self.molmask2N.Selection()
        else:
            molsel2N = [0 for i in range(natom)]

        for i in range(natom):
            if sel1[i] and not molsel1[i]:
                self.err = 'ERROR: scmask1 must be a subset of mol1mask'
                return None
            if sel2[i] and not molsel2[i]:
                self.err = 'ERROR: scmask2 must be a subset of mol2mask'
                return None
            if molsel1[i] and molsel2[i]:
                self.err = 'ERROR: mol1mask can not overlap with mol2mask.'
                return None

        if not hasattr(self.parm, 'coords'):
            self.err = 'ERROR: Load coordinates before merging topology.'
            return None
      
        # we now have enough info to remap the atom indicies if an atom in
        # molsel2 has no overlap (dihedrals, angles, bonds) with sel2 then we
        # can just delete it (it is redundant).

        keep_mask = [0 for i in range(natom)]

        for i in range(natom):
            if molsel2[i]:
                for j in range(natom):
                    if sel2[j]:
                        atm1 = self.parm.atom_list[i]
                        atm2 = self.parm.atom_list[j]

                        if (atm1 in atm2.bond_partners or
                            atm1 in atm2.angle_partners or
                            atm1 in atm2.dihedral_partners):
                            keep_mask[i] = 1

        nremove = sum(molsel2) - sum(sel2)
      
        #We use an ambmask here to minimize changes to the code and
        #not introduce extra maintenance issues
        remove_mask = []
        #remove_map[old_atm_idx] = new_atm_idx
        remove_map = [0 for i in range(natom)] 

        new_atm_idx = 0
        for i in range(natom):
            if molsel2[i] == 1 and sel2[i] == 0:
                remove_mask.append('%d' % (i+1))
            else:
                remove_map[i] = new_atm_idx
                new_atm_idx += 1

        remove_str = '@' + ','.join(remove_mask)

        #However, if there is overlap, we need to re-index the atoms involved
        #Create a map from molsel2 to molsel1 excluding sel2/sel1 respectively.

        mol1common = []
        mol2common = []

        for i in range(natom):
            if molsel1[i] == 1 and sel1[i] == 0:                  
                mol1common.append(i)

        for i in range(natom):
            if molsel2[i] == 1 and sel2[i] == 0:                  
                mol2common.append(i)

        if len(mol1common) != len(mol2common):
            self.err = ('ERROR: The number of nonsoftcore atoms in mol1mask '
                        'and mol2mask must be the same.')
            return None

        mol2common_sort = []
        #reorder mol2common so that it matches mol1common
        for i in range(len(mol1common)):
            atm_i = mol1common[i]
            for j in range(len(mol2common)):
                atm_j = mol2common[j]
                diff_count = 0
                for k in range(3):
                    diff = (self.parm.coords[3*atm_i + k] - 
                            self.parm.coords[3*atm_j + k])
                    if abs(diff) < self.tol:
                        diff_count += 1

                if diff_count == 3:
                    mol2common_sort.append(atm_j)

        mol2common = mol2common_sort

        #check again if we didn't match all coords
        if len(mol1common) != len(mol2common):
            self.err = ('ERROR: The number of nonsoftcore atoms in mol1mask '
                        'and mol2mask must be the same.')
            return None

        for i in range(len(mol1common)):
            atm_i = mol1common[i]
            atm_j = mol2common[i]               
            for k in range(3):
                diff = (self.parm.coords[3*atm_i + k] - 
                        self.parm.coords[3*atm_j + k])
                if abs(diff) > self.tol:
                    self.err = ('ERROR: Common (nonsoftcore) atoms must have'
                                'the same coordinates.')
                    return None
      
        for j in range(natom):
            if keep_mask[j] == 1 and sel2[j] == 0:
                atm = self.parm.atom_list[j]
                idx = mol1common[mol2common.index(j)]
                atm_new = self.parm.atom_list[idx]

                for k in range(natom):
                    if sel2[k]:
                        atm2 = self.parm.atom_list[k]
                        #update partners -- the exclusion list will be updated 
                        #when the file is written out
                        if atm in atm2.bond_partners:
                            atm.bond_partners.remove(atm2)
                            atm2.bond_partners.remove(atm)
                            atm2.bond_to(atm_new)
                            atm_new.bond_to(atm2)

                        if atm in atm2.angle_partners:
                            atm.angle_partners.remove(atm2)
                            atm2.angle_partners.remove(atm)
                            atm2.angle_to(atm_new)
                            atm_new.angle_to(atm2)

                        if atm in atm2.dihedral_partners:
                            atm.dihedral_partners.remove(atm2)
                            atm2.dihedral_partners.remove(atm)
                            atm2.dihedral_to(atm_new)
                            atm_new.dihedral_to(atm2)

                        #Now go through each array re-indexing the atoms
                        #Check to make sure that this is a bond/angle/dihed 
                        #involving the common atom j and the softcore atom k
                      
                        for bond in (self.parm.bonds_inc_h, 
                                     self.parm.bonds_without_h):
                            for i in range(len(bond)):
                                holder = bond[i]
                                if (holder.atom1.starting_index == j and 
                                    holder.atom2.starting_index == k):
                                    holder.atom1 = atm_new
                                elif (holder.atom2.starting_index == j and 
                                      holder.atom1.starting_index == k):
                                    holder.atom2 = atm_new

                        for angle in (self.parm.angles_inc_h, 
                                      self.parm.angles_without_h):
                            for i in range(len(angle)):
                                holder = angle[i]
                                if holder.atom1.starting_index == j:
                                    if (holder.atom2.starting_index == k or 
                                        holder.atom3.starting_index == k):
                                        holder.atom1 = atm_new
                                elif holder.atom2.starting_index == j:
                                    if (holder.atom1.starting_index == k or 
                                        holder.atom3.starting_index == k):
                                        holder.atom2 = atm_new
                                elif holder.atom3.starting_index == j:
                                    if (holder.atom1.starting_index == k or 
                                        holder.atom2.starting_index == k):
                                        holder.atom3 = atm_new
                  
                        for dihed in (self.parm.dihedrals_inc_h, 
                                      self.parm.dihedrals_without_h):
                            for i in range(len(dihed)):
                                holder = dihed[i]
                                if holder.atom1.starting_index == j:
                                    if (holder.atom2.starting_index == k or 
                                        holder.atom3.starting_index == k or 
                                        holder.atom4.starting_index == k):
                                        holder.atom1 = atm_new
                                elif holder.atom2.starting_index == j:
                                    if (holder.atom1.starting_index == k or 
                                        holder.atom3.starting_index == k or 
                                        holder.atom4.starting_index == k):
                                        holder.atom2 = atm_new
                                elif holder.atom3.starting_index == j:
                                    if (holder.atom1.starting_index == k or 
                                        holder.atom2.starting_index == k or 
                                        holder.atom4.starting_index == k):
                                        holder.atom3 = atm_new
                                elif holder.atom4.starting_index == j:
                                    if (holder.atom1.starting_index == k or 
                                        holder.atom2.starting_index == k or 
                                        holder.atom3.starting_index == k):
                                        holder.atom4 = atm_new               

        self.parm.atom_list.changed = True

        if nremove > 0:
            self.parm.delete_mask(remove_str)

        new_sc_atm1 = []
        new_sc_atm2 = []
        new_sc_atm1_int = []
        new_sc_atm2_int = []
        for i in range(natom):
            if sel1[i] or molsel1N[i]:
                new_sc_atm1_int.append(remove_map[i])
                new_sc_atm1.append('%d' % (remove_map[i]+1))
            elif sel2[i] or molsel2N[i]:
                new_sc_atm2_int.append(remove_map[i])
                new_sc_atm2.append('%d' % (remove_map[i]+1))

        # Do not allow definition where dihedrals cross through the softcore 
        # region. This is generally breaking a ring (and possibly other cases), 
        # and can cause problems with the 1-4 nonbonded calculations.
        # This can be worked-around: 
        # Define your softcore region so that it includes the ring.        
        for dihed in (self.parm.dihedrals_inc_h,self.parm.dihedrals_without_h):
            for i in range(len(dihed)):
                holder = dihed[i]
                # skip impropers, these are not used to define 1-4 interactions
                # so these can cross through the softcore region
                if holder.signs[1] < 0: continue
                atmi = holder.atom1.starting_index
                atmj = holder.atom2.starting_index
                atmk = holder.atom3.starting_index
                atml = holder.atom4.starting_index
                if (atmj in new_sc_atm1_int or 
                    atmk in new_sc_atm1_int or 
                    atmj in new_sc_atm2_int or 
                    atmk in new_sc_atm2_int): #dihedral includes sc atoms 
               
                    #endpoint atoms are not softcore
                    #we are crossing through the softcore region
                    if (atmi not in new_sc_atm1_int and 
                        atmi not in new_sc_atm2_int and 
                        atml not in new_sc_atm1_int and 
                        atml not in new_sc_atm2_int):
                        self.err = (
                                'ERROR: Can not have dihedral cross through '
                                'softcore region. (DIHED : %d %d %d %d)\n '
                                'Usually this means you have defined the '
                                'softcore region in a way that breaks a '
                                'ring.\n Try redefining your softcore region '
                                'to include the ring or at least three '
                                'consecutive atoms.' %
                                (atmi+1, atmj+1, (atmk+1)*holder.signs[0], 
                                 (atml+1) * holder.signs[1])
                        )
                        return
                     
        self.sc_mask1 = '@' + ','.join(new_sc_atm1)
        self.sc_mask2 = '@' + ','.join(new_sc_atm2)

        if self.err is not None:
            return self.err
        else:
            ret_str = ("Merging molecules %s and %s into the same molecule.\n"
                       % (self.molmask1, self.molmask2))
            ret_str2 = ("Use softcore mask:\ntimask1=\'%s\',\ntimask2=\'%s\',"
                        % (self.sc_mask1, self.sc_mask2))
            ret_str3 = ("\nscmask1=\'%s\',\nscmask2=\'%s\',"
                        % (self.sc_mask1, self.sc_mask2))
            return ret_str + ret_str2 + ret_str3

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class source(Action):
    """
    Sources a file with a list of parmed commands
    """
    needs_parm = False
    def init(self, arg_list):
        self.filename = arg_list.get_next_string()

    def __str__(self):
        return 'Sourcing %s' % self.filename

    def execute(self):
        """
        This is a no-op, since a separate command interpreter for this file is
        launched inside parmed_cmd.py
        """
        pass

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class parm(Action):
    """
    Either adds a new parm to the list of available parms to edit in ParmEd, or
    it sets a new 'active' parm to edit by default with new commands
    """
    needs_parm = False # We don't need a parm for this action
    def init(self, arg_list):
        from glob import glob
        self.new_active_parm = arg_list.get_key_string('select', None)
        self.copied_parm = arg_list.get_key_string('copy', None)
        self.new_parm = None
        # Add as many parms as we want with one command, and support globbing
        new_parm = len(arg_list.unmarked()) or None
        # Make sure we specified only one operating mode
        if (self.new_active_parm is None and self.copied_parm is None and
                new_parm is None):
            raise ParmError('Improper usage of `parm\' command')
        opts = (self.new_active_parm, self.copied_parm, new_parm)
        nnone = 0
        for opt in opts:
            if opt is None:
                nnone += 1
        if nnone != 2:
            raise ParmError('Improper usage of `parm\' -- choose one behavior')
        # Now get our indexes for copied and new parms if they are integers,
        # otherwise we index based on filename
        if self.new_active_parm is not None:
            try:
                self.new_active_parm = int(self.new_active_parm)
            except ValueError:
                pass
        elif self.copied_parm is not None:
            try:
                self.copied_parm = int(self.copied_parm)
            except ValueError:
                pass
        else:
            self.new_parm = []
            new_parm = arg_list.get_next_string(optional=True)
            while new_parm is not None:
                listparms = glob(new_parm)
                if not listparms:
                    warnings.warn('No files matching %s' % new_parm,
                                  NonexistentParmWarning)
                    continue
                self.new_parm.extend(glob(new_parm))
                new_parm = arg_list.get_next_string(optional=True)
            if not self.new_parm:
                raise NonexistentParm('No matching parm files')

    def __str__(self):
        if self.new_active_parm is not None:
            try:
                idx = self.parm_list.index(self.new_active_parm)
            except IndexError:
                return ('%s not in parm list. Doing nothing' %
                        self.new_active_parm)
            return 'Setting new active parm [%s]' % self.parm_list[idx]
        elif self.new_parm is not None:
            return ('Adding prmtop %s to parm list. %s is the active parm.' %
                    (', '.join(self.new_parm), self.new_parm[-1]))
        elif self.copied_parm is not None:
            try:
                idx = self.parm_list.index(self.copied_parm)
            except IndexError:
                return ('%s not in parm list. Doing nothing' %
                        self.new_active_parm)
            return ("Copying prmtop %s to parm list. %s's copy is the active "
                    "parm." % (self.parm_list[idx], self.parm_list[idx]))
        return 'Internal error!' # should never reach here
   
    def execute(self):
        """ Either set the new active parm or add the new parm """
        from copy import copy
        if self.new_parm is not None:
            # Add the new parm
            for new_parm in self.new_parm:
                try:
                    self.parm_list.add_parm(new_parm)
                except IOError:
                    warnings.warn('Could not open %s for reading' % new_parm,
                                  SeriousParmWarning)
        elif self.new_active_parm is not None:
            try:
                self.parm_list.set_new_active(self.new_active_parm)
            except IndexError:
                warnings.warn('%s is not in the parm list!' %
                              self.new_active_parm, SeriousParmWarning)
        elif self.copied_parm is not None:
            try:
                self.parm_list.add_parm(copy(self.parm_list[self.copied_parm]))
            except IndexError:
                warnings.warn('%s is not in the parm list!' % self.copied_parm,
                              SeriousParmWarning)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class ls(Action):
    """
    Lists directory contents. Like UNIX 'ls'
    """
    needs_parm = False
    def init(self, arg_list):
        from glob import glob
        self.args = []
        # Process the argument list to mimic the real ls as much as possible
        while True:
            try:
                arg = arg_list.get_next_string()
                if not arg.startswith('-'):
                    # Glob this argument
                    globarg = glob(arg)
                    if len(globarg) > 0:
                        self.args.extend(globarg)
                    else:
                        self.args.append(arg)
                else:
                    self.args.append(arg)
            except NoArgument:
                break

    def __str__(self):
        from subprocess import Popen, PIPE
        process = Popen(['/bin/ls', '-C'] + self.args, stdout=PIPE, stderr=PIPE)
        out, err = process.communicate('')
        process.wait()
        return out + err

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class cd(Action):
    """
    Changes to a new directory like UNIX 'cd'
    """
    needs_parm = False
    def init(self, arg_list):
        from glob import glob
        from os.path import expanduser, expandvars
        mydir = expanduser(expandvars(arg_list.get_next_string()))
        self.directory = glob(mydir)

    def __str__(self):
        if len(self.directory) != 1:
            return 'Change directory failed'
        if not os.path.isdir(self.directory[0]):
            return '%s does not exist. cd failed.' % self.directory[0]
        return 'New working directory: %s' % self.directory[0]

    def execute(self):
        if len(self.directory) < 1:
            warnings.warn('No recognized directories given to cd',
                          SeriousParmWarning)
            return
        elif len(self.directory) > 1:
            warnings.warn('More than one file/directory given to cd',
                          SeriousParmWarning)
            return
        if not os.path.isdir(self.directory[0]):
            warnings.warn('%s is not a directory' % self.directory[0],
                          SeriousParmWarning)
            return
        # If we've gotten this far, go ahead and change to the directory
        os.chdir(self.directory[0])

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class listparms(Action):
    """
    Lists all of the loaded topology files
    """
    needs_parm = False

    def init(self, arg_list):
        pass

    def __str__(self):
        if self.parm_list.empty():
            return "No topology files are loaded"

        retstr = 'Loaded topology files:'
        for i, parm in enumerate(self.parm_list):
            retstr += '\n[%d]\t%s' % (i, parm)
            if parm is self.parm_list.parm:
                retstr += ' (active)'

        return retstr

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class interpolate(Action):
    """
    Interpolates between two topology files (VDW and electrostatic terms only).
    If [eleconly] is present, only the charges will be interpolated. <nparm> is
    the number of 'interpolated' topology files you want (in addition to the two
    end-points). The second prmtop must be specified if there are more than 2 in
    the list.
    """
    supported_classes = ('AmberParm', 'ChamberParm')
    def init(self, arg_list):
        # Make sure we have at least 2 prmtops
        if len(self.parm_list) < 2:
            raise NonexistentParm('Must have 2 topology files to interpolate!')
        parm2 = arg_list.get_key_string('parm2', None)
        if parm2 is None and len(self.parm_list) == 2:
            self.parm2 = self.parm_list[1]
            if self.parm_list[0] != self.parm:
                self.parm2 = self.parm_list[0]
        elif parm2 is None:
            raise AmbiguousParmError('You must identify parm2 if more than 2 '
                                     'parm instances exist!')
        else:
            try:
                parm2 = int(parm2)
            except ValueError:
                pass
            self.parm2 = self.parm_list[parm2]

        self.startnum = arg_list.get_key_int('startnum', 1)
        self.prefix = arg_list.get_key_string('prefix', str(self.parm))
        self.eleconly = arg_list.has_key('eleconly')
        self.nparm = arg_list.get_next_int()
        if self.nparm <= 0:
            raise ArgumentError('Must have >= 1 prmtop')
        self.diff_vdw = False
        self._check_parms()

    def __str__(self):
        extra = ''
        if self.eleconly and self.diff_vdw:
            extra = ' [only interpolating charges]'
        return 'Creating %d interpolated prmtops between %s and %s' % (
                                    self.nparm, self.parm, self.parm2) + extra

    def _check_parms(self):
        """ Makes sure that the atoms in both parms are all the same """
        parm1, parm2 = self.parm, self.parm2
        if parm1.ptr('natom') != parm2.ptr('natom'):
            raise IncompatibleParmsError('%s and %s have different #s of '
                                         'atoms!' % (parm1, parm2))
        ndiff = 0
        for atom1, atom2 in zip(parm1.atom_list, parm2.atom_list):
            if atom1.atname != atom2.atname:
                ndiff += 1
        if ndiff > 0:
            warnings.warn('%d atoms have different names b/w %s and %s' %
                          (ndiff, parm1, parm2), SeriousParmWarning)
        for atm1, atm2 in zip(parm1.atom_list, parm2.atom_list):
            if ((parm1.LJ_radius[atm1.nb_idx-1] !=
                 parm2.LJ_radius[atm2.nb_idx-1]) or
                (parm1.LJ_depth[atm1.nb_idx-1] !=
                 parm2.LJ_radius[atm2.nb_idx-1])):
                self.diff_vdw = True

    def execute(self):
        """ Interpolates the prmtops """
        from ParmedTools.arraytools import NumberArray
        if self.diff_vdw and not self.eleconly:
            raise NotImplemented('No support for scaling vdW parameters yet!')

        parm1, parm2 = self.parm, self.parm2
        # Original charges for parm 1
        orig_chg1 = parm1.parm_data['CHARGE']
        chg1 = NumberArray(parm1.parm_data['CHARGE'])
        chg2 = NumberArray(parm2.parm_data['CHARGE'])
        diff = chg2 - chg1
        diff *= 1 / (self.nparm + 1)
        for i in range(self.nparm):
            new_chg = chg1 + diff * (i + 1)
            parm1.parm_data['CHARGE'] = [c for c in new_chg]
            newname = '%s.%d' % (self.prefix, i+self.startnum)
            print 'Printing %s' % newname
            parm1.writeParm(newname)
        # Restore the original charges
        parm1.parm_data['CHARGE'] = orig_chg1

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class summary(Action):
    """
    Prints out a summary of prmtop contents
    """

    nucleic = ['A', 'G', 'C', 'U', 'DA', 'DG', 'DC', 'DT', 'AP', 'CP', 'DAP',
               'DCP', 'AE', 'CE', 'DCE', 'DAE', 'GE', 'DGE', 'DTE', 'UE']

    amino = ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'AS4', 'CYM', 'CYS', 'CYX',
             'GLH', 'GLN', 'GLU', 'GLY', 'GL4', 'HID', 'HIE', 'HIP', 'HYP',
             'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR',
             'TRP', 'TYR', 'VAL']

    anions = ['Cl-', 'Br-', 'F-', 'I-']

    cations = ['Na+', 'Li+', 'Mg+', 'Rb+', 'MG', 'Cs+']

    solvent = ['WAT', 'HOH']

    def init(self, arg_list):
        pass

    def __str__(self):
        """ Collect statistics """
        nnuc = namin = ncion = naion = nwat = nunk = 0
        for i, res in enumerate(self.parm.parm_data['RESIDUE_LABEL']):
            if res in summary.nucleic:
                nnuc += 1
            elif res in summary.amino:
                namin += 1
            elif res in summary.solvent:
                nwat += 1
            elif res in summary.anions:
                naion += 1
            elif res in summary.cations:
                ncion += 1
            else:
                nunk += 1
      
        tmass = sum(self.parm.parm_data['MASS'])
        tchg = sum(self.parm.parm_data['CHARGE'])

        retval = ('Amino Acid Residues:   %d\n'
                  'Nucleic Acid Residues: %d\n'
                  'Number of cations:     %d\n'
                  'Number of anions:      %d\n'
                  'Num. of solvent mols:  %d\n' 
                  'Num. of unknown atoms: %d\n'
                  'Total charge (e-):     %.4f\n'
                  'Total mass (amu):      %.4f\n'
                  'Number of atoms:       %d\n'
                  'Number of residues:    %d\n' %
                  (namin, nnuc, ncion, naion, nwat, nunk, tchg, tmass,
                   self.parm.ptr('natom'), self.parm.ptr('nres'))
        )

        if self.parm.ptr('ifbox') == 1:
            if hasattr(self.parm, 'box'):
                a, b, c = self.parm.box[:3]
            else:
                a, b, c = self.parm.parm_data['BOX_DIMENSIONS'][1:]
            v = a * b * c
            # Get the total volume (and density) of orthorhombic box
            retval += ('System volume (ang^3): %f\n' 
                       'System density (g/mL): %f\n' %
                       (v, sum(self.parm.parm_data['MASS']) / (v * 0.602204))
            )
        elif self.parm.ptr('ifbox') == 2:
            # General triclinic cell
            if hasattr(self.parm, 'box'):
                a, b, c, alpha, beta, gamma = self.parm.box[:]
            else:
                a, b, c = self.parm.parm_data['BOX_DIMENSIONS'][1:]
                alpha = beta = gamma = self.parm.parm_data['BOX_DIMENSIONS'][0]
            # Convert to radians
            cosa = math.cos(alpha * math.pi / 180)
            cosb = math.cos(beta * math.pi / 180)
            cosg = math.cos(gamma * math.pi / 180)
            v = a * b * c * math.sqrt(1 - cosa*cosa - cosb*cosb - cosg*cosg +
                                      2 * cosa*cosb*cosg)
            retval += ('System volume (ang^3): %f\n' 
                       'System density (g/mL): %f\n' %
                       (v, sum(self.parm.parm_data['MASS']) / (v * 0.602204))
            )
        return retval

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class scale(Action):
    """
    Multiplies all values in a particular section of the prmtop by a scalar
    factor
    """

    def init(self, arg_list):
        self.flag = arg_list.get_next_string().upper()
        self.factor = arg_list.get_next_float()
        if not self.flag in self.parm.flag_list:
            raise ArgumentError('%s is not a valid prmtop flag!' % self.flag)

    def __str__(self):
        return 'Scaling data in %s by %f' % (self.flag, self.factor)

    def execute(self):
        try:
            for i in range(len(self.parm.parm_data[self.flag])):
                self.parm.parm_data[self.flag][i] *= self.factor
            self.parm.flush_data_changes()
        except TypeError:
            raise ArgumentError('Cannot scale data in %%FLAG %s' % self.flag)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class lmod(Action):
    """
    Adjusts Lennard Jones parameters to work with the LMOD code in Amber
    (changes LJ A coefficients that are 0 to 1000).
    """
    supported_classes = ('AmberParm', 'ChamberParm')

    def init(self, arg_list):
        pass
   
    def __str__(self):
        return ('Making adjustments for LMOD (LJ A-coef. for H atoms bonded '
                'to O)')

    def execute(self):
        for i, val in enumerate(self.parm.parm_data['LENNARD_JONES_ACOEF']):
            self.parm.parm_data['LENNARD_JONES_ACOEF'][i] = val or 1000.0

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class addpdb(Action):
    """
    Adds PDB information to new flags in the topology file to enable analyses
    based on the original residue information in the PDB file, <filename>. It
    adds the flags:
        RESIDUE_CHAINID: The chain ID of each residue (* if LEaP added it)
        RESIDUE_ICODE: Insertion code, if it exists
        RESIDUE_NUMBER: Original residue number in the PDB
        ATOM_ELEMENT: Atomic element (redundant now, not printed by default)

    The 'strict' keyword turns residue mismatches (NOT solvent) into errors
    The 'elem' keyword will force printing of the element names.
    The 'allicodes' keyword forces insertion codes to be printed, even if every
    one will be blank (so that parsers can count on that section existing).

    Residues _not_ in the PDB will be assigned a CHAINID of '*' and
    RESIDUE_NUMBER of 0.

    Historical info:
        This action is based on, and reproduces the key results of, the
        historical 'add_pdb' program. It is a bit more flexible, though.
    """

    def init(self, arg_list):
        self.pdbfile = arg_list.get_next_string()
        self.elements = arg_list.has_key('elem')
        self.printicodes = arg_list.has_key('allicodes')
        if arg_list.has_key('strict'):
            warnings.filterwarnings('error', category=AddPDBWarning)
        else:
            warnings.filterwarnings('always', category=AddPDBWarning)
        self.pdbpresent = ('RESIDUE_NUMBER' in self.parm.flag_list or
                           'RESIDUE_CHAINID' in self.parm.flag_list or
                           'RESIDUE_ICODE' in self.parm.flag_list or
                           'ATOM_ELEMENT' in self.parm.flag_list
        )

    def __str__(self):
        if self.pdbpresent:
            return 'PDB information already present in %s. Doing nothing'
        retstr = 'Adding PDB information from %s' % self.pdbfile
        if self.elements: retstr += '\n\t[printing elements from prmtop]'
        return retstr

    def execute(self):
        from chemistry.system import ChemicalSystem
        if self.pdbpresent: return
        pdb = ChemicalSystem.load_from_pdb(self.pdbfile)
        resnums = [0 for i in range(self.parm.ptr('nres'))]
        chainids = ['*' for i in range(self.parm.ptr('nres'))]
        icodes = ['' for i in range(self.parm.ptr('nres'))]
        for i, res in enumerate(pdb):
            try:
                reslab = self.parm.parm_data['RESIDUE_LABEL'][i].strip()
                resname = res.name.strip()
                if resname != reslab:
                    if (not resname in ('WAT', 'HOH') or 
                        not reslab in ('WAT', 'HOH')):
                        # Allow for 3's and 5's in terminal nucleic acid
                        # residues
                        if resname in ('A', 'C', 'G', 'U', 'DA',
                                       'DG', 'DC', 'DT'):
                            if reslab[-1] in '35':
                                reslab = reslab[:-1]
                        if reslab != resname:
                            warnings.warn('Residue name mismatch [#%d] %s vs. '
                                          '%s' % (i+1, resname, reslab),
                                          AddPDBWarning)
                resnums[i] = res.number
                chainids[i] = res.chain
                icodes[i] = res.insertion_code.strip()
            except IndexError:
                raise AddPDBError('PDB %s has more residues than prmtop %s' %
                                (self.pdbfile, self.parm))
        ncmts = ['Residue number (resSeq) read from PDB file; DIMENSION(NRES)']
        if self.printicodes or any(icodes):
            haveicodes = True
            ncmts += ['Residue insertion codes (iCode) present in %FLAG '
                      'RESIDUE_ICODE']
        else:
            haveicodes = False
            ncmts += ['Residue insertion code (iCode) not present in PDB file',
                      'If present: %FLAG RESIDUE_ICODE, %FORMAT(20a4)']
        self.parm.addFlag('RESIDUE_NUMBER', '20I4', data=resnums,
                          comments=ncmts)
        self.parm.addFlag('RESIDUE_CHAINID', '20a4', data=chainids,
                          comments=['Residue chain ID (chainId) read from PDB '
                                    'file; DIMENSION(NRES)']
        )
        if haveicodes:
            self.parm.addFlag('RESIDUE_ICODE', '20a4', data=icodes,
                comments=['Residue insertion code (iCode) read from PDB file; '
                        'DIMENSION(NRES)']
            )
        if self.elements:
            self.parm.addFlag('ATOM_ELEMENT', '20a4',
                data=['%2s' % (_Element[atm.atomic_number].upper()) 
                            for atm in self.parm.atom_list
                ], comments=['Atom element name read from topology file']
            )

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class deletepdb(Action):
    """
    This action deletes the flags added by addPDB if they are present.
    """

    def init(self, arg_list):
        self.pdbpresent = ('RESIDUE_NUMBER' in self.parm.flag_list or
                           'RESIDUE_CHAINID' in self.parm.flag_list or
                           'RESIDUE_ICODE' in self.parm.flag_list or
                           'ATOM_ELEMENT' in self.parm.flag_list
        )

    def __str__(self):
        if self.pdbpresent:
            return 'Deleting PDB info from %s' % self.parm
        return 'No PDB information present in %s. Doing nothing' % self.parm

    def execute(self):
        if not self.pdbpresent: return
        self.parm.deleteFlag('RESIDUE_NUMBER')
        self.parm.deleteFlag('RESIDUE_CHAINID')
        self.parm.deleteFlag('ATOM_ELEMENT')
        self.parm.deleteFlag('RESIDUE_ICODE')

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class add12_6_4(Action):
    """
    Adds the LENNARD_JONES_CCOEF term for the new divalent metal ion 12-6-4
    Lennard-Jones potential term. If provided, the mask will allow you to
    specify which ions are divalent. The C4 parameter between the metal ion and
    water can either be taken from Ref. [1] for the requested [watermodel]
    (TIP3P, TIP4PEW, or SPCE) or provided in the file specified by the c4file
    keyword.  The polarizabilities must be present in the the polfile file. The
    chemical symbol of the element will be used to determine the atom type.
    Parameters are expected in a file with 2 columns:
        <atom type>   <parameter>

    All defaults come from Ref. [1]

    [1] Pengfei Li and Kenneth M. Merz, J. Chem. Theory Comput., 2014, 10,
        289-297
    """
    supported_classes = ('AmberParm', 'ChamberParm')
   
    _supported_wms = ('TIP3P', 'TIP4PEW', 'SPCE')

    def init(self, arg_list):
        import os
        self.mask = AmberMask(self.parm,
                        arg_list.get_next_mask(optional=True, default=':ZN'))
        self.c4file = arg_list.get_key_string('c4file', None)
        self.watermodel = arg_list.get_key_string('watermodel', None)
        self.polfile = arg_list.get_key_string('polfile',
                            os.path.join(os.getenv('AMBERHOME'), 'dat', 'leap',
                            'parm', 'lj_1264_pol.dat'))
        self.tunfactor = arg_list.get_key_float('tunfactor', None)

        if self.c4file is None:
            if self.watermodel is None:
                self.watermodel = 'TIP3P'
        else:
            self.watermodel = self.watermodel.upper()
            if not self.watermodel in self._supported_wms:
                raise LJ12_6_4Error('Defaults exist for water models ' +
                                    ', '.join(self._supported_wms))
            else:
                if self.watermodel is not None:
                    raise LJ12_6_4Error('c4file and watermodel are mutually '
                                        'exclusive')

    def __str__(self):
        retstr = ('Adding Lennard-Jones C-coefficient for 12-6-4 potential. '
                  'Polarizabilities read from [%s]. ' % self.polfile)
        if self.c4file is None:
            retstr += ('Using default C4 parameters for water model [%s].' %
                       self.watermodel)
        else:
            retstr += 'Reading C4 parameters from [%s].' % self.c4file

        return retstr

    def execute(self):
        from ParmedTools.add1264 import params1264 as params
        if 'LENNARD_JONES_CCOEF' in self.parm.flag_list:
            self.parm.deleteFlag('LENNARD_JONES_CCOEF')
        self.parm.addFlag('LENNARD_JONES_CCOEF', '5E16.8',
                num_items=len(self.parm.parm_data['LENNARD_JONES_ACOEF']),
                comments=['For 12-6-4 potential used for divalent metal ions'])
        for i, param in enumerate(params(self.parm, self.mask, self.c4file,
                                         self.watermodel, self.polfile,
                                         self.tunfactor)):
            self.parm.parm_data['LENNARD_JONES_CCOEF'][i] = param

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class hmassrepartition(Action):
    """
    This action implements hydrogen mass repartitioning in the system by
    changing the mass of each hydrogen to the desired value (the default new
    hydrogen mass is 3.024 daltons) and adjusting the mass of the atom to which
    it is bonded by the amount required to leave the total mass unchanged. By
    default, water hydrogen masses are unchanged (the SETTLE algorithm for
    implementing constraints on water molecules is analytical). Water masses can
    be repartitioned as well with the 'dowater' keyword.
    """

    def init(self, arg_list):
        self.changewater = arg_list.has_key('dowater')
        self.new_h_mass = arg_list.get_next_float(optional=True, default=3.024)

    def __str__(self):
        retstr = ('Repartitioning hydrogen masses to %s daltons. ' %
                  self.new_h_mass)
        if self.changewater:
            return retstr + 'Also changing water hydrogen masses.'
        return retstr + 'Not changing water hydrogen masses.'

    def execute(self):
        # Back up the masses in case something goes wrong
        original_masses = self.parm.parm_data['MASS'][:]
        for i, atom in enumerate(self.parm.atom_list):
            if atom.atomic_number != 1: continue
            if not self.changewater and atom.residue.resname in ('WAT', 'HOH'):
                continue
            heteroatom = None
            heteroidx = 0
            bondeds = list(atom.bond_partners)
            while heteroidx < len(bondeds):
                if bondeds[heteroidx].atomic_number != 1:
                    heteroatom = bondeds[heteroidx]
                    break
                heteroidx += 1
            if heteroatom is None:
                # Only bonded to other hydrogens. Weird, but do not repartition
                warnings.warn('H atom detected not bound to heteroatom. '
                              'Ignoring.', ParmWarning)
                continue
            transfermass = self.new_h_mass - self.parm.parm_data['MASS'][i]
            oi = heteroatom.starting_index
            self.parm.parm_data['MASS'][i] = self.new_h_mass
            self.parm.parm_data['MASS'][oi] -= transfermass

        # Now make sure that all masses are positive, or revert masses and
        # raise an exception
        for i, mass in enumerate(self.parm.parm_data['MASS']):
            if mass <= 0 and self.parm.atom_list[i].atomic_number > 0:
                self.parm.parm_data['MASS'] = original_masses
                raise HMassRepartitionError('Too much mass removed from atom '
                                            '%d. Hydrogen masses must be '
                                            'smaller.' % i)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class openmm(Action):
    """
    This class will read a sander/pmemd input file and run an equivalent
    simulation using the OpenMM Python API. It uses a topology file that has
    been defined (and/or modified) here. The command-line flags are identical to
    those used for sander and pmemd (described in the Amber manual). You can
    additionally specify a computational platform available through OpenMM such
    as CUDA, OpenCL, Reference, and CPU.

    The default prmtop will be the 'active' parm. The prmtop can be changed
    using either the 'parm' keyword or the '-p' flag, but it must resolve to one
    of the stored topology files. Any parm given with the '-p' flag has
    precedence. If present, the 'dcd' keyword triggers writting DCD-formatted
    trajectory files instead of NetCDF when ioutfm=1. The DCD trajectory writing
    offers binary trajectory support without requiring a NetCDF-Python library.
   
    The 'progress' keyword triggers printing of ParmEd's progress in setting up
    and starting the calculation.

    The 'script' keyword provides a file to write an OpenMM script that performs
    the same calculation requested by this ParmEd command. The 'norun' command
    will prevent anything from actually being run and can only be used when a
    script file is provided.
    """
    supported_classes = ('AmberParm', 'ChamberParm')

    def init(self, arg_list):
        parm = arg_list.get_key_string('-p', default=None)
        if parm is not None:
            self.parm = self.parm_list[parm]
        # This consumes the remaining arguments
        self.arg_list = ArgumentList(arg_list)

    def __str__(self):
        return ("Running Amber-style simulation with OpenMM using the command "
                "string:\n\t%s\nThis may take awhile..." % self.arg_list)

    def execute(self):
        """ Runs the OpenMM simulation """
        from ParmedTools.simulations.openmm import simulate, HAS_OPENMM
        if not HAS_OPENMM:
            raise SimulationError('OpenMM could not be imported. Skipping.')
        # First try to load a restart file if it was supplied
        inptraj = self.arg_list.has_key('-y', mark=False)
        has_inpcrd = self.arg_list.has_key('-c', mark=False)
        if not hasattr(self.parm, 'rst7') and not inptraj and not has_inpcrd:
            raise SimulationError('No input coordinates provided.')
        # Eliminate some incompatibilities that are easy to catch now
        if self.parm.ptr('ifbox') > 1:
            raise SimulationError('OpenMM only supports orthorhombic boxes '
                                  'currently.')

        simulate(self.parm, self.arg_list)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class energy(Action):
    """
    This action will use OpenMM to compute a single-point energy for the loaded
    structure (you must use 'loadRestart' prior to this command to load
    coordinates). The following options and keywords are supported:

    The options are:

    platform <platform> : OpenMM compute platform to use. Options are CUDA,
            OpenCL, Reference, and CPU. Consult the OpenMM manual for details
    precision <precision model> : Precision model to use. Options are single,
         double, and mixed. Reference platform is always double and CPU platform
         is always single. Mixed (default) uses single precision for
         calculations and double for accumulation
    cutoff <cut> : The size of the non-bonded cutoff, in Angstroms. Default 8 A
         for periodic systems or infinite for nonperiodic systems

    For systems with no periodic box information:

    igb <IGB> : An integer corresponding to the desired GB model. May be 1, 2,
         5, 7, or 8 as described by the sander and pmemd manual. Default 5
    saltcon <conc> : The salt concentration, in mol/L, modeled by a Debye
         screening parameter. Default 0.0

    For periodic systems:

    Ewald : Use an Ewald sum to compute long-range electrostatics instead of PME
    nodisper : Do not use a long-range vdW dispersion correction

    Other options:

    decompose : Print bond, angle, dihedral, and nonbonded energies separately
    applayer : Use OpenMM's class to compute the energy
    """
    supported_classes = ('AmberParm', 'ChamberParm')

    output = sys.stdout

    def init(self, arg_list):
        self.arg_list = ArgumentList(arg_list)

    def __str__(self):
        return 'Computing a single-point energy for %s' % self.parm

    def execute(self):
        from ParmedTools.simulations.openmm import energy, HAS_OPENMM
        if not HAS_OPENMM:
            raise SimulationError('OpenMM could not be imported. Skipping.')

        energy(self.parm, self.arg_list, self.output)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class deletebond(Action):
    """
    This action deletes any bonds that occur between the atoms in two masks.

    <mask1> : Amber mask defining one of the atoms in a bond
    <mask2> : Amber mask defining the other atom in the bond
    [verbose] : Print out every bond that is deleted as well as the number of
                other valence terms that were eliminated.

    All bonds will be matched in which one atom comes from <mask1> and the other
    atom comes from <mask2>. This action will also delete any other valence term
    (like angles and dihedrals) that would get severed by the deletion of one of
    the bonds.
    """
    supported_classes = ('AmberParm', 'ChamberParm')

    def init(self, arg_list):
        self.mask1 = AmberMask(self.parm, arg_list.get_next_mask())
        self.mask2 = AmberMask(self.parm, arg_list.get_next_mask())
        self.verbose = arg_list.has_key('verbose')
        # Go through each atom in mask1 and see if it forms a bond with any atom
        # in mask2.
        self.del_h_bonds = set()
        self.del_noh_bonds = set()
        deleted_bond_list = set()
        for i in self.mask1.Selected():
            ai = self.parm.atom_list[i]
            for j in self.mask2.Selected():
                aj = self.parm.atom_list[j]
                # Skip if these two atoms are identical
                if ai is aj: continue
                if aj in ai.bonds():
                    # A bond exists here. Now find the bond in the appropriate
                    # list (is there a hydrogen or not?) and add its index to
                    # the relevant list.
                    if 1 in (ai.atomic_number, aj.atomic_number):
                        for ii, bond in enumerate(self.parm.bonds_inc_h):
                            if ai in bond and aj in bond:
                                self.del_h_bonds.add(ii)
                                deleted_bond_list.add(bond)
                                break
                    else:
                        for ii, bond in enumerate(self.parm.bonds_without_h):
                            if ai in bond and aj in bond:
                                self.del_noh_bonds.add(ii)
                                deleted_bond_list.add(bond)
                                break
        # Now go through all of our other valence terms and collect the terms
        # we need to delete.
        if not deleted_bond_list:
            # Nothing to delete
            return
        self.del_h_angles = set()
        self.del_noh_angles = set()
        self.del_h_dihedrals = set()
        self.del_noh_dihedrals = set()
        if self.parm.chamber:
            self.del_impropers = set()
            self.del_ureybrad = set()
            self.del_cmap = set()
        for bond in deleted_bond_list:
            for i, angle in enumerate(self.parm.angles_inc_h):
                if bond in angle:
                    self.del_h_angles.add(i)
            for i, angle in enumerate(self.parm.angles_without_h):
                if bond in angle:
                    self.del_noh_angles.add(i)
            for i, dihed in enumerate(self.parm.dihedrals_inc_h):
                if bond in dihed:
                    self.del_h_dihedrals.add(i)
            for i, dihed in enumerate(self.parm.dihedrals_without_h):
                if bond in dihed:
                    self.del_noh_dihedrals.add(i)
            if self.parm.chamber:
                for i, ureybrad in enumerate(self.parm.urey_bradley):
                    if bond in ureybrad:
                        self.del_ureybrad.add(i)
                for i, improper in enumerate(self.parm.improper):
                    if bond in improper:
                        self.del_impropers.add(i)
                if hasattr(self.parm, 'cmap'):
                    for i, cmap in enumerate(self.parm.cmap):
                        if bond in cmap:
                            self.del_cmap.add(i)

    def __str__(self):
        if not self.del_h_bonds and not self.del_noh_bonds:
            return 'No bonds to delete'
        if not self.verbose:
            return 'Deleting the %d bonds found between %s and %s' % (
                    len(self.del_h_bonds) + len(self.del_noh_bonds),
                    self.mask1, self.mask2)
        # Now we want full statistics
        retstr = 'Deleting %d bonds between %s and %s:\n' % (
                    len(self.del_h_bonds) + len(self.del_noh_bonds),
                    self.mask1, self.mask2)
        for i in sorted(list(self.del_h_bonds)):
            a1 = self.parm.bonds_inc_h[i].atom1
            a2 = self.parm.bonds_inc_h[i].atom2
            retstr += '\t%d [%s %d] %s --- %d [%s %d] %s\n' % (
                    a1.starting_index+1, a1.residue.resname, a1.residue.idx,
                    a1.atname, a2.starting_index+1, a2.residue.resname,
                    a2.residue.idx, a2.atname)
        for i in sorted(list(self.del_noh_bonds)):
            a1 = self.parm.bonds_without_h[i].atom1
            a2 = self.parm.bonds_without_h[i].atom2
            retstr += '\t%d [%s %d] %s --- %d [%s %d] %s\n' % (
                    a1.starting_index+1, a1.residue.resname, a1.residue.idx,
                    a1.atname, a2.starting_index+1, a2.residue.resname,
                    a2.residue.idx, a2.atname)
        retstr += 'Deleting %d angles, ' % (len(self.del_h_angles) +
                                            len(self.del_noh_angles))
        if self.parm.chamber:
            retstr += ('%d Urey-Bradleys, %d impropers,\n         %d dihedrals '
                        'and %d CMAPs' % (
                        len(self.del_ureybrad), len(self.del_impropers),
                        len(self.del_h_dihedrals) + len(self.del_noh_dihedrals),
                        len(self.del_cmap))
            )
        else:
            retstr += ' and %d dihedrals' % (len(self.del_h_dihedrals) +
                                            len(self.del_noh_dihedrals))
        return retstr

    @staticmethod
    def _dfl(selection, mylist):
        """ Delete From List """
        if not selection: return
        for i in reversed(list(selection)):
            del mylist[i]

    def execute(self):
        if not self.del_h_bonds and not self.del_noh_bonds:
            # Nothing to do...
            return
        self._dfl(self.del_h_bonds, self.parm.bonds_inc_h)
        self._dfl(self.del_noh_bonds, self.parm.bonds_without_h)
        self._dfl(self.del_h_angles, self.parm.angles_inc_h)
        self._dfl(self.del_noh_angles, self.parm.angles_without_h)
        self._dfl(self.del_h_dihedrals, self.parm.dihedrals_inc_h)
        self._dfl(self.del_noh_dihedrals, self.parm.dihedrals_without_h)
        # Now get rid of chamber sections
        if self.parm.chamber:
            self._dfl(self.del_ureybrad, self.parm.urey_bradley)
            self._dfl(self.del_impropers, self.parm.improper)
            if self.del_cmap:
                self._dfl(self.del_cmap, self.parm.cmap)
        # If we had anything to do, remake the parm
        self.parm.remake_parm()

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class fixtopology(Action):
    """
    This action will look through the lists of bonds, angles, and dihedrals both
    with and without hydrogen and make sure that all of the parameters are
    assigned to the correct lists. LEaP is known in certain cases to misassign
    angles and torsions if dummy atoms are involved. This action will correct
    the parameter assignments if any mistake was made
    """
    supported_classes = ('ChamberParm', 'AmberParm')

    def init(self, arg_list):
        # First check to see if there is any problem
        self.needs_fixing = False
        # Bond Fix, Angle Fix, Dihedral Fix -- variables to determine which
        # specific lists need fixing
        self.bf = self.af = self.df = False
        for bnd in self.parm.bonds_inc_h:
            if bnd.atom1.element != 1 and bnd.atom2.element != 1:
                self.needs_fixing = self.bf = True
                break
        for bnd in self.parm.bonds_without_h:
            if bnd.atom1.element == 1 or bnd.atom2.element == 1:
                self.needs_fixing = self.bf = True
                break
        for ang in self.parm.angles_inc_h:
            if (ang.atom1.element != 1 and ang.atom2.element != 1 and
                ang.atom3.element != 1):
                self.needs_fixing = self.af = True
                break
        for ang in self.parm.angles_without_h:
            if (ang.atom1.element == 1 or ang.atom2.element == 1 or
                ang.atom3.element == 1):
                self.needs_fixing = self.af = True
                break
        for dih in self.parm.dihedrals_inc_h:
            if (dih.atom1.element != 1 and dih.atom2.element != 1 and
                dih.atom3.element != 1 and dih.atom4.element != 1):
                self.needs_fixing = self.df = True
                break
        for dih in self.parm.dihedrals_without_h:
            if (dih.atom1.element == 1 or dih.atom2.element == 1 or
                dih.atom3.element == 1 or dih.atom4.element == 1):
                self.needs_fixing = self.df = True
                break

    def __str__(self):
        if self.needs_fixing:
            return 'Fixing bond/angle/dihedral list assignments'
        return 'No bond/angle/dihedral list problems detected. Doing nothing.'

    def execute(self):
        if not self.needs_fixing: return
        # This is the tracked list type we're using
        listtype = type(self.parm.bonds_inc_h)
        if self.bf:
            # Need to fix bonds
            bonds_inc_h = listtype()
            bonds_without_h = listtype()
            for bnd in self.parm.bonds_inc_h + self.parm.bonds_without_h:
                if bnd.atom1.element == 1 or bnd.atom2.element == 1:
                    bonds_inc_h.append(bnd)
                else:
                    bonds_without_h.append(bnd)
            self.parm.bonds_inc_h = bonds_inc_h
            self.parm.bonds_without_h = bonds_without_h
        if self.af:
            # Need to fix angles
            angles_inc_h = listtype()
            angles_without_h = listtype()
            for ang in self.parm.angles_inc_h + self.parm.angles_without_h:
                if (ang.atom1.element == 1 or ang.atom2.element == 1 or
                    ang.atom3.element == 1):
                    angles_inc_h.append(ang)
                else:
                    angles_without_h.append(ang)
            self.parm.angles_inc_h = angles_inc_h
            self.parm.angles_without_h = angles_without_h
        if self.df:
            # Need to fix dihedrals
            dihedrals_inc_h = listtype()
            dihedrals_without_h = listtype()
            for dih in self.parm.dihedrals_inc_h+self.parm.dihedrals_without_h:
                if (dih.atom1.element == 1 or dih.atom2.element == 1 or
                    dih.atom3.element == 1 or dih.atom4.element == 1):
                    dihedrals_inc_h.append(dih)
                else:
                    dihedrals_without_h.append(dih)
            self.parm.dihedrals_inc_h = dihedrals_inc_h
            self.parm.dihedrals_without_h = dihedrals_without_h
        self.parm.remake_parm()

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class chamber(Action):
    """
    This action will read CHARMM parameter, topology (RTF), and stream (STR)
    files and apply those parameters to a structure defined in a CHARMM PSF
    (protein structure file). XPLOR, CHARMM, and VMD-generated PSF files are all
    supported. You may specify -top, -param, and -str as many times as you want
    to provide multiple parameter files. All topology files are read first (in
    the order they are specified), followed by all parameter files, and finally
    all stream files are read in.  Topology files are only necessary if the
    parameter files do not define the atom types (newer CHARMM force fields
    define atom types directly in the parameter file.

    Options:
        -top        CHARMM topology, or Residue Topology File (RTF) file
        -param      CHARMM parameter file
        -str        CHARMM stream file. Only RTF and parameter sections are read
        -psf        CHARMM PSF file
        -crd        Input coordinate file (must be PDB file currently)
        -nocmap     Do not use any CMAP parameters
        usechamber  Use the 'chamber' program to write a topology file instead
        -box        Box dimensions. If no angles are defined, they are assumed
                    to be 90 degrees (orthorhombic box)

    If the PDB file has a CRYST1 record, the box information will be set from
    there. Any box info given on the command-line will override the box found in
    the crd file.
    """
    needs_parm = False

    def init(self, arg_list):
        from os.path import expandvars, expanduser
        # First get all of the topology files
        self.topfiles, self.paramfiles, self.streamfiles = [], [], []
        while arg_list.has_key('-top', mark=False):
            topfile = expanduser(arg_list.get_key_string('-top', None))
            topfile = expandvars(topfile)
            if not os.path.exists(topfile):
                raise FileDoesNotExist('CHARMM RTF file %s does not exist' %
                                       topfile)
            self.topfiles.append(topfile)
        while arg_list.has_key('-param', mark=False):
            param = expanduser(arg_list.get_key_string('-param', None))
            param = expandvars(param)
            if not os.path.exists(param):
                raise FileDoesNotExist('CHARMM parameter file %s does not exist'
                                       % param)
            self.paramfiles.append(param)
        while arg_list.has_key('-str', mark=False):
            stream = expanduser(arg_list.get_key_string('-str', None))
            stream = expandvars(stream)
            if not os.path.exists(stream):
                raise FileDoesNotExist('CHARMM stream file %s does not exist' %
                                       stream)
            self.streamfiles.append(stream)
        crdfile = arg_list.get_key_string('-crd', None)
        if crdfile is not None:
            self.crdfile = expanduser(expandvars(crdfile))
            if not os.path.exists(self.crdfile):
                raise FileDoesNotExist('Coordinate file %s does not exist' %
                                       self.crdfile)
        else:
            self.crdfile = None
        self.cmap = not arg_list.has_key('-nocmap')
        self.usechamber = arg_list.has_key('usechamber')
        psf = arg_list.get_key_string('-psf', None)
        if psf is None:
            raise InputError('chamber requires a PSF file')
        self.psf = expanduser(expandvars(psf))
        if not os.path.exists(self.psf):
            raise InputError('chamber PSF file %s cannot be found' % self.psf)
        box = arg_list.get_key_string('-box', None)

        if box is not None:
            try:
                self.box = [float(w) for w in box.replace(',', ' ').split()]
            except ValueError:
                raise InputError('Box info must be comma-delimited floats')
            if len(self.box) == 3:
                self.box += [90.0, 90.0, 90.0]
            elif len(self.box) != 6:
                raise InputError('Box must be 3 lengths or 3 lengths and 3 '
                                 'angles!')
        else:
            self.box = None

        # Make sure we have legal input
        if not self.paramfiles and not self.streamfiles:
            raise InputError('No parameter files were provided!')
        if self.usechamber:
            # We need a topfile, paramfile, and crd file for chamber
            if not self.topfiles or not self.paramfiles:
                raise InputError('The chamber program requires both RTF '
                                 '(topology) and PAR (parameter) files.')
            if not self.crdfile:
                raise InputError('The chamber program requires a CRD file.')

    def __str__(self):
        retstr = 'Creating chamber topology file from PSF %s, ' % self.psf
        if self.topfiles:
            retstr += 'RTF files [%s], ' % (', '.join(self.topfiles))
            if not self.streamfiles or not self.paramfiles:
                retstr += 'and '
        if self.paramfiles:
            retstr += 'PAR files [%s]' % (', '.join(self.paramfiles))
            if self.streamfiles:
                retstr += ', and '
        if self.streamfiles:
            retstr += 'STR files [%s].' % (', '.join(self.streamfiles))
        if self.crdfile is not None:
            retstr += ' Coords from %s.' % self.crdfile
        if self.cmap:
            retstr += ' Using CMAP.'
        else:
            retstr += ' NO CMAP.'
        if self.box is not None:
            retstr += ' Box info %s.' % (self.box)
        if self.usechamber:
            retstr += ' Using chamber program.'
        return retstr

    def execute(self):
        from chemistry.amber._chamberparm import ConvertFromPSF
        from chemistry.charmm.parameters import CharmmParameterSet
        from chemistry.charmm.psf import CharmmPsfFile
        from chemistry.system import ChemicalSystem
        from subprocess import Popen
        import tempfile
        if self.usechamber:
            tmpprm = tempfile.mktemp(suffix='.parm7')
            tmpcrd = tempfile.mktemp(suffix='.rst7')
            args = ['chamber', '-top', self.topfiles[0], '-param',
                    self.paramfiles[0], '-psf', self.psf, '-crd', self.crdfile,
                    '-str']
            args += self.topfiles[1:] + self.paramfiles[1:]
            args += self.streamfiles[:]
            args += ['-p', tmpprm, '-inpcrd', tmpcrd]
            if self.cmap:
                args.append('-cmap')
            else:
                args.append('-nocmap')
            if self.box is not None:
                args += ['-box'] + [str(x) for x in self.box[:3]]
            process = Popen(args)
            if process.wait() != 0:
                raise ChamberError('chamber failed with command\n\t' +
                                   ' '.join(args))
            # Now we have created a topology file
            self.parm_list.add_parm(tmpprm, tmpcrd)
            return
        # We're not using chamber, do the conversion in-house
        try:
            parmset = CharmmParameterSet()
            for tfile in self.topfiles:
                parmset.read_topology_file(tfile)
            for pfile in self.paramfiles:
                parmset.read_parameter_file(pfile)
            for sfile in self.streamfiles:
                parmset.read_stream_file(sfile)
            # All parameters are loaded, now condense the types
            parmset.condense()
        except ChemError:
            raise ChamberError('Problem reading CHARMM parameter sets')

        # Now read the PSF
        try:
            psf = CharmmPsfFile(self.psf)
        except ChemError:
            raise ChamberError('Problem reading CHARMM PSF')

        # Read the PDB and set the box information
        if self.crdfile is not None:
            crd = ChemicalSystem.load_from_pdb(self.crdfile)
            # Set the box info from self.box if set
            if self.box is None and crd.box is not None:
                psf.set_box(crd.box[:])
            elif self.box is not None:
                psf.set_box(self.box[:])
            else:
                psf.set_box(None)
        else:
            # Set the box information
            psf.set_box(self.box)

        nsets = len(parmset.parametersets)
        frcfield = '%2d' % nsets
        frcfield += ('\n%2d' % nsets).join(parmset.parametersets)
        # Delete the CMAP list if requested
        if not self.cmap:
            try:
                del psf.cmap_list
            except AttributeError:
                pass
        # Now load the parameters
        try:
            psf.load_parameters(parmset)
        except ChemError:
            raise ChamberError('Problem assigning parameters to PSF')
        parm = ConvertFromPSF(psf, frcfield).view(ChamberParm)
        if self.crdfile is not None:
            parm.load_coordinates(crd.positions())
        parm.prm_name = self.psf
        self.parm_list.add_parm(parm)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class minimize(Action):
    """
    This action takes a structure and minimizes the energy using OpenMM.
    Following this action, the coordinates stored in the topology will be the
    minimized structure
    
    General Options:
        cutoff <cutoff>     Nonbonded cutoff in Angstroms
        restrain <mask>     Selection of atoms to restrain
        weight <k>          Weight of positional restraints (kcal/mol/A^2)
        norun               Do not run the calculation
        script <script>     Name of the Python script to write to run this
                            calculation

    Implicit Solvent options:
        igb <IGB>           GB model to use (0=No GB, 1,2,5,7,8 correspond to
                            Amber models)
        saltcon <conc>      Salt concentration for GB in Molarity

    The implicit solvent options will be ignored for systems with periodic
    boundary conditions
    """
    supported_classes = ('ChamberParm', 'AmberParm')

    def init(self, arg_list):
        self.cutoff = arg_list.get_key_float('cutoff', None)
        mask = arg_list.get_key_mask('restrain', None)
        self.igb = arg_list.get_key_int('igb', 0)
        self.saltcon = arg_list.get_key_float('saltcon', 0.0)
        self.weight = arg_list.get_key_float('weight', 0.0)
        self.norun = arg_list.has_key('norun')
        self.script = arg_list.get_key_string('script', None)
        self.platform = arg_list.get_key_string('platform', None)
        self.precision = arg_list.get_key_string('precision', 'mixed')
        # Check for legal values
        if self.parm.ptr('ifbox') == 0:
            if self.cutoff is None or self.cutoff > 500:
                self.cutoff = None # no cutoff
        elif self.cutoff is None:
            self.cutoff = 8.0
        elif self.cutoff < 7:
            raise InputError('Cutoff unreasonably small.')
        if self.parm.ptr('ifbox') == 0 and self.saltcon < 0:
            raise InputError('Salt concentration must be non-negative')
        if mask is not None:
            self.restrain = AmberMask(self.parm, mask)
            if self.weight <= 0:
                raise InputError('Restraints require a restraint stiffness '
                                 'larger than 0 kcal/mol/A^2')
        else:
            self.restrain = None
        if self.platform not in ('CUDA', 'OpenCL', 'CPU', 'Reference', None):
            raise InputError('platform must be CUDA, OpenCL, CPU or Reference '
                             '(NOT %s)' % self.platform)
        if self.precision not in ('mixed', 'single', 'double'):
            raise InputError('precision must be single, double, or mixed.')
        if self.parm.ptr('ifbox') == 0 and not self.igb in (0, 1, 2, 5, 7, 8):
            raise InputError('Illegal igb value (%d). Must be 0, 1, 2, 5, 7, '
                             'or 8' % self.igb)

    def __str__(self):
        retstr = 'Minimizing %s ' % self.parm
        if self.parm.ptr('ifbox'):
            retstr += 'with PME '
        elif self.igb == 0:
            retstr += 'in gas phase '
        elif self.igb == 1:
            retstr += 'with GB(HCT) '
        elif self.igb == 2:
            retstr += 'with GB(OBC1) '
        elif self.igb == 5:
            retstr += 'with GB(OBC2) '
        elif self.igb == 7:
            retstr += 'with GB(GBneck) '
        elif self.igb == 8:
            retstr += 'with GB(GBneck2) '
        if self.cutoff is None:
            retstr += 'and no cutoff. '
        else:
            retstr += 'and a cutoff of %.2f Angstroms. ' % self.cutoff
        if self.restrain is not None:
            retstr += 'Restraining %s with weights %f. ' % (self.restrain,
                                                            self.weight)
        return retstr.strip()

    def execute(self):
        from ParmedTools.simulations.openmm import minimize, HAS_OPENMM
        if not HAS_OPENMM:
            raise SimulationError('OpenMM could not be imported. Skipping.')

        if not hasattr(self.parm, 'coords'):
            raise SimulationError('You must load coordinates before "minimize"')
        minimize(self.parm, self.igb, self.saltcon, self.cutoff,
                 self.restrain, self.weight, self.script, self.platform,
                 self.precision, self.norun)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

