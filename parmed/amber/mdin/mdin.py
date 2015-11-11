"""
This module will create an amber mdin file for either sander or      
pmemd (or others). The program specification loads the appropriate   
dictionaries with default values, etc. It can read and write mdins. 

                           GPL LICENSE INFO                             

Copyright (C) 2009 - 2014 Dwight Mcgee, Bill Miller III, and Jason Swails

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""

# This module will create and read a sander/pmemd input
from parmed.amber.mdin.cntrl import cntrl
from parmed.amber.mdin.ewald import ewald
from parmed.amber.mdin.pb import pb
from parmed.amber.mdin.qmmm import qmmm
from parmed.exceptions import InputError
from parmed.utils.six import string_types
from parmed.utils.six.moves import range

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def addOn(line, string, file):
    if len(line.strip()) == 0:
        return line + string
    elif len(line) + len(string) > 40:
        file.write(line + '\n')
        return ' ' + string
    else:
        return line + string

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Mdin(object):
   
   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def __init__(self, program = 'sander', verbosity = 1):
        # define instance data
        self.program = program   # which program we're creating the input file for
        self.cntrl_obj = cntrl() # object with cntrl namelist vars in a dictionary
        self.ewald_obj = ewald() # object with ewald namelist vars in a dictionary
        self.pb_obj = pb()       # object with pb namelist vars in a dictionary
        self.qmmm_obj = qmmm()   # object with qmmm namelist vars in a dictionary
        self.verbosity = 0       # verbosity level: 0 -- print nothing
                                 #                  1 -- print errors
                                 #                  2 -- 1 + warnings
                                 #                  3 -- 2 + notes
        self.cards = []          # array that has all of the input cards that come
                                 # after namelists
        self.cntrl_nml = {}      # dictionary with cntrl namelist vars
        self.cntrl_nml_defaults = {} # dictionary with default cntrl namelist vars
        self.ewald_nml = {}          # dictionary with ewald namelist vars
        self.ewald_nml_defaults = {} # dictionary with default ewald namelist vars
        self.pb_nml = {}             # dictionary with pb namelist vars
        self.pb_nml_defaults = {}    # dictionary with default pb namelist vars
        self.qmmm_nml = {}           # dictionary with qmmm namelist vars
        self.qmmm_nml_defaults = {}  # dictionary with default qmmm namelist vars
        self.valid_namelists = []    # array with valid namelists for each program
        self.title = 'mdin prepared by mdin.py'   # title for the mdin file


        if self.program == "sander":
            self.cntrl_nml = self.cntrl_obj.sander
            self.ewald_nml = self.ewald_obj.sander
            self.pb_nml = self.pb_obj.sander
            self.qmmm_nml = self.qmmm_obj.sander
            self.valid_namelists = ['cntrl','ewald','qmmm','pb']
        elif self.program == "sander.APBS":
            self.cntrl_nml = self.cntrl_obj.sander
            self.pb_nml = self.pb_obj.sanderAPBS
            self.valid_namelists = ['cntrl','apbs']
        elif self.program == "pmemd":
            self.cntrl_nml = self.cntrl_obj.pmemd
            self.ewald_nml = self.ewald_obj.pmemd
            self.valid_namelists = ['cntrl','ewald']
        else:
            raise InputError('Error: program (%s) unrecognized!' %
                                   self.program)

        self.cntrl_nml_defaults = self.cntrl_nml.copy()
        self.ewald_nml_defaults = self.ewald_nml.copy()
        self.pb_nml_defaults = self.pb_nml.copy()
        self.qmmm_nml_defaults = self.qmmm_nml.copy()

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def write(self, filename = 'mdin'):
      
        # open the file for writing and write the header and &cntrl namelist
        file = open(filename,'w')
        file.write(self.title + '\n')
        file.write('&cntrl\n')
        # automatic indent of single space
        line = ' '
        # add any variable that is different from the default to the mdin file
        for var in self.cntrl_nml.keys():
            if self.cntrl_nml[var] == self.cntrl_nml_defaults[var]: continue
            if isinstance(self.cntrl_nml[var], string_types):
                line = addOn(line, "%s='%s', " % (var, self.cntrl_nml[var]), file)
            else:
                line = addOn(line, "%s=%s, " % (var, self.cntrl_nml[var]), file)

        # flush any remaining items that haven't yet been printed to the mdin file
        if len(line.strip()) != 0:
            file.write(line + '\n')

        # end the namelist
        file.write('/\n')

        # print the ewald namelist if any variables differ from the default
        line = ' '
        has_been_printed = False  # keep track if this namelist has been printed
        for var in self.ewald_nml.keys():
            if self.ewald_nml[var] == self.ewald_nml_defaults[var]: continue
            if not has_been_printed:
                file.write('&ewald\n')
                has_been_printed = True
            if isinstance(self.ewald_nml_defaults[var], string_types):
                line = addOn(line, "%s='%s', " % (var, self.ewald_nml[var]), file)
            else:
                line = addOn(line, "%s='%s', " % (var, self.ewald_nml[var]), file)

        # flush any remaining items that haven't been printed to the mdin file
        if len(line.strip()) != 0:
            file.write(line + '\n')

        # end the namelist
        if has_been_printed:
            file.write('/\n')

        # print the pb namelist if any variables differ from the original
        line = ' '
        has_been_printed = False # keep track if this namelist has been printed
        for var in self.pb_nml.keys():
            if self.pb_nml[var] == self.pb_nml_defaults[var]: continue
            if not has_been_printed:
                if self.program == 'sander.APBS':
                    file.write('&apbs\n')
                else:
                    file.write('&pb\n')
                has_been_printed = True
            if isinstance(self.pb_nml[var], string_types):
                line = addOn(line,"%s='%s', " % (var, self.pb_nml[var]), file)
            else:
                line = addOn(line,"%s=%s, " % (var, self.pb_nml[var]), file)

        # flush any remaining items that haven't been printed to the mdin file
        if len(line.strip()) != 0:
            file.write(line + '\n')

        # end the namelist
        if has_been_printed:
            file.write('/\n')

        # print the qmmm namelist if any variables differ from the original
        line = ' '
        has_been_printed = False # keep track if this namelist has been printed
        if self.cntrl_nml['ifqnt'] == 1:
            for var in self.qmmm_nml.keys():
                if self.qmmm_nml[var] == self.qmmm_nml_defaults[var]: continue
                if not has_been_printed:
                    file.write('&qmmm\n')
                    has_been_printed = True
                if isinstance(self.qmmm_nml_defaults[var], string_types):
                    line = addOn(line, "%s='%s', " % (var, self.qmmm_nml[var]),
                                 file)
                else:
                    line = addOn(line, "%s=%s, " % (var, self.qmmm_nml[var]), file)
      
            # flush any remaining items that haven't been printed to the mdin file
            if len(line.strip()) != 0:
                file.write(line + '\n')

            # end the namelist
            if has_been_printed:
                file.write('/\n')

        # Write the cards to the input file
        for i in range(len(self.cards)):
            file.write(self.cards[i].strip() + '\n')
        if len(self.cards) != 0:
            file.write('END\n')
      
        file.close()

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
   
    def read(self, filename = 'mdin'):
        lines = open(filename, 'r').readlines()

        # split up input file into separate fields by comma
        blocks = []  # namelists in the order they appear
        block_fields = [] # array of arrays that correspond to entries in 
                          # namelists found in "blocks" above
        inblock = False
        lead_comment = True
        for i in range(len(lines)):
            if not inblock and not lines[i].strip().startswith('&') and \
                    lead_comment:
                continue
            elif not inblock and not lines[i].strip().startswith('&') and \
                    not lead_comment:
                final_ended = True
                for j in range(i,len(lines)):
                    if lines[j].strip().startswith('&'):
                        final_ended = False
                if final_ended and len(lines[i].strip()) != 0:
                    self.cards.append(lines[i])
            elif not inblock and lines[i].strip().startswith('&'):
                lead_comment = False
                inblock = True
                block = lines[i].strip()[1:].lower()
                blocks.append(block)    # add the name of the namelist to "blocks"
                block_fields.append([]) # add empty array to be filled with entries 
                                        # for given namelist
                if not block in self.valid_namelists:
                    raise InputError('Invalid namelist (%s) in input '
                            'file (%s) for %s' % (lines[i].strip(), filename,
                            self.program))
            elif inblock and (lines[i].strip() == '/' or 
                            lines[i].strip() == '&end'):
                inblock = False
            elif inblock and lines[i].strip().startswith('&'):
                raise InputError('Invalid input file (%s). Terminate each '
                                'namelist before another is started' % filename)
            elif inblock:
                items = lines[i].strip().split(',')
                j = 0
                while j < len(items):
                    items[j] = items[j].strip()
                    if len(items[j]) == 0:
                        items.pop(j)
                    else:
                        j += 1
                block_fields[len(block_fields)-1].extend(items)

        # take out the last END in the cards if it's there
        if len(self.cards) != 0 and \
                self.cards[len(self.cards)-1].strip().upper() == 'END':
            self.cards.pop()

        # combine any multi-element fields: e.g. rstwt=1,2,3,
        begin_field = -1
        for i in range(len(block_fields)):
            for j in range(len(block_fields[i])):
                if not '=' in block_fields[i][j]:
                    if begin_field == -1:
                        raise InputError('Invalid input file (%s).' %
                                               filename)
                    else:
                        block_fields[i][begin_field] += ',' + block_fields[i][j]
                else:
                    begin_field = j

        # now parse through the options and add them to the dictionaries
        for i in range(len(block_fields)):
            for j in range(len(block_fields[i])):
                if not '=' in block_fields[i][j]:
                    continue
                else:
                    var = block_fields[i][j].split('=')
                    self.change(blocks[i], var[0].strip(), var[1].strip())

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def change(self, namelist, variable, value):
        """ Change the value of a variable without adding a new key-pair """
      
        variable = variable.lower()
        if isinstance(value, string_types):
            if (value.startswith('"') and value.endswith('"')) or (
                   value.startswith("'") and value.endswith("'")):
                value = value[1:-1]
      
        if namelist == "cntrl":
            if variable in self.cntrl_nml.keys():
                mytype = type(self.cntrl_nml_defaults[variable])
                self.cntrl_nml[variable] = mytype(value)
            else:
                raise InputError('Unknown variable (%s) in &cntrl!' % variable)
        elif namelist == 'ewald': 
            if variable in self.ewald_nml.keys():
                mytype = type(self.ewald_nml_defaults[variable])
                self.ewald_nml[variable] = mytype(value)
            else:
                raise InputError('Unknown variable (%s) in &ewald!' % variable)
        elif namelist == 'pb' or namelist == 'apbs':
            if variable in self.pb_nml.keys():
                mytype = type(self.pb_nml_defaults[variable])
                self.pb_nml[variable] = mytype(value)
            else:
                raise InputError('Unknown variable (%s) in &%s!' % 
                                (variable, namelist))
        elif namelist == 'qmmm':
            if variable in self.qmmm_nml.keys():
                mytype = type(self.qmmm_nml_defaults[variable])
                self.qmmm_nml[variable] = mytype(value)
            else:
                raise InputError('Unknown variable (%s) in &qmmm' % variable)
        else:
            raise InputError('Unknown namelist (%s)!' % namelist)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def check(self):
        return True

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def SHAKE(self):
        self.change('cntrl','ntf', 2)
        self.change('cntrl','ntc', 2)
        self.change('cntrl','dt', 0.002)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def constPressure(self, press=1.0, taup=1.0):
        self.change('cntrl','ntb', 2)
        self.change('cntrl','ntp', 1)
        self.change('cntrl','pres0', press)
        self.change('cntrl','taup', taup)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def constVolume(self):
        self.change('cntrl','ntb', 1)
        self.change('cntrl','ntp', 0)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def constTemp(self, ntt=3, temp=300.0, gamma_ln=2.0, ig=-1, tautp=1.0):
        self.change('cntrl','ntt', ntt)
        self.change('cntrl','temp0', temp)
        self.change('cntrl','tempi', temp)
        self.change('cntrl','gamma_ln', gamma_ln if ntt==3 else 0)
        self.change('cntrl','ig', ig)
        self.change('cntrl','tautp', tautp)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def constpH(self, solvph=7.0, igb=2, ntcnstph=10):
        self.change('cntrl','icnstph', 1)
        self.change('cntrl','solvph', solvph)
        self.change('cntrl','ntcnstph', ntcnstph)
        self.change('cntrl','igb', igb)
        self.change('cntrl','ntb', 0)
        self.change('cntrl','saltcon', 0.1)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def restrainHeavyAtoms(self, restraint_wt=0.0):
        self.change('cntrl','ntr', 1)
        self.change('cntrl','restraint_wt', restraint_wt)
        self.change('cntrl','restraintmask', '!@H=')

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def restrainBackbone(self, restraint_wt=0.0):
        self.change('cntrl','ntr', 1)
        self.change('cntrl','restraint_wt', restraint_wt)
        self.change('cntrl','restraintmask', '@N,CA,C')

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def genBorn(self, igb=5, rgbmax=25.0):
        self.change('cntrl','igb', igb)
        self.change('cntrl','ntb', 0)
        self.change('cntrl','ntp', 0)
        self.change('cntrl','rgbmax', rgbmax)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def time(self, time=1000.0, dt=-1): # time in ps
        if dt == -1:
            if self.cntrl_nml['ntc'] == 2 and self.cntrl_nml['ntf'] == 2:
                dt = 0.002
            else:
                dt = 0.001
        time = int(time / dt)

        self.change('cntrl','dt', dt)
        self.change('cntrl','nstlim', time)
        self.change('cntrl','imin', 0)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def heat(self, tempi=0.0, temp0=300.0, ntt=3, tautp=5.0, ig=-1, gamma_ln=5.0):
        self.constVolume()
        self.change('cntrl','tempi', tempi)
        self.change('cntrl','temp0', temp0)
        self.change('cntrl','ntt', ntt)
        self.change('cntrl','tautp', tautp)
        self.change('cntrl','ig', ig)
        self.change('cntrl','gamma_ln', gamma_ln if ntt==3 else 0)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def restart(self,ntx=5):
        self.change('cntrl','irest',1)
        self.change('cntrl','ntx',ntx)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def TI(self, clambda=0.0):
        self.change('cntrl','clambda', clambda)
        self.change('cntrl','icfe',1)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def softcore_TI(self, scalpha=0.5, scmask='', crgmask='', logdvdl=0):
        self.change('cntrl','icfe',1)
        self.change('cntrl','ifsc',1)
        self.change('cntrl','scalpha',scalpha)
        self.change('cntrl','scmask',scmask)
        self.change('cntrl','crgmask',crgmask)
        self.change('cntrl','logdvdl',logdvdl)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def minimization(self, imin=1, maxcyc=1, ncyc=10, ntmin=1):
        self.change('cntrl','imin', imin)
        self.change('cntrl','maxcyc', maxcyc)
        self.change('cntrl','ncyc', ncyc)
        self.change('cntrl','ntmin', ntmin)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def AddCard(self, title='Residues in card', cardString='RES 1'):
        self.cards.append('%s\n%s\nEND' % (title, cardString))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
