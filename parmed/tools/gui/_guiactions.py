"""
This module contains all of the methods to invoke specific ParmEd actions listed
in that module. Each method will establish the necessary variables (StringVar),
create whatever window it needs to get information from the user (or give
information to the user), wait for that window to close, then dispatch that
Action to the class in ParmEd actions.

Follow the general trend if you wish to add your method to the GUI. Note, any
method that you want accessible through the GUI must have an action method put
here with the same name as the class found in ParmEd actions.
"""

from parmed.utils.six.moves import range
from parmed.utils.six.moves.tkinter import *
from parmed.utils.six.moves.tkinter_messagebox import (
        askyesno, showinfo, showerror)
from parmed.tools import actions
from parmed.tools.argumentlist import ArgumentList
from parmed.tools.gui.guifiletools import file_chooser, save_file_chooser
from parmed.tools.gui import _guiwidgets

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def writefrcmod(root, amber_prmtop, messages):
    """ Dumps an frcmod file to a given filename """
    fname = save_file_chooser('frcmod', '.frcmod')
    if fname: 
        action = actions.writeFrcmod(amber_prmtop, fname)
        action.execute()
        messages.write('%s\n' % action)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def loadrestrt(root, amber_prmtop, messages):
    """ Finds a file to load as a restart """
    type_list = [('Inpcrd', '*.inpcrd'), ('Inpcrd', '*.crd'),
                 ('Restart', '*.restrt'), ('Restart', '*.rst7'),
                 ('All Files', '*')]
    fname = file_chooser('Amber Coordinate File', type_list)
    if fname: 
        action = actions.loadRestrt(amber_prmtop, ArgumentList(fname))
        messages.write('%s\n' % action)
        action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def writeoff(root, amber_prmtop, messages):
    """ Dumps an OFF library to a given filename """
    fname = save_file_chooser('OFF library', '.lib')
    if fname: 
        action = actions.writeOFF(amber_prmtop, ArgumentList(fname))
        messages.write('%s\n' % action)
        action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def changeradii(root, amber_prmtop, messages):
    """ Allows users to change the GB Radius set """
    title = 'Choose a Radius Set'
    desc = ('Select a radius set for implicit solvent\n' +
            'calculations. This has the same effect as\n' +
            '"set default PBRadii <value>" in tleap')
    radius_selection = StringVar()
    namelist = ['bondi', 'mbondi', 'mbondi2', 'mbondi3', 'amber6']
    cmd_window = _guiwidgets.RadioButtonWindow(
                    amber_prmtop, title, desc, radius_selection, namelist)
    # Wait until the window is destroyed, then get the variable and pass it
    # over to the class
    cmd_window.wait_window()
    sel = str(radius_selection.get())
    if sel:
        action = actions.changeRadii(amber_prmtop, ArgumentList(sel))
        messages.write('%s\n' % action)
        action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def changeljpair(root, amber_prmtop, messages):
    """ Changes a pair-wise LJ interaction for pre-combined epsilon/Rmin """
    # The variables we need for changeljpair
    widget_list = [('MaskEntry', 'Atom(s) Type 1 Mask'),
                   ('MaskEntry', 'Atom(s) Type 2 Mask'),
                   ('Entry', 'Combined Radius'),
                   ('Entry', 'Combined Well Depth')]
    # Variable list -- we need 2 masks and 2 floats
    var_list = [StringVar(), StringVar(), StringVar(), StringVar()]
    # description
    description = ' '.join(actions.changeLJPair.__doc__.split())
    cmd_window = _guiwidgets.ActionWindow('changeLJPair', amber_prmtop,
                                widget_list, var_list, description)
    cmd_window.wait_window()
    # Make sure we didn't cancel (or just press OK with no input), or just leave
    vars_exist = True in [bool(v.get()) for v in var_list]
    if not vars_exist: return
    # Now that we did something, do it
    var_list = [v.get() for v in var_list]
    try:
        action = actions.changeLJPair(amber_prmtop,ArgumentList(var_list))
        messages.write('%s\n' % action)
        action.execute()
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def outparm(root, amber_prmtop, messages):
    """ Output a final topology file """
    fname = [save_file_chooser('prmtop', '.prmtop')]
    if amber_prmtop.parm.coords is not None and fname[0]:
        fname.append(save_file_chooser('inpcrd', '.inpcrd'))
    if fname[0]:
        action = actions.outparm(amber_prmtop, ArgumentList(fname))
        messages.write('%s\n' % action)
        action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printflags(root, amber_prmtop, messages):
    """ Prints all of the flags in the topology file """
    window = Toplevel(root)
    window.resizable(True, True)
    window.title('%%FLAG list in %s' % amber_prmtop.parm)
    text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                         pady=5, width=80, height=20)
    text.pack(fill=BOTH, expand=1)
    action = actions.printFlags(amber_prmtop, ArgumentList(''))
    text.write(action)
    window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printpointers(root, amber_prmtop, messages):
    """ Prints all of the flags in the topology file """
    window = Toplevel(root)
    window.resizable(True, True)
    window.title('POINTER list in %s' % amber_prmtop.parm)
    text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                         pady=5, width=80, height=20)
    text.pack(fill=BOTH, expand=1)
    action = actions.printPointers(amber_prmtop, ArgumentList(''))
    text.write(action)
    window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def changelj14pair(root, amber_prmtop, messages):
    """ Changes specific 1-4 Lennard Jones pairs """
    # Only good for chamber topologies
    if not amber_prmtop.parm.chamber:
        showerror('Incompatible',
                  'changeLJ14Pair is only valid for chamber topologies!')
        return
    # variables we need for changelj14pair
    widget_list = [('MaskEntry', 'Atom(s) Type 1 Mask'),
                   ('MaskEntry', 'Atom(s) Type 2 Mask'),
                   ('Entry', 'Combined Radius'),
                   ('Entry', 'Combined Well Depth')]
    # Variable list -- we need 2 masks and 2 floats
    var_list = [StringVar(), StringVar(), StringVar(), StringVar()]
    # description
    description = ' '.join(actions.changeLJ14Pair.__doc__.split())
    cmd_window = _guiwidgets.ActionWindow('changeLJ14Pair', amber_prmtop,
                                widget_list, var_list, description)
    cmd_window.wait_window()
    # Make sure we didn't cancel (or just press OK with no input), or just leave
    vars_exist = True in [bool(v.get()) for v in var_list]
    if not vars_exist: return
    # Now that we did something, do it
    var_list = [v.get() for v in var_list]
    try:
        action=actions.changeLJ14Pair(amber_prmtop,ArgumentList(var_list))
        messages.write('%s\n' % action)
        action.execute()
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def checkvalidity(root, amber_prmtop, messages):
    """ Basic validity checks """
    # Create our Info window
    window = Toplevel(root)
    window.resizable(True, True)
    window.title('Problems with %s' % amber_prmtop.parm)
    text = _guiwidgets.ExitingScrollText(window, None, spacing3=2, padx=5,
                                         pady=5, width=81, height=30)
    text.pack(fill=BOTH, expand=1)
    # Set this text to catch the output of our action
    actions.checkValidity.output = text
    # Initialize our action
    action = actions.checkValidity(amber_prmtop, ArgumentList(''))
    messages.write('%s\n' % action)
    action.execute()
    text.write(action)

    window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def change(root, amber_prmtop, messages):
    """ Allows us to change a specific atomic property """
    # The spinbox is sent with the Spinbox, label, and then a list of all of the
    # values to give to it
    widget_list = [('Spinbox', 'Property to change', 'CHARGE', 'MASS',
                    'RADII', 'SCREEN', 'ATOM_NAME', 'AMBER_ATOM_TYPE',
                    'ATOM_TYPE_INDEX', 'ATOMIC_NUMBER'),
                    ('MaskEntry', 'Atoms to change'),
                    ('Entry', 'New Value for Property')]
    # We need 3 string variables, then get the description
    var_list = [StringVar(), StringVar(), StringVar()]
    description = 'Changes the property of given atoms to a new value'
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('change', amber_prmtop,
                        widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    vars_found = True in [bool(v.get()) for v in var_list]
    if not vars_found: return
    # If we did, store them and pass it to the class
    var_list = [v.get() for v in var_list]
    try:
        action = actions.change(amber_prmtop, ArgumentList(var_list))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    action.execute()
    messages.write('%s\n' % action)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printinfo(root, amber_prmtop, messages):
    """ Prints all of the info in a given FLAG """
    # Set up the window
    # variables we need for printInfo
    widget_list = [('Entry', '%FLAG you want info from')]
    # Variable list -- we need a single string
    var_list = [StringVar()]
    # description
    description = ' '.join(actions.printInfo.__doc__.split())
    cmd_window = _guiwidgets.ActionWindow('printInfo', amber_prmtop,
                                widget_list, var_list, description)
    cmd_window.wait_window()
    # Make sure we didn't cancel (or just press OK with no input), or just leave
    var = var_list[0].get()
    if not var: return
    # Now that we did something, do it
    action = actions.printInfo(amber_prmtop, ArgumentList(var))
    if not action.found:
        showerror('Not Found!', '%%FLAG %s not found!' % var.upper())
        return

    window = Toplevel(root)
    window.resizable(True, True)
    window.title('%%FLAG %s Info in %s' % (var.upper(), amber_prmtop))
    text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                            pady=5, width=82, height=20)
    text.write(action)
    text.pack(fill=BOTH, expand=1)
    window.wait_window()

    messages.write('Wrote info for flag %s\n' % var.upper())

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def addljtype(root, amber_prmtop, messages):
    """ Turns given mask into a new LJ atom type """
    # We need a mask, new radius, new epsilon, and for chamber topologies,
    # a new radius-1-4 and epsilon-1-4. Don't add the latter ones until we
    # know if we have a chamber prmtop or not.
    widget_list = [('MaskEntry', 'Atoms to make new LJ Type'),
                   ('Entry', 'New Radius (default old radius)'),
                   ('Entry', 'New Depth (default old depth)')]
    # We need 5 string variables, then get the description. 
    var_list = [StringVar(), StringVar(), StringVar(), StringVar(), StringVar()]
    description=('Turns given mask into a new LJ atom type. Uses the radius\n'
                 'and well depth from the first atom type in <mask> if none\n'
                 'are provided.')
    if amber_prmtop.parm.chamber:
        widget_list += [('Entry', 'New Radius for 1-4 Terms'),
                        ('Entry', 'New Depth for 1-4 Terms')]
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('addLJType', amber_prmtop,
                        widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    vars_found = True in [bool(v.get()) for v in var_list]
    if not vars_found: return
    # addljtype expects any _non_specified variables to be None
    var_list = [v.get().strip() for v in var_list]
    kw_var_list = [var_list[0]]
    if var_list[1]: kw_var_list.extend(['radius', var_list[1]])
    if var_list[2]: kw_var_list.extend(['epsilon', var_list[2]])
    if amber_prmtop.parm.chamber and var_list[3]:
        kw_var_list.extend(['radius_14', var_list[3]])
    if amber_prmtop.parm.chamber and var_list[4]:
        kw_var_list.extend(['epsilon_14', var_list[4]])
    try:
        action = actions.addLJType(amber_prmtop,ArgumentList(kw_var_list))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    action.execute()
    messages.write('%s\n' % action)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def setmolecules(root, amber_prmtop, messages):
    """ Sets the molecules. Asks users if they want ions to be solute or not """
    title = 'Set Molecularity'
    desc = ('This will determine the molecular topology (ATOMS_PER_MOLECULE\n'
            'and SOLVENT_POINTERS). Do you want the ions to be considered\n'
            'part of the solute?')
    solute_ions = StringVar()
    namelist = ['Yes', 'No']
    cmd_window = _guiwidgets.RadioButtonWindow(
                            amber_prmtop, title, desc, solute_ions, namelist)
    # Wait until the window is destroyed, then get the variable and pass it
    # over to the class
    cmd_window.wait_window()
    sel = str(solute_ions.get())
    if not sel: return

    if sel == 'Yes':
        sel = 'True'
    else:
        sel = 'False'
    action = actions.setMolecules(amber_prmtop, 
                                        ArgumentList('solute_ions '+sel))
    messages.write('%s\n' % action)
    action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printdetails(root, amber_prmtop, messages):
    """ Prints details about a given Amber mask """
    title = 'Print Details'
    mask = StringVar()
    cmd_window = Toplevel(root)
    cmd_window.title(title)
    mask_entry = _guiwidgets.MaskEntry(cmd_window, amber_prmtop,
                                'Input an Amber Mask', mask, cmd_window)
    mask_entry.config(pady=10)
    mask_entry.grid(row=0, column=0, sticky=N+E+S+W)
    button = Button(cmd_window, text='OK / Quit', command=cmd_window.destroy)
    button.grid(row=1, column=0, sticky=N+E+S+W)
    cmd_window.wait_window()
    if not mask.get(): return
    # Print our mask
    window = Toplevel(root)
    window.resizable(True, True)
    window.title('Atom information for mask %s' % mask.get())
    text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                         pady=5, width=100, height=20)
    text.pack(fill=BOTH, expand=1)
    action = actions.printDetails(amber_prmtop, mask.get())
    text.write(action)
    messages.write('Printed Amber Mask details on [%s]\n' % mask.get())
    window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def netcharge(root, amber_prmtop, messages):
    """ Calculates the net charge and shows its value """
    mask = StringVar()
    cmd_window = Toplevel(root)
    mask_entry = _guiwidgets.MaskEntry(cmd_window, amber_prmtop,
                        'Mask from which to calculate charge', mask, cmd_window)
    mask_entry.config(pady=10)
    mask_entry.grid(row=0, column=0, sticky=N+E+S+W)
    button = Button(cmd_window, text='OK / Quit', command=cmd_window.destroy)
    button.grid(row=1, column=0, sticky=N+E+S+W)
    cmd_window.wait_window()
    if not mask.get(): return
    action = actions.netCharge(amber_prmtop, ArgumentList(mask.get()))
    chg = action.execute()
    showinfo('Net Charge', 'The net charge of [%s] is %.4f' % (mask.get(), chg))

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def strip(root, amber_prmtop, messages):
    """ Strips a mask from the topology file """
    # We need a mask, new radius, new epsilon, and for chamber topologies,
    # a new radius-1-4 and epsilon-1-4. Don't add the latter ones until we
    # know if we have a chamber prmtop or not.
    widget_list = [('MaskEntry', 'Atoms to strip from topology')]
    # We need 5 string variables, then get the description. 
    var_list = [StringVar()]
    description=('Strips the selected atoms from the topology file. All\n'
                 'remaining atoms and parameters remain unchanged. Any\n'
                 'parameters associated with stripped atoms are removed.')
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('strip', amber_prmtop,
                                        widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    var = var_list[0]
    if not var.get(): return
    try:
        action = actions.strip(amber_prmtop, ArgumentList(var.get()))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    messages.write('%s\n' % action)
    action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def addexclusions(root, amber_prmtop, messages):
    """ Adds atoms to other atoms' exclusion list """
    # We need 2 masks
    widget_list = [('MaskEntry', 'Atoms to add excluded atoms to'),
                   ('MaskEntry', 'Atoms to exclude from other mask')]
    # We have 2 mask variables
    var_list = [StringVar(), StringVar()]
    # Description
    description = ('Allows you to add arbitrary excluded atoms to exclusion\n'
                   'lists. This omits all non-bonded interactions in the '
                   'direct-space\ncalculation, but does not omit interactions '
                   'from adjacent cells in\nperiodic simulations')
    cmd_window = _guiwidgets.ActionWindow('addExclusions', amber_prmtop,
                                          widget_list, var_list, description)
    cmd_window.wait_window()
   
    # Bail out if we didn't get any variables
    if not True in [bool(v.get()) for v in var_list]: return
   
    var_list = [v.get() for v in var_list]
    try:
        act = actions.addExclusions(amber_prmtop, ArgumentList(var_list))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    act.execute()
    messages.write('%s\n' % act)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def changeprotstate(root, amber_prmtop, messages):
    # We need a mask and a state
    widget_list = [('MaskEntry', 'Residue to change protonation state'),
                   ('Entry', 'Protonation state to change to')]
    var_list = [StringVar(), StringVar()]
    description=('Changes the protonation state of a pH-active titratable residue that\n'
                 'can be treated with constant pH MD in Amber.')
    cmd_window = _guiwidgets.ActionWindow('changeProtState', amber_prmtop,
                        widget_list, var_list, description)
    cmd_window.wait_window()
   
    # Bail out if we didn't get any variables
    if not True in [bool(v.get()) for v in var_list]: return
   
    var_list = [v.get() for v in var_list]
    try:
        action = actions.changeProtState(amber_prmtop, 
                                               ArgumentList(var_list))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    action.execute()
    messages.write('%s\n' % action)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def changeredustate(root, amber_prmtop, messages):
    # We need a mask and a state
    widget_list = [('MaskEntry', 'Residue to change reduction state'),
                   ('Entry', 'Reduction state to change to')]
    var_list = [StringVar(), StringVar()]
    description=('Changes the reduction state of a redox-active titratable residue that\n'
                 'can be treated with constant redox potential MD in Amber.')
    cmd_window = _guiwidgets.ActionWindow('changeRedoxState', amber_prmtop,
                        widget_list, var_list, description)
    cmd_window.wait_window()
   
    # Bail out if we didn't get any variables
    if not True in [bool(v.get()) for v in var_list]: return
   
    var_list = [v.get() for v in var_list]
    try:
        action = actions.changeRedoxState(amber_prmtop, 
                                               ArgumentList(var_list))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    action.execute()
    messages.write('%s\n' % action)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def scee(root, amber_prmtop, messages):
    # We need a value
    widget_list = [('Entry', '1-4 Electrostatic Scaling Factor')]
    var_list = [StringVar()]
    description = 'Adjust the scaling factor for 1-4 electrostatic interactions'

    cmd_window = _guiwidgets.ActionWindow('scee', amber_prmtop,
                                          widget_list, var_list, description)
    cmd_window.wait_window()
   
    var = var_list[0].get()

    # Bail out if we didn't get any variables
    if not var: return

    try:
        action = actions.scee(amber_prmtop, ArgumentList(var))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    messages.write('%s\n' % action)
    action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def scnb(root, amber_prmtop, messages):
    # We need a value
    widget_list = [('Entry', '1-4 van der Waals Scaling Factor')]
    var_list = [StringVar()]
    description = 'Adjust the scaling factor for 1-4 van der Waals interactions'

    cmd_window = _guiwidgets.ActionWindow('scee', amber_prmtop,
                                          widget_list, var_list, description)
    cmd_window.wait_window()
   
    var = var_list[0].get()

    # Bail out if we didn't get any variables
    if not var: return

    try:
        action = actions.scnb(amber_prmtop, ArgumentList(var))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    messages.write('%s\n' % action)
    action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printljtypes(root, amber_prmtop, messages):
    """
    Prints all of the atoms that have the same LJ type as atoms in a given mask
    """
    # We need a mask
    widget_list = [('MaskEntry', 'Atom Mask or Atom Type Name')]
    var_list = [StringVar()]
    description = ('Get a list of all atoms that share a common LJ type')
    cmd_window = _guiwidgets.ActionWindow('printLJTypes', amber_prmtop,
                                          widget_list, var_list, description)
    cmd_window.wait_window()
   
    var = var_list[0].get()

    # Bail out if we didn't get any variables
    if not var: return
   
    # Instantiate our action
    try:
        action = actions.printLJTypes(amber_prmtop, ArgumentList(var))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return

    # Create our Info window
    window = Toplevel(root)
    window.resizable(True, True)
    window.title('LJ Type List')
    text = _guiwidgets.ExitingScrollText(window, None, spacing3=2, padx=5,
                                         pady=5, width=45, height=30)
    text.pack(fill=BOTH, expand=1)
    text.write(action)
    messages.write('Printed LJ types for [%s]\n' % var)
    window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def changeljsingletype(root, amber_prmtop, messages):
    """ Changes radius/well depth of a single LJ type given by the mask """
    # We need a mask, radius, and well depth
    widget_list = [('MaskEntry', 'Mask to change LJ Type'),
                   ('Entry', 'New LJ Radius'), ('Entry', 'New LJ Depth')]
    var_list = [StringVar(), StringVar(), StringVar()]
    description = "Change a given atom type's LJ radius and well depth"
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('changeLJSingleType', amber_prmtop,
                                          widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    vars_found = True in [bool(v.get()) for v in var_list]
    if not vars_found: return
    # addljtype expects any _non_specified variables to be None
    var_list = [v.get() for v in var_list]
    for i, v in enumerate(var_list):
        if not v: var_list[i] = None
    # If we did, store them and pass it to the class
    try:
        action = actions.changeLJSingleType(amber_prmtop,
                                                  ArgumentList(var_list))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    messages.write('%s\n' % action)
    action.execute()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def addcoarsegrain(root, amber_prmtop, messages):
    """ Adds coarse graining to topology file via Lula's algo """
    # This implementation doesn't exist anywhere yet, so disable it
    showerror('Warning', 'This functionality is not implemented in Amber yet!')
    return

    # We need a file
    fname = file_chooser('Coarse Grain Parameter', 
                         [('Coarse Grain Parameters', '*.cgparm'),
                          ('All Files', '*')]
    )
   
    if not fname: return

    try:
        action = actions.addCoarseGrain(amber_prmtop, ArgumentList(fname))
        messages.write('%s\n' % action)
        action.execute()
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
   
#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def definesolvent(root, amber_prmtop, messages):
    """ Allows you to define what you consider to be solvent molecules """
    # We need a molecule #
    widget_list = [('Entry', 'List of solvent residue names')]
    var_list = [StringVar()]
    description =('Tell ParmEd that these residue names should be considered\n'
                  'to be solvent residues. This should be a comma-delimited '
                  'list')
    cmd_window = _guiwidgets.ActionWindow('defineSolvent', amber_prmtop,
                                          widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    var = var_list[0].get()
    if not var: return
    # addljtype expects any _non_specified variables to be None
    try:
        action = actions.defineSolvent(amber_prmtop, ArgumentList(var))
        messages.write('%s\n' % action)
        action.execute()
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printbonds(root, amber_prmtop, messages):
    """ Prints bonds containing atoms from a given mask """
    widget_list = [('MaskEntry', 'Atoms 1 or Atoms to analyze for bonds'),
                   ('MaskEntry', 'Atoms 2 (optional)')]
    # We need 5 string variables, then get the description. 
    var_list = [StringVar(), StringVar()]
    description='Prints all bonds containing at least 1 atom in the given mask'
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('printBonds', amber_prmtop,
                                          widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    var = ' '.join([x.get() for x in var_list]).strip()
    if not var: return
    try:
        action = actions.printBonds(amber_prmtop, ArgumentList(var))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    messages.write('Printed BONDs for %s\n' % var)
    # Now make the text window
    window = Toplevel(root)
    window.resizable(True, True)
    window.title('BOND list in %s' % amber_prmtop.parm)
    text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                         pady=5, width=65, height=20)
    text.pack(fill=BOTH, expand=1)
    text.write(action)
    window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printangles(root, amber_prmtop, messages):
    """ Prints angles containing atoms from a given mask """
    widget_list = [('MaskEntry', 'Atoms 1 or Atoms to analyze for angles'),
                   ('MaskEntry', 'Atoms 2 (optional)'),
                   ('MaskEntry', 'Atoms 3 (optional)')]
    var_list = [StringVar(), StringVar(), StringVar()]
    description='Prints all angles containing at least 1 atom in the given mask'
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('printAngles', amber_prmtop,
                                        widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    var = ' '.join([x.get() for x in var_list]).strip()
    if not var: return
    try:
        action = actions.printAngles(amber_prmtop, ArgumentList(var))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    messages.write('Printed ANGLEs for %s\n' % var)
    # Now make the text window
    window = Toplevel(root)
    window.resizable(True, True)
    window.title('ANGLE list in %s' % amber_prmtop.parm)
    text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                          pady=5, width=90, height=20)
    text.pack(fill=BOTH, expand=1)
    text.write(action)
    window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printdihedrals(root, amber_prmtop, messages):
    """ Prints dihedrals containing atoms from a given mask """
    widget_list = [('MaskEntry', 'Atoms 1 or Atoms to analyze for dihedrals'),
                   ('MaskEntry', 'Atoms 2 (optional)'),
                   ('MaskEntry', 'Atoms 3 (optional)'),
                   ('MaskEntry', 'Atoms 4 (optional)')]
    var_list = [StringVar(), StringVar(), StringVar(), StringVar()]
    desc = 'Prints all dihedrals containing at least 1 atom in the given mask'
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('printDihedrals', amber_prmtop,
                                          widget_list, var_list, desc)
    cmd_window.wait_window()
    # See if we got any variables back
    var = ' '.join([x.get() for x in var_list]).strip()
    if not var: return
    try:
        action = actions.printDihedrals(amber_prmtop, ArgumentList(var))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    messages.write('Printed DIHEDRALs for %s\n' % var)
    # Now make the text window
    window = Toplevel(root)
    window.resizable(True, True)
    window.title('DIHEDRAL list in %s' % amber_prmtop.parm)
    text = _guiwidgets.ExitingScrollText(window, None, spacing3=5, padx=5,
                                         pady=5, width=140, height=20)
    text.pack(fill=BOTH, expand=1)
    text.write(action)
    window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def setbond(root, amber_prmtop, messages):
    """ Sets (adds or changes) a bond in the topology file """
    # We need 2 masks, a force constant, and an equilibrium distance
    widget_list = [('MaskEntry', 'First atom in bond'),
                   ('MaskEntry', 'Second atom in bond'),
                   ('Entry', 'Force constant (kcal/mol Ang**2)'),
                   ('Entry', 'Equilibrium Distance (Ang)')]
    # We need 4 variables
    var_list = [StringVar(), StringVar(), StringVar(), StringVar()]
    description = ('Sets a bond in the topology file with the given Force '
                   'constant in kcal/mol/Ang**2\nand the given equilibrium '
                   'distance in Angstroms. Both masks must specify only a\n'
                   'single atom. If the bond exists, it will be replaced. If '
                   'it doesn\'t, it will be added.')
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('setBond', amber_prmtop, 
                                          widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    vars_found = True in [bool(v.get()) for v in var_list]
    if not vars_found: return
    # If we did, pass them through
    var_list = [v.get() for v in var_list]
    try:
        action = actions.setBond(amber_prmtop, ArgumentList(var_list))
        messages.write('%s\n' % action)
        action.execute()
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def setangle(root, amber_prmtop, messages):
    """ Sets (adds or changes) an angle in the topology file """
    # We need 3 masks, a force constant, and an equilibrium angle
    widget_list = [('MaskEntry', 'First atom in angle'),
                   ('MaskEntry', 'Second (middle) atom in angle'),
                   ('MaskEntry', 'Third atom in angle'),
                   ('Entry', 'Force constant (kcal/mol rad**2)'),
                   ('Entry', 'Equilibrium Angle (Degrees)')]
    # We need 5 variables
    var_list = [StringVar(), StringVar(), StringVar(), StringVar(), StringVar()]
    description = ('Sets an angle in the topology file with the given Force '
                   'constant in kcal/mol/rad**2\nand the given equilibrium '
                   'angle in Degrees. All three masks must specify only a\n'
                   'single atom. If the angle exists, it will be replaced. If '
                   'it doesn\'t, it will be added.')
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('setAngle', amber_prmtop, 
                                          widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    vars_found = True in [bool(v.get()) for v in var_list]
    if not vars_found: return
    # If we did, pass them through
    var_list = [v.get() for v in var_list]
    try:
        action = actions.setAngle(amber_prmtop, ArgumentList(var_list))
        messages.write('%s\n' % action)
        action.execute()
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def adddihedral(root, amber_prmtop, messages):
    """ Adds a dihedral (improper, multiterm, or normal) to the prmtop """
    # We need 4 masks, phi_k, periodicity, phase, scee/scnb, dihedral type
    widget_list = [
            ('MaskEntry', 'First (end) atom in dihedral'),
            ('MaskEntry', 'Second (middle) atom in dihedral'),
            ('MaskEntry', 'Third (middle) atom in dihedral'),
            ('MaskEntry', 'Fourth (end) atom in dihedral'),
            ('Entry', 'Phi Force constant (kcal/mol)'),
            ('Entry', 'Periodicity'),
            ('Entry', 'Phase (Degrees)'),
            ('Entry', 'EEL scaling factor'),
            ('Entry', 'VDW scaling factor'),
            ('Spinbox', 'Dihedral type', 'normal', 'improper')
    ]
    # We need 10 variables
    var_list = [StringVar() for i in range(10)]
    description = ('Adds a dihedral in the topology file with the given Phi '
                   'Force constant in kcal/mol the\ngiven phase in Degrees '
                   'and the given periodicity. All masks must specify only \n'
                   'a single atom. The default Amber values for SCEE/SCNB are '
                   '1.2 and 2.0, respectively.\nSee the Amber manual for '
                   'details about normal, multiterm, and improper dihedrals')
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('setDihedral', amber_prmtop, 
                                          widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    vars_found = True in [bool(v.get()) for v in var_list]
    if not vars_found: return
    # If we did, pass them through
    var_list = [v.get() for v in var_list]
    # Fill scee/scnb in with default values
    if var_list[7] is None: var_list[7] = '1.2'
    if var_list[8] is None: var_list[8] = '2.0'
    # The last argument is a keyword, so append that, then swap the last 2 args
    var_list.insert(9, 'type')
    try:
        action = actions.addDihedral(amber_prmtop, ArgumentList(var_list))
        messages.write('%s\n' % action)
        action.execute()
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    
#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def addatomicnumber(root, amber_prmtop, messages):
    """ Asks the user if they want to add ATOMIC_NUMBER to the prmtop """
    response = askyesno('addAtomicNumber',
                        'Do you want to add the ATOMIC_NUMBER section to %s?' % 
                        amber_prmtop.parm)
    if response:
        action = actions.addAtomicNumber(amber_prmtop, ArgumentList(''))
        action.execute()
        messages.write('%s\n' % action)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def deletedihedral(root, amber_prmtop, messages):
    """ Deletes a dihedral between 4 given atoms """
    # We need 4 masks
    widget_list = [('MaskEntry', 'First (end) atom in dihedral'),
                   ('MaskEntry', 'Second (middle) atom in dihedral'),
                   ('MaskEntry', 'Third (middle) atom in dihedral'),
                   ('MaskEntry', 'Fourth (end) atom in dihedral')]
    # We need 4 variables
    var_list = [StringVar() for i in range(4)]
    description = ('Deletes dihedrals between the atoms specified in mask1, ' +
                   'mask2, mask3, and mask4.\nIt will try to match dihedrals ' +
                   'only in the order given or reverse order. Dihedrals are\n' +
                   'specified by atom N in mask1, mask2, mask3, and mask4,' +
                   'where N is the selected\natom in that mask.')
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('deleteDihedral', amber_prmtop, 
                                          widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    vars_found = True in [bool(v.get()) for v in var_list]
    if not vars_found: return
    # If we did, pass them through
    var_list = [v.get() for v in var_list]
    try:
        action = actions.deleteDihedral(amber_prmtop,
                                              ArgumentList(var_list))
        messages.write('%s\n' % action)
        action.execute()
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def printljmatrix(root, amber_prmtop, messages):
    """ Prints A/B matrix coefficients """
    # We need 1 mask
    widget_list = [('MaskEntry', 'Atoms To Find LJ Interactions')]
    # We need 1 variable
    var_list = [StringVar()]
    description = ('Prints all A- and B-Coefficient elements between the atom '
                   'types specified in <mask> with every other atom type')
    # Create the window, open it, and wait for it to close
    cmd_window = _guiwidgets.ActionWindow('printLJMatrix', amber_prmtop, 
                                          widget_list, var_list, description)
    cmd_window.wait_window()
    # Bail out if we cancelled or something
    var = var_list[0].get()
    if not var: return
    try:
        action = actions.printLJMatrix(amber_prmtop, ArgumentList(var))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return

    # Create our Info window
    window = Toplevel(root)
    window.resizable(True, True)
    window.title('LJ Type List')
    text = _guiwidgets.ExitingScrollText(window, None, spacing3=2, padx=5,
                                         pady=5, width=45, height=30)
    text.pack(fill=BOTH, expand=1)
    text.write(action)
    messages.write('Printed LJ matrix for types [%s]\n' % var)
    window.wait_window()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def timerge(root, amber_prmtop, messages):
    """ Merges a topology file with 2 molecules into a single prmtop for TI """
    # We need coordinates, so check for that
    if amber_prmtop.parm.coords is None:
        showerror(root, 'tiMerge requires you to load coordinates first!')
        return
    # We need 2 masks, a force constant, and an equilibrium distance
    widget_list = [('MaskEntry', 'Mask Unique to Lambda = 0'),
                   ('MaskEntry', 'Mask Unique to Lambda = 1'),
                   ('MaskEntry', 'Softcore Mask in Lambda = 0'),
                   ('MaskEntry', 'Softcore Mask in Lambda = 1'),
                   ('MaskEntry', 'Softcore Mols NOT to Merge (Lambda=0)'),
                   ('MaskEntry', 'Softcore Mols NOT to Merge (Lambda=1)'),
                   ('Entry', 'Tolerance')]
    # We need 4 variables
    var_list = [StringVar() for i in range(7)]
    var_list[6].set('0.0001')
    description = 'Merges 2 molecules inside a prmtop for use with softcore TI'
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('tiMerge', amber_prmtop, 
                                          widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    vars_found = True in [bool(v.get()) for v in var_list]
    if not vars_found: return
    # If we did, pass them through
    var_list = [v.get() for v in var_list]
    if not var_list[6].strip():
        var_list[6] = '0.0001'
    var_list.insert(6, 'tol')
    actions.tiMerge.output = messages
    try:
        action = actions.tiMerge(amber_prmtop, ArgumentList(var_list))
        messages.write('%s\n' % action)
        action.execute()
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def addpdb(root, amber_prmtop, messages):
    """ Adds PDB information to a topology file """

    # Get the name of the file
    type_list = [('PDB', '*.pdb'), ('All Files', '*')]
    fname = file_chooser('PDB File', type_list)
    if not fname: return

    widget_list = [('Checkbutton', 'Require all residue names be consistent.'),
                   ('Checkbutton', 'Print element names in %FLAG ELEMENT'),
                   ('Checkbutton', 'Print all insertion codes, even if all are '
                                   'blank.')]
   
    var_list = [StringVar(value='no') for i in range(3)]
    description = ('Adds information from the PDB file (e.g., CHAIN IDs) to '
                   'the topology file')
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('addPDB', amber_prmtop, widget_list,
                                          var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    vars_found = True in [bool(v.get()) for v in var_list]
    if not vars_found: return
    newvars = [fname]
    # If we did, pass them through
    if var_list[0] == 'yes': newvars.append('strict')
    if var_list[1] == 'yes': newvars.append('elem')
    if var_list[1] == 'yes': newvars.append('allicodes')
    try:
        action = actions.addPDB(amber_prmtop, ArgumentList(newvars))
        messages.write('%s\n' % action)
        action.execute()
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def deletepdb(root, amber_prmtop, messages):
    """ Asks the user if they want to delete the sections added by addPDB """
    response = askyesno('deletePDB',
                        'Do you want to delete the addPDB flags from %s?' %
                        amber_prmtop.parm)
    if response:
        action = actions.deletePDB(amber_prmtop, ArgumentList(''))
        action.execute()
        messages.write('%s\n' % action)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def lmod(root, amber_prmtop, messages):
    """ Prepares a prmtop file for use with lmod """
    response = askyesno('lmod', 'Do you want to adjust LJ A-coefficients for '
                                'use with lmod?')
    if response:
        action = actions.lmod(amber_prmtop, ArgumentList(''))
        action.execute()
        messages.write('%s\n' % action)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def scale(root, amber_prmtop, messages):
    """ Allows us to scale every value in a particular prmtop section """
    # The spinbox is sent with the Spinbox, label, and then a list of all of the
    # values to give to it
    changeable_properties = ['Spinbox', 'Property to scale']
    for flag in amber_prmtop.parm.flag_list:
        if amber_prmtop.parm.formats[flag].type is float:
            changeable_properties.append(flag)
    widget_list = [changeable_properties, ('Entry', 'Value to scale by')]
    # We need 2 string variables, then get the description
    var_list = [StringVar(), StringVar()]
    description = 'Scales all values in a given prmtop section by a given value'
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('scale', amber_prmtop,
                                          widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    vars_found = var_list[0].get() or var_list[1].get()
    if not vars_found: return
    # If we did, store them and pass it to the class
    var_list = [v.get() for v in var_list]
    try:
        action = actions.scale(amber_prmtop, ArgumentList(var_list))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    action.execute()
    messages.write('%s\n' % action)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def summary(root, amber_prmtop, messages):
    """ Prints a summary of the topology file to the messages """
    summary = actions.summary(amber_prmtop, ArgumentList(''))
    messages.write('%s\n' % summary)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def hmassrepartition(root, amber_prmtop, messages):
    """
    Repartitions mass by making H-atoms heavier and attached heteroatoms
    lighter by the same amount to preserve the same total mass
    """
    # The spinbox is sent with the Spinbox, label, and then a list of all of the
    # values to give to it
    widget_list = [('Entry', 'New hydrogen mass (daltons)'),
                   ('Checkbutton', 'Repartition waters?')]
    # We need 2 string variables, then get the description
    var_list = [StringVar(value='3.024'), StringVar(value='no')]
    description = ('Repartitions atomic masses to keep total mass the same '
                   'while increasing H-atom masses to allow using a larger '
                   'time step.')
    # Create the window, open it, then wait for it to close
    cmd_window = _guiwidgets.ActionWindow('HMassChange', amber_prmtop,
                                          widget_list, var_list, description)
    cmd_window.wait_window()
    # See if we got any variables back
    vars_found = var_list[0].get() or var_list[1].get()
    if not vars_found: return
    # If we did, store them and pass it to the class
    var_list = [v.get() for v in var_list]
    # Convert the check button into the correct keyword
    if var_list[1] == 'yes':
        var_list[1] = 'dowater'
    else:
        var_list[1] = ''
    try:
        action = actions.HMassRepartition(amber_prmtop,
                                                ArgumentList(var_list))
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
        return
    action.execute()
    messages.write('%s\n' % action)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def add12_6_4(root, amber_prmtop, messages):
    """
    Adds 12-6-4 potential energy terms
    """
    # The variables we need for changeljpair
    widget_list = [('MaskEntry', 'Divalent ion mask'),
                   ('FileSelector', 'C4 Parameter File'),
                   ('Entry', 'Water model (instead of C4 Params)'),
                   ('FileSelector', 'Pol. Param File'),
                   ('Entry', 'tunfactor')]
    # Variable list -- we need 2 masks and 2 floats
    var_list = [StringVar(), StringVar(), StringVar(), StringVar(), StringVar()]
    var_list[1].set('Pick C4 File')
    var_list[3].set('Pick Pol. File')
    # description
    description = ('Add the r^-4 Lennard-Jones parameter for the 12-6-4 term\n'
                   'used typically for multi-valent ion parameters')
    cmd_window = _guiwidgets.ActionWindow('add12_6_4', amber_prmtop,
                                          widget_list, var_list, description)
    cmd_window.wait_window()
    # Make sure we didn't cancel (or just press OK with no input), or just leave
    vars_exist = [bool(v.get()) for v in var_list]
    if var_list[1].get() == 'Pick C4 File': vars_exist[1] = False
    if var_list[3].get() == 'Pick Pol. File': vars_exist[3] = False
    if vars_exist[1] and vars_exist[2]:
        showerror('Cannot select both C4 parameter file AND water model')
        return
    vars_exist = True in vars_exist
    if not vars_exist: return
    # Now that we did something, do it
    var_list = [v.get() for v in var_list]
    if var_list[1] == 'Pick C4 File': var_list[1] = ''
    if var_list[3] == 'Pick Pol. File': var_list[3] = ''
    args = [var_list[0]]
    if var_list[1]: args.extend(['c4file', var_list[1]])
    if var_list[2]: args.extend(['watermodel', var_list[2]])
    if var_list[3]: args.extend(['polfile', var_list[3]])
    if var_list[4]: args.extend(['tunfactor', var_list[4]])
    try:
        action = actions.add12_6_4(amber_prmtop, ArgumentList(args))
        messages.write('%s\n' % action)
        action.execute()
    except Exception as err:
        showerror('Unexpected Error!', '%s: %s' % (type(err).__name__, err))
