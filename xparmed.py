#!/usr/bin/env python

"""
This program is the GUI incarnation of parmed.py.  It depends on the Tkinter
package, which is the Tcl/Tk interface with Python, as the GUI toolkit. This is
the most common GUI toolkit in Python (and is included in the stdlib).

All it really does is give users a point-and-click, guided tour through 
everything parmed can do
"""

from chemistry import load_file
from optparse import OptionParser
from os.path import exists, split
from ParmedTools.exceptions import ParmError
from ParmedTools import ParmedActions
from ParmedTools import __version__
from ParmedTools.gui.guitools import ParmedApp
from ParmedTools.gui.guifiletools import file_chooser
from ParmedTools.logos import Logo
from ParmedTools.ParmedActions import Action
from ParmedTools.parmlist import ParmList
import Tkinter as tk
from tkMessageBox import showerror
import sys

debug = False

def excepthook(exception_type, exception_value, tb):
    """ Default exception handler """
    import traceback
    if debug: traceback.print_tb(tb)
    showerror('Fatal Error','%s: %s' % (exception_type.__name__,
                                        exception_value))
    sys.exit(1)

def main():
    """ The main function """
    global excepthook, debug
    # Launch the root window
    root = tk.Tk()
    root.resizable(True, True)

    # Replace the default excepthook with mine
    sys.excepthook = excepthook

    # See if we were provided a topology file on the command-line
    parser = OptionParser(usage = '%prog [<prmtop>]')
    parser.add_option('-d', '--debug', dest='debug', default=False, 
                      action='store_true', help='Show detailed tracebacks ' +
                      'when an error is detected.')
    opt, args = parser.parse_args()

    debug = opt.debug

    # If the user provided a CL argument, that is the prmtop_name. Otherwise,
    # open up a file choosing dialog box to get the input from the user
    if len(args) == 0:
        prmtop_name = file_chooser('Topology')
    elif len(args) == 1:
        prmtop_name = args[0]
    else:
        sys.stderr.write('Unknown command-line options. Ignoring\n')
        prmtop_name = file_chooser('Topology')

    # If we chose no prmtop file, 
    if not prmtop_name: raise ParmError('No prmtop chosen!')

    # Load the amber prmtop and check for errors
    amber_prmtop = ParmList()
    parm = load_file(prmtop_name)
    amber_prmtop.add_parm(parm)

    # Make this overwritable -- the all of the file save boxes will ask the user
    # for verification before saving over an existing file. There's no need for
    # the AmberParm security layer.
    Action.overwrite = True

    fname = split(prmtop_name)[1]
    root.title('xParmED: Editing/viewing [%s] Choose an operation' % fname)

    # Now build the action list on root
    app = ParmedApp(root, amber_prmtop)
    app.pack(fill=tk.BOTH, expand=1)
    root.mainloop()

    print('Thank you for using xParmEd!\n%s' % Logo())

if __name__ == '__main__': main()
