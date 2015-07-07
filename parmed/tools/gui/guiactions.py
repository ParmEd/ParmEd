""" 
Dispatcher for GUI Action buttons. It calls all of the functions from the
_guiactions, so that's where all of the class-specific methods should be put
"""

from parmed.utils.six.moves.tkinter_messagebox import showerror
from parmed.tools.gui import _guiactions

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def gui_action_dispatcher(root, amber_prmtop, action_name, messages):
    """ 
    Dispatches all of the GUI actions given the action name. This is the only
    externally accessible *thing* in this module. All of the action-specific
    methods are in _guiactions with the same name as those found in
    parmed.tools.actions 
    """

    if not hasattr(_guiactions, action_name.lower()):
        showerror('Not Implemented.', action_name + 
                  ' is not implemented in xParmEd! Use parmed instead.')
        return None

    # Call the function to establish our action
    getattr(_guiactions, action_name.lower())(root, amber_prmtop, messages)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~
