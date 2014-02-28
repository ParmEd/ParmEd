""" 
Dispatcher for GUI Action buttons. It calls all of the functions from the
_guiactions, so that's where all of the class-specific methods should be put
"""

from ParmedTools.gui import _guiactions
from tkMessageBox import showerror

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def gui_action_dispatcher(root, amber_prmtop, action_name, messages):
   """ 
   Dispatches all of the GUI actions given the action name. This is the only
   externally accessible *thing* in this module. All of the action-specific
   methods are in _guiactions with the same name as those found in ParmedActions
   """

   if not hasattr(_guiactions, action_name.lower()):
      showerror('Not Implemented.', action_name + 
                ' is not implemented in xParmEd! Use parmed.py instead.')
      return None

   # Call the function to establish our action
   getattr(_guiactions, action_name.lower())(root, amber_prmtop, messages)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~
