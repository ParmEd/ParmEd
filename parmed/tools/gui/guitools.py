"""
The GUI components of xparmed
"""
from __future__ import division
from parmed.utils.six import iteritems
from parmed.utils.six.moves import range
from parmed.utils.six.moves.tkinter import *
import parmed.utils.six.moves.tkinter_messagebox as tkMessageBox
from parmed.tools.gui.guiactions import gui_action_dispatcher
from parmed.tools.gui._guiwidgets import MessageWindow

BPR = 3 # Buttons Per Row

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class ParmedApp(Frame):
    """ Main Parmed App Window """
    skipped_actions = ('parmout', 'go', 'quit', 'help', 'setoverwrite', 'ls',
                       'interpolate', 'listparms', 'parm', 'source', 'cd')
    def __init__(self, master, amber_prmtop):
        """ Instantiates the window """
        import math
        from parmed.tools.actions import COMMANDMAP
        global BPR
        Frame.__init__(self, master)
        self.parm = amber_prmtop
        actions = {}
        for action, actioncls in iteritems(COMMANDMAP):
            action_name = actioncls.__name__
            # Skip actions that don't get buttons
            if action in self.skipped_actions: continue
            # Save outparm for later
            elif action == 'outparm':
                outparm_action = _Description(action_name, actioncls.__doc__)
                continue
            actions[action] = _Description(action_name, actioncls.__doc__)
        self.message_frame = MessageWindow(self)
        # Now make the action button frame and put a title on it
        action_frame = Frame(self)
        label = Label(action_frame,text='ParmEd Actions (right-click for help)')
        action_keys = actions.keys(); action_keys.sort()
        # Number of Items Per Column
        IPC = int(math.ceil(len(action_keys) / BPR))
        for i, action in enumerate(action_keys):
            button = ActionButton(action_frame, self.parm, self.message_frame,
                                  actions[action])
            button.grid(column=i//IPC, row=i%IPC+1, padx=10, pady=5, 
                        sticky=N+S+E+W)
        action_frame.grid(column=0, row=0, pady=10, padx=5, sticky=N+S+E+W)
        self.message_frame.grid(column=1,row=0,pady=10,padx=5,sticky=N+S+E+W)
        label.grid(column=0, row=0, columnspan=BPR, sticky=N+S+E+W)
        # Add the two final buttons
        self.outparm = OutparmButton(self, self.parm, self.message_frame,
                                     outparm_action)
        self.outparm.grid(column=0, row=1, padx=10, sticky=N+S+E+W)
        self.exit = Button(self,text='Exit xParmEd',command=self.master.destroy)
        self.exit.grid(column=1, row=1, padx=10, sticky=N+S+E+W)
        # Control expandability of root window
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=5)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=0)
        # Control expandability of action frame
        for i in range(int(math.ceil(len(action_keys)/BPR)) + 1):
            action_frame.rowconfigure(i+1, weight=1)
        for i in range(BPR):
            action_frame.columnconfigure(i, weight=1)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class _Description(object):
    " This formats a description so it can be displayed in different places "

    def __init__(self, action_name, desc_string):
        self.desc_string = desc_string
        self.action_name = action_name


    def info_box(self, event=None):
        """ Display an info button """
        tkMessageBox.showinfo(title=self.action_name, 
                        message=str(' '.join(self.desc_string.split())).strip())

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class ActionButton(Button):
    " A button to dispatch the command to a particular function in guiactions "
    def __init__(self, root, amber_prmtop, messages, action_desc):
        self.amber_prmtop = amber_prmtop
        self.action_desc = action_desc
        self.messages = messages
        Button.__init__(self, root, text=action_desc.action_name,
                        command=self.execute)
        self.bind("<Button-3>", self.action_desc.info_box)

    def execute(self):
        """ Sends the action name to the GUI dispatcher """
        gui_action_dispatcher(self.master, self.amber_prmtop, 
                              self.action_desc.action_name, self.messages)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class OutparmButton(ActionButton):
    " An ActionButton specifically for outparm "
    def __init__(self, root, amber_prmtop, messages, action_desc):
        ActionButton.__init__(self, root, amber_prmtop, messages, action_desc)
        self.config(text='Write a New Topology File Now (can be used any '
                         'number of times)')

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~
