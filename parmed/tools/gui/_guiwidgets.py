"""
This is a module that contains all of the necessary widgets that xParmEd will
use
"""
# Pull in all of Tkinter into our top-level namespace. This is common practice
# for Tkinter
from parmed.utils.six.moves.tkinter import *
from parmed.utils.six.moves.tkinter_messagebox import showerror, showinfo
from parmed.tools.argumentlist import ArgumentList
from parmed.tools.gui.guifiletools import file_chooser

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class MessageWindow(Frame):
    """ 
    This is a window that just contains all of the messages printed out as 
    actions are completed.
    """
    def __init__(self, root):
        Frame.__init__(self, root)
        # Create the scrollable text area that cannot be modified
        self.text = Text(self, spacing3=5, width=60, wrap=WORD)
        scroller = Scrollbar(self, orient=VERTICAL)
        self.text.config(yscrollcommand=scroller.set)
        scroller.config(command=self.text.yview)
        self.text.config(state=DISABLED)
        self.text.grid(column=0, row=0, sticky=N+S+E+W)
        scroller.grid(column=1, row=0, sticky=N+S+E+W)
        # Make sure our window expands the way we want it to
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.rowconfigure(0, weight=1)
        # Add the title
        self.write('Message log from ParmEd Actions\n')

    def add_line(self, text):
        self.write(text+'\n')

    def write(self, text):
        self.text.config(state=NORMAL)
        self.text.insert(END, text)
        self.text.config(state=DISABLED)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class BaseParmedWindow(Toplevel):
    """ 
    The base xParmEd window. It sits on top of a parent window and disables it
    until it is destroyed itself (in which case it re-enables it)
    """
    def __init__(self, amber_prmtop, title='Parmed Action Window'):
        Toplevel.__init__(self)
        self.amber_prmtop = amber_prmtop
        self.title(title)
        # Bind all events to this window
        self.grab_set()
   
#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class RadioButtonWindow(BaseParmedWindow):
    """ 
    This is a window when the user can specify ONE option of a group of options
    with a single set of radio buttons
    """
    def __init__(self, amber_prmtop, title, desc, variable, button_namelist):
        BaseParmedWindow.__init__(self, amber_prmtop, title)
        self.variable = variable
        # Add a label describing what we're doing
        label = Label(self, text=desc, justify=LEFT, pady=5, padx=10)
        label.grid(column=0, row=0)
        frame = Frame(self, pady=5, padx=10)
        for i, name in enumerate(button_namelist):
            radio = Radiobutton(frame, variable=variable, text=name, value=name)
            radio.grid(column=0, row=i, sticky=W)
        frame.grid(column=0, row=1, sticky=W)
        # Give users a way out
        button = Button(self, text='OK', command=self.ok_destroy)
        button.grid(column=0, row=2, sticky=N+S+E+W)
        button = Button(self, text='Cancel', command=self.cancel_destroy)
        button.grid(column=0, row=3, sticky=N+S+E+W)

    def cancel_destroy(self):
        self.variable.set('')
        self.destroy()

    def ok_destroy(self):
        if not self.variable.get():
            showinfo('Missing Variable', "You did not choose an option.")
        self.destroy()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class MaskEntry(Frame):
    """
    This class is used to accept a Mask from a user. It gives them the
    option of inspecting their mask to see what atoms they will be selecting
    """
    def __init__(self, root, amber_prmtop, desc, var, restore_to):
        self.root = root
        self.restore_to = restore_to
        self.var = var
        self.amber_prmtop = amber_prmtop
        Frame.__init__(self, root)
        # Put a label for the entry
        label = Label(self, text=desc)
        label.grid(column=0, row=0, sticky=N+S+E+W)
        # Put an entry field in the frame
        self.entry = Entry(self, width=40, textvariable=var)
        self.entry.grid(column=0, row=1, sticky=N+S+E+W)
        # Put a button underneath to evaluate the mask
        self.button = Button(self, text='Evaluate Current Mask', 
                             command=self.evaluate_me)
        self.button.grid(column=0, row=2, sticky=N+S+E+W)
        # Now make a toplevel window for the evaluation, but keep it withdrawn
        # until it's needed
        self.window = Toplevel(root)
        self.text = ScrollText(self.window, restore_to, spacing3=5, padx=5,
                               pady=5, width=100, height=20)
        self.text.pack(fill=BOTH, expand=1)
        self.window.withdraw()
        self.window.resizable(True, True)
        self.window.transient(self)

    def evaluate_me(self):
        """ Evaluates the mask and dumps the contents into a Text box """
        from parmed.tools import printDetails
        from parmed.exceptions import MaskError
        mask = self.var.get().strip()
        if not mask: return
        # Clear out the text
        self.window.deiconify()
        self.text.clear()
        try: 
            maskstr = printDetails(self.amber_prmtop, ArgumentList(mask))
            self.text.write(str(maskstr))
        except MaskError:
            self.text.write('Bad mask string [%s]' % mask)
        self.window.grab_set()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class ActionWindow(BaseParmedWindow):
    """ A window with a bunch of user-entry fields
    
    Parameters
    ----------
    title : str
        Title that should be placed on the top bar of the window
    amber_prmtop : :class:`ParmList`
        List of topology files being modified
    widget_list : list of (str, str)
        Widget class names, descriptions of widgets
    var_list : list of XyzVar
        List of Var classes (StringVar, IntVar, FloatVar) corresponding to the
        entries in the widget_list
    description : str
        Description of the action being performed
    """
    def __init__(self, title, amber_prmtop, widget_list, var_list, description):
        self.cancelled = False
        self.var_list = var_list
        # widget_list: list of tuples having widget class names and description
        # title: title of the window (action name)
        # var_list: list of variables that go in each widget
        # description: Description to print out on the window
        BaseParmedWindow.__init__(self, amber_prmtop, title)
        # Add description label
        main_desc = Label(self, text=description, justify=CENTER,
                          pady=5, padx=10)
        main_desc.grid(column=0, row=0, sticky=N+S+E+W)
        # Make a frame to put all of our widgets in
        main_frame = Frame(self, padx=10, pady=5)
        # Now add all of our widgets to our 
        for i, wtemp in enumerate(widget_list):
            wname, wdesc = wtemp[0], wtemp[1]
            var = var_list[i]
            local_frame = Frame(main_frame)
            if wname == 'MaskEntry':
                mywidget = MaskEntry(local_frame, self.amber_prmtop,
                                     wdesc, var, self)
                mywidget.grid(row=0, column=0, sticky=N+S+E+W)
            elif wname == 'Entry':
                wlab = Label(local_frame, text=wdesc)
                wlab.grid(row=0, column=0, sticky=N+S+E+W)
                mywidget = Entry(local_frame, width=40, textvariable=var)
                mywidget.grid(row=1, column=0, sticky=N+S+E+W)
            elif wname == 'Spinbox':
                wlab = Label(local_frame, text=wdesc)
                wlab.grid(row=0, column=0, sticky=N+S+E+W)
                mywidget = Spinbox(local_frame, textvariable=var,
                                   values=wtemp[2:])
                mywidget.grid(row=1, column=0, sticky=N+S+E+W)
            elif wname == 'Checkbutton':
                mywidget = Checkbutton(local_frame, text=wdesc, variable=var,
                                    onvalue='yes', offvalue='no')
                mywidget.grid(row=0, column=0, sticky=N+S+E+W)
            elif wname == 'FileSelector':
                def callback():
                    file_chooser('C4 Parameter File',
                            extensions=[('All files', '*')], set_var=var)
                mywidget = Button(local_frame, textvariable=var,
                                  command=callback)
                wlab = Label(local_frame, text=wdesc)
                wlab.grid(row=0, column=0, sticky=N+S+E+W)
                mywidget.grid(row=1, column=0, sticky=N+S+E+W)
            else:
                showerror('Error!', '%s not implemented yet!' % wname)
                self.destroy()
            local_frame.grid(row=i//2, column=i%2, sticky=N+S+E+W)

        # Now pack the main frame
        main_frame.grid(column=0, row=1, sticky=N+S+E+W)
        # Now make the OK and Cancel buttons
        ok_can_frame = Frame(self, padx=10, pady=5)
        ok_can_frame.grid(column=0, row=2)
        button = Button(ok_can_frame, text='OK', command=self.destroy)
        button.grid(row=0, column=0, sticky=N+S+E+W, padx=10)
        button = Button(ok_can_frame, text='Cancel', command=self.cancel)
        button.grid(row=0, column=1, sticky=N+S+E+W, padx=10)

    def cancel(self):
        """ Clears out the variables and bails out """
        self.cancelled = True
        for var in self.var_list: var.set('')
        self.destroy()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class ScrollText(Frame):
    """ A scrollable text (in the vertical direction) """
    ok_button_label = 'OK / Go Back'
    def __init__(self, root, restore_to, **options):
        self.restore_to = restore_to
        Frame.__init__(self, root)
        self.text = Text(self, **options)
        self.text.grid(column=0, row=0, sticky=N+E+S+W)
        self.scroller = Scrollbar(self, orient=VERTICAL,
                                  command=self.text.yview)
        self.scroller.grid(column=1, row=0, sticky=N+S+E+W)
        self.text.config(yscrollcommand=self.scroller.set)
        self.button = Button(self, text=self.ok_button_label, command=self.kill)
        self.button.grid(column=0, row=1, sticky=E+W, columnspan=2)
        root.protocol('WM_DELETE_WINDOW', self.kill)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.rowconfigure(0, weight=1)
        self.text.config(state=DISABLED)

    def write(self, text):
        """ Writes text to the text window """
        self.text.config(state=NORMAL)
        self.text.insert(END, str(text))
        self.text.config(state=DISABLED)

    def clear(self):
        """ Clears all text from the window """
        self.text.config(state=NORMAL)
        self.text.delete(1.0, END)
        self.text.config(state=DISABLED)

    def kill(self):
        self.master.withdraw()
        # Have the old window grab the focus
        if self.restore_to: 
            self.restore_to.grab_set()
        else: 
            self.grab_release()

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class ExitingScrollText(ScrollText):
    """ A ScrollText that will quit after its OK button is clicked """
    ok_button_label = 'OK / Quit'
    def kill(self):
        self.master.destroy()
