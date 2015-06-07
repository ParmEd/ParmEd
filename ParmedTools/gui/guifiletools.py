""" Tools for getting/giving file names easily for opening and saving """

import parmed.utils.six.moves.tkinter_tkfiledialog as tkFileDialog

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def file_chooser(ftype, extensions=None, set_var=None):
    """ Opens a file dialog to choose a topology file """
    if extensions is None:
        extensions = [('Amber Prmtop', '*.prmtop'), ('Amber Prmtop', '*.parm7'),
                      ('Amber Prmtop', '*.top'), ('All Files', '*')]
    fname = tkFileDialog.askopenfilename(filetypes=extensions,
                            title='Select an Amber %s file to load' % ftype)
    if set_var is not None:
        set_var.set(fname)
    return fname

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def save_file_chooser(ftype, ext, ext_opts=None):
    """ Opens a file dialog to choose a file name to save a new file as """
    if ext_opts is None:
        ext_opts = [('All files', '*')]
    return tkFileDialog.asksaveasfilename(filetypes=ext_opts,
                        defaultextension=ext, title='Save %s file as' % ftype)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~
