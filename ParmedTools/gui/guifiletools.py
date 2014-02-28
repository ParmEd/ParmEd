""" Tools for getting/giving file names easily for opening and saving """

import tkFileDialog

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def file_chooser(ftype, extensions = [('Amber Prmtop', '*.prmtop'),
                                      ('Amber Prmtop', '*.parm7'),
                                      ('Amber Prmtop', '*.top'),
                                      ('All Files', '*')]):
   """ Opens a file dialog to choose a topology file """
   return tkFileDialog.askopenfilename(filetypes=extensions,
                       title='Select an Amber %s file to load' % ftype)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

def save_file_chooser(ftype, ext, ext_opts=[('All files', '*')]):
   """ Opens a file dialog to choose a file name to save a new file as """
   return tkFileDialog.asksaveasfilename(filetypes=ext_opts,
                defaultextension=ext, title='Save %s file as' % ftype)

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

