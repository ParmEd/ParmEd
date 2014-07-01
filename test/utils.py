"""
Useful functions for the test cases
"""
from os.path import join, split, abspath
import sys

def get_fn(filename):
    """
    Gets the full path of the file name for a particular test file

    Parameters
    ----------
    filename : str
        Name of the file to get

    Returns
    -------
    str
        Name of the test file with the full path location
    """
    return join(split(abspath(__file__))[0], 'files', filename)
