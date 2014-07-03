"""
Useful functions for the test cases
"""
import os
from os.path import join, split, abspath
import sys
import unittest

# Patches for older Ambers.

if not hasattr(unittest.TestCase, 'assertIsInstance'):
    class TestCase(unittest.TestCase):
        
        def assertIsInstance(self, thing, type):
            if not isinstance(thing, type):
                standardMsg = '%s is not an instance of %r' % (obj, type)
                self.fail(self._formatMessage(msg, standardMsg))

    unittest.TestCase = TestCase

def get_fn(filename, written=False):
    """
    Gets the full path of the file name for a particular test file

    Parameters
    ----------
    filename : str
        Name of the file to get
    written : bool=False
        Was this file written by a test? (If so, it is put in a different
        location)

    Returns
    -------
    str
        Name of the test file with the full path location
    """
    if written:
        return join(split(abspath(__file__))[0], 'files', 'writes', filename)
    else:
        return join(split(abspath(__file__))[0], 'files', filename)

def has_scipy():
    try:
        import scipy.io.netcdf as nc
        return True
    except ImportError:
        return False

def has_netcdf4():
    try:
        import netCDF4
        return True
    except ImportError:
        return False

def has_scientific():
    try:
        from Scientific.IO.NetCDF import NetCDFFile
        return True
    except ImportError:
        return False

def has_pynetcdf():
    try:
        import pynetcdf
        return True
    except ImportError:
        return False

def has_numpy():
    try:
        import numpy as np
        return True
    except ImportError:
        return False
