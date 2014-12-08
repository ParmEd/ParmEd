"""
Useful functions for the test cases
"""
import os
from os.path import join, split, abspath
import sys
import unittest
try:
    import numpy
except ImportError:
    numpy = None

# Patches for older Pythons.

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

def get_saved_fn(filename):
    """
    Gets the full path of a file name of a saved test case that is used for
    comparison with a generated file

    Parameters
    ----------
    filename : str
        Name of the file to get
    
    Returns
    -------
    str
        Name of the test file with the full path location
    """
    return join(split(abspath(__file__))[0], 'files', 'saved', filename)

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
    return numpy is not None

def diff_files(file1, file2, ignore_whitespace=True,
               absolute_error=None, relative_error=None):
    """
    Compares 2 files line-by-line

    Parameters
    ----------
    file1 : str or file-like
        Name of the first file to compare or first file object to compare
    file2 : str or file-like
        Name of the second file to compare or second file object to compare
    ignore_whitespace : bool=True
        If true, ignores differences in leading and trailing whitespace
    absolute_error : float=None
        If set, only differences greater than the absolute error will trigger
        failures. Cannot be used with relative_error
    relative_error : float=None
        If set, only relative differences greater than the relative error will
        trigger failures. Cannot be used with absolute_error

    Returns
    -------
    bool
        True if files match. False if they do not or one file does not exist
    
    Notes
    -----
    This routine is not protected against bad types of input. AttributeError may
    be raised if readline() is not implemented on any file-like object passed in
    """
    if absolute_error is not None and relative_error is not None:
        raise ValueError('Cannot specify absolute_error AND relative_error')
    if absolute_error is not None: absolute_error = float(absolute_error)
    if relative_error is not None: relative_error = float(relative_error)
    if isinstance(file1, str):
        try:
            f1 = open(file1, 'r')
        except IOError:
            print('Could not find %s' % file1)
            return False
    else:
        f1 = file1
        file1 = str(file1)
    if isinstance(file2, str):
        try:
            f2 = open(file2, 'r')
        except IOError:
            print('Could not find %s' % file2)
            return False
    else:
        f2 = file2
        file2 = str(file2)
    try:
        l1 = f1.readline()
        l2 = f2.readline()
        i = 1
        same = True
        if ignore_whitespace:
            while l1 and l2:
                if l1.strip() != l2.strip():
                    if l1.startswith('%VERSION') and l2.startswith('%VERSION'):
                        l1 = f1.readline()
                        l2 = f2.readline()
                        continue
                    if not detailed_diff(l1,l2,absolute_error,relative_error):
                        same = False
                        record_diffs(i, file1, file2, l1, l2)
                l1 = f1.readline()
                l2 = f2.readline()
                i += 1
        else:
            while l1 and l2:
                if l1 != l2:
                    if l1.startswith('%VERSION') and l2.startswith('%VERSION'):
                        l1 = f1.readline()
                        l2 = f2.readline()
                        continue
                    if not detailed_diff(l1,l2,absolute_error,relative_error):
                        same = False
                        record_diffs(i, file1, file2, l1, l2)
                l1 = f1.readline()
                l2 = f2.readline()
                i += 1
    finally:
        f1.close()
        f2.close()

    return same

def record_diffs(i, f1, f2, l1, l2):
    if not os.path.isdir(get_fn('diffs')):
        os.makedirs(get_fn('diffs'))
    f = open(os.path.join(get_fn('diffs'), 'TEST_FAILURES.diff'), 'a')
    f.write('# diff %s %s [line %d]\n' % (f1, f2, i))
    f.write('< %s> %s' % (l1, l2))
    f.close()

def detailed_diff(l1, l2, absolute_error=None, relative_error=None):
    """
    Check individual fields to make sure numbers are numerically equal if the
    lines differ. Also ignore fields that appear to be a file name, since those
    will be system-dependent
    """
    fdir = os.path.split(get_fn('writes'))[0]
    w1 = l1.split()
    w2 = l2.split()
    if len(w1) != len(w2): return False
    for wx, wy in zip(w1, w2):
        try:
            wx = float(wx)
            wy = float(wy)
        except ValueError:
            if isinstance(wx, float) or (wx != wy and not
                    (wx.startswith(fdir) or wy.startswith(fdir))):
                return False
        else:
            if wx != wy:
                if absolute_error is not None and abs(wx-wy) > absolute_error:
                    return False
                elif relative_error is not None:
                    if wx == 0 or wy == 0 and abs(wx-wy) > relative_error:
                        return False
                    if abs((wx / wy) - 1) > relative_error:
                        return False
    return True
