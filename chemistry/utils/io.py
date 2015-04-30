"""
Tools to aid in input/output within the chemistry package
"""
from __future__ import print_function, division, absolute_import

__all__ = ['genopen']

try:
    import bz2
except ImportError:
    bz2 = None
try:
    import gzip
except ImportError:
    gzip = None
from io import TextIOWrapper
from chemistry.utils.six import PY2

def genopen(name, mode='r', buffering=None):
    """
    Opens a file, automatically detecting compression schemes by filename
    extension. Note, these files are opened in binary format, so all lines need
    to be decoded after reading or encoded before writing to make sure that the
    code is Py2-Py3 compatible.

    Parameters
    ----------
    name : str
        Name of the file to open
    mode : str, optional
        Whether to open the file to 'r'ead or 'w'rite. Default is 'r'
    buffering : int, optional
        The buffer size to use. You are suggested to leave this at the default

    Returns
    -------
    file : file-like
        A file-like object in the requested mode
    """
    if mode not in ['w', 'r']:
        raise ValueError('open mode must be "w" or "r"')

    mode += 'b'

    if name.endswith('.bz2'):
        if bz2 is None:
            raise ImportError('bz2 unavailable; cannot read %s' % name)
        if buffering is not None:
            if PY2:
                return bz2.BZ2File(name, mode, buffering)
            else:
                return TextIOWrapper(bz2.BZ2File(name, mode, buffering))
        else:
            if PY2:
                return bz2.BZ2File(name, mode)
            else:
                return TextIOWrapper(bz2.BZ2File(name, mode))
    elif name.endswith('.gz'):
        if gzip is None:
            raise ImportError('gzip is unavailable; cannot read %s' % name)
        if buffering is not None:
            if PY2:
                return gzip.open(name, mode, buffering)
            else:
                return TextIOWrapper(gzip.open(name, mode, buffering))
        else:
            if PY2:
                return gzip.open(name, mode, buffering)
            else:
                return TextIOWrapper(gzip.open(name, mode))

    if buffering is not None:
        return open(name, mode, buffering)
    return open(name, mode)
