"""
Tools to aid in input/output within the parmed package
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
from io import TextIOWrapper, BytesIO
from parmed.utils.six import PY2
from parmed.utils.six.moves.urllib.request import urlopen

def genopen(name, mode='r'):
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

    Returns
    -------
    file : file-like
        A file-like object in the requested mode
    """
    if mode not in ['w', 'r', 'a']:
        raise ValueError('open mode must be "w", "r", or "a"')

    # Handle arbitrary online files. file:// is just an alias for a local file
    is_url = False
    if name.startswith('file:///'):
        name = name[7:]
    elif name.startswith('http://') or name.startswith('https://'):
        is_url = True
        if mode in ['w', 'a']:
            raise ValueError('Cannot write or append a webpage')

    if name.endswith('.bz2'):
        if mode == 'a':
            raise ValueError('Cannot open Bzipped files in append mode')
        if bz2 is None:
            raise ImportError('bz2 unavailable; cannot read %s' % name)
        # BZ2File does not have a way of taking an arbitrary file-like object,
        # so we have to read everything into memory, decompress it, and then
        # pass it back as a BytesIO object wrapped with TextIOWrapper if it is a
        # URL
        if is_url:
            fileobj = BytesIO()
            fileobj.write(bz2.decompress(urlopen(name).read()))
            fileobj.seek(0)
            return TextIOWrapper(fileobj)
        # Not a URL -- handle like a regular file
        if PY2:
            return bz2.BZ2File(name, mode+'b')
        else:
            return TextIOWrapper(bz2.BZ2File(name, mode+'b'))
    elif name.endswith('.gz'):
        if gzip is None:
            raise ImportError('gzip is unavailable; cannot read %s' % name)
        if PY2:
            if is_url:
                # addinfourl in Python 2 does not have a "tell" attribute, so we
                # need to take the same approach for BZ2File above with the
                # BytesIO object... sigh. Yet another reason to migrate to
                # Python 3
                fileobj = BytesIO()
                fileobj.write(urlopen(name).read())
                fileobj.seek(0)
                return gzip.GzipFile(fileobj=fileobj, mode='r')
            else:
                return gzip.open(name, mode+'b')
        else:
            if is_url:
                return TextIOWrapper(
                        gzip.GzipFile(fileobj=urlopen(name), mode='r')
                )
            else:
                return TextIOWrapper(gzip.open(name, mode+'b'))

    if is_url:
        if PY2:
            return urlopen(name)
        else:
            return TextIOWrapper(urlopen(name))
    else:
        return open(name, mode)
