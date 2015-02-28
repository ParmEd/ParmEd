"""
Tools to aid in input/output within the chemistry package
"""

__all__ = ['genopen']

try:
    import gzip
except ImportError:
    gzip = None

try:
    import bz2
except ImportError:
    bz2 = None

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
            return bz2.BZ2File(name, mode, buffering)
        else:
            return bz2.BZ2File(name, mode)
    elif name.endswith('.gz'):
        if gzip is None:
            raise ImportError('gzip is unavailable; cannot read %s' % name)
        if buffering is not None:
            return gzip.open(name, mode, buffering)
        else:
            return gzip.open(name, mode)

    if buffering is not None:
        return open(name, mode, buffering)
    return open(name, mode)

class TextToBinaryFile(object):
    """ Allows you to write text to a file open only for bytes in Python 3 """
    def __init__(self, fileobj):
        self._handle = fileobj

    def write(self, stuff):
        try:
            self._handle.write(stuff.encode())
        except (AttributeError, TypeError):
            self._handle.write(stuff)

    def read(self, *args, **kwargs):
        stuff = self._handle.read(*args, **kwargs)
        try:
            return stuff.decode()
        except AttributeError:
            return stuff

    def readline(self):
        line = self._handle.readline()
        try:
            return line.decode()
        except AttributeError:
            return line

    def readlines(self):
        lines = self._handle.readlines()
        for i, line in enumerate(lines):
            try:
                lines[i] = line.decode()
            except AttributeError:
                pass

    def __iter__(self):
        for line in self._handle:
            try:
                yield line.decode()
            except AttributeError:
                yield line

    def close(self):
        self._handle.close()

    def flush(self):
        self._handle.flush()
