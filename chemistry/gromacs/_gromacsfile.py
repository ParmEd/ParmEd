"""
Provides a class for reading GROMACS-style files. The key component to these
files is that the ; character is a comment character and everything after ; is
ignored.
"""
from chemistry.exceptions import GromacsFileError

class GromacsFile(object):
    """
    A GROMACS file that recognizes the ";" character as a 'comment' token. It
    can be iterated over and generally treated like a file object, but only
    spits out strings that have been truncated at its first comment character.
    
    There is currently no way to recognize a ; as a _non_ comment character,
    since allowing an escape character does not seem to be common practice and
    would likely introduce negative performance implications.
    """

    def __init__(self, fname, mode='r'):
        if mode not in ('r', 'w'):
            raise ValueError('Cannot open GromacsFile with mode "%s"' % mode)
        if mode == 'r':
            self.status = 'OLD'
        else:
            self.status = 'NEW'
        try:
            self._handle = open(fname, mode)
        except IOError, e:
            raise GromacsFileError(str(e))
        self.closed = False
        self.line_number = 0

    def write(self, *args, **kwargs):
        return self._handle.write(*args, **kwargs)

    def __iter__(self):
        # Iterate over the file
        for line in self._handle:
            try:
                idx = line.index(';')
                end = '\n'
            except ValueError:
                # There is no comment...
                idx = None
                end = ''
            yield '%s%s' % (line[:idx], end)

    def readline(self):
        self.line_number += 1
        line = self._handle.readline()
        try:
            idx = line.index(';')
            end = '\n'
        except ValueError:
            idx = None
            end = ''
        return '%s%s' % (line[:idx], end)

    def readlines(self):
        return [line for line in self]

    def read(self):
        return ''.join(self.readlines())

    def close(self):
        self._handle.close()
        self.closed = True

    def rewind(self):
        """ Return to the beginning of the file """
        self._handle.seek(0)

    def __del__(self):
        try:
            self.closed or self._handle.close()
        except AttributeError:
            # It didn't make it out of the constructor
            pass
