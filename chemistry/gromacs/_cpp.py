"""
A little utiltity for performing C-like preprocessing using standard CPP
directives like #if, #ifdef, and #define.

Written by Jason Swails
"""
from chemistry.exceptions import PreProcessorError, PreProcessorWarning
from chemistry.utils.six import string_types, iteritems, wraps
from os import path
import re
import warnings

ppre = re.compile(r'#\s*(ifdef|ifndef|if|else|elif|endif|define|undef|include)'
                  r'\s*(.+)?')
ppcomments = re.compile(r'(?://.+|/\*(?:.*)\*/)')
includere = re.compile(r'[<"](.+)[>"]')

def _strip_pp_comments(func):
    """
    Decorator to apply to functions that will strip out C-style comments before
    calling the preprocessor function on the resulting arguments to the
    preprocessor directive
    """
    @wraps(func)
    def wrapper(self, args):
        args = ppcomments.sub('', args)
        return func(self, args)
    return wrapper

# To track where in the "if-elif-else" block each conditional is
_IN_IF = 'in if'
_IN_ELIF = 'in elif'
_IN_ELSE = 'in else'

class CPreProcessor(object):
    """
    Steps through a file line-by-line, yielding only the preprocessed contents

    Parameters
    ----------
    fname : str or file-like
        The file name of the file to open or the open file object
    defines : dict{str : value}, optional
        The dict of defines to apply to this preprocessed file if any;
        equivalent to the list of -DXYZ=ABC in the standard preprocessor
    includes : list of str, optional
        List of directory names to search for included files, if any
    notfound_fatal : bool, optional
        If True, include files not found are fatal. If False, they will simply
        be skipped (with a warning emitted). Default True

    Notes
    -----
    If ``fname`` is a file name, the directory containing that file name is the
    first directory searched for include files. If ``fname`` is a file-like
    object, then the current directory is the first searched.
    """

    def __init__(self, fname, defines=None, includes=None, notfound_fatal=True):
        if isinstance(fname, string_types):
            self._fileobj = open(fname, 'r')
            self._ownhandle = True
            curpath = path.abspath(path.split(fname)[0])
        else:
            self._fileobj = fname
            self._ownhandle = False
            curpath = path.abspath(path.curdir)
        if includes is None:
            self._includes = [curpath]
        else:
            self._includes = [curpath] + list(includes)
        if defines is None:
            self.defines = {}
        else:
            # Convert every define to a string
            self.defines = dict()
            for define, value in iteritems(defines):
                self.defines[define] = str(value)
        self._notfound_fatal = notfound_fatal

        # Now to keep track of other basic logic stuff
        self._ifstack = []
        self._elsestack = []
        self._satisfiedstack = []
        self._num_ignoring_if = 0
        self._includefile = None

    def __del__(self):
        self.close()

    def close(self):
        if self._ownhandle:
            self._fileobj.close()

    def readline(self):
        return next(iter(self))

    def readlines(self):
        return [line for line in self]

    def read(self):
        return ''.join(self.readlines())

    def tell(self):
        return self._fileobj.tell()

    def seek(self, value):
        raise NotImplementedError('Cannot seek through a preprocessed file')

    def __iter__(self):
        """ Step through the preprocessed file line-by-line """
        for line in self._fileobj:
            rematch = ppre.match(line)
            if rematch:
                cmd, args = rematch.groups()
                args = args or ''
                self._ppcmdmap[cmd](self, args)
                # If we defined an include file, step through it
                if self._includefile is not None:
                    for line in self._includefile:
                        yield line
                    self._includefile.close()
                    # We have to pass our defines back to our caller
                    self.defines = self._includefile.defines
                    self._includefile = None
                continue

            if self._satisfiedstack and not self._satisfiedstack[-1]:
                # We are inside an unsatisfied conditional
                continue
            for define, value in iteritems(self.defines):
                line = line.replace(define, value)

            yield line

    @_strip_pp_comments
    def _pp_if(self, args):
        if self._satisfiedstack and not self._satisfiedstack[-1]:
            self._num_ignoring_if += 1
            return
        args = args.strip()
        if args.strip() in ('0', '1'):
            satisfied = bool(int(args))
            self._ifstack.append(str(bool(satisfied)))
            self._satisfiedstack.append(bool(satisfied))
            self._elsestack.append(_IN_IF)
            return
        raise NotImplementedError('Only "#if 0|1" is currently supported')

    @_strip_pp_comments
    def _pp_ifdef(self, args):
        if self._satisfiedstack and not self._satisfiedstack[-1]:
            self._num_ignoring_if += 1
            return
        words = args.split()
        if len(words) == 0:
            raise PreProcessorError('Bad #ifdef syntax: "#ifdef %s"' % args)
        elif len(words) > 1:
            warnings.warn('Ignored tokens in #ifdef: %s' % ','.join(words[1:]),
                          PreProcessorWarning)
        self._ifstack.append('%s in self.defines' % words[0])
        self._elsestack.append(_IN_IF)
        self._satisfiedstack.append(words[0] in self.defines)

    @_strip_pp_comments
    def _pp_ifndef(self, args):
        if self._satisfiedstack and not self._satisfiedstack[-1]:
            self._num_ignoring_if += 1
            return
        words = args.split()
        if len(words) == 0:
            raise PreProcessorError('Bad #ifndef syntax: "#ifndef %s"' % args)
        elif len(words) > 1:
            warnings.warn('Ignored tokens in #ifndef: %s' % ','.join(words[1:]),
                          PreProcessorWarning)
        self._ifstack.append('%s not in self.defines' % words[0])
        self._elsestack.append(_IN_IF)
        self._satisfiedstack.append(words[0] not in self.defines)

    @_strip_pp_comments
    def _pp_elif(self, args):
        raise NotImplementedError('#elif conditionals are not yet implemented.')

    @_strip_pp_comments
    def _pp_else(self, args):
        if self._num_ignoring_if > 0:
            return
        if not self._ifstack:
            raise PreProcessorError('#else missing #if(def)')
        words = args.split()
        if len(words) > 0:
            warnings.warn('Ignored tokens in #else: %s' % ','.join(words[1:]),
                          PreProcessorWarning)
        if self._elsestack[-1] == _IN_ELSE:
            raise PreProcessorError('#else following #else')
        self._elsestack[-1] = _IN_ELSE
        self._satisfiedstack[-1] = not self._satisfiedstack[-1]

    @_strip_pp_comments
    def _pp_endif(self, args):
        if self._num_ignoring_if > 0:
            self._num_ignoring_if -= 1
            return
        if not self._ifstack:
            raise PreProcessorError('#endif missing #if(def)')
        self._ifstack.pop()
        self._satisfiedstack.pop()
        self._elsestack.pop()

    @_strip_pp_comments
    def _pp_include(self, args):
        if self._satisfiedstack and not self._satisfiedstack[-1]:
            return
        # Locate the include file
        rematch = includere.match(args)
        if not rematch:
            raise PreProcessorError('Bad #include syntax')
        includefile = rematch.groups()[0]
        for folder in self._includes:
            testfile = path.join(folder, includefile)
            if path.isfile(testfile):
                break
        else:
            if self._notfound_fatal:
                raise PreProcessorError('Could not find %s' % includefile)
            warnings.warn('Could not find %s; skipping' % includefile,
                          PreProcessorWarning)
        self._includefile = CPreProcessor(testfile,
                                          defines=self.defines,
                                          includes=self._includes,
                                          notfound_fatal=self._notfound_fatal)

    @_strip_pp_comments
    def _pp_define(self, args):
        if self._satisfiedstack and not self._satisfiedstack[-1]:
            return
        # Define a new variable
        words = args.split()
        if len(words) == 1:
            self.defines[words[0]] = '1'
        elif len(words) >= 2:
            self.defines[words[0]] = args[len(words[0]):].strip()
        elif len(words) == 0:
            raise PreProcessorError('Nothing defined in #define')

    @_strip_pp_comments
    def _pp_undef(self, args):
        if self._satisfiedstack and not self._satisfiedstack[-1]:
            return
        # Define a new variable
        words = args.split()
        if len(words) == 1:
            try:
                del self.defines[words[0]]
            except KeyError:
                pass
        elif len(words) > 1:
            raise PreProcessorError("Too many tokens found in #undef: %s" %
                                    ','.join(words))
        elif len(words) == 0:
            raise PreProcessorError('Nothing defined in #undef')

    _ppcmdmap = {'if' : _pp_if, 'elif' : _pp_elif, 'ifdef' : _pp_ifdef,
                 'else' : _pp_else, 'define' : _pp_define, 'undef' : _pp_undef,
                 'include' : _pp_include, 'endif' : _pp_endif,
                 'ifndef' : _pp_ifndef}
