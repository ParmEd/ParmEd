"""
This module is designed to implement features not found in Python 2.4 so they
can be used while still supporting Python 2.4. Importing this module has
limits, however, and will not introduce new syntax that was introduced
post-Python 2.4 (like the 'as' keyword in the "except Exception as e"
construct).  It only supports additions that can be written in pure python,
like "any" and "all".  As such, the performance of some functions implemented
here may be a bit lower than what you would find in the native implementations
in later Python versions.
"""

# This can be used as "from compat24 import *" to implement the full
# compatibility layer. The relevant licenses for the various parts of the
# compatibility layer are shown below.

#  -------------- For OrderedDict --------------
#Copyright (c) 2009 Raymond Hettinger
#
#Permission is hereby granted, free of charge, to any person
#obtaining a copy of this software and associated documentation files
#(the "Software"), to deal in the Software without restriction,
#including without limitation the rights to use, copy, modify, merge,
#publish, distribute, sublicense, and/or sell copies of the Software,
#and to permit persons to whom the Software is furnished to do so,
#subject to the following conditions:
#
#    The above copyright notice and this permission notice shall be
#    included in all copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#    OTHER DEALINGS IN THE SOFTWARE.
#  ---------------------------------------------


__all__ = ['any', 'all', 'property', 'wraps']

import collections as _collections
try:
    from functools import wraps
except ImportError:
    wraps = None # Define this later
try:
    import __builtin__ as __builtins__
except ImportError:
    pass

# Support "any" and "all" functions
if not hasattr(__builtins__, 'any'):
    def any(iterable):
        for it in iterable:
            if it: return True
        return False
    def all(iterable):
        for it in iterable:
            if not it: return False
        return True
else:
    any = __builtins__.any
    all = __builtins__.all

# Support property.setter
if not hasattr(__builtins__.property, 'setter'):
    # Taken from https://gist.github.com/romuald/1104222
    class property(__builtins__.property):
        __metaclass__ = type

        def setter(self, method):
            return property(self.fget, method, self.fdel)

        def deleter(self, method):
            return property(self.fget, self.fset, method)

        @__builtins__.property
        def __doc__(self):
            """ Set doc correctly for subclass """
            return self.fget.__doc__
else:
    property = __builtins__.property

# Support collections.OrderedDict
if 'OrderedDict' not in dir(_collections):
    # This class is taken from Raymond Hettinger's ordereddict package on PyPI
    # Clever use of sentinels. It maintains dict lookup performance with a
    # slight increase in required memory and marginally lower performance for
    # some basic tasks. The OrderedDict implementation in Python 2.7 is also
    # pure-python, so it has similar performance as the one implemented here.
    from UserDict import DictMixin
    class OrderedDict(dict, DictMixin):

        def __init__(self, *args, **kwargs):
            if len(args) > 1:
                raise TypeError('expected at most 1 arguments, got %d' %
                                len(args))
            try:
                self.__end
            except AttributeError:
                self.clear()
            self.update(*args, **kwargs)

        def clear(self):
            self.__end = end = []
            end += [None, end, end]     # sentinel node for doubly linked list
            self.__map = {}             # key --> [key, prev, next]
            dict.clear(self)

        def __setitem__(self, key, value):
            if key not in self:
                end = self.__end
                curr = end[1]
                curr[2] = end[1] = self.__map[key] = [key, curr, end]
            dict.__setitem__(self, key, value)
      
        def __delitem__(self, key):
            dict.__delitem__(self, key)
            key, prev, next = self.__map.pop(key)
            prev[2] = next
            next[1] = prev

        def __iter__(self):
            end = self.__end
            curr = end[2]
            while curr is not end:
                yield curr[0]
                curr = curr[2]

        def __reversed__(self):
            end = self.__end
            curr = end[1]
            while curr is not end:
                yield curr[0]
                curr = curr[1]

        def popitem(self, last=True):
            if not self:
                raise KeyError('dictionary is empty')
            if last:
                key = reversed(self).next()
            else:
                key = iter(self).next()
            value = self.pop(key)
            return key, value

        def __reduce__(self):
            items = [[k, self[k]] for k in self]
            tmp = self.__map, self.__end
            del self.__map, self.__end
            inst_dict = vars(self).copy()
            self.__map, self.__end = tmp
            if inst_dict:
                return (self.__class__, (items,), inst_dict)
            return self.__class__, (items,)

        def keys(self):
            return list(self)

        setdefault = DictMixin.setdefault
        update = DictMixin.update
        pop = DictMixin.pop
        values = DictMixin.values
        items = DictMixin.items
        iterkeys = DictMixin.iterkeys
        itervalues = DictMixin.itervalues
        iteritems = DictMixin.iteritems

        def __repr__(self):
            if not self:
                return '%s()' % self.__class__.__name__
            return '%s(%r)' % (self.__class__.__name__, self.items())

        def copy(self):
            return self.__class__(self)

        @classmethod
        def fromkeys(cls, iterable, value=None):
            d = cls()
            for key in iterable:
                d[key] = value
            return d

        def __eq__(self, other):
            if isinstance(other, OrderedDict):
                if len(self) != len(other):
                    return False
                for p, q, in zip(self.items(), other.items()):
                    if p != q:
                        return False
                    return True
            return dict.__eq__(self, other)

        def __ne__(self, other):
            return not self == other

    _collections.OrderedDict = OrderedDict
    del OrderedDict, DictMixin

# Define the "wraps" function in a Python 2.4-compatible way. Take the Python
# 2.5 code and implement a pure-python _functools.partial taken from
# StackOverflow
if wraps is None:
    # Python module wrapper for _functools C module
    # to allow utilities written in Python to be added
    # to the functools module.
    # Written by Nick Coghlan <ncoghlan at gmail.com>
    #   Copyright (C) 2006 Python Software Foundation.
    # See C source code for _functools credits/copyright

    def partial(func, *args, **kwargs):
        """ Emulate Python 2.6's functools.partial """
        def newfunc(*fargs, **fkwargs):
            nkw = kwargs.copy()
            nkw.update(fkwargs)
            nar = args + fargs
            return func(*nar, **nkw)

        newfunc.func = func
        newfunc.args = args
        newfunc.keywords = kwargs

        return newfunc

    # update_wrapper() and wraps() are tools to help write
    # wrapper functions that can handle naive introspection

    WRAPPER_ASSIGNMENTS = ('__module__', '__name__', '__doc__')
    WRAPPER_UPDATES = ('__dict__',)
    def update_wrapper(wrapper,
                       wrapped,
                       assigned = WRAPPER_ASSIGNMENTS,
                       updated = WRAPPER_UPDATES):
        """Update a wrapper function to look like the wrapped function
    
           wrapper is the function to be updated
           wrapped is the original function
           assigned is a tuple naming the attributes assigned directly
           from the wrapped function to the wrapper function (defaults to
           functools.WRAPPER_ASSIGNMENTS)
           updated is a tuple naming the attributes off the wrapper that
           are updated with the corresponding attribute from the wrapped
           function (defaults to functools.WRAPPER_UPDATES)
        """
        for attr in assigned:
            setattr(wrapper, attr, getattr(wrapped, attr))
        for attr in updated:
            getattr(wrapper, attr).update(getattr(wrapped, attr, {}))
        # Return the wrapper so this can be used as a decorator via partial()
        return wrapper

    def wraps(wrapped,
              assigned=WRAPPER_ASSIGNMENTS,
              updated=WRAPPER_UPDATES):
        """Decorator factory to apply update_wrapper() to a wrapper function
    
           Returns a decorator that invokes update_wrapper() with the decorated
           function as the wrapper argument and the arguments to wraps() as the
           remaining arguments. Default arguments are as for update_wrapper().
           This is a convenience function to simplify applying partial() to
           update_wrapper().
        """
        return partial(update_wrapper, wrapped=wrapped,
                       assigned=assigned, updated=updated)
