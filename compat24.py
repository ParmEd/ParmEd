# This module is designed to implement features not found in Python 2.4 so they
# can be used while still supporting Python 2.4

# This can be used as "from compat24 import *" to implement the full
# compatibility layer

__all__ = ['any', 'all', 'property']

import __builtin__

# Support "any" and "all" functions
if not hasattr(__builtin__, 'any'):
   def any(iterable):
      for it in iterable:
         if it: return True
      return False
   def all(iterable):
      for it in iterable:
         if not it: return False
      return True
else:
   any = __builtin__.any
   all = __builtin__.all

# Support property.setter
if not hasattr(__builtin__.property, 'setter'):
   # Taken from https://gist.github.com/romuald/1104222
   class property(__builtin__.property):
      __metaclass__ = type

      def setter(self, method):
         return property(self.fget, method, self.fdel)

      def deleter(self, method):
         return property(self.fget, self.fset, method)

      @__builtin__.property
      def __doc__(self):
         """ Set doc correctly for subclass """
         return self.fget.__doc__
else:
   property = __builtin__.property
