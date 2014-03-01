"""
General utilities for cpinutil.py
"""
from cpinutils.exceptions import *

def process_arglist(arglist, argtype):
   """
   This processes an argument list with an arbitrary number of arguments that
   may or may not be comma-delimited (with any number of spaces)
   """
   # If the argument list is not set, just return None
   if arglist is None:
      return None
   # Otherwise, process the arguments
   processed_args = []
   for arg in arglist:
      # Delete any spaces, split on commas, and add this to processed arg list
      for arg in arg.replace(' ', '').split(','):
         if not arg: continue
         try:
            processed_args.append(argtype(arg))
         except ValueError:
            raise CpinInputError('Expected type %s for argument. Got %s' % 
                                 (argtype.__name__, arg))
   
   return processed_args
