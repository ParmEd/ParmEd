"""
List of exceptions used in cpin creation
"""
import sys

class CpinError(Exception):
   " Generic error with Cpin file generation "

class CpinWarning(Warning):
   " Generic warning during Cpin file generation "

class CpinResidueError(CpinError):
   " Error adding a residue to the CPIN file "

class CpinChargeWarning(CpinWarning):
   " Bad charge definitions that are not consistent with protonation states "

class CpinRefEneWarning(CpinWarning):
   " If not all reference energies are properly pKa-adjusted "

class CpinInputWarning(CpinWarning):
   """ If there is a non-fatal problem with the input variables """

class CpinInputError(CpinError):
   " If the user provides bad input "

def replace_excepthook(debug):
   " This function replaces sys.excepthook with one that suppresses tracebacks "
   def excepthook(exception_type, exception_value, tb):
      import traceback
      sys.stderr.write('%s: %s\n' % (exception_type.__name__, exception_value))
      sys.exit(1)
   if not debug:
      sys.excepthook = excepthook
