""" A list of all of the tools used by parmed """

__version__ = '14.0'
__author__ = 'Jason Swails'
__all__ = [] # This is populated with the ParmEd Actions below

# For the purposes of the API, all of the actions from ParmedActions will be
# imported here and renamed according to the camelCase convention used in
# ParmedActions.Usages. Importing directly from ParmedActions will cease to be
# documented. This should clarify things.

import ParmedTools.ParmedActions as _PA

for key in _PA.Usages:
    # Skip actions that only make sense for the ParmEd interpreter
    if key in ('help', 'go', 'quit', 'source', 'ls', 'cd'): continue
    name = _PA.Usages[key].split()[0]
    exec('%s = _PA.%s' % (name, key))
    __all__.append(name)

del _PA
