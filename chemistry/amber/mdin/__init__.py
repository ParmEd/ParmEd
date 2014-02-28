from chemistry import __version__ as _chemistry_version

__version__ = _chemistry_version
__author__ = "Jason Swails <jason.swails@gmail.com>"
__license__ = 'GPL v3'

del _chemistry_version

from chemistry.amber.mdin.mdin import Mdin

# For backwards-compatibility

mdin = Mdin

__all__ = ['mdin', 'Mdin']
