"""
This package contains code for reading TINKER-style parameter and system-setup
files (like output from "analyze" and the xyz-format file).
"""

__all__ = ['system', 'parameterfile', 'tinkerfiles', 'XyzFile']
__authors__ = 'Jason Swails'

from parmed.tinker.tinkerfiles import XyzFile
