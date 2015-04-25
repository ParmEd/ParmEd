"""
A package dealing with PyRosetta integration.
"""

from chemistry import __version__ as _chemistry_version

__version__ = _chemistry_version
__author__ = "Carlos Xavier Hernandez <cxh@stanford.edu>"

__all__ = ['RosettaPose']

from chemistry.rosetta.pose import RosettaPose

# Now let's modify structure.Structure and add our write methods from our
# various formats
# from chemistry.structure import Structure as _Structure
# _Structure.dump_pose = RosettaPose.dump
