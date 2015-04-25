"""
A package dealing with PyRosetta integration.
"""

__all__ = ['RosettaPose']

from chemistry.rosetta import RosettaPose

# Now let's modify structure.Structure and add our write methods from our
# various formats
# from chemistry.structure import Structure as _Structure
# _Structure.dump_pose = RosettaPose.dump
