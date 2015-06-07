"""
A package dealing with PyRosetta integration.
"""

__author__ = "Carlos Xavier Hernandez <cxh@stanford.edu>"

__all__ = ['RosettaPose']

from parmed.rosetta.pose import RosettaPose

# Now let's modify structure.Structure and add our write methods from our
# various formats
# from parmed.structure import Structure as _Structure
# _Structure.dump_pose = RosettaPose.dump
