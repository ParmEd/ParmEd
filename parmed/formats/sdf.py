"""
This module contains parsers for sdf file format 
extension described at [TODO] 
"""
from __future__ import print_function, division, absolute_import

from parmed.formats.registry import FileFormatType
from parmed.utils.io import genopen
from parmed.utils.six import add_metaclass
from parmed import rdkit

@add_metaclass(FileFormatType)
class SDFFile(object):
    """ Class to read SDF file """

    @staticmethod
    def id_format(filename):
        """ Identify the file as a SDF file format or not

        Parameters
        ----------
        filename : str
            Name of the file to test whether or not it is a sdf file

        Returns
        -------
        is_fmt : bool
            True if it is a sdf file, False otherwise
        """
        f = genopen(filename, 'r')
        try:
            for line in f:
                if not line.strip(): continue
                if line.startswith('M  END') or line.startswith('$$$$'):
                    return True
                else:
                    continue
            return False
        finally:
            f.close()

    @staticmethod
    def parse(filename, structure=False):
        return rdkit.from_sdf(filename, structure=structure)
