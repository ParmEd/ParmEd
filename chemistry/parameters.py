"""
This module contains classes for parsing and processing CHARMM parameter,
topology, and stream files. It only extracts atom properties from the
topology files and extracts all parameters from the parameter files

Author: Jason M. Swails
"""
from chemistry import (AtomType, BondType, AngleType, DihedralType,
                       ImproperType, CmapType)
from chemistry.utils.six.moves import range
from collections import OrderedDict

class ParameterSet(object):
    """
    Stores a parameter set defining a force field

    Attributes
    ----------
    atom_types_str : dict(str:AtomType)
        Dictionary mapping the names of the atom types to the corresponding
        AtomType instances
    atom_types_int : dict(int:AtomType)
        Dictionary mapping the serial indexes of the atom types to the
        corresponding AtomType instances
    atom_types_double : dict((str,int):AtomType)
        Dictionary mapping the (name,number) tuple of the atom types to the
        corresponding AtomType instances
    bond_types : dict((str,str):AtomType)
        Dictionary mapping the 2-element tuple of the names of the two atom
        types involved in the bond to the BondType instances
    angle_types : dict((str,str,str):AngleType)
        Dictionary mapping the 3-element tuple of the names of the three atom
        types involved in the angle to the AngleType instances
    urey_bradley_types : dict((str,str,str):BondType)
        Dictionary mapping the 3-element tuple of the names of the three atom
        types involved in the angle to the BondType instances of the
        Urey-Bradley terms
    dihedral_types : dict((str,str,str,str):list(DihedralType))
        Dictionary mapping the 4-element tuple of the names of the four atom
        types involved in the dihedral to the DihedralType instances. Since each
        torsion term can be a multiterm expansion, each item corresponding to a
        key in this dict is a list of `DihedralType`s for each term in the
        expansion
    improper_types : dict((str,str,str,str):ImproperType)
        Dictionary mapping the 4-element tuple of the names of the four atom
        types involved in the improper torsion to the ImproperType instances
    cmap_types : dict((str,str,str,str,str):CmapType)
        Dictionary mapping the 5-element tuple of the names of the five atom
        types involved in the correction map to the CmapType instances
    nbfix_types : dict((str,str):(float,float))
        Dictionary mapping the 2-element tuple of the names of the two atom
        types whose LJ terms are modified to the tuple of the (epsilon,rmin)
        terms for that off-diagonal term
    """

    def __init__(self, *args):
        # Instantiate the list types
        self.atom_types = self.atom_types_str = OrderedDict()
        self.atom_types_int = OrderedDict()
        self.atom_types_tuple = OrderedDict()
        self.bond_types = OrderedDict()
        self.angle_types = OrderedDict()
        self.urey_bradley_types = OrderedDict()
        self.dihedral_types = OrderedDict()
        self.improper_types = OrderedDict()
        self.rb_torsion_types = OrderedDict()
        self.cmap_types = OrderedDict()
        self.nbfix_types = OrderedDict()
        self.parametersets = []

    def condense(self, do_dihedrals=True):
        """
        This function goes through each of the parameter type dicts and
        eliminates duplicate types. After calling this function, every unique
        bond, angle, dihedral, improper, or cmap type will pair with EVERY key
        in the type mapping dictionaries that points to the equivalent type

        Parameters
        ----------
        do_dihedrals : bool=True
            Dihedrals can take the longest time to compress since testing their
            equality takes the longest (this is complicated by the existence of
            multi-term torsions). This flag will allow you to *skip* condensing
            the dihedral parameter types (for large parameter sets, this can cut
            the compression time in half)

        Returns
        -------
        self
            The instance that is being condensed

        Notes
        -----
        The return value allows you to condense the types at construction time.

        Example
        -------
        >>> params = ParameterSet().condense()
        >>> params
        <chemistry.parameters.ParameterSet at 0x7f88757de090>
        """
        # First scan through all of the bond types
        self._condense_types(self.bond_types)
        self._condense_types(self.angle_types)
        self._condense_types(self.urey_bradley_types)
        if do_dihedrals:
            self._condense_types(self.dihedral_types)
            self._condense_types(self.rb_torsion_types)
        self._condense_types(self.improper_types)
        self._condense_types(self.cmap_types)
        return self

    @staticmethod
    def _condense_types(typedict):
        """
        Loops through the given dict and condenses all types.

        Parameter
        ---------
        typedict : dict
            Type dictionary to condense
        """
        keylist = typedict.keys()
        for i in range(len(keylist) - 1):
            key1 = keylist[i]
            for j in range(i+1, len(keylist)):
                key2 = keylist[j]
                if typedict[key1] == typedict[key2]:
                    typedict[key2] = typedict[key1]
