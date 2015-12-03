"""
This module contains the class for storing and creating/converting/writing
OpenMM-style ffxml files defining a force field

Author(s): Jason Swails
"""
from __future__ import absolute_import, print_function, division

from copy import copy as _copy
from parmed.formats.registry import FileFormatType
from parmed.parameters import ParameterSet
from parmed.topologyobjects import NoUreyBradley
from parmed.utils.six import add_metaclass, string_types

@add_metaclass(FileFormatType)
class OpenMMParameterSet(ParameterSet):
    """ Class storing parameters from an OpenMM parameter set

    Parameters
    ----------
    filenames : str, list of str, file-like, or list of file-like; optional
        Either the name of a file or a list of filenames from which parameters
        should be parsed.

    Notes
    -----
    Order is important in the list of files provided. The parameters are loaded
    in the order they are provided, and any parameters that are specified in
    multiple places are overwritten (that is, the *last* occurrence is the
    parameter type that is used)

    See Also
    --------
    :class:`parmed.parameters.ParameterSet`
    """

    @staticmethod
    def id_format(filename):
        """
        Identifies the file type as either an Amber-style frcmod or parm.dat
        file.

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is an Amber-style parameter file. False otherwise.
        """
        # Not implemented yet
        return False

    def __init__(self, *filenames):
        super(OpenMMParameterSet, self).__init__()
        if filenames:
            raise NotImplementedError('Cannot yet read OpenMM Parameter sets')

    @classmethod
    def from_parameterset(cls, params, copy=False):
        """
        Instantiates a CharmmParameterSet from another ParameterSet (or
        subclass). The main thing this feature is responsible for is converting
        lower-case atom type names into all upper-case and decorating the name
        to ensure each atom type name is unique.

        Parameters
        ----------
        params : :class:`parmed.parameters.ParameterSet`
            ParameterSet containing the list of parameters to be converted to a
            CHARMM-compatible set
        copy : bool, optional
            If True, the returned parameter set is a deep copy of ``params``. If
            False, the returned parameter set is a shallow copy. Default False.

        Returns
        -------
        new_params : OpenMMParameterSet
            OpenMMParameterSet with the same parameters as that defined in the
            input parameter set
        """
        new_params = cls()
        if copy:
            # Make a copy so we don't modify the original
            params = _copy(params)
        new_params.atom_types = new_params.atom_types_str = params.atom_types
        new_params.atom_types_int = params.atom_types_int
        new_params.atom_types_tuple = params.atom_types_tuple
        new_params.bond_types = params.bond_types
        new_params.angle_types = params.angle_types
        new_params.urey_bradley_types = params.urey_bradley_types
        new_params.dihedral_types = params.dihedral_types
        new_params.improper_types = params.improper_types
        new_params.improper_periodic_types = params.improper_periodic_types
        new_params.rb_torsion_types = params.rb_torsion_types
        new_params.cmap_types = params.cmap_types
        new_params.nbfix_types = params.nbfix_types
        new_params.pair_types = params.pair_types
        new_params.parametersets = params.parametersets
        new_params._combining_rule = params.combining_rule
        new_params.residues = params.residues

        return new_params

    def write(self, dest, provenance=None):
        """ Write the parameter set to an XML file for use with OpenMM

        Parameters
        ----------
        dest : str or file-like
            The name of the file or the file-like object (with a ``write``
            attribute) to which the XML file will be written
        provenance : str, optional
            Provenance information for the force field being converted
        """
        if isinstance(dest, string_types):
            dest = genopen(dest, 'w')
            own_handle = True
        else:
            own_handle = False

        try:
            if provenance is not None:
                self._write_omm_provenance(provenance)
            self._write_omm_residues(dest)
            self._write_omm_bonds(dest)
            self._write_omm_angles(dest)
            self._write_omm_ubs(dest)
            self._write_omm_dihedrals(dest)
            self._write_omm_periodic_impropers(dest)
            self._write_omm_impropers(dest)
            self._write_omm_rb_torsions(dest)
            self._write_omm_cmaps(dest)
            self._write_omm_nonbonded(dest)
        finally:
            if own_handle:
                dest.close()
