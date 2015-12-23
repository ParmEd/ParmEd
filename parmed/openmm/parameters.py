"""
This module contains the class for storing and creating/converting/writing
OpenMM-style ffxml files defining a force field

Author(s): Jason Swails
"""
from __future__ import absolute_import, print_function, division

from copy import copy as _copy
import datetime
from parmed.formats.registry import FileFormatType
from parmed.modeller.residue import ResidueTemplate
from parmed.parameters import ParameterSet
from parmed.periodic_table import Element, Mass
#from parmed.topologyobjects import NoUreyBradley
from parmed import unit as u
from parmed.utils.io import genopen
from parmed.utils.six import add_metaclass, string_types, iteritems
from parmed.utils.six.moves import range

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
        provenance : dict, optional
            If present, the XML file will be tagged with the available fields.
            The keys of this dictionary are turned into the XML tags, and the
            values become the contents of that tag. Default is no provenance

        Notes
        -----
        The generated XML file will have the XML tag ``DateGenerated`` added to
        the provenance information set to the current date. Therefore, you
        should not provide this information in ``provenance`` (it will be
        removed if it is provided).
        """
        if isinstance(dest, string_types):
            dest = genopen(dest, 'w')
            own_handle = True
        else:
            own_handle = False

        try:
            dest.write('<ForceField>\n')
            self._write_omm_provenance(dest, provenance)
            self._write_omm_atom_types(dest)
            self._write_omm_residues(dest)
            self._write_omm_bonds(dest)
            self._write_omm_angles(dest)
            self._write_omm_dihedrals(dest)
            self._write_omm_impropers(dest)
#           self._write_omm_rb_torsions(dest)
            self._write_omm_cmaps(dest)
#           self._write_omm_scripts(dest)
#           self._write_omm_nonbonded(dest)
        finally:
            dest.write('</ForceField>\n')
            if own_handle:
                dest.close()

    def _write_omm_provenance(self, dest, provenance):
        dest.write(' <Info>\n')
        dest.write('  <DateGenerated>%02d-%02d-%02d</DateGenerated>\n' %
                   datetime.datetime.now().timetuple()[:3])
        provenance = provenance if provenance is not None else {}
        for item, key in iteritems(provenance):
            if item == 'DateGenerated': continue
            dest.write('  <%s>%s</%s>\n' % (item, key, item))
        dest.write(' </Info>\n')

    def _write_omm_atom_types(self, dest):
        dest.write(' <AtomTypes>\n')
        for name, atom_type in iteritems(self.atom_types):
            assert atom_type.atomic_number >= 0, 'Atomic number not set!'
            element = Element[atom_type.atomic_number]
            dest.write('  <Type name="%s" element="%s" mass="%f"/>\n'
                % (name, element, atom_type.mass or Mass[element])
            )
        dest.write(' </AtomTypes>\n')

    def _write_omm_residues(self, dest):
        self.typeify_templates()
        dest.write(' <Residues>\n')
        for name, residue in iteritems(self.residues):
            if not isinstance(residue, ResidueTemplate):
                continue
            dest.write('  <Residue name="%s">\n' % residue.name)
            for atom in residue.atoms:
                dest.write('   <Atom name="%s" type="%s" charge="%f"/>\n' %
                           (atom.name, atom.type, atom.charge))
            for bond in residue.bonds:
                dest.write('   <Bond atomName1="%s" atomName2="%s" />\n' %
                           (bond.atom1.name, bond.atom2.name))
            if residue.head is not None:
                dest.write('   <ExternalBond atomName="%s">\n' %
                           residue.head.name)
            if residue.tail is not None:
                dest.write('   <ExternalBond atomName="%s">\n' %
                           residue.tail.name)
        dest.write(' </Residues>\n')

    def _write_omm_bonds(self, dest):
        if not self.bond_types: return
        dest.write(' <HarmonicBondForce>\n')
        bonds_done = set()
        lconv = u.angstroms.conversion_factor_to(u.nanometers)
        kconv = u.kilocalorie.conversion_factor_to(u.kilojoule) / lconv**2 * 2
        for (a1, a2), bond in iteritems(self.bond_types):
            if (a1, a2) in bonds_done: continue
            bonds_done.add((a1, a2))
            bonds_done.add((a2, a1))
            dest.write('  <Bond type1="%s" type2="%s", length="%f", k="%f"/>\n'
                       % (a1, a2, bond.req*lconv, bond.k*kconv))
        dest.write(' </HarmonicBondForce>\n')

    def _write_omm_angles(self, dest):
        if not self.angle_types: return
        dest.write(' <HarmonicAngleForce>\n')
        angles_done = set()
        tconv = u.degree.conversion_factor_to(u.radians)
        kconv = u.kilocalorie.conversion_factor_to(u.kilojoule) * 2
        for (a1, a2, a3), angle in iteritems(self.angle_types):
            if (a1, a2, a3) in angles_done: continue
            angles_done.add((a1, a2, a3))
            angles_done.add((a3, a2, a1))
            dest.write('  <Angle type1="%s" type2="%s" type3="%s" '
                       'angle="%s" k="%s"/>\n' %
                       (a1, a2, a3, angle.theteq*tconv, angle.k*kconv))
        dest.write(' </HarmonicAngleForce>\n')

    def _write_omm_dihedrals(self, dest):
        if not self.dihedral_types: return
        # In ParameterSet, dihedral_types is *always* of type DihedralTypeList.
        # The from_structure method ensures that, even if the containing
        # Structure has separate dihedral entries for each torsion
        dest.write(' <PeriodicTorsionForce>\n')
        diheds_done = set()
        pconv = u.degree.conversion_factor_to(u.radians)
        kconv = u.kilocalorie.conversion_factor_to(u.kilojoule)
        def nowild(name):
            return name if name != 'X' else ''
        for (a1, a2, a3, a4), dihed in iteritems(self.dihedral_types):
            if (a1, a2, a3, a4) in diheds_done: continue
            diheds_done.add((a1, a2, a3, a4))
            diheds_done.add((a4, a3, a2, a1))
            dest.write('  <Proper type1="%s" type2="%s" type3="%s" type4="%s"'
                       % (nowild(a1), a2, a3, nowild(a4)))
            for i, term in enumerate(dihed):
                i += 1
                dest.write(' periodicity%d="%d" phase%d="%f" k%d="%f"' %
                           (i, term.per, i, term.phase*pconv, i,
                            term.phi_k*kconv))
            dest.write('/>\n')
        # Now do the periodic impropers. OpenMM expects the central atom to be
        # listed first. ParameterSet goes out of its way to list it third
        # (consistent with Amber) except in instances where order is random (as
        # in CHARMM parameter files). But CHARMM parameter files don't have
        # periodic impropers, so we don't have to worry about that here.
        for (a2, a3, a1, a4), improp in iteritems(self.improper_periodic_types):
            # Try to make the wild-cards in the middle
            if a4 == 'X':
                if a2 != 'X':
                    a2, a4 = a4, a2
                elif a3 != 'X':
                    a3, a4 = a4, a3
            if a2 != 'X' and a3 == 'X':
                # Single wild-card entries put the wild-card in position 2
                a2, a3 = a3, a2
            dest.write('  <Improper type1="%s" type2="%s" type3="%s" '
                       'type4="%s" periodicity1="%d" phase1="%f" k1="%f"/>\n' %
                       (a1, nowild(a2), nowild(a3), nowild(a4), improp.per,
                        improp.phase*pconv, improp.phi_k*kconv)
            )
        dest.write(' </PeriodicTorsionForce>\n')

    def _write_omm_impropers(self, dest):
        if not self.improper_types: return
        dest.write(' <CustomTorsionForce energy="k*(theta-theta0)^2">\n')
        dest.write('  <PerTorsionParameter name="k"/>\n')
        dest.write('  <PerTorsionParameter name="theta0"/>\n')
        kconv = u.kilocalorie.conversion_factor_to(u.kilojoule)
        tconv = u.degree.conversion_factor_to(u.radian)
        def nowild(name):
            return name if name != 'X' else ''
        for (a1, a2, a3, a4), improp in iteritems(self.improper_types):
            dest.write('  <Improper type1="%s" type2="%s" type3="%s" type4="%s"'
                       ' k="%f" theta0="%f"/>\n' %
                       (nowild(a1), nowild(a2), nowild(a3), nowild(a4),
                       improp.psi_k*kconv, improp.psi_eq*tconv)
            )
        dest.write(' </CustomTorsionForce>\n')

    def _write_omm_cmaps(self, dest):
        if not self.cmap_types: return
        dest.write(' <CmapTorsionForce>\n')
        maps = dict()
        counter = 0
        econv = u.kilocalorie.conversion_factor_to(u.kilojoule)
        for _, cmap in iteritems(self.cmap_types):
            if id(cmap) in maps: continue
            maps[id(cmap)] = counter
            counter += 1
            dest.write('  <Map>\n')
            grid = cmap.grid.switch_range().T
            for i in range(cmap.resolution):
                dest.write('  ')
                base = i * cmap.resolution
                for j in range(cmap.resolution):
                    dest.write(' %f' % (grid[base+j]*econv))
                dest.write('\n')
            dest.write('  </Map>\n')
        used_torsions = set()
        for (a1, a2, a3, a4, _, _, _, a5), cmap in iteritems(self.cmap_types):
            if (a1, a2, a3, a4, a5) in used_torsions: continue
            used_torsions.add((a1, a2, a3, a4, a5))
            used_torsions.add((a5, a4, a3, a2, a1))
            dest.write('   <Torsion map="%d" type1="%s" type=2"%s" type3="%s" '
                       'type4="%s" type5="%s"/>\n' %
                       (maps[id(cmap)], a1, a2, a3, a4, a5)
            )
        dest.write(' </CmapTorsionForce>\n')
