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
from parmed.periodic_table import Element
from parmed.topologyobjects import NoUreyBradley
from parmed import unit as u
from parmed.utils.io import genopen
from parmed.utils.six import add_metaclass, string_types, iteritems
from parmed.utils.six.moves import range
import warnings
from parmed.exceptions import ParameterWarning

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
        new_params.default_scee = params.default_scee
        new_params.default_scnb = params.default_scnb
        # add only ResidueTemplate instances (no ResidueTemplateContainers)
        for name, residue in iteritems(params.residues):
            if isinstance(residue, ResidueTemplate):
                new_params.residues[name] = residue

        return new_params

    def write(self, dest, provenance=None, write_unused=True, separate_ljforce=False):
        """ Write the parameter set to an XML file for use with OpenMM

        Parameters
        ----------
        dest : str or file-like
            The name of the file or the file-like object (with a ``write``
            attribute) to which the XML file will be written
        provenance : dict, optional
            If present, the XML file will be tagged with the available fields.
            Keys of the dictionary become XML element tags, the values of the
            dictionary must be instances of any of:
            - str / unicode (Py2) or str (Py3) - one XML element with this
            content is written
            - list - one XML element per each item of the list is written, all
            these XML elements use the same tag (key in provenance dict)
            - dict - one of the keys of this dict must be the same as the key of
            of the provenance dict under which this dict is nested. The value
            under this key becomes the content of the XML element. Remaining keys
            and their values are used to construct attributes of the XML element.
            Note that OrderedDict's should be used to ensure appropriate order
            of the XML elements and their attributes.
            Default is no provenance.
            Example (unordered):
            provenance = {'Reference' : ['Nature', 'Cell'],
                          'Source' : {'Source': 'leaprc.ff14SB', sourcePackage :
                          'AmberTools', sourcePackageVersion : '15'},
                          'User' : 'Mark'}
        write_unused : bool
            If False: a) residue templates using unavailable atom types will not
            be written, b) atom types that are not used in any of the residue
            templates remaining and parameters including those atom types will
            not be written. A ParameterWarning is issued if any such residues are
            found in a).
        separate_ljforce : bool
            If True will use a separate LennardJonesForce to create a
            CostumNonbondedForce to compute L-J interactions. It will set sigma
            to 1 and epsilon to 0 in the NonbondedForce so that the
            NonbondedForce  only calculates the electrostatic contribution. It
            should be set to True when converting a CHARMM force field file that
            doesn't have pair-specific  L-J modifications (NBFIX in CHARMM) so
            that the ffxml conversion is compatible with the main charmm36.xml file.
            Note:
            ----
            When pair-specific L-J modifications are present (NBFIX in CHARMM), this
            behavior is always present and this flag is ignored.

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
        if not write_unused:
            skip_residues = self._find_unused_residues()
            skip_types = self._find_unused_types(skip_residues)
            if skip_residues:
                warnings.warn('Some residue templates using unavailable AtomTypes '
                              'were found. They will not be written to the ffxml '
                              'as write_unused is set to False', ParameterWarning)
        else:
            skip_residues = set()
            skip_types = set()
        if self.atom_types:
            try:
                self.typeify_templates()
            except KeyError:
                warnings.warn('Some residue templates are using unavailable '
                              'AtomTypes', ParameterWarning)
        try:
            dest.write('<ForceField>\n')
            self._write_omm_provenance(dest, provenance)
            self._write_omm_atom_types(dest, skip_types)
            self._write_omm_residues(dest, skip_residues)
            self._write_omm_bonds(dest, skip_types)
            self._write_omm_angles(dest, skip_types)
            self._write_omm_urey_bradley(dest, skip_types)
            self._write_omm_dihedrals(dest, skip_types)
            self._write_omm_impropers(dest, skip_types)
#           self._write_omm_rb_torsions(dest, skip_types)
            self._write_omm_cmaps(dest, skip_types)
            self._write_omm_scripts(dest, skip_types)
            self._write_omm_nonbonded(dest, skip_types, separate_ljforce)
            self._write_omm_LennardJonesForce(dest, skip_types, separate_ljforce)
        finally:
            dest.write('</ForceField>\n')
            if own_handle:
                dest.close()

    def _find_unused_residues(self):
        skip_residues = set()
        for name, residue in iteritems(self.residues):
            if any((atom.type not in self.atom_types for atom in residue.atoms)):
                skip_residues.add(name)
        return skip_residues

    def _find_unused_types(self, skip_residues):
        keep_types = set()
        for name, residue in iteritems(self.residues):
            if name not in skip_residues:
                for atom in residue.atoms:
                    keep_types.add(atom.type)
        return {typ for typ in self.atom_types if typ not in keep_types}

    @staticmethod
    def _templhasher(residue):
        if len(residue.atoms) == 1:
            atom = residue.atoms[0]
            return hash((atom.atomic_number, atom.type, atom.charge))
        # TODO implement hash for polyatomic residues
        return id(residue)

    def _write_omm_provenance(self, dest, provenance):
        dest.write(' <Info>\n')
        dest.write('  <DateGenerated>%02d-%02d-%02d</DateGenerated>\n' %
                   datetime.datetime.now().timetuple()[:3])
        provenance = provenance if provenance is not None else {}
        for tag, content in iteritems(provenance):
            if tag == 'DateGenerated': continue
            if not isinstance(content, list):
                content = [content]
            for sub_content in content:
                if isinstance(sub_content, string_types):
                    dest.write('  <%s>%s</%s>\n' % (tag, sub_content, tag))
                elif isinstance(sub_content, dict):
                    if tag not in sub_content:
                        raise KeyError('Content of an attribute-containing element '
                                       'specified incorrectly.')
                    attributes = [key for key in sub_content if key != tag]
                    element_content = sub_content[tag]
                    dest.write('  <%s' % tag)
                    for attribute in attributes:
                        dest.write(' %s="%s"' % (attribute, sub_content[attribute]))
                    dest.write('>%s</%s>\n' % (element_content, tag))
                else:
                    raise TypeError('Incorrect type of the %s element content' % tag)
        dest.write(' </Info>\n')

    def _write_omm_atom_types(self, dest, skip_types):
        if not self.atom_types: return
        dest.write(' <AtomTypes>\n')
        for name, atom_type in iteritems(self.atom_types):
            if name in skip_types: continue
            assert atom_type.atomic_number >= 0, 'Atomic number not set!'
            if atom_type.atomic_number == 0:
                dest.write('  <Type name="%s" class="%s" mass="%s"/>\n'
                           % (name, name, atom_type.mass)
                           )
            else:
                element = Element[atom_type.atomic_number]
                dest.write('  <Type name="%s" class="%s" element="%s" mass="%s"/>\n'
                           % (name, name, element, atom_type.mass)
                          )
        dest.write(' </AtomTypes>\n')

    def _write_omm_residues(self, dest, skip_residues):
        if not self.residues: return
        written_residues = set()
        dest.write(' <Residues>\n')
        for name, residue in iteritems(self.residues):
            if name in skip_residues: continue
            templhash = OpenMMParameterSet._templhasher(residue)
            if templhash in written_residues: continue
            written_residues.add(templhash)
            if residue.override_level == 0:
                dest.write('  <Residue name="%s">\n' % residue.name)
            else:
                dest.write('  <Residue name="%s" override="%d">\n' % (residue.name,
                           residue.override_level))
            for atom in residue.atoms:
                dest.write('   <Atom name="%s" type="%s" charge="%s"/>\n' %
                           (atom.name, atom.type, atom.charge))
            for bond in residue.bonds:
                dest.write('   <Bond atomName1="%s" atomName2="%s"/>\n' %
                           (bond.atom1.name, bond.atom2.name))
            if residue.head is not None:
                dest.write('   <ExternalBond atomName="%s"/>\n' %
                           residue.head.name)
            if residue.tail is not None and residue.tail is not residue.head:
                dest.write('   <ExternalBond atomName="%s"/>\n' %
                           residue.tail.name)
            dest.write('  </Residue>\n')
        dest.write(' </Residues>\n')

    def _write_omm_bonds(self, dest, skip_types):
        if not self.bond_types: return
        dest.write(' <HarmonicBondForce>\n')
        bonds_done = set()
        lconv = u.angstroms.conversion_factor_to(u.nanometers)
        kconv = u.kilocalorie.conversion_factor_to(u.kilojoule) / lconv**2 * 2
        for (a1, a2), bond in iteritems(self.bond_types):
            if any((a in skip_types for a in (a1, a2))): continue
            if (a1, a2) in bonds_done: continue
            bonds_done.add((a1, a2))
            bonds_done.add((a2, a1))
            dest.write('  <Bond type1="%s" type2="%s" length="%s" k="%s"/>\n'
                       % (a1, a2, bond.req*lconv, bond.k*kconv))
        dest.write(' </HarmonicBondForce>\n')

    def _write_omm_angles(self, dest, skip_types):
        if not self.angle_types: return
        dest.write(' <HarmonicAngleForce>\n')
        angles_done = set()
        tconv = u.degree.conversion_factor_to(u.radians)
        kconv = u.kilocalorie.conversion_factor_to(u.kilojoule) * 2
        for (a1, a2, a3), angle in iteritems(self.angle_types):
            if any((a in skip_types for a in (a1, a2, a3))): continue
            if (a1, a2, a3) in angles_done: continue
            angles_done.add((a1, a2, a3))
            angles_done.add((a3, a2, a1))
            dest.write('  <Angle type1="%s" type2="%s" type3="%s" '
                       'angle="%s" k="%s"/>\n' %
                       (a1, a2, a3, angle.theteq*tconv, angle.k*kconv))
        dest.write(' </HarmonicAngleForce>\n')

    def _write_omm_dihedrals(self, dest, skip_types):
        if not self.dihedral_types and not self.improper_periodic_types: return
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
            if any((a in skip_types for a in (a1, a2, a3, a4))): continue
            if (a1, a2, a3, a4) in diheds_done: continue
            diheds_done.add((a1, a2, a3, a4))
            diheds_done.add((a4, a3, a2, a1))
            dest.write('  <Proper type1="%s" type2="%s" type3="%s" '
                       'type4="%s"' % (nowild(a1), a2, a3, nowild(a4)))
            for i, term in enumerate(dihed):
                i += 1
                dest.write(' periodicity%d="%d" phase%d="%s" k%d="%s"' %
                           (i, term.per, i, term.phase*pconv, i,
                            term.phi_k*kconv))
            dest.write('/>\n')
        # Now do the periodic impropers. OpenMM expects the central atom to be
        # listed first. ParameterSet goes out of its way to list it third
        # (consistent with Amber) except in instances where order is random (as
        # in CHARMM parameter files). But CHARMM parameter files don't have
        # periodic impropers, so we don't have to worry about that here.
        for (a2, a3, a1, a4), improp in iteritems(self.improper_periodic_types):
            if any((a in skip_types for a in (a1, a2, a3, a4))): continue
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
                       'type4="%s" periodicity1="%d" phase1="%s" k1="%s"/>\n' %
                       (a1, nowild(a2), nowild(a3), nowild(a4), improp.per,
                        improp.phase*pconv, improp.phi_k*kconv)
            )
        dest.write(' </PeriodicTorsionForce>\n')

    def _write_omm_impropers(self, dest, skip_types):
        if not self.improper_types: return
        dest.write(' <CustomTorsionForce energy="k*(theta-theta0)^2">\n')
        dest.write('  <PerTorsionParameter name="k"/>\n')
        dest.write('  <PerTorsionParameter name="theta0"/>\n')
        kconv = u.kilocalorie.conversion_factor_to(u.kilojoule)
        tconv = u.degree.conversion_factor_to(u.radian)
        def nowild(name):
            return name if name != 'X' else ''
        for (a1, a2, a3, a4), improp in iteritems(self.improper_types):
            if any((a in skip_types for a in (a1, a2, a3, a4))): continue
            dest.write('  <Improper type1="%s" type2="%s" type3="%s" '
                       'type4="%s" k="%s" theta0="%s"/>\n' %
                       (nowild(a1), nowild(a2), nowild(a3), nowild(a4),
                       improp.psi_k*kconv, improp.psi_eq*tconv)
            )
        dest.write(' </CustomTorsionForce>\n')

    def _write_omm_urey_bradley(self, dest, skip_types):
        if not self.urey_bradley_types: return None
        dest.write(' <!-- Urey-Bradley terms -->\n')
        dest.write(' <AmoebaUreyBradleyForce>\n')
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        _ambfrc = u.kilocalorie_per_mole/u.angstrom**2
        _ommfrc = u.kilojoule_per_mole/u.nanometer**2
        frc_conv = _ambfrc.conversion_factor_to(_ommfrc)
        ureys_done = set()
        for (a1, a2, a3), urey in iteritems(self.urey_bradley_types):
            if any((a in skip_types for a in (a1, a2, a3))): continue
            if (a1, a2, a3) in ureys_done: continue
            if urey == NoUreyBradley: continue
            dest.write('  <UreyBradley type1="%s" type2="%s" type3="%s" d="%s" k="%s"/>\n'
                       % (a1, a2, a3, urey.req*length_conv, urey.k*frc_conv))

        dest.write(' </AmoebaUreyBradleyForce>\n')

    def _write_omm_cmaps(self, dest, skip_types):
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
                    dest.write(' %s' % (grid[base+j]*econv))
                dest.write('\n')
            dest.write('  </Map>\n')
        used_torsions = set()
        for (a1, a2, a3, a4, _, _, _, a5), cmap in iteritems(self.cmap_types):
            if any((a in skip_types for a in (a1, a2, a3, a4, a5))): continue
            if (a1, a2, a3, a4, a5) in used_torsions: continue
            used_torsions.add((a1, a2, a3, a4, a5))
            used_torsions.add((a5, a4, a3, a2, a1))
            dest.write('   <Torsion map="%d" type1="%s" type2="%s" '
                       'type3="%s" type4="%s" type5="%s"/>\n' %
                       (maps[id(cmap)], a1, a2, a3, a4, a5)
            )
        dest.write(' </CmapTorsionForce>\n')

    def _write_omm_nonbonded(self, dest, skip_types, separate_ljforce):
        if not self.atom_types: return
        # Compute conversion factors for writing in natrual OpenMM units.
        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        ene_conv = u.kilocalories.conversion_factor_to(u.kilojoules)

        # Get the 1-4 scaling factors from the torsion list
        scee, scnb = set(), set()
        for key in self.dihedral_types:
            dt = self.dihedral_types[key]
            for t in dt:
                if t.scee: scee.add(t.scee)
                if t.scnb: scnb.add(t.scnb)
        if len(scee) > 1:
            raise NotImplementedError('Cannot currently handle mixed 1-4 '
                    'scaling: Elec. Scaling factors %s detected' %
                    (', '.join([str(x) for x in scee])))
        if len(scnb) > 1:
            raise NotImplementedError('Cannot currently handle mixed 1-4 '
                    'scaling: L-J Scaling factors %s detected' %
                    (', '.join([str(x) for x in scnb])))
        if len(scee) > 0:
            coulomb14scale = 1.0 / scee.pop()
        else:
            coulomb14scale = 1.0 / self.default_scee
        if len(scnb) > 0:
            lj14scale = 1.0 / scnb.pop()
        else:
            lj14scale = 1.0 / self.default_scnb

        # Write NonbondedForce records.
        dest.write(' <NonbondedForce coulomb14scale="%s" lj14scale="%s">\n' %
                   (coulomb14scale, lj14scale))
        dest.write('  <UseAttributeFromResidue name="charge"/>\n')
        for name, atom_type in iteritems(self.atom_types):
            if name in skip_types: continue
            if (atom_type.rmin is not None) and (atom_type.epsilon is not None):
                sigma = atom_type.sigma * length_conv  # in md_unit_system
                epsilon = atom_type.epsilon * ene_conv # in md_unit_system
            else:
                # Dummy atom
                sigma = 1.0
                epsilon = 0.0

            if self.nbfix_types or separate_ljforce:
                # turn off L-J. Will use LennardJonesForce to use CostumNonbondedForce to compute L-J interactions
                sigma = 1.0
                epsilon = 0.0
            # Ensure we don't have sigma = 0
            if (sigma == 0.0):
                if (epsilon == 0.0):
                    sigma = 1.0 # reset sigma = 1
                else:
                    raise ValueError("For atom type '%s', sigma = 0 but "
                                     "epsilon != 0." % name)

            dest.write('  <Atom type="%s" sigma="%s" epsilon="%s"/>\n' %
                       (name, sigma, abs(epsilon)))
        dest.write(' </NonbondedForce>\n')

    def _write_omm_LennardJonesForce(self, dest, skip_types, separate_ljforce):
        if not self.nbfix_types and not separate_ljforce: return
        # Convert Conversion factors for writing in natural OpenMM units
        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        ene_conv = u.kilocalories.conversion_factor_to(u.kilojoules)

        scnb = set()
        for key in self.dihedral_types:
            dt = self.dihedral_types[key]
            for t in dt:
                if t.scnb: scnb.add(t.scnb)
        if len(scnb) > 1:
            raise NotImplementedError('Cannot currently handle mixed 1-4 '
                    'scaling: L-J Scaling factors %s detected' %
                    (', '.join([str(x) for x in scnb])))
        if len(scnb) > 0:
            lj14scale = 1.0 / scnb.pop()
        else:
            lj14scale = 1.0 / self.default_scnb

        # write L-J records
        dest.write(' <LennardJonesForce lj14scale="%s">\n' % lj14scale)
        for name, atom_type in iteritems(self.atom_types):
            if name in skip_types: continue
            if (atom_type.rmin is not None) and (atom_type.epsilon is not None):
                sigma = atom_type.sigma * length_conv  # in md_unit_system
                epsilon = atom_type.epsilon * ene_conv # in md_unit_system
            else:
                # Dummy atom
                sigma = 1.0
                epsilon = 0.0

            # Ensure we don't have sigma = 0
            if (sigma == 0.0):
                if (epsilon == 0.0):
                    sigma = 1.0 # reset sigma = 1
                else:
                    raise ValueError("For atom type '%s', sigma = 0 but "
                                     "epsilon != 0." % name)

            dest.write('  <Atom type="%s" sigma="%s" epsilon="%s"/>\n' %
                       (name, sigma, abs(epsilon)))

        # write NBFIX records
        for (atom_types, value) in iteritems(self.nbfix_types):
            emin = value[0] * ene_conv
            rmin = value[1] * length_conv
            # convert to sigma
            sigma = 2 * rmin/(2**(1.0/6))
            dest.write('  <NBFixPair type1="%s" type2="%s" sigma="%s" epsilon="%s"/>\n' %
                       (atom_types[0], atom_types[1], sigma, emin))
        dest.write(' </LennardJonesForce>\n')


    def _write_omm_scripts(self, dest, skip_types):
        # Not currently implemented, so throw an exception if any unsupported
        # options are specified
        if self.combining_rule == 'geometric':
            raise NotImplementedError('Geometric combining rule not currently '
                                      'supported.')
