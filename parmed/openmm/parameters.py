"""
This module contains the class for storing and creating/converting/writing
OpenMM-style ffxml files defining a force field

Author(s): Jason Swails
"""
from __future__ import absolute_import, print_function, division

from copy import copy as _copy
import math
from functools import wraps
from contextlib import closing
import datetime
from parmed.constants import DEFAULT_ENCODING
from parmed.formats.registry import FileFormatType
from parmed.modeller.residue import ResidueTemplate, PatchTemplate
from parmed.parameters import ParameterSet
from parmed.periodic_table import Element
from parmed.topologyobjects import NoUreyBradley
from parmed import unit as u
from parmed.utils.io import genopen
from parmed.utils.six import add_metaclass, string_types, iteritems
from parmed.utils.six.moves import range
import warnings
from parmed.exceptions import ParameterWarning, IncompatiblePatchError
import itertools

from collections import OrderedDict

try:
    from lxml import etree
except ImportError:
    etree = None

import logging
LOGGER = logging.getLogger(__name__)

def needs_lxml(func):
    """
    Decorator to raise an ImportError if a function requires lxml but it is not
    present
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        if etree is None:
            raise ImportError('required package lxml could not be found')
        return func(*args, **kwargs)
    return wrapper

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
    def _remediate_residue_template(cls, params, residue):
        """
        Modify non-compliant residue templates to conform with OpenMM requirements.

        * To correctly detect waters, OpenMM ffxml water models must not contain
          non-chemical bond constraints. Theses are removed when importing
          foreign parameter sets (e.g., CHARMM) into OpenMM parameter sets,
          and not restored on conversion from OpenMM to other formats

        """
        # Populate atomic numbers in residue template
        # TODO: This can be removed if the parameter readers are guaranteed to populate this correctly
        for atom in residue.atoms:
            if atom.type not in params.atom_types_str:
                warnings.warn('Residue {} contains atom type {} not found in parameter set and will be dropped.'.format(residue.name, atom.type))
                return False
            atom.atomic_number = params.atom_types_str[atom.type].atomic_number

        # Check waters
        if residue.empirical_chemical_formula == 'H2O':
            # Remove any H-H bonds if they are present
            for bond in list(residue.bonds):
                if (bond.atom1.element_name == 'H') and (bond.atom2.element_name == 'H'):
                    # Remove nonphysical H-H bonds
                    LOGGER.debug('Deleting H-H bond from water residue {}'.format(residue.name))
                    residue.delete_bond(bond)
                elif (bond.atom1.atomic_number == 0) or (bond.atom2.atomic_number == 0):
                    LOGGER.debug('Deleting bonds to virtual sites in residue {}'.format(residue.name))
                    residue.delete_bond(bond)
        return True

    @classmethod
    def from_parameterset(cls, params, copy=False):
        """
        Instantiates a CharmmParameterSet from another ParameterSet (or
        subclass). The main thing this feature is responsible for is converting
        lower-case atom type names into all upper-case and decorating the name
        to ensure each atom type name is unique.

        Warning
        -------
        Converting parameter sets to OpenMM can be lossy, and can modify the
        original parameter set unless ``copy=True``:
        * To correctly detect waters, OpenMM ffxml water models must not contain
          non-chemical bond constraints. Theses are removed when importing
          foreign parameter sets (e.g., CHARMM) into OpenMM parameter sets,
          and not restored on conversion from OpenMM to other formats.

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
        # Add only ResidueTemplate instances (no ResidueTemplateContainers)
        # Maintain original residue ordering
        remediated_residues = list()
        for name, residue in iteritems(params.residues):
            if isinstance(residue, ResidueTemplate):
                if OpenMMParameterSet._remediate_residue_template(new_params, residue):
                    remediated_residues.append(residue)
            new_params.residues = OrderedDict()
            for residue in remediated_residues:
                new_params.residues[residue.name] = residue
        for name, patch in iteritems(params.patches):
            if isinstance(patch, PatchTemplate):
                new_params.patches[name] = patch

        return new_params

    @needs_lxml
    def write(self, dest, provenance=None, write_unused=True, separate_ljforce=False,
              improper_dihedrals_ordering='default', charmm_imp=False):
        """ Write the parameter set to an XML file for use with OpenMM

        Parameters
        ----------
        dest : str or file-like
            The name of the file or the file-like object (with a ``write``
            attribute) to which the XML file will be written
        provenance : dict, optional
            If present, the XML file will be tagged with the available fields.
            Keys of the dictionary become XML etree.Element tags, the values of the
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
        improper_dihedrals_ordering : str
            The ordering to use when assigning improper torsions in OpenMM.  Default is 'default',
            other option is 'amber'
        charmm_imp: bool
            If True, will check for existence of IMPR in each residue and patch template,
            and write out the explicit improper definition without wildcards in the ffxml file.

        Notes
        -----
        The generated XML file will have the XML tag ``DateGenerated`` added to
        the provenance information set to the current date. Therefore, you
        should not provide this information in ``provenance`` (it will be
        removed if it is provided).
        """
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

        [valid_residues_for_patch, valid_patches_for_residue] = self._determine_valid_patch_combinations(skip_residues)
        LOGGER.debug('Valid patch combinations:')
        for patch_name in self.patches:
            LOGGER.debug('%8s : %s', patch_name, valid_residues_for_patch[patch_name])

        if charmm_imp:
            self._find_explicit_impropers()

        self._compress_impropers()

        root = etree.Element('ForceField')
        self._write_omm_provenance(root, provenance)
        self._write_omm_atom_types(root, skip_types)
        self._write_omm_residues(root, skip_residues, valid_patches_for_residue)
        self._write_omm_patches(root, valid_residues_for_patch)
        self._write_omm_bonds(root, skip_types)
        self._write_omm_angles(root, skip_types)
        self._write_omm_urey_bradley(root, skip_types)
        self._write_omm_dihedrals(root, skip_types, improper_dihedrals_ordering)
        self._write_omm_impropers(root, skip_types)
        #self._write_omm_rb_torsions(root, skip_types)
        self._write_omm_cmaps(root, skip_types)
        self._write_omm_scripts(root, skip_types)
        self._write_omm_nonbonded(root, skip_types, separate_ljforce)
        self._write_omm_LennardJonesForce(root, skip_types, separate_ljforce)

        tree = etree.ElementTree(root)

        xml = etree.tostring(tree, encoding=DEFAULT_ENCODING, pretty_print=True).decode('utf-8')

        if isinstance(dest, string_types):
            with closing(genopen(dest, 'w')) as f:
                f.write(xml)
        else:
            dest.write(xml)

    def _find_explicit_impropers(self):
        improper_harmonic = {}
        improper_periodic = {}

        def get_types(residue, atomname):
            """Return list of atom type(s) that match the given atom name.
            # TODO: Find a more general way to discover all the types that need to be expanded for +N and -C
            """
            C_types = ['CC', 'CD', 'C'] # atom types associated with '-C'
            N_types = ['NH1', 'NH2', 'NH3', 'N', 'NP'] # atom types associated with '+N'

            a_names = [a.name for a in residue.atoms]
            a_types = [a.type for a in residue.atoms]

            if atomname == '-C':
                return C_types
            elif atomname == '+N':
                return N_types
            elif atomname[0] in ['-', '+']:
                raise ValueError('Unknown atom name %s' % atomname)
            else:
                return [ a_types[a_names.index(atomname)] ]

        for name, residue in iteritems(self.residues):
            for impr in residue._impr:
                types = [ get_types(residue, atomname) for atomname in impr ]
                for (t1, t2, t3, t4) in itertools.product(*types):
                    MATCH = False
                    key = tuple(sorted((t1, t2, t3, t4)))
                    altkeys1 = (t1, t2, t3, t4)
                    altkeys2 = (t4, t3, t2, t1)
                    if key in self.improper_types:
                        improper_harmonic[altkeys1] = self.improper_types[key]
                        MATCH = True
                    elif key in self.improper_periodic_types:
                        improper_periodic[altkeys1] = self.improper_periodic_types[key]
                        MATCH = True
                    elif altkeys1 in self.improper_periodic_types:
                        improper_periodic[altkeys1] = self.improper_periodic_types[altkeys1]
                        MATCH = True
                    elif altkeys2 in self.improper_periodic_types:
                        improper_periodic[altkeys1] = self.improper_periodic_types[altkeys2]
                        MATCH = True
                    else:
                        # Check for wildcards
                        key_placeholder = None
                        for anchor in itertools.combinations([t1, t2, t3, t4], 2):
                            key = tuple(sorted([anchor[0], anchor[1], 'X', 'X']))
                            if key in self.improper_types:
                                if MATCH and key != key_placeholder:
                                    flag = (altkeys1[0], altkeys1[-1])
                                    if flag[0] == key_placeholder[0] and flag[1] == key_placeholder[1]:
                                        # Match was already found.
                                        warnings.warn("{} and {} match improper {}. Using {}".format(key, key_placeholder,
                                                      altkeys1, key_placeholder))
                                        break

                                    if flag[0] == key[0] and flag[1] == key[1]:
                                        improper_harmonic[altkeys1] = self.improper_types[key]
                                        warnings.warn("{} and {} match improper {}. Using {}".format(key, key_placeholder,
                                                        altkeys1, key), ParameterWarning)

                                MATCH = True
                                key_placeholder = key
                                improper_harmonic[altkeys1] = self.improper_types[key]
                            if key not in self.improper_types:
                                for anchor in itertools.combinations([t1, t2, t3, t4], 2):
                                    key = tuple(sorted([anchor[0], anchor[1], 'X', 'X']))
                                    if key in self.improper_periodic_types:
                                        if MATCH and key != key_placeholder:
                                            flag = (altkeys1[0], altkeys1[-1])
                                            if flag[0] == key_placeholder[0] and flag[1] == key_placeholder[1]:
                                                # Match already found.
                                                warnings.warn("{} and {} match improper {}. Using {}".format(key,
                                                              key_placeholder, altkeys1, key_placeholder), ParameterWarning)
                                                break

                                            if flag[0] == key[0] and flag[1] == key[1]:
                                                warnings.warn("More than one improper matches for {}. Using {}".format(
                                                        altkeys1, key), ParameterWarning)
                                                improper_periodic[altkeys1] = self.improper_periodic_types[key]

                                        MATCH = True
                                        key_placeholder = key
                                        improper_periodic[altkeys1] = self.improper_periodic_types[key]
                    if not MATCH:
                        warnings.warn("No improper parameter found for {}".format(altkeys1), ParameterWarning)

        self.improper_periodic_types = improper_periodic
        self.improper_types = improper_harmonic

    def _compress_impropers(self):
        """
        OpenMM's ForceField cannot handle impropers that match the same four atoms
        in more than one order, so Peter Eastman wants us to compress duplicates
        and increment the spring constant accordingly.

        """
        if not self.improper_types: return

        unique_keys = OrderedDict() # unique_keys[key] is the key to retrieve the improper from improper_types
        improper_types = OrderedDict() # replacement for self.improper_types with compressed impropers
        for atoms, improper in iteritems(self.improper_types):
            # Compute a unique key
            unique_key = tuple(sorted(atoms))
            if unique_key in unique_keys:
                # Accumulate spring constant, discarding this contribution
                # TODO: Do we need to check if `psi_eq` is the same?
                atoms2 = unique_keys[unique_key]
                improper_types[atoms2].psi_k += improper.psi_k
            else:
                # Store this improper
                unique_keys[unique_key] = atoms
                improper_types[atoms] = improper

        self.improper_types = improper_types

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

    @needs_lxml
    def _write_omm_provenance(self, root, provenance):
        info = etree.SubElement(root, 'Info')

        date_generated = etree.SubElement(info, "DateGenerated")
        date_generated.text = '%02d-%02d-%02d' % datetime.datetime.now().timetuple()[:3]

        provenance = provenance or dict()
        for tag, content in iteritems(provenance):
            if tag == 'DateGenerated': continue
            if not isinstance(content, list):
                content = [content]
            for sub_content in content:
                if isinstance(sub_content, string_types):
                    item = etree.Element(tag)
                    item.text = sub_content
                    info.append(item)
                elif isinstance(sub_content, dict):
                    if tag not in sub_content:
                        raise KeyError('Content of an attribute-containing element '
                                       'specified incorrectly.')
                    attributes = [key for key in sub_content if key != tag]
                    element_content = sub_content[tag]
                    attributes = { k : str(v) for (k,v) in sub_content.items() }
                    item = etree.SubElement(info, tag, **attributes)
                    item.text = str(element_content)
                else:
                    raise TypeError('Incorrect type of the %s element content' % tag)

    @needs_lxml
    def _write_omm_atom_types(self, xml_root, skip_types):
        if not self.atom_types: return
        xml_section = etree.SubElement(xml_root, "AtomTypes")
        for name, atom_type in iteritems(self.atom_types):
            if name in skip_types: continue
            assert atom_type.atomic_number >= 0, 'Atomic number not set!'
            properties = { 'name' : name, 'class' : name, 'mass' : str(atom_type.mass) }
            if atom_type.atomic_number == 0:
                etree.SubElement(xml_section, 'Type', **properties)
            else:
                element = Element[atom_type.atomic_number]
                etree.SubElement(xml_section, 'Type', element=str(element), **properties)

    @needs_lxml
    def _write_omm_residues(self, xml_root, skip_residues, valid_patches_for_residue=None):
        if not self.residues: return
        if valid_patches_for_residue is None:
            valid_patches_for_residue = dict()
        written_residues = set()
        xml_section = etree.SubElement(xml_root, 'Residues')
        for name, residue in iteritems(self.residues):
            if name in skip_residues: continue
            templhash = OpenMMParameterSet._templhasher(residue)
            if templhash in written_residues: continue
            written_residues.add(templhash)
            # Write residue
            if residue.override_level == 0:
                xml_residue = etree.SubElement(xml_section, 'Residue', name=residue.name)
            else:
                xml_residue = etree.SubElement(xml_section, 'Residue', name=residue.name, override=str(residue.override_level))
            # Write residue contents
            for atom in residue.atoms:
                etree.SubElement(xml_residue, 'Atom', name=atom.name, type=atom.type, charge=str(atom.charge))
            for bond in residue.bonds:
                etree.SubElement(xml_residue, 'Bond', atomName1=bond.atom1.name, atomName2=bond.atom2.name)
            for (index, lonepair) in enumerate(residue.lonepairs):
                (lptype, a1, a2, a3, a4, r, theta, phi) = lonepair
                if lptype == 'relative':
                    xweights = [-1.0, 0.0, 1.0]
                elif lptype == 'bisector':
                    xweights = [-1.0, 0.5, 0.5]
                else:
                    raise ValueError('Unknown lonepair type: '+lptype)
                r /= 10.0 # convert to nanometers
                theta *= math.pi / 180.0 # convert to radians
                phi = (180 - phi) * math.pi / 180.0 # convert to radians
                p = [r*math.cos(theta), r*math.sin(theta)*math.cos(phi), r*math.sin(theta)*math.sin(phi)]
                p = [x if abs(x) > 1e-10 else 0 for x in p] # Avoid tiny numbers caused by roundoff error
                etree.SubElement(xml_residue, 'VirtualSite', type="localCoords", index=str(index),
                    siteName=a1, atomName1=a2, atomName2=a3, atomName3=a4,
                    wo1="1", wo2="0", wo3="0",
                    wx1=str(xweights[0]), wx2=str(xweights[1]), wx3=str(xweights[2]),
                    wy1="0", wy2="-1", wy3="1",
                    p1=str(p[0]), p2=str(p[1]), p3=str(p[2]))
            for atom in residue.connections:
                etree.SubElement(xml_residue, 'ExternalBond', atomName=atom.name)
            if residue.head is not None:
                etree.SubElement(xml_residue, 'ExternalBond', atomName=residue.head.name)
            if residue.tail is not None and residue.tail is not residue.head:
                etree.SubElement(xml_residue, 'ExternalBond', atomName=residue.tail.name)
            if residue.name in valid_patches_for_residue:
                for patch_name in valid_patches_for_residue[residue.name]:
                    etree.SubElement(xml_residue, 'AllowPatch', name=patch_name)

    def _determine_valid_patch_combinations(self, skip_residues):
        """
        Determine valid (permissible) combinations of patches with residues that
        lead to integral net charges.

        Parameters
        ----------
        skip_residues : set of ResidueTemplate
            List of residues to skip

        Returns
        -------
        valid_residues_for_patch : dict
            valid_residues_for_patch[patch] is a list of residues compatible with that patch
        valid_patches_for_residue : dict
            valid_patches_for_residue[residue] is a list of patches compatible with that residue

        """
        # Attempt to patch every residue, recording only valid combinations.
        valid_residues_for_patch = OrderedDict()
        for patch in self.patches.values():
            valid_residues_for_patch[patch.name] = list()
        valid_patches_for_residue = OrderedDict()
        for residue in self.residues.values():
            valid_patches_for_residue[residue.name] = list()

        for patch in self.patches.values():
            for residue in self.residues.values():
                if residue in skip_residues: continue
                # Attempt to patch the residue.
                try:
                    residue.apply_patch(patch)
                except IncompatiblePatchError as e:
                    # Patching failed; continue to next patch
                    LOGGER.debug('%8s x %8s : %s', patch.name, residue.name, e)
                    continue

                valid_residues_for_patch[patch.name].append(residue.name)
                valid_patches_for_residue[residue.name].append(patch.name)

        return [valid_residues_for_patch, valid_patches_for_residue]

    @needs_lxml
    def _write_omm_patches(self, xml_root, valid_residues_for_patch, write_apply_to_residue=False):
        """
        Write patch definitions for OpenMM ForceField

        Parameters
        ----------
        xml_root : lxml.etree.Element
            The XML Element write the <Patches> section to.
        valid_residues_for_patch : dict of str : str
            valid_residues_for_patch[patch_name] lists the residue names valid for this patch
        write_apply_to_residue : bool, optional, default=False
            If True, will write <ApplyToResidue> tags.

        """
        if not self.patches: return
        written_patches = set()
        xml_patches = etree.SubElement(xml_root, 'Patches')
        for name, patch in iteritems(self.patches):
            # Require that at least one valid patch combination exists for this patch
            if (name not in valid_residues_for_patch) or (len(valid_residues_for_patch[name])==0):
                continue

            templhash = OpenMMParameterSet._templhasher(patch)
            if templhash in written_patches: continue
            written_patches.add(templhash)
            if patch.override_level == 0:
                patch_xml = etree.SubElement(xml_patches, 'Patch', name=patch.name)
            else:
                patch_xml = etree.SubElement(xml_patches, 'Patch', name=patch.name, override=str(patch.override_level))

            # Construct an example patched residue
            # TODO: We should also ensure that *all* compatible residues have the
            # same added/changed atoms and deleted bonds, just to be safe.
            residue_name = valid_residues_for_patch[name][0]
            try:
                residue = self.residues[residue_name]
            except KeyError as e:
                msg =  'Compatible residue not found in self.residues\n'
                msg += '   patch name: %s\n' % name
                msg += '   valid patch combinations: %s\n' % valid_residues_for_patch[name]
                msg += '   residue name: %s\n' % residue_name
                msg += str(e)
                raise(msg)
            patched_residue = residue.apply_patch(patch)
            #patch_xml.append( etree.Comment('Using residue %s as reference template' % residue.name) )

            for atom in patch.atoms:
                if atom.name not in residue:
                    etree.SubElement(patch_xml, 'AddAtom', name=atom.name, type=atom.type, charge=str(atom.charge))
                else:
                    etree.SubElement(patch_xml, 'ChangeAtom', name=atom.name, type=atom.type, charge=str(atom.charge))

            for atom_name in patch.delete_atoms:
                etree.SubElement(patch_xml, 'RemoveAtom', name=atom_name)

            for bond in patch.bonds:
                etree.SubElement(patch_xml, 'RemoveBond', atomName1=bond.atom1.name, atomName2=bond.atom2.name)

            for bond in patched_residue.bonds:
                if (bond.atom1.name not in residue) or (bond.atom2.name not in residue):
                    etree.SubElement(patch_xml, 'AddBond', atomName1=bond.atom1.name, atomName2=bond.atom2.name)
            for bond in residue.bonds:
                if (bond.atom1.name not in patched_residue) or (bond.atom2.name not in patched_residue):
                    etree.SubElement(patch_xml, 'RemoveBond', atomName1=bond.atom1.name, atomName2=bond.atom2.name)

            if (residue.head is not None) and (patched_residue.head is None):
                etree.SubElement(patch_xml, 'RemoveExternalBond', atomName=residue.head.name)
            if (residue.tail is not None) and (patched_residue.tail is None):
                etree.SubElement(patch_xml, 'RemoveExternalBond', atomName=residue.tail.name)

            if (residue.head is None) and (patched_residue.head is not None):
                etree.SubElement(patch_xml, 'AddExternalBond', atomName=patched_residue.head.name)
            if (residue.tail is None) and (patched_residue.tail is not None):
                etree.SubElement(patch_xml, 'AddExternalBond', atomName=patched_residue.tail.name)

            if write_apply_to_residue:
                for residue_name in valid_residues_for_patch[patch.name]:
                    etree.SubElement(patch_xml, 'ApplyToResidue', name=residue_name)

    @needs_lxml
    def _write_omm_bonds(self, xml_root, skip_types):
        if not self.bond_types: return
        xml_force = etree.SubElement(xml_root, 'HarmonicBondForce')
        bonds_done = set()
        lconv = u.angstroms.conversion_factor_to(u.nanometers)
        kconv = u.kilocalorie.conversion_factor_to(u.kilojoule) / lconv**2 * 2
        for (a1, a2), bond in iteritems(self.bond_types):
            if any((a in skip_types for a in (a1, a2))): continue
            if (a1, a2) in bonds_done: continue
            bonds_done.add((a1, a2))
            bonds_done.add((a2, a1))
            etree.SubElement(xml_force, 'Bond', type1=a1, type2=a2, length=str(bond.req*lconv), k=str(bond.k*kconv))

    @needs_lxml
    def _write_omm_angles(self, xml_root, skip_types):
        if not self.angle_types: return
        xml_force = etree.SubElement(xml_root, 'HarmonicAngleForce')
        angles_done = set()
        tconv = u.degree.conversion_factor_to(u.radians)
        kconv = u.kilocalorie.conversion_factor_to(u.kilojoule) * 2
        for (a1, a2, a3), angle in iteritems(self.angle_types):
            if any((a in skip_types for a in (a1, a2, a3))): continue
            if (a1, a2, a3) in angles_done: continue
            angles_done.add((a1, a2, a3))
            angles_done.add((a3, a2, a1))
            etree.SubElement(xml_force, 'Angle', type1=a1, type2=a2, type3=a3, angle=str(angle.theteq*tconv), k=str(angle.k*kconv))

    @needs_lxml
    def _write_omm_dihedrals(self, xml_root, skip_types, improper_dihedrals_ordering):
        if not self.dihedral_types and not self.improper_periodic_types: return
        # In ParameterSet, dihedral_types is *always* of type DihedralTypeList.
        # The from_structure method ensures that, even if the containing
        # Structure has separate dihedral entries for each torsion
        if improper_dihedrals_ordering == 'default':
            xml_force = etree.SubElement(xml_root, 'PeriodicTorsionForce')
        else:
            xml_force = etree.SubElement(xml_root, 'PeriodicTorsionForce', ordering=improper_dihedrals_ordering)
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
            terms = dict()
            for i, term in enumerate(dihed):
                i += 1
                terms['periodicity%d' % i] = str(term.per)
                terms['phase%d' % i] = str(term.phase*pconv)
                terms['k%d' % i] = str(term.phi_k*kconv)
            etree.SubElement(xml_force, 'Proper', type1=nowild(a1), type2=a2, type3=a3, type4=nowild(a4), **terms)
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
            etree.SubElement(xml_force, 'Improper', type1=a1, type2=nowild(a2), type3=nowild(a3), type4=nowild(a4),
                       periodicity1=str(improp.per), phase1=str(improp.phase*pconv), k1=str(improp.phi_k*kconv))

    @needs_lxml
    def _write_omm_impropers(self, xml_root, skip_types):
        if not self.improper_types: return
        xml_force = etree.SubElement(xml_root, 'CustomTorsionForce', energy="k*(theta-theta0)^2")
        etree.SubElement(xml_force, 'PerTorsionParameter', name="k")
        etree.SubElement(xml_force, 'PerTorsionParameter', name="theta0")
        kconv = u.kilocalorie.conversion_factor_to(u.kilojoule)
        tconv = u.degree.conversion_factor_to(u.radian)
        def nowild(name):
            return name if name != 'X' else ''
        for (a1, a2, a3, a4), improp in iteritems(self.improper_types):
            if any((a in skip_types for a in (a1, a2, a3, a4))): continue
            etree.SubElement(xml_force, 'Improper', type1=nowild(a1), type2=nowild(a2), type3=nowild(a3), type4=nowild(a4),
                       k=str(improp.psi_k*kconv), theta0=str(improp.psi_eq*tconv))

    @needs_lxml
    def _write_omm_urey_bradley(self, xml_root, skip_types):
        if not self.urey_bradley_types: return None
        xml_root.append( etree.Comment("Urey-Bradley terms") )
        xml_force = etree.SubElement(xml_root, 'AmoebaUreyBradleyForce')
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        _ambfrc = u.kilocalorie_per_mole/u.angstrom**2
        _ommfrc = u.kilojoule_per_mole/u.nanometer**2
        frc_conv = _ambfrc.conversion_factor_to(_ommfrc)
        ureys_done = set()
        for (a1, a2, a3), urey in iteritems(self.urey_bradley_types):
            if any((a in skip_types for a in (a1, a2, a3))): continue
            if (a1, a2, a3) in ureys_done: continue
            if urey == NoUreyBradley: continue
            etree.SubElement(xml_force, 'UreyBradley', type1=a1, type2=a2, type3=a3, d=str(urey.req*length_conv), k=str(urey.k*frc_conv))

    @needs_lxml
    def _write_omm_cmaps(self, xml_root, skip_types):
        if not self.cmap_types: return
        xml_force = etree.SubElement(xml_root, 'CMAPTorsionForce')
        maps = dict()
        counter = 0
        econv = u.kilocalorie.conversion_factor_to(u.kilojoule)
        for _, cmap in iteritems(self.cmap_types):
            if id(cmap) in maps: continue
            maps[id(cmap)] = counter
            counter += 1
            xml_map = etree.SubElement(xml_force, 'Map')
            grid = cmap.grid.switch_range().T
            map_string = ''
            for i in range(cmap.resolution):
                base = i * cmap.resolution
                for j in range(cmap.resolution):
                    map_string += ' %s' % (grid[base+j]*econv)
                map_string += '\n'
            xml_map.text = map_string
        used_torsions = set()
        for (a1, a2, a3, a4, _, _, _, a5), cmap in iteritems(self.cmap_types):
            if any((a in skip_types for a in (a1, a2, a3, a4, a5))): continue
            if (a1, a2, a3, a4, a5) in used_torsions: continue
            used_torsions.add((a1, a2, a3, a4, a5))
            used_torsions.add((a5, a4, a3, a2, a1))
            etree.SubElement(xml_force, 'Torsion', map=str(maps[id(cmap)]),
                       type1=a1, type2=a2, type3=a3, type4=a4, type5=a5)

    @needs_lxml
    def _write_omm_nonbonded(self, xml_root, skip_types, separate_ljforce):
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
        xml_force = etree.SubElement(xml_root, 'NonbondedForce', coulomb14scale=str(coulomb14scale), lj14scale=str(lj14scale))
        etree.SubElement(xml_force, 'UseAttributeFromResidue', name="charge")
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
            else:
                # NonbondedForce cannot handle distinct 14 parameters
                # We need to use a separate LennardJonesForce instead
                # TODO: Can we autodetect this and switch on separate_ljforce earlier?
                if (atom_type.rmin_14 != atom_type.rmin) or (atom_type.epsilon_14 != atom_type.epsilon):
                    raise NotImplementedError('OpenMM <NonbondedForce> cannot handle '
                        'distinct 1-4 sigma and epsilon parameters; '
                        'use separate_ljforce=True instead')

            # Ensure we don't have sigma = 0
            if (sigma == 0.0):
                if (epsilon == 0.0):
                    sigma = 1.0 # reset sigma = 1
                else:
                    raise ValueError("For atom type '%s', sigma = 0 but "
                                     "epsilon != 0." % name)

            etree.SubElement(xml_force, 'Atom', type=name, sigma=str(sigma), epsilon=str(abs(epsilon)))

    @needs_lxml
    def _write_omm_LennardJonesForce(self, xml_root, skip_types, separate_ljforce):
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
        xml_force = etree.SubElement(xml_root, 'LennardJonesForce', lj14scale=str(lj14scale))
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

            # Handle special values used for 14 interactions
            if (atom_type.rmin_14 != atom_type.rmin) or (atom_type.epsilon_14 != atom_type.epsilon):
                sigma14 = atom_type.sigma_14 * length_conv  # in md_unit_system
                epsilon14 = atom_type.epsilon_14 * ene_conv # in md_unit_system

                # Ensure we don't have sigma = 0
                if sigma14 == 0.0:
                    if (epsilon14 == 0.0):
                        sigma14 = 1.0 # reset sigma = 1
                    else:
                        raise ValueError("For atom type '%s', sigma_14 = 0 but "
                                        "epsilon_14 != 0." % name)
            else:
                sigma14 = None
                epsilon14 = None

            attributes = { 'type' : name, 'sigma' : str(sigma), 'epsilon' : str(abs(epsilon)) }
            if epsilon14 is not None:
                attributes['epsilon14'] = str(abs(epsilon14))
            if sigma14 is not None:
                attributes['sigma14'] = str(sigma14)
            etree.SubElement(xml_force, 'Atom', **attributes)

        # write NBFIX records
        for (atom_types, value) in iteritems(self.nbfix_types):
            emin = value[0] * ene_conv
            rmin = value[1] * length_conv
            # convert to sigma; note that NBFIX types are not rmin/2 but rmin
            sigma = rmin/(2**(1.0/6))
            etree.SubElement(xml_force, 'NBFixPair', type1=atom_types[0], type2=atom_types[1], sigma=str(sigma), epsilon=str(emin))

    def _write_omm_scripts(self, dest, skip_types):
        # Not currently implemented, so throw an exception if any unsupported
        # options are specified
        if self.combining_rule == 'geometric':
            raise NotImplementedError('Geometric combining rule not currently supported.')
