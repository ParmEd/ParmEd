"""
This module contains classes for parsing and processing CHARMM parameter,
topology, and stream files. It only extracts atom properties from the
topology files and extracts all parameters from the parameter files

Author: Jason M. Swails
Contributors:
Date: Apr. 13, 2015
"""
from __future__ import division
from parmed import (Atom, AtomType, BondType, AngleType, DihedralType,
                       DihedralTypeList, ImproperType, CmapType, NoUreyBradley)
from parmed.charmm._charmmfile import CharmmFile, CharmmStreamFile
from parmed.exceptions import CharmmError
from parmed.modeller import ResidueTemplate, PatchTemplate
from parmed.parameters import ParameterSet
from parmed.periodic_table import AtomicNum, element_by_mass
from parmed.utils.six import iteritems
from parmed.utils.six.moves import zip
import os
import re
import warnings

_penaltyre = re.compile(r'penalty\s*=\s*([\d\.]+)')

class _EmptyStringIterator(object):
    """ Always yields an empty string """
    def __iter__(self):
        while True:
            yield ''
    def __getitem__(self, idx):
        return ''

class CharmmParameterSet(ParameterSet):
    """
    Stores a parameter set defined by CHARMM files. It stores the equivalent of
    the information found in the MASS section of the CHARMM topology file
    (TOP/RTF) and all of the information in the parameter files (PAR)

    Parameters
    ----------
    filenames : list of str
        The list of topology, parameter, and stream files to load into the
        parameter set. The following file type suffixes are recognized:
            .rtf, .top -- Residue topology file
            .par, .prm -- Parameter file
            .str -- Stream file
            .inp -- If "par" is in the file name, it is a parameter file. If
                    "top" is in the file name, it is a topology file.
                    Otherwise, ValueError is raised.

    Attributes
    ----------
    atom_types : dict(str:AtomType)
        Dictionary mapping the names of the atom types to the corresponding
        AtomType instances
    atom_types_str : dict(str:AtomType)
        alias for atom_types
    atom_types_int : dict(int:AtomType)
        Dictionary mapping the serial indexes of the atom types to the
        corresponding AtomType instances
    atom_types_tuple : dict((str,int):AtomType)
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

    @staticmethod
    def _convert(data, type, msg=''):
        """
        Converts a data type to a desired type, raising CharmmError if it
        fails
        """
        try:
            return type(data)
        except ValueError:
            raise CharmmError('Could not convert %s to %s' % (msg, type))

    def __init__(self, *args):
        # Instantiate the list types
        super(CharmmParameterSet, self).__init__(self)
        self.parametersets = []
        self.residues = dict()
        self.patches = dict()

        # Load all of the files
        tops, pars, strs = [], [], []
        for arg in args:
            if arg.endswith('.rtf') or arg.endswith('.top'):
                tops.append(arg)
            elif arg.endswith('.par') or arg.endswith('.prm'):
                pars.append(arg)
            elif arg.endswith('.str'):
                strs.append(arg)
            elif arg.endswith('.inp'):
                # Only consider the file name (since the directory is likely
                # "toppar" and will screw up file type detection)
                fname = os.path.split(arg)[1]
                if 'par' in fname:
                    pars.append(arg)
                elif 'top' in fname:
                    tops.append(arg)
                else:
                    raise ValueError('Unrecognized file type: %s' % arg)
            else:
                raise ValueError('Unrecognized file type: %s' % arg)
        for top in tops: self.read_topology_file(top)
        for par in pars: self.read_parameter_file(par)
        for strf in strs: self.read_stream_file(strf)

    @classmethod
    def load_set(cls, tfile=None, pfile=None, sfiles=None):
        """
        Instantiates a CharmmParameterSet from a Topology file and a Parameter
        file (or just a Parameter file if it has all information)

        Parameters
        ----------
        tfile : str
            The name of the Topology (RTF/TOP) file to parse
        pfile : str
            The name of the Parameter (PAR/PRM) file to parse
        sfiles : list(str)
            Iterable of stream (STR) file names

        Returns
        -------
        New CharmmParameterSet populated with parameters found in the provided
        files
            
        Notes
        -----
        The RTF file is read first (if provided), followed by the PAR file,
        followed by the list of stream files (in the order they are provided).
        Parameters in each stream file will overwrite those that came before (or
        simply append to the existing set if they are different)
        """
        inst = cls()
        if tfile is not None:
            inst.read_topology_file(tfile)
        if pfile is not None:
            inst.read_parameter_file(pfile)
        if isinstance(sfiles, str):
            # The API docstring requests a list, but allow for users to pass a
            # string with a single filename instead
            inst.read_stream_file(sfiles)
        elif sfiles is not None:
            for sfile in sfiles:
                inst.read_stream_file(sfile)
        return inst

    def read_parameter_file(self, pfile, comments=None):
        """
        Reads all of the parameters from a parameter file. Versions 36 and
        later of the CHARMM force field files have an ATOMS section defining
        all of the atom types.  Older versions need to load this information
        from the RTF/TOP files.

        Parameters
        ----------
        pfile : str or list of lines
            Name of the CHARMM parameter file to read or list of lines to parse
            as a file
        comments : list of str, optional
            List of comments on each of the pfile lines (if pfile is a list of
            lines)

        Notes
        -----
        The atom types must all be loaded by the end of this routine.  Either
        supply a PAR file with atom definitions in them or read in a RTF/TOP
        file first. Failure to do so will result in a raised RuntimeError.
        """
        conv = CharmmParameterSet._convert
        if isinstance(pfile, str):
            own_handle = True
            f = CharmmFile(pfile)
        else:
            own_handle = False
            f = pfile
            if not isinstance(f, CharmmFile) and comments is None:
                comments = _EmptyStringIterator()
        # What section are we parsing?
        section = None
        # The current cmap we are building (these span multiple lines)
        current_cmap = None
        current_cmap2 = None
        current_cmap_data = []
        current_cmap_res = 0
        nonbonded_types = dict() # Holder
        parameterset = None
        read_first_nonbonded = False
        for i, line in enumerate(f):
            line = line.strip()
            try:
                comment = f.comment
            except AttributeError:
                comment = comments[i]
            if not line:
                # This is a blank line
                continue
            if parameterset is None and line.strip().startswith('*>>'):
                parameterset = line.strip()[1:78]
                continue
            # Set section if this is a section header
            if line.startswith('ATOMS'):
                section = 'ATOMS'
                continue
            if line.startswith('BONDS'):
                section = 'BONDS'
                continue
            if line.startswith('ANGLES'):
                section = 'ANGLES'
                continue
            if line.startswith('DIHEDRALS'):
                section = 'DIHEDRALS'
                continue
            if line.startswith('IMPROPER'):
                section = 'IMPROPER'
                continue
            if line.startswith('CMAP'):
                section = 'CMAP'
                continue
            if line.startswith('NONBONDED'):
                read_first_nonbonded = False
                section = 'NONBONDED'
                continue
            if line.startswith('NBFIX'):
                section = 'NBFIX'
                continue
            if line.startswith('HBOND'):
                section = None
                continue
            # It seems like files? sections? can be terminated with 'END'
            if line.startswith('END'): # should this be case-insensitive?
                section = None
                continue
            # If we have no section, skip
            if section is None: continue
            # See if our comments define a penalty for this line
            pens = _penaltyre.findall(comment)
            if len(pens) == 1:
                penalty = float(pens[0])
            else:
                penalty = None
            # Now handle each section specifically
            if section == 'ATOMS':
                if not line.startswith('MASS'): continue # Should this happen?
                words = line.split()
                try:
                    idx = conv(words[1], int, 'atom type')
                    name = words[2]
                    mass = conv(words[3], float, 'atom mass')
                except IndexError:
                    raise CharmmError('Could not parse MASS section.')
                # The parameter file might or might not have an element name
                try:
                    elem = words[4]
                    if len(elem) == 2:
                        elem = elem[0] + elem[1].lower()
                    atomic_number = AtomicNum[elem]
                except (IndexError, KeyError):
                    # Figure it out from the mass
                    atomic_number = AtomicNum[element_by_mass(mass)]
                atype = AtomType(name=name, number=idx, mass=mass,
                                 atomic_number=atomic_number)
                self.atom_types_str[atype.name] = atype
                self.atom_types_int[atype.number] = atype
                self.atom_types_tuple[(atype.name, atype.number)] = atype
                continue
            if section == 'BONDS':
                words = line.split()
                try:
                    type1 = words[0]
                    type2 = words[1]
                    k = conv(words[2], float, 'bond force constant')
                    req = conv(words[3], float, 'bond equilibrium dist')
                except IndexError:
                    raise CharmmError('Could not parse bonds.')
                key = (min(type1, type2), max(type1, type2))
                bond_type = BondType(k, req)
                self.bond_types[(type1, type2)] = bond_type
                self.bond_types[(type2, type1)] = bond_type
                bond_type.penalty = penalty
                continue
            if section == 'ANGLES':
                words = line.split()
                try:
                    type1 = words[0]
                    type2 = words[1]
                    type3 = words[2]
                    k = conv(words[3], float, 'angle force constant')
                    theteq = conv(words[4], float, 'angle equilibrium value')
                except IndexError:
                    raise CharmmError('Could not parse angles.')
                angle_type = AngleType(k, theteq)
                self.angle_types[(type1, type2, type3)] = angle_type
                self.angle_types[(type3, type2, type1)] = angle_type
                # See if we have a urey-bradley
                try:
                    ubk = conv(words[5], float, 'Urey-Bradley force constant')
                    ubeq = conv(words[6], float, 'Urey-Bradley equil. value')
                    ubtype = BondType(ubk, ubeq)
                    ubtype.penalty = penalty
                except IndexError:
                    ubtype = NoUreyBradley
                self.urey_bradley_types[(type1, type2, type3)] = ubtype
                self.urey_bradley_types[(type3, type2, type1)] = ubtype
                angle_type.penalty = penalty
                continue
            if section == 'DIHEDRALS':
                words = line.split()
                try:
                    type1 = words[0]
                    type2 = words[1]
                    type3 = words[2]
                    type4 = words[3]
                    k = conv(words[4], float, 'dihedral force constant')
                    n = conv(words[5], float, 'dihedral periodicity')
                    phase = conv(words[6], float, 'dihedral phase')
                except IndexError:
                    raise CharmmError('Could not parse dihedrals.')
                key = (type1, type2, type3, type4)
                # See if this is a second (or more) term of the dihedral group
                # that's already present.
                dihedral = DihedralType(k, n, phase)
                dihedral.penalty = penalty
                if key in self.dihedral_types:
                    # See if the existing dihedral type list has a term with
                    # the same periodicity -- If so, replace it
                    replaced = False
                    for i, dtype in enumerate(self.dihedral_types[key]):
                        if dtype.per == dihedral.per:
                            # Replace. Warn if they are different
                            if dtype != dihedral:
                                warnings.warn('Replacing dihedral %r with %r' % 
                                              (dtype, dihedral))
                            self.dihedral_types[key][i] = dihedral
                            replaced = True
                            break
                    if not replaced:
                        self.dihedral_types[key].append(dihedral)
                    # Now do the other order
                    replaced = False
                    key = (type4, type3, type2, type1)
                    for i, dtype in enumerate(self.dihedral_types[key]):
                        if dtype.per == dihedral.per:
                            self.dihedral_types[key][i] = dihedral
                            replaced = True
                            break
                    if not replaced:
                        self.dihedral_types[key].append(dihedral)
                else: # key not present
                    dtl = DihedralTypeList()
                    dtl.append(dihedral)
                    self.dihedral_types[(type1,type2,type3,type4)] = dtl
                    self.dihedral_types[(type4,type3,type2,type1)] = dtl
                continue
            if section == 'IMPROPER':
                words = line.split()
                try:
                    type1 = words[0]
                    type2 = words[1]
                    type3 = words[2]
                    type4 = words[3]
                    k = conv(words[4], float, 'improper force constant')
                    theteq = conv(words[5], float, 'improper equil. value')
                except IndexError:
                    raise CharmmError('Could not parse dihedrals.')
                # If we have a 7th column, that is the real psi0 (and the 6th
                # is just a dummy 0)
                try:
                    tmp = conv(words[6], float, 'improper equil. value')
                    theteq = tmp
                except IndexError:
                    pass # Do nothing
                # Improper types seem not to have the central atom defined in
                # the first place, so just have the key a fully sorted list. We
                # still depend on the PSF having properly ordered improper atoms
                key = tuple(sorted([type1, type2, type3, type4]))
                improp = ImproperType(k, theteq)
                self.improper_types[key] = improp
                improp.penalty = penalty
                continue
            if section == 'CMAP':
                # This is the most complicated part, since cmap parameters span
                # many lines. We won't do much error catching here.
                words = line.split()
                try:
                    holder = [float(w) for w in words]
                    current_cmap_data.extend(holder)
                except ValueError:
                    # We assume this is a definition of a new CMAP, so
                    # terminate the last CMAP if applicable
                    if current_cmap is not None:
                        # We have a map to terminate
                        ty = CmapType(current_cmap_res, current_cmap_data)
                        self.cmap_types[current_cmap] = ty
                        self.cmap_types[current_cmap2] = ty
                    try:
                        type1 = words[0]
                        type2 = words[1]
                        type3 = words[2]
                        type4 = words[3]
                        type5 = words[4]
                        type6 = words[5]
                        type7 = words[6]
                        type8 = words[7]
                        res = conv(words[8], int, 'CMAP resolution')
                    except IndexError:
                        raise CharmmError('Could not parse CMAP data.')
                    # order the torsions independently
                    k1 = [type1,type2,type3,type4,type5,type6,type7,type8]
                    k2 = [type8,type7,type6,type5,type4,type3,type2,type1]
                    current_cmap = tuple(min(k1, k2))
                    current_cmap2 = tuple(max(k1, k2))
                    current_cmap_res = res
                    current_cmap_data = []
                continue
            if section == 'NONBONDED':
                # Now get the nonbonded values
                words = line.split()
                try:
                    atype = words[0]
                    # 1st column is ignored
                    epsilon = conv(words[2], float, 'vdW epsilon term')
                    rmin = conv(words[3], float, 'vdW Rmin/2 term')
                except IndexError:
                    # If we haven't read our first nonbonded term yet, we may
                    # just be parsing the settings that should be used. So
                    # soldier on
                    if not read_first_nonbonded: continue
                    raise CharmmError('Could not parse nonbonded terms.')
                except CharmmError:
                    if not read_first_nonbonded: continue
                    raise
                else:
                    # OK, we've read our first nonbonded section for sure now
                    read_first_nonbonded = True
                # See if we have 1-4 parameters
                try:
                    # 4th column is ignored
                    eps14 = conv(words[5], float, '1-4 vdW epsilon term')
                    rmin14 = conv(words[6], float, '1-4 vdW Rmin/2 term')
                except IndexError:
                    eps14 = rmin14 = None
                nonbonded_types[atype] = [epsilon, rmin, eps14, rmin14]
                continue
            if section == 'NBFIX':
                words = line.split()
                try:
                    at1 = words[0]
                    at2 = words[1]
                    emin = abs(conv(words[2], float, 'NBFIX Emin'))
                    rmin = conv(words[3], float, 'NBFIX Rmin')
                    try:
                        emin14 = abs(conv(words[4], float, 'NBFIX Emin 1-4'))
                        rmin14 = conv(words[5], float, 'NBFIX Rmin 1-4')
                    except IndexError:
                        emin14 = rmin14 = None
                    try:
                        self.atom_types_str[at1].add_nbfix(at2, rmin, emin,
                                                           rmin14, emin14)
                        self.atom_types_str[at2].add_nbfix(at1, rmin, emin,
                                                           rmin14, emin14)
                    except KeyError:
                        # Some stream files define NBFIX terms with an atom that
                        # is defined in another toppar file that does not
                        # necessarily have to be loaded. As a result, not every
                        # NBFIX found here will necessarily need to be applied.
                        # If we can't find a particular atom type, don't bother
                        # adding that nbfix and press on
                        pass
                except IndexError:
                    raise CharmmError('Could not parse NBFIX terms.')
                self.nbfix_types[(min(at1, at2), max(at1, at2))] = (emin, rmin)
        # If we had any CMAP terms, then the last one will not have been added
        # yet. Add it here
        if current_cmap is not None:
            ty = CmapType(current_cmap_res, current_cmap_data)
            self.cmap_types[current_cmap] = ty
            self.cmap_types[current_cmap2] = ty
        # Now we're done. Load the nonbonded types into the relevant AtomType
        # instances. In order for this to work, all keys in nonbonded_types
        # must be in the self.atom_types_str dict. Raise a RuntimeError if this
        # is not satisfied
        try:
            for key in nonbonded_types:
                self.atom_types_str[key].set_lj_params(*nonbonded_types[key])
        except KeyError:
            raise RuntimeError('Atom type %s not present in AtomType list' %
                               key)
        if parameterset is not None: self.parametersets.append(parameterset)
        if own_handle: f.close()

    def read_topology_file(self, tfile):
        """
        Reads _only_ the atom type definitions from a topology file. This is
        unnecessary for versions 36 and later of the CHARMM force field.

        Parameters
        ----------
        tfile : str
            Name of the CHARMM topology file to read
        """
        conv = CharmmParameterSet._convert
        if isinstance(tfile, str):
            own_handle = True
            f = iter(CharmmFile(tfile))
        else:
            own_handle = False
            f = tfile
        hpatch = tpatch = None # default Head and Tail patches
        residues = dict()
        patches = dict()
        hpatches = dict()
        tpatches = dict()
        line = next(f)
        try:
            while line:
                line = line.strip()
                if line[:4] == 'MASS':
                    words = line.split()
                    try:
                        idx = conv(words[1], int, 'atom type')
                        name = words[2]
                        mass = conv(words[3], float, 'atom mass')
                    except IndexError:
                        raise CharmmError('Could not parse MASS section of %s' %
                                          tfile)
                    # The parameter file might or might not have an element name
                    try:
                        elem = words[4]
                        if len(elem) == 2:
                            elem = elem[0] + elem[1].lower()
                        atomic_number = AtomicNum[elem]
                    except (IndexError, KeyError):
                        # Figure it out from the mass
                        atomic_number = AtomicNum[element_by_mass(mass)]
                    atype = AtomType(name=name, number=idx, mass=mass,
                                     atomic_number=atomic_number)
                    self.atom_types_str[atype.name] = atype
                    self.atom_types_int[atype.number] = atype
                    self.atom_types_tuple[(atype.name, atype.number)] = atype
                elif line[:4] == 'DECL':
                    pass # Not really sure what this means
                elif line[:4] == 'DEFA':
                    words = line.split()
                    if len(words) < 5:
                        warnings.warn('DEFA line has %d tokens; expected 5' %
                                      len(words))
                    else:
                        it = iter(words[1:5])
                        for tok, val in zip(it, it):
                            if val.upper() == 'NONE':
                                val = None
                            if tok.upper().startswith('FIRS'):
                                hpatch = val
                            elif tok.upper() == 'LAST':
                                tpatch = val
                            else:
                                warnings.warn('DEFA patch %s unknown' % val)
                elif line[:4].upper() in ('RESI', 'PRES'):
                    restype = line[:4].upper()
                    # Get the residue definition
                    words = line.split()
                    resname = words[1]
                    # Assign default patches
                    hpatches[resname] = hpatch
                    tpatches[resname] = tpatch
                    try:
                        charge = float(words[2])
                    except (IndexError, ValueError):
                        warnings.warn('No charge for %s' % resname)
                    if restype == 'RESI':
                        res = ResidueTemplate(resname)
                    elif restype == 'PRES':
                        res = PatchTemplate(resname)
                    else:
                        assert False, 'restype != RESI or PRES'
                    line = next(f)
                    group = []
                    ictable = []
                    while line:
                        if line[:5].upper() == 'GROUP':
                            if group:
                                res.groups.append(group)
                            group = []
                        elif line[:4].upper() == 'ATOM':
                            words = line.split()
                            name, type, charge = words[1:4]
                            charge = float(charge)
                            atom = Atom(name=name, type=type, charge=charge)
                            group.append(atom)
                            res.add_atom(atom)
                        elif line.strip() and line.split()[0].upper() in \
                                ('BOND', 'DOUBLE'):
                            words = line.split()[1:]
                            it = iter(words)
                            for a1, a2 in zip(it, it):
                                if a1.startswith('-'):
                                    res.head = res[a2]
                                    continue
                                if a2.startswith('-'):
                                    res.head = res[a1]
                                    continue
                                if a1.startswith('+'):
                                    res.tail = res[a2]
                                    continue
                                if a2.startswith('+'):
                                    res.tail = res[a1]
                                    continue
                                # Apparently PRES objects do not need to put +
                                # or - in front of atoms that belong to adjacent
                                # residues
                                if restype == 'PRES' and (a1 not in res or
                                                          a2 not in res):
                                    continue
                                res.add_bond(a1, a2)
                        elif line[:4].upper() == 'CMAP':
                            pass
                        elif line[:5].upper() == 'DONOR':
                            pass
                        elif line[:6].upper() == 'ACCEPT':
                            pass
                        elif line[:2].upper() == 'IC':
                            w = line.split()[1:]
                            ictable.append((w[:4], [float(x) for x in w[4:]]))
                        elif line[:3].upper() == 'END':
                            break
                        elif line[:5].upper() == 'PATCH':
                            it = iter(line.split()[1:])
                            for tok, val in zip(it, it):
                                if val.upper() == 'NONE': val = None
                                if tok.upper().startswith('FIRS'):
                                    hpatches[resname] = val
                                elif tok.upper().startswith('LAST'):
                                    tpatches[resname] = val
                        elif line[:5].upper() == 'DELETE':
                            pass
                        elif line[:4].upper() == 'IMPR':
                            it = iter(line.split()[1:])
                            for a1, a2, a3, a4 in zip(it, it, it, it):
                                if a2[0] == '-' or a3[0] == '-' or a4 == '-':
                                    res.head = res[a1]
                        elif line[:4].upper() in ('RESI', 'PRES', 'MASS'):
                            # Back up a line and bail
                            break
                        line = next(f)
                    if group: res.groups.append(group)
                    _fit_IC_table(res, ictable)
                    if restype == 'RESI':
                        residues[resname] = res
                    elif restype == 'PRES':
                        patches[resname] = res
                    else:
                        assert False, 'restype != RESI or PRES'
                    # We parsed a line we need to look at. So don't update the
                    # iterator
                    continue
                # Get the next line and cycle through
                line = next(f)
        except StopIteration:
            pass

        # Go through the patches and add the appropriate one
        for resname, res in iteritems(residues):
            if hpatches[resname] is not None:
                try:
                    res.first_patch = patches[hpatches[resname]]
                except KeyError:
                    warnings.warn('Patch %s not found' % hpatches[resname])
            if tpatches[resname] is not None:
                try:
                    res.last_patch = patches[tpatches[resname]]
                except KeyError:
                    warnings.warn('Patch %s not found' % tpatches[resname])
        # Now update the residues and patches with the ones we parsed here
        self.residues.update(residues)
        self.patches.update(patches)

        if own_handle: f.close()

    def read_stream_file(self, sfile):
        """
        Reads RTF and PAR sections from a stream file and dispatches the
        sections to read_topology_file or read_parameter_file

        Parameters
        ----------
        sfile : str or CharmmStreamFile
            Stream file to parse
        """
        if isinstance(sfile, CharmmStreamFile):
            f = sfile
        else:
            f = CharmmStreamFile(sfile)

        title, section, comments = f.next_section()
        while title is not None and section is not None:
            words = title.lower().split()
            if words[1] == 'rtf':
                # This is a Residue Topology File section.
                self.read_topology_file(iter(section))
            elif words[1].startswith('para'):
                # This is a Parameter file section
                self.read_parameter_file(section, comments)
            title, section, comments = f.next_section()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _fit_IC_table(res, ictable):
    """
    Determines cartesian coordinates from an internal coordinate table stored in
    CHARMM residue topology files

    Parameters
    ----------
    res : ResidueTemplate
        The residue template for which coordinates are being determined
    ictable : list[tuple(atoms, measurements)]
        The internal coordinate table

    Notes
    -----
    This method assigns an xx, xy, and xz attribute to ``res``. For the time
    being, this is just a placeholder, as its functionality has not yet been
    implemented (CHARMM does not use a 'traditional' Z-matrix, and I don't know
    of any existing code that will compute the proper cartesian coordinates from
    the form used by CHARMM)
    """
    for atom in res:
        atom.xx = atom.xy = atom.xz = 0.0
