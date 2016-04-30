"""
This module contains parsers for mol2-format files (with support for the mol3
extension described at http://q4md-forcefieldtools.org/Tutorial/leap-mol3.php
"""
from __future__ import print_function, division, absolute_import

import copy
from parmed.exceptions import Mol2Error, ParameterWarning
from parmed.formats.registry import FileFormatType
from parmed.modeller.residue import ResidueTemplate, ResidueTemplateContainer
from parmed.periodic_table import element_by_name, AtomicNum
from parmed.residue import AminoAcidResidue, RNAResidue, DNAResidue
from parmed.structure import Structure
from parmed.topologyobjects import Atom, Bond
from parmed.utils.io import genopen
from parmed.utils.six import add_metaclass, string_types
import warnings

@add_metaclass(FileFormatType)
class Mol2File(object):
    """ Class to read and write TRIPOS Mol2 files """

    BOND_ORDER_MAP = dict(ar=1.5, am=1.25)
    REVERSE_BOND_ORDER_MAP = {1.25 : 'am', 1.5 : 'ar'}

    #===================================================

    @staticmethod
    def id_format(filename):
        """ Identify the file as a Mol2 (or Mol3) file format or not

        Parameters
        ----------
        filename : str
            Name of the file to test whether or not it is a mol2 file

        Returns
        -------
        is_fmt : bool
            True if it is a mol2 (or mol3) file, False otherwise
        """
        f = genopen(filename, 'r')
        try:
            for line in f:
                if line.startswith('#'): continue
                if not line.strip(): continue
                return line.startswith('@<TRIPOS>')
            return False
        finally:
            f.close()

    #===================================================

    @staticmethod
    def parse(filename, structure=False):
        """ Parses a mol2 file (or mol3) file

        Parameters
        ----------
        filename : str or file-like
            Name of the file to parse or file-like object to parse from
        structure : bool, optional
            If True, the return value is a :class:`Structure` instance. If
            False, it is either a :class:`ResidueTemplate` or
            :class:`ResidueTemplateContainter` instance, depending on whether
            there is one or more than one residue defined in it. Default is
            False

        Returns
        -------
        molecule : :class:`Structure`, :class:`ResidueTemplate`, or
                   :class:`ResidueTemplateContainer`
            The molecule defined by this mol2 file

        Raises
        ------
        Mol2Error
            If the file format is not recognized or non-numeric values are
            present where integers or floating point numbers are expected. Also
            raises Mol2Error if you try to parse a mol2 file that has multiple
            @<MOLECULE> entries with ``structure=True``.
        """
        if isinstance(filename, string_types):
            f = genopen(filename, 'r')
            own_handle = True
        else:
            f = filename
            own_handle = False
        rescont = ResidueTemplateContainer()
        struct = Structure()
        restemp = ResidueTemplate()
        mol_info = []
        multires_structure = False
        try:
            section = None
            last_residue = None
            headtail = 'head'
            molecule_number = 0
            for line in f:
                if line.startswith('#'): continue
                if not line.strip() and section is None: continue
                if line.startswith('@<TRIPOS>'):
                    section = line[9:].strip()
                    if section == 'MOLECULE' and (restemp.atoms or rescont):
                        if structure:
                            raise Mol2Error('Cannot convert MOL2 with multiple '
                                            '@<MOLECULE>s to a Structure')
                        # Set the residue name from the MOL2 title if the
                        # molecule had only 1 residue and it was given a name in
                        # the title
                        if not multires_structure and mol_info[0]:
                            restemp.name = mol_info[0]
                        multires_structure = False
                        rescont.append(restemp)
                        restemp = ResidueTemplate()
                        struct = Structure()
                        last_residue = None
                        molecule_number += 1
                        mol_info = []
                    continue
                if section is None:
                    raise Mol2Error('Bad mol2 file format')
                if section == 'MOLECULE':
                    # Section formatted as follows:
                    #   mol_name
                    #   num_atoms [num_bonds [num_substr [num_feat [num_sets]]]]
                    #   mol_type
                    #   charge_type
                    #   [status_bits]
                    #   [mol_comment]
                    # TODO: Do something with the name.
                    if len(mol_info) == 0:
                        mol_info.append(line.strip())
                    elif len(mol_info) == 1:
                        mol_info.append([int(x) for x in line.split()])
                    elif len(mol_info) == 2:
                        mol_info.append(line.strip())
                    elif len(mol_info) == 3:
                        mol_info.append(line.strip())
                    # Ignore the rest
                    continue
                if section == 'ATOM':
                    # Section formatted as follows:
                    #   atom_id -- serial number of atom
                    #   atom_name -- name of the atom
                    #   x -- X-coordinate of the atom
                    #   y -- Y-coordinate of the atom
                    #   z -- Z-coordinate of the atom
                    #   atom_type -- type of the atom
                    #   subst_id -- Residue serial number
                    #   subst_name -- Residue name
                    #   charge -- partial atomic charge
                    #   status_bit -- ignored
                    words = line.split()
                    id = int(words[0])
                    name = words[1]
                    x = float(words[2])
                    y = float(words[3])
                    z = float(words[4])
                    typ = words[5]
                    try:
                        resid = int(words[6])
                    except IndexError:
                        resid = 0
                    try:
                        resname = words[7]
                    except IndexError:
                        resname = 'UNK'
                    if 'NO_CHARGES' not in mol_info:
                        try:
                            charge = float(words[8])
                        except IndexError:
                            charge = 0
                    else:
                        charge = 0
                    if last_residue is None:
                        last_residue = (resid, resname)
                        restemp.name = resname
                    atom = Atom(name=name, type=typ, number=id, charge=charge)
                    atom.xx, atom.xy, atom.xz = x, y, z
                    struct.add_atom(atom, resname, resid)
                    if last_residue != (resid, resname):
                        rescont.append(restemp)
                        restemp = ResidueTemplate()
                        restemp.name = resname
                        last_residue = (resid, resname)
                        multires_structure = True
                    try:
                        restemp.add_atom(copy.copy(atom))
                    except ValueError:
                        # Allow mol2 files being parsed as a Structure to have
                        # duplicate atom names
                        if not structure:
                            raise
                    continue
                if section == 'BOND':
                    # Section formatted as follows:
                    #   bond_id -- serial number of bond (ignored)
                    #   origin_atom_id -- serial number of first atom in bond
                    #   target_atom_id -- serial number of other atom in bond
                    #   bond_type -- string describing bond type
                    #   status_bits -- ignored
                    words = line.split()
                    int(words[0]) # Bond serial number... redundant and ignored
                    a1 = int(words[1])
                    a2 = int(words[2])
                    try:
                        order = words[3]
                    except IndexError:
                        order = 1.0
                    if order in Mol2File.BOND_ORDER_MAP:
                        order = Mol2File.BOND_ORDER_MAP[order]
                    else:
                        try:
                            order = float(order)
                        except ValueError:
                            warnings.warn('Mol2 bond order not recognized: %s' %
                                          order, ParameterWarning)
                            order = 1.0
                    atom1 = struct.atoms.find_original_index(a1)
                    atom2 = struct.atoms.find_original_index(a2)
                    struct.bonds.append(Bond(atom1, atom2, order=order))
                    # Now add it to our residue container
                    # See if it's a head/tail connection
                    if atom1.residue is not atom2.residue:
                        if atom1.residue.idx == len(rescont):
                            res1 = restemp
                        elif atom1.residue.idx < len(rescont):
                            res1 = rescont[atom1.residue.idx]
                        assert atom.residue.idx <= len(rescont), 'Bad bond!'
                        if atom2.residue.idx == len(rescont):
                            res2 = restemp
                        elif atom2.residue.idx < len(rescont):
                            res2 = rescont[atom2.residue.idx]
                        assert atom.residue.idx <= len(rescont), 'Bad bond!'
                        assert res1 is not res2, 'BAD identical residues'
                        idx1 = atom1.idx - atom1.residue[0].idx
                        idx2 = atom2.idx - atom2.residue[0].idx
                        if atom1.residue.idx < atom2.residue.idx:
                            res1.tail = res1[idx1]
                            res2.head = res2[idx2]
                        else:
                            res1.head = res1[idx1]
                            res2.tail = res2[idx2]
                    elif not multires_structure:
                        if not structure:
                            restemp.add_bond(a1-1, a2-1, order)
                    else:
                        # Same residue, add the bond
                        offset = atom1.residue[0].idx
                        if atom1.residue.idx == len(rescont):
                            res = restemp
                        else:
                            res = rescont[atom1.residue.idx]
                        res.add_bond(atom1.idx-offset, atom2.idx-offset, order)
                    continue
                if section == 'CRYSIN':
                    # Section formatted as follows:
                    #   a -- length of first unit cell vector
                    #   b -- length of second unit cell vector
                    #   c -- length of third unit cell vector
                    #   alpha -- angle b/w b and c
                    #   beta -- angle b/w a and c
                    #   gamma -- angle b/w a and b
                    #   space group -- number of space group (ignored)
                    #   space group setting -- ignored
                    words = line.split()
                    box = [float(w) for w in words[:6]]
                    if len(box) != 6:
                        raise ValueError('%d box dimensions found; needed 6' %
                                         len(box))
                    struct.box = copy.copy(box)
                    rescont.box = copy.copy(box)
                    continue
                if section == 'SUBSTRUCTURE':
                    # Section formatted as follows:
                    #   subst_id -- residue number
                    #   subst_name -- residue name
                    #   root_atom -- first atom of residue
                    #   subst_type -- ignored (usually 'RESIDUE')
                    #   dict_type -- type of substructure (ignored)
                    #   chain -- chain ID of residue
                    #   sub_type -- type of the chain
                    #   inter_bonds -- # of inter-substructure bonds
                    #   status -- ignored
                    #   comment -- ignored
                    words = line.split()
                    if not words: continue
                    id = int(words[0])
                    resname = words[1]
                    try:
                        chain = words[5]
                    except IndexError:
                        chain = ''
                    # Set the chain ID
                    for res in struct.residues:
                        if res.number == id and res.name == resname:
                            res.chain = chain
                    continue
                # MOL3 sections
                if section == 'HEADTAIL':
                    atname, residx = line.split()
                    residx = int(residx)
                    if residx in (0, 1) or residx - 1 == len(rescont):
                        res = restemp
                    elif residx - 1 < len(rescont):
                        res = rescont[residx-1]
                    else:
                        raise Mol2Error('Residue out of range in head/tail')
                    for atom in res:
                        if atom.name == atname:
                            if headtail == 'head':
                                res.head = atom
                                headtail = 'tail'
                            else:
                                res.tail = atom
                                headtail = 'head'
                            break
                    else:
                        if headtail == 'head':
                            headtail = 'tail'
                        else:
                            headtail = 'head'
                    continue
                if section == 'RESIDUECONNECT':
                    words = line.split()
                    residx = int(words[0])
                    if residx - 1 == len(rescont):
                        res = restemp
                    elif residx - 1 < len(rescont):
                        res = rescont[residx-1]
                    else:
                        raise Mol2Error('Residue out of range in '
                                        'residueconnect')
                    for a in words[3:]:
                        if a == '0': continue
                        for atom in res:
                            if atom.name == a:
                                res.connections.append(atom)
                                break
                        else:
                            raise Mol2Error('Residue connection atom %s not '
                                            'found in residue %d' % (a, residx))
            if structure:
                for atom in struct.atoms:
                    anum = _guess_atomic_number(atom.name, restemp)
                    if anum == 0:
                        anum = _guess_atomic_number(atom.type, restemp)
                    atom.atomic_number = anum
                return struct
            elif len(rescont) > 0:
                if not multires_structure and mol_info[0]:
                    restemp.name = mol_info[0]
                rescont.append(restemp)
                for res in rescont:
                    for atom in res.atoms:
                        anum = _guess_atomic_number(atom.name, restemp)
                        if anum == 0:
                            anum = _guess_atomic_number(atom.type, restemp)
                        atom.atomic_number = anum
                return rescont
            else:
                for atom in restemp.atoms:
                    anum = _guess_atomic_number(atom.name, restemp)
                    if anum == 0:
                        anum = _guess_atomic_number(atom.type, restemp)
                    atom.atomic_number = anum
                return restemp
        except ValueError as e:
            raise Mol2Error('String conversion trouble: %s' % e)
        finally:
            if own_handle: f.close()

    #===================================================

    @staticmethod
    def write(struct, dest, mol3=False, split=False, compress_whitespace=False):
        """ Writes a mol2 file from a structure or residue template

        Parameters
        ----------
        struct : :class:`Structure` or :class:`ResidueTemplate` or
                 :class:`ResidueTemplateContainer`
            The input structure to write the mol2 file from
        dest : str or file-like obj
            Name of the file to write or open file handle to write to
        mol3 : bool, optional
            If True and ``struct`` is a ResidueTemplate or container, write
            HEAD/TAIL sections. Default is False
        split : bool, optional
            If True and ``struct`` is a ResidueTemplateContainer or a Structure
            with multiple residues, each residue is printed in a separate
            @<MOLECULE> section that appear sequentially in the output file
        compress_whitespace : bool, optional
            If True, seprate fields on one line with a single space instead of
            aligning them with whitespace. This is useful for parsers that
            truncate lines at 80 characters (e.g., some versions of OpenEye).
            However, it will not look as "neat" upon visual inspection in a text
            editor. Default is False.
        """
        own_handle = False
        if not hasattr(dest, 'write'):
            own_handle = True
            dest = genopen(dest, 'w')
        if split:
            # Write sequentially if it is a multi-residue container or Structure
            if isinstance(struct, ResidueTemplateContainer):
                try:
                    for res in struct:
                        Mol2File.write(res, dest, mol3,
                                       compress_whitespace=compress_whitespace)
                finally:
                    if own_handle: dest.close()
                return
            elif isinstance(struct, Structure) and len(struct.residues) > 1:
                try:
                    for res in ResidueTemplateContainer.from_structure(struct):
                        Mol2File.write(res, dest, mol3,
                                       compress_whitespace=compress_whitespace)
                finally:
                    if own_handle: dest.close()
                return
        try:
            if isinstance(struct, ResidueTemplateContainer):
                natom = sum([len(c) for c in struct])
                # To find the number of bonds, we need to total number of bonds
                # + the number of bonds that would be formed by "stitching"
                # together residues via their head and tail
                bonds = []
                charges = []
                bases = [1 for res in struct]
                for i, res in enumerate(struct):
                    if i < len(struct) - 1:
                        bases[i+1] = bases[i] + len(res)
                for i, res in enumerate(struct):
                    for bond in res.bonds:
                        bonds.append((bond.atom1.idx+bases[i],
                                      bond.atom2.idx+bases[i],
                                      bond.order))
                    if i < len(struct)-1 and (res.tail is not None and
                            struct[i+1].head is not None):
                        bonds.append((res.tail.idx+bases[i],
                                      struct[i+1].head.idx+bases[i+1],
                                      bond.order))
                    charges.extend([a.charge for a in res])
                residues = struct
                name = struct.name or struct[0].name
            else:
                natom = len(struct.atoms)
                bonds = [(b.atom1.idx+1, b.atom2.idx+1, b.order)
                            for b in struct.bonds]
                if isinstance(struct, ResidueTemplate):
                    residues = [struct]
                    name = struct.name
                else:
                    residues = struct.residues
                    name = struct.residues[0].name
                charges = [a.charge for a in struct.atoms]
            dest.write('@<TRIPOS>MOLECULE\n')
            dest.write('%s\n' % name)
            dest.write('%d %d %d 0 1\n' % (natom, len(bonds), len(residues)))
            if len(residues) == 1:
                dest.write('SMALL\n')
            else:
                for residue in residues:
                    if AminoAcidResidue.has(residue.name):
                        dest.write('PROTEIN\n')
                        break
                    if (RNAResidue.has(residue.name) or
                            DNAResidue.has(residue.name)):
                        dest.write('NUCLEIC\n')
                        break
                else:
                    dest.write('BIOPOLYMER\n')
            if not any(charges):
                dest.write('NO_CHARGES\n')
                printchg = False
            else:
                dest.write('USER_CHARGES\n')
                printchg = True
            # See if we want to print box info
            if hasattr(struct, 'box') and struct.box is not None:
                box = struct.box
                dest.write('@<TRIPOS>CRYSIN\n')
                if compress_whitespace:
                    fmt = '%.4f %.4f %.4f %.4f %.4f %.4f 1 1\n'
                else:
                    fmt = '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f  1  1\n'
                dest.write(fmt % tuple(box))
            # Now do ATOM section
            dest.write('@<TRIPOS>ATOM\n')
            j = 1
            for i, res in enumerate(residues):
                for atom in res:
                    try:
                        x = atom.xx
                    except AttributeError:
                        x = 0
                    try:
                        y = atom.xy
                    except AttributeError:
                        y = 0
                    try:
                        z = atom.xz
                    except AttributeError:
                        z = 0
                    if compress_whitespace:
                        fmt = '%d %s %.4f %.4f %.4f %s %d %s'
                    else:
                        fmt = '%8d %-8s %10.4f %10.4f %10.4f %-8s %6d %-8s'
                    dest.write(fmt % (j, atom.name, x, y, z,
                               atom.type.strip() or atom.name, i+1, res.name))
                    if printchg:
                        if compress_whitespace:
                            fmt = ' %.6f\n'
                        else:
                            fmt = ' %10.6f\n'
                        dest.write(fmt % atom.charge)
                    else:
                        dest.write('\n')
                    j += 1
            dest.write('@<TRIPOS>BOND\n')
            for i, bond in enumerate(bonds):
                if bond[2] in Mol2File.REVERSE_BOND_ORDER_MAP:
                    order = Mol2File.REVERSE_BOND_ORDER_MAP[bond[2]]
                else:
                    order = int(bond[2])
                if compress_whitespace:
                    fmt = '%d %d %d %s\n'
                else:
                    fmt = '%8d %8d %8d %s\n'
                dest.write(fmt % (i+1, bond[0], bond[1], order))
            dest.write('@<TRIPOS>SUBSTRUCTURE\n')
            first_atom = 0
            for i, res in enumerate(residues):
                if not hasattr(res, 'chain') or not res.chain:
                    chain = '****'
                else:
                    chain = res.chain
                intresbonds = 0
                if isinstance(res, ResidueTemplate):
                    if i != len(residues)-1 and (res.tail is not None and
                            residues[i+1].head is not None):
                        intresbonds += 1
                    if i != 0 and (res.head is not None and residues[i-1].tail
                            is not None):
                        intresbonds += 1
                else:
                    for atom in res:
                        for a2 in atom.bond_partners:
                            if a2.residue is not res:
                                intresbonds += 1
                if compress_whitespace:
                    fmt = '%d %s %d RESIDUE %d %s ROOT %d\n'
                else:
                    fmt = '%8d %-8s %8d RESIDUE %4d %-4s ROOT %6d\n'
                dest.write(fmt % (i+1, res.name, first_atom+1, 0, chain[:4],
                           intresbonds))
                first_atom += len(res)
            if mol3:
                dest.write('@<TRIPOS>HEADTAIL\n')
                for i, res in enumerate(residues):
                    if isinstance(res, ResidueTemplate):
                        if res.head is not None:
                            dest.write('%s %d\n' % (res.head.name, i+1))
                        else:
                            dest.write('0 0\n')
                        if res.tail is not None:
                            dest.write('%s %d\n' % (res.tail.name, i+1))
                        else:
                            dest.write('0 0\n')
                    else:
                        head = tail = None
                        for atom in res:
                            for a2 in atom.bond_partners:
                                if a2.residue.idx == res.idx - 1:
                                    head = atom
                                if a2.residue.idx == res.idx + 1:
                                    tail = atom
                        if head is not None:
                            dest.write('%s %d\n' % (head.name, i+1))
                        else:
                            dest.write('0 0\n')
                        if tail is not None:
                            dest.write('%s %d\n' % (tail.name, i+1))
                        else:
                            dest.write('0 0\n')
                dest.write('@<TRIPOS>RESIDUECONNECT\n')
                for i, res in enumerate(residues):
                    if isinstance(res, ResidueTemplate):
                        con = [res.head, res.tail, None, None, None, None]
                        for i, a in enumerate(res.connections):
                            con[i+2] = a
                    else:
                        con = [None, None, None, None, None, None]
                        ncon = 2
                        for atom in res:
                            for a2 in atom.bond_partners:
                                if a2.residue.idx == res.idx - 1:
                                    con[0] = atom
                                elif a2.residue.idx == res.idx + 1:
                                    con[1] = atom
                                elif a2.residue.idx != res.idx:
                                    con[ncon] = atom
                                    ncon += 1
                    dest.write('%d' % (i+1))
                    for a in con:
                        if a is not None:
                            dest.write(' %s' % a.name)
                        else:
                            dest.write(' 0')
                    dest.write('\n')
        finally:
            if own_handle: dest.close()

def _guess_atomic_number(name, residue=None):
    """ Guesses the atomic number """
    # Special-case single-atom residues, which are almost always ions
    name = ''.join(c for c in name if c.isalpha())
    if residue is None or len(residue.atoms) == 1:
        if len(name) > 1:
            try:
                return AtomicNum[name[0].upper() + name[1].lower()]
            except KeyError:
                return AtomicNum[element_by_name(name)]
    return AtomicNum[element_by_name(name)]
