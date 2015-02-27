"""
This module contains parsers for mol2-format files (with support for the mol3
extension described at http://q4md-forcefieldtools.org/Tutorial/leap-mol3.php
"""
from chemistry.exceptions import Mol2Error
from chemistry.formats.io import genopen, TextToBinaryFile
from chemistry.formats.registry import FileFormatType
from chemistry.modeller import ResidueTemplate, ResidueTemplateContainer
from chemistry.structure import Structure
from chemistry.topologyobjects import Atom, Bond
import copy

class Mol2File(object):
    """ Class to read and write TRIPOS Mol2 files """
    __metaclass__ = FileFormatType

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
        f = TextFileToBinary(genopen(filename, 'r'))
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
        filename : str
            Name of the file to parse
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
            present where integers or floating point numbers are expected
        """
        f = TextToBinaryFile(genopen(filename, 'r'))
        rescont = ResidueTemplateContainer()
        struct = Structure()
        restemp = ResidueTemplate()
        mol_info = []
        try:
            section = None
            last_residue = None
            headtail = 'head'
            for line in f:
                if line.startswith('#'): continue
                if not line.strip() and section is None: continue
                if line.startswith('@<TRIPOS>'):
                    section = line[9:].strip()
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
                    restemp.add_atom(copy.copy(atom))
                    continue
                if section == 'BOND':
                    # Section formatted as follows:
                    #   bond_id -- serial number of bond (ignored)
                    #   origin_atom_id -- serial number of first atom in bond
                    #   target_atom_id -- serial number of other atom in bond
                    #   bond_type -- string describing bond type (ignored)
                    #   status_bits -- ignored
                    words = line.split()
                    int(words[0]) # Bond serial number... redundant and ignored
                    a1 = int(words[1])
                    a2 = int(words[2])
                    atom1 = struct.atoms.find_original_index(a1)
                    atom2 = struct.atoms.find_original_index(a2)
                    struct.bonds.append(Bond(atom1, atom2))
                    # Now add it to our residue container
                    # See if it's a head/tail connection
                    if atom1.residue is not atom2.residue:
                        if atom1.residue.idx == len(rescont):
                            res1 = restemp
                        elif atom1.residue.idx < len(rescont):
                            res1 = rescont[atom1.residue.idx]
                        else:
                            raise Mol2Error('Bad bonding pattern detected')
                        if atom2.residue.idx == len(rescont):
                            res2 = restemp
                        elif atom1.residue.idx < len(rescont):
                            res2 = rescont[atom2.residue.idx]
                        else:
                            raise Mol2Error('Bad bonding pattern detected')
                        assert res1 is not res2, 'BAD identical residues'
                        idx1 = atom1.idx - atom1.residue[0].idx
                        idx2 = atom2.idx - atom2.residue[0].idx
                        if atom1.residue.idx < atom2.residue.idx:
                            res1.tail = res1[idx1]
                            res2.head = res2[idx2]
                        else:
                            res1.head = res1[idx1]
                            res2.tail = res2[idx2]
                    else:
                        # Same residue, add the bond
                        offset = atom1.residue[0].idx
                        if atom1.residue.idx == len(rescont):
                            res = restemp
                        else:
                            res = rescont[atom1.residue.idx]
                        res.add_bond(atom1.idx-offset, atom2.idx-offset)
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
                    box = [float(x) for x in words[:6]]
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
                    id = int(words[0])
                    resname = words[1]
                    root_atom = int(words[2])
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
                    if residx - 1 == len(rescont):
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
                        # Did not find the atom
                        raise Mol2Error('Could not find %s atom %s in residue '
                                        '%d' % (headtail, atname, residx))
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
                                atom.connections.append(atom)
                                break
                        else:
                            raise Mol2Error('Residue connection atom %s not '
                                            'found in residue %d' % (a, residx))
            if structure:
                return struct
            elif len(rescont) > 0:
                rescont.append(restemp)
                return rescont
            else:
                return restemp
        except ValueError, e:
            raise Mol2Error('String conversion trouble: %s' % e)
        finally:
            f.close()

    #===================================================

    @staticmethod
    def write(struct, filename, mol3=False):
        """ Writes a mol2 file from a structure or residue template

        Parameters
        ----------
        struct : :class:`Structure` or :class:`ResidueTemplate` or
                 :class:`ResidueTemplateContainer`
            The input structure to write the mol2 file from
        filename : str or file-like obj
            Name of the file to write
        mol3 : bool, optional
            If True and ``struct`` is a ResidueTemplate or container, write
            HEAD/TAIL sections. Default is False
        """
        raise NotImplementedError('Not implemented yet')
