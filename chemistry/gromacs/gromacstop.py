"""
This module contains functionality relevant to loading a GROMACS topology file
and building a Structure from it
"""
from chemistry.exceptions import GromacsTopologyError, GromacsTopologyWarning
from chemistry.formats.io import genopen, TextToBinaryFile
from chemistry.formats.registry import FileFormatType
from chemistry.gromacs._gromacsfile import GromacsFile
from chemistry.gromacs import _cpp as cpp
from chemistry.structure import Structure
from chemistry.topologyobjects import (Atom, Bond, Angle, Dihedral, TrackedList,
            NonbondedException)
from chemistry.periodic_table import element_by_mass, AtomicNum
from chemistry import unit as u
import warnings

_ppre = re.compile(r'#([A-Za-z]+) *(\S+)?')

class GromacsTopologyFile(Structure):
    """
    Loads a GROMACS topology file

    Parameters
    ----------
    fname : str=None
        Name of the GROMACS topology file to parse, if any

    Attributes
    ----------
    parameterset : ParameterSet
        The set of parameters defining a force field
    """
    __metaclass__ = FileFormatType

    #===================================================

    @staticmethod
    def id_format(filename):
        """ Identifies the file as a GROMACS topology file

        Parameters
        ----------
        filename : str
            Name of the file to check if it is a gromacs topology file

        Returns
        -------
        is_fmt : bool
            If it is identified as a gromacs topology, return True. False
            otherwise
        """
        f = TextToBinaryFile(genopen(filename))
        try:
            for line in f:
                if line.startswith(';'): continue
                if line.startswith('#'):
                    if line.startswith('#if'): continue
                    if line.startswith('#define'): continue
                    if line.startswith('#include'): continue
                    if line.startswith('#undef'): continue
                    return False
                if line.strip() == '[ moleculetype ]': return True
                return False
            return False
        finally:
            f.close()

    #===================================================

    def __init__(self, fname=None, defines=None):
        super(GromacsTopologyFile, self).__init__()
        self.name = fname
        self.parameterset = ParameterSet()
        # This protects us against using topologies defining a functional form
        # that I have no idea how to deal with
        self.unknown_functional = False
        if fname is not None:
            self.rdparm(fname, defines)

    #===================================================

    def rdparm(self, fname, defines=None):
        """
        Reads the GROMACS topology file

        Parameters
        ----------
        fname : str
            The name of the file to read
        defines : list of str=None
            If specified, this is the set of defines to use when parsing the
            topology file
        """
        params = inst.parameterset
        f = GromacsFile(fname)
        if defines is None:
            defines = []
        try:
            current_section = None
            current_define = None
            negated = False # Is the current define negated?
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # See if we ignore this line based on the current define
                ignore_line = define is not None and (
                        (current_define in defines and not negated) or
                        (current_define not in defines and negated)
                )
                elif line[0] == '#':
                    rematch = _ppre.match(line)
                    if not rematch:
                        raise GromacsTopologyError('Trouble matching command '
                                                   'directive [%s]' % line)
                    cmd, arg = rematch.groups()
                    if cmd == 'else':
                        if current_define is None:
                            raise GromacsTopologyError('Found #else without '
                                                       'matching #if')
                        negated = True
                    elif cmd == 'endif':
                        if arg:
                            warnings.warn('Extra tokens after #endif ignored')
                        current_define = None
                    elif cmd == 'ifdef':
                        if not arg:
                            raise GromacsTopologyError('#ifdef what??')
                        current_define = 
                    elif cmd == 'include':
                        if ignore_line: continue # We don't want to include this
                        incfname = arg.replace('"', '').replace("'", '').strip()
                        inst.rdparm(inst, _get_absolute_path(incfname), defines)
                    else:
                        raise GromacsTopologyError('#%s directive unsupported'
                                                    % cmd)
                elif line[0] == '[' and not ignore_line:
                    current_section = line.replace('[', '').replace(']', '')
                    current_section = current_section.strip()
                elif current_section == 'moleculetype':
                    inst.molecule_type = line.split()
                elif current_section == 'atoms':
                    words = line.split()
                    if len(words) < 8:
                        mass = -1
                        atomic_number = -1
                    else:
                        mass = float(words[7])
                        atomic_number = AtomicNum[element_by_mass(mass)]
                    if len(words) < 7:
                        charge = None
                    else:
                        charge = float(words[6])
                    atom = Atom(atomic_number=atomic_number, name=words[4],
                                type=words[1], charge=charge)
                    inst.residues.add_atom(atom, words[3], int(words[2]))
                    inst.atoms.append(atom)
                elif current_section == 'bonds':
                    words = line.split()
                    i, j = int(words[0]), int(words[1])
                    if words[2] != '1':
                        warnings.warn('bond funct != 1; unknown functional',
                                      GromacsTopologyWarning)
                        inst.unknown_functional = True
                    inst.bonds.append(Bond(inst.atoms[i-1], inst.atoms[j-1]))
                elif current_section == 'pairs':
                    words = line.split()
                    i, j = int(words[0]), int(words[1])
                    if words[2] != '1':
                        warnings.warn('pairs funct != 1; unknown functional',
                                      GromacsTopologyWarning)
                        inst.unknown_functional = True
                    inst.adjusts.append(NonbondedException(inst.atoms[i-1],
                                                           inst.atoms[j-1]))
                elif current_section == 'angles':
                    words = line.split()
                    i, j, k = int(words[0]), int(words[1]), int(words[2])
                    if words[3] != '1':
                        warnings.warn('angles funct != 1; unknown functional',
                                      GromacsTopologyWarning)
                        inst.unknown_functional = True
                    inst.angles.append(Angle(inst.atoms[i-1], inst.atoms[j-1],
                                             inst.atoms[k-1])
                    )
                elif current_section == 'dihedrals':
                    words = line.split()
                    i, j, k, l = [int(x) for x in words[:4]]
                    if words[4] == '1':
                        # Normal dihedral
                        dih = Dihedral(inst.atoms[i-1], inst.atoms[j-1],
                                       inst.atoms[k-1], inst.atoms[l-1])
                        inst.dihedrals.append(dih)
                    elif words[4] == '2':
                        # Improper
                        imp = Improper(inst.atoms[i-1], inst.atoms[j-1],
                                       inst.atoms[k-1], inst.atoms[l-1])
                        inst.impropers.append(imp)
                    else:
                        # ??? unknown
                        warnings.warn('dihedrals funct != 1 or 2; unknown '
                                      'functional', GromacsTopologyWarning)
                        dih = Dihedral(inst.atoms[i-1], inst.atoms[j-1],
                                       inst.atoms[k-1], inst.atoms[l-1])
                        inst.dihedrals.append(dih)
                elif current_section == 'system':
                    inst.title = line
                elif current_section == 'defaults':
                    words = line.split()
                    if len(words) < 4:
                        raise GromacsTopologyError('Too few fields in '
                                                   '[ defaults ]')
                    if fields[0] != '1':
                        warnings.warn('Unsupported nonbonded type; unknown '
                                      'functional', GromacsTopologyWarning)
                        inst.unknown_functional = True
                    if fields[1] != '2':
                        warnings.warn('Unsupported combining rule',
                                      GromacsTopologyWarning)
                        inst.unknown_functional = True
                    if fields[2].lower() == 'no':
                        warnings.warn('gen_pairs=no is not supported')
                        inst.unknown_functional = True
                    inst._fudgeLJ = float(words[3])
                    inst._fudgeQQ = float(words[4])
                elif current_section == 'bondtypes':
                    words = line.split()
                    r = float(words[3]) * u.nanometers
                    k = (float(words[4]) / 2) * (
                            u.kilojoules_per_mole / u.nanometers**2)
                    if words[2] != '1':
                        warnings.warn('bondtypes funct != 1; unknown '
                                      'functional', GromacsTopologyWarning)
                        unst.unknown_functional = True
                    ptype = BondType(k, r)
                    params.bond_types[(words[0], words[1])] = ptype
                    params.bond_types[(words[1], words[0])] = ptype
                elif current_section == 'angletypes':
                    words = line.split()
                    theta = float(words[4]) * u.degrees
                    k = (float(words[5]) / 2) * (
                            u.kilojoules_per_mole / u.radians**2)
                    if words[2] != '1' and words[2] != '5':
                        warnings.warn('angletypes funct != 1; unknown '
                                      'functional', GromacsTopologyWarning)
                        inst.unknown_functional = True
                    if words[2] == '5':
                        # Contains the angle with urey-bradley
                        ub0 = float(words[6])
                        cub = float(words[7])
                        if cub == 0:
                            ub = NoUreyBradley
                        else:
                            ub0 *= u.nanometers
                            cub *= u.kilojoules_per_mole / u.nanometers**2
                            ub = BondType(cub, ub0)
                            params.urey_bradley_types[(words[0], words[2])] = ub
                            params.urey_bradley_types[(words[2], words[0])] = ub
                    ptype = AngleType(k, theta)
                    params.angle_types[(words[0], words[1], words[2])] = ptype
                    params.angle_types[(words[2], words[1], words[0])] = ptype
                elif current_section == 'dihedraltypes':
                    words = line.split()
                    replace = False
                    dtype = 'normal'
                    a1, a2, a3, a4 = words[:4]
                    if words[4] == '1':
                        pass
                    if words[4] == '4':
                        replace = True
                    elif words[4] == '9':
                        pass
                    elif words[4] == '2':
                        replace = True
                        dtype = 'improper'
                    elif words[4] == '5':
                        dtype = 'rbtorsion'
                    else:
                        warnings.warn('dihedraltypes funct not supported',
                                      GromacsTopologyWarning)
                        inst.unknown_functional = True
                    # Do the proper types
                    if dtype == 'normal':
                        phase = float(words[5]) * u.degrees
                        phi_k = float(words[6]) * u.kilojoules_per_mole
                        per = int(words[7])
                        dt = DihedralType(phi_k, per, phase)
                        key = (words[0], words[1], words[2], words[3])
                        rkey = (words[3], words[2], words[1], words[0])
                        if replace or not key in params.dihedral_types:
                            dtl = DihedralTypeList()
                            dtl.append(dt)
                            params.dihedral_types[key] = dtl
                            params.dihedral_types[rkey] = dtl
                        else:
                            params.dihedral_types[key].append(dt)
                    elif dtype == 'improper':
                        theta = float(words[5])*u.degrees
                        k = float(words[6])*u.kilojoules_per_mole/u.radians**2
                        a1, a2, a3, a4 = words[:4]
                        ptype = ImproperType(k, theta)
                        params.improper_types[(a1, a2, a3, a4)] = ptype
                    elif dtype == 'rbtorsion':
                        a1, a2, a3, a4 = words[:4]
                        c0, c1, c2, c3, c4, c5 = [float(x) for x in words[5:11]]
                        ptype = RBTorsionType(c0, c1, c2, c3, c4, c5)
                        params.rb_torsion_types[(a1, a2, a3, a4)] = ptype
                        params.rb_torsion_types[(a4, a3, a2, a1)] = ptype
        finally:
            f.close()
