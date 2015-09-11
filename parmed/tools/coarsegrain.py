"""
This module handles adding coarse graining data to topology files set up by
Luciana Capece. It has been implemented here in parmed by Jason Swails, but all
ideas were Luciana's. Here is a description of what this section does to
topology files:

ANGLE_FORCE_CONSTANT has been rigged to provide an integer index to the specific
parameter found in the given parameter file. This will be used to set up 4
different (new) prmtop sections: ANGLE_COEF_A, ANGLE_COEF_B, ANGLE_COEF_C, and
ANGLE_COEF_D. ANGLE_FORCE_CONSTANT now becomes useless, so I will add a comment
to the topology file saying such.

DIHEDRAL_FORCE_CONSTANT has been rigged in much the same way, creating 8 new
sections: DIHEDRAL_AMPLITUDE_1,2,3,4 and DIHEDRAL_PHASE_1,2,3,4.
DIHEDRAL_FORCE_CONSTANT, DIHEDRAL_PHASE, and DIHEDRAL_PERIODICITY become useless
and are commented as such (though not removed).

This module sets up ANGLE and DIHEDRAL classes to make printing them easier,
then parses the parameter file, then sets the respective sections of the prmtop
file.
"""
from parmed.utils.six.moves import range
import warnings

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Angle(object):
    """ Angle between 3 bonded atoms """
    def __init__(self, atom1, atom2, atom3, acoef, bcoef, ccoef, dcoef):
        self.atom1 = atom1.strip()
        self.atom2 = atom2.strip()
        self.atom3 = atom3.strip()
        self.acoef = float(acoef)
        self.bcoef = float(bcoef)
        self.ccoef = float(ccoef)
        self.dcoef = float(dcoef)

    def __eq__(self, other):
        """
        Defines if 2 angles are equal: the middle atom has to be the same, and
        the 2 outer ones have to be the same in either order, or they can be
        wildcards (X)
        """
        # First see if we are looking to match the angles in forward or reverse
        reversed_angles = False
        if self.atom2 != other.atom2: return False
        if self.atom1 != other.atom1:
            if self.atom1 != other.atom3:
                if self.atom1 != 'X' and other.atom1 != 'X':
                    if other.atom3 == 'X':
                        reversed_angles = True
                    else:
                        return False
                # end if self.atom1 != 'X' and other.atom1 != 'X'
            else:
                reversed_angles = True
        # Now check the rest (non-reversed case first)
        if not reversed_angles:
            if (self.atom3 != other.atom3 and self.atom3 != 'X' and
                other.atom3 != 'X'):
                return False
        else: # reversed angles
            if (self.atom3 != other.atom1 and self.atom3 != 'X' and
                other.atom1 != 'X'):
                return False
        # If we passed all these tests so far, then we must be the same!
        return True

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Dihedral(object):
    """ Dihedral between 4 bonded atoms """
    def __init__(self, atom1, atom2, atom3, atom4, 
                 ampl1, ampl2, ampl3, ampl4,
                 phase1, phase2, phase3, phase4):
        self.atom1, self.atom2 = atom1, atom2
        self.atom3, self.atom4 = atom3, atom4
        self.ampl1, self.ampl2 = ampl1, ampl2
        self.ampl3, self.ampl4 = ampl3, ampl4
        self.phase1, self.phase2 = phase1, phase2
        self.phase3, self.phase4 = phase3, phase4

    def __eq__(self, other):
        """
        2 Dihedrals are equal if they go in forward/reverse order and all atom
        types either match or we have wild cards (which can only be on the ends)
        """
        reversed_angles = False
        if self.atom2 != other.atom2 or self.atom3 != other.atom3:
            if self.atom2 != other.atom3 or self.atom3 != other.atom2:
                return False
            else:
                reversed_angles = True
        if not reversed_angles:
            if self.atom1 != other.atom1:
                if self.atom1 != 'X' and other.atom1 != 'X':
                    return False
            if self.atom4 != other.atom4:
                if self.atom4 != 'X' and other.atom4 != 'X':
                    return False
        else:
            if self.atom1 != other.atom4:
                if self.atom1 != 'X' and other.atom4 != 'X':
                    return False
            if self.atom4 != other.atom1:
                if self.atom4 != 'X' and other.atom4 != 'X':
                    return False
        # There must be no problems if we made it this far
        return True

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def addCoarseGrain(parm, param_file):
    from parmed.tools.exceptions import CoarseGrainError
    """ Adds coarse graining sections to the prmtop """
    # Angles and Dihedrals dictionaries for storing parameters parsed from
    # parameter file
    angle_params = {}
    dihedral_params = {}
    # First parse the params file, filling dictionaries of angles and dihedrals
    # whose keys are their integer index values
    params = open(param_file, 'r')
    reading_angles = False
    reading_dihedrals = False
    file_line = params.readline()
    while file_line:
        # Filter out comments
        if '#' in file_line:
            file_line = file_line[:file_line.index('#')]
        # Filter out blank lines
        file_line = file_line.strip()
        if len(file_line) == 0: 
            file_line = params.readline()
            continue
        # See if we start reading ANGLEs yet:
        if file_line[:5].upper() == 'ANGLE':
            reading_angles = True
            reading_dihedrals = False
            file_line = params.readline()
            continue
        # See if we start reading DIHEDRALs yet
        if file_line[:8].upper() == 'DIHEDRAL':
            reading_dihedrals = True
            reading_angles = False
            file_line = params.readline()
            continue
        # Pass over any INDEX header lines -- they're just there for readability
        if (reading_dihedrals or reading_angles) and file_line[:5] == 'INDEX':
            file_line = params.readline()
            continue
        if reading_angles:
            line_parts = file_line.split()
            try:
                idx = int(line_parts[0])
                acoef = float(line_parts[1])
                bcoef = float(line_parts[2])
                ccoef = float(line_parts[3])
                dcoef = float(line_parts[4])
                atom_types = line_parts[5].split('-')
                angle_params[idx] = Angle(atom_types[0], atom_types[1],
                                          atom_types[2], acoef, bcoef, ccoef,
                                          dcoef)
            except ValueError as err:
                raise CoarseGrainError(
                        'Unexpected format in Coarse Grain angles. Expected '
                        'different data type: %s. See format specification'
                        % err)
            except IndexError:
                raise CoarseGrainError(
                        'Unexpected format in Coarse Grain parameter file. '
                        'Expected more data fields on the line. See format '
                        'specification')
            file_line = params.readline()
            continue
        if reading_dihedrals:
            line_parts = file_line.split()
            try:
                idx = int(line_parts[0])
                ampl1 = float(line_parts[1])
                ampl2 = float(line_parts[2])
                ampl3 = float(line_parts[3])
                ampl4 = float(line_parts[4])
                atom_types = line_parts[5].split('-')
                file_line = params.readline().strip()
                line_parts = file_line.split()
                phase1 = float(line_parts[0])
                phase2 = float(line_parts[1])
                phase3 = float(line_parts[2])
                phase4 = float(line_parts[3])
                dihedral_params[idx] = Dihedral(atom_types[0], atom_types[1],
                                                atom_types[2], atom_types[3],
                                                ampl1, ampl2, ampl3, ampl4,
                                                phase1, phase2, phase3, phase4)
            except ValueError as err:
                raise CoarseGrainError(
                        'Unexpected format in Coarse Grain dihedrals. Expected '
                        'different data type: %s. See format specification'
                        % err)
            except IndexError:
                raise CoarseGrainError(
                        'Unexpected format in Coarse Grain parameter file. '
                        'Expected more data fields on the line. See format '
                        'specification')
            file_line = params.readline()
            continue
        warnings.warn('Line (%s) ignored in Coarse Grain parameter file' %
                      file_line)
        file_line = params.readline()
    # End reading while

    # Now let's add comments to our topology file!
    parm.parm_comments['ANGLE_FORCE_CONSTANT'].append(
                'This section is only an index and is not used for Coarse '
                'grained topologies')
    parm.parm_comments['DIHEDRAL_FORCE_CONSTANT'].append(
                'This section is only an index and is not used for Coarse '
                'grained topologies')
    parm.parm_comments['DIHEDRAL_PERIODICITY'].append(
                'This section is not used for Coarse grained topologies')
    parm.parm_comments['DIHEDRAL_PHASE'].append(
                'This section is not used for Coarse grained topologies')
   
    # Now let's add our new sections
    parm.add_flag('ANGLE_COEF_A','5E16.8',parm.ptr('numang'),
                 comments='A Coefficient for Coarse grained force field')
    parm.add_flag('ANGLE_COEF_B','5E16.8',parm.ptr('numang'),
                 comments='B Coefficient for Coarse grained force field')
    parm.add_flag('ANGLE_COEF_C','5E16.8',parm.ptr('numang'),
                 comments='C Coefficient for Coarse grained force field')
    parm.add_flag('ANGLE_COEF_D','5E16.8',parm.ptr('numang'),
                 comments='D Coefficient for Coarse grained force field')
    parm.add_flag('DIHEDRAL_AMPLITUDE_1','5E16.8',parm.ptr('nptra'),
                 comments='1st Dihedral Amplitude for coarse grained force '
                          'field')
    parm.add_flag('DIHEDRAL_AMPLITUDE_2','5E16.8',parm.ptr('nptra'),
                 comments='2nd Dihedral Amplitude for coarse grained force '
                          'field'
    )
    parm.add_flag('DIHEDRAL_AMPLITUDE_3','5E16.8',parm.ptr('nptra'),
                 comments='3rd Dihedral Amplitude for coarse grained force '
                          'field'
    )
    parm.add_flag('DIHEDRAL_AMPLITUDE_4','5E16.8',parm.ptr('nptra'),
                 comments='4th Dihedral Amplitude for coarse grained force '
                          'field'
    )
    parm.add_flag('DIHEDRAL_PHASE_1','5E16.8',parm.ptr('nptra'),
                 comments='1st Dihedral Phase for coarse grained force field')
    parm.add_flag('DIHEDRAL_PHASE_2','5E16.8',parm.ptr('nptra'),
                 comments='2nd Dihedral Phase for coarse grained force field')
    parm.add_flag('DIHEDRAL_PHASE_3','5E16.8',parm.ptr('nptra'),
                 comments='3rd Dihedral Phase for coarse grained force field')
    parm.add_flag('DIHEDRAL_PHASE_4','5E16.8',parm.ptr('nptra'),
                 comments='4th Dihedral Phase for coarse grained force field')

    for i in range(len(parm.parm_data['ANGLE_FORCE_CONSTANT'])):
        try: 
            index = int(parm.parm_data['ANGLE_FORCE_CONSTANT'][i])
            angl = angle_params[index]
        except KeyError:
            raise CoarseGrainError('Missing angle parameters for index %d' %
                                   int(index))
        parm.parm_data['ANGLE_COEF_A'][i] = angl.acoef
        parm.parm_data['ANGLE_COEF_B'][i] = angl.bcoef
        parm.parm_data['ANGLE_COEF_C'][i] = angl.ccoef
        parm.parm_data['ANGLE_COEF_D'][i] = angl.dcoef
   
    for i in range(len(parm.parm_data['DIHEDRAL_FORCE_CONSTANT'])):
        try:
            index = int(parm.parm_data['DIHEDRAL_FORCE_CONSTANT'][i])
            dihe = dihedral_params[index]
        except KeyError:
            raise CoarseGrainError('Missing dihedral parameters for index %d' %
                                   int(index))
        parm.parm_data['DIHEDRAL_AMPLITUDE_1'][i] = dihe.ampl1
        parm.parm_data['DIHEDRAL_AMPLITUDE_2'][i] = dihe.ampl2
        parm.parm_data['DIHEDRAL_AMPLITUDE_3'][i] = dihe.ampl3
        parm.parm_data['DIHEDRAL_AMPLITUDE_4'][i] = dihe.ampl4
        parm.parm_data['DIHEDRAL_PHASE_1'][i] = dihe.phase1
        parm.parm_data['DIHEDRAL_PHASE_2'][i] = dihe.phase2
        parm.parm_data['DIHEDRAL_PHASE_3'][i] = dihe.phase3
        parm.parm_data['DIHEDRAL_PHASE_4'][i] = dihe.phase4
