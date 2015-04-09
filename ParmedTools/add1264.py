"""
This module contains the necessary machinery to add the LENNARD_JONES_CCOEF to a
topology file given a list of atomic polarizabilities
"""
from __future__ import division

from ParmedTools.exceptions import LJ12_6_4Error, DuplicateParamWarning
import warnings

WATER_POL = 1.444 # Polarizability of water

DEFAULT_C4_PARAMS = {
        'TIP3P' : {'Li1': 27.0, 'Na1': 0.0, 'K1': 8.0, 'Rb1': 0.0, 'Cs1': 2.0, 'Tl1': 50.0,
                   'Cu1': 7.0, 'Ag1': 83.0, 'F-1': -27.0, 'Cl-1': -38.0, 'Br-1': -39.0,
                   'I-1': -45.0, 'Be2': 186.5, 'Cu2': 290.9, 'Ni2': 212.8, 'Zn2': 231.6,
                   'Co2': 209.7, 'Cr2': 136.8, 'Fe2': 163.0, 'Mg2': 132.9, 'V2': 195.7,
                   'Mn2': 146.1, 'Hg2': 288.8, 'Cd2': 185.6, 'Ca2': 87.3, 'Sn2': 187.9,
                   'Sr2': 82.7, 'Ba2': 71.9, 'Al3': 399.0, 'Fe3': 428.0, 'Cr3': 258.0,
                   'In3': 347.0, 'Tl3': 456.0, 'Y3': 216.0, 'La3': 152.0, 'Ce3': 230.0,
                   'Pr3': 264.0, 'Nd3': 213.0, 'Sm3': 230.0, 'Eu3': 259.0, 'Gd3': 198.0,
                   'Tb3': 235.0, 'Dy3': 207.0, 'Er3': 251.0, 'Tm3': 282.0, 'Lu3': 249.0,
                   'Hf4': 827.0, 'Zr4': 761.0, 'Ce4': 706.0, 'U4': 1034.0, 'Pu4': 828.0,
                   'Th4': 512.0},
        'TIP4PEW' : {'Li1': 36.0, 'Na1': 9.0, 'K1': 24.0, 'Rb1': 13.0, 'Cs1': 16.0, 'Tl1': 65.0,
                     'Cu1': 21.0, 'Ag1': 94.0, 'F-1': -67.0, 'Cl-1': -66.0, 'Br-1': -68.0,
                     'I-1': -62.0, 'Be2': 228.5, 'Cu2': 339.2, 'Ni2': 259.2, 'Zn2': 272.3,
                     'Co2': 252.8, 'Cr2': 177.4, 'Fe2': 201.1, 'Mg2': 180.5, 'V2': 244.8,
                     'Mn2': 192.3, 'Hg2': 335.2, 'Cd2': 233.7, 'Ca2': 128.0, 'Sn2': 231.4,
                     'Sr2': 118.9, 'Ba2': 112.5, 'Al3': 488.0, 'Fe3': 519.0, 'Cr3': 322.0,
                     'In3': 425.0, 'Tl3': 535.0, 'Y3': 294.0, 'La3': 243.0, 'Ce3': 315.0,
                     'Pr3': 348.0, 'Nd3': 297.0, 'Sm3': 314.0, 'Eu3': 345.0, 'Gd3': 280.0,
                     'Tb3': 313.0, 'Dy3': 298.0, 'Er3': 328.0, 'Tm3': 356.0, 'Lu3': 331.0,
                     'Hf4': 956.0, 'Zr4': 895.0, 'Ce4': 835.0, 'U4': 1183.0, 'Pu4': 972.0,
                     'Th4': 625.0},
        'SPCE' :  {'Li1': 33.0, 'Na1': 6.0, 'K1': 19.0, 'Rb1': 7.0, 'Cs1': 12.0, 'Tl': 61.0,
                  'Cu1': 9.0, 'Ag1': 92.0, 'F-1': -53.0, 'Cl-1': -55.0, 'Br-1': -51.0,
                  'I-1': -51.0, 'Be2': 188.1, 'Cu2': 304.4, 'Ni2': 205.2, 'Zn2': 231.2,
                  'Co2': 209.2, 'Cr2': 131.2, 'Fe2': 155.4, 'Mg2': 122.2, 'V2': 206.6,
                  'Mn2': 154.9, 'Hg2': 300.2, 'Cd2': 198.8, 'Ca2': 89.0, 'Sn2': 201.1,
                  'Sr2': 96.3, 'Ba2': 85.8, 'Al3': 406.0, 'Fe3': 442.0, 'Cr3': 254.0,
                  'In3': 349.0, 'Tl3': 455.0, 'Y3': 209.0, 'La3': 165.0, 'Ce3': 242.0,
                  'Pr3': 272.0, 'Nd3': 235.0, 'Sm3': 224.0, 'Eu3': 273.0, 'Gd3': 186.0,
                  'Tb3': 227.0, 'Dy3': 206.0, 'Er3': 247.0, 'Tm3': 262.0, 'Lu3': 247.0,
                  'Hf4': 810.0, 'Zr4': 760.0, 'Ce4': 694.0, 'U4': 1043.0, 'Pu4': 828.0,
                  'Th4': 513.0}
}

def params1264(parm, mask, c4file, watermodel, polfile, tunfactor):
   
    from chemistry import periodic_table as pt

    global DEFAULT_C4_PARAMS, WATER_POL

    try:
        pollist = _get_params(polfile)
    except ValueError:
        raise LJ12_6_4Error('Bad polarizability file %s. Expected a file with '
                            '2 columns: <Amber Atom Type> <Polarizability>' %
                            polfile)
   
    if c4file is None:
        c4list = DEFAULT_C4_PARAMS[watermodel]
    else:
        try:
            c4list = _get_params(c4file)
        except ValueError:
            raise LJ12_6_4Error('Bad C4 parameter file %s. Expected a file '
                                'with 2 columns: <Atom Element> '
                                '<C4 Parameter>' % c4file)


    print "***********************************************************"
    # Determine which atom type was treated as the center metal ion
    typelist = dict()
    for i in mask.Selected():
        typeidx = parm.parm_data['ATOM_TYPE_INDEX'][i]
        chgi = parm.parm_data['CHARGE'][i]
        if typeidx in typelist: continue
        typelist[typeidx] = (parm.atoms[i].atomic_number, int(chgi))
        print "The selected metal ion is", \
                pt.Element[parm.atoms[i].atomic_number]
    types = typelist.keys()
    types.sort()

    # Determine the C4 term between the atom type of center metal ion 
    # and every atom type in the prmtop file 

    # 1. Get the list of AMBER_ATOM_TYPE and ATOM_TYPE_INDEX for all 
    # the atoms in prmtop file

    amberatomtypelist1 = parm.parm_data['AMBER_ATOM_TYPE']
    atomtypeindexlist1 = parm.parm_data['ATOM_TYPE_INDEX']

    # 2. Have the representative AMBER_ATOM_TYPE for the each certain
    # ATOM_TYPE_INDEX

    amberatomtypelist2 = []
    atomtypeindexlist2 = []
    for i in xrange(parm.pointers['NATOM']):
        if not atomtypeindexlist1[i] in atomtypeindexlist2:
            amberatomtypelist2.append(amberatomtypelist1[i])
            atomtypeindexlist2.append(atomtypeindexlist1[i])

    #3.Generate the C4 term for each atom type
    result = [0.0 for i in xrange(len(parm.parm_data['LENNARD_JONES_ACOEF']))]

 
    ntypes = parm.pointers['NTYPES']
    for typ in types:
#       print "The biggest ATOM_TYPE_INDEX =", ntypes
#       print "The selected +2 metal ion has the ATOM_TYPE_INDEX =",typ
        i = typ - 1

#       print "***********************************************************"
#       print "Here are the atom types which have been added C4 term:"
        # for the situation j < i
        for j in xrange(typ):
            atomtypej = amberatomtypelist2[j]
            jj = atomtypeindexlist2[j] - 1
            idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*jj+i]-1
            try:
                c4 = c4list[ pt.Element[typelist[typ][0]] + str(typelist[typ][1]) ]
                pol = pollist[atomtypej]
#               print ('ATOM_TYPE_INDEX = %d; AMBER_ATOM_TYPE=%s; '
#                      'Polarizability=%s' %(jj+1, atomtypej, pol))
            except KeyError:
                raise LJ12_6_4Error("Could not find parameters for "
                                    "ATOM_TYPE %s" % atomtypej )

            if (atomtypej=="OW"):
                result[idx] = c4 
            else:   
                result[idx] = c4 / WATER_POL * pol * tunfactor 
               
            # for the situation i =< j 
        for j in xrange(typ, parm.ptr('ntypes')):
            atomtypej = amberatomtypelist2[j]
            jj = atomtypeindexlist2[j] - 1
            idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+jj] - 1
            try:
                c4 = c4list[ pt.Element[typelist[typ][0]] + str(typelist[typ][1]) ]
                pol = pollist[atomtypej]
#               print ('ATOM_TYPE_INDEX = %d; AMBER_ATOM_TYPE= %s; '
#                      'Polarizability=%s' %(jj+1, atomtypej, pol))
            except KeyError:
                raise LJ12_6_4Error("Could not find parameters for "
                                    "ATOM_TYPE %s" % atomtypej)

            if (atomtypej=="OW"):
                result[idx] = c4
            else:
                result[idx] = c4 / WATER_POL * pol * tunfactor
#       print "***********************************************************"
    return result


def _get_params(fname):
    params = dict()
    f = open(fname, 'r')
    for line in f:
        atomtype, param = line.split()[:2]
        param = float(param)
        if atomtype in params and abs(param - params[atomtype]) > 0.0001:
            warnings.warn('Atom type %s has multiple parameters in %s.' %
                          (atomtype, fname), DuplicateParamWarning)
        params[atomtype] = param

    return params
