"""
This module contains the necessary machinery to add the LENNARD_JONES_CCOEF to a
topology file given a list of atomic polarizabilities
"""
from __future__ import division

from ParmedTools.exceptions import LJ12_6_4Error, DuplicateParamWarning
import warnings

WATER_POL = 1.444 # Polarizability of water

DEFAULT_C4_PARAMS = dict(
        TIP3P = dict(Be=186.5, Cu=290.9, Ni=212.8, Zn=231.6, Co=209.7,
                     Cr=136.8, Fe=163.0, Mg=132.9, V=195.7, Mn=146.1,
                     Hg=288.8, Cd=185.6, Ca=87.3, Sn=187.9, Sr=82.7,
                     Ba=71.9),
        TIP4PEW = dict(Be=228.5, Cu=339.2, Ni=259.2, Zn=272.3, Co=252.8,
                       Cr=177.4, Fe=201.1, Mg=180.5, V=244.8, Mn=192.3,
                       Hg=335.2, Cd=233.7, Ca=128.0, Sn=231.4, Sr=118.9,
                       Ba=112.5),
        SPCE = dict(Be=188.1, Cu=304.4, Ni=205.2, Zn=231.2, Co=209.2,
                    Cr=131.2, Fe=155.4, Mg=122.2, V=206.6, Mn=154.9,
                    Hg=300.2, Cd=198.8, Ca=89.0, Sn=201.1, Sr=96.3,
                    Ba=85.8),
)

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

    if tunfactor is None:
        tunfactor = 1.0

    print "***********************************************************"
    # Determine which atom type was treated as the center metal ion
    typelist = dict()
    for i in mask.Selected():
        typeidx = parm.parm_data['ATOM_TYPE_INDEX'][i]
        if typeidx in typelist: continue
        typelist[typeidx] = parm.atom_list[i].atomic_number
        print "The selected metal ion is", \
                pt.Element[parm.atom_list[i].atomic_number]
    types = typelist.keys()
    types.sort()

    # Determine the C4 term between the atom type of center metal ion 
    # and every atom type in the prmtop file 

    # 1. Get the list of AMBER_ATOM_TYPE and ATOM_TYPE_INDEX for all 
    # the atoms in prmtop file

    amberatomtypelist1 = []
    atomtypeindexlist1 = [] 
    for i in range(0, parm.pointers['NATOM']):
        amberatomtypelist1.append(parm.parm_data['AMBER_ATOM_TYPE'][i])
        atomtypeindexlist1.append(parm.parm_data['ATOM_TYPE_INDEX'][i])

    # 2. Have the represetative AMBER_ATOM_TYPE for the each certain 
    # ATOM_TYPE_INDEX

    amberatomtypelist2 = []
    atomtypeindexlist2 = []
    for i in range(0, parm.pointers['NATOM']):
        if not atomtypeindexlist1[i] in atomtypeindexlist2:
            amberatomtypelist2.append(amberatomtypelist1[i])
            atomtypeindexlist2.append(atomtypeindexlist1[i])

    #3.Generate the C4 term for each atom type
    result = [0.0 for i in range(len(parm.parm_data['LENNARD_JONES_ACOEF']))]

 
    ntypes = parm.pointers['NTYPES']
    for typ in types:
#       print "The biggest ATOM_TYPE_INDEX =", ntypes
#       print "The selected +2 metal ion has the ATOM_TYPE_INDEX =",typ
        i = typ - 1

#       print "***********************************************************"
#       print "Here are the atom types which have been added C4 term:"
        # for the situation j < i
        for j in range(0 , typ):
            atomtypej = amberatomtypelist2[j]
            idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*j+i]-1
            try:
                c4 = c4list[ pt.Element[typelist[typ]] ]
                pol = pollist[atomtypej]
#               print ('ATOM_TYPE_INDEX = %d; AMBER_ATOM_TYPE=%s; '
#                      'Polarizability=%s' %(j+1, atomtypej, pol))
            except KeyError:
                raise LJ12_6_4Error("Could not find parameters for "
                                    "ATOM_TYPE %s" % atomtypej )

            if (atomtypej=="OW"):
                result[idx] = c4 
            else:   
                result[idx] = c4 / WATER_POL * pol * tunfactor 
               
            # for the situation i =< j 
        for j in range(typ, parm.ptr('ntypes')):
            atomtypej = amberatomtypelist2[j]
            idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
            try:
                c4 = c4list[ pt.Element[typelist[typ]] ]
                pol = pollist[atomtypej]
#               print ('ATOM_TYPE_INDEX = %d; AMBER_ATOM_TYPE= %s; '
#                      'Polarizability=%s' %(j+1, atomtypej, pol))
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
