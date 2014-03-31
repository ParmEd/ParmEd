""" Changes the LJ coefficient terms only for a specific type-pair """

def ChLJPair(parm, atom_1, atom_2, rmin, eps, one_4=False):
   
    if one_4:
        key = 'LENNARD_JONES_14'
    else:
        key = 'LENNARD_JONES'

    # Make sure that atom type 1 comes first
    a1, a2 = sorted([atom_1, atom_2])
    ntypes = parm.pointers['NTYPES']
   
    # Find the atom1 - atom2 interaction (adjusting for indexing from 0)
    term_idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*(a1-1)+a2-1] - 1
   
    # Now change the ACOEF and BCOEF arrays, assuming pre-combined values
    parm.parm_data['%s_ACOEF' % key][term_idx] = eps * rmin**12
    parm.parm_data['%s_BCOEF' % key][term_idx] = 2 * eps * rmin**6
