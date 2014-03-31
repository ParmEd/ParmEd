from ParmedTools.exceptions import ChangeRadiiError

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def bondi(parm):
    """ Sets the bondi radii """
    for i, atom in enumerate(parm.atom_list):
        # Radius of C atom depends on what type it is
        if atom.atomic_number == 6:
            if atom.attype.startswith('C1') and atom.mass > 13.0:
                parm.parm_data['RADII'][i] = 2.2
            if atom.attype.startswith('C2') and atom.mass > 14.0:
                parm.parm_data['RADII'][i] = 2.2
            if atom.attype.startswith('C3') and atom.mass > 15.0:
                parm.parm_data['RADII'][i] = 2.2
            else:
                parm.parm_data['RADII'][i] = 1.7
        # All other elements have fixed radii for all types/partners
        elif atom.atomic_number == 1: parm.parm_data['RADII'][i] = 1.2
        elif atom.atomic_number == 7: parm.parm_data['RADII'][i] = 1.55
        elif atom.atomic_number == 8: parm.parm_data['RADII'][i] = 1.5
        elif atom.atomic_number == 9: parm.parm_data['RADII'][i] = 1.5
        elif atom.atomic_number == 14: parm.parm_data['RADII'][i] = 2.1
        elif atom.atomic_number == 15: parm.parm_data['RADII'][i] = 1.85
        elif atom.atomic_number == 16: parm.parm_data['RADII'][i] = 1.8
        elif atom.atomic_number == 17: parm.parm_data['RADII'][i] = 1.5
        else: parm.parm_data['RADII'][i] = 1.5

    parm.parm_data['RADIUS_SET'][0] = 'Bondi radii (bondi)'
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def amber6(parm):
    """ Sets the amber6 radii """
    for i, atom in enumerate(parm.atom_list):
        # Radius of H atom depends on element it is bonded to
        bondeds = list(atom.bond_partners)
        if atom.atomic_number == 1:
            if bondeds[0].atomic_number == 6: # carbon
                parm.parm_data['RADII'][i] = 1.3
            elif bondeds[0].atomic_number in (8, 16): # oxygen or sulfur
                parm.parm_data['RADII'][i] = 0.8
            else: # anything else
                parm.parm_data['RADII'][i] = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if atom.attype.startswith('C1') and atom.mass > 13.0:
                parm.parm_data['RADII'][i] = 2.2
            if atom.attype.startswith('C2') and atom.mass > 14.0:
                parm.parm_data['RADII'][i] = 2.2
            if atom.attype.startswith('C3') and atom.mass > 15.0:
                parm.parm_data['RADII'][i] = 2.2
            else:
                parm.parm_data['RADII'][i] = 1.7
        # All other elements have fixed radii for all types/partners
        elif atom.atomic_number == 7: parm.parm_data['RADII'][i] = 1.55
        elif atom.atomic_number == 8: parm.parm_data['RADII'][i] = 1.5
        elif atom.atomic_number == 9: parm.parm_data['RADII'][i] = 1.5
        elif atom.atomic_number == 14: parm.parm_data['RADII'][i] = 2.1
        elif atom.atomic_number == 15: parm.parm_data['RADII'][i] = 1.85
        elif atom.atomic_number == 16: parm.parm_data['RADII'][i] = 1.8
        elif atom.atomic_number == 17: parm.parm_data['RADII'][i] = 1.5
        else: parm.parm_data['RADII'][i] = 1.5
   
    parm.parm_data['RADIUS_SET'][0] = 'amber6 modified Bondi radii (amber6)'
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi(parm):
    """ Sets the mbondi radii """
    for i, atom in enumerate(parm.atom_list):
        # Radius of H atom depends on element it is bonded to
        if atom.atomic_number == 1:
            bondeds = list(atom.bond_partners)
            if bondeds[0].atomic_number in (6, 7): # C or N
                parm.parm_data['RADII'][i] = 1.3
            elif bondeds[0].atomic_number in (8, 16): # O or S
                parm.parm_data['RADII'][i] = 0.8
            else:
                parm.parm_data['RADII'][i] = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if atom.attype.startswith('C1') and atom.mass > 13.0:
                parm.parm_data['RADII'][i] = 2.2
            if atom.attype.startswith('C2') and atom.mass > 14.0:
                parm.parm_data['RADII'][i] = 2.2
            if atom.attype.startswith('C3') and atom.mass > 15.0:
                parm.parm_data['RADII'][i] = 2.2
            else:
                parm.parm_data['RADII'][i] = 1.7
        # All other elements have fixed radii for all types/partners
        elif atom.atomic_number == 7: parm.parm_data['RADII'][i] = 1.55
        elif atom.atomic_number == 8:  parm.parm_data['RADII'][i] = 1.5
        elif atom.atomic_number == 9:  parm.parm_data['RADII'][i] = 1.5
        elif atom.atomic_number == 14: parm.parm_data['RADII'][i] = 2.1
        elif atom.atomic_number == 15: parm.parm_data['RADII'][i] = 1.85
        elif atom.atomic_number == 16: parm.parm_data['RADII'][i] = 1.8
        elif atom.atomic_number == 17: parm.parm_data['RADII'][i] = 1.5
        else: parm.parm_data['RADII'][i] = 1.5

    parm.parm_data['RADIUS_SET'][0] = 'modified Bondi radii (mbondi)'
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi2(parm):
    """ Sets the mbondi2 radii """
    for i, atom in enumerate(parm.atom_list):
        # Radius of H atom depends on element it is bonded to
        if atom.atomic_number == 1:
            if atom.bond_partners[0].atomic_number == 7:
                parm.parm_data['RADII'][i] = 1.3
            else:
                parm.parm_data['RADII'][i] = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if atom.attype.startswith('C1') and atom.mass > 13.0:
                parm.parm_data['RADII'][i] = 2.2
            if atom.attype.startswith('C2') and atom.mass > 14.0:
                parm.parm_data['RADII'][i] = 2.2
            if atom.attype.startswith('C3') and atom.mass > 15.0:
                parm.parm_data['RADII'][i] = 2.2
            else:
                parm.parm_data['RADII'][i] = 1.7
        # All other elements have fixed radii for all types/partners
        elif atom.atomic_number == 7: parm.parm_data['RADII'][i] = 1.55
        elif atom.atomic_number == 8: parm.parm_data['RADII'][i] = 1.5
        elif atom.atomic_number == 9: parm.parm_data['RADII'][i] = 1.5
        elif atom.atomic_number == 14: parm.parm_data['RADII'][i] = 2.1
        elif atom.atomic_number == 15: parm.parm_data['RADII'][i] = 1.85
        elif atom.atomic_number == 16: parm.parm_data['RADII'][i] = 1.8
        elif atom.atomic_number == 17: parm.parm_data['RADII'][i] = 1.5
        else: parm.parm_data['RADII'][i] = 1.5

    parm.parm_data['RADIUS_SET'][0] = 'H(N)-modified Bondi radii (mbondi2)'
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi3(parm):
    """ Sets mbondi3 radii """
    mbondi2(parm) # start from mbondi2 radii
    for i, atom in enumerate(parm.atom_list):
        # Adjust OE (GLU), OD (ASP), and HH/HE (ARG)
        if atom.residue.resname in ('GLU', 'ASP', 'GL4', 'AS4'):
            if parm.parm_data['ATOM_NAME'][i].startswith('OE') or \
                parm.parm_data['ATOM_NAME'][i].startswith('OD'):
                parm.parm_data['RADII'][i] = 1.4
        elif atom.residue.resname == 'ARG':
            if parm.parm_data['ATOM_NAME'][i].startswith('HH') or \
                parm.parm_data['ATOM_NAME'][i].startswith('HE'):
                parm.parm_data['RADII'][i] = 1.17
        # Adjust carboxylate O radii on C-Termini. Don't just do the end
        # residue, since we can have C-termini in the middle as well
        # (i.e., 2-chain dimers)
        if parm.parm_data['ATOM_NAME'][i] == 'OXT':
            parm.parm_data['RADII'][i] = 1.4
            parm.parm_data['RADII'][i-1] = 1.4

    parm.parm_data['RADIUS_SET'][0] = \
                'ArgH and AspGluO modified Bondi2 radii (mbondi3)'
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _screen1(parm):
    """ Applies the first set of screening parameters found in tleap source """
    for i, atom in enumerate(parm.atom_list):
        if atom.atomic_number == 1:
            parm.parm_data['SCREEN'][i] = 0.85
        elif atom.atomic_number == 6:
            parm.parm_data['SCREEN'][i] = 0.72
        elif atom.atomic_number == 7:
            parm.parm_data['SCREEN'][i] = 0.79
        elif atom.atomic_number == 8:
            parm.parm_data['SCREEN'][i] = 0.85
        elif atom.atomic_number == 9:
            parm.parm_data['SCREEN'][i] = 0.88
        elif atom.atomic_number == 15:
            parm.parm_data['SCREEN'][i] = 0.86
        elif atom.atomic_number == 16:
            parm.parm_data['SCREEN'][i] = 0.96
        else:
            parm.parm_data['SCREEN'][i] = 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _screen2(parm):
    """
    Applies the second set of screening parameters found in tleap source
    (unused as far as I can tell)
    """
    for i, atom in enumerate(parm.atom_list):
        if atom.atomic_number == 1:
            parm.parm_data['SCREEN'][i] = 0.8461
        elif atom.atomic_number == 6:
            parm.parm_data['SCREEN'][i] = 0.9615
        elif atom.atomic_number == 7:
            parm.parm_data['SCREEN'][i] = 0.9343
        elif atom.atomic_number == 8:
            parm.parm_data['SCREEN'][i] = 1.0088
        elif atom.atomic_number == 11:
            parm.parm_data['SCREEN'][i] = 1.0000
        elif atom.atomic_number == 12:
            parm.parm_data['SCREEN'][i] = 1.0000
        elif atom.atomic_number == 15:
            parm.parm_data['SCREEN'][i] = 1.0700
        elif atom.atomic_number == 16:
            parm.parm_data['SCREEN'][i] = 1.1733
        else:
            parm.parm_data['SCREEN'][i] = 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _screen3(parm):
    """
    Applies the third and final set of screening parameters found in tleap
    source (unused as far as I can tell)
    """
    for i, atom in enumerate(parm.atom_list):
        if atom.atomic_number == 1:
            parm.parm_data['SCREEN'][i] = 0.8846
        elif atom.atomic_number == 6:
            parm.parm_data['SCREEN'][i] = 0.9186
        elif atom.atomic_number == 7:
            parm.parm_data['SCREEN'][i] = 0.8733
        elif atom.atomic_number == 8:
            parm.parm_data['SCREEN'][i] = 0.8836
        elif atom.atomic_number == 11:
            parm.parm_data['SCREEN'][i] = 1.0000
        elif atom.atomic_number == 12:
            parm.parm_data['SCREEN'][i] = 1.0000
        elif atom.atomic_number == 15:
            parm.parm_data['SCREEN'][i] = 0.9604
        elif atom.atomic_number == 16:
            parm.parm_data['SCREEN'][i] = 0.9323
        else:
            parm.parm_data['SCREEN'][i] = 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def ChRad(parm, radii_set):
    if not radii_set in ["mbondi", "mbondi2", "mbondi3", "bondi", "amber6"]:
        raise ChangeRadiiError("You must choose from mbondi, mbondi2, mbondi3,"
                               " bondi, or amber6 radii sets!")
   
    eval("%s(parm)" % radii_set)
