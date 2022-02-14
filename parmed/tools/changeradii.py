from parmed.tools.exceptions import ChangeRadiiError

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def bondi(parm):
    """ Sets the bondi radii """
    for i, atom in enumerate(parm.atoms):
        # Radius of C atom depends on what type it is
        if atom.atomic_number == 6:
            if atom.type.startswith('C1') and atom.mass > 13.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C2') and atom.mass > 14.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C3') and atom.mass > 15.0:
                atom.solvent_radius = 2.2
            else:
                atom.solvent_radius = 1.7
        # All other elements have fixed radii for all types/partners
        elif atom.atomic_number == 1:
            atom.solvent_radius = 1.2
        elif atom.atomic_number == 7:
            atom.solvent_radius = 1.55
        elif atom.atomic_number == 8:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 9:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 14:
            atom.solvent_radius = 2.1
        elif atom.atomic_number == 15:
            atom.solvent_radius = 1.85
        elif atom.atomic_number == 16:
            atom.solvent_radius = 1.8
        elif atom.atomic_number == 17:
            atom.solvent_radius = 1.7
        else:
            atom.solvent_radius = 1.5

    try:
        parm.parm_data['RADIUS_SET'][0] = 'Bondi radii (bondi)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def amber6(parm):
    """ Sets the amber6 radii """
    for i, atom in enumerate(parm.atoms):
        # Radius of H atom depends on element it is bonded to
        bondeds = list(atom.bond_partners)
        if atom.atomic_number == 1:
            if bondeds[0].atomic_number == 6: # carbon
                atom.solvent_radius = 1.3
            elif bondeds[0].atomic_number in (8, 16): # oxygen or sulfur
                atom.solvent_radius = 0.8
            else: # anything else
                atom.solvent_radius = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if atom.type.startswith('C1') and atom.mass > 13.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C2') and atom.mass > 14.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C3') and atom.mass > 15.0:
                atom.solvent_radius = 2.2
            else:
                atom.solvent_radius = 1.7
        # All other elements have fixed solvent_radius for all types/partners
        elif atom.atomic_number == 7:
            atom.solvent_radius = 1.55
        elif atom.atomic_number == 8:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 9:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 14:
            atom.solvent_radius = 2.1
        elif atom.atomic_number == 15:
            atom.solvent_radius = 1.85
        elif atom.atomic_number == 16:
            atom.solvent_radius = 1.8
        elif atom.atomic_number == 17:
            atom.solvent_radius = 1.7
        else:
            atom.solvent_radius = 1.5
    try:
        parm.parm_data['RADIUS_SET'][0] = 'amber6 modified Bondi radii (amber6)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi(parm):
    """ Sets the mbondi radii """
    for i, atom in enumerate(parm.atoms):
        # Radius of H atom depends on element it is bonded to
        if atom.atomic_number == 1:
            bondeds = list(atom.bond_partners)
            if bondeds[0].atomic_number in (6, 7): # C or N
                atom.solvent_radius = 1.3
            elif bondeds[0].atomic_number in (8, 16): # O or S
                atom.solvent_radius = 0.8
            else:
                atom.solvent_radius = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if atom.type.startswith('C1') and atom.mass > 13.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C2') and atom.mass > 14.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C3') and atom.mass > 15.0:
                atom.solvent_radius = 2.2
            else:
                atom.solvent_radius = 1.7
        # All other elements have fixed solvent_radius for all types/partners
        elif atom.atomic_number == 7:
            atom.solvent_radius = 1.55
        elif atom.atomic_number == 8:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 9:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 14:
            atom.solvent_radius = 2.1
        elif atom.atomic_number == 15:
            atom.solvent_radius = 1.85
        elif atom.atomic_number == 16:
            atom.solvent_radius = 1.8
        elif atom.atomic_number == 17:
            atom.solvent_radius = 1.7
        else:
            atom.solvent_radius = 1.5
    try:
        parm.parm_data['RADIUS_SET'][0] = 'modified Bondi radii (mbondi)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi_pb2(parm):
    """
    Sets the mbondi_pb2 radii

    Use ropt values for halogens in Table 3 (https://pubs.acs.org/doi/full/10.1021/acs.jctc.9b00106)

                     pb2
    halogen	 rstd	 MAEstd	 MAEopt (ropt)
    Cl	     1.70    1.544	 1.529  (1.76)
    Br	     1.85    0.792	 0.743  (1.97)
    I	     1.98    0.905	 0.858  (2.09)

    These ropt values have been obtained using the "pb2" routine described
    in https://pubs.acs.org/doi/full/10.1021/acs.jctc.9b00106

    setup	radii set(s)	ΔGnonpolar	      γ	        β        inp value in gmx_MMPBSA
    pb2	    mbondi	        γ * SASA + β	  0.00720	0.0000   1

    """
    for i, atom in enumerate(parm.atoms):
        # Radius of H atom depends on element it is bonded to
        if atom.atomic_number == 1:
            bondeds = list(atom.bond_partners)
            if bondeds[0].atomic_number in (6, 7): # C or N
                atom.solvent_radius = 1.3
            elif bondeds[0].atomic_number in (8, 16): # O or S
                atom.solvent_radius = 0.8
            else:
                atom.solvent_radius = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if atom.type.startswith('C1') and atom.mass > 13.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C2') and atom.mass > 14.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C3') and atom.mass > 15.0:
                atom.solvent_radius = 2.2
            else:
                atom.solvent_radius = 1.7
        # All other elements have fixed solvent_radius for all types/partners
        elif atom.atomic_number == 7:
            atom.solvent_radius = 1.55
        elif atom.atomic_number == 8:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 9:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 14:
            atom.solvent_radius = 2.1
        elif atom.atomic_number == 15:
            atom.solvent_radius = 1.85
        elif atom.atomic_number == 16:
            atom.solvent_radius = 1.8
        # Cl
        elif atom.atomic_number == 17:
            atom.solvent_radius = 1.76
        # Br
        elif atom.atomic_number == 35:
            atom.solvent_radius = 1.97
        # I
        elif atom.atomic_number == 53:
            atom.solvent_radius = 2.09
        else:
            atom.solvent_radius = 1.5
    try:
        parm.parm_data['RADIUS_SET'][0] = 'modified Bondi radii with ropt for halogens (mbondi_pb2)'
    except AttributeError:
        pass
    _screen1(parm)


def mbondi_pb3(parm):
    """
    Sets the mbondi_pb3 radii

    Use ropt values for halogens in Table 3 (https://pubs.acs.org/doi/full/10.1021/acs.jctc.9b00106)

                     pb3
    halogen	 rstd	 MAEstd	 MAEopt (ropt)
    Cl	     1.70    1.358	 1.236 (2.20)
    Br	     1.85    0.770	 0.681 (2.04)
    I	     1.98    0.769	 0.697 (2.19)

    These ropt values have been obtained using the "pb3" routine described
    in https://pubs.acs.org/doi/full/10.1021/acs.jctc.9b00106

    setup	radii set(s)	ΔGnonpolar	                γ	     β        inp value in gmx_MMPBSA
    pb3	    mbondi + Rmin	γ * SAV + β + ΔGdispersion	0.03780	 –0.5692  2

    The mbondi radii set is used in the ΔGpolar calculation, while the value corresponding to Rmin in GAFF is
    used for the calculation of ΔGnonpolar

    """
    for i, atom in enumerate(parm.atoms):
        # Radius of H atom depends on element it is bonded to
        if atom.atomic_number == 1:
            bondeds = list(atom.bond_partners)
            if bondeds[0].atomic_number in (6, 7): # C or N
                atom.solvent_radius = 1.3
            elif bondeds[0].atomic_number in (8, 16): # O or S
                atom.solvent_radius = 0.8
            else:
                atom.solvent_radius = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if atom.type.startswith('C1') and atom.mass > 13.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C2') and atom.mass > 14.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C3') and atom.mass > 15.0:
                atom.solvent_radius = 2.2
            else:
                atom.solvent_radius = 1.7
        # All other elements have fixed solvent_radius for all types/partners
        elif atom.atomic_number == 7:
            atom.solvent_radius = 1.55
        elif atom.atomic_number == 8:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 9:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 14:
            atom.solvent_radius = 2.1
        elif atom.atomic_number == 15:
            atom.solvent_radius = 1.85
        elif atom.atomic_number == 16:
            atom.solvent_radius = 1.8
        # Cl
        elif atom.atomic_number == 17:
            atom.solvent_radius = 2.20
        # Br
        elif atom.atomic_number == 35:
            atom.solvent_radius = 2.04
        # I
        elif atom.atomic_number == 53:
            atom.solvent_radius = 2.19
        else:
            atom.solvent_radius = 1.5
    try:
        parm.parm_data['RADIUS_SET'][0] = 'modified Bondi radii with ropt for halogens (mbondi_pb3)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi2(parm):
    """ Sets the mbondi2 radii """
    for i, atom in enumerate(parm.atoms):
        # Radius of H atom depends on element it is bonded to
        if atom.atomic_number == 1:
            if atom.bond_partners[0].atomic_number == 7:
                atom.solvent_radius = 1.3
            else:
                atom.solvent_radius = 1.2
        # Radius of C atom depends on what type it is
        elif atom.atomic_number == 6:
            if atom.type.startswith('C1') and atom.mass > 13.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C2') and atom.mass > 14.0:
                atom.solvent_radius = 2.2
            if atom.type.startswith('C3') and atom.mass > 15.0:
                atom.solvent_radius = 2.2
            else:
                atom.solvent_radius = 1.7
        # All other elements have fixed solvent_radius for all types/partners
        elif atom.atomic_number == 7:
            atom.solvent_radius = 1.55
        elif atom.atomic_number == 8:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 9:
            atom.solvent_radius = 1.5
        elif atom.atomic_number == 14:
            atom.solvent_radius = 2.1
        elif atom.atomic_number == 15:
            atom.solvent_radius = 1.85
        elif atom.atomic_number == 16:
            atom.solvent_radius = 1.8
        elif atom.atomic_number == 17:
            atom.solvent_radius = 1.7
        else:
            atom.solvent_radius = 1.5
    try:
        parm.parm_data['RADIUS_SET'][0] = 'H(N)-modified Bondi radii (mbondi2)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mbondi3(parm):
    """ Sets mbondi3 radii """
    mbondi2(parm) # start from mbondi2 radii
    for i, atom in enumerate(parm.atoms):
        # Adjust OE (GLU), OD (ASP), and HH/HE (ARG)
        if atom.residue.name in ('GLU', 'ASP', 'GL4', 'AS4'):
            if atom.name.startswith('OE') or atom.name.startswith('OD'):
                atom.solvent_radius = 1.4
        elif atom.residue.name == 'ARG':
            if atom.name.startswith('HH') or atom.name.startswith('HE'):
                atom.solvent_radius = 1.17
        # Adjust carboxylate O radii on C-Termini. Don't just do the end
        # residue, since we can have C-termini in the middle as well
        # (i.e., 2-chain dimers)
        if atom.name == 'OXT':
            atom.solvent_radius = 1.4
            parm.atoms[i-1].solvent_radius = 1.4

    try:
        parm.parm_data['RADIUS_SET'][0] = \
                'ArgH and AspGluO modified Bondi2 radii (mbondi3)'
    except AttributeError:
        pass
    _screen1(parm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _screen1(parm):
    """ Applies the first set of screening parameters found in tleap source """
    for i, atom in enumerate(parm.atoms):
        if atom.atomic_number == 1:
            atom.screen = 0.85
        elif atom.atomic_number == 6:
            atom.screen = 0.72
        elif atom.atomic_number == 7:
            atom.screen = 0.79
        elif atom.atomic_number == 8:
            atom.screen = 0.85
        elif atom.atomic_number == 9:
            atom.screen = 0.88
        elif atom.atomic_number == 15:
            atom.screen = 0.86
        elif atom.atomic_number == 16:
            atom.screen = 0.96
        else:
            atom.screen = 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _screen2(parm):
    """
    Applies the second set of screening parameters found in tleap source
    (unused as far as I can tell)
    """
    for i, atom in enumerate(parm.atoms):
        if atom.atomic_number == 1:
            atom.screen = 0.8461
        elif atom.atomic_number == 6:
            atom.screen = 0.9615
        elif atom.atomic_number == 7:
            atom.screen = 0.9343
        elif atom.atomic_number == 8:
            atom.screen = 1.0088
        elif atom.atomic_number == 11:
            atom.screen = 1.0000
        elif atom.atomic_number == 12:
            atom.screen = 1.0000
        elif atom.atomic_number == 15:
            atom.screen = 1.0700
        elif atom.atomic_number == 16:
            atom.screen = 1.1733
        else:
            atom.screen = 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _screen3(parm):
    """
    Applies the third and final set of screening parameters found in tleap
    source (unused as far as I can tell)
    """
    for i, atom in enumerate(parm.atoms):
        if atom.atomic_number == 1:
            atom.screen = 0.8846
        elif atom.atomic_number == 6:
            atom.screen = 0.9186
        elif atom.atomic_number == 7:
            atom.screen = 0.8733
        elif atom.atomic_number == 8:
            atom.screen = 0.8836
        elif atom.atomic_number == 11:
            atom.screen = 1.0000
        elif atom.atomic_number == 12:
            atom.screen = 1.0000
        elif atom.atomic_number == 15:
            atom.screen = 0.9604
        elif atom.atomic_number == 16:
            atom.screen = 0.9323
        else:
            atom.screen = 0.8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_call_method = dict(bondi=bondi, mbondi=mbondi, mbondi_pb2=mbondi_pb2, mbondi_pb3=mbondi_pb3, mbondi2=mbondi2,
                    mbondi3=mbondi3, amber6=amber6, charmm_radii=charmm_radii)

def ChRad(parm, radii_set):
    global _call_method
    if radii_set not in _call_method:
        raise ChangeRadiiError("You must choose from %s radii sets" %
                               ', '.join(_call_method.keys()))
    _call_method[radii_set](parm)
