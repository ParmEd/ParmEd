#!/usr/bin/env python

import openmoltools.utils

in_prmtop = 'toluene_cyclohexane_10_500.prmtop'
in_crd = 'toluene_cyclohexane_10_500.inpcrd'

openmoltools.utils.convert_via_acpype( 'mixture', in_prmtop, in_crd, out_top = 'toluene_cyclohexane_10_500_acpype.top', out_gro = 'toluene_cyclohexane_10_500_acpype.gro' )
openmoltools.utils.amber_to_gromacs( 'mixture', in_prmtop, in_crd, out_top = 'toluene_cyclohexane_10_500_parmed.top', out_gro = 'toluene_cyclohexane_10_500_parmed.gro' )
