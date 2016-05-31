""" Standard residue templates for biomolecular residues """

import os
from parmed.amber.offlib import AmberOFFLibrary as _parser

StandardBiomolecularResidues = _parser.parse(
        os.path.join(os.path.split(__file__)[0], 'data',
                     'standard_residues.lib')
)
