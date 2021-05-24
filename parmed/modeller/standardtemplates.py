""" Standard residue templates for biomolecular residues """

import os
from ..amber.offlib import AmberOFFLibrary

__all__ = ['StandardBiomolecularResidues']

StandardBiomolecularResidues = AmberOFFLibrary.parse(
    os.path.join(os.path.split(__file__)[0], 'data', 'standard_residues.lib')
)
