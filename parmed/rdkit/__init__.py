"""
A package dealing with RDKit integration.
"""
from parmed.rdkit.rdkit import RDKit

__all__ = ['load_rdkit', 'from_smiles', 'from_sdf']

load_rdkit = RDKit.load
from_smiles = RDKit.from_smiles
from_sdf = RDKit.from_sdf
