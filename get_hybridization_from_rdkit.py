"""
The purpose of this script is to convert all of the CCD-generated molecules to RDKit in order to
extract the hybridization for each templated residue
"""
import json
from pathlib import Path
from rdkit import Chem
from pathlib import Path
from parmed.rdkit import RDKit
from parmed.modeller.standardtemplates import get_standard_residue_template_library

residues = get_standard_residue_template_library()

rdkit_map = dict()

for name, template in residues.items():
    mol = RDKit.to_mol(template)
    Chem.SanitizeMol(mol, Chem.SANITIZE_SETHYBRIDIZATION)
    rdkit_map[name] = mol

with Path("ccd_residue_templates.json").open("r") as f:
    data = json.load(f)

for res in data:
    if res["name"] not in rdkit_map:
        print(f"Residue {res['name']} is not in the RDKit map")
        continue
    assert len(res["atoms"]) == rdkit_map[res["name"]].GetNumAtoms()
    for atom_data, atom in zip(res["atoms"], rdkit_map[res["name"]].GetAtoms()):
        atom_data["hybridization"] = atom.GetHybridization()

with Path("ccd_residue_templates_with_hybridization.json").open("w") as f:
    json.dump(data, f)
