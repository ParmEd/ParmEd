""" 
Code to modify the molecule section of the topology file. It will combine
adjacent molecules. The code to determine molecularity is in readparm.py
"""
def combineMolecules(parm, molnum):
   """ Combines 2 molecules into a single one """
   if molnum >= parm.parm_data['SOLVENT_POINTERS'][2]:
      parm.parm_data['SOLVENT_POINTERS'][1] -= 1
      parm.parm_data['ATOMS_PER_MOLECULE'][molnum-1] += \
         parm.parm_data['ATOMS_PER_MOLECULE'][molnum]
      del parm.parm_data['ATOMS_PER_MOLECULE'][molnum]
   elif molnum < parm.parm_data['SOLVENT_POINTERS'][2]:
      parm.parm_data['SOLVENT_POINTERS'][1] -= 1
      parm.parm_data['SOLVENT_POINTERS'][2] -= 1
      parm.parm_data['ATOMS_PER_MOLECULE'][molnum-1] += \
         parm.parm_data['ATOMS_PER_MOLECULE'][molnum]
      del parm.parm_data['ATOMS_PER_MOLECULE'][molnum]
